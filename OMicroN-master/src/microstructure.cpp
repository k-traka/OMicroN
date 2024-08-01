/* microstructure.cpp
   OMicroN (optimising microstructures numerically) simulation program
   cpp-file containing Microstructure class implementation
*/

#include "microstructure.h"
#include <cstdlib> // For exit()
#include <iomanip> // For std::setprecision
#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem> // For std::filesystem
#include <numeric>
#include <algorithm>
#include <cmath> // For std::abs, std::exp etc.
#include <optional>
#include <tiffio.h> // For exporting TIFF images
#include "application.h"
#include <eigen3/Eigen/Dense>
#include <cstdint>

/**
 * @brief Constructor for the Microstructure class.
 *
 * Initializes the microstructure simulation with default or user-specified parameters. This function
 * sets up necessary data structures and prepares the class for simulation steps.
 *
 * @param inputFile The path to the input file containing user-defined settings.
 */

Microstructure::Microstructure(const UserSettings &userSettings)
    : mvUserSettings(userSettings), mvpStateVars(nullptr), mvpTCK(nullptr), mvOrientations()
{
    LOG_F(INFO, "Initializing microstructure for simulation");

    // Initialize internal variables of the class which change per simulation step
    mvStillSimulating = true;
    mvSimulationStep = 0;
    mvTime = 0.0;
    mvCheckICellsTime = 0;
    mvTimePassedFromPreviousOutput = 0;
    mvMaxDiffusivityInTimeStep = 0.0;
    mvStoreSoluteEvolutionTime = 0.;

    printMicrostructureParameters();

    if (mvUserSettings.mvTimeStep > 1.e-15)
    {
        mvTimeStep = mvUserSettings.mvTimeStep;
        LOG_F(INFO, "Setting time step from user inputs equal to %g", mvTimeStep);
    }

    // Initialize Orientations
    mvOrientations.SetOrientationParameters(mvUserSettings);

    if (mvUserSettings.mvIsRexAndGG == 0 && mvUserSettings.mvHasSoluteDiffusion == 0)
    {
        std::cerr << "You have chosen neither a rex&GG simulation nor a carbon/transformation related simulation. Check your input file." << std::endl;
        exit(-1);
    }

    // Initialize StateVars
    mvpStateVars = new StateVars(
        mvUserSettings.mvNx, mvUserSettings.mvNy, mvUserSettings.mvNz, mvUserSettings.mvDx,
        mvUserSettings.mvGridType, mvUserSettings.mvHasSoluteDiffusion, mvUserSettings.mvIsRexAndGG,
        mvUserSettings.mvAllowInterfaceMovementDuringPartitioning, !mvUserSettings.mvIsPartitioningHappeningInDiffusionStep,
        !mvUserSettings.mvIsSoluteSegregationHappeningInDiffusionStep);

    // Initialize ThermChemKin
    mvpTCK = new ThermChemKin(
        mvUserSettings.mvStartTemperature, mvUserSettings.mvAverageCarbon,
        mvUserSettings.mvGrainBoundaryMobilityProExponential, mvUserSettings.mvGrainBoundaryMobilityActivationEnergy,
        mvUserSettings.mvGrainBoundaryEnergy, mvUserSettings.mvPhaseBoundaryMobilityProExponential,
        mvUserSettings.mvPhaseBoundaryMobilityActivationEnergy, mvUserSettings.mvPhaseBoundaryEnergy,
        mvUserSettings.mvSoluteDiffusivityPreFactor, mvUserSettings.mvSoluteDiffusivityActivationEnergy);

    if (mvUserSettings.mvIncludeCarbonInterphasePartitioning)
    {
        if (!mvUserSettings.mvThermodynamicsDataFilename.empty())
        {
            mvpTCK->LoadChemicalPotentialsAndLocalEquilibriumXC(mvUserSettings.mvThermodynamicsDataFilename.c_str(), mvUserSettings.mvConsiderConstantIronAtomsPerCell);
        }
        else
        {
            LOG_F(INFO, "Didn't provide Gibbs energies file. Using default parameters A, B, C, D.");
        }

        double paramA, paramB, paramC, paramD;
        mvpTCK->GetEqMuCarbonRelatedParametersABCD(&paramA, &paramB, &paramC, &paramD);
        LOG_F(INFO, "Parameters for XcEqBCC and XCEqFCC: A = %g, B = %g, C = %g, D = %g", paramA, paramB, paramC, paramD);
    }

    ReadMicrostructureFile(mvUserSettings.mvStartFileName);
    ValidateOrientationSettingsAndCheckNonIndexedPixels();
    SetStateVariablesBasedOnSimAndGiveInfoForMicrostructuralState();
    InitInterfaceAndInterphaseCells();
    int MapType = 0;
    ExportColorCodedTIFF(mvUserSettings.mvOutputFolderPath + "/LatticeStructure_Init_Map.tiff", MapType);

    if (mvUserSettings.mvIsRexAndGG)
        WriteMicrostructureRX(mvUserSettings.mvOutputFolderPath + "/MicroData_Init.txt");
    else
    {
        MapType = 1;
        ExportColorCodedTIFF(mvUserSettings.mvOutputFolderPath + "/Carbon_Init_Map.tiff", MapType);
        ExportColorbar(mvUserSettings.mvOutputFolderPath + "/Carbon_Init_Colorbar.tiff");

        if (mvUserSettings.mvSoluteSegregationDislocations)
        {
            MapType = 2;
            ExportColorCodedTIFF(mvUserSettings.mvOutputFolderPath + "/CarbonTrapped_Init_Map.tiff", MapType);
            ExportColorbar(mvUserSettings.mvOutputFolderPath + "/CarbonTrapped_Init_Colorbar.tiff");
        }
    }
}

// Print Microstructure Parameters
void Microstructure::printMicrostructureParameters() const
{
    LOG_F(INFO, "I will read microstructure file = %s", mvUserSettings.mvStartFileName.c_str());
    LOG_F(INFO, "I will read thermodynamic data from file = %s", mvUserSettings.mvThermodynamicsDataFilename.c_str());
    LOG_F(INFO, "I will write all simulation outputs in = %s", mvUserSettings.mvOutputFolderPath.c_str());
}

// Find Index
std::optional<int> Microstructure::findIndex(const std::vector<int> &vectorToCheck, int id)
{
    auto it = std::find(vectorToCheck.begin(), vectorToCheck.end(), id);
    if (it != vectorToCheck.end())
    {
        return std::distance(vectorToCheck.begin(), it);
    }
    return std::nullopt;
}

void Microstructure::ReadMicrostructureFile(const std::string &startFileName)
{
    Eigen::Vector3f ea;
    float CI;
    FILE *fp;
    float x, y, z, deltax, rho, id;
    int i, j, k, index, NCellsToRead, OriId, min_id, max_id, RX, IsPixelToBeConsidered, IsCellGammaFiber, IsCellThetaFiber, IsCellComponent;
    std::vector<int> cellIdsRead;

    int phase;
    bool is2d;
    bool CheckedState = false;
    const float DDFactor = 1.e12;
    const float DxFactor = 1.e6;

    int IsOnlyAnOrientationSubsetConsidered = 0;

    fp = fopen(startFileName.c_str(), "r");
    if (!fp)
    {
        LOG_F(ERROR, "Error reading microstructure file %s", startFileName.c_str());
        exit(EXIT_FAILURE);
    }

    // read number of cells in file
    if (fscanf(fp, "%d \n", &NCellsToRead) != 1)
    {
        LOG_F(ERROR, "Error reading microstructure file at line 1: %s", startFileName.c_str());
        exit(EXIT_FAILURE);
    }

    LOG_F(INFO, "NCellsToRead %d", NCellsToRead);

    if (NCellsToRead != mvpStateVars->mvN)
    {
        LOG_F(ERROR, "In ReadMicrostructure %s number of cells in file not equal to user supplied system dimensions, user %d file %d", startFileName.c_str(), mvpStateVars->mvN, NCellsToRead);
        exit(EXIT_FAILURE);
    }

    LOG_F(INFO, "Reading %s data points z %s expected %d data points from %s", std::to_string(mvpStateVars->mvN).c_str(), std::to_string(mvpStateVars->mvNz).c_str(), std::to_string(NCellsToRead).c_str(), startFileName.c_str());

    if (IsOnlyAnOrientationSubsetConsidered != 1 && IsOnlyAnOrientationSubsetConsidered != 2 && IsOnlyAnOrientationSubsetConsidered != 3 && IsOnlyAnOrientationSubsetConsidered != 4 && IsOnlyAnOrientationSubsetConsidered != 5 && IsOnlyAnOrientationSubsetConsidered != 6)
        LOG_F(INFO, "All cells (based on orientation) participate in the simulation");

    // allocate grains
    is2d = (mvpStateVars->mvNz == 1);

    int nOris = mvpStateVars->mvN;

    // read cells in file
    deltax = DxFactor * mvpStateVars->mvDx;
    float deltay = deltax;
    if (mvpStateVars->mvGridType == 1)
        deltay *= pow(0.75, 0.5);

    i = j = k = 0;
    index = 0;

    for (int line = 0; line < mvpStateVars->mvN; line++)
    {
        if (IsOnlyAnOrientationSubsetConsidered == 1 || IsOnlyAnOrientationSubsetConsidered == 2 || IsOnlyAnOrientationSubsetConsidered == 3 || IsOnlyAnOrientationSubsetConsidered == 4 || IsOnlyAnOrientationSubsetConsidered == 5 || IsOnlyAnOrientationSubsetConsidered == 6)
        {
            if (fscanf(fp, "%f %f %f %d %f %f %f %f %f %d %d %d %d \n",
                       &x, &y, &z, &phase, &rho, &ea[0], &ea[1], &ea[2], &CI, &RX, &IsCellGammaFiber, &IsCellThetaFiber, &IsCellComponent) != 13)
            {
                LOG_F(ERROR, "Error reading microstructure file %s at line %d", startFileName.c_str(), line + 1);
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            if (fscanf(fp, "%f %f %f %d %f %f %f %f %f %d \n",
                       &x, &y, &z, &phase, &rho, &ea[0], &ea[1], &ea[2], &CI, &RX) != 10)
            {
                LOG_F(ERROR, "Error reading microstructure file %s at line %d", startFileName.c_str(), line + 1);
                exit(EXIT_FAILURE);
            }
        }

        if (abs(round(x / deltax) - i) > 1 || abs(round(y / deltay) - j) > 1 || (!is2d && abs(round(z / deltax) - k) > 1))
        {
            LOG_F(ERROR, "Wrong coordinates in line %d. read: %d %d ; should be: %d %d", line + 2, int(round(x / deltax)), int(round(y / deltay)), i, j);
            exit(EXIT_FAILURE);
        }

        mvpStateVars->SetLatticeIdOfCell(index, phase);
        // LOG_F(INFO, "read %d and %d", phase, mvpStateVars->GetLatticeIdOfCell(index));

        id = index;
        // Add orientation and get its id
        OriId = mvOrientations.addOriAndReturnId(ea);

        // Set orientation id of the cell
        mvpStateVars->SetOriIdOfCell(index, OriId);

        // Debug print
        if (index == 0)
            LOG_F(INFO, "Line: %d, Index: %d, Orientation ID: %d", line, index, OriId);

        // Debug print
        if (index == 0)
        {
            Eigen::Vector3f eaAssigned = mvOrientations.GetEulerAngles(OriId);
            LOG_F(INFO, "ea phi1: %f Phi: %f phi2 %f", eaAssigned[0], eaAssigned[1], eaAssigned[2]);
        }

        if (!CheckedState && RX != 0)
        {
            LOG_F(WARNING, "Keep in mind you read recrystallized cell - if you were not intending to then stop simulation and change last column in microstructure file in 0 values");
            CheckedState = true;
        }
        mvpStateVars->SetEBSDRelatedStuffForRexAndGG(index, CI, RX, rho * DDFactor);

        if (++i == mvpStateVars->mvNx)
        {
            if (++j == mvpStateVars->mvNy)
            {
                k++;
                j = 0;
            }
            i = 0;
        };

        index++;
    }

    LOG_F(INFO, "Finished reading microstructure");

    // Debug print orientations
    // mvOrientations.DebugPrintOrientations();
}

void Microstructure::WriteMicrostructureRX(const std::string &fullPath)
{

    std::ofstream outputFile(fullPath);
    if (!outputFile.is_open())
    {
        LOG_F(ERROR, "Error opening output file: %s", fullPath.c_str());
        LOG_F(ERROR, "Current path: %s", std::filesystem::current_path().c_str());

        return;
    }

    // Write data to the file
    int N = mvpStateVars->mvN;
    Eigen::Vector3f eaOfCell;

    outputFile << " x " << " y " << " z " << " phi1 " << " Phi " << " phi2 " << " oriId " << " IsRex " << " MaxAnglePassed " << std::endl;

    for (int index = 0; index < N; index++)
    {
        int oriId = mvpStateVars->GetOriIdOfCell(index);
        eaOfCell = mvOrientations.GetEulerAngles(oriId);
        float x = mvpStateVars->GetXFromIndex(index) * 1.e6;
        float y = mvpStateVars->GetYFromIndex(index) * 1.e6;
        float z = mvpStateVars->GetZFromIndex(index) * 1.e6;
        outputFile << x << " " << y << " " << z << " "
                   << eaOfCell[0] << " " << eaOfCell[1] << " " << eaOfCell[2] << " "
                   << oriId << " " << mvpStateVars->mvpRX[index] << " " << mvpStateVars->mvpMaxGbPassed[index] << std::endl;
    }

    LOG_F(INFO, "Successfully wrote data to: %s", fullPath.c_str());

    outputFile.close();
}

// Helper function to generate color from value
void Microstructure::GetColorFromValue(float value, float minValue, float maxValue, uint8_t &r, uint8_t &g, uint8_t &b)
{
    float ratio = 2 * (value - minValue) / (maxValue - minValue);
    r = std::max(0.0f, 255 * (ratio - 1));
    g = std::max(0.0f, 255 * (1 - std::abs(ratio - 1)));
    b = std::max(0.0f, 255 * (1 - ratio));
}

// // Test function to validate the IPF color coding
// void Microstructure::TestIPFColorCoding() {
//     Eigen::Vector3f ea;
//     uint8_t r, g, b;

//     // Test case 1: should be reddish
//     ea = Eigen::Vector3f(5.0f, 5.0f, 45.0f); // ea[1] is arbitrary
//         // Define the principal directions
//     Eigen::Vector3f xDir(1, 0, 0);
//     Eigen::Vector3f yDir(0, 1, 0);
//     Eigen::Vector3f zDir(0, 0, 1);

//     // Compute IPF colors for each direction
//     computeIPFColor(eulerAngles, xDir, rX, gX, bX);
//     computeIPFColor(eulerAngles, yDir, rY, gY, bY);
//     computeIPFColor(eulerAngles, zDir, rZ, gZ, bZ);
//     mvOrientations.GetIPFColorFromOrientation(ea, r, g, b);
//     std::cout << "Test case 1: R=" << int(r) << " G=" << int(g) << " B=" << int(b) << std::endl;

//     // Test case 2: ea[1]=54.7 degrees and ea[2]=45 degrees should be blueish
//     ea = Eigen::Vector3f(30.0f, 54.7f, 45.0f); // ea[0] is arbitrary
//     mvOrientations.GetIPFColorFromOrientation(ea, r, g, b);
//     std::cout << "Test case 2: R=" << int(r) << " G=" << int(g) << " B=" << int(b) << std::endl;

//     // Test case 3: ea[0]=85 degrees and ea[2]=45 degrees should be greenish
//     ea = Eigen::Vector3f(85.0f, 30.0f, 45.0f); // ea[1] is arbitrary
//     mvOrientations.GetIPFColorFromOrientation(ea, r, g, b);
//     std::cout << "Test case 3: R=" << int(r) << " G=" << int(g) << " B=" << int(b) << std::endl;
// }

void Microstructure::GetColorCoding(int width, std::vector<uint8_t> &raster, int MapType)
{
    int i, j, k;
    float minCarbon = std::numeric_limits<float>::max();
    float maxCarbon = std::numeric_limits<float>::lowest();

    if (MapType == 1 || MapType == 2)
    {
        for (int index = 0; index < mvpStateVars->mvN; ++index)
        {
            float carbon = (MapType == 1) ? (*mvpStateVars->mvpXC)[index] : (*mvpStateVars->mvpXCTrapped)[index];
            minCarbon = std::min(minCarbon, carbon);
            maxCarbon = std::max(maxCarbon, carbon);
        }
    }

    for (int index = 0; index < mvpStateVars->mvN; ++index)
    {
        mvpStateVars->IndexToIJK(index, &i, &j, &k);
        int pixelIndex = j * width + i;

        uint8_t r = 0, g = 0, b = 0;
        if (MapType == 1 || MapType == 2)
        {
            float carbon = (MapType == 1) ? (*mvpStateVars->mvpXC)[index] : (*mvpStateVars->mvpXCTrapped)[index];
            GetColorFromValue(carbon, minCarbon, maxCarbon, r, g, b);
        }
        else
        {
            if (mvpStateVars->IsCellFCC(index))
            {
                r = 255;
                g = 153;
                b = 255; // pink
            }
            else if (mvpStateVars->IsCellBCC(index))
            {
                r = 51;
                g = 204;
                b = 255; // light blue
            }
            else if (mvpStateVars->IsCellHCP(index))
            {
                r = 0;
                g = 204;
                b = 0; // green
            }
            if (mvpStateVars->IsCellInterface(index))
            {
                float DarkFactor = GetCellsMisorientationDarkeningFactor(index);
                r -= (r * DarkFactor);
                g -= (g * DarkFactor);
                b -= (b * DarkFactor);
            }
        }

        raster[pixelIndex * 3] = r;
        raster[pixelIndex * 3 + 1] = g;
        raster[pixelIndex * 3 + 2] = b;
    }
}

void Microstructure::ExportColorCodedTIFF(const std::string &fullPath, int MapType)
{
    int width = mvUserSettings.mvNx;
    int height = mvUserSettings.mvNy;
    std::vector<uint8_t> raster(width * height * 3);

    TIFF *out = TIFFOpen(fullPath.c_str(), "w");
    if (!out)
    {
        LOG_F(ERROR, "Could not open %s for writing", fullPath.c_str());
        return;
    }

    GetColorCoding(width, raster, MapType);

    float FactorToScale = (mvpStateVars->mvGridType == 1) ? 0.866f : 1.0f;
    int scaledHeight = static_cast<int>(height * FactorToScale);

    std::vector<uint8_t> scaledRaster(width * scaledHeight * 3, 0);
    if (FactorToScale >= 0.99f)
    {
        scaledRaster = raster;
    }
    else
    {
        for (int y = 0; y < scaledHeight; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                int originalY = static_cast<int>(y / FactorToScale);
                if (originalY < height)
                {
                    int originalIndex = (originalY * width + x) * 3;
                    int scaledIndex = (y * width + x) * 3;

                    scaledRaster[scaledIndex] = raster[originalIndex];
                    scaledRaster[scaledIndex + 1] = raster[originalIndex + 1];
                    scaledRaster[scaledIndex + 2] = raster[originalIndex + 2];
                }
            }
        }
    }

    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, scaledHeight);
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

    TIFFWriteEncodedStrip(out, 0, scaledRaster.data(), scaledRaster.size());
    TIFFClose(out);
    scaledRaster.clear();
    raster.clear();
}

void Microstructure::ExportColorbar(const std::string &fullPath)
{
    int width = 256; // Width of the colorbar
    int height = 20; // Height of the colorbar
    std::vector<uint8_t> colorbar(width * height * 3);

    float minCarbon = std::numeric_limits<float>::max();
    float maxCarbon = std::numeric_limits<float>::lowest();

    for (int index = 0; index < mvpStateVars->mvN; ++index)
    {
        float carbon = (*mvpStateVars->mvpXC)[index];
        minCarbon = std::min(minCarbon, carbon);
        maxCarbon = std::max(maxCarbon, carbon);
    }

    for (int x = 0; x < width; ++x)
    {
        float value = minCarbon + (maxCarbon - minCarbon) * x / (width - 1);
        uint8_t r, g, b;
        GetColorFromValue(value, minCarbon, maxCarbon, r, g, b);

        for (int y = 0; y < height; ++y)
        {
            int pixelIndex = (y * width + x) * 3;
            colorbar[pixelIndex] = r;
            colorbar[pixelIndex + 1] = g;
            colorbar[pixelIndex + 2] = b;
        }
    }

    TIFF *out = TIFFOpen(fullPath.c_str(), "w");
    if (!out)
    {
        LOG_F(ERROR, "Could not open %s for writing", fullPath.c_str());
        return;
    }

    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

    TIFFWriteEncodedStrip(out, 0, colorbar.data(), colorbar.size());
    TIFFClose(out);
}

bool Microstructure::IsTheNumberFinite(double numberToCheck)
{
    if (numberToCheck > std::numeric_limits<double>::max() || numberToCheck < std::numeric_limits<double>::lowest())
        return false;
    return true;
}

bool Microstructure::IsTheNumberFinite(float numberToCheck)
{
    if (numberToCheck > std::numeric_limits<float>::max() || numberToCheck < std::numeric_limits<float>::lowest())
        return false;
    return true;
}

double Microstructure::MisorientationBetweenCells(int index1, int index2, bool MaxHagb, bool CheckAlsoForSameOrientation, bool CheckIfCSL, int *CSLRelationship)
{
    double mis;
    *CSLRelationship = 0;
    Eigen::Vector3f ea1, ea2;

    int OriId1 = mvpStateVars->GetOriIdOfCell(index1);
    int OriId2 = mvpStateVars->GetOriIdOfCell(index2);

    if (!CheckAlsoForSameOrientation && OriId1 == OriId2)
        return 0.;

    mis = mvOrientations.CalculateMisorientationBetweenTwoOriIds(OriId1, OriId2);
    if (!IsTheNumberFinite(mis))
    { // isnan is undefined now in latest cmath (from math) and cannot be reloaded. But we need to check because if the ea is identical ro the ea1 the misorientation is nan
        // std::cout<<" ea[0]  "<<ea[0] <<" ea[1] "<<ea[1]<<" ea[2] "<<ea[2]<<" ea1[0]  "<<ea1[0] <<" ea1[1] "<<ea1[1]<<" ea1[2] "<<ea1[2]<<std::endl;
        mis = 0.;
    }
    if (mis < 0. || mis > 64.)
    {
        ea1 = mvOrientations.GetEulerAngles(OriId1);
        ea2 = mvOrientations.GetEulerAngles(OriId2);

        std::cerr << " Wrong misorientation value in MisorientationBetweenCells: phi1 " << ea1[0] << " Phi " << ea1[1] << " phi2 " << ea1[2] << " phi1 " << ea2[0] << " Phi " << ea2[1] << " phi2 " << ea2[2] << " gave mis " << mis << std::endl;
    }
    if (CheckIfCSL)
    {
        ea1 = mvOrientations.GetEulerAngles(OriId1);
        ea2 = mvOrientations.GetEulerAngles(OriId2);
        float deltaMisorientation = mvOrientations.DistanceOfBoundaryFromCSL19a(ea1, ea2);
        if (deltaMisorientation < 15.)
            *CSLRelationship = 19;
        // LOG_F(INFO,"CSL19a between orientations %e, %e,%e, and %e, %e,%e, as the distance is  %e ", ea[0] ,ea[1] ,ea[2],ea1[0] ,ea1[1] ,ea1[2], deltaMisorientation);
    }

    if (mis > mvUserSettings.mvHAGB && MaxHagb)
        return mvUserSettings.mvHAGB; // if we want the misorientation only to see if it is hagb (e.g. mobility calculations) we return max{mvUserSettings.mvHAGB,mis}, Otherwise (e.g. for plotting purposes) we return actual value (MaxHagb==false in that case).
    return mis;
}

void Microstructure::ValidateOrientationSettingsAndCheckNonIndexedPixels(void)
{
    int CSLRelationship;
    float maxMis = 0.;
    int nbrListTemp[26];
    int NeighNonInd = 0;
    int NonInd = 0;

    for (int i = 0; i < mvpStateVars->mvN; i++)
    {
        if (mvUserSettings.mvIsTestRVE)
            mvpStateVars->mvpCI[i] = 1.;
        mvpStateVars->SetIfCellHasNonIndNeighbours(i, 0);
        // misorientation from itself - to check if error is smaller than MisTol
        double mis = MisorientationBetweenCells(i, i, true, true, false, &CSLRelationship);
        if (mis > maxMis)
            maxMis = mis;
        if (maxMis > mvUserSettings.mvLowerMisorientationCutOff)
        {
            std::cerr << "change LowerMisorientationCutOff in Settings file, you have super low cutoff angle, i.e. " << mvUserSettings.mvLowerMisorientationCutOff << " and max misorientation (error) between same orientation  " << maxMis << std::endl;
        }
    }
    if (mvUserSettings.mvHasSoluteDiffusion)
    {
        for (int i = 0; i < mvpStateVars->mvN; i++)
        {
            mvpStateVars->GetNeighbourCellOffsets(i, nbrListTemp);
            int SameLatticeTypeNeighbours = 0;
            for (int nbr = 0; nbr < mvpStateVars->mvNAll; nbr++)
            {
                if (mvpStateVars->mvpLatticeId[nbrListTemp[nbr]] == mvpStateVars->mvpLatticeId[i])
                {
                    SameLatticeTypeNeighbours++;
                    break;
                }
            }
            if (SameLatticeTypeNeighbours < 1)
                mvpStateVars->mvpCI[i] = 0.5 * mvUserSettings.mvMinimumCI;
        }
    }

    for (int i = 0; i < mvpStateVars->mvN; i++)
    {
        if (mvpStateVars->mvpCI[i] < mvUserSettings.mvMinimumCI)
            NonInd++;
        mvpStateVars->GetNeighbourCellOffsets(i, nbrListTemp);
        for (int nbr = 0; nbr < mvpStateVars->mvNAll; nbr++)
        {
            if (mvpStateVars->mvpCI[nbrListTemp[nbr]] < mvUserSettings.mvMinimumCI)
            {
                mvpStateVars->SetIfCellHasNonIndNeighbours(i, 1);
                NeighNonInd++;
                break;
            }
        }
    }

    LOG_F(INFO, "Finished checking MaxMisorientationError, equal to %g and you have set tol %g", maxMis, mvUserSettings.mvLowerMisorientationCutOff);
    LOG_F(INFO, "Finished checking orientation map, you have non Indexed %d and cells with non indexed neighbours non Indexed %d", NonInd, NeighNonInd);
}

void Microstructure::InitInterfaceAndInterphaseCells()
{
    int nbrList1Temp[26];
    int nBoundaries = 0;
    int nInterphases = 0;
    int indexNbr;
    int CSLRelationship;
    float mis;
    for (int index = 0; index < mvpStateVars->mvN; index++)
    {
        mvpStateVars->GetNeighbourCellOffsets(index, nbrList1Temp);
        for (int nbr = 0; nbr < mvpStateVars->mvNAll; nbr++)
        {
            indexNbr = nbrList1Temp[nbr];
            if (mvpStateVars->GetLatticeIdOfCell(index) != mvpStateVars->GetLatticeIdOfCell(indexNbr))
            {
                nInterphases++;
                mvpStateVars->SetCellAsInterphase(index);
                AddInterphaseIndex(index);
                break;
            }
        }
    }
    for (int index = 0; index < mvpStateVars->mvN; index++)
    {
        mvpStateVars->GetNeighbourCellOffsets(index, nbrList1Temp);
        for (int nbr = 0; nbr < mvpStateVars->mvNAll; nbr++)
        {
            indexNbr = nbrList1Temp[nbr];
            if (mvpStateVars->GetLatticeIdOfCell(index) == mvpStateVars->GetLatticeIdOfCell(indexNbr))
            {
                mis = MisorientationBetweenCells(index, indexNbr, true, false, false, &CSLRelationship);
                if (mis >= mvUserSettings.mvLowerMisorientationCutOff)
                {
                    nBoundaries++;
                    mvpStateVars->SetCellAsInterface(index);
                    AddInterfaceIndex(index);
                    break;
                }
            }
        }
    }
    LOG_F(INFO, "Initialized (sub)grain and phase boundary cells: we have %d phase boundary cells and %d (sub)grain boundary cells out of %d cells",
          nInterphases, nBoundaries, mvpStateVars->mvN);
}

void Microstructure::CheckInterfaceAndInterphaseCells()
{
    LOG_F(INFO, " Checking boundary and interphase cells ");
    int nbrList1Temp[26];
    int nBoundaries = 0;
    int nInterphases = 0;
    int nBoundariesRemoved = 0;
    int nInterphasesRemoved = 0;
    int indexNbr;
    int CSLRelationship;
    float mis;
    for (int index = 0; index < mvpStateVars->mvN; index++)
    {
        bool IsIPhase = IsCellIndeedInterphase(index);
        if (IsIPhase)
            nInterphases++;
        if (IsIPhase && !mvpStateVars->IsCellInterphase(index))
            ABORT_F(" A cell is found with dissimilar (phase) neighbour - YET THIS CELL IS NOT STORED AS PHASE BOUNDARY CELL. The program will exit and check what happened ");
        if (IsIPhase && mvpStateVars->IsCellInterphase(index))
            continue;
        if (!IsIPhase && mvpStateVars->IsCellInterphase(index))
        {
            RemoveInterphaseIndex(index);
            nInterphasesRemoved++;
        }
        bool IsBoundary = IsCellIndeedInterface(index);
        if (IsBoundary)
            nBoundaries++;

        if (IsBoundary && !mvpStateVars->IsCellInterface(index))
            ABORT_F(" A cell is found with dissimilar (orientation) neighbour - YET THIS CELL IS NOT STORED AS A (SUB)GRAIN BOUNDARY CELL. The program will exit and check what happened ");
        if (!IsBoundary && mvpStateVars->IsCellInterface(index))
        {
            RemoveInterfaceIndex(index);
            nBoundariesRemoved++;
        }
    }
    LOG_F(INFO, " All fine, we checked AND PASSED the check for (sub)grain and phase boundary cells: we have %d interphases and %d boundaries out of %d cells", nInterphases, nBoundaries, mvpStateVars->mvN);
    LOG_IF_F(WARNING, nBoundariesRemoved > 0 || nInterphasesRemoved > 0, "It is not wrong for the simulation, but maybe check your Switch/Updating cells functions. Removed %d interphases and %d boundaries", nInterphasesRemoved, nBoundariesRemoved);
    LOG_F(INFO, " We check again after %d simulation steps", mvEveryThatManySimStepsToCheckICells);
}

float Microstructure::SetStateVariablesBasedOnSimAndGiveInfoForMicrostructuralState()
{
    NcellsForCarbonPart = mvpStateVars->mvNNearest;
    NcellsForIntMig = mvpStateVars->mvNAll;
    LOG_F(INFO, "Setting state variables of full-field simulations, i.e. GG, Rex, Or some simulation types for carbon long range Diffusion");
    BetaFactorForInterfaceMigrationRates = 1.;

    double tolerance = 1e-6;
    std::vector<double> validParameters = {0.5, 0.7, 1.0};

    auto isValidParameter = [&](double param)
    {
        return std::any_of(validParameters.begin(), validParameters.end(), [&](double v)
                           { return std::abs(param - v) < tolerance; });
    };

    if (!isValidParameter(mvUserSettings.mvRexGGParameterC))
    {
        ABORT_F("Cannot simulate full-field Rex or Grain Growth with the value of %g you gave for parameter C. "
                "Go back to input settings and put C either 0.5, 0.7 or 1., but preferably 0.7",
                mvUserSettings.mvRexGGParameterC);
    }
    else
    {
        if (mvUserSettings.mvGridType == 1 && std::abs(mvUserSettings.mvRexGGParameterC - 0.5) < tolerance)
            BetaFactorForInterfaceMigrationRates = 0.66;
        if (mvUserSettings.mvGridType == 1 && std::abs(mvUserSettings.mvRexGGParameterC - 0.7) < tolerance)
            BetaFactorForInterfaceMigrationRates = 0.47;
        if (mvUserSettings.mvGridType == 1 && std::abs(mvUserSettings.mvRexGGParameterC - 1.0) < tolerance)
            BetaFactorForInterfaceMigrationRates = 0.37;
        if (mvUserSettings.mvGridType == 0 && mvpStateVars->mvIs2DSimulation && std::abs(mvUserSettings.mvRexGGParameterC - 0.5) < tolerance)
            BetaFactorForInterfaceMigrationRates = 0.58;
        if (mvUserSettings.mvGridType == 0 && mvpStateVars->mvIs2DSimulation && std::abs(mvUserSettings.mvRexGGParameterC - 0.7) < tolerance)
            BetaFactorForInterfaceMigrationRates = 0.49;
        if (mvUserSettings.mvGridType == 0 && mvpStateVars->mvIs2DSimulation && std::abs(mvUserSettings.mvRexGGParameterC - 1.0) < tolerance)
            BetaFactorForInterfaceMigrationRates = 0.33;
        if (mvUserSettings.mvGridType == 0 && !mvpStateVars->mvIs2DSimulation && std::abs(mvUserSettings.mvRexGGParameterC - 0.5) < tolerance)
            BetaFactorForInterfaceMigrationRates = 0.37;
        if (mvUserSettings.mvGridType == 0 && !mvpStateVars->mvIs2DSimulation && std::abs(mvUserSettings.mvRexGGParameterC - 0.7) < tolerance)
            BetaFactorForInterfaceMigrationRates = 0.32;
        if (mvUserSettings.mvGridType == 0 && !mvpStateVars->mvIs2DSimulation && std::abs(mvUserSettings.mvRexGGParameterC - 1.0) < tolerance)
            BetaFactorForInterfaceMigrationRates = 0.22;
        if (BetaFactorForInterfaceMigrationRates > 0.99)
            LOG_F(INFO, "Check grid and parameter C, because somehow your kinetic factor is equal to %g", BetaFactorForInterfaceMigrationRates);
        else
            LOG_F(INFO, "For grid type %d and parameter C %g your boundary migration rate now matches the theory, after scaling by a factor of %g your mobility/boundary energy input",
                  mvUserSettings.mvGridType,
                  mvUserSettings.mvRexGGParameterC,
                  BetaFactorForInterfaceMigrationRates);
    }

    int nbrList1Temp[26];
    int CSLRelationship;
    float maxmis = 0.;
    double mis, BoundaryEnergyInGrid, IncSigma;
    int BCells = 0;
    cells_rx = 0;
    frx = 0.;
    Eigen::Vector3f ea1, ea2;

    //    Eigen::Vector3f CSL19a,CSL21a,CSL27a;
    Eigen::Quaternionf q;

    if (mvUserSettings.mvHasSoluteDiffusion && !mvpStateVars->mvpDiffusivityC_Ratio)
    {
        LOG_F(WARNING, " strange - you have full diffusion yet the mvpDiffusivityC_Ratio is not initialized - i initialize it now but TAKE CARE (IT IS PART OF Ax=B system in CGsolver)");
        mvpStateVars->mvpDiffusivityC_Ratio = new double[mvpStateVars->mvN * mvpStateVars->mvNNearest]; // although it only matters for having different diffusivities (e.g. Agren, or solving martensite-austenite), the way the Ax=b code is written now  needs the Diff matrix (i.e. it is part of A matrix, and takes values of 1 if all diffusivities are equal (i.e. equal to the maxdiffusivity))
        std::memset(mvpStateVars->mvpDiffusivityC_Ratio, 1., mvpStateVars->mvN * mvpStateVars->mvNNearest);
    }
    if (mvUserSettings.mvHasSoluteDiffusion && !mvpStateVars->mvpKappaFactorForCTrappedDefects)
    {
        LOG_F(WARNING, " strange - you have full diffusion yet the kappaFactorForCdefects is not initialized - i initialize it now but TAKE CARE (IT IS PART OF Ax=B system in CGsolver)");
        mvpStateVars->mvpKappaFactorForCTrappedDefects = new double[mvpStateVars->mvN]; // although it only matters for Csegregation to defects, the way the code is writtem now whole Cdiffusion step is solved
        std::memset(mvpStateVars->mvpKappaFactorForCTrappedDefects, 0., mvpStateVars->mvN);
    }

    for (int i = 0; i < mvpStateVars->mvN; i++)
    {
        if (!mvUserSettings.mvReadDislocationDensity)
            mvpStateVars->mvpRho[i] = 1.e12;

        if (mvUserSettings.mvIsTestRVE)
        {
            mvpStateVars->mvpHasNeighNonInd[i] = 0;
            mvpStateVars->mvpCI[i] = 1.;
            mvpStateVars->mvpRX[i] = 1.;
            mvpStateVars->mvpKAM[i] = 1.;
        }
        mvpStateVars->mvpKAM[i] = 0.;

        if (mvUserSettings.mvIsRexAndGG || mvUserSettings.mvAllowInterfaceMovementDuringPartitioning)
        {

            if (mvpStateVars->mvpBoundaryCell[i] && !mvpStateVars->mvpRX[i])
            {
                BCells++;
                mvpStateVars->mvpCI[i] = 0.5 * mvUserSettings.mvMinimumCI; // CI affects whether rex is allowed to initiate there,
            }
            if (mvUserSettings.mvIncludeCSL19FastGrowth)
            {
                if (mvpStateVars->mvpCSLRelationshipNow == nullptr)
                    mvpStateVars->mvpCSLRelationshipNow = new int[mvpStateVars->mvN];
                mvpStateVars->mvpCSLRelationshipNow[i] = 0;
            }
        }
    }

    LOG_F(INFO, "bcells %d", BCells);
    BoundaryEnergyInGrid = 0.;

    float WeightedBoundaryArea;

    double meanKAM_FCC = 0.;
    double meanKAM_BCC = 0.;
    int NFCC = 0;
    int NBCC = 0;
    int numberInterphases = 0;

    // This loop will determine which cell's orientation (and neighborhood's misorientations) are
    // trustful enough (useful for EBSD inputs) to allow growth or rex to initiate there (also boundary cells cannot - there is no HAGB in reality)
    // Also, some quantities for determining later start here
    for (int i = 0; i < mvpStateVars->mvN; i++)
    {
        mvpStateVars->mvpKAM[i] = 0.;
        if (mvpStateVars->IsCellBCC(i))
            NBCC++;
        if (mvpStateVars->IsCellFCC(i))
            NFCC++;
        if (!mvUserSettings.mvAllowPeriodicBoundConditions && mvpStateVars->mvpBoundaryCell[i])
            continue;

        mvpStateVars->GetNeighbourCellOffsets(i, nbrList1Temp);
        float NNeighs = 0.;
        bool CellIsHAGB = false;
        for (int nbr = 0; nbr < mvpStateVars->mvNAll; nbr++)
        {

            // if (mvpStateVars->mvpLatticeId[nbrList1Temp[nbr]] != mvpStateVars->mvpLatticeId[i]) // Interphase mis doesnt contribute to average mis
            //     continue;
            int OriId1 = mvpStateVars->GetOriIdOfCell(i);
            ea1 = mvOrientations.GetEulerAngles(OriId1);
            int OriId2 = mvpStateVars->GetOriIdOfCell(nbrList1Temp[nbr]);
            ea2 = mvOrientations.GetEulerAngles(OriId2);

            if (ea1[0] == ea2[0] && ea1[1] == ea2[1] && ea1[2] == ea2[2])
                continue;
            if (mvpStateVars->mvpCI[nbrList1Temp[nbr]] < mvUserSettings.mvMinimumCI || mvpStateVars->mvpCI[i] < mvUserSettings.mvMinimumCI)
                continue;

            if (mvpStateVars->mvpLatticeId[nbrList1Temp[nbr]] != mvpStateVars->mvpLatticeId[i])
            {
                // Interphase mis doesnt contribute to average mis
                numberInterphases++;
                continue;
            }

            mis = MisorientationBetweenCells(i, nbrList1Temp[nbr], true, false, mvUserSettings.mvIncludeCSL19FastGrowth, &CSLRelationship);
            if (CSLRelationship > 0)
            {
                mvpStateVars->mvpCSLRelationshipNow[i] = 19;
                mvpStateVars->mvpCSLRelationshipNow[nbrList1Temp[nbr]] = 19;
            }

            if (mis >= mvUserSettings.mvLowerMisorientationCutOff)
            {
                mvpStateVars->mvpKAM[i] += mis;
                if (mis >= mvUserSettings.mvHAGB)
                    CellIsHAGB = true;
            }

            NNeighs += 1.;

            IncSigma = 0.;
            if (mis >= mvUserSettings.mvLowerMisorientationCutOff && mvpStateVars->mvpCI[nbrList1Temp[nbr]] >= mvUserSettings.mvMinimumCI && mvpStateVars->mvpCI[i] >= mvUserSettings.mvMinimumCI)
                IncSigma = GetBoundaryEnergyForThisMisorientation(mis);
            if (mvUserSettings.mvGridType == 0)
            {
                WeightedBoundaryArea = mvpStateVars->mvDx / mvpStateVars->GetDistanceBetweenCells(i, nbrList1Temp[nbr]);
                IncSigma *= WeightedBoundaryArea;
            }
            BoundaryEnergyInGrid += (IncSigma);
        }

        if (NNeighs > 0.)
        {
            mvpStateVars->mvpKAM[i] /= NNeighs;
            if (CellIsHAGB)
                mvpStateVars->mvpKAM[i] = mvUserSettings.mvHAGB;
        }

        if (NNeighs < 0.01)
            mvpStateVars->mvpKAM[i] = 0.;
        if (mvpStateVars->IsCellFCC(i))
        {
            meanKAM_FCC += mvpStateVars->mvpKAM[i];
        }
        if (mvpStateVars->IsCellBCC(i))
        {
            meanKAM_BCC += mvpStateVars->mvpKAM[i];
        }
    }

    if (NBCC > 0)
        LOG_IF_F(INFO, NBCC > 0, "Calculated BCC: %d cells in, we have ave kam %e cumulative kam %e ", NBCC, meanKAM_BCC / NBCC, meanKAM_BCC);

    if (NFCC > 0)
        LOG_IF_F(INFO, NFCC > 0, "Calculated FCC: %d cells in, we have ave kam %e cumulative kam %e ", NFCC, meanKAM_FCC / NFCC, meanKAM_FCC);

    meanKAM_BCC /= NBCC;
    meanKAM_FCC /= NFCC;
    // std::cout<<" all fine 11 "<<std::endl;

    double DDMaxBCCScaled = 0.;
    if (mvUserSettings.mvHasSoluteDiffusion)
    {
        for (int i = 0; i < mvpStateVars->mvN; i++)
        {
            (*mvpStateVars->mvpXC)[i] = mvUserSettings.mvAverageCarbon;
            // std::cout<<" all fine 111 "<<std::endl;

            if (mvUserSettings.mvSoluteSegregationDislocations)
            {
                (*mvpStateVars->mvpXCTrapped)[i] = 0.;
                mvpStateVars->mvpKappaFactorForCTrappedDefects[i] = 0.;
            }
            // std::cout<<" all fine 1111 "<<std::endl;

            if (mvUserSettings.mvAllowInterfaceMovementDuringPartitioning)
            {
                grewFCC = 0.;
                grewBCC = 0.;
            }
            // std::cout<<" all fine 11111 "<<std::endl;

            // std::cout<<" all fine 1 "<<std::endl;
            // This is just a way to make articial strips of high defect density in microstructure to test Ctrapping to defects
            if (mvUserSettings.mvIsTestRVE && mvUserSettings.mvSoluteSegregationDislocations)
            {
                int distKAM2 = mvpStateVars->mvNx / 5.;
                int distKAM10 = mvpStateVars->mvNx / 4.;
                int widthKAMZones = mvpStateVars->mvNx / 20.;
                mvpStateVars->mvpKAM[i] = 0.;
                mvpStateVars->mvpRho[i] = 0.;
                int x, y, z;
                // prevent double counting of same corner
                mvpStateVars->IndexToIJK(i, &x, &y, &z);
                if ((x % distKAM2 < widthKAMZones / 2.) || (distKAM2 - (x % distKAM2) < widthKAMZones / 2.))
                {
                    mvpStateVars->mvpKAM[i] = 1.;
                    mvpStateVars->mvpRho[i] = 1.e15;
                }
                if ((x % distKAM10 < widthKAMZones / 2.) || (distKAM10 - (x % distKAM10) < widthKAMZones / 2.))
                {
                    mvpStateVars->mvpKAM[i] = 5.;
                    mvpStateVars->mvpRho[i] = 5.e15;
                }
            }
            // std::cout<<" all fine 2 "<<std::endl;

            // if orientations are per cell then kam is going to be recalculated later and so dislospacing and number. Otherwise we keep the KAM as read.
            if (mvUserSettings.mvSoluteSegregationDislocations)
            {
                double rho = 1.e12;
                if (!mvUserSettings.mvReadDislocationDensity && !mvUserSettings.mvIsTestRVE)
                {
                    float b = (mvpStateVars->IsCellFCC(i) ? mvUserSettings.mvBurgersVectorFCC : mvUserSettings.mvBurgersVectorBCC);
                    b *= 1.e6; // IN MICROMETER
                    float disloTypeCoeff = 3.;
                    double rho = disloTypeCoeff * mvpStateVars->mvpKAM[i] * 0.0175 / b / mvpStateVars->mvDx / 1.e6; // IN 1/SQUARE MICROMETER
                    mvpStateVars->mvpRho[i] = rho * 1.e12;                                                          // IN 1/SQUARE METER
                }
                mvpStateVars->mvpKappaFactorForCTrappedDefects[i] = mvUserSettings.mvKappaFactorForCsegInClusterArea * mvpStateVars->mvpRho[i];
                if (!mvUserSettings.mvAllowDislocationsSegregationInAustenite && mvpStateVars->IsCellFCC(i))
                    mvpStateVars->mvpKappaFactorForCTrappedDefects[i] = 0.;
            }
        }
        // std::cout<<" all fine 3 "<<std::endl;

        if (mvUserSettings.mvSoluteSegregationDislocations)
        {
            mvpStateVars->SetCarbonTrappedAndCarbonFreeForGivenXcTot();
            // std::cout<<" all fine 4 "<<std::endl;

            meanKAM_FCC = 0.;
            meanKAM_BCC = 0.;
            double meanRhoBCC = 0.;
            double meanRhoFCC = 0.;
            double meanXCTrappedFCC = 0.;
            double meanXCTrappedBCC = 0.;
            double meanKappaFCC = 0.;
            double meanKappaBCC = 0.;
            for (int i = 0; i < mvpStateVars->mvN; i++)
            {
                if (mvpStateVars->IsCellFCC(i))
                {
                    if (mvpStateVars->mvpKAM[i] > 0.)
                        meanKAM_FCC += mvpStateVars->mvpKAM[i];
                    if (mvpStateVars->mvpRho[i] > 0.)
                        meanRhoFCC += (mvpStateVars->mvpRho[i] * 1.e-12);
                    if (mvpStateVars->mvpKappaFactorForCTrappedDefects[i] > 0.)
                    {
                        meanXCTrappedFCC += (*mvpStateVars->mvpXCTrapped)[i];
                        meanKappaFCC += mvpStateVars->mvpKappaFactorForCTrappedDefects[i];
                    }
                }
                if (mvpStateVars->IsCellBCC(i))
                {
                    if (mvpStateVars->mvpKAM[i] > 0.)
                        meanKAM_BCC += mvpStateVars->mvpKAM[i];
                    if (mvpStateVars->mvpRho[i] > 0.)
                        meanRhoBCC += (mvpStateVars->mvpRho[i] * 1.e-12);
                    if (mvpStateVars->mvpKappaFactorForCTrappedDefects[i] > 0.)
                    {
                        meanXCTrappedBCC += (*mvpStateVars->mvpXCTrapped)[i];
                        meanKappaBCC += mvpStateVars->mvpKappaFactorForCTrappedDefects[i];
                    }
                }
            }
            // std::cout<<" all fine 5 "<<std::endl;

            meanKAM_BCC /= NBCC;
            meanKAM_FCC /= NFCC;
            meanXCTrappedBCC /= NBCC;
            meanXCTrappedFCC /= NFCC;
            meanRhoBCC /= NBCC;
            meanRhoFCC /= NFCC;
            meanKappaBCC /= NBCC;
            meanKappaFCC /= NFCC;
            meanRhoBCC *= 1.e12;
            meanRhoFCC *= 1.e12;

            LOG_IF_F(INFO, NBCC > 0, "Just to check again, martensite, we have %d cells, and we have ave kam %e, cumulative kam %e, ave xcTrapped %e, cumulative xcTrapped %e : ", NBCC, meanKAM_BCC, meanKAM_BCC * NBCC, meanXCTrappedBCC, meanXCTrappedBCC * NBCC);
            LOG_IF_F(INFO, NFCC > 0, "Just to check again, FCC, we have %d cells, and we have ave kam %e, cumulative kam %e, ave xcTrapped %e, cumulative xcTrapped %e : ", NFCC, meanKAM_FCC, meanKAM_FCC * NFCC, meanXCTrappedFCC, meanXCTrappedFCC * NFCC);
            LOG_IF_F(INFO, (NBCC > 0 && mvUserSettings.mvReadDislocationDensity), "Just to check again, we read dislocation density from micro file, we have martensite, we have %d cells, and we have ave kam %e, cumulative kam %e, ave xcTrapped %e, cumulative xcTrapped %e, ave Rho %e, ave Kappa %e: ", NBCC, meanKAM_BCC, meanKAM_BCC * NBCC, meanXCTrappedBCC, meanXCTrappedBCC * NBCC, meanRhoBCC, meanKappaBCC);
            LOG_IF_F(INFO, (NBCC > 0 && !mvUserSettings.mvReadDislocationDensity), "Just to check again, we calculated here dislodensity from kam, we have martensite, we have %d cells, and we have ave kam %e, cumulative kam %e, ave xcTrapped %e, cumulative xcTrapped %e, ave Rho %e, ave Kappa %e: ", NBCC, meanKAM_BCC, meanKAM_BCC * NBCC, meanXCTrappedBCC, meanXCTrappedBCC * NBCC, meanRhoBCC, meanKappaBCC);
            LOG_IF_F(INFO, NFCC > 0, "Just to check again, austenite, we have %d cells, and we have ave kam %e, cumulative kam %e, ave xcTrapped %e, cumulative xcTrapped %e, ave Rho %e, ave Kappa %e: ", NFCC, meanKAM_FCC, meanKAM_FCC * NFCC, meanXCTrappedFCC, meanXCTrappedFCC * NFCC, meanRhoFCC, meanKappaFCC);
        }
    }

    BoundaryEnergyInGrid *= 0.5; // to correct double counting
    float BoundaryEnergyInGridPerCell = BoundaryEnergyInGrid / mvpStateVars->mvN;
    LOG_F(INFO, "Ncells %d A/v %g", mvpStateVars->mvN, mvpStateVars->mvMaxBoundaryAreaPerVol);
    LOG_F(INFO, "sigma %g TotalBoundEnergyInGrid %g and per cell %g maxmis %g", mvUserSettings.mvGrainBoundaryEnergy, BoundaryEnergyInGrid * mvpStateVars->mvMaxBoundaryAreaPerVol, BoundaryEnergyInGridPerCell * mvpStateVars->mvMaxBoundaryAreaPerVol, maxmis);
    LOG_F(INFO, "Finished settings properties for rex, gg, diffusion");
    mvpStateVars->eraseKAM();
    return BoundaryEnergyInGridPerCell;
}

double Microstructure::GetBoundaryEnergyForThisMisorientation(double mis)
{
    double sigma_hagb = mvUserSettings.mvGrainBoundaryEnergy;
    double r = mis / mvUserSettings.mvHAGB;
    if (r > 0.99)
        return sigma_hagb;
    if (r < mvUserSettings.mvLowerMisorientationCutOff / mvUserSettings.mvHAGB)
        return 0.; // <MisTol means no Boundary
    return sigma_hagb * r * (1. - log(r));
}

double Microstructure::GetMolarVolumeOfCell(int index)
{
    double Vm = std::numeric_limits<double>::quiet_NaN(); // Initialize Vm with NaN
    if (mvpStateVars->IsCellFCC(index))
        Vm = mvpTCK->AusteniteMolarVolume();
    if (mvpStateVars->IsCellBCC(index))
        Vm = mvpTCK->FerriteMolarVolume();
    if (std::isnan(Vm))
    { // Check if Vm is still NaN
        ABORT_F(" For now only iron BCC and iron FCC are handled for simulations - unknown phase so also molar volume ");
    }

    return Vm;
}

double Microstructure::GetMolarVolumeOfLattice(int latticeId)
{
    double Vm = std::numeric_limits<double>::quiet_NaN(); // Initialize Vm with NaN

    if (mvpStateVars->IsLatticeIdFCC(latticeId))
    {
        Vm = mvpTCK->AusteniteMolarVolume();
    }
    if (mvpStateVars->IsLatticeIdBCC(latticeId))
    {
        Vm = mvpTCK->FerriteMolarVolume();
    }

    if (std::isnan(Vm))
    { // Check if Vm is still NaN
        ABORT_F(" For now only iron BCC and iron FCC are handled for simulations - unknown phase so also molar volume ");
    }

    return Vm;
}
double Microstructure::GetMobilityForThisMisorientation(int latticeId, double mis)
{
    double HumphreysB = 5; ///< Parameter B in Humphreys' mobility equation (Acta Mater 45, 1997)  e.g. HumphreysB=5.
    double HumphreysN = 4; ///< Parameter N in Humphreys' mobility equation (Acta Mater 45, 1997)  e.g. HumphreysN=4.

    double MobilityOfMaxMisorientation = mvpTCK->GetGrainBoundaryMobility();
    // double MobilityOfMaxMisorientation = (mvUserSettings.UseSoluteDrag) ? mvpTCK->GetMobilityWithSD(phase, 1., 2, true) : mvpTCK->GetMobilityPure(phase); // [mol m / J sec]
    double Vm = GetMolarVolumeOfLattice(latticeId);
    double FactorForUnitsConversion = Vm;
    MobilityOfMaxMisorientation *= FactorForUnitsConversion; // [m4 / J sec]
    double r = mis / mvUserSettings.mvHAGB;
    double M;
    M = MobilityOfMaxMisorientation;
    if (r >= 0.99)
        return M;
    if (r < mvUserSettings.mvLowerMisorientationCutOff / mvUserSettings.mvHAGB)
        return 0.; // <MisTol means no Boundary
    M *= (1. - std::exp(-HumphreysB * pow(r, HumphreysN)));
    return M;
}

void Microstructure::CalculateDirectionOfForceAtThisBoundary(int index1, int index2, double *ForceUnitVectorToAxisX, double *ForceUnitVectorToAxisY, double *ForceUnitVectorToAxisZ)
{
    float dx, dy, dz;
    double unitVectorToAxisX, unitVectorToAxisY, unitVectorToAxisZ;
    dx = fabs(mvpStateVars->GetXFromIndex(index1) - mvpStateVars->GetXFromIndex(index2));
    dy = fabs(mvpStateVars->GetYFromIndex(index1) - mvpStateVars->GetYFromIndex(index2));
    dz = fabs(mvpStateVars->GetZFromIndex(index1) - mvpStateVars->GetZFromIndex(index2));
    if (mvUserSettings.mvGridType == 1)
    {
        if (dx > 0.2 * mvpStateVars->mvDx && dy < 0.2 * mvpStateVars->mvDx)
        {
            unitVectorToAxisX = 1.;
            unitVectorToAxisY = 0.;
            unitVectorToAxisZ = 0.;
        }
        else
        { // in 2d hexagonal all von Neumann boundaries (besides two first) have an angle of 60deg to horizontal axis
            unitVectorToAxisX = 1. / 2.;
            unitVectorToAxisY = sqrt(3.) / 2.;
            unitVectorToAxisZ = 0.;
        }
    }
    else
    {
        if (dx > 0.2 * mvpStateVars->mvDx && dy < 0.2 * mvpStateVars->mvDx && (mvpStateVars->mvIs2DSimulation || dz < 0.5 * mvpStateVars->mvDx))
        {
            unitVectorToAxisX = 1.;
            unitVectorToAxisY = 0.;
            unitVectorToAxisZ = 0.;
        }
        if (dx < 0.5 * mvpStateVars->mvDx && dy > 0.5 * mvpStateVars->mvDx && (mvpStateVars->mvIs2DSimulation || dz < 0.5 * mvpStateVars->mvDx))
        {
            unitVectorToAxisX = 0.;
            unitVectorToAxisY = 1.;
            unitVectorToAxisZ = 0.;
        }
        if (dx < 0.5 * mvpStateVars->mvDx && dy < 0.5 * mvpStateVars->mvDx && (!mvpStateVars->mvIs2DSimulation && dz > 0.5 * mvpStateVars->mvDx))
        {
            unitVectorToAxisX = 0.;
            unitVectorToAxisY = 0.;
            unitVectorToAxisZ = 1.;
        }
        if (dx > 0.5 * mvpStateVars->mvDx && dy > 0.5 * mvpStateVars->mvDx && (mvpStateVars->mvIs2DSimulation || dz < 0.5 * mvpStateVars->mvDx))
        {
            unitVectorToAxisX = 1. / sqrt(2.);
            unitVectorToAxisY = 1. / sqrt(2.);
            unitVectorToAxisZ = 0.;
        }
        if (dx > 0.5 * mvpStateVars->mvDx && dy < 0.5 * mvpStateVars->mvDx && (!mvpStateVars->mvIs2DSimulation && dz > 0.5 * mvpStateVars->mvDx))
        {
            unitVectorToAxisX = 1. / sqrt(2.);
            unitVectorToAxisY = 0.;
            unitVectorToAxisZ = 1. / sqrt(2.);
        }
        if (dx < 0.5 * mvpStateVars->mvDx && dy > 0.5 * mvpStateVars->mvDx && (!mvpStateVars->mvIs2DSimulation && dz > 0.5 * mvpStateVars->mvDx))
        {
            unitVectorToAxisX = 0.;
            unitVectorToAxisY = 1. / sqrt(2.);
            unitVectorToAxisZ = 1. / sqrt(2.);
        }
        if (dx > 0.5 * mvpStateVars->mvDx && dy > 0.5 * mvpStateVars->mvDx && (!mvpStateVars->mvIs2DSimulation && dz > 0.5 * mvpStateVars->mvDx))
        {
            unitVectorToAxisX = 1. / sqrt(3.);
            unitVectorToAxisY = 1. / sqrt(3.);
            unitVectorToAxisZ = 1. / sqrt(3.);
        }
    }
    *ForceUnitVectorToAxisX = unitVectorToAxisX;
    *ForceUnitVectorToAxisY = unitVectorToAxisY;
    *ForceUnitVectorToAxisZ = unitVectorToAxisZ;
}

double Microstructure::GetAreaPerVolumeOfBoundary(int index1, int index2)
{
    double AreaPerVolForThisBoundary;
    // FIRST CHECK IF WE ARE IN TEST MODE, AND PARTICULARLY FOR THE "WRONG" V-N GRID
    if (mvUserSettings.mvGridType == 0 && mvpStateVars->mvIs2DSimulation && mvUserSettings.mvTestVonNeumannInSquare)
        return 1. / mvpStateVars->mvDx;
    // NOW WE ARE THEN USING THE CORRECT ONES
    if (mvUserSettings.mvGridType == 1)
        AreaPerVolForThisBoundary = mvpStateVars->mvMaxBoundaryAreaPerVol;
    else
    {
        double WeightedBoundaryArea = mvpStateVars->mvDx / mvpStateVars->GetDistanceBetweenCells(index1, index2);
        AreaPerVolForThisBoundary = WeightedBoundaryArea * mvpStateVars->mvMaxBoundaryAreaPerVol;
    }
    return AreaPerVolForThisBoundary;
}

/** @brief Performs the computation of all cells' Traka's reorientation rates to transform due to recrystallization and/or grain growth
 * Traka, K, et al, Comp. Mat. Science, 2021.
 *
 * @param latticeId the crystal structure simulated for rex&gg
 * @return the max re-orientation rate among all cells for this simulation step
 */
float Microstructure::CalcReorientationRatesForFullFieldRexOrGG(int latticeId)
{
    int index, indexNbrToGrow, indexNbrOther, indexDragged;
    int lall[26];
    int lall2[26];
    int OriIdToConsume; // ori id of (per loop) current cell (index) being examined for its potential re-orientation
    int OriIdToGrow;    // ori id of candidate cell to grow into cell index
    int OriIdNbr;       // ori id of any neighbour cell of index

    LOG_F(INFO, "calculating reorientation rates, number of cells so far where growth/rex initiated are %d", cells_rx);
    double MisTol = mvUserSettings.mvLowerMisorientationCutOff; // misorientation tolerance to consider two orientations the same (in degrees)

    double AreaPerVolumeForThisBoundary, IncSigma;                                                                                          // calculated per boundary i.e. to calculate DG, but also to calculate force (the area)
    double Component_X_OfUnitVectorOfForceAtBoundary, Component_Y_OfUnitVectorOfForceAtBoundary, Component_Z_OfUnitVectorOfForceAtBoundary; // calculated per boundary (only for migrating boundaries of each possible reorientation)
    double CurrentStateEnergyPerUnitVolume, EnergyPerVolForThisBoundaryOfCurrentState;                                                      // used when calculating current state
    double NextStateEnergyPerUnitVolume, EnergyPerVolForThisBoundaryOfNextState;                                                            // used when calculating next state
    double TotalForceX_PerVol, TotalForceY_PerVol, TotalForceZ_PerVol;                                                                      // components of reorientation force, calculated gradually (i.e. frst unit vectors of all migrating boundaries times areas, then multiply with driving pressure etc )
    double DG_ForReorientation, MobilityOfBoundaries, TotalForceForReOrientationPerVol, fDot;                                               // used when calculating finally reorientation rate of cell.

    float fDotMax = 0.;
    float fDotMaxICellsGroup = 0.;

    int CSLRelationship;
    int cellsMoving = 0;
    int iCells = 0;

    int NNeighbours = mvpStateVars->mvNAll;
    if (mvUserSettings.mvGridType == 0 && mvpStateVars->mvIs2DSimulation && mvUserSettings.mvTestVonNeumannInSquare)
    {
        NNeighbours = mvpStateVars->mvNNearest;
        LOG_F(WARNING, " KEEP IN MIND, YOU ARE USING GRAIN GROWTH, SQUARE GRID AND CONSIDERING ONLY NEAREST NEIGHBOURS. I HOPE YOU ARE JUST TESTING IF EVERYTHING WORKS AS IT SHOULD (i.e. NOT ACCEPTABLE GG DUE TO HIGH ANISOTROPY). OTHERWISE CHECK/CHANGE USER SETTINGS FOR TestVonNeumannInSquare OPTION TO ZERO.");
    }
    double mis[NNeighbours][NNeighbours]; // misorientation between pairs of neighbours, mis[i][i] between i-th neighbour and reorienting grain
    fDotMaxICellsGroup = 0.;

    for (int index : interfaceIndices)
    {
        if (mvpStateVars->GetLatticeIdOfCell(index) != latticeId)
            continue;
        iCells++;
        if (!mvUserSettings.mvAllowPeriodicBoundConditions && mvpStateVars->mvpBoundaryCell[index])
            continue;
        OriIdToConsume = mvpStateVars->GetOriIdOfCell(index);
        mvpStateVars->mvpOneOfTheNeighboursGrowingInto[index] = -1;
        mvpStateVars->mvpConsumptionRate[index] = 0;
        if (mvUserSettings.mvIncludeCSL19FastGrowth && mvpStateVars->mvpCSLRelationshipNow == nullptr)
            mvpStateVars->mvpCSLRelationshipNow = new int[mvpStateVars->mvN];
        if (mvUserSettings.mvIncludeCSL19FastGrowth)
            mvpStateVars->mvpCSLRelationshipNow[index] = 0;

        mvpStateVars->GetNeighbourCellOffsets(index, lall);
        CurrentStateEnergyPerUnitVolume = 0.; // boundary energy of cell i at time t (now)
        bool canGetConsumed = false;          // this is for now. To avoid calculating the while mis[nbr][nbr] in case index cell is not eligible to be consumed due to boundaries less than tolerance

        for (int nbr = 0; nbr < NNeighbours; nbr++)
        {
            indexNbrToGrow = lall[nbr];
            int OriIdToGrow = mvpStateVars->GetOriIdOfCell(indexNbrToGrow);
            mis[nbr][nbr] = 0.;
            if (!mvUserSettings.mvAllowPeriodicBoundConditions && (fabs(mvpStateVars->GetXFromIndex(index) - mvpStateVars->GetXFromIndex(lall[nbr])) > 5. * mvpStateVars->mvDx || fabs(mvpStateVars->GetYFromIndex(index) - mvpStateVars->GetYFromIndex(lall[nbr])) > 5. * mvpStateVars->mvDx))
                continue;

            if (OriIdToConsume != OriIdToGrow)
                mis[nbr][nbr] = MisorientationBetweenCells(index, indexNbrToGrow, true, false, false, &CSLRelationship);
            if (mis[nbr][nbr] >= MisTol)
                canGetConsumed = true;
        }
        if (!canGetConsumed)
            continue;

        for (int nbr = 0; nbr < NNeighbours; nbr++)
        {
            indexNbrToGrow = lall[nbr];
            OriIdToGrow = mvpStateVars->GetOriIdOfCell(indexNbrToGrow);
            mis[nbr][nbr] = 0.;
            if (OriIdToConsume != OriIdToGrow)
            {
                mis[nbr][nbr] = MisorientationBetweenCells(index, indexNbrToGrow, true, false, mvUserSettings.mvIncludeCSL19FastGrowth, &CSLRelationship);
                if (CSLRelationship > 0)
                {
                    mvpStateVars->mvpCSLRelationshipNow[index] = 19;
                    mvpStateVars->mvpCSLRelationshipNow[indexNbrToGrow] = 19;
                }
            }

            if (!mvUserSettings.mvAllowPeriodicBoundConditions && (fabs(mvpStateVars->GetXFromIndex(index) - mvpStateVars->GetXFromIndex(indexNbrToGrow)) > 5. * mvpStateVars->mvDx || fabs(mvpStateVars->GetYFromIndex(index) - mvpStateVars->GetYFromIndex(indexNbrToGrow)) > 5. * mvpStateVars->mvDx))
                mis[nbr][nbr] = 0.;

            for (int nbr1 = 0; nbr1 < NNeighbours; nbr1++)
            {
                if (nbr1 == nbr)
                    continue;
                indexNbrOther = lall[nbr1];

                if (!mvUserSettings.mvAllowPeriodicBoundConditions && (fabs(mvpStateVars->GetXFromIndex(indexNbrOther) - mvpStateVars->GetXFromIndex(indexNbrToGrow)) > 5. * mvpStateVars->mvDx || fabs(mvpStateVars->GetYFromIndex(indexNbrOther) - mvpStateVars->GetYFromIndex(indexNbrToGrow)) > 5. * mvpStateVars->mvDx))
                    mis[nbr][nbr1] = 0.;
                OriIdNbr = mvpStateVars->GetOriIdOfCell(indexNbrOther);
                if (nbr1 < nbr)
                {
                    mis[nbr][nbr1] = mis[nbr1][nbr];
                    continue;
                }
                if (OriIdToGrow == OriIdNbr)
                {
                    mis[nbr][nbr1] = 0.;
                    continue;
                }
                mis[nbr][nbr1] = MisorientationBetweenCells(indexNbrToGrow, indexNbrOther, true, false, false, &CSLRelationship);
            }

            if (mis[nbr][nbr] < MisTol)
                continue; // almost same orientation, still not a boundary
            IncSigma = GetBoundaryEnergyForThisMisorientation(mis[nbr][nbr]);
            AreaPerVolumeForThisBoundary = GetAreaPerVolumeOfBoundary(index, indexNbrToGrow);
            EnergyPerVolForThisBoundaryOfCurrentState = AreaPerVolumeForThisBoundary * IncSigma;
            CurrentStateEnergyPerUnitVolume += EnergyPerVolForThisBoundaryOfCurrentState;
        }

        for (int nbr = 0; nbr < NNeighbours; nbr++)
        {
            NextStateEnergyPerUnitVolume = 0.;
            DG_ForReorientation = 0.;
            TotalForceX_PerVol = 0.;
            TotalForceY_PerVol = 0.;
            TotalForceZ_PerVol = 0.;

            indexNbrToGrow = lall[nbr];
            if (!mvUserSettings.mvAllowPeriodicBoundConditions && (fabs(mvpStateVars->GetXFromIndex(index) - mvpStateVars->GetXFromIndex(indexNbrToGrow)) > 5. * mvpStateVars->mvDx || fabs(mvpStateVars->GetYFromIndex(index) - mvpStateVars->GetYFromIndex(indexNbrToGrow)) > 5. * mvpStateVars->mvDx))
                continue;

            if (mvpStateVars->mvpRX[indexNbrToGrow] == 0 && (mvpStateVars->mvpCI[index] < mvUserSettings.mvMinimumCI || mvpStateVars->mvpCI[indexNbrToGrow] < mvUserSettings.mvMinimumCI || mvpStateVars->mvpHasNeighNonInd[index] == 1 || mvpStateVars->mvpHasNeighNonInd[indexNbrToGrow] == 1))
                continue;
            if (mis[nbr][nbr] < MisTol)
                continue; // not a candidate

            AreaPerVolumeForThisBoundary = GetAreaPerVolumeOfBoundary(index, indexNbrToGrow);
            CalculateDirectionOfForceAtThisBoundary(index, indexNbrToGrow, &Component_X_OfUnitVectorOfForceAtBoundary, &Component_Y_OfUnitVectorOfForceAtBoundary, &Component_Z_OfUnitVectorOfForceAtBoundary);
            TotalForceX_PerVol += (Component_X_OfUnitVectorOfForceAtBoundary * AreaPerVolumeForThisBoundary);
            TotalForceY_PerVol += (Component_Y_OfUnitVectorOfForceAtBoundary * AreaPerVolumeForThisBoundary);
            TotalForceZ_PerVol += (Component_Z_OfUnitVectorOfForceAtBoundary * AreaPerVolumeForThisBoundary);

            OriIdToGrow = mvpStateVars->GetOriIdOfCell(indexNbrToGrow);

            for (int nbr1 = 0; nbr1 < NNeighbours; nbr1++)
            {
                if (nbr1 == nbr)
                    continue;
                indexNbrOther = lall[nbr1];

                if (mis[nbr][nbr1] < MisTol)
                {
                    if (mvpStateVars->mvpRX[indexNbrOther] == 1 || (mvpStateVars->mvpHasNeighNonInd[indexNbrOther] == 0 && mvpStateVars->mvpCI[indexNbrOther] >= mvUserSettings.mvMinimumCI))
                    {
                        AreaPerVolumeForThisBoundary = GetAreaPerVolumeOfBoundary(index, indexNbrOther);
                        CalculateDirectionOfForceAtThisBoundary(index, indexNbrOther, &Component_X_OfUnitVectorOfForceAtBoundary, &Component_Y_OfUnitVectorOfForceAtBoundary, &Component_Z_OfUnitVectorOfForceAtBoundary);
                        TotalForceX_PerVol += (Component_X_OfUnitVectorOfForceAtBoundary * AreaPerVolumeForThisBoundary);
                        TotalForceY_PerVol += (Component_Y_OfUnitVectorOfForceAtBoundary * AreaPerVolumeForThisBoundary);
                        TotalForceZ_PerVol += (Component_Z_OfUnitVectorOfForceAtBoundary * AreaPerVolumeForThisBoundary);
                    }
                }
                else
                {
                    IncSigma = GetBoundaryEnergyForThisMisorientation(mis[nbr][nbr1]);
                    AreaPerVolumeForThisBoundary = GetAreaPerVolumeOfBoundary(index, indexNbrOther);
                    EnergyPerVolForThisBoundaryOfNextState = AreaPerVolumeForThisBoundary * IncSigma;

                    mvpStateVars->GetNeighbourCellOffsets(indexNbrOther, lall2);
                    indexDragged = lall2[nbr];

                    if (mvpStateVars->mvpRX[indexDragged] == 1 || (mvpStateVars->mvpCI[indexNbrOther] >= mvUserSettings.mvMinimumCI && mvpStateVars->mvpCI[indexDragged] >= mvUserSettings.mvMinimumCI && mvpStateVars->mvpHasNeighNonInd[indexNbrOther] == 0 && mvpStateVars->mvpHasNeighNonInd[indexDragged] == 0))
                    {
                        if (fabs(MisorientationBetweenCells(indexNbrToGrow, indexDragged, true, false, false, &CSLRelationship)) < MisTol && fabs(MisorientationBetweenCells(indexNbrOther, indexDragged, true, false, false, &CSLRelationship)) >= MisTol)
                        {
                            EnergyPerVolForThisBoundaryOfNextState *= (1. - mvUserSettings.mvRexGGParameterC);
                        }
                    }
                    NextStateEnergyPerUnitVolume += EnergyPerVolForThisBoundaryOfNextState;
                }
            }

            DG_ForReorientation = CurrentStateEnergyPerUnitVolume - NextStateEnergyPerUnitVolume;
            if (!IsTheNumberFinite(DG_ForReorientation))
                LOG_F(ERROR, "DG_ForReorientation %d with nbr %d is nan, DG", index, nbr);
            if (!IsTheNumberFinite(TotalForceX_PerVol) || !IsTheNumberFinite(TotalForceY_PerVol) || !IsTheNumberFinite(TotalForceZ_PerVol))
                LOG_F(ERROR, "Forces also have NaN values");

            if (DG_ForReorientation <= 0.)
                continue;
            TotalForceX_PerVol *= DG_ForReorientation;
            TotalForceY_PerVol *= DG_ForReorientation;
            TotalForceZ_PerVol *= DG_ForReorientation;
            TotalForceForReOrientationPerVol = sqrt(pow(TotalForceX_PerVol, 2.) + pow(TotalForceY_PerVol, 2.) + pow(TotalForceZ_PerVol, 2.));

            MobilityOfBoundaries = GetMobilityForThisMisorientation(mvpStateVars->GetLatticeIdOfCell(index), mis[nbr][nbr]);

            fDot = TotalForceForReOrientationPerVol * MobilityOfBoundaries;
            fDot *= BetaFactorForInterfaceMigrationRates;
            if (mvUserSettings.mvIncludeCSL19FastGrowth && mvpStateVars->mvpCSLRelationshipNow[index] == 0)
                fDot *= 0.1;

            if (fDot <= 0.)
                continue;
            if (!IsTheNumberFinite(fDot))
                LOG_F(ERROR, "fdot of %d with nbr %d is nan, mob %f force %f Enow %f Enext %f",
                      index, nbr, MobilityOfBoundaries, TotalForceForReOrientationPerVol,
                      CurrentStateEnergyPerUnitVolume, CurrentStateEnergyPerUnitVolume);
            if (fDot > 0. && mvpStateVars->mvpConsumptionRate[index] <= 0.)
                cellsMoving++;

            if (fDot <= mvpStateVars->mvpConsumptionRate[index])
                continue;
            mvpStateVars->mvpConsumptionRate[index] = fDot;
            mvpStateVars->mvpOneOfTheNeighboursGrowingInto[index] = indexNbrToGrow;
        }

        if (mvpStateVars->mvpConsumptionRate[index] > fDotMaxICellsGroup)
            fDotMaxICellsGroup = mvpStateVars->mvpConsumptionRate[index];
    }

    LOG_IF_F(INFO, iCells > 0, "number of iCells %d for lattice id %d for max mob in [meter4  / J s]  %g sigma Max %g and kinetic factor beta %g fDotMax %g A/V %g cells moving %d",
             iCells,
             latticeId,
             GetMobilityForThisMisorientation(latticeId, mvUserSettings.mvHAGB),
             GetBoundaryEnergyForThisMisorientation(mvUserSettings.mvHAGB),
             BetaFactorForInterfaceMigrationRates,
             fDotMaxICellsGroup,
             mvpStateVars->mvMaxBoundaryAreaPerVol,
             cellsMoving);

    return fDotMaxICellsGroup;
}

int Microstructure::ReOrientCellsForFullFieldRexOrGG(int latticeId)
{
    int isize = 0;
    int moved = 0;
    int reoriented = 0;

    // Loop over interface cells
    for (int index : interfaceIndices)
    {
        if (mvpStateVars->GetLatticeIdOfCell(index) != latticeId)
            continue;
        isize++;
        if (mvpStateVars->mvpConsumptionRate[index] <= 0.)
            continue;
        moved++;
        mvpStateVars->mvpConsumedFraction[index] += mvpStateVars->mvpConsumptionRate[index] * mvTimeStep;

        if (mvpStateVars->mvpConsumedFraction[index] >= 1.)
        {
            if (mvpStateVars->mvpRX[index] == 0)
                frx++;
            int indexThatGrew = mvpStateVars->mvpOneOfTheNeighboursGrowingInto[index];
            SwitchCell(index, indexThatGrew, true);
            reoriented++;
            ConsumedCellsInStepBySamePhase.push_back(index);

            // Debug log
            // LOG_F(INFO, "Cell reoriented: index %d, indexThatGrew %d, reoriented %d", index, indexThatGrew, reoriented);
        }
    }

    // Log the summary of the step
    LOG_F(INFO, "mvTime %g steps %d dt %g frx %g and icells %d moved %d reoriented %d",
          mvTime, mvSimulationStep, mvTimeStep, frx, isize, moved, reoriented);

    return reoriented;
}

int Microstructure::TransformCellsFromInterfaceMigration()
{
    int isize = 0;
    int moved = 0;
    int consumed = 0;
    int BCC_grew = 0;
    int FCC_grew = 0;
    int transformedIntoBCC = 0;
    int transformedIntoFCC = 0;
    // Loop over interface cells
    for (int index : interphaseIndices)
    {
        isize++;
        if (mvpStateVars->mvpConsumptionRate[index] <= 0.)
            continue;
        moved++;
        mvpStateVars->mvpConsumedFraction[index] += mvpStateVars->mvpConsumptionRate[index] * mvTimeStep;
        if (mvpStateVars->IsCellFCC(mvpStateVars->mvpOneOfTheNeighboursGrowingInto[index]))
            FCC_grew++;
        else
            BCC_grew++;

        if (mvpStateVars->mvpConsumedFraction[index] >= 1.)
        {
            if (mvpStateVars->mvpRX[index] == 0)
                frx++;
            int indexThatGrew = mvpStateVars->mvpOneOfTheNeighboursGrowingInto[index];
            SwitchCell(index, indexThatGrew, false);
            consumed++;
            ConsumedCellsInStepByDifferentPhase.push_back(index);
            if (mvpStateVars->IsCellFCC(mvpStateVars->mvpOneOfTheNeighboursGrowingInto[index]))
                transformedIntoFCC++;
            else
                transformedIntoBCC++;
            // Debug log
            // LOG_F(INFO, "Cell reoriented: index %d, indexThatGrew %d, reoriented %d", index, indexThatGrew, reoriented);
        }
    }

    grewFCC += transformedIntoFCC;
    grewBCC += transformedIntoBCC;

    LOG_F(INFO, "mvTime %g dt %g with %d cells being consumed by BCC, and %d cells just transformed into BCC, with %d cells being consumed by FCC, and %d cells just transformed into FCC ",
          mvTime,
          mvTimeStep,
          BCC_grew,
          transformedIntoBCC,
          FCC_grew,
          transformedIntoFCC);
    return consumed;
}

double Microstructure::CalculateDiffusivityMatrixAndReturnMaxDiffusivityInStep()
{
    double maxDiffusivity = 0.;
    double maxDcFCC = 0.;
    int phaseOfMaxD;
    double Xc1OfMaxD, Xc2OfMaxD;
    int nbrList[26];
    double XcFreeIndex, XcFreeIndexNbr;
    bool isInterphase;
    double XcMaxAllowedToIncreaseAgrenDiffusivityInFCC = 0.2;
    if (!mvUserSettings.mvConcentrationDependentDiffusivityInAustenite)
    {
        // No need to recalculate D matrix (it has the value of 1 for all pairs of cells, BUT we still need to know the max diffusivity (e.g. temperature changed from previous step)
        maxDiffusivity = mvpTCK->GetSoluteDiffusivityFromUser();
        LOG_F(INFO, "Max diffusivity in step is %e , for all phases (well FCC only matters here) and for all cells ", maxDiffusivity);
        return maxDiffusivity;
    }
    else
    {
        if (mvUserSettings.mvConcentrationDependentDiffusivityInAustenite)
        {
            LOG_IF_F(INFO, mvUserSettings.mvXcToUseConstantDiffusivityAgrenInAustenite >= 0.00001, "You have mart diffusivity equal to %e, and you also consider austenite diffusivity being xC dependent, but you take an effective value everywhere - i.e. Agren D for xC = %e ", mvpTCK->GetCDiffusivityMartensite(), mvUserSettings.mvXcToUseConstantDiffusivityAgrenInAustenite);
            LOG_IF_F(INFO, mvUserSettings.mvXcToUseConstantDiffusivityAgrenInAustenite < 0.00001, "You have mart diffusivity equal to %e, and you also consider austenite diffusivity being xC dependent from Agren (up to xC = %e at.frac.), and so max Dc in FCC is %e", mvpTCK->GetCDiffusivityMartensite(), XcMaxAllowedToIncreaseAgrenDiffusivityInFCC, mvpTCK->GetCDiffusivityAgren(XcMaxAllowedToIncreaseAgrenDiffusivityInFCC));
        }
        else
            LOG_F(INFO, "You have mart diffusivity equal to %e, and you also consider austenite diffusivity being different (but same for all FCC cells) ", mvpTCK->GetCDiffusivityMartensite());

        for (int index = 0; index < mvpStateVars->mvN; index++)
        {

            LOG_IF_F(ERROR, (*mvpStateVars->mvpXC)[index] < 0., " In Global, now that we are calculating D, I found cell of index %d, having XCTot %e, and mvpKappaFactorForCTrappedDefects %e ", index, (*mvpStateVars->mvpXC)[index], mvpStateVars->mvpKappaFactorForCTrappedDefects[index]);

            mvpStateVars->GetNeighbourCellOffsets(index, nbrList);
            for (int nbr = 0; nbr < mvpStateVars->mvNNearest; nbr++)
            {
                int indexNbr = nbrList[nbr];

                if (mvUserSettings.mvSoluteSegregationDislocations)
                {
                    XcFreeIndex = (*mvpStateVars->mvpXC)[index] - (*mvpStateVars->mvpXCTrapped)[index];
                    XcFreeIndexNbr = (*mvpStateVars->mvpXC)[indexNbr] - (*mvpStateVars->mvpXCTrapped)[indexNbr];
                }
                else
                {
                    XcFreeIndex = (*mvpStateVars->mvpXC)[index];
                    XcFreeIndexNbr = (*mvpStateVars->mvpXC)[indexNbr];
                }

                mvpStateVars->mvpDiffusivityC_Ratio[mvpStateVars->mvNNearest * index + nbr] = 1.; // for i and nbr in dissimilar phases it doesnt matter what we set in rhe Diffusivity matrix (they wont be considered - AX=b is only solved for phase-same elements)

                double PairDiffusivity = 0.;
                if (mvUserSettings.mvXcToUseConstantDiffusivityAgrenInAustenite < 0.00001)
                {
                    if (mvpStateVars->IsCellFCC(index))
                        PairDiffusivity += mvpTCK->GetCDiffusivityAgrenWithRestrictionForXC(XcFreeIndex, XcMaxAllowedToIncreaseAgrenDiffusivityInFCC);
                    else
                        PairDiffusivity += mvpTCK->GetCDiffusivityMartensite();
                    if (mvpStateVars->IsCellFCC(indexNbr))
                        PairDiffusivity += mvpTCK->GetCDiffusivityAgrenWithRestrictionForXC(XcFreeIndexNbr, XcMaxAllowedToIncreaseAgrenDiffusivityInFCC);
                    else
                        PairDiffusivity += mvpTCK->GetCDiffusivityMartensite();
                }
                else
                {
                    if (mvpStateVars->IsCellFCC(index))
                        PairDiffusivity += mvpTCK->GetCDiffusivityAgren(mvUserSettings.mvXcToUseConstantDiffusivityAgrenInAustenite);
                    else
                        PairDiffusivity += mvpTCK->GetCDiffusivityMartensite();
                    if (mvpStateVars->IsCellFCC(indexNbr))
                        PairDiffusivity += mvpTCK->GetCDiffusivityAgren(mvUserSettings.mvXcToUseConstantDiffusivityAgrenInAustenite);
                    else
                        PairDiffusivity += mvpTCK->GetCDiffusivityMartensite();
                }

                mvpStateVars->mvpDiffusivityC_Ratio[mvpStateVars->mvNNearest * index + nbr] = 0.5 * PairDiffusivity;
                // not yet ratio, it actually has the absolute diff value of the pair

                if (mvpStateVars->IsCellFCC(index) && mvpStateVars->IsCellFCC(indexNbr) && maxDcFCC < mvpStateVars->mvpDiffusivityC_Ratio[mvpStateVars->mvNNearest * index + nbr])
                    maxDcFCC = mvpStateVars->mvpDiffusivityC_Ratio[mvpStateVars->mvNNearest * index + nbr];

                if (mvpStateVars->mvpDiffusivityC_Ratio[mvpStateVars->mvNNearest * index + nbr] > maxDiffusivity)
                {
                    maxDiffusivity = mvpStateVars->mvpDiffusivityC_Ratio[mvpStateVars->mvNNearest * index + nbr];
                    phaseOfMaxD = mvpStateVars->GetLatticeIdOfCell(index);
                    Xc1OfMaxD = XcFreeIndex;
                    Xc2OfMaxD = XcFreeIndexNbr;
                    if (mvpStateVars->GetLatticeIdOfCell(index) != mvpStateVars->GetLatticeIdOfCell(indexNbr))
                        isInterphase = true;
                    else
                        isInterphase = false;
                }
                //            std::cout<<" In CarbonDiffusionStep Aust diffusivity is "<<PairDiffusivity<<" so ratio "<<mvpDiffusivityC_Ratio[mvNNearest*i+ nbr]<<std::endl;
            }
        }
    }
    LOG_F(INFO, "Max diffusivity in FCC %e", maxDcFCC);
    LOG_IF_F(INFO, !isInterphase, "Max diffusivity in step is %e , for phase  %d, and it corresponds to a pair of XCfree1 = %e, and XCfree2 = %e", maxDiffusivity, phaseOfMaxD, Xc1OfMaxD, Xc2OfMaxD);
    LOG_IF_F(INFO, isInterphase, "Max diffusivity in step is %e , for interphase cell, and it corresponds to a pair of XCfree1 = %e, and XCfree2 = %e", maxDiffusivity, Xc1OfMaxD, Xc2OfMaxD);
    LOG_F(INFO, "This means that max Rate factor is %e", maxDiffusivity * mvUserSettings.mvTimeStep / mvpStateVars->mvDx / mvpStateVars->mvDx);
    LOG_IF_F(WARNING, maxDiffusivity * mvUserSettings.mvTimeStep / mvpStateVars->mvDx / mvpStateVars->mvDx > 1., "Careful with your rate factor  >1 ,i.e.is %e", maxDiffusivity * mvUserSettings.mvTimeStep / mvpStateVars->mvDx / mvpStateVars->mvDx);

    // NOW go calclate the ratio
    for (int index = 0; index < mvpStateVars->mvN; index++)
    {
        mvpStateVars->GetNeighbourCellOffsets(index, nbrList);
        for (int nbr = 0; nbr < mvpStateVars->mvNNearest; nbr++)
        {
            mvpStateVars->mvpDiffusivityC_Ratio[mvpStateVars->mvNNearest * index + nbr] /= maxDiffusivity;
        }
    }
    LOG_F(INFO, "Finished setting diffusivity ratio matrix, max diffusivity is equal to %g", maxDiffusivity);
    return maxDiffusivity;
}

float Microstructure::CalcConsumptionRatesInterfaceMigrationDuringPartitioning()
{
    int index, indexNeigh, TypeInt, Nitypes = 4;
    int lall[26];
    LOG_F(INFO, "calculating consumption rates  (phase transformation)");
    float fDotMax = 0.;
    float fDotMaxICellsGroup = 0.;
    int cellsMoving = 0;
    int iCells = 0;
    fDotMax = 0.;
    double DG_Chem_InMartCells = 0.;
    double DG_Chem_InAustCells = 0.;
    double MartCellsUnderConsumption = 0.;
    double AustCellsUnderConsumption = 0.;
    double DeltaEnergyChemPerMole, DeltaEnergyChemPerVol, AreaPerVolumeForThisBoundary;                                                     // calculated per boundary to calculate force (the area)
    double Component_X_OfUnitVectorOfForceAtBoundary, Component_Y_OfUnitVectorOfForceAtBoundary, Component_Z_OfUnitVectorOfForceAtBoundary; // calculated per boundary (only for migrating boundaries of each possible reorientation)
    double TotalForceX_PerVol, TotalForceY_PerVol, TotalForceZ_PerVol;                                                                      // components of reorientation force, calculated gradually (i.e. frst unit vectors of all migrating boundaries times areas, then multiply with driving pressure etc )
    double DG_ForTransformationOfCell, MobilityOfBoundaries, TotalForceForTransformationPerVol, fDot;                                       // used when calculating finally reorientation rate of cell.
    int IndexToConsiderIfToBeConsumed, GrainIdToGiveIfToBeConsumed;
    double Vm;

    for (int index : interphaseIndices)
    {
        iCells++;
        if (!mvUserSettings.mvAllowPeriodicBoundConditions && mvpStateVars->mvpBoundaryCell[index])
            continue;

        fDot = 0.;
        DG_ForTransformationOfCell = 0.;
        TotalForceX_PerVol = 0.;
        TotalForceY_PerVol = 0.;
        TotalForceZ_PerVol = 0.;
        IndexToConsiderIfToBeConsumed = -1;
        GrainIdToGiveIfToBeConsumed = -1;
        mvpStateVars->mvpConsumptionRate[index] = -1.;
        double CurrentXC = (*mvpStateVars->mvpXC)[index];

        int NeighNumber = 0;
        double DeltaEnergyChemPerMoleForThisPair;
        double NeighXC;
        double BoundarySigma;
        DeltaEnergyChemPerMole = 0.;
        mvpStateVars->GetNeighbourCellOffsets(index, lall);
        double BoundaryEnergyPerVolOfCurrentState = 0.;
        double BoundaryEnergyPerVolOfNextState = 0.;
        for (int nbr = 0; nbr < NcellsForIntMig; nbr++)
        {
            indexNeigh = lall[nbr];
            if (mvpStateVars->IsCellFCC(indexNeigh) && mvpStateVars->mvpCI[indexNeigh] < 0.)
                continue;
            double Weight = mvpStateVars->mvDx / mvpStateVars->GetDistanceBetweenCells(index, indexNeigh);

            BoundarySigma = mvUserSettings.mvPhaseBoundaryEnergy;
            AreaPerVolumeForThisBoundary = GetAreaPerVolumeOfBoundary(index, indexNeigh);
            if (mvpStateVars->GetLatticeIdOfCell(index) == mvpStateVars->GetLatticeIdOfCell(indexNeigh))
                BoundaryEnergyPerVolOfNextState += AreaPerVolumeForThisBoundary * BoundarySigma;
            else
                BoundaryEnergyPerVolOfCurrentState += AreaPerVolumeForThisBoundary * BoundarySigma;

            if (mvpStateVars->GetLatticeIdOfCell(index) == mvpStateVars->GetLatticeIdOfCell(indexNeigh))
                continue;
            NeighNumber++;
            NeighXC = (*mvpStateVars->mvpXC)[indexNeigh];

            if (mvpStateVars->IsCellFCC(index))
                DeltaEnergyChemPerMoleForThisPair = mvpTCK->GetMuSubstitutionalInBCC(NeighXC) - mvpTCK->GetMuSubstitutionalInFCC(CurrentXC);
            else
                DeltaEnergyChemPerMoleForThisPair = mvpTCK->GetMuSubstitutionalInFCC(NeighXC) - mvpTCK->GetMuSubstitutionalInBCC(CurrentXC);

            DeltaEnergyChemPerMole += (Weight * DeltaEnergyChemPerMoleForThisPair);
        }
        if (NeighNumber > 0)
            DeltaEnergyChemPerMole = DeltaEnergyChemPerMole / NeighNumber;
        else
            continue;

        if (DeltaEnergyChemPerMole >= 0)
            continue; // cell only wants to transform if EphaseOther-EphaseNow<0

        Vm = (mvpStateVars->IsCellFCC(index)) ? mvpTCK->AusteniteMolarVolume() : mvpTCK->FerriteMolarVolume();
        // so convert in J/m3
        DeltaEnergyChemPerVol = DeltaEnergyChemPerMole / Vm;
        // Also for convention make DeltaEnergy positive since we know we have energy release;
        DeltaEnergyChemPerVol = fabs(DeltaEnergyChemPerVol);
        DG_ForTransformationOfCell = DeltaEnergyChemPerVol; // J/m3

        // std::cout<<" DeltaEnergyChemPerVol "<<DeltaEnergyChemPerVol<<" StrainEnergyInsidePerVol "<<StrainEnergyInsidePerVol<<std::endl;

        // no lets find which boundaries will exert force to migrate and the respective directions/magnitudes of force (depending on their area per vol)
        for (int nbr = 0; nbr < NcellsForIntMig; nbr++)
        {
            indexNeigh = lall[nbr];
            if (mvpStateVars->GetLatticeIdOfCell(index) == mvpStateVars->GetLatticeIdOfCell(indexNeigh))
                continue;
            if (mvpStateVars->IsCellFCC(indexNeigh) && mvpStateVars->mvpCI[indexNeigh] < 0.)
                continue;
            if (IndexToConsiderIfToBeConsumed < 0)
                IndexToConsiderIfToBeConsumed = indexNeigh;

            AreaPerVolumeForThisBoundary = GetAreaPerVolumeOfBoundary(index, indexNeigh);
            CalculateDirectionOfForceAtThisBoundary(index, indexNeigh, &Component_X_OfUnitVectorOfForceAtBoundary, &Component_Y_OfUnitVectorOfForceAtBoundary, &Component_Z_OfUnitVectorOfForceAtBoundary);
            TotalForceX_PerVol += (Component_X_OfUnitVectorOfForceAtBoundary * AreaPerVolumeForThisBoundary); // 1/m
            TotalForceY_PerVol += (Component_Y_OfUnitVectorOfForceAtBoundary * AreaPerVolumeForThisBoundary);
            TotalForceZ_PerVol += (Component_Z_OfUnitVectorOfForceAtBoundary * AreaPerVolumeForThisBoundary);
        }

        TotalForceX_PerVol *= DG_ForTransformationOfCell; // now finally we have N/m3 since we did J/m3 * 1/m
        TotalForceY_PerVol *= DG_ForTransformationOfCell;
        TotalForceZ_PerVol *= DG_ForTransformationOfCell;

        if (mvpStateVars->IsCellFCC(index))
        {
            DG_Chem_InAustCells += DeltaEnergyChemPerVol;
            AustCellsUnderConsumption += 1.;
        }
        else
        {
            DG_Chem_InMartCells += DeltaEnergyChemPerVol;
            MartCellsUnderConsumption += 1.;
        }

        TotalForceForTransformationPerVol = sqrt(pow(TotalForceX_PerVol, 2.) + pow(TotalForceY_PerVol, 2.) + pow(TotalForceZ_PerVol, 2.)); // N/m3

        MobilityOfBoundaries = mvpTCK->GetPhaseBoundaryMobility(); // mol m/J sec
        MobilityOfBoundaries *= Vm;                                // convert directly to m4/Jsec
        if (IsTheNumberFinite(TotalForceForTransformationPerVol))
            fDot = TotalForceForTransformationPerVol * MobilityOfBoundaries; // (N/m3) * (m4/Jsec) = 1/sec
        else
            fDot = 0.;

        if (fDot <= 0.)
            continue;

        cellsMoving++; // nothing, just to display during simulation
        if (fDot <= mvpStateVars->mvpConsumptionRate[index])
            continue;
        mvpStateVars->mvpConsumptionRate[index] = fDot;
        mvpStateVars->mvpOneOfTheNeighboursGrowingInto[index] = IndexToConsiderIfToBeConsumed;
        if (mvpStateVars->mvpConsumptionRate[index] > fDotMaxICellsGroup)
            fDotMaxICellsGroup = mvpStateVars->mvpConsumptionRate[index];
    }

    LOG_IF_F(INFO, iCells > 0, "number of iCells %d for max mob [j*m/(mol*s)] %g max mob [j*m4/s] %g and kinetic factor beta %g fDotMax %g A/V %g cells moving %d",
             iCells,
             mvpTCK->GetPhaseBoundaryMobility(),
             mvpTCK->GetPhaseBoundaryMobility() * Vm,
             BetaFactorForInterfaceMigrationRates,
             fDotMaxICellsGroup,
             mvpStateVars->mvMaxBoundaryAreaPerVol,
             cellsMoving);

    if (AustCellsUnderConsumption > 0.)
        Mean_DG_Chem_InFCC_Cells[Mean_DG_Chem_InFCC_Cells.size() - 1] = DG_Chem_InAustCells / AustCellsUnderConsumption;
    if (MartCellsUnderConsumption > 0.)
        Mean_DG_Chem_InBCC_Cells[Mean_DG_Chem_InBCC_Cells.size() - 1] = DG_Chem_InMartCells / MartCellsUnderConsumption;

    if (fDotMaxICellsGroup > fDotMax)
        fDotMax = fDotMaxICellsGroup;
    return fDotMax;
}

void Microstructure::SwitchCell(int index, int indexThatGrew, bool IsRexOrGG)
{
    int CSLRelationship;
    mvpStateVars->SetOriIdOfCell(index, mvpStateVars->GetOriIdOfCell(indexThatGrew));
    mvpStateVars->SetLatticeIdOfCell(index, mvpStateVars->GetLatticeIdOfCell(indexThatGrew));
    mvpStateVars->mvpConsumedFraction[index] = 0.;
    mvpStateVars->mvpOneOfTheNeighboursGrowingInto[index] = -1;
    if (IsRexOrGG)
    {
        mvpStateVars->mvpRX[index] = 1;
        double gb = MisorientationBetweenCells(index, indexThatGrew, false, false, false, &CSLRelationship);
        if (mvpStateVars->mvpMaxGbPassed[index] < gb)
            mvpStateVars->mvpMaxGbPassed[index] = gb;
    }
}

void Microstructure::UpdateInterfaceCells()
{
    for (int i = 0; i < ConsumedCellsInStepBySamePhase.size(); i++)
        UpdateThisInterfaceCellAndNeighbours(ConsumedCellsInStepBySamePhase[i]);
}

void Microstructure::UpdateInterphaseCells()
{
    for (int i = 0; i < ConsumedCellsInStepByDifferentPhase.size(); i++)
        UpdateThisInterphaseCellAndNeighbours(ConsumedCellsInStepByDifferentPhase[i]);
}

void Microstructure::UpdateThisInterfaceCellAndNeighbours(int indexThatGotConsumed)
{
    bool IsICell;
    int indexNbr;
    // check first own cell (it can happen if the consumed cell was itself a subgrain/grain which is now eliminated)
    IsICell = IsCellIndeedInterface(indexThatGotConsumed);
    if (!IsICell && mvpStateVars->IsCellInterface(indexThatGotConsumed))
        RemoveInterfaceIndex(indexThatGotConsumed);

    // now check all its neighborhood
    int NbrList[26];
    mvpStateVars->GetNeighbourCellOffsets(indexThatGotConsumed, NbrList);
    for (int nbr = 0; nbr < mvpStateVars->mvNAll; nbr++)
    {
        indexNbr = NbrList[nbr];
        IsICell = IsCellIndeedInterface(indexNbr);
        if (!IsICell && mvpStateVars->IsCellInterface(indexNbr))
            RemoveInterfaceIndex(indexNbr);
        if (IsICell && !mvpStateVars->IsCellInterface(indexNbr))
            AddInterfaceIndex(indexNbr);
    }
}

void Microstructure::UpdateThisInterphaseCellAndNeighbours(int indexThatGotConsumed)
{
    bool IsICell;
    int indexNbr;
    // check first own cell (it can happen if the consumed cell was itself a subgrain/grain which is now eliminated)
    IsICell = IsCellIndeedInterphase(indexThatGotConsumed);
    if (!IsICell && mvpStateVars->IsCellInterphase(indexThatGotConsumed))
        RemoveInterphaseIndex(indexThatGotConsumed);

    // now check all its neighborhood
    int NbrList[26];
    mvpStateVars->GetNeighbourCellOffsets(indexThatGotConsumed, NbrList);
    for (int nbr = 0; nbr < mvpStateVars->mvNAll; nbr++)
    {
        indexNbr = NbrList[nbr];
        IsICell = IsCellIndeedInterphase(indexNbr);
        if (!IsICell && mvpStateVars->IsCellInterphase(indexNbr))
            RemoveInterphaseIndex(indexNbr);
        if (IsICell && !mvpStateVars->IsCellInterphase(indexNbr))
            AddInterphaseIndex(indexNbr);
    }
}

float Microstructure::GetCellsMisorientationDarkeningFactor(int index)
{
    int indexNbr;
    int NbrList[26];
    int CSLRelationship;
    mvpStateVars->GetNeighbourCellOffsets(index, NbrList);
    float MaxMis = 0.;
    for (int nbr = 0; nbr < mvpStateVars->mvNAll; nbr++)
    {
        indexNbr = NbrList[nbr];
        if (mvpStateVars->GetOriIdOfCell(index) == mvpStateVars->GetOriIdOfCell(indexNbr))
            continue;
        float mis = MisorientationBetweenCells(index, indexNbr, true, false, false, &CSLRelationship);
        if (mis > mvUserSettings.mvLowerMisorientationCutOff && mis > MaxMis)
            MaxMis = mis;
        if (MaxMis > mvUserSettings.mvHAGB)
            break; // no reason to continue - we already have max misorentation
    }
    float DarkFactor = (MaxMis - mvUserSettings.mvLowerMisorientationCutOff) / (mvUserSettings.mvHAGB - mvUserSettings.mvLowerMisorientationCutOff);

    return DarkFactor;
}

bool Microstructure::IsCellIndeedInterface(int index)
{
    int indexNbr;
    int NbrList[26];
    int CSLRelationship;
    mvpStateVars->GetNeighbourCellOffsets(index, NbrList);
    for (int nbr = 0; nbr < mvpStateVars->mvNAll; nbr++)
    {
        indexNbr = NbrList[nbr];
        if (mvpStateVars->GetOriIdOfCell(index) == mvpStateVars->GetOriIdOfCell(indexNbr))
            continue;
        float mis = MisorientationBetweenCells(index, indexNbr, true, false, false, &CSLRelationship);
        if (mis > mvUserSettings.mvLowerMisorientationCutOff)
            return true;
    }
    return false;
}

bool Microstructure::IsCellIndeedInterphase(int index)
{
    int indexNbr;
    int NbrList[26];
    mvpStateVars->GetNeighbourCellOffsets(index, NbrList);
    for (int nbr = 0; nbr < mvpStateVars->mvNAll; nbr++)
    {
        indexNbr = NbrList[nbr];
        if (mvpStateVars->GetLatticeIdOfCell(index) == mvpStateVars->GetLatticeIdOfCell(indexNbr))
            continue;
        else
            return true;
    }
    return false;
}

/**
 * @brief Updates set of (sub)grain boundary cells by adding the cell of index, and updates cell's (index) identity regarding whether it is interface
 *
 * @param index the cell's index that will now be recognized as (sub)grain boundary cell
 */
void Microstructure::AddInterfaceIndex(int index)
{
    interfaceIndices.insert(index);
    mvpStateVars->SetCellAsInterface(index);
}

/**
 * @brief Updates set of (sub)grain boundary cells by removing the cell of index, and updates cell's (index) identity regarding whether it is interface
 *
 * @param index the cell's index that will no longer be recognized as (sub)grain boundary cell
 */
void Microstructure::RemoveInterfaceIndex(int index)
{
    interfaceIndices.erase(index);
    mvpStateVars->SetCellAsNonInterface(index);
}

/**
 * @brief Updates set of phase boundary cells by adding the cell of index, and updates cell's (index) identity regarding whether it is interface
 *
 * @param index the cell's index that will now be recognized as phase cell
 */
void Microstructure::AddInterphaseIndex(int index)
{
    interphaseIndices.insert(index);
    mvpStateVars->SetCellAsInterphase(index);
}

/**
 * @brief Updates set of phase boundary cells by removing the cell of index, and updates cell's (index) identity regarding whether it is interface
 *
 * @param index the cell's index that will now no longer recognized as phase cell
 */

void Microstructure::RemoveInterphaseIndex(int index)
{
    interphaseIndices.erase(index);
    mvpStateVars->SetCellAsNonInterphase(index);
}

/**
 * @brief Calculates the average carbon evolution behaviour in the different lattices as well as in interfaces, in defects, etc..
 * The purpose of this function is to save in the log file the quantities as well as to save them for later export
 */

void Microstructure::CarbonTemporalOverallQuantitiesForLaterAnalysis()
{

    if (mvUserSettings.mvAllowInterfaceMovementDuringPartitioning || mvUserSettings.mvIncludeCarbonInterphasePartitioning)
    {
        if (!mvUserSettings.mvThermodynamicsDataFilename.empty() && mvTime == 0)
        {
            LOG_F(INFO, "xC G_FCC G_BCC muSubsFCC muSubsBCC : ");
            for (int i = 0; i < 70; i++)
            {
                LOG_F(INFO, " %g %g %g", i * 0.01, mvpTCK->GetMuSubstitutionalInFCC(i * 0.01), mvpTCK->GetMuSubstitutionalInBCC(i * 0.01));
            }
            LOG_F(INFO, "And example: %g %g %g ", 0.124148, mvpTCK->GetMuSubstitutionalInFCC(0.124148), mvpTCK->GetMuSubstitutionalInBCC(0.124148));
        }
    }
    TimeEvol.push_back(mvTime);
    Mean_DG_Chem_InFCC_Cells.push_back(0.);
    Mean_DG_Chem_InBCC_Cells.push_back(0.);
    BCC_Growth.push_back(0.);
    FCC_Growth.push_back(0.);
    FCC_XC.push_back(0.);
    BCC_XC.push_back(0.);
    FCC_XC_Dispersion.push_back(0.);
    BCC_XC_Dispersion.push_back(0.);
    FCC_IcellsXC.push_back(0.);
    BCC_IcellsXC.push_back(0.);
    FCC_XCTrappedKappaFactor.push_back(0.);
    BCC_XCTrappedKappaFactor.push_back(0.);
    FCC_TrappedXC.push_back(0.);
    BCC_TrappedXC.push_back(0.);
    double XcTot = 0.;
    double XcTotTrapped = 0.;
    double XcKappaFactor = 0.;
    double XcTotTrappedBCC = 0.;
    double XcTotTrappedFCC = 0.;
    double XcKappaFactorBCC = 0.;
    double XcKappaFactorFCC = 0.;
    double XcTotBCC = 0.;
    double XcTotFCC = 0.;
    int CellsBCC = 0;
    int CellsFCC = 0;
    int NBCCInt = 0;
    int NFCCInt = 0;
    double TotMoleFrac = 0.;
    double TotMoleFracFCC = 0.;
    double TotMoleFracBCC = 0.;
    double TotMoleFracFCCICells = 0.;
    double TotMoleFracBCCICells = 0.;

    double DispersionCarbonFCC = 0.;
    double DispersionCarbonBCC = 0.;

    double XcBCCInt = 0.;
    double XcFCCInt = 0.;
    double MoleFracOfCell = 1.;

    for (int i = 0; i < mvpStateVars->mvN; i++)
    {
        MoleFracOfCell = 1.;
        if (mvUserSettings.mvConsiderConstantIronAtomsPerCell)
            MoleFracOfCell = (1. - mvUserSettings.mvAverageCarbon) / (1. - (*mvpStateVars->mvpXC)[i]);

        XcTot += MoleFracOfCell * (*mvpStateVars->mvpXC)[i];
        TotMoleFrac += MoleFracOfCell;
        if (mvUserSettings.mvSoluteSegregationDislocations)
        {
            XcKappaFactor += mvpStateVars->GetKappaFactorForCTrappedInCell(i);
            XcTotTrapped += MoleFracOfCell * mvpStateVars->GetXCTrappedInCell(i);
        }
        //                std::cout<<"time "<<mvTime<<" mvpStateVars->mvpXcDislo[i] "<<mvpStateVars->mvpXcDislo[i]<<" mvpStateVars->mvpCellFractionClustered[i] "<<mvpStateVars->mvpCellFractionClustered[i]<<std::endl;
        if (mvpStateVars->IsCellBCC(i))
        {
            CellsBCC++;
            TotMoleFracBCC += MoleFracOfCell;
            XcTotBCC += MoleFracOfCell * (*mvpStateVars->mvpXC)[i];
            if (mvUserSettings.mvSoluteSegregationDislocations)
            {
                XcKappaFactorBCC += mvpStateVars->GetKappaFactorForCTrappedInCell(i);
                XcTotTrappedBCC += MoleFracOfCell * mvpStateVars->GetXCTrappedInCell(i);
            }
        }
        if (mvpStateVars->IsCellFCC(i))
        {
            CellsFCC++;
            TotMoleFracFCC += MoleFracOfCell;
            XcTotFCC += MoleFracOfCell * (*mvpStateVars->mvpXC)[i];
            if (mvUserSettings.mvSoluteSegregationDislocations)
            {
                XcKappaFactorFCC += mvpStateVars->GetKappaFactorForCTrappedInCell(i);
                XcTotTrappedFCC += MoleFracOfCell * mvpStateVars->GetXCTrappedInCell(i);
            }
        }
        if (mvpStateVars->IsCellBCC(i))
        {
            if (mvpStateVars->IsCellInterphase(i))
            {
                NBCCInt++;
                TotMoleFracBCCICells += MoleFracOfCell;
                XcBCCInt += MoleFracOfCell * (*mvpStateVars->mvpXC)[i];
            }
        }
        if (mvpStateVars->IsCellFCC(i))
        {
            if (mvpStateVars->IsCellInterphase(i))
            {
                NFCCInt++;
                TotMoleFracFCCICells += MoleFracOfCell;
                XcFCCInt += MoleFracOfCell * (*mvpStateVars->mvpXC)[i];
            }
        }
    }

    for (int i = 0; i < mvpStateVars->mvN; i++)
    {
        MoleFracOfCell = 1.;
        if (mvUserSettings.mvConsiderConstantIronAtomsPerCell)
            MoleFracOfCell = (1. - mvUserSettings.mvAverageCarbon) / (1. - (*mvpStateVars->mvpXC)[i]);
        if (mvpStateVars->IsCellFCC(i))
            DispersionCarbonFCC += pow(MoleFracOfCell * (*mvpStateVars->mvpXC)[i] - XcTotFCC / TotMoleFracFCC, 2.);
        else
            DispersionCarbonBCC += pow(MoleFracOfCell * (*mvpStateVars->mvpXC)[i] - XcTotBCC / TotMoleFracBCC, 2.);
    }

    DispersionCarbonFCC /= TotMoleFracFCC;
    DispersionCarbonBCC /= TotMoleFracBCC;

    XcTot /= TotMoleFrac;
    XcTotTrapped /= TotMoleFrac;
    XcKappaFactor /= mvpStateVars->mvN;

    XcTotFCC /= TotMoleFracFCC;
    XcTotBCC /= TotMoleFracBCC;
    XcFCCInt /= TotMoleFracFCCICells;
    XcBCCInt /= TotMoleFracBCCICells;
    XcKappaFactorFCC /= CellsFCC;
    XcKappaFactorBCC /= CellsBCC;
    XcTotTrappedFCC /= TotMoleFracFCC;
    XcTotTrappedBCC /= TotMoleFracBCC;

    FCC_XC[FCC_XC.size() - 1] = XcTotFCC;
    BCC_XC[BCC_XC.size() - 1] = XcTotBCC;
    FCC_IcellsXC[FCC_IcellsXC.size() - 1] = XcFCCInt;
    BCC_IcellsXC[BCC_IcellsXC.size() - 1] = XcBCCInt;
    FCC_XCTrappedKappaFactor[FCC_XCTrappedKappaFactor.size() - 1] = XcKappaFactorFCC;
    BCC_XCTrappedKappaFactor[BCC_XCTrappedKappaFactor.size() - 1] = XcKappaFactorBCC;
    FCC_TrappedXC[FCC_TrappedXC.size() - 1] = XcTotTrappedFCC;
    BCC_TrappedXC[BCC_TrappedXC.size() - 1] = XcTotTrappedBCC;
    BCC_Growth[BCC_Growth.size() - 1] = grewBCC;
    FCC_Growth[FCC_Growth.size() - 1] = grewFCC;

    FCC_XC_Dispersion[FCC_XC_Dispersion.size() - 1] = DispersionCarbonFCC;
    BCC_XC_Dispersion[BCC_XC_Dispersion.size() - 1] = DispersionCarbonBCC;

    LOG_F(INFO, " Right before C diffusion step, at time %e seconds, I have average C %e , in my %d BCC ICells, and I have average C %e, in my %d FCC ICells ", mvTimeStep, XcBCCInt, NBCCInt, XcFCCInt, NFCCInt);
    LOG_F(INFO, " Now for %d BCC cells I have average C %e , and dispersion XC %e ", CellsBCC, XcTotBCC, DispersionCarbonBCC);
    LOG_F(INFO, " Now for %d FCC cells I have average C %e , and dispersion XC %e ", CellsFCC, XcTotFCC, DispersionCarbonFCC);
    LOG_F(INFO, " I also have in the whole RVE average C %e , average XC kappa factor  %e , and total Xc trapped from overall %e ", XcTot, XcKappaFactor, XcTotTrapped);
    LOG_F(INFO, " Now for %d BCC cells I have average C %e , average XC kappa factor %e , and total Xc trapped from overall %e ", CellsBCC, XcTotBCC, XcKappaFactorBCC, XcTotTrappedBCC);
    LOG_F(INFO, " Now for %d FCC cells have average C %e , average XC kappa factor %e , and total Xc trapped from overall %e ", CellsFCC, XcTotFCC, XcKappaFactorFCC, XcTotTrappedFCC);
    LOG_F(INFO, " Total BCC cells that have grown are %d ", grewBCC);
    LOG_F(INFO, " Total FCC cells that have grown are %d ", grewFCC);
}

/**
 * @brief Updates temperature based on simulation time and input settings, and informs thermodynamics object about it (which should NECESSARILY ALWAYS BE UPDATED)
 *
 */
void Microstructure::SetTemperatureAndInformInThermo()
{
    float timeFraction = mvTime / mvUserSettings.mvTimeTotal;
    float TemperatureRangeCovered = (mvUserSettings.mvEndTemperature - mvUserSettings.mvStartTemperature) * timeFraction;
    float CurrentTemperature = mvUserSettings.mvStartTemperature + TemperatureRangeCovered;
    mvpTCK->SetTemperatureForThisStep(CurrentTemperature);
}

void Microstructure::SimulationStep()
{

    ConsumedCellsInStepBySamePhase.clear();
    ConsumedCellsInStepByDifferentPhase.clear();
    SetTemperatureAndInformInThermo();
    LOG_F(INFO, "Simulation step %d time is %f temperature is %f", mvSimulationStep, mvTime, mvpTCK->GetTemperatureForThisStep());

    if (mvUserSettings.mvHasSoluteDiffusion)
    {
        mvTimeStep = mvUserSettings.mvTimeStep;
        mvMaxDiffusivityInTimeStep = CalculateDiffusivityMatrixAndReturnMaxDiffusivityInStep();
        // LOG_F(INFO, "mvMaxDiffusivityInTimeStep %g",mvMaxDiffusivityInTimeStep);
        mvpStateVars->SoluteDiffusionStep(mvpTCK, mvTimeStep, mvUserSettings.mvSoluteSegregationDislocations, mvMaxDiffusivityInTimeStep, mvUserSettings.mvIsPartitioningHappeningInDiffusionStep, mvUserSettings.mvIsSoluteSegregationHappeningInDiffusionStep);
    }

    if (mvUserSettings.mvHasSoluteDiffusion)
    {
        if ((mvSimulationStep == 0 || mvSimulationStep == 1) && mvStoreSoluteEvolutionTime < mvStoreSoluteEvolutionTimeInterval)
            CarbonTemporalOverallQuantitiesForLaterAnalysis();
        if (mvStoreSoluteEvolutionTime >= mvStoreSoluteEvolutionTimeInterval)
        {
            CarbonTemporalOverallQuantitiesForLaterAnalysis();
            mvStoreSoluteEvolutionTime = 0.;
        }
    }

    if (mvUserSettings.mvIsRexAndGG)
    {
        float mvMaxReRate = CalcReorientationRatesForFullFieldRexOrGG(1);
        mvTimeStep = 1.0 / mvMaxReRate;
        LOG_F(INFO, "Rex & GG max rate is %f and time step %f", mvMaxReRate, mvTimeStep);
        int NReoriented = ReOrientCellsForFullFieldRexOrGG(1);
        if (NReoriented)
        {
            UpdateInterfaceCells();
        }
    }
    if (mvUserSettings.mvAllowInterfaceMovementDuringPartitioning)
    {
        float mvMaxReRate = CalcConsumptionRatesInterfaceMigrationDuringPartitioning();
        LOG_F(INFO, "Phase transformation max rate is %f and time step %f", mvMaxReRate, mvTimeStep);
        int NReoriented = TransformCellsFromInterfaceMigration();
        if (NReoriented)
        {
            UpdateInterphaseCells();
        }
    }

    mvSimulationStep++;
    mvTime += mvTimeStep;
    mvCheckICellsTime++;
    mvTimePassedFromPreviousOutput += mvTimeStep;
    mvStoreSoluteEvolutionTime += mvTimeStep;

    if (mvCheckICellsTime >= mvEveryThatManySimStepsToCheckICells)
    {
        mvCheckICellsTime = 0;
        CheckInterfaceAndInterphaseCells();
    }

    if (mvTimePassedFromPreviousOutput >= mvUserSettings.mvTimeRangeToWriteOutput)
    {
        int MapType = 0;
        ExportColorCodedTIFF(mvUserSettings.mvOutputFolderPath + "/LatticeStructure_Time_" + std::to_string(int(mvTime)) + "_Map.tiff", MapType);

        if (mvUserSettings.mvIsRexAndGG)
            WriteMicrostructureRX(mvUserSettings.mvOutputFolderPath + "/MicroData_Time_" + std::to_string(int(mvTime)) + ".txt");
        else
        {
            MapType = 1;
            ExportColorCodedTIFF(mvUserSettings.mvOutputFolderPath + "/Carbon_Time_" + std::to_string(int(mvTime)) + "_Map.tiff", MapType);
            ExportColorbar(mvUserSettings.mvOutputFolderPath + "/Carbon_Time_" + std::to_string(int(mvTime)) + "_Colorbar.tiff");
            if (mvUserSettings.mvSoluteSegregationDislocations)
            {
                MapType = 2;
                ExportColorCodedTIFF(mvUserSettings.mvOutputFolderPath + "/CarbonTrapped_Time_" + std::to_string(int(mvTime)) + "_Map.tiff", MapType);
                ExportColorbar(mvUserSettings.mvOutputFolderPath + "/CarbonTrapped_Time_" + std::to_string(int(mvTime)) + "_Colorbar.tiff");
            }
        }

        mvTimePassedFromPreviousOutput = 0.0;
    }
}
