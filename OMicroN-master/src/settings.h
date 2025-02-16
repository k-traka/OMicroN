/* Settings.h
   OMicroN (optimising microstructures numerically) simulation program
   Header file containing the UserSettings class definitions and implementation
*/

#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem> // For creating directories
#include <stdexcept>  // For std::runtime_error
#include <cstdio>     // For std::remove

/** @brief Handles all the simulation parameters and sets them according to the input file and/or default values.
 *
 */
class UserSettings
{
public:
    /**
     * @brief Constructor that initializes settings from the input file.
     * @param inputFile Path to the input file containing simulation parameters.
     */
    UserSettings(const std::string &inputFile);

    /**
     * @brief Sets default values for all parameters. If a parameter is found in the input file, its value will be used instead.
     */
    void setDefaultValues();

    /**
     * @brief Writes the settings to the output file, including default values.
     * @param outputFile Path to the output file where settings will be saved.
     */
    void writeToFile(const std::string &outputFile);

    /**
     * @brief Prints the parameters as read from the input file.
     */
    void printParameters();

    /**
     * @brief Returns the path to the input (starting) microstructure file.
     * @return Path to the input microstructure file.
     */
    std::string getStartFileName() const
    {
        return mvStartFileName;
    }

    /**
     * @brief Returns the path to the folder where simulation outputs are stored.
     * @return Path to the output folder.
     */
    std::string getOutputFolderPath() const
    {
        return mvOutputFolderPath;
    }

    /**
     * @brief Returns the path to the file containing thermodynamic data.
     * @return Path to the thermodynamic data file.
     */
    std::string getThermodynamicDataFile() const
    {
        return mvThermodynamicsDataFilename;
    }

    /**
     * @brief Desired time step for solute redistribution.
     * Note: When recrystallization and grain growth are simulated, the time step will adapt to the maximum reorientation rate.
     * Note: When simulating solute redistribution the user must set it according to the highest expected diffusion rate, which in the numerical system is dt * maxDiffusivityInTimeStep / (mvDx * mvDx). So dt must be chosen such as dt * maxDiffusivityInTimeStep / (mvDx * mvDx) is lower than 1.
     */
    double mvTimeStep;

    /**
     * @brief Temperature at the beginning of the applied treatment (in K).
     */
    double mvStartTemperature;

    /**
     * @brief Temperature at the end of the applied treatment (in K).
     */
    double mvEndTemperature;

    /**
     * @brief Total simulation time (s).
     */
    double mvTimeTotal;

    /**
     * @brief Number of elements in the x direction.
     */
    int mvNx;

    /**
     * @brief Number of elements in the y direction.
     */
    int mvNy;

    /**
     * @brief Number of elements in the z direction.
     */
    int mvNz;

    /**
     * @brief Grid spacing (m).
     */
    double mvDx;

    /**
     * @brief Integer determining the type of grid used.
     * 0 is regular (square/cubic) and 1 is hexagonal.
     */
    int mvGridType;

    /**
     * @brief Flag indicating whether solute redistribution (e.g., diffusion, partitioning, trapping) should be simulated.
     */
    int mvHasSoluteDiffusion;


    /**
     * @brief Flag indicating whether solidification and solute redistribution (e.g., diffusion, partitioning, trapping) should be simulated.
     */
    int mvIsSolidification;

    /**
     * @brief Flag indicating whether deformation should be simulated.
     */
    int mvIsDeformation;
    /**
     * @brief Flag indicating whether solidification does not consider surface energy
     */
    int mvUseBoundaryEnergyIncrement;

    /**
     * @brief Flag indicating whether solidification does not consider anisotropic surface energy
     */
    int mvUseAnisotropicBoundaryEnergy;
    /**
     * @brief Flag indicating whether the interstitial redistribution and defect trapping is hydrogen (otherwise it is carbon) 
     *      */

int mvIsHydrogen;
    /**
     * @brief Flag indicating whether the interstitial redistribution and defect trapping is hydrogen refers to charging
     *      */

int mvIsHydrogenSource;
    /**
     * @brief Flag indicating whether solidification considers Solute Drag
     *      */
    int mvUseSoluteDrag;

    /**
     * @brief Flag indicating whether interface migration in the same phase (e.g., recrystallization and grain growth) should be simulated.
     */
    int mvIsRexAndGG;

    /**
     * @brief Misorientation below which any misorientation between pixels is considered noise and no orientation gradient is considered (e.g., 0.4 degrees).
     */
    float mvLowerMisorientationCutOff;

    /**
     * @brief A cell is not considered as an Icell if it has smaller misorientation from its neigbours than  (i.e. it CANNOT grow/be consumed) but may still consider as misoriented from its neighbours (i.e. if  mvLowerMisorientationCutOff< mis < mvLowerMisorientationForLAGB then the cell is deformed but it moves so slowly so it is not iterated for being itself consumed.
     */
    float mvLowerMisorientationForLAGB;

    /**
     * @brief Misorientation above which a high-angle grain boundary is considered. Energy density and mobility are constant above this value (e.g., 15 degrees).
     */
    float mvHAGB;

    /**
     * @brief Path to the input (starting) microstructure file.
     */
    std::string mvStartFileName;

    /**
     * @brief Path to the folder where simulation output files are written.
     */
    std::string mvOutputFolderPath;

    /**
     * @brief Path to the HDF5 file containing thermodynamic inputs for solute partitioning and/or phase transformations.
     */
    std::string mvThermodynamicsDataFilename;

    /**
     * @brief Flag indicating whether the input microstructure file also has Ori Id 
     */
    int mvReadOriId;
    /**
     * @brief Flag indicating whether the simulation involves microstructure evolution (considering imported state variables) or code testing.
     */
    int mvIsTestRVE;

    /**
     * @brief Flag indicating whether the code is tested for grain growth. Turn on means modeling grain growth in an unacceptable grid, useful for testing.
     */
    int mvTestVonNeumannInSquare;

    /**
     * @brief Flag indicating whether boundary conditions are enabled. Turning off is useful for test simulations of specific cases or recrystallization in very small RVEs.
     */
    int mvAllowPeriodicBoundConditions;

    /**
     * @brief Minimum confidence index for valid crystal orientations imported (useful when importing uncleaned EBSD microstructures).
     */
    float mvMinimumCI;

    /**
     * @brief Parameter C that calibrates the driving force for grain growth and recrystallization. Accepted values are 0.5, 0.7, 1.0.
     */
    float mvRexGGParameterC;
    /**
     * @brief Flag indicating that Rex or GG will be simulated in faster way (which may have error - though negligible)
     */
    int mvFasterRexAndGG;
    /**
     * @brief Flag indicating whether the dislocation density provided by the input file should be used in the simulation. Turn off means calculating local dislocation density based on orientation gradients.
     */
    int mvReadDislocationDensity;

    /**
     * @brief Flag indicating whether CSL19 boundaries (Ibe and Lucke) grow faster.
     */
    int mvIncludeCSL19FastGrowth;

    /**
     * @brief Flag indicating whether only a subset of orientations is considered in the microstructure evolution.
     */
    int mvIsOnlyAnOrientationSubsetConsidered;

    /**
     * @brief Flag indicating whether mass conservation for solutes refers to C/Fe or C/(Fe+C). C/Fe is more accurate, while C/(Fe+C) is commonly used.
     */
    int mvConsiderConstantIronAtomsPerCell;

    /**
     * @brief Flag indicating whether solute trapping to defects should be considered.
     */
    int mvSoluteSegregationDislocations;

    /**
     * @brief Mesoscale enrichment ratio at defects (e.g., k = 7*10^-15).
     */
    double mvKappaFactorForCsegInClusterArea;

    /**
     * @brief Flag indicating whether solute trapping to defects in the soft phase (e.g., austenite) should be considered.
     */
    int mvAllowDislocationsSegregationInAustenite;

    /**
     * @brief Flag indicating whether solute phase partitioning is to be simulated.
     */
    int mvIncludeCarbonInterphasePartitioning;

    /**
     * @brief Flag indicating whether solute trapping to defects is part of the numerical solution of diffusion. Recommended for mesh- and time step-independent simulation.
     */
    int mvIsSoluteSegregationHappeningInDiffusionStep;

    /**
     * @brief Flag indicating whether solute phase partitioning is part of the numerical solution of diffusion. Recommended for mesh- and time step-independent simulation.
     */
    int mvIsPartitioningHappeningInDiffusionStep;

    /**
     * @brief Flag indicating whether solute phase partitioning follows thermodynamic inputs for local concentrations. Turn off means using pre-defined analytical expressions.
     */
    int mvCarbonPartitioningFromInterpolatedSolutions;

    /**
     * @brief Flag indicating whether phase transformation simulation is enabled.
     */
    int mvAllowInterfaceMovementDuringPartitioning;

    /**
     * @brief Flag indicating whether diffusivity follows Agren's concentration-dependent expression (useful when simulating carbon).
     */
    int mvConcentrationDependentDiffusivityInAustenite;

    /**
     * @brief Carbon concentration that is used as "effective" value at which Agren's concentration-dependent expression gives diffusivity (useful when the simulation leads to high localized carbon contents which would then require very high time steps).
     */
    double mvXcToUseConstantDiffusivityAgrenInAustenite;

    /**
     * @brief Concentration of solute (carbon or other) that redistributes through the simulation (atomic fraction).
     */
    double mvAverageCarbon;

    /**
     * @brief Concentration of solute (carbon or other) that redistributes through the simulation (atomic fraction).
     */
    double mvAverageCopper;
        /**
     * @brief Concentration of solute (carbon or other) that redistributes through the simulation (atomic fraction).
     */
    double mvAverageTin;

    /**
     * @brief BCC lattice Burgers vector for the material investigated (m).
     */
    double mvBurgersVectorBCC;

    /**
     * @brief FCC lattice Burgers vector for the material investigated (m).
     */
    double mvBurgersVectorFCC;

    /**
     * @brief HCP lattice Burgers vector for the material investigated (m).
     */
    double mvBurgersVectorHCP;

    /**
     * @brief Time range to write output (s).
     */
    double mvTimeRangeToWriteOutput;

    /**
     * @brief Activation energy for solute diffusion (J/mol).
     */
    double mvSoluteDiffusivityActivationEnergy;

    /**
     * @brief Pre-exponential factor for solute diffusion (J/mol).
     */
    double mvSoluteDiffusivityPreFactor;

    /**
     * @brief Pre-exponential factor for grain boundary migration (J/mol).
     */
    double mvGrainBoundaryMobilityProExponential;

    /**
     * @brief Activation energy for grain boundary migration (J/mol).
     */
    double mvGrainBoundaryMobilityActivationEnergy;

    /**
     * @brief Grain boundary energy density (J/m^2).
     */
    double mvGrainBoundaryEnergy;

    /**
     * @brief Pre-exponential factor for phase boundary migration.
     */
    double mvPhaseBoundaryMobilityProExponential;

    /**
     * @brief Activation energy for phase boundary migration (J/mol).
     */
    double mvPhaseBoundaryMobilityActivationEnergy;

    /**
     * @brief Phase boundary energy density (J/m^2).
     */
    double mvPhaseBoundaryEnergy;

private:
    void readFromFile(const std::string &inputFile);
    void readStartFile(const std::string &startFileName);
    void createOutputFolder(const std::string &folderPath);
    void clearOutputFolder(const std::string &folderPath);
};

#endif // SETTINGS_H
