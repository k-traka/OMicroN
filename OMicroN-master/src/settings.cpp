#include "settings.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>

UserSettings::UserSettings(const std::string& inputFile) {
    setDefaultValues();
    readFromFile(inputFile);
    readStartFile(mvStartFileName); // Read start file using mvStartFileName from input file
    
    // Create the output folder if it doesn't exist
    createOutputFolder(mvOutputFolderPath);
}

void UserSettings::setDefaultValues() {
    mvTimeStep = 0.1;
    mvStartTemperature = 300.0;
    mvEndTemperature = 300.0;
    mvTimeTotal = 1000.0;
    mvNx = 10;
    mvNy = 10;
    mvNz = 10;
    mvDx = 1.0;
    mvGridType = 0;
    mvHasSoluteDiffusion = 0;
    mvIsRexAndGG = 1;
    mvLowerMisorientationCutOff = 0.3;
    mvHAGB = 15.;
    mvIsTestRVE = 0;
    mvTestVonNeumannInSquare = 0;
    mvAllowPeriodicBoundConditions = 1;
    mvMinimumCI = 0.1;
    mvRexGGParameterC = 0.7;
    mvReadDislocationDensity = 0;
    mvIncludeCSL19FastGrowth = 0;
    mvIsOnlyAnOrientationSubsetConsidered = 0;
    mvConsiderConstantIronAtomsPerCell = 0;
    mvSoluteSegregationDislocations = 0;
    mvKappaFactorForCsegInClusterArea = 0.0;
    mvAllowDislocationsSegregationInAustenite = 0;
    mvIncludeCarbonInterphasePartitioning = 0;
    mvIsSoluteSegregationHappeningInDiffusionStep = 0;
    mvIsPartitioningHappeningInDiffusionStep = 0;
    mvCarbonPartitioningFromInterpolatedSolutions = 0;
    mvAllowInterfaceMovementDuringPartitioning = 0;
    mvConcentrationDependentDiffusivityInAustenite = 0;
    mvXcToUseConstantDiffusivityAgrenInAustenite = 0.0;
    mvAverageCarbon = 0.0;
    mvGrainBoundaryEnergy = 0.0;
    mvPhaseBoundaryEnergy = 0.0;
    mvGrainBoundaryMobilityProExponential = 0.0;
    mvPhaseBoundaryMobilityProExponential = 0.0;
    mvGrainBoundaryMobilityActivationEnergy = 0.0;
    mvPhaseBoundaryMobilityActivationEnergy = 0.0;
    mvBurgersVectorBCC = 0.0;
    mvBurgersVectorFCC = 0.0;
    mvBurgersVectorHCP = 0.0;
    mvTimeRangeToWriteOutput = 10;
    mvSoluteDiffusivityPreFactor = 0.002;
    mvSoluteDiffusivityActivationEnergy = 125000;
}

void UserSettings::writeToFile(const std::string& outputFileName) {
    std::string outputFile = mvOutputFolderPath + "/" +  outputFileName;
    std::ofstream out(outputFile);
    if (!out.is_open()) {
        std::cerr << "Unable to open file: " << outputFile << std::endl;
        return;
    }

    out << "TimeStep=" << mvTimeStep << std::endl;
    out << "StartTemperature=" << mvStartTemperature << std::endl;
    out << "EndTemperature=" << mvEndTemperature << std::endl;
    out << "TimeTotal=" << mvTimeTotal << std::endl;
    out << "Nx=" << mvNx << std::endl;
    out << "Ny=" << mvNy << std::endl;
    out << "Nz=" << mvNz << std::endl;
    out << "Dx=" << mvDx << std::endl;
    out << "GridType=" << mvGridType << std::endl;
    out << "IsRexAndGG=" << mvIsRexAndGG << std::endl;
    out << "HasSoluteDiffusion=" << mvHasSoluteDiffusion << std::endl;
    out << "LowerMisorientationCutOff=" << mvLowerMisorientationCutOff << std::endl;
    out << "HAGB=" << mvHAGB << std::endl;
    out << "StartFileName=" << mvStartFileName << std::endl;
    out << "ThermodynamicsDataFilename=" << mvThermodynamicsDataFilename << std::endl;
    out << "OutputFolderPath=" << mvOutputFolderPath << std::endl;
    out << "IsTestRVE=" << mvIsTestRVE << std::endl;
    out << "TestVonNeumannInSquare=" << mvTestVonNeumannInSquare << std::endl;
    out << "AllowPeriodicBoundConditions=" << mvAllowPeriodicBoundConditions << std::endl;
    out << "MinimumCI=" << mvMinimumCI << std::endl;
    out << "RexGGParameterC=" << mvRexGGParameterC << std::endl;
    out << "ReadDislocationDensity=" << mvReadDislocationDensity << std::endl;
    out << "IncludeCSL19FastGrowth=" << mvIncludeCSL19FastGrowth << std::endl;
    out << "IsOnlyAnOrientationSubsetConsidered=" << mvIsOnlyAnOrientationSubsetConsidered << std::endl;
    out << "ConsiderConstantIronAtomsPerCell=" << mvConsiderConstantIronAtomsPerCell << std::endl;
    out << "SoluteSegregationDislocations=" << mvSoluteSegregationDislocations << std::endl;
    out << "KappaFactorForCsegInClusterArea=" << mvKappaFactorForCsegInClusterArea << std::endl;
    out << "AllowDislocationsSegregationInAustenite=" << mvAllowDislocationsSegregationInAustenite << std::endl;
    out << "IncludeCarbonInterphasePartitioning=" << mvIncludeCarbonInterphasePartitioning << std::endl;
    out << "IsSoluteSegregationHappeningInDiffusionStep=" << mvIsSoluteSegregationHappeningInDiffusionStep << std::endl;
    out << "IsPartitioningHappeningInDiffusionStep=" << mvIsPartitioningHappeningInDiffusionStep << std::endl;
    out << "CarbonPartitioningFromInterpolatedSolutions=" << mvCarbonPartitioningFromInterpolatedSolutions << std::endl;
    out << "AllowInterfaceMovementDuringPartitioning=" << mvAllowInterfaceMovementDuringPartitioning << std::endl;
    out << "ConcentrationDependentDiffusivityInAustenite=" << mvConcentrationDependentDiffusivityInAustenite << std::endl;
    out << "XcToUseConstantDiffusivityAgrenInAustenite=" << mvXcToUseConstantDiffusivityAgrenInAustenite << std::endl;
    out << "AverageCarbon=" << mvAverageCarbon << std::endl;
    out << "GrainBoundaryEnergy=" << mvGrainBoundaryEnergy << std::endl;
    out << "PhaseBoundaryEnergy=" << mvPhaseBoundaryEnergy << std::endl;
    out << "GrainBoundaryMobilityProExponential=" << mvGrainBoundaryMobilityProExponential << std::endl;
    out << "PhaseBoundaryMobilityProExponential=" << mvPhaseBoundaryMobilityProExponential << std::endl;
    out << "GrainBoundaryMobilityActivationEnergy=" << mvGrainBoundaryMobilityActivationEnergy << std::endl;
    out << "PhaseBoundaryMobilityActivationEnergy=" << mvPhaseBoundaryMobilityActivationEnergy << std::endl;
    out << "BurgersVectorBCC=" << mvBurgersVectorBCC << std::endl;
    out << "BurgersVectorFCC=" << mvBurgersVectorFCC << std::endl;
    out << "BurgersVectorHCP=" << mvBurgersVectorHCP << std::endl;
    out << "TimeRangeToWriteOutput=" << mvTimeRangeToWriteOutput << std::endl;
    out << "SoluteDiffusivityPreFactor=" << mvSoluteDiffusivityPreFactor << std::endl;
    out << "SoluteDiffusivityActivationEnergy=" << mvSoluteDiffusivityActivationEnergy << std::endl;
    out.close();
}

void UserSettings::printParameters() {
    std::cout << "TimeStep=" << mvTimeStep << std::endl;
    std::cout << "StartTemperature=" << mvStartTemperature << std::endl;
    std::cout << "EndTemperature=" << mvEndTemperature << std::endl;    std::cout << "TimeTotal=" << mvTimeTotal << std::endl;
    std::cout << "Nx=" << mvNx << std::endl;
    std::cout << "Ny=" << mvNy << std::endl;
    std::cout << "Nz=" << mvNz << std::endl;
    std::cout << "Dx=" << mvDx << std::endl;
    std::cout << "GridType=" << mvGridType << std::endl;
    std::cout << "IsRexAndGG=" << mvIsRexAndGG << std::endl;
    std::cout << "HasSoluteDiffusion=" << mvHasSoluteDiffusion << std::endl;
    std::cout << "LowerMisorientationCutOff=" << mvLowerMisorientationCutOff << std::endl;
    std::cout << "HAGB=" << mvHAGB << std::endl;
    std::cout << "StartFileName=" << mvStartFileName << std::endl;
    std::cout << "ThermodynamicsDataFilename=" << mvThermodynamicsDataFilename << std::endl;
    std::cout << "OutputFolderPath=" << mvOutputFolderPath << std::endl;
    std::cout << "IsTestRVE=" << mvIsTestRVE << std::endl;
    std::cout << "TestVonNeumannInSquare=" << mvTestVonNeumannInSquare << std::endl;
    std::cout << "AllowPeriodicBoundConditions=" << mvAllowPeriodicBoundConditions << std::endl;
    std::cout << "MinimumCI=" << mvMinimumCI << std::endl;
    std::cout << "RexGGParameterC=" << mvRexGGParameterC << std::endl;
    std::cout << "ReadDislocationDensity=" << mvReadDislocationDensity << std::endl;
    std::cout << "IncludeCSL19FastGrowth=" << mvIncludeCSL19FastGrowth << std::endl;
    std::cout << "IsOnlyAnOrientationSubsetConsidered=" << mvIsOnlyAnOrientationSubsetConsidered << std::endl;
    std::cout << "ConsiderConstantIronAtomsPerCell=" << mvConsiderConstantIronAtomsPerCell << std::endl;
    std::cout << "SoluteSegregationDislocations=" << mvSoluteSegregationDislocations << std::endl;
    std::cout << "KappaFactorForCsegInClusterArea=" << mvKappaFactorForCsegInClusterArea << std::endl;
    std::cout << "AllowDislocationsSegregationInAustenite=" << mvAllowDislocationsSegregationInAustenite << std::endl;
    std::cout << "IncludeCarbonInterphasePartitioning=" << mvIncludeCarbonInterphasePartitioning << std::endl;
    std::cout << "IsSoluteSegregationHappeningInDiffusionStep=" << mvIsSoluteSegregationHappeningInDiffusionStep << std::endl;
    std::cout << "IsPartitioningHappeningInDiffusionStep=" << mvIsPartitioningHappeningInDiffusionStep << std::endl;
    std::cout << "CarbonPartitioningFromInterpolatedSolutions=" << mvCarbonPartitioningFromInterpolatedSolutions << std::endl;
    std::cout << "AllowInterfaceMovementDuringPartitioning=" << mvAllowInterfaceMovementDuringPartitioning << std::endl;
    std::cout << "ConcentrationDependentDiffusivityInAustenite=" << mvConcentrationDependentDiffusivityInAustenite << std::endl;
    std::cout << "XcToUseConstantDiffusivityAgrenInAustenite=" << mvXcToUseConstantDiffusivityAgrenInAustenite << std::endl;
    std::cout << "AverageCarbon=" << mvAverageCarbon << std::endl;
    std::cout << "GrainBoundaryEnergy=" << mvGrainBoundaryEnergy << std::endl;
    std::cout << "PhaseBoundaryEnergy=" << mvPhaseBoundaryEnergy << std::endl;
    std::cout << "GrainBoundaryMobilityProExponential=" << mvGrainBoundaryMobilityProExponential << std::endl;
    std::cout << "PhaseBoundaryMobilityProExponential=" << mvPhaseBoundaryMobilityProExponential << std::endl;
    std::cout << "GrainBoundaryMobilityActivationEnergy=" << mvGrainBoundaryMobilityActivationEnergy << std::endl;
    std::cout << "PhaseBoundaryMobilityActivationEnergy=" << mvPhaseBoundaryMobilityActivationEnergy << std::endl;
    std::cout << "BurgersVectorBCC=" << mvBurgersVectorBCC << std::endl;
    std::cout << "BurgersVectorFCC=" << mvBurgersVectorFCC << std::endl;
    std::cout << "BurgersVectorHCP=" << mvBurgersVectorHCP << std::endl;
    std::cout << "TimeRangeToWriteOutput=" << mvTimeRangeToWriteOutput << std::endl;
    std::cout << "SoluteDiffusivityPreFactor=" << mvSoluteDiffusivityPreFactor << std::endl;
    std::cout << "SoluteDiffusivityActivationEnergy=" << mvSoluteDiffusivityActivationEnergy << std::endl;

    
}

void UserSettings::readFromFile(const std::string& inputFile) {
    std::ifstream in(inputFile);
    if (!in.is_open()) {
        std::cerr << "Unable to open file: " << inputFile << std::endl;
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            std::cout << "Key: " << key << ", Value: " << value << std::endl; // Debug print
            if (key == "TimeStep") {
                mvTimeStep = std::stod(value);
            } else if (key == "StartTemperature") {
                mvStartTemperature = std::stod(value);
            } else if (key == "EndTemperature") {
                mvEndTemperature = std::stod(value);
            } else if (key == "TimeTotal") {
                mvTimeTotal = std::stod(value);
            } else if (key == "Nx") {
                mvNx = std::stoi(value);
            } else if (key == "Ny") {
                mvNy = std::stoi(value);
            } else if (key == "Nz") {
                mvNz = std::stoi(value);
            } else if (key == "Dx") {
                mvDx = std::stod(value);
            } else if (key == "GridType") {
                mvGridType = std::stoi(value);
            } else if (key == "IsRexAndGG") {
                mvIsRexAndGG = std::stoi(value);
            } else if (key == "HasSoluteDiffusion") {
                mvHasSoluteDiffusion = std::stoi(value);
            } else if (key == "LowerMisorientationCutOff") {
                mvLowerMisorientationCutOff = std::stod(value);
            } else if (key == "HAGB") {
                mvHAGB = std::stod(value);
            } else if (key == "StartFileName") {
                mvStartFileName = value;
            } else if (key == "ThermodynamicsDataFilename"){
                mvThermodynamicsDataFilename = value;
            } else if (key == "OutputFolderPath") {
                mvOutputFolderPath = value;
            } else if (key == "IsTestRVE") {
                mvIsTestRVE = std::stoi(value);
            } else if (key == "TestVonNeumannInSquare") {
                mvTestVonNeumannInSquare = std::stoi(value);
            } else if (key == "AllowPeriodicBoundConditions") {
                mvAllowPeriodicBoundConditions = std::stoi(value);
            } else if (key == "MinimumCI") {
                mvMinimumCI = std::stod(value);
            } else if (key == "RexGGParameterC") {
                mvRexGGParameterC = std::stod(value);
            } else if (key == "ReadDislocationDensity") {
                mvReadDislocationDensity = std::stoi(value);
            } else if (key == "IncludeCSL19FastGrowth") {
                mvIncludeCSL19FastGrowth = std::stoi(value);
            } else if (key == "IsOnlyAnOrientationSubsetConsidered") {
                mvIsOnlyAnOrientationSubsetConsidered = std::stoi(value);
            } else if (key == "ConsiderConstantIronAtomsPerCell") {
                mvConsiderConstantIronAtomsPerCell = std::stoi(value);
            } else if (key == "SoluteSegregationDislocations") {
                mvSoluteSegregationDislocations = std::stoi(value);
            } else if (key == "KappaFactorForCsegInClusterArea") {
                mvKappaFactorForCsegInClusterArea = std::stod(value);
            } else if (key == "AllowDislocationsSegregationInAustenite") {
                mvAllowDislocationsSegregationInAustenite = std::stoi(value);
            } else if (key == "IncludeCarbonInterphasePartitioning") {
                mvIncludeCarbonInterphasePartitioning = std::stoi(value);
            } else if (key == "IsSoluteSegregationHappeningInDiffusionStep") {
                mvIsSoluteSegregationHappeningInDiffusionStep = std::stoi(value);
            } else if (key == "IsPartitioningHappeningInDiffusionStep") {
                mvIsPartitioningHappeningInDiffusionStep = std::stoi(value);
            } else if (key == "CarbonPartitioningFromInterpolatedSolutions") {
                mvCarbonPartitioningFromInterpolatedSolutions = std::stoi(value);
            } else if (key == "AllowInterfaceMovementDuringPartitioning") {
                mvAllowInterfaceMovementDuringPartitioning = std::stoi(value);
            } else if (key == "ConcentrationDependentDiffusivityInAustenite") {
                mvConcentrationDependentDiffusivityInAustenite = std::stoi(value);
            } else if (key == "XcToUseConstantDiffusivityAgrenInAustenite") {
                mvXcToUseConstantDiffusivityAgrenInAustenite = std::stod(value);
            } else if (key == "AverageCarbon") {
                mvAverageCarbon = std::stod(value);
            } else if (key == "GrainBoundaryEnergy") {
                mvGrainBoundaryEnergy = std::stod(value);
            } else if (key == "PhaseBoundaryEnergy") {
                mvPhaseBoundaryEnergy = std::stod(value);
            } else if (key == "GrainBoundaryMobilityProExponential") {
                mvGrainBoundaryMobilityProExponential = std::stod(value);
            } else if (key == "PhaseBoundaryMobilityProExponential") {
                mvPhaseBoundaryMobilityProExponential = std::stod(value);
            } else if (key == "GrainBoundaryMobilityActivationEnergy") {
                mvGrainBoundaryMobilityActivationEnergy = std::stod(value);
            } else if (key == "PhaseBoundaryMobilityActivationEnergy") {
                mvPhaseBoundaryMobilityActivationEnergy = std::stod(value);
            } else if (key == "BurgersVectorBCC") {
                mvBurgersVectorBCC = std::stod(value);
            } else if (key == "BurgersVectorFCC") {
                mvBurgersVectorFCC = std::stod(value);
            } else if (key == "BurgersVectorHCP") {
                mvBurgersVectorHCP = std::stod(value);
            } else if (key == "TimeRangeToWriteOutput") {
                mvTimeRangeToWriteOutput = std::stod(value);
            } else if (key == "SoluteDiffusivityPreFactor") {
                mvSoluteDiffusivityPreFactor = std::stod(value);
            } else if (key == "SoluteDiffusivityActivationEnergy") {
                mvSoluteDiffusivityActivationEnergy = std::stod(value);
            } else {
                std::cerr << "Unknown key: " << key << std::endl; // Debug for unknown keys
            }
        }
    }

    in.close();
}


void UserSettings::readStartFile(const std::string& startFileName) {
    std::ifstream in(startFileName);
    if (!in.is_open()) {
        std::cerr << "Unable to open file: " << startFileName << std::endl;
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        // Process start file content as needed
    }

    in.close();
}

void UserSettings::createOutputFolder(const std::string& folderPath) {
    try {
        if (!std::filesystem::exists(folderPath)) {
            std::filesystem::create_directories(folderPath); // Creates all non-existing directories in the path
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }
}

void UserSettings::clearOutputFolder(const std::string& folderPath) {
    std::filesystem::path path(folderPath);
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        std::filesystem::remove_all(entry.path());
    }
}
