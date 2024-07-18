#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem> // For creating directories
#include <stdexcept> // For std::runtime_error
#include <cstdio> // For std::remove

class UserSettings {
public:
    UserSettings(const std::string& inputFile);

    void setDefaultValues();
    void writeToFile(const std::string& outputFile);
    void printParameters();

    // Public getters for parameters
    std::string getStartFileName() const {
        return mvStartFileName;
    }

    std::string getOutputFolderPath() const {
        return mvOutputFolderPath;
    }

    std::string getThermodynamicDataFile() const {
        return mvThermodynamicsDataFilename;
    }

    double mvTimeStep;
    double mvStartTemperature;
    double mvEndTemperature;
    double mvTimeTotal;
    int mvNx, mvNy, mvNz;
    double mvDx;
    int mvGridType;
    int mvHasSoluteDiffusion;
    int mvIsRexAndGG;
    float mvLowerMisorientationCutOff;
    float mvHAGB;
    std::string mvStartFileName;
    std::string mvOutputFolderPath; 
    std::string mvThermodynamicsDataFilename; 

    int mvIsTestRVE;
    int mvTestVonNeumannInSquare;
    int mvAllowPeriodicBoundConditions;
    float mvMinimumCI;
    float mvRexGGParameterC;
    int mvReadDislocationDensity;
    int mvIncludeCSL19FastGrowth;
    int mvIsOnlyAnOrientationSubsetConsidered;
    int mvConsiderConstantIronAtomsPerCell;
    int mvSoluteSegregationDislocations;
    double mvKappaFactorForCsegInClusterArea;
    int mvAllowDislocationsSegregationInAustenite;
    int mvIncludeCarbonInterphasePartitioning;
    int mvIsSoluteSegregationHappeningInDiffusionStep;
    int mvIsPartitioningHappeningInDiffusionStep;
    int mvCarbonPartitioningFromInterpolatedSolutions;
    int mvAllowInterfaceMovementDuringPartitioning;
    int mvConcentrationDependentDiffusivityInAustenite;
    double mvXcToUseConstantDiffusivityAgrenInAustenite;
    double mvAverageCarbon;
    double mvGrainBoundaryEnergy;
    double mvPhaseBoundaryEnergy;
    double mvGrainBoundaryMobilityProExponential; ///< Grain boundary mobility pre-factor (mol m/J s) e.g. M0AA=1.95
    double mvPhaseBoundaryMobilityProExponential; ///< Phase boundary mobility pre-factor (mol m/J s) e.g. M0AA=1.95
    double mvGrainBoundaryMobilityActivationEnergy;   ///<  Grain boundary mobility activation energy (J/mol) e.g.  QgAA=140000.0
    double mvPhaseBoundaryMobilityActivationEnergy;  ///<  Phase boundary mobility activation energy (J/mol) e.g.  QgAA=140000.0
    double mvBurgersVectorBCC;
    double mvBurgersVectorFCC;
    double mvBurgersVectorHCP;
    double mvTimeRangeToWriteOutput;
    double mvSoluteDiffusivityPreFactor;
    double mvSoluteDiffusivityActivationEnergy;

private:
    void readFromFile(const std::string& inputFile);
    void readStartFile(const std::string& startFileName);
    void createOutputFolder(const std::string& folderPath);
    void clearOutputFolder(const std::string& folderPath);
};

#endif // SETTINGS_H
