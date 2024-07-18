/* microstructure.h
   OMicroN (optimising microstructures numerically) simulation program
   header file containing Microstructure class definitions and implementation
*/

#ifndef Microstructure_H
#define Microstructure_H

#include <vector>
#include <string>
#include <unordered_set>
#include <optional> // Add this line for std::optional
#include "settings.h"
#include "stateVars.h"
#include "orientation.h"
#include "ThermoChemKinetics.h"
#include <loguru.hpp>

/**
 * @class Microstructure
 * @brief Manages the microstructure simulation for solid state transformations (phase transformations,
 * recrystallization, grain growth) and solute redistribution (partitioning, diffusion, trapping to defects).
 *
 * This class handles the initialization, simulation steps, and state management
 * for the microstructure in the OMicroN simulation program.
 * 
 * It is the only (and most main) class that communicates with all the others, setting (with function) state variables and parameters to other classes,
 * and getting information from all instances.
 * 
 * Finally, it is the main class in terms of writing and exporting simulation outputs.
 */


class Microstructure {
public:
    /**
     * @brief Constructs a Microstructure object with user settings.
     * 
     * @param userSettings Reference to user settings.
     */
    Microstructure(const UserSettings& userSettings);

    /**
     * @brief Checks if the simulation is still running.
     * 
     * @return True if the simulation is still running, otherwise false.
     */
    bool IsStillSimulating() { return mvStillSimulating; }

    /**
     * @brief Prints the microstructure parameters.
     */
    void printMicrostructureParameters();

    /**
     * @brief Executes a simulation step.
     */
    void SimulationStep();

    /**
     * @brief Adds an index to the set of interface indices.
     * 
     * @param index Index to add.
     */
    void AddInterfaceIndex(int index);

    /**
     * @brief Removes an index from the set of interface indices.
     * 
     * @param index Index to remove.
     */
    void RemoveInterfaceIndex(int index);

    /**
     * @brief Adds an index to the set of interphase indices.
     * 
     * @param index Index to add.
     */
    void AddInterphaseIndex(int index);

    /**
     * @brief Removes an index from the set of interphase indices.
     * 
     * @param index Index to remove.
     */
    void RemoveInterphaseIndex(int index);

    /**
     * @brief Gets the set of interface indices.
     * 
     * @return Const reference to the set of interface indices.
     */
    const std::unordered_set<int>& GetInterfaceIndices() const { return interfaceIndices; }

    /**
     * @brief Gets the set of interphase indices.
     * 
     * @return Const reference to the set of interphase indices.
     */
    const std::unordered_set<int>& GetInterphaseIndices() const { return interphaseIndices; }



private:
   // Data members containing information of newly transformed cells. 
    std::vector<int> ConsumedCellsInStepBySamePhase;
        std::vector<int> ConsumedCellsInStepByDifferentPhase;

           //Internal parameters 
float mvStoreSoluteEvolutionTimeInterval = 0.1;
    int mvEveryThatManySimStepsToCheckICells = 10; ///< Frequency of checking interface cells.
    float mvStoreSoluteEvolutionTime;
    int mvCheckICellsTime; ///< Time to check interface cells.
    float mvTimePassedFromPreviousOutput; ///< Time that been passed from previous output time
double mvMaxDiffusivityInTimeStep;
    bool mvStillSimulating; ///< Flag to indicate if the simulation is still running.
    float mvTime; ///< Current simulation time.
    int mvSimulationStep; ///< Current simulation step.
    int NcellsForCarbonPart; ///< Number of cells for carbon partitioning.
    int NcellsForIntMig; ///< Number of cells for interface migration.
    float BetaFactorForInterfaceMigrationRates; ///< Beta factor for interface migration rates.
    float mvTimeStep; ///< Simulation time step.
    int cells_rx; ///< Number of recrystallized cells.
    float frx; ///< Fraction of recrystallized cells.
    int grewFCC; ///< Number of FCC cells that grew.
    int grewBCC; ///< Number of BCC cells that grew.


   // Data members for internally storing evolution data
    std::vector<float> TimeEvol;
    std::vector<float> Mean_DG_Chem_InFCCCells;
    std::vector<float> Mean_DG_SE_InFCCCells;
    std::vector<float> Mean_DG_Chem_InBCCCells;
    std::vector<float> Mean_DG_SE_InBCCCells;
    std::vector<float> BCC_Growth;
    std::vector<float> FCC_Growth;
    std::vector<float> FCC_XC;
    std::vector<float> BCC_XC;
    std::vector<float> FCC_XC_Dispersion;
    std::vector<float> BCC_XC_Dispersion;
    std::vector<float> FCC_IcellsXC;
    std::vector<float> BCC_IcellsXC;
    std::vector<float> FCC_XCTrappedKappaFactor;
    std::vector<float> BCC_XCTrappedKappaFactor;
    std::vector<float> FCC_TrappedXC;
    std::vector<float> BCC_TrappedXC;

    const UserSettings& mvUserSettings; ///< Reference to user settings.
    StateVars* mvpStateVars; ///< Pointer to state variables.
    ThermChemKin* mvpTCK; ///< Pointer to thermodynamic / chemical / kinetic data.
    Orientations mvOrientations; ///< Object to manage orientations.

    std::unordered_set<int> interfaceIndices; ///< Set of interface indices.
    std::unordered_set<int> interphaseIndices; ///< Set of interphase indices.

    /// @name Functions for initializing, checking, and setting state variables based on application, by using UserSettings and informing also appropriately classes Cells, ThermChemKin
    /// @{
    void InitInterfaceAndInterphaseCells();
    void CheckInterfaceAndInterphaseCells();
    void ValidateOrientationSettingsAndCheckNonIndexedPixels();
    float SetStateVariablesBasedOnSimAndGiveInfoForMicrostructuralState();
    /// @}

   /// @name Functions for analyzing and saving in log file the solute redistribution average behaviour 
    /// @{
    void CarbonTemporalOverallQuantitiesForLaterAnalysis();
    /// @}


    /// @name Functions for writing output files containing pre-defined state variables based on application
    /// @{
    void ReadMicrostructureFile(const std::string& startFileName);
    void WriteMicrostructureRX(const std::string& writeFileName);
    void GetColorCoding(int width, std::vector<uint8_t>& raster, int MapType);
    void ExportColorCodedTIFF(const std::string& fullPath, int MapType);
    void ExportColorbar(const std::string& fullPath);

    /// @}

    /// @name Functions for connecting simulation conditions (application) to the CA settings per simulation step
    /// @{
    void SetTemperatureAndInformInThermo();
    /// @}

    /// @name Functions used for interface motion (e.g., Recrystallization, grain growth, phase transformations)
    /// @{
    void CalculateDirectionOfForceAtThisBoundary(int index1, int index2, double* ForceUnitVectorToAxisX, double* ForceUnitVectorToAxisY, double* ForceUnitVectorToAxisZ);
    double GetAreaPerVolumeOfBoundary(int index1, int index2);
    float CalcReorientationRatesForFullFieldRexOrGG(int latticeId);
    int ReOrientCellsForFullFieldRexOrGG(int latticeId);
    void SwitchCell(int index, int indexThatGrew, bool IsRexOrGG);
    void UpdateInterfaceCells();
    void UpdateInterphaseCells();
    bool IsCellIndeedInterface(int index);
    bool IsCellIndeedInterphase(int index);
        void UpdateThisInterfaceCellAndNeighbours(int indexThatGotConsumed);
    void UpdateThisInterphaseCellAndNeighbours(int indexThatGotConsumed);

    /// @}

  /// @name Functions used for solute redistribution (e.g. carbon partitioning / trapping to defects)
    /// @{

    double CalculateDiffusivityMatrixAndReturnMaxDiffusivityInStep();


    /// @name Auxiliary functions for repetitive transitions within CA state variables
    /// @{
    double GetMolarVolumeOfLattice(int latticeId);
    double GetMolarVolumeOfCell(int index);
    double GetMobilityForThisMisorientation(int latticeId, double mis);
    double MisorientationBetweenCells(int index1, int index2, bool MaxHagb, bool CheckAlsoForSameOrientation, bool CheckIfCSL, int* CSLRelationship);
    double GetBoundaryEnergyForThisMisorientation(double mis);
    /// @}

    /// @name Auxiliary functions for computational operations
    /// @{
    bool IsTheNumberFinite(float numberToCheck);
    bool IsTheNumberFinite(double numberToCheck);
    std::optional<int> findIndex(const std::vector<int>& vectorToCheck, int id);
    /// @}
};

#endif // Microstructure_H
