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
#include "nodes.h"  // Include the Node definition



/**
 * @class Microstructure
 * @brief Manages the microstructure simulation for solid state transformations (phase transformations,
 * recrystallization, grain growth) and solute redistribution (partitioning, diffusion, trapping to defects).
 *
 * This class handles the initialization, simulation steps, and state management
 * for the microstructure in the OMicroN simulation program. It interacts with all other classes, setting
 * state variables and parameters, and retrieving information from instances.
 *
 * Additionally, it is the main class responsible for writing and exporting simulation outputs.
 */
class Microstructure
{
public:
    /**
     * @brief Constructs a Microstructure object with user settings.
     *
     * @param userSettings Reference to user settings.
     */
    Microstructure(const UserSettings &userSettings);

    /**
     * @brief Checks if the simulation is still running.
     *
     * @return True if the simulation is still running, otherwise false.
     */
    bool IsStillSimulating() const { return mvStillSimulating; }

    /**
     * @brief Prints the microstructure parameters.
     */
    void printMicrostructureParameters() const;

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
    const std::unordered_set<int> &GetInterfaceIndices() const { return interfaceIndices; }

    /**
     * @brief Gets the set of interphase indices.
     *
     * @return Const reference to the set of interphase indices.
     */
    const std::unordered_set<int> &GetInterphaseIndices() const { return interphaseIndices; }

    Eigen::SparseMatrix<float> K;  // Global stiffness matrix
    Eigen::VectorXf F;             // Force vector
    Eigen::VectorXf U; 
void initializeNodes();
    void solveFEM();
void assembleStiffnessMatrix();
float computeConditionNumber(const Eigen::SparseMatrix<float>& mat);

void initializeForceAndStiffnessMatrix();
void applyBoundaryConditions();
     float computeStrainStressForCell(int cellIndex);
     void applyUniaxialTension(float deformationRate);
bool isMatrixSymmetric(const Eigen::SparseMatrix<float>& K);
void applyDisplacementBoundaryConditions(
    Eigen::SparseMatrix<float>& K,
    Eigen::VectorXf& F,
    Eigen::VectorXf& U);
private:

    std::vector<Node> nodes; // Manage nodes for FEM


    // Data members containing information of newly transformed cells.
    std::vector<int> ConsumedCellsInStepBySamePhase;
    std::vector<int> ConsumedCellsInStepByDifferentPhase;

    // Internal parameters
    float mvStoreSoluteEvolutionTimeInterval = 0.5;
        float mvStoreSoluteEvolutionTimeStepInterval = 10;

    int mvEveryThatManySimStepsToCheckICells = 100; ///< Frequency of checking interface cells.
    int mvEveryThatManySimStepsExport = 1000; ///
    float mvStoreSoluteEvolutionTime;
    int mvStoreSoluteEvolutionTimeStep;
    int mvCheckICellsTime;                ///< Time to check interface cells.
    float mvTimePassedFromPreviousOutput; ///< Time that has passed from previous output time.
    int mvTimeStepsPassedFromPreviousOutput; ///< Time steps that have passed from previous output time.
    double mvMaxDiffusivityInTimeStep;
    bool mvStillSimulating = true;              ///< Flag to indicate if the simulation is still running.
    float mvTime;                               ///< Current simulation time.
    int mvSimulationStep;                       ///< Current simulation step.
    int NcellsForCarbonPart;                    ///< Number of neighbour cells for carbon partitioning.
    int NcellsForIntMig;                        ///< Number of neighbour cells for interface migration.
    float BetaFactorForInterfaceMigrationRates; ///< Beta factor for interface migration rates.
    float mvTimeStep;                           ///< Simulation time step.
    int cells_rx;                               ///< Number of recrystallized cells.
    float frx;                                  ///< Fraction of recrystallized cells.
    int grewFCC;                                ///< Number of FCC cells that grew.
    int grewBCC;                                ///< Number of BCC cells that grew.
    int grewHCP;                                ///< Number of HCP cells that grew.

    int NCellsFCC;                              ///< Number of FCC cells in RVE.
    int NCellsBCC;                              ///< Number of BCC cells in RVE.
    int NCellsHCP;                              ///< Number of HCP cells in RVE.
    int NCellsLiquid;                           ///< Number of liquid cells in RVE.

    // Data members for internally storing evolution data
    std::vector<float> TimeEvol;
    std::vector<float> Mean_DG_Chem_InFCC_Cells;
    std::vector<float> Mean_DG_Chem_InBCC_Cells;
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

    const UserSettings &mvUserSettings; ///< Reference to user settings.
    StateVars *mvpStateVars;            ///< Pointer to state variables.
    ThermChemKin *mvpTCK;               ///< Pointer to thermodynamic / chemical / kinetic data.
    Orientations mvOrientations;        ///< Object to manage orientations.

    std::unordered_set<int> interfaceIndices;  ///< Set of interface indices.
    std::unordered_set<int> interphaseIndices; ///< Set of interphase indices.

    /// @name Functions for initializing, checking, and setting state variables based on application, using UserSettings and informing classes Cells, ThermChemKin
    /// @{
    void InitInterfaceAndInterphaseCells(); ///< (Sub)grain boundary (based on local misorientation) and phase boundary cells (based on dissimilar phase neighbours) are determined for the first time and initialized.
    void CheckInterfaceAndInterphaseCells(); ///< (Sub)grain boundary (based on local misorientation) and phase boundary cells (based on dissimilar phase neighbours) are checked during the simulation for debugging purposes and to ensure validity when developing new models.
    void ValidateOrientationSettingsAndCheckNonIndexedPixels(); ///< Checkes the state of the initial microstructure that is read (e.g. the orientations are valid, how many cells are non indexed etc.)
    float SetStateVariablesBasedOnSimAndGiveInfoForMicrostructuralState(); ///< Sets values on the state variables of each cell depending on the desired simulation and application.
    void SetStateVariablesForSolidification();
    /// @}

    /// @name Functions for analyzing and saving in log file the solute redistribution average behavior
    /// @{
    void CarbonTemporalOverallQuantitiesForLaterAnalysis();
    /// @}

    
    /// @name Functions for analyzing phase fractions
    /// @{
    void CalculateNumbersOfCells();
    /// @}

    /// @name Functions for writing output files containing pre-defined state variables based on application
    /// @{
    void ReadMicrostructureFile(const std::string &startFileName); ///< Determines the format of microstructure file 
    void ReadMicrostructureTextFile(const std::string &startFileName); ///< Reads microstructure file. The file should contain one header line with the number of pixels/voxels. After this the columns are:  x (um), y (um),  z (um), latticeId (0 for FCC, 1 for BCC, 2 for HCP), disl.density (1/um2),  euler phi1 (degrees), euler Phi (degrees), euler phi2 (degrees), confidence index,  rex or not (0 or 1)
    void ReadMicrostructureHDF5File(const std::string &startFileName); ///< Reads microstructure file. The file should contain one header line with the number of pixels/voxels. After this the columns are:  x (um), y (um),  z (um), latticeId (0 for FCC, 1 for BCC, 2 for HCP), disl.density (1/um2),  euler phi1 (degrees), euler Phi (degrees), euler phi2 (degrees), confidence index,  rex or not (0 or 1)
    void MakeMicrostructure();
    void MakeMicrostructureTestDeformation();
    void MakeMicrostructureRolled();
   void WriteNodesPositions(const std::string &fullPath);
void WriteMicrostructureQP(const std::string &fullPath);
    void WriteMicrostructureRX(const std::string &writeFileName);
void ExportColorCodedTIFF(const std::string &fullPath, int MapType);
void exportNodePositionsAsTiff(const std::string& filename);
void exportDeformedTiff(const std::string& filename, const std::vector<float>& data, int nx, int ny, int nz);
    void ExportColorbar(const std::string &fullPath);
    void GetColorFromValue(float value, float minValue, float maxValue, uint8_t &r, uint8_t &g, uint8_t &b);
    float GetCellsMisorientationDarkeningFactor(int index);
std::vector<float> GetNDCrystalDirection(int index, int crystalStructure);
void GetColorCoding(int width, std::vector<uint8_t> &raster, int MapType);
std::vector<std::vector<std::vector<float>>> GetOrientationMatrices();
std::vector<int> GetCrystalStructures();
void GetIPFColor(const std::vector<float>& crystalDirection, int crystalStructure, uint8_t &r, uint8_t &g, uint8_t &b);
std::vector<float> ReduceToFundamentalTriangle(const std::vector<float>& direction);

// void GetIPFColor(const std::vector<float>& crystalDirection, int crystalStructure, uint8_t &r, uint8_t &g, uint8_t &b);
std::vector<float> Normalize(const std::vector<float>& vec);


    /// @}

    /// @name Functions for connecting simulation conditions (application) to the CA settings per simulation step
    /// @{
    void SetTemperatureAndInformInThermo();
    /// @}

    /**
     * @name Functions used for interface motion (e.g., Recrystallization, grain growth, phase transformations)
     * @{
     */
    /**
     * @name Common functions
     * @{
     */
    void CalculateDirectionOfForceAtThisBoundary(int index1, int index2, double *ForceUnitVectorToAxisX, double *ForceUnitVectorToAxisY, double *ForceUnitVectorToAxisZ);
    double GetAreaPerVolumeOfBoundary(int index1, int index2);
    void SwitchCell(int index, int indexThatGrew, bool IsRexOrGG);
    void UpdateNumbersOfCells(int newPhase, int oldPhase);
    /** @} */

    /**
     * @name Auxiliary functions
     * @{
     */
    bool IsCellIndeedInterphase(int index);
    bool IsCellIndeedInterface(int index);
    
    /** @} */

    /**
     * @name Recrystallization and grain growth
     * @{
     */
    float CalcReorientationRatesForFullFieldRexOrGG(int latticeId);
    int ReOrientCellsForFullFieldRexOrGG(int latticeId);
    void UpdateInterfaceCells(void);
    void UpdateThisInterfaceCellAndNeighbours(int indexThatGotConsumed);
    /** @} */

    /**
     * @name Phase transformations
     * @{
     */
    float CalcConsumptionRatesInterfaceMigrationDuringPartitioning(void);
    int TransformCellsFromInterfaceMigration(void);
    void UpdateInterphaseCells(void);
    void UpdateThisInterphaseCellAndNeighbours(int indexThatGotConsumed);
    double CalcConsumptionRatesDuringSolidification();
    int TransformCellsFromSolidification();
    /** @} */



double GetAnisotropicBoundaryEnergyFactor(int indexGrowing,int indexConsumed);
double ComputeAngle(const Eigen::Vector3d& vec1, const Eigen::Vector3d& vec2);
Eigen::Vector3d TransformToGlobal(const Eigen::Vector3d& crystalDirection, const Eigen::Matrix3d& rotationMatrix);
Eigen::Vector3d FindClosestPlaneNormal(const Eigen::Vector3d& dirCrystal, const std::vector<Eigen::Vector3d>& crystalPlanes);
Eigen::Vector3d CalculatePlaneNormal(const Eigen::Vector3d& dirCartesian, const Eigen::Matrix3d& rotationMatrix);
    /** @} */



    /**
     * @name Functions used for solute redistribution (e.g., carbon partitioning / trapping to defects)
     * @{
     */
    double CalculateDiffusivityMatrixAndReturnMaxDiffusivityInStep(); /**< Calculate diffusivity matrix and return the maximum diffusivity in the current step. */
        void ChargeTheRVE();

    /** @} */
 /**
     * @name Functions used for solute redistribution during solidification (e.g., liquid partitioning / trapping to defects)
     * @{
     */
void CalculateDiffusivityMatrixForImpurities(double maxDiffusivityOfElement);


    /**
     * @name Auxiliary functions for repetitive transitions within CA state variables
     * @{
     */

    /**
     * @brief Gets the molar volume of the cell with the specified index.
     *
     * @param index The index of the cell.
     * @return The molar volume of the cell.
     */
    double GetMolarVolumeOfCell(int index);

    /**
     * @brief Gets the molar volume of the cell with the lattice id.
     *
     * @param latticeId The lattice id cell.
     * @return The molar volume of the cell.
     */
    double GetMolarVolumeOfLattice(int latticeId);

    /**
     * @brief Gets the mobility for the given misorientation in the specified lattice.
     *
     * @param latticeId The ID of the lattice.
     * @param mis The misorientation value.
     * @return The mobility associated with the given misorientation.
     */
    double GetMobilityForThisMisorientation(int latticeId, double mis);

    /**
     * @brief Computes the misorientation between two cells and optionally checks for CSL relationship.
     *
     * @param index1 The index of the first cell.
     * @param index2 The index of the second cell.
     * @param MaxHagb Flag indicating whether to consider maximum Hagb.
     * @param CheckAlsoForSameOrientation Flag indicating whether to check for the same orientation.
     * @param CheckIfCSL Flag indicating whether to check for CSL (Coincidence Site Lattice) relationship.
     * @param CSLRelationship Pointer to an integer to store the CSL relationship, if applicable.
     * @return The computed misorientation between the two cells.
     */
    double MisorientationBetweenCells(int index1, int index2, bool MaxHagb, bool CheckAlsoForSameOrientation, bool CheckIfCSL, int *CSLRelationship);

    /**
     * @brief Gets the boundary energy for the specified misorientation.
     *
     * @param mis The misorientation value.
     * @return The boundary energy associated with the given misorientation.
     */
    double GetBoundaryEnergyForThisMisorientation(double mis);
    /** @} */

    /**
     * @name Auxiliary functions for computational operations
     * @{
     */

    /**
     * @brief Checks if the given float number is finite.
     *
     * @param numberToCheck The float number to check.
     * @return True if the number is finite, otherwise false.
     */
    bool IsTheNumberFinite(float numberToCheck);

    /**
     * @brief Checks if the given double number is finite.
     *
     * @param numberToCheck The double number to check.
     * @return True if the number is finite, otherwise false.
     */
    bool IsTheNumberFinite(double numberToCheck);

    /**
     * @brief Finds the index of the given ID in the vector, if it exists.
     *
     * @param vectorToCheck The vector to search in.
     * @param id The ID to find.
     * @return An optional containing the index if found, otherwise std::nullopt.
     */
    std::optional<int> findIndex(const std::vector<int> &vectorToCheck, int id);
        /**
     * @brief Generates "nNumbers" numbers (integers) from a range [nMin,nMax]
     *
     * @param nIntegers the number of integers to generate
     * @param nMax the minimum number 
     * @param nMax the maximum number
     * @return the list of randomly generated integers
     */
    std::vector<int> generateRandomIntegers(int nIntegers, int nMin, int nMax);

    /** @} */
};

#endif // Microstructure_H
