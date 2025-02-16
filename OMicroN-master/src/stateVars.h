
#ifndef StateVars_H
#define StateVars_H

#include <cstdlib>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseCore>
#include <set>
#include <vector>
#include <iostream>
#include <loguru.hpp>
#include <ThermoChemKinetics.h>
#include <eigen3/Eigen/Dense> // operations between pointers

class StateVars;

/// [Periodic] boundaries of the grid
enum BOUNDARIES
{
    X_LOWER = 1 << 0, ///< X-axis lower boundary
    X_UPPER = 1 << 1, ///< X-axis upper boundary
    Y_LOWER = 1 << 2, ///< Y-axis lower boundary
    Y_UPPER = 1 << 3, ///< Y-axis upper boundary
    Z_LOWER = 1 << 4, ///< Z-axis lower boundary
    Z_UPPER = 1 << 5, ///< Z-axis upper boundary
};




class StateVars {
private:
    /// @name Cellular automata geometry data
    /// @{
    float MinCFree = 1.e-15;
    bool mvKeepConstantIronAtomsPerCell; ///<Internal flag to keep in mind that SoluteDiffusionStep (and partitioning) solve XC equations with XC [carbon atoms/iron atoms]. In this case, while other functions (and classes thermo, Microstructure etc.) have carbon in at. frac., we convert before solving the diffusion step (and then go back to at.frac.)
    bool mvIs2DSimulation;                                ///< Internal flag for checking if we are running a 2D simulation (i.e. Nz==1)

    bool mvIncludeSoluteDiffusion;              ///<Flag to inform the constructor whether we have solute diffusion in simulation, in order to initalize relevant pointers
    bool mvIncludeRexAndGG;                     ///<Flag to inform the constructor whether we have recrystallization and/or grain growth in simulation, in order to initalize relevant pointers
    bool mvAllowIntMigration;                   ///<Flag to inform the constructor whether we have phase transformations in simulation, in order to initalize relevant pointers
    int mvN;                                    ///< Total number of material elements (#mvN = #mvNx * #mvNy * #mvNz)
    int mvNx;                                   ///< Number of material elements along the X direction
    int mvNy;                                   ///< Number of material elements along the Y direction
    int mvNz;                                   ///< Number of material elements along the Z direction
    int mvNxy;                                  ///< #mvNx * #mvNy
    float mvDx;                                 ///< Grid spacing (material element dimension) in m
    float mvDx2;                                ///< #mvDx squared
    float mvDx3;                                ///< #mvDx cubed
    int mvGridType;                ///< GridType: 0 for Regular 1 for hexagonal
    float mvMaxBoundaryAreaPerVol; ///< Ratio of one boundary's area divided by material elements's volume: It depends on the type grid



    int mvNNearest;          ///< Number of nearest material elements in the Moore neighbourhood (2D: 4, 3D: 6)
    int mvNNext;             ///< Number of second nearest material elements in the Moore neighbourhood (2D: 4, 3D: 12)
    int mvNNextNext;         ///< Number of third nearest material elements in the Moore neighbourhood (2D: 0, 3D: 8)
    int mvNAll;              ///< Total number of material elements in the Moore Neighbourhood neighbours (2D: 8, 3D: 26)
    int mvNNeighboursMatter; ///< Total number of material elements that we care about depending on grid and application, e.g. if fulldiffusion we care about only nearest, if full-field grain growth we consider only acceptable grid settings (because of GB energy anisotropy)

    int mvpAllNeighbours[26]; ///< Relative indices (1D array) of all (first, second, and third) neighbours
float *mvpVonMisesStress = nullptr;

    /// @name The pointers to most basid (and always initialized during simulation) state variables of material elements information
    /// @{
    int *mvpLatticeId = nullptr;           ///< LatticeType type of the cell
    unsigned char *mvpBoundaryCell = nullptr; ///< Flag that determines if a cell is on a system (periodic) boundary
    std::vector<int> mvpOriId;  ///< Orientation ID of the cell
    int *mvpIsInterface = nullptr;           ///< Flag that determines if a cell is a (sub)grain boundary cell
    int *mvpIsInterphase = nullptr;        ///< Flag that determines if a cell is a phase boundary cell
    /// @}

    /// @name Local composition data
    /// @{
// NECESSARY FUNCTIONS
    Eigen::VectorXd* mvpXC  = nullptr;        ///< Carbon concentration per cell
    Eigen::VectorXd* mvpXCTrapped = nullptr;        ///< Carbon concentration trapped per cell
// NOT NECESSARY (AND NOT ADVISABLE TO USE) FUNCTIONS - partitioning and trapping should ideally take place inside the diffusion step solver (they should be part of Ax=b system)
    Eigen::VectorXd *mvpXCEqAtInterphase = nullptr; 
   Eigen::VectorXd *mvpXCTrappedEqWithNbr = nullptr; 
    /// @}


    Eigen::VectorXd* mvpXCu  = nullptr;        ///< Copper concentration per cell
    Eigen::VectorXd* mvpXSn  = nullptr;        ///< Tin concentration per cell

    /// @name Carbon diffusion due to interface partitioning. Also carbon segragation
    double *mvpKappaFactorForSoluteTrappedDefects = nullptr;
    double *mvpDiffusivitySolute_Ratio = nullptr;    ///< ratio of actual concentration depenent diffusivity (Agren) to the max carbon diffusivity. Written as ratio to directly go in the matrix A when solving Ax=b.

    int *mvpOneOfTheNeighboursGrowingInto = nullptr; ///< This is the material element's (or one of the material elements) growing to the material element with the orientation of grain mvpIdToGive. It is useful only for rexInitiation to take info from cell or processing
     float *mvpConsumptionRate = nullptr;           ///< A cell's re-orientation rate into the orientation of one neighbours
    float *mvpConsumedFraction = nullptr;          ///< A cell's re-orientation fraction (divided by boundary area) into the orientation of one of its neighbours
    /// @name Grid related data

    /// @name GG and/or static recrystallisation data
   
    int *mvpCSLRelationshipNow = nullptr;

    /// @name Input mmicrostructure related data i.e. EBSD
    float *mvpCI = nullptr;           ///< Cell's CI: Only useful when EBSD maps are used (for diminishing processing artifacts)
    int *mvpHasNeighNonInd = nullptr; ///< If cell has neighbour with non-reliable CI: Only useful when EBSD maps are used (for diminishing processing artifacts)
    float *mvpKAM = nullptr;  
    float *mvpRho = nullptr;  

    /// @name Info for later processing if REX
    int *mvpRX = nullptr;            ///< Is the cell RX (passed by any boundary - lagb or hagb)
    float *mvpMaxGbPassed = nullptr; ///< Max misorientation angle that boundary had when surpassing this cell (processing purposes to analyze extent of recovery, rex)


    /// @name Diffusion calculation step data
    /// @{
    float mvAlphaDt;                     ///< Factor used for solving carbon diffusion time step, mvAlphaDt = dt * maxDiffusivityInTimeStep / (mvDx * mvDx);

    /// @}

   /// @name Initialization functions
    /// @{
    void InitializeBasicStateVariables(void); 
    void InitializeSoluteRedistributionAndDefectsVariables(void);
    void InitializeStateVariablesRelatedToInterfaceMotion(void);
    void InitBoundaryCells(void);
    void allocateCarbonVectors(bool initXCEqInt, bool initXCEqDef);
    void InitDiffusivitiesVector(void);
    void InitEBSDRelatedStuff();
    void InitStateVariablesRelatedToRexAndGG();
    void InitImpurities();

    /// @}

    /// @name "Check" functions
    /// @{
    void validateAllNeighbourCells(void);
   /// @}

  
    /// @name Diffusion calculation functions
    /// @{
/** @brief Numerical solution of diffusion step based on conjugate gradient method
 * @param b the b=Ax (which is equal to current concentratons initially)
 * @param TCK_p pointer to ThermChemKin instance
 * @param maxit maximum allowed iterations
 * @param tolerance tolerance for convergence
 * @param initial_guess the first x values guessed (set here as current solute concentrations)
 * @param IsPartitioningHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 * @param IsSoluteSegregationHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 */
Eigen::VectorXd NumericalSolverCG(const Eigen::VectorXd &b, ThermChemKin *TCK_p, int maxit, double tolerance, const Eigen::VectorXd &initial_guess, bool IsPartitioningHappeningHere, bool IsSoluteSegregationHappeningHere,bool UseChemPot, int highSolPhase);
  
  /** @brief Performs one computation of the system Ax=b for Traka's diffusion / partitioning / trapping to defects
 * Traka, 2024, Acta Materialia.
 *
 * @param x the solute concentrations vector which remains here as such
 * @param TCK_p pointer to ThermChemKin instance
 * @param b the solute concentrations vector which here changes according to the reduction of Ax
 * @param p_dotAp dot product of new b
 * @param IsPartitioningHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 * @param IsSoluteSegregationHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 */

    void SystemATimesX(const Eigen::VectorXd &x, ThermChemKin *th_p, Eigen::VectorXd &b, double &p_dotAp, bool IsPartitioningHappeningHere, bool IsSoluteSegregationHappeningHere, int HighSolPhase);

    ThermChemKin* mvpTCK; ///< Pointer to thermodynamic data.

bool IsTheNumberFinite(double numberToCheck);
   
     /** @brief Performs one computation of the system Ax=b for Traka's diffusion / partitioning / trapping to defects
 * Traka, 2024, Acta Materialia. Here it is for solidification.
 *
 * @param x the solute concentrations vector which remains here as such
 * @param TCK_p pointer to ThermChemKin instance
 * @param b the solute concentrations vector which here changes according to the reduction of Ax
 * @param p_dotAp dot product of new b
 * @param IsPartitioningHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 * @param IsSoluteSegregationHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 */
void SystemATimesXSolidification(const Eigen::VectorXd &x, ThermChemKin *TCK_p, Eigen::VectorXd &b, double &p_dotAp, bool IsPartitioningHappeningHere, bool IsSoluteSegregationHappeningHere);



public:

    /// @name Constructors and destructor
    /// @{

    /// @brief Default constructor
    StateVars(void)
        : mvN(-1)
    {
    }
    StateVars(int nx, int ny, int nz, float dx, int grid, bool IncludeSoluteDiffusion, bool IncludeRexAndGG, bool AllowIntMigration, bool initXCEqInt, bool initXCEqDef, bool isSolidification, bool isDeformation);
    ~StateVars(void);
    /// @}

   
void SetVonMisesStress(int index, float stress){mvpVonMisesStress[index] = stress;};
float GetVonMisesStress(int index){return mvpVonMisesStress[index];};


     /// @name Functions for setting state variables of cells
    /// @{

    void SetEBSDRelatedStuffForRexAndGG(int index,float CI, int RX, float rho);
    void SetStateVariablesRelatedToRexAndGG(int index,float ReRate, int ReFraction, float MaxAnglePassed);
  
    void SetCellAsInterface(int index){ mvpIsInterface[index] = 1;}
    void SetCellAsInterphase(int index){ mvpIsInterphase[index] = 1;}
    void SetCellAsNonInterface(int index){ mvpIsInterface[index] = 0;}
    void SetCellAsNonInterphase(int index){ mvpIsInterphase[index] = 0;}

    int IsCellInterface(int index){return mvpIsInterface[index];}
    int IsCellInterphase(int index){return mvpIsInterphase[index];}

    void SetOriIdOfCell(int index, int OriId) { mvpOriId[index] = OriId; }
    void SetIfCellHasNonIndNeighbours(int index, int HasNeighNonInd) { mvpHasNeighNonInd[index] = HasNeighNonInd; }
    void SetLatticeIdOfCell(int index, int latticeValue){mvpLatticeId[index]=latticeValue;}
///@}
    /// @name Function for deleting KAM - it is called after initialing all pointers with simulation-necessary state variables (KAM is not useful after setting dislocation density (if not read)).
    /// @{
void eraseKAM(){   
     if (mvpKAM)
    {
        delete[] mvpKAM;
        mvpKAM = nullptr;
    }
    }
    /// @name Functions for retrieving state variables of cells as they are stored in CA (orientation id, lattice id)
    /// @{
    int GetLatticeIdOfCell(int index) { return mvpLatticeId[index]; }
    int GetOriIdOfCell(int index) { return mvpOriId[index]; }

    /// @name Functions for retrieving lattice name based on lattice Id stored in CA
    /// @{
bool IsCellFCC(int index)
{
 return (mvpLatticeId[index] == 0);
     }

bool IsCellBCC(int index)
{
 return (mvpLatticeId[index] == 1);
     }

bool IsCellHCP(int index)
{
 return (mvpLatticeId[index] == 2);
     }

bool IsCellLiquid(int index)
{
 return (mvpLatticeId[index] == 3);
     }


  
    /// @name Functions for retrieving lattice name based on lattice Id stored in CA
    /// @{
bool IsLatticeIdFCC(int latticeId)
{
 return (latticeId == 0);
     }
bool IsLatticeIdBCC(int latticeId)
{
 return (latticeId == 1);
     }
     bool IsLatticeIdHCP(int latticeId)
{
 return (latticeId == 2);
     }
     bool IsLatticeIdLiquid(int latticeId)
{
 return (latticeId == 3);
     }


    /// @name Functions for retrieving geometry data
    /// @{
    std::vector<int> GetAllNeighbourCells(int index) const;
    void GetNeighbourCellOffsets(int index, int *all_p) const;


    int GetGridType(void) const { return mvGridType; }




    // i, j, k -> index
    int IJKToIndex(int i, int j, int k) const;
    /// Alias of IJK2Index
    int GetIndexFromIJK(int i, int j, int k) const { return IJKToIndex(i, j, k); };

    // index -> i, j, k (grid positions, dimensionless)
    void IndexToIJK(int index, int *i, int *j, int *k) const;
    void IndexToIJK(int index, float *i, float *j, float *k) const;
    std::array<int, 3> GetIJKFromIndex(int index) const;

    // index -> x, y, z (grid coordinates, in m)
    float GetXFromIndex(int index) const;
    float GetYFromIndex(int index) const;
    float GetZFromIndex(int index) const;
    std::array<float, 3> GetXYZFromIndex(float index) const;

    /// @brief Returns the distance between two cells (in m)
    float GetDistanceBetweenCells(int indexA, int indexB) const;


    /** @brief Returns (optionally used) local carbon concentration
     *
     * @param index Cell index
     *
     * @return Local carbon concentration of cell in at. fraction is available,
     * otherwise, returns -1
     */

    double GetXCOfCell(int index) const { return (mvpXC != nullptr) ? (*mvpXC)[index] : -1; }
    double GetKappaFactorForCTrappedInCell(int index) const { return (mvpKappaFactorForSoluteTrappedDefects != nullptr) ? mvpKappaFactorForSoluteTrappedDefects[index] : 0; }
    double GetXCTrappedInCell(int index) const { return (mvpXCTrapped != nullptr) ? (*mvpXCTrapped)[index] : 0; }
    void ConvertTotCellConcentrationsInCarbonPerIron() const;
    void ConvertTotCellConcentrationsInAtFraction() const;

    /** @brief Returns the value of mvpBoundaryCell. This value is associated to the position
     * of the cell on the grid. If it belongs to the boundary, it assumes a value > 0.
     *
     * @param index Cell index
     *
     * @return The value of mvpBoundaryCell, which is > 0 if the cell belongs to the periodic
     * boundary. 0 otherwise.
     */
    unsigned char GetBoundaryCell(int index) const { return mvpBoundaryCell[index]; }

    /** @brief Returns flag that determines if a cell is on a system (periodic) boundary
     *
     * @param index Cell index
     *
     * @return A boolean: true if cell belongs to periodic boundary, false otherwise.
     */
    bool IsCellOnBoundary(int index) const { return GetBoundaryCell(index) > 0; }


    /// @name Diffusion calculation
    /// @{
    
/** @brief Solves one diffusion step
 *
 * @param TCK_p pointer to ThermChemKin instance
 * @param dt predetermined time step
 * @param AllowSoluteSegregation bool on whether we have solute trapping in general during simulation
 * @param maxDiffusivityInTimeStep the maximum value of diffusivity between neighbour cells calculated in microstructure class
 * @param IsPartitioningHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 * @param IsSoluteSegregationHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 */

void SoluteDiffusionStep(ThermChemKin *TCK_p, double dt, bool AllowSoluteSegregation, double maxDiffusivityInTimeStep, bool IsPartitioningHappeningHere, bool IsSoluteSegregationHappeningHere, int ElementNumber, int LatticeOfHighSol);
    void PutBackCarbonTrappedAndCalculateNewCarbonFreeCarbonTrapped();
    void SetCarbonTrappedAndCarbonFreeForGivenXcTot();
    /// @}

    // Node-Cell Connectivity
    std::vector<int> getNodesForCell(int cellIndex) const;

    // Node Position
    Eigen::Vector3f getNodePosition(int nodeIndex) const;

    // FEM-related methods
    Eigen::MatrixXf computeElementStiffness(int cellIndex) const;


    friend class Microstructure;
};

#endif