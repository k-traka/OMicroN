/* stateVars.cpp
   OMicroN (optimising microstructures numerically) simulation program
   cpp-file containing stateVars class implementation
*/

#include <algorithm>
#include <cmath>
#include <cstring>
#include <set>
#include <vector>
#include "stateVars.h"
#include <iostream>

/** @brief Constructs a new instance of StateVars, initialing the grid and the basic state variables, as well as the state variables related to the simulation (i.e. froom user defined settings)
 *
 * @param nx Number of cells along X direction
 * @param ny Number of cells along Y direction
 * @param nz Number of cells along Z direction
 * @param dx Grid spacing in m
 * @param grid Type of grid (e.g. 0 for square/cubic and 1 for hexagonal)

 * @param IncludeSoluteDiffusion Flag signalling if solute diffusion will be simulated
 * @param IncludeRexAndGG Flag signalling if recrystallization and/or grain growth will be simulated
 * @param AllowIntMigration Flag signalling if phase transformations (interface migration between dissimilar phases) simulated

 * @param initXCEqInt Flag signaling whether the carbon equilibrium at local interphase partitioning takes place before or during the diffusion step. If before (ALTHOUGH THIS IS NOT RECOMMENDED) then initXCEqInt = true and should be initialized here
 * @param initXCEqDef Flag signaling whether the carbon equilibrium between free lattice and local defects takes place before or during the diffusion step. If before (ALTHOUGH THIS IS NOT RECOMMENDED) then initXCEqInt = true and should be initialized here
 */

StateVars::StateVars(int nx, int ny, int nz, float dx, int grid, bool IncludeSoluteDiffusion, bool IncludeRexAndGG, bool AllowIntMigration, bool initXCEqInt, bool initXCEqDef)
    : StateVars()
{
    LOG_F(INFO, "Initializing Cells");

    LOG_F(INFO, "Grid of: nx=%d, ny=%d, nz=%d, dx=%g", nx, ny, nz, dx);

    int i, j, k, mindex;
    int nearest[6];   // Relative indices (1D array) of nearest neighbours
    int next[12];     // Relative indices (1D array) of next nearest neighbours
    int next_next[8]; // Relative indices (1D array) of 3rd neighbours
    mvNx = nx;
    mvNy = ny;
    mvNz = nz;
    mvNxy = nx * ny;
    mvN = mvNxy * nz;

    mvGridType = grid;
    mvIncludeRexAndGG = IncludeRexAndGG;
    mvIncludeSoluteDiffusion = IncludeSoluteDiffusion;
    mvAllowIntMigration = AllowIntMigration;

    if (mvGridType == 0)
    {
        if (mvNz == 1)
        { // We are doing a 2D simulation
            mvDx = dx;
            mvDx2 = mvDx;        // This means that areas are treated as lengths
            mvDx3 = mvDx * mvDx; // And volumes as areas!!
            mvIs2DSimulation = true;
            mvNNearest = 4;
            mvNNext = 4;
            mvNNextNext = 0;
            mvNAll = 8;
            mvMaxBoundaryAreaPerVol = floor(1. / mvDx / 2.); // in Sq(1,2) grid the max boundary area is half the side. We assume directly we are dealing with all neighbours, i.e. Sq(1,2),  (for grain growth the square von-neumann grid is bad)
        }
        else
        {
            mvDx = dx;
            mvDx2 = mvDx * mvDx;
            mvDx3 = mvDx2 * mvDx;
            mvIs2DSimulation = false;
            mvNNearest = 6;
            mvNNext = 12;
            mvNNextNext = 8;
            mvNAll = 26;
            mvMaxBoundaryAreaPerVol = floor(1. / mvDx / 4.); // in Cub(1,2) grid the max boundary area is a quarter of the side, and it corresponds to the nearest neighbours.
        }
    }
    else
    {
        if (mvGridType != 1)
        {
            ABORT_F("Unknown grid type : we have only regular (0) and hexagonal (1). CHECK INPUT FILE SETTINGS GridType.");
        }
        if (mvNz != 1)
        {
            ABORT_F("3D Simulation in hexagonal grid not supported. CHECK INPUT FILE SETTINGS GridType and Nz");
        }
        // we are then in 2D hexagonal
        mvDx = dx;
        mvDx2 = mvDx;        // This means that areas are treated as lengths
        mvDx3 = mvDx * mvDx; // And volumes as areas!!
        mvIs2DSimulation = true;
        mvNNearest = 6;
        mvNNext = 0;
        mvNNextNext = 0;
        mvNAll = 6;
        mvMaxBoundaryAreaPerVol = floor(1. / mvDx / 1.5); // in Hex(1,1) grid all boundaries have are per cell volume equal to side/1.5 (regular hexagons)
    }

    // init neighbour lists
    if (mvGridType == 0)
    {
        if (this->mvIs2DSimulation)
        {
            nearest[0] = -1;  /* x - 1 */
            nearest[1] = 1;   /* x + 1 */
            nearest[2] = nx;  /* y + 1 */
            nearest[3] = -nx; /* y - 1 */

            next[0] = -1 - nx; /* x - 1 y - 1*/
            next[1] = 1 - nx;  /* x + 1 y - 1*/
            next[2] = -1 + nx; /* x - 1 y + 1*/
            next[3] = 1 + nx;  /* x + 1 y + 1*/
        }
        else
        {                        // 3D simulation
            nearest[0] = -1;     /* x - 1 */
            nearest[1] = 1;      /* x + 1 */
            nearest[2] = nx;     /* y + 1 */
            nearest[3] = -nx;    /* y - 1 */
            nearest[4] = mvNxy;  /* z + 1 */
            nearest[5] = -mvNxy; /* z - 1 */

            next[0] = -1 - nx; /* x - 1 y - 1*/
            next[1] = 1 - nx;  /* x + 1 y - 1*/
            next[2] = -1 + nx; /* x - 1 y + 1*/
            next[3] = 1 + nx;  /* x + 1 y + 1*/
            next[4] = -1 - mvNxy;
            next[5] = 1 - mvNxy;
            next[6] = -1 + mvNxy;
            next[7] = 1 + mvNxy;
            next[8] = -nx - mvNxy;
            next[9] = nx - mvNxy;
            next[10] = -nx + mvNxy;
            next[11] = nx + mvNxy; /* y + 1 z + 1*/

            next_next[0] = 1 + nx + mvNxy; /* x + 1 y + 1 z + 1*/
            next_next[1] = -1 + nx + mvNxy;
            next_next[2] = 1 - nx + mvNxy;
            next_next[3] = -1 - nx + mvNxy; /* x - 1 y - 1 z + 1*/
            next_next[4] = 1 + nx - mvNxy;
            next_next[5] = -1 + nx - mvNxy;
            next_next[6] = 1 - nx - mvNxy;
            next_next[7] = -1 - nx - mvNxy; /* x - 1 y - 1 z - 1*/
        }
    }
    else
    {
        nearest[0] = -1;      /* x - 1 */
        nearest[1] = 1;       /* x + 1 */
        nearest[2] = nx - 1;  /* y + 1 */
        nearest[3] = nx;      /* y - 1 */
        nearest[4] = -nx - 1; /* y - 1 */
        nearest[5] = -nx;     /* y - 1 */

        next[0] = -2 - nx; /* x - 1 y - 1*/
        next[1] = 1 - nx;  /* x + 1 y - 1*/
        next[2] = -2 + nx; /* x - 1 y + 1*/
        next[3] = 1 + nx;  /* x + 1 y + 1*/
        next[4] = -2 * nx; /* x - 1 y + 1*/
        next[5] = 2 * nx;  /* x + 1 y + 1*/
    }

    /*Set relative neighbour indices (1D array)*/
    for (i = 0; i < mvNNearest; i++)
        mvpAllNeighbours[i] = nearest[i];
    for (i = 0; i < mvNNext; i++)
        mvpAllNeighbours[mvNNearest + i] = next[i];
    for (i = 0; i < mvNNextNext; i++)
        mvpAllNeighbours[mvNNearest + mvNNext + i] = next_next[i];

    // initialize state variables related to all simulations
    InitEBSDRelatedStuff();
    InitializeBasicStateVariables();

    InitBoundaryCells();
    std::cout << " fine 5" << std::endl;

    validateAllNeighbourCells();

    if (mvIncludeRexAndGG || AllowIntMigration)
    {
        InitializeStateVariablesRelatedToInterfaceMotion();
    }

    std::cout << " fine 5" << std::endl;

    if (IncludeSoluteDiffusion)
    {
        allocateCarbonVectors(initXCEqInt, initXCEqDef);
        InitializeSoluteRedistributionAndDefectsVariables();
        InitDiffusivitiesVector();
    }
}

/** @brief Destructor
 *
 * Takes care of freeing memory assigned to pointers
 */
StateVars::~StateVars()
{
    // Check and delete raw pointers
    if (mvpXC)
    {
        delete mvpXC;
        mvpXC = nullptr;
    }
    if (mvpXCTrapped)
    {
        delete mvpXCTrapped;
        mvpXCTrapped = nullptr;
    }
    if (mvpLatticeId)
    {
        delete[] mvpLatticeId;
        mvpLatticeId = nullptr;
    }
    if (mvpBoundaryCell)
    {
        delete[] mvpBoundaryCell;
        mvpBoundaryCell = nullptr;
    }
    if (mvpIsInterface)
    {
        delete[] mvpIsInterface;
        mvpIsInterface = nullptr;
    }
    if (mvpIsInterphase)
    {
        delete[] mvpIsInterphase;
        mvpIsInterphase = nullptr;
    }
    if (mvpKappaFactorForCTrappedDefects)
    {
        delete[] mvpKappaFactorForCTrappedDefects;
        mvpKappaFactorForCTrappedDefects = nullptr;
    }
    if (mvpDiffusivityC_Ratio)
    {
        delete[] mvpDiffusivityC_Ratio;
        mvpDiffusivityC_Ratio = nullptr;
    }
    if (mvpOneOfTheNeighboursGrowingInto)
    {
        delete[] mvpOneOfTheNeighboursGrowingInto;
        mvpOneOfTheNeighboursGrowingInto = nullptr;
    }
    if (mvpConsumptionRate)
    {
        delete[] mvpConsumptionRate;
        mvpConsumptionRate = nullptr;
    }
    if (mvpConsumedFraction)
    {
        delete[] mvpConsumedFraction;
        mvpConsumedFraction = nullptr;
    }
    if (mvpCSLRelationshipNow)
    {
        delete[] mvpCSLRelationshipNow;
        mvpCSLRelationshipNow = nullptr;
    }
    if (mvpCI)
    {
        delete[] mvpCI;
        mvpCI = nullptr;
    }
    if (mvpHasNeighNonInd)
    {
        delete[] mvpHasNeighNonInd;
        mvpHasNeighNonInd = nullptr;
    }
    if (mvpKAM)
    {
        delete[] mvpKAM;
        mvpKAM = nullptr;
    }
    if (mvpRho)
    {
        delete[] mvpRho;
        mvpRho = nullptr;
    }
    if (mvpRX)
    {
        delete[] mvpRX;
        mvpRX = nullptr;
    }
    if (mvpMaxGbPassed)
    {
        delete[] mvpMaxGbPassed;
        mvpMaxGbPassed = nullptr;
    }
}


/** @brief Initializes pointers with state variables related any type of simulation.
 */
void StateVars::InitializeBasicStateVariables()
{
    mvpBoundaryCell = new unsigned char[mvN];

    if (!mvpLatticeId)
    {
        mvpLatticeId = new int[mvN];
    }

    mvpOriId.resize(mvN);
    std::fill(mvpOriId.begin(), mvpOriId.end(), -1);

    if (!mvpIsInterface)
    {
        mvpIsInterface = new int[mvN];
        std::fill(mvpIsInterface, mvpIsInterface + mvN, 0); // Fill with zeros
    }
    if (!mvpIsInterphase)
    {
        mvpIsInterphase = new int[mvN];
        std::fill(mvpIsInterphase, mvpIsInterphase + mvN, 0); // Fill with zeros
    }
}

/** @brief Initializes pointers with state variables related to simulating interface migration (i.e. recrystallization / grain growth / phase transformation).
 */
void StateVars::InitializeStateVariablesRelatedToInterfaceMotion()
{
    mvpConsumptionRate = new float[mvN];
    std::fill(mvpConsumptionRate, mvpConsumptionRate + mvN, 0.); // Fill with zeros
    mvpConsumedFraction = new float[mvN];
    std::fill(mvpConsumedFraction, mvpConsumedFraction + mvN, 0.); // Fill with zeros
    mvpOneOfTheNeighboursGrowingInto = new int[mvN];
    std::fill(mvpOneOfTheNeighboursGrowingInto, mvpOneOfTheNeighboursGrowingInto + mvN, -1); // Fill with -1
    //  the following for analysis after simulation
    mvpMaxGbPassed = new float[mvN];
    std::fill(mvpMaxGbPassed, mvpMaxGbPassed + mvN, 0.); // Fill with 0
}

/** @brief Initializes pointer containing information on whether a cell lies on the outer edges (or faces for 3D) of the RVE
 */
void StateVars::InitBoundaryCells()
{
    /*Init boundary cells*/
    unsigned char btype;
    int mindex = 0;

    for (int k = 0; k < mvNz; k++)
    {
        for (int j = 0; j < mvNy; j++)
        {
            for (int i = 0; i < mvNx; i++)
            {
                btype = 0;
                if (i == 0)
                    btype |= BOUNDARIES::X_LOWER;
                if (i == mvNx - 1)
                    btype |= BOUNDARIES::X_UPPER;
                if (j == 0)
                    btype |= BOUNDARIES::Y_LOWER;
                if (j == mvNy - 1)
                    btype |= BOUNDARIES::Y_UPPER;
                if (!mvIs2DSimulation)
                {
                    if (k == 0)
                        btype |= BOUNDARIES::Z_LOWER;
                    if (k == mvNz - 1)
                        btype |= BOUNDARIES::Z_UPPER;
                }
                if (btype)
                    mvpBoundaryCell[mindex] = btype;
                mindex++;
            }
        }
    }
}

/** @brief Initializes pointers with state variables which are necessary when reading the starting microstructure file and especially if it comes from EBSD
 */
void StateVars::InitEBSDRelatedStuff()
{
    if (mvpKAM == nullptr)
        mvpKAM = new float[mvN];
    if (mvpCI == nullptr)
        mvpCI = new float[mvN];
    if (mvpRX == nullptr)
    {
        mvpRX = new int[mvN];
    }
    if (mvpRho == nullptr)
        mvpRho = new float[mvN];
    std::memset(mvpRho, 1.e12, mvN);
    if (mvpHasNeighNonInd == nullptr)
        mvpHasNeighNonInd = new int[mvN];
}

/** @brief Initialises measurement-related quantities when microstructure comes from EBSD
 * *
 * @param index the index of cell for which the state variables will be set
 * @param CI the confidence index of the measured pixel
 * @param RX state variable regarding the state of the pixel (normally non-recrystallized for deformed inputs) - unless we read partially recrystallized state or state after previous simulation
 * @param rho the dislocation density measured / calculated / simulated from other software

 */
void StateVars::SetEBSDRelatedStuffForRexAndGG(int index, float CI, int RX, float rho)
{
    mvpCI[index] = CI;
    mvpRX[index] = RX;
    mvpRho[index] = rho;
}

/** @brief Allocates memory for carbon concentration of all cells based
 *
 * Use this function to reserve memory for the carbon concentration of all cells.
 *
 * @param *mvpXC Pointer to Eigen vector to store carbon concentration information
 * @param *mvpXC Trapped Pointer to Eigen vector to store carbon concentration trapped
 * @param *mvpXCEqAtInterphase (USE ONLY IF THERE IS ANY PROBLEM ENCOUNTERED WITH DIFFUSION EQUATION - otherwise local equilibrium partitioning takes place later) Pointer to Eigen vector to store carbon concentration equilibrium for phase partitioning at interphase cells
 * @param *mvpXCTrappedEqWithNbr (USE ONLY IF THERE IS ANY PROBLEM ENCOUNTERED WITH DIFFUSION EQUATION - otherwise local equilibrium trapping/free takes place later) Trapped Pointer to Eigen vector to store carbon concentration trapped equilibrium
 */
void StateVars::allocateCarbonVectors(bool initXCEqInt, bool initXCEqDef)
{
    /*First check if already allocated */
    if (mvpXC && mvpXCTrapped)
        return;
    if (!mvpXC)
    {
        mvpXC = new Eigen::VectorXd(mvN);
        mvpXC->setZero(mvN);
    }
    if (!mvpXCTrapped)
    {
        mvpXCTrapped = new Eigen::VectorXd(mvN);
        mvpXCTrapped->setZero(mvN);
    }
    if (initXCEqInt)
    {
        mvpXCEqAtInterphase = new Eigen::VectorXd(mvN);
        mvpXCEqAtInterphase->setZero(mvN);
        LOG_F(WARNING, " It is not advisable to do partitioning outside the diffusion solver, as this can cause (minor) mesh dependencies. Anyway i initialize the relevant variables, but make sure this is what you want. If not change the input file parameter IsPartitioningHappeningInDiffusionStep into 1");
    }
    if (initXCEqDef)
    {
        mvpXCTrappedEqWithNbr = new Eigen::VectorXd(mvN);
        mvpXCTrappedEqWithNbr->setZero(mvN);
        LOG_F(WARNING, " It is not advisable to do trapping to defects equilibrium outside the diffusion solver, as this can cause (minor) mesh dependencies. Anyway i initialize the relevant variables, but make sure this is what you want. If not change the input file parameter IsSoluteSegregationHappeningInDiffusionStep into 1");
    }
}

/** @brief Initializes pointers with state variables which are necessary for solute redistrubution (i.e. diffusivities per pair of cells,  and trapping to defects per cell).
 */

void StateVars::InitializeSoluteRedistributionAndDefectsVariables()
{
    mvpKappaFactorForCTrappedDefects = new double[mvN];   // although it only matters for Csegregation to defects, the way the code is writtem now whole Cdiffusion step is solved
    mvpDiffusivityC_Ratio = new double[mvN * mvNNearest]; // although it only matters for having different diffusivities (e.g. Agren, or solving martensite-austenite), the way the Ax=b code is written now  needs the Diff matrix (i.e. it is part of A matrix, and takes values of 1 if all diffusivitties are equal (i.e. equal to the maxdiffusivity))
    std::memset(mvpKappaFactorForCTrappedDefects, 0., mvN);
    std::memset(mvpDiffusivityC_Ratio, 1., mvN * mvNNearest);
}








/** @brief Calculates in place the pointer to the array with indices of all neighbour cells, depending on the grid and settings used in the simulation
 *
 * @param index Index of cell for which to return the neighbour cells
 * @param *all_p Pointer to array of offsets for all neighbours
 */

void StateVars::GetNeighbourCellOffsets(const int index, int *all_p) const
{
    int i;
    int ci, cj, ck;
    unsigned char btype;

    for (i = 0; i < mvNAll; i++)
        all_p[i] = index + mvpAllNeighbours[i];
    /* take care of periodic boundaries*/
    btype = mvpBoundaryCell[index];

    IndexToIJK(index, &ci, &cj, &ck);
    // do this before returning - in hex grid
    if (mvGridType == 1 && cj % 2 != 0)
    {
        all_p[2] += 1;
        all_p[3] += 1;
        all_p[4] += 1;
        all_p[5] += 1;
    }

    if (!btype)
        return;

    if (mvGridType == 0)
    {
        if (mvIs2DSimulation)
        {
            if (btype & BOUNDARIES::X_LOWER)
            {
                all_p[0] += mvNx;
                all_p[4] += mvNx;
                all_p[6] += mvNx;
            }
            if (btype & BOUNDARIES::X_UPPER)
            {
                all_p[1] -= mvNx;
                all_p[5] -= mvNx;
                all_p[7] -= mvNx;
            }
            if (btype & BOUNDARIES::Y_UPPER)
            {
                all_p[2] -= mvNxy;
                all_p[6] -= mvNxy;
                all_p[7] -= mvNxy;
            }
            if (btype & BOUNDARIES::Y_LOWER)
            {
                all_p[3] += mvNxy;
                all_p[4] += mvNxy;
                all_p[5] += mvNxy;
            }
        }
        else
        {
            if (btype & BOUNDARIES::X_LOWER)
            {
                all_p[0] += mvNx;
                all_p[6 + 0] += mvNx;
                all_p[6 + 2] += mvNx;
                all_p[6 + 4] += mvNx;
                all_p[6 + 6] += mvNx;
                all_p[6 + 12 + 1] += mvNx;
                all_p[6 + 12 + 3] += mvNx;
                all_p[6 + 12 + 5] += mvNx;
                all_p[6 + 12 + 7] += mvNx;
            }
            if (btype & BOUNDARIES::X_UPPER)
            {
                all_p[1] -= mvNx;
                all_p[6 + 1] -= mvNx;
                all_p[6 + 3] -= mvNx;
                all_p[6 + 5] -= mvNx;
                all_p[6 + 7] -= mvNx;
                all_p[6 + 12 + 0] -= mvNx;
                all_p[6 + 12 + 2] -= mvNx;
                all_p[6 + 12 + 4] -= mvNx;
                all_p[6 + 12 + 6] -= mvNx;
            }
            if (btype & BOUNDARIES::Y_UPPER)
            {
                all_p[2] -= mvNxy;
                all_p[6 + 2] -= mvNxy;
                all_p[6 + 3] -= mvNxy;
                all_p[6 + 9] -= mvNxy;
                all_p[6 + 11] -= mvNxy;
                all_p[6 + 12 + 0] -= mvNxy;
                all_p[6 + 12 + 1] -= mvNxy;
                all_p[6 + 12 + 4] -= mvNxy;
                all_p[6 + 12 + 5] -= mvNxy;
            }
            if (btype & BOUNDARIES::Y_LOWER)
            {
                all_p[3] += mvNxy;
                all_p[6 + 0] += mvNxy;
                all_p[6 + 1] += mvNxy;
                all_p[6 + 8] += mvNxy;
                all_p[6 + 10] += mvNxy;
                all_p[6 + 12 + 2] += mvNxy;
                all_p[6 + 12 + 6] += mvNxy;
                all_p[6 + 12 + 3] += mvNxy;
                all_p[6 + 12 + 7] += mvNxy;
            }
            if (btype & BOUNDARIES::Z_UPPER)
            {
                all_p[4] -= mvN;
                all_p[6 + 6] -= mvN;
                all_p[6 + 7] -= mvN;
                all_p[6 + 10] -= mvN;
                all_p[6 + 11] -= mvN;
                all_p[6 + 12 + 0] -= mvN;
                all_p[6 + 12 + 1] -= mvN;
                all_p[6 + 12 + 2] -= mvN;
                all_p[6 + 12 + 3] -= mvN;
            }
            if (btype & BOUNDARIES::Z_LOWER)
            {
                all_p[5] += mvN;
                all_p[6 + 4] += mvN;
                all_p[6 + 5] += mvN;
                all_p[6 + 8] += mvN;
                all_p[6 + 9] += mvN;
                all_p[6 + 12 + 4] += mvN;
                all_p[6 + 12 + 5] += mvN;
                all_p[6 + 12 + 6] += mvN;
                all_p[6 + 12 + 7] += mvN;
            }
        }
    }
    else
    {
        if (btype & BOUNDARIES::X_LOWER)
        {
            all_p[0] += mvNx;
            if (cj % 2 == 0)
            {
                all_p[2] += mvNx;
                all_p[4] += mvNx;
            }
        }
        if (btype & BOUNDARIES::X_UPPER)
        {
            all_p[1] -= mvNx;
            if (cj % 2 != 0)
            {
                all_p[3] -= mvNx;
                all_p[5] -= mvNx;
            }
        }
        if (btype & BOUNDARIES::Y_LOWER)
        {
            all_p[4] += mvNxy;
            all_p[5] += mvNxy;
        }
        if (btype & BOUNDARIES::Y_UPPER)
        {
            all_p[2] -= mvNxy;
            all_p[3] -= mvNxy;
        }
    }
}

/** @brief Returns a vector with indices of all neighbour cells, depending on the grid and settings used in the simulation
 *
 * @param index Index of cell for which to return the neighbour cells
 *
 * @return A vector with all neighbours
 */
std::vector<int> StateVars::GetAllNeighbourCells(int index) const
{
    int i;
    int ci, cj, ck;
    unsigned char btype;
    std::vector<int> nbrVec;

    for (i = 0; i < mvNAll; i++)
        nbrVec.push_back(index + mvpAllNeighbours[i]);
    /* take care of periodic boundaries*/

    IndexToIJK(index, &ci, &cj, &ck);
    if (mvGridType == 1 && cj % 2 != 0)
    {
        nbrVec[2] += 1;
        nbrVec[3] += 1;
        nbrVec[4] += 1;
        nbrVec[5] += 1;
    }

    btype = mvpBoundaryCell[index];
    if (!btype)
        return nbrVec;

    if (mvGridType == 0)
    {
        if (mvIs2DSimulation)
        {
            if (btype & BOUNDARIES::X_LOWER)
            {
                nbrVec[0] += mvNx;
                nbrVec[4] += mvNx;
                nbrVec[6] += mvNx;
            }
            if (btype & BOUNDARIES::X_UPPER)
            {
                nbrVec[1] -= mvNx;
                nbrVec[5] -= mvNx;
                nbrVec[7] -= mvNx;
            }
            if (btype & BOUNDARIES::Y_UPPER)
            {
                nbrVec[2] -= mvNxy;
                nbrVec[6] -= mvNxy;
                nbrVec[7] -= mvNxy;
            }
            if (btype & BOUNDARIES::Y_LOWER)
            {
                nbrVec[3] += mvNxy;
                nbrVec[4] += mvNxy;
                nbrVec[5] += mvNxy;
            }
        }
        else
        {
            if (btype & BOUNDARIES::X_LOWER)
            {
                nbrVec[0] += mvNx;
                nbrVec[6 + 0] += mvNx;
                nbrVec[6 + 2] += mvNx;
                nbrVec[6 + 4] += mvNx;
                nbrVec[6 + 6] += mvNx;
                nbrVec[6 + 12 + 1] += mvNx;
                nbrVec[6 + 12 + 3] += mvNx;
                nbrVec[6 + 12 + 5] += mvNx;
                nbrVec[6 + 12 + 7] += mvNx;
            }
            if (btype & BOUNDARIES::X_UPPER)
            {
                nbrVec[1] -= mvNx;
                nbrVec[6 + 1] -= mvNx;
                nbrVec[6 + 3] -= mvNx;
                nbrVec[6 + 5] -= mvNx;
                nbrVec[6 + 7] -= mvNx;
                nbrVec[6 + 12 + 0] -= mvNx;
                nbrVec[6 + 12 + 2] -= mvNx;
                nbrVec[6 + 12 + 4] -= mvNx;
                nbrVec[6 + 12 + 6] -= mvNx;
            }
            if (btype & BOUNDARIES::Y_UPPER)
            {
                nbrVec[2] -= mvNxy;
                nbrVec[6 + 2] -= mvNxy;
                nbrVec[6 + 3] -= mvNxy;
                nbrVec[6 + 9] -= mvNxy;
                nbrVec[6 + 11] -= mvNxy;
                nbrVec[6 + 12 + 0] -= mvNxy;
                nbrVec[6 + 12 + 1] -= mvNxy;
                nbrVec[6 + 12 + 4] -= mvNxy;
                nbrVec[6 + 12 + 5] -= mvNxy;
            }
            if (btype & BOUNDARIES::Y_LOWER)
            {
                nbrVec[3] += mvNxy;
                nbrVec[6 + 0] += mvNxy;
                nbrVec[6 + 1] += mvNxy;
                nbrVec[6 + 8] += mvNxy;
                nbrVec[6 + 10] += mvNxy;
                nbrVec[6 + 12 + 2] += mvNxy;
                nbrVec[6 + 12 + 6] += mvNxy;
                nbrVec[6 + 12 + 3] += mvNxy;
                nbrVec[6 + 12 + 7] += mvNxy;
            }
            if (btype & BOUNDARIES::Z_UPPER)
            {
                nbrVec[4] -= mvN;
                nbrVec[6 + 6] -= mvN;
                nbrVec[6 + 7] -= mvN;
                nbrVec[6 + 10] -= mvN;
                nbrVec[6 + 11] -= mvN;
                nbrVec[6 + 12 + 0] -= mvN;
                nbrVec[6 + 12 + 1] -= mvN;
                nbrVec[6 + 12 + 2] -= mvN;
                nbrVec[6 + 12 + 3] -= mvN;
            }
            if (btype & BOUNDARIES::Z_LOWER)
            {
                nbrVec[5] += mvN;
                nbrVec[6 + 4] += mvN;
                nbrVec[6 + 5] += mvN;
                nbrVec[6 + 8] += mvN;
                nbrVec[6 + 9] += mvN;
                nbrVec[6 + 12 + 4] += mvN;
                nbrVec[6 + 12 + 5] += mvN;
                nbrVec[6 + 12 + 6] += mvN;
                nbrVec[6 + 12 + 7] += mvN;
            }
        }
    }
    else
    {
        if (btype & BOUNDARIES::X_LOWER)
        {
            nbrVec[0] += mvNx;
            if (cj % 2 == 0)
            {
                nbrVec[2] += mvNx;
                nbrVec[4] += mvNx;
            }
        }
        if (btype & BOUNDARIES::X_UPPER)
        {
            nbrVec[1] -= mvNx;
            if (cj % 2 != 0)
            {
                nbrVec[3] -= mvNx;
                nbrVec[5] -= mvNx;
            }
        }
        if (btype & BOUNDARIES::Y_LOWER)
        {
            nbrVec[4] += mvNxy;
            nbrVec[5] += mvNxy;
        }
        if (btype & BOUNDARIES::Y_UPPER)
        {
            nbrVec[2] -= mvNxy;
            nbrVec[3] -= mvNxy;
        }
    }
    return nbrVec;
}



/** @brief Converts X, Y, Z grid positions (i, j, k) to cell index
 *
 * @param i X coordinate as grid position (dimensionless)
 * @param j Y coordinate as grid position (dimensionless)
 * @param k Z coordinate as grid position (dimensionless)
 *
 * @return Cell index
 */
int StateVars::IJKToIndex(int i, int j, int k) const
{
    return i + j * mvNx + k * mvNxy;
}

/** @brief Converts cell index to X, Y, Z grid positions (i, j, k)
 *
 * @param index Cell index
 * @param *i X coordinate as grid position returned by pointer value
 * @param *j Y coordinate as grid position returned by pointer value
 * @param *k Z coordinate as grid position returned by pointer value
 */
void StateVars::IndexToIJK(int index, int *i, int *j, int *k) const
{
    *k = index / mvNxy;
    *j = (index - (*k) * mvNxy) / mvNx;
    *i = index % mvNx;
}

/** @brief Similar to Index2IJK, but returns i, j, k as std:array
 *
 * @param index Cell index
 *
 * @return X, Y, Z coordinates as grid positions as std::array
 */
std::array<int, 3> StateVars::GetIJKFromIndex(int index) const
{
    int i, j, k;
    IndexToIJK(index, &i, &j, &k);
    return {i, j, k};
}

/** @brief Gets the X coordinate (in m) from the cell index
 *
 * @param index Cell index
 *
 * @return The X coordinate in m
 */
float StateVars::GetXFromIndex(int index) const
{
    int k = index / mvNxy;
    int j = (index - k * mvNxy) / mvNx;
    int i = index - k * mvNxy - j * mvNx;
    if (mvGridType == 0 || j % 2 == 0)
        return (i * mvDx);
    else
        return (i * mvDx + mvDx / 2.); // in Hexagonal grid odd numbered rows (starting from 0) are shifted to the right by dx/2}
}

/** @brief Gets the Y coordinate (in m) from the cell index
 *
 * @param index Cell index
 *
 * @return The Y coordinate in m
 */
float StateVars::GetYFromIndex(int index) const
{
    int k = index / mvNxy;
    int j = (index - k * mvNxy) / mvNx;
    if (mvGridType == 0)
        return (j * mvDx);
    return (j * mvDx * pow(0.75, 0.5));
}

/** @brief Gets the Z coordinate (in m) from the cell index
 *
 * @param index Cell index
 *
 * @return The Z coordinate in m
 */
float StateVars::GetZFromIndex(int index) const
{
    int k = index / mvNxy;
    return (k * mvDx);
}

/** @brief Gets the XYZ coordinates (in m) from the cell index
 *
 * @param index Cell index
 *
 * @return The XYZ coordinates in m as an std::array
 */
std::array<float, 3> StateVars::GetXYZFromIndex(float index) const
{
    return {GetXFromIndex(index), GetYFromIndex(index), GetZFromIndex(index)};
}

void StateVars::validateAllNeighbourCells(void)
{
    int nbr, index, indexo, i, j, k;
    int lall[mvNAll];

    for (index = 0; index < mvN; index++)
    {
        GetNeighbourCellOffsets(index, lall);
        IndexToIJK(index, &i, &j, &k);
        for (nbr = 0; nbr < mvNAll; nbr++)
        {
            indexo = lall[nbr];

            if (mvIs2DSimulation && indexo < 0 && !(i <= 1 || j <= 1 || i >= mvNx - 2 || j >= mvNy - 2))
            {
                LOG_F(WARNING, "neighbour: %d of %d too low: %d; i %d j %d k %d", nbr, index, indexo, i, j, k);
            }

            if (mvIs2DSimulation && indexo >= mvN && !mvpBoundaryCell[index] && !(i <= 1 || j <= 1 || i >= mvNx - 2 || j >= mvNy - 2))
            {
                LOG_F(WARNING, "neighbour: %d of %d too low: %d; i %d j %d k %d", nbr, index, indexo, i, j, k);
            }

            if (!mvIs2DSimulation && indexo < 0 && !mvpBoundaryCell[index] /**(x<=1 || y<=1 || x>=mvNx-2 || y>=mvNy-2 || z<=1 || z>=mvNz-1)*/)
            {
                LOG_F(WARNING, "neighbour: %d of %d too high: %d max: %d; i %d j %d k %d", nbr, index, indexo, mvN, i, j, k);
            }

            if (!mvIs2DSimulation && indexo >= mvN && !mvpBoundaryCell[index] /*&& !(x<=1 || y<=1 || x>=mvNx-2 || y>=mvNy-2 || z<=1 || z>=mvNz-1)*/)
            {
                LOG_F(WARNING, "neighbour: %d of %d too high: %d max: %d; i %d j %d k %d", nbr, index, indexo, mvN, i, j, k);
            }
        }
    }
}

/// @brief Initializes #mvpDiffusivityC_Ratio containing diffusivities per pair of neighbour cells through which solute diffuses, in matrix form
void StateVars::InitDiffusivitiesVector()
{
    mvpDiffusivityC_Ratio = new double[mvN * mvNNearest];
    int i, j;
    for (i = 0; i < mvN; ++i)
    {
        for (j = 0; j < mvNNearest; ++j)
        {
            mvpDiffusivityC_Ratio[mvNNearest * i + j] = 1.; // if you dont want variability in diffusivities the values will remain 1. But still, initialize it pls even because (for now) when solving Ax=b, A is matrix considering diffusivities for each neighbour pair.
        }
    }
}

/** @brief Gets the distance between two cells in m (taking care of periodic boundaries)
 *
 * @param indexA index of first cell
 * @param indexB index of second cell
 *
 * @return distance between the two cells in m
 */

// Changed to do directly di dj dk for xyz insteaf of ijk (and not only in return) - it is easier in case we have also hex grid

float StateVars::GetDistanceBetweenCells(int indexA, int indexB) const
{
    float di, dj, dk;
    float i, j, k, iB, jB, kB;

    i = GetXFromIndex(indexA);
    j = GetYFromIndex(indexA);
    k = GetZFromIndex(indexA);
    iB = GetXFromIndex(indexB);
    jB = GetYFromIndex(indexB);
    kB = GetZFromIndex(indexB);
    di = i - iB;
    dj = j - jB;
    dk = k - kB;

    if (di > mvDx * 0.5 * mvNx)
        di -= mvDx * mvNx;
    else if (di < -mvDx * 0.5 * mvNx)
        di += mvDx * mvNx;
    if (mvGridType == 0)
    {
        if (dj > mvDx * 0.5 * mvNy)
            dj -= mvDx * mvNy;
        else if (dj < -mvDx * 0.5 * mvNy)
            dj += mvDx * mvNy;
        if (dk > mvDx * 0.5 * mvNz)
            dk -= mvDx * mvNz;
        else if (dk < -mvDx * 0.5 * mvNz)
            dk += mvDx * mvNz;
    }
    else
    {
        if (dj > mvDx * 0.5 * mvNy)
            dj -= pow(0.75, 0.5) * mvDx * mvNy;
        else if (dj < -mvDx * 0.5 * mvNy)
            dj += pow(0.75, 0.5) * mvDx * mvNy;
        dk = 0.;
    }
    return sqrt(di * di + dj * dj + dk * dk);
}

/** @brief Solves one diffusion step
 *
 * @param TCK_p pointer to ThermChemKin instance
 * @param dt predetermined time step
 * @param AllowSoluteSegregation bool on whether we have solute trapping in general during simulation
 * @param maxDiffusivityInTimeStep the maximum value of diffusivity between neighbour cells calculated in microstructure class
 * @param IsPartitioningHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 * @param IsSoluteSegregationHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 */

void StateVars::SoluteDiffusionStep(ThermChemKin *TCK_p, double dt, bool AllowSoluteSegregation, double maxDiffusivityInTimeStep, bool IsPartitioningHappeningHere, bool IsSoluteSegregationHappeningHere)
{
    LOG_F(INFO, "Solving solute diffusion mvAlphaDt %g, maxDiffusivityInTimeStep %g, dt %g, dx %g", mvAlphaDt, maxDiffusivityInTimeStep, dt, mvDx);
    LOG_IF_F(INFO, IsPartitioningHappeningHere, "Including carbon partitioning in Ax = b");
    LOG_IF_F(INFO, !IsPartitioningHappeningHere, "Carbon partitioning already solved (not part of Ax = b)");
    LOG_IF_F(INFO, IsSoluteSegregationHappeningHere && AllowSoluteSegregation, "Including carbon trapping in Ax = b");
    LOG_IF_F(INFO, !IsSoluteSegregationHappeningHere && AllowSoluteSegregation, "Carbon trapping already solved (not part of Ax = b)");
    LOG_IF_F(INFO, !mvKeepConstantIronAtomsPerCell, "Constant number of moles per cell; XCgradient means XC [at.fraction.]");
    LOG_IF_F(INFO, mvKeepConstantIronAtomsPerCell, "Constant number of Fe atoms per cell; XCgradient means XC [carbon / Fe]");
    int maxit = 500;
    double tolerance = 5e-5; // relative tolerance: iterative linear solver stops at iteration k when
    mvAlphaDt = dt * maxDiffusivityInTimeStep / (mvDx * mvDx);

    if (AllowSoluteSegregation && !IsSoluteSegregationHappeningHere)
    {
        *mvpXC -= *mvpXCTrapped;
    }

    bool IsSystemAsBSolvedWithCarbonTrapping = (AllowSoluteSegregation && IsSoluteSegregationHappeningHere);

    if (mvKeepConstantIronAtomsPerCell)
        ConvertTotCellConcentrationsInCarbonPerIron();

    *mvpXC = NumericalSolverCG(*mvpXC, TCK_p, 500, 5e-5, *mvpXC, IsPartitioningHappeningHere, IsSystemAsBSolvedWithCarbonTrapping);
    if (mvKeepConstantIronAtomsPerCell)
        ConvertTotCellConcentrationsInAtFraction();
    if (AllowSoluteSegregation)
    {
        if (!IsSoluteSegregationHappeningHere)
            PutBackCarbonTrappedAndCalculateNewCarbonFreeCarbonTrapped();
        else
            SetCarbonTrappedAndCarbonFreeForGivenXcTot();
    }
}

void StateVars::PutBackCarbonTrappedAndCalculateNewCarbonFreeCarbonTrapped()
{
    *mvpXC += *mvpXCTrapped;
    SetCarbonTrappedAndCarbonFreeForGivenXcTot();
}

void StateVars::SetCarbonTrappedAndCarbonFreeForGivenXcTot()
{
    for (int i = 0; i < mvN; i++)
    {
        if (mvpKappaFactorForCTrappedDefects[i] > 0.)
        {
            double FactorTot = mvpKappaFactorForCTrappedDefects[i] + 1.;
            (*mvpXCTrapped)[i] = (*mvpXC)[i] * mvpKappaFactorForCTrappedDefects[i] / FactorTot;
        }
        else
        {
            (*mvpXCTrapped)[i] = 0.;
        }
        if ((*mvpXC)[i] < 0. || (*mvpXCTrapped)[i] < 0. || (*mvpXC)[i] - (*mvpXCTrapped)[i] < 0.)
            LOG_F(ERROR, " wrong calculations check cell %d, has XCTot %e, XCFree %e, mvpXCTrapped %e, mvpKappaFactorForCTrappedDefects %e, mvpKappaFactorForCTrappedDefects*rho %e ", i, (*mvpXC)[i], (*mvpXC)[i] - (*mvpXCTrapped)[i], (*mvpXCTrapped)[i], mvpKappaFactorForCTrappedDefects[i], mvpKappaFactorForCTrappedDefects[i] / mvpRho[i] / 1.e12);
    }
}

/** @brief Numerical solution of diffusion step based on conjugate gradient method
 * @param b the b=Ax (which is equal to current concentratons initially)
 * @param TCK_p pointer to ThermChemKin instance
 * @param maxit maximum allowed iterations
 * @param tolerance tolerance for convergence
 * @param initial_guess the first x values guessed (set here as current solute concentrations)
 * @param IsPartitioningHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 * @param IsSoluteSegregationHappeningHere bool on whether the diffusion equation includes the solute trapping equilibrium (i.e. if trapping is part of the numerical system Ax=b)
 */

Eigen::VectorXd StateVars::NumericalSolverCG(const Eigen::VectorXd &b, ThermChemKin *TCK_p, int maxit, double tolerance, const Eigen::VectorXd &initial_guess, bool IsPartitioningHappeningHere, bool IsSoluteSegregationHappeningHere)
{
    Eigen::VectorXd x = initial_guess;
    Eigen::VectorXd r(b.size());
    Eigen::VectorXd p(b.size());
    Eigen::VectorXd Ap(b.size());
    double rsold, alpha, rsnew, p_dotAp;

    // Compute the initial residual r = b - A*x for the initial guess
    SystemATimesX(x, TCK_p, Ap, p_dotAp, IsPartitioningHappeningHere, IsSoluteSegregationHappeningHere);
    r = b - Ap; // Initial residual
    p = r;
    rsold = r.dot(r);

    if (rsold < tolerance * tolerance)
    {
        LOG_F(INFO, "Initial guess is already good enough.");
        return x;
    }

    for (int i = 0; i < maxit; ++i)
    {

        SystemATimesX(p, TCK_p, Ap, p_dotAp, IsPartitioningHappeningHere, IsSoluteSegregationHappeningHere);
        alpha = rsold / p_dotAp;
        x += alpha * p;
        r -= alpha * Ap;
        rsnew = r.dot(r);

        // LOG_F(INFO, "Iteration %d: residual = %e", i, sqrt(rsnew));

        if (sqrt(rsnew) < tolerance)
        {
            // LOG_F(INFO, "Converged after %d iterations.", i);
            break;
        }

        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    }

    return x;
}

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

void StateVars::SystemATimesX(const Eigen::VectorXd &x, ThermChemKin *TCK_p, Eigen::VectorXd &b, double &p_dotAp, bool IsPartitioningHappeningHere, bool IsSoluteSegregationHappeningHere)
{
    int lall[mvNAll];
    double XcAlphaEq, XcGammaEq;
    // LOG_F(INFO," mvAlphaDt %g",mvAlphaDt);
    for (int i = 0; i < mvN; i++)
    {
        int j;
        float NbrCounter;
        double RatioXcFreeToXcTotForIndex = 1.;
        if (IsSoluteSegregationHappeningHere)
            RatioXcFreeToXcTotForIndex /= (1. + mvpKappaFactorForCTrappedDefects[i]);

        b[i] = 0.0;
        GetNeighbourCellOffsets(i, lall);
        NbrCounter = 0;
        for (int nbr = 0; nbr < mvNNearest; ++nbr)
        {
            j = lall[nbr];
            double RatioXcFreeToXcTotForNbr = 1.;
            if (IsSoluteSegregationHappeningHere)
                RatioXcFreeToXcTotForNbr /= (1. + mvpKappaFactorForCTrappedDefects[j]);

            if (mvpLatticeId[j] != mvpLatticeId[i])
            {
                if (!IsPartitioningHappeningHere)
                    continue;
                double XcFCC = (IsCellFCC(i)) ? RatioXcFreeToXcTotForIndex * x[i] : RatioXcFreeToXcTotForNbr * x[j];
                double XcBCC = (IsCellFCC(i)) ? RatioXcFreeToXcTotForNbr * x[j] : RatioXcFreeToXcTotForIndex * x[i];
                double XcTot = 0.5 * (XcBCC + XcFCC);
                XcGammaEq = TCK_p->GetEqXcFCCAtThisInterface(XcBCC, XcFCC);
                XcAlphaEq = 2. * XcTot - XcGammaEq;
                double XcPartitionedAtBoundary = (IsCellFCC(i)) ? XcGammaEq : XcAlphaEq;
                if (!IsTheNumberFinite(XcGammaEq) || !IsTheNumberFinite(XcAlphaEq))
                    ABORT_F(" XcGammaEq %g and XcAlphaEq %g", XcGammaEq, XcAlphaEq);
                // LOG_F(INFO," partitioning solution FCC %g and BCC %g xcPart %g", XcGammaEq,XcAlphaEq, XcPartitionedAtBoundary);

                b[i] -= (mvpDiffusivityC_Ratio[mvNNearest * i + nbr] * XcPartitionedAtBoundary);
                NbrCounter += mvpDiffusivityC_Ratio[mvNNearest * i + nbr];
                continue;
            }
            b[i] -= (mvpDiffusivityC_Ratio[mvNNearest * i + nbr] * x[j] * RatioXcFreeToXcTotForNbr);
            NbrCounter += mvpDiffusivityC_Ratio[mvNNearest * i + nbr];
        }
        b[i] *= mvAlphaDt;
        b[i] += (NbrCounter * mvAlphaDt * RatioXcFreeToXcTotForIndex + 1) * x[i];
        if (!IsTheNumberFinite(b[i]))
            ABORT_F(" b[i] %g ", b[i]);
    }
    p_dotAp = b.dot(x);
}

bool StateVars::IsTheNumberFinite(double numberToCheck)
{
    if (numberToCheck > std::numeric_limits<double>::max() || numberToCheck < std::numeric_limits<double>::lowest())
        return false;
    return true;
}
