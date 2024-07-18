/* ThermoChemKinetics.cpp
   OMicroN (optimising microstructures numerically) simulation program
   header file containing ThermChemKin class definitions and implementation
*/


#ifndef ThermoChemKinetics_H
#define ThermoChemKinetics_H

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "hdf5_utils.h"
#include <cmath>



/// @name Helper constants
/// @{
#define K_BOLTZMANN 1.3806504e-23 ///< Boltzmann constant (J/K)
#define R_GAS 8.314472 ///< Gas constant [m3 Pa K-1 mol-1]
/// @}

class ThermChemKin {
public:
    // Constructor
    ThermChemKin(double userXCo,double userStartTemperature, double userGB_Mo,double userGB_Qg,double userGB_E,double userPB_Mo,double userPB_Qg,double userPB_E, double user_DiffPreFactor, double user_DiffQg);

    // Destructor
    ~ThermChemKin();

    // Member functions
    void SetTemperatureForThisStep(double temp){mvT = temp;}
    double GetTemperatureForThisStep(){return mvT;}

    double GetGrainBoundaryMobility() const { return (mvGrainBoundaryM0 * std::exp(-mvGrainBoundaryQg / (R_GAS * mvT))); }
    double GetPhaseBoundaryMobility() const { return (mvPhaseBoundaryM0 * std::exp(-mvPhaseBoundaryQg / (R_GAS * mvT))); }
    double GetGrainBoundaryEnergy() const { return (mvGrainBoundaryEnergy); }
    double GetPhaseBoundaryEnergy() const { return (mvPhaseBoundaryEnergy); }

// Helper functions
    double CarbonWtPercentToAtFraction(const double XcInWt);
    double CarbonAtFractionToWtPercent(const double XcInAtFraction);

// partitioning for unknown chemical potentials
     void CalculateXcEqNextToInterphase(const double T, double InitialXcAlphaInAt, double InitialXcGammaInAt, double *XcAlphaEq, double *XcGammaEq);

// Diffusivities of solute / carbon based on temperature (internally stored always, mvT) and possibly other state variables
    double GetCDiffusivityAgren(const double xC) const;
    double GetCDiffusivityMartensite() const;
    double GetCDiffusivityAgrenWithRestrictionForXC(const double xC, const double MaxAllowedXcToAffectDiffusivity) const;
    double GetSoluteDiffusivityFromUser() const;

    void LoadChemicalPotentialsAndLocalEquilibriumXC(const char *Thermofilename, bool KeepConstantIronAtomsPerCell);
    void MakeTablesYForXCEqVSCarbonInterface(const double* XCLocalEqFCC_AsRead, const double* XCLocalEqBCC_AsRead);
    void MakeTablesYForGibbsVSCarbon(const double* GibbsFCC_AsRead, const double* GibbsBCC_AsRead, const double* MuSubstitutionalFCC_AsRead, const double* MuSubstitutionalBCC_AsRead);
    void GetEqMuCarbonRelatedParametersABCD(double *A, double *B, double *C, double *D);
    double GetMuSubstitutionalInBCC(const double xC) const;
    double GetMuSubstitutionalInFCC(const double xC) const;
    double GetEqXcFCCAtThisInterface(const double xCBCCNow, const double xCFCCNow) const;
    double GetValueLocalEqXC_Clamped(const double X, const double* Y) const;
    double GetValueGibbs_Clamped(const double X, const double* Y) const;

    

/** Returns the Molar Volume for iron BCC
 * TODO: read the lattice data 
 */
double FerriteMolarVolume(){ return (7.1052e-6);}

/** Returns the Molar Volume for iron FCC
 *  * TODO: read the lattice data 
 */
double AusteniteMolarVolume(){   return (0.0000073713716);}


private:
    // Member variables
    double mvT; // temperature that SHOULD BE RESET/CHANGED IF DESIRED at every simulation step according to applied temperature time profile
    double mvMuCEqParamA, mvMuCEqParamB, mvMuCEqParamC, mvMuCEqParamD;
    double mvXC0;
    double mvSoluteDiffQg,mvSoluteDiffPreFactor;
    double mvGrainBoundaryM0, mvGrainBoundaryQg, mvGrainBoundaryEnergy;
    double mvPhaseBoundaryM0, mvPhaseBoundaryQg, mvPhaseBoundaryEnergy;
    double* mvpGibbsFCC_TableXC = nullptr;
    double* mvpGibbsBCC_TableXC = nullptr;
    double* mvpMuSubstitutionalFCC_TableXC = nullptr;
    double* mvpMuSubstitutionalBCC_TableXC = nullptr;
    double* mvpXCLocalEqFCC_TableXC = nullptr;
    double* mvpXCLocalEqBCC_TableXC = nullptr;
    double mvXCmaxInChemTables, mvXCminInChemTables, mvXCStepForChemTables;
    int mvNumberOfStepsForChemTables;
    double mvXCmaxInLocalEqRelationship, mvXCminInLocalEqRelationship, mvXCStepForLocalEqRelationship;
    int mvNumberOfStepsForLocalEqRelationship;
};

#endif // THERMODYNAMICS_H