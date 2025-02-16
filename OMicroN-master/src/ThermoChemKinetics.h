/* ThermoChemKinetics.h
   OMicroN (optimising microstructures numerically) simulation program
   Header file containing the ThermChemKin class definitions and implementation
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

/// @name Constants
/// @{
#define K_BOLTZMANN 1.3806504e-23 ///< Boltzmann constant (J/K)
#define R_GAS 8.314472            ///< Gas constant (m^3 Pa K^-1 mol^-1)
/// @}

/**
 * @class ThermChemKin
 * @brief Manages all the physical quantities that affect thermodynamics and kinetics on the microstructure simulation.
 * Stores all thermodynamic, chemical, kinetic properties in relation to lattice, temperature, solute concentration, defect density, etc..
 * These affect all simulated processes recrystallization, grain growth), solute redistribution (partitioning, diffusion, trapping to defects), phase transformations.
 */

class ThermChemKin
{
public:
    /**
     * @brief Constructor for ThermChemKin class.
     *
     * Initializes the ThermChemKin object with user-defined parameters.
     *
     * @param userXCo Initial carbon concentration (atomic fraction).
     * @param userStartTemperature Initial temperature (K).
     * @param userGB_Mo Pre-exponential factor for grain boundary migration.
     * @param userGB_Qg Activation energy for grain boundary migration (J/mol).
     * @param userGB_E Grain boundary energy density (J/m^2).
     * @param userPB_Mo Pre-exponential factor for phase boundary migration.
     * @param userPB_Qg Activation energy for phase boundary migration (J/mol).
     * @param userPB_E Phase boundary energy density (J/m^2).
     * @param user_DiffPreFactor Pre-exponential factor for solute diffusion.
     * @param user_DiffQg Activation energy for solute diffusion (J/mol).
     */
    ThermChemKin(double userXCo, double userStartTemperature, double userGB_Mo, double userGB_Qg, double userGB_E, double userPB_Mo, double userPB_Qg, double userPB_E, double user_DiffPreFactor, double user_DiffQg, double userXCu, double userXSn);

    /// Destructor
    ~ThermChemKin();

    /**
     * @brief Sets the temperature for the current simulation step.
     *
     * @param temp Temperature (K).
     */
    void SetTemperatureForThisStep(double temp) { mvT = temp; }

    /**
     * @brief Gets the temperature for the current simulation step.
     *
     * @return Temperature (K).
     */
    double GetTemperatureForThisStep() const { return mvT; }

    /**
     * @brief Calculates the grain boundary mobility.
     *
     * @return Grain boundary mobility.
     */
    double GetGrainBoundaryMobility() const { return (mvGrainBoundaryM0 * std::exp(-mvGrainBoundaryQg / (R_GAS * mvT))); }

    /**
     * @brief Calculates the phase boundary mobility.
     *
     * @return Phase boundary mobility.
     */
    double GetPhaseBoundaryMobility() const { return (mvPhaseBoundaryM0 * std::exp(-mvPhaseBoundaryQg / (R_GAS * mvT))); }

    /**
     * @brief Gets the grain boundary energy density.
     *
     * @return Grain boundary energy density (J/m^2).
     */
    double GetGrainBoundaryEnergy() const { return mvGrainBoundaryEnergy; }

    /**
     * @brief Gets the phase boundary energy density.
     *
     * @return Phase boundary energy density (J/m^2).
     */
    double GetPhaseBoundaryEnergy() const { return mvPhaseBoundaryEnergy; }

    /// @name Helper Functions
    /// @{

    /**
     * @brief Converts carbon weight percent to atomic fraction.
     *
     * @param XcInWt Carbon weight percent.
     * @return Carbon atomic fraction.
     */
    double CarbonWtPercentToAtFraction(const double XcInWt);

    /**
     * @brief Converts carbon atomic fraction to weight percent.
     *
     * @param XcInAtFraction Carbon atomic fraction.
     * @return Carbon weight percent.
     */
    double CarbonAtFractionToWtPercent(const double XcInAtFraction);

    /**
     * @brief Calculates the equilibrium carbon concentration next to the interphase.
     *
     * @param T Temperature (K).
     * @param InitialXcAlphaInAt Initial carbon concentration in alpha phase (atomic fraction).
     * @param InitialXcGammaInAt Initial carbon concentration in gamma phase (atomic fraction).
     * @param XcAlphaEq Pointer to store equilibrium carbon concentration in alpha phase.
     * @param XcGammaEq Pointer to store equilibrium carbon concentration in gamma phase.
     */
    void CalculateXcEqNextToInterphase(const double T, double InitialXcAlphaInAt, double InitialXcGammaInAt, double *XcAlphaEq, double *XcGammaEq);

    /**
     * @brief Gets the carbon diffusivity based on temperature using the Agren model.
     *
     * @param xC Carbon atomic fraction.
     * @return Carbon diffusivity.
     */
    double GetCDiffusivityAgren(const double xC) const;

    /**
     * @brief Gets the carbon diffusivity for martensite.
     *
     * @return Carbon diffusivity.
     */
    double GetCDiffusivityMartensite() const;

    /**
     * @brief Gets the carbon diffusivity with restriction for maximum allowed carbon concentration.
     *
     * @param xC Carbon atomic fraction.
     * @param MaxAllowedXcToAffectDiffusivity Maximum allowed carbon atomic fraction to affect diffusivity.
     * @return Carbon diffusivity.
     */
    double GetCDiffusivityAgrenWithRestrictionForXC(const double xC, const double MaxAllowedXcToAffectDiffusivity) const;

    /**
     * @brief Gets the solute diffusivity from user-defined parameters.
     *
     * @return Solute diffusivity.
     */
    double GetSoluteDiffusivityFromUser() const;

    /// @name Chemical Potentials and Thermodynamic Equilibrium
    /// @{

    /**
     * @brief Loads chemical potentials and local equilibrium data from a file.
     *
     * @param Thermofilename Filename containing the thermodynamic data.
     * @param KeepConstantIronAtomsPerCell Boolean flag to keep constant iron atoms per cell.
     */
    void LoadChemicalPotentialsAndLocalEquilibriumXC(const char *Thermofilename, bool KeepConstantIronAtomsPerCell);

    /**
     * @brief Creates tables for Y vs. equilibrium carbon concentration at the interface.
     *
     * @param XCLocalEqFCC_AsRead Pointer to data for FCC.
     * @param XCLocalEqBCC_AsRead Pointer to data for BCC.
     */
    void MakeTablesYForXCEqVSCarbonInterface(const double *XCLocalEqFCC_AsRead, const double *XCLocalEqBCC_AsRead);

    /**
     * @brief Creates tables for Gibbs energy vs. carbon concentration.
     *
     * @param GibbsFCC_AsRead Pointer to Gibbs energy data for FCC.
     * @param GibbsBCC_AsRead Pointer to Gibbs energy data for BCC.
     * @param MuSubstitutionalFCC_AsRead Pointer to chemical potential data for FCC.
     * @param MuSubstitutionalBCC_AsRead Pointer to chemical potential data for BCC.
     */
    void MakeTablesYForGibbsVSCarbon(const double *GibbsFCC_AsRead, const double *GibbsBCC_AsRead, const double *MuSubstitutionalFCC_AsRead, const double *MuSubstitutionalBCC_AsRead);

    /**
     * @brief Gets the parameters A, B, C, D for equilibrium carbon concentration.
     *
     * @param A Pointer to store parameter A.
     * @param B Pointer to store parameter B.
     * @param C Pointer to store parameter C.
     * @param D Pointer to store parameter D.
     */
    void GetEqMuCarbonRelatedParametersABCD(double *A, double *B, double *C, double *D);

    /**
     * @brief Gets the substitutional chemical potential in BCC.
     *
     * @param xC Carbon atomic fraction.
     * @return Substitutional chemical potential (J/mol).
     */
    double GetMuSubstitutionalInBCC(const double xC) const;

    /**
     * @brief Gets the substitutional chemical potential in FCC.
     *
     * @param xC Carbon atomic fraction.
     * @return Substitutional chemical potential (J/mol).
     */
    double GetMuSubstitutionalInFCC(const double xC) const;

    /**
     * @brief Gets the equilibrium carbon concentration in FCC at the interface.
     *
     * @param xCBCCNow Carbon concentration in BCC (atomic fraction).
     * @param xCFCCNow Carbon concentration in FCC (atomic fraction).
     * @return Equilibrium carbon concentration in FCC (atomic fraction).
     */
    double GetEqXcFCCAtThisInterface(const double xCBCCNow, const double xCFCCNow) const;

    /**
     * @brief Gets the local equilibrium carbon concentration, clamped to table limits.
     *
     * @param X Carbon concentration.
     * @param Y Pointer to the data table.
     * @return Local equilibrium carbon concentration.
     */
    double GetValueLocalEqXC_Clamped(const double X, const double *Y) const;

    /**
     * @brief Gets the Gibbs energy, clamped to table limits.
     *
     * @param X Carbon concentration.
     * @param Y Pointer to the data table.
     * @return Gibbs energy.
     */
    double GetValueGibbs_Clamped(const double X, const double *Y) const;

    /**
     * @brief Returns the molar volume for iron BCC.
     *
     * @return Molar volume for iron BCC (m^3/mol).
     */
    double FerriteMolarVolume() const { return 7.1052e-6; }

    /**
     * @brief Returns the molar volume for iron FCC.
     *
     * @return Molar volume for iron FCC (m^3/mol).
     */
    double AusteniteMolarVolume() const { return 0.0000073713716; }


/** @brief Returns the equilibrium concentration at liquid based on partitioning coefficient as a function of carbon
     *  @param xSolute the current (before further partitioning) solute concentration at interface
     *  @param xCarbon the current (before further partitioning) carbon concentration at interface
     *  @return k * xSolute
 */
double GetEqSoluteAtLiquidAtThisInterface(const double xSolute, const double xCarbon) const;
 double GetAverageCopper(){return mvXCu;}
 double GetAverageTin(){return mvXSn;}

private:
    /**
     * @brief Temperature for the simulation step.
     *
     * This variable stores the temperature (in Kelvin) that is used in calculations
     * for the current simulation step. It should be updated as per the applied temperature profile.
     */
    double mvT;

    /**
     * @brief Parameters A, B, C, D for partitioning (normally unused).
     *
     * These parameters are used for partitioning calculations, if tables are not provided.
     * - A: Parameter A for partitioning.
     * - B: Parameter B for partitioning.
     * - C: Parameter C for partitioning.
     * - D: Parameter D for partitioning.
     */
    double mvMuCEqParamA, mvMuCEqParamB, mvMuCEqParamC, mvMuCEqParamD;

    /**
     * @brief Average carbon concentration (atomic fraction).
     *
     * This variable represents the average carbon concentration in the material,
     * expressed as an atomic fraction.
     */
    double mvXC0;

    /**
     * @brief Average copper concentration (atomic fraction).
     *
     * This variable represents the average copper concentration in the material,
     * expressed as an atomic fraction.
     */
    double mvXCu;

        /**
     * @brief Average tin concentration (atomic fraction).
     *
     * This variable represents the average tin concentration in the material,
     * expressed as an atomic fraction.
     */
    double mvXSn;

    /**
     * @brief Activation energy for solute diffusion (J/mol).
     *
     * This variable holds the activation energy required for solute diffusion
     * within the material, measured in Joules per mole.
     */
    double mvSoluteDiffQg;

    /**
     * @brief Pre-exponential factor for solute diffusion.
     *
     * This variable represents the pre-exponential factor used in the Arrhenius equation
     * for solute diffusion. It is a coefficient that scales the diffusion rate.
     */
    double mvSoluteDiffPreFactor;

    /**
     * @brief Pre-exponential factor for grain boundary migration.
     *
     * This variable denotes the pre-exponential factor used in the Arrhenius equation
     * for grain boundary migration.
     */
    double mvGrainBoundaryM0;

    /**
     * @brief Activation energy for grain boundary migration (J/mol).
     *
     * This variable contains the activation energy required for grain boundary migration,
     * measured in Joules per mole.
     */
    double mvGrainBoundaryQg;

    /**
     * @brief Grain boundary energy density (J/m^2).
     *
     * This variable holds the energy density of grain boundaries, expressed in Joules per
     * square meter. It quantifies the energy associated with grain boundaries.
     */
    double mvGrainBoundaryEnergy;

    /**
     * @brief Pre-exponential factor for phase boundary migration.
     *
     * This variable represents the pre-exponential factor used in the Arrhenius equation
     * for phase boundary migration.
     */
    double mvPhaseBoundaryM0;

    /**
     * @brief Activation energy for phase boundary migration (J/mol).
     *
     * This variable holds the activation energy required for phase boundary migration,
     * measured in Joules per mole.
     */
    double mvPhaseBoundaryQg;

    /**
     * @brief Phase boundary energy density (J/m^2).
     *
     * This variable contains the energy density of phase boundaries, expressed in Joules per
     * square meter. It quantifies the energy associated with phase boundaries.
     */
    double mvPhaseBoundaryEnergy;

    /**
     * @brief Pointer to table containing Gibbs energy (J/mol) in FCC vs. carbon concentration (atomic fraction).
     *
     * This pointer points to a table that provides the Gibbs energy values for FCC as a function
     * of carbon concentration, expressed in atomic fraction. The table is used for thermodynamic calculations.
     */
    double *mvpGibbsFCC_TableXC = nullptr;

    /**
     * @brief Pointer to table containing Gibbs energy (J/mol) in BCC vs. carbon concentration (atomic fraction).
     *
     * This pointer points to a table that provides the Gibbs energy values for BCC as a function
     * of carbon concentration, expressed in atomic fraction. The table is used for thermodynamic calculations.
     */
    double *mvpGibbsBCC_TableXC = nullptr;

    /**
     * @brief Pointer to table containing chemical potential of substitutional atoms (J/mol) in FCC vs. carbon concentration (atomic fraction).
     *
     * This pointer points to a table that provides the chemical potential of substitutional atoms
     * in FCC as a function of carbon concentration, expressed in atomic fraction. The table is used
     * for calculating chemical potentials.
     */
    double *mvpMuSubstitutionalFCC_TableXC = nullptr;

    /**
     * @brief Pointer to table containing chemical potential of substitutional atoms (J/mol) in BCC vs. carbon concentration (atomic fraction).
     *
     * This pointer points to a table that provides the chemical potential of substitutional atoms
     * in BCC as a function of carbon concentration, expressed in atomic fraction. The table is used
     * for calculating chemical potentials.
     */
    double *mvpMuSubstitutionalBCC_TableXC = nullptr;

    /**
     * @brief Pointer to table containing local equilibrium carbon concentration (atomic fraction) in FCC.
     *
     * This pointer points to a table that provides the local equilibrium carbon concentration in FCC
     * as a function of the total interphase carbon concentration, expressed in atomic fraction.
     */
    double *mvpXCLocalEqFCC_TableXC = nullptr;

    /**
     * @brief Pointer to table containing local equilibrium carbon concentration (atomic fraction) in BCC.
     *
     * This pointer points to a table that provides the local equilibrium carbon concentration in BCC
     * as a function of the total interphase carbon concentration, expressed in atomic fraction.
     */
    double *mvpXCLocalEqBCC_TableXC = nullptr;

    /**
     * @brief Table limits and step size for chemical tables.
     *
     * These variables define the limits and step size for the chemical tables used in calculations:
     * - mvXCmaxInChemTables: Maximum carbon concentration in the chemical tables.
     * - mvXCminInChemTables: Minimum carbon concentration in the chemical tables.
     * - mvXCStepForChemTables: Step size between concentration values in the chemical tables.
     */
    double mvXCmaxInChemTables, mvXCminInChemTables, mvXCStepForChemTables;

    /**
     * @brief Number of steps for chemical tables.
     *
     * This variable specifies the number of discrete steps used in the chemical tables for interpolation.
     */
    int mvNumberOfStepsForChemTables;

    /**
     * @brief Table limits and step size for local equilibrium relationship.
     *
     * These variables define the limits and step size for the local equilibrium relationship tables:
     * - mvXCmaxInLocalEqRelationship: Maximum carbon concentration in the local equilibrium relationship tables.
     * - mvXCminInLocalEqRelationship: Minimum carbon concentration in the local equilibrium relationship tables.
     * - mvXCStepForLocalEqRelationship: Step size between concentration values in the local equilibrium relationship tables.
     */
    double mvXCmaxInLocalEqRelationship, mvXCminInLocalEqRelationship, mvXCStepForLocalEqRelationship;

    /**
     * @brief Number of steps for local equilibrium relationship.
     *
     * This variable specifies the number of discrete steps used in the local equilibrium relationship tables for interpolation.
     */
    int mvNumberOfStepsForLocalEqRelationship;
};

#endif // ThermoChemKinetics_H
