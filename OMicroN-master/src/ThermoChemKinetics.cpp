/* ThermoChemKinetics.cpp
   OMicroN (optimising microstructures numerically) simulation program
   cpp-file containing ThermChemKin class implementation
*/

#include "ThermoChemKinetics.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "hdf5_utils.h"
#include <algorithm>
#include <H5Cpp.h>




ThermChemKin::ThermChemKin(double userStartTemperature, double userXCo,double userGB_Mo,double userGB_Qg,double userGB_E,double userPB_Mo,double userPB_Qg,double userPB_E, double user_DiffPreFactor, double user_DiffQg, double userXCu, double userXSn) {
    mvT = userStartTemperature; // only in initialization mvT is this temperature. Later it is set per simulation step here
    mvMuCEqParamA = -90885;
    mvMuCEqParamB = 61;
    mvMuCEqParamC = 9837;
    mvMuCEqParamD = -5.7;

    mvXC0 = userXCo; // acerage carbon concentration (at. fraction)
mvXCu = userXCu;
mvXSn = userXSn;

    mvGrainBoundaryM0 = userGB_Mo;
    mvGrainBoundaryQg = userGB_Qg;
    mvGrainBoundaryEnergy = userGB_E;

    mvPhaseBoundaryM0 = userPB_Mo;
    mvPhaseBoundaryQg = userPB_Qg;
    mvPhaseBoundaryEnergy = userPB_E;

    mvSoluteDiffPreFactor = user_DiffPreFactor;
    mvSoluteDiffQg = user_DiffQg;
    
}

ThermChemKin::~ThermChemKin() {
    delete[] mvpGibbsFCC_TableXC;
    delete[] mvpGibbsBCC_TableXC;
    delete[] mvpMuSubstitutionalFCC_TableXC;
    delete[] mvpMuSubstitutionalBCC_TableXC;
    delete[] mvpXCLocalEqFCC_TableXC;
    delete[] mvpXCLocalEqBCC_TableXC;
}


/** Returns the converted wt pct to at frac for carbon - iron
 *  *   @param XcInWt concentration of carbon in wt pct
 *  * TODO: read the lattice data 
 */
double ThermChemKin::CarbonWtPercentToAtFraction(const double XcInWt) {
    double Na = 6.02214076 * 1.e23;
    double CarbonAtoms = XcInWt * Na / 12.01;
    double FeAtoms = (100. - XcInWt) * Na / 55.845;
    CarbonAtoms *= 0.01;
    FeAtoms *= 0.01;
    double Tot = CarbonAtoms + FeAtoms;
    return CarbonAtoms / Tot;
}


/** Returns the converted at frac to wt pct for carbon - iron
 *   @param XcInAtFraction concentration of carbon in at fraction
 *  *  * TODO: read the lattice data 
 */
double ThermChemKin::CarbonAtFractionToWtPercent(const double XcInAtFraction) {
    double Na = 6.02214076 * 1.e23;
    double weightCarbon = 12.01 * XcInAtFraction * 100. / Na;
    double FractionFe = (1. - XcInAtFraction);
    double weightFe = 55.845 * 100. * FractionFe / Na;
    return (weightCarbon * 100.) / (weightCarbon + weightFe);
}
/** @brief The purpose of this function is to partition the solute ONLY IN CASE the chemical potentials are known
/** Calculates local equilibrium (partitioning) carbon concentrations based on analytical function for iron-carbon that is solved numerically
 *  *   @param T the current temperature
 *      @param InitialXcAlphaInAt the current carbon concentration in BCC next to interphase
 *      @param InitialXcGammaInAt the current carbon concentration in FCC next to interphase
 *      @return XcAlphaEq the partitioned (and local equilibrium) carbon concentration in BCC next to interphase
 *  *   @return XcGammaEq the partitioned (and local equilibrium) carbon concentration in FCC next to interphase
 *  * TODO: read the lattice data 
 */
void ThermChemKin::CalculateXcEqNextToInterphase(const double T, double InitialXcAlphaInAt, double InitialXcGammaInAt, double *XcAlphaEq, double *XcGammaEq) {
    double fAlpha = 0.5;
    double fGamma = 0.5;
    double mvT = T;
    double A = mvMuCEqParamA;
    double B = mvMuCEqParamB;
    double C = mvMuCEqParamC;
    double D = mvMuCEqParamD;
    double XcTot = fGamma * InitialXcGammaInAt + fAlpha * InitialXcAlphaInAt;
    double XcGammaTrial = InitialXcGammaInAt;
    double XcAlphaTrial = InitialXcAlphaInAt;
    double XcGammaTrialPre = XcGammaTrial;
    int iter;
    int maxIter = 100;
    double WeightForTol = 1.e-5;
    double relTol;

    for (iter = 0; iter < maxIter; iter++) {
        XcGammaTrial = XcGammaTrialPre;
        double InExp = ((A + B * mvT) + (C + D * mvT) * CarbonAtFractionToWtPercent(XcGammaTrial)) / (R_GAS * mvT);
        double f = (XcTot - fGamma * XcGammaTrial) / fAlpha - XcGammaTrial * std::exp(InExp);
        double fPrime = -fGamma / fAlpha - XcGammaTrial * std::exp(InExp) * (C + D * mvT) / (R_GAS * mvT) - std::exp(InExp);
        XcGammaTrial -= (f / fPrime);
        relTol = XcGammaTrial - XcGammaTrialPre;
        if (relTol < 0.) relTol = -relTol;
        if (relTol < XcGammaTrialPre * WeightForTol && iter > 10) break;
        XcGammaTrialPre = XcGammaTrial;
    }

    if (relTol >= XcGammaTrialPre * WeightForTol || iter > maxIter - 3) {
        // LOG_F(ERROR, "check here %g from original xcgammma %g and alpha i got relTol %g and XcGammaPre %g and iter %d",
        //       XcGammaTrial, InitialXcGammaInAt, relTol, XcGammaTrialPre, iter);
        std::cerr << "check here " << XcGammaTrial << " from original xcgammma " << InitialXcGammaInAt << " and alpha i got relTol " << relTol << " and XcGammaPre " << XcGammaTrialPre << " and iter " << iter << std::endl;
    }

    XcAlphaTrial = (XcTot - fGamma * XcGammaTrial) / fAlpha;

    if (iter == maxIter) {
        // LOG_F(ERROR, "problem iter %d with XcGammaTrial %g check InitialXcGamma wtpct %g InitialXcAlpha wtpct %g",
        //       iter, XcGammaTrial, CarbonAtFractionToWtPercent(InitialXcGammaInAt), CarbonAtFractionToWtPercent(InitialXcAlphaInAt));
        std::cerr << "problem iter " << iter << " with XcGammaTrial " << XcGammaTrial << " check InitialXcGamma wtpct " << CarbonAtFractionToWtPercent(InitialXcGammaInAt) << " InitialXcAlpha wtpct " << CarbonAtFractionToWtPercent(InitialXcAlphaInAt) << std::endl;
    }

    if (XcGammaTrial < 0. || XcAlphaTrial < 0.) {
        // LOG_F(ERROR, "check XcGamma negative %g XcAlphaTrial %g InitialXcGamma %g InitialXcAlpha %g",
        //       XcGammaTrial, XcAlphaTrial, InitialXcGammaInAt, InitialXcAlphaInAt);
        std::cerr << "check XcGamma negative " << XcGammaTrial << " XcAlphaTrial " << XcAlphaTrial << " InitialXcGamma " << InitialXcGammaInAt << " InitialXcAlpha " << InitialXcAlphaInAt << std::endl;
    }

    *XcGammaEq = XcGammaTrial;
    *XcAlphaEq = XcAlphaTrial;
}



       /** @brief Gets solute diffusivity using Arrhenius type equation at
     * temperature #mT given user defined input values 
     *
     * @return Diffusivity in \f$m^2 s^{-1}\f$
     */
    double ThermChemKin::GetSoluteDiffusivityFromUser() const 
    { 
        return (mvSoluteDiffPreFactor * exp(-mvSoluteDiffQg / (R_GAS * mvT)));
         }


/** @brief Gets Ågren's composition dependent austenite carbon diffusivity at
 * temperature #mVT
 * (J. Ågren, Scr. Metall. 20 (1986) 1507–1510)
 *
 * @param Carbon concentration in at. fraction
 *
 * @return Diffusivity in \f$m^2 s^{-1}\f$
 */
double ThermChemKin::GetCDiffusivityAgren(const double xC) const 
{
    double yC = xC / (1.0 - xC);
    double D0 = 4.53e-7 * (1.0 + yC * (1.0 - yC) * 8339.9 / mvT);
    return D0 * std::exp(-(1. / mvT - 2.221e-4) * (17767.0 - yC * 26436.0));
}


/** @brief Gets Ågren's  martensite austenite carbon diffusivity at
 * temperature #mVT
 * (J. Ågren, Journal of Physics and Chemistry of Solids (1982)) *
 * @return Diffusivity in \f$m^2 s^{-1}\f$
 */

double ThermChemKin::GetCDiffusivityMartensite() const {
    double D0 = 2.e-6;
    double FirstExp = std::exp(-10115. / mvT);
    double SecondExp = std::exp(0.5898 * (1. + 2. * atan(1.48985 - 15309. / mvT) / M_PI));
    return D0 * FirstExp * SecondExp;
}

    /** @brief This function is useful if there are tiny FCC grains that will enrich so much that diffusivity becomes very high and time steps reqired for diffusion are then very small.
     *  Gets carbon diffusivity using Agren's expression (for temperature mvT) but with a max allowed value. 
     *  @param xC the local carbon concentration
     *  @param MaxAllowedXcToAffectDiffusivity the max allowed value
     * @return Diffusivity in \f$m^2 s^{-1}\f$
     */
double ThermChemKin::GetCDiffusivityAgrenWithRestrictionForXC(const double xC, const double MaxAllowedXcToAffectDiffusivity) const {
    double yC = (xC < MaxAllowedXcToAffectDiffusivity) ? xC / (1.0 - xC) : MaxAllowedXcToAffectDiffusivity / (1.0 - MaxAllowedXcToAffectDiffusivity);
    double D0 = 4.53e-7 * (1.0 + yC * (1.0 - yC) * 8339.9 / mvT);
    return D0 * std::exp(-(1. / mvT - 2.221e-4) * (17767.0 - yC * 26436.0));
}


/** @brief Loads the hdf5 file that contains chemical potentials of substitutional lattice, the equilibrium (paritioning) solute concentration for various interphase compositions
     *  @param filename the name of the file with the data
     *  @param KeepConstantIronAtomsPerCell user defined flag signaling whether equilibrium requires constant at. fractions or substitional-to-total fraction 
 */

void ThermChemKin::LoadChemicalPotentialsAndLocalEquilibriumXC(const char* filename, bool KeepConstantIronAtomsPerCell) {
    H5_IO file(filename, "r");

    // Read datasets using H5_IO class
    std::vector<double> xC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/xC");
    std::vector<double> GibbsJoulePerMoleFCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/GibbsJoulePerMoleFCC");
    std::vector<double> GibbsJoulePerMoleBCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/GibbsJoulePerMoleBCC");
    std::vector<double> GibbsJoulePerMoleSubsFCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/GibbsJoulePerMoleSubsFCC");
    std::vector<double> GibbsJoulePerMoleSubsBCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/GibbsJoulePerMoleSubsBCC");
    std::vector<double> muCarbonFCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/muCarbonFCC");
    std::vector<double> muCarbonBCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/muCarbonBCC");
    std::vector<double> muFeFCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/muFeFCC");
    std::vector<double> muFeBCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/muFeBCC");
    std::vector<double> muSubsFCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/muSubsFCC");
    std::vector<double> muSubsBCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/muSubsBCC");

    std::vector<double> LocalXcPossible, LocalXcEqFCC, LocalXcEqBCC;
    if (KeepConstantIronAtomsPerCell) {
        LocalXcPossible = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/XcTotLocalModifiedMoleFracToCarbonPerIron");
        LocalXcEqFCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/XcEqLocalFCC_ModifiedMoleFracToCarbonPerIron");
        LocalXcEqBCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/XcEqLocalBCC_ModifiedMoleFracToCarbonPerIron");
    } else {
        LocalXcPossible = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/XcTotLocal");
        LocalXcEqFCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/XcEqLocalFCC");
        LocalXcEqBCC = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/XcEqLocalBCC");
    }

    for (int i = 0; i < 20; i++) {
        std::cout << "xcTotLocal: " << LocalXcPossible[i] << std::endl;
        std::cout << "xceqfcc: " << LocalXcEqFCC[i] << std::endl;
        std::cout << "xceqBcc: " << LocalXcEqBCC[i] << std::endl;
    }

    std::vector<double> A = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/EqMuCarbonParamA");
    std::vector<double> B = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/EqMuCarbonParamB");
    std::vector<double> C = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/EqMuCarbonParamC");
    std::vector<double> D = file.ReadDataset<double>("/tables/thermodynamicsPartitioning/EqMuCarbonParamD");
    mvMuCEqParamA = A[0];
    mvMuCEqParamB = B[0];
    mvMuCEqParamC = C[0];
    mvMuCEqParamD = D[0];

    std::cout << "Set muC parameters A = " << mvMuCEqParamA
              << " , B = " << mvMuCEqParamB
              << " , C = " << mvMuCEqParamC
              << " , D = " << mvMuCEqParamD << std::endl;

    double minCarbon = *std::min_element(xC.begin(), xC.end());
    double maxCarbon = *std::max_element(xC.begin(), xC.end());
    double minCarbonForLocalEq = *std::min_element(LocalXcPossible.begin(), LocalXcPossible.end());
    double maxCarbonForLocalEq = *std::max_element(LocalXcPossible.begin(), LocalXcPossible.end());

    mvXCmaxInChemTables = maxCarbon;
    mvXCminInChemTables = minCarbon;
    mvXCStepForChemTables = (mvXCmaxInChemTables - mvXCminInChemTables) / xC.size();
    mvNumberOfStepsForChemTables = xC.size();
    mvXCminInLocalEqRelationship = minCarbonForLocalEq;
    mvXCmaxInLocalEqRelationship = maxCarbonForLocalEq;
    mvXCStepForLocalEqRelationship = (mvXCmaxInLocalEqRelationship - mvXCminInLocalEqRelationship) / LocalXcPossible.size();
    mvNumberOfStepsForLocalEqRelationship = LocalXcPossible.size();

    MakeTablesYForXCEqVSCarbonInterface(LocalXcEqFCC.data(), LocalXcEqBCC.data());
    MakeTablesYForGibbsVSCarbon(GibbsJoulePerMoleFCC.data(), GibbsJoulePerMoleBCC.data(), muSubsFCC.data(), muSubsBCC.data());
}


/** @brief Makes tables for carbon local equilibrium due to partitioning
     *  @param XCLocalEqFCC_AsRead the FCC eq. concentration for given (AND EQUIDISTANT) total carbon at interphase
     *  @param XCLocalEqBCC_AsRead the BCC eq. concentration for given (AND EQUIDISTANT) total carbon at interphase
 */
void ThermChemKin::MakeTablesYForXCEqVSCarbonInterface(const double* XCLocalEqFCC_AsRead, const double* XCLocalEqBCC_AsRead) {
    mvpXCLocalEqFCC_TableXC = new double[mvNumberOfStepsForLocalEqRelationship];
    mvpXCLocalEqBCC_TableXC = new double[mvNumberOfStepsForLocalEqRelationship];
    for (int i = 0; i < mvNumberOfStepsForLocalEqRelationship; i++) {
        mvpXCLocalEqFCC_TableXC[i] = XCLocalEqFCC_AsRead[i];
        mvpXCLocalEqBCC_TableXC[i] = XCLocalEqBCC_AsRead[i];
    }
}

/** @brief Makes tables for carbon local equilibrium due to partitioning
     *  @param GibbsFCC_AsRead the Gibbs energy in FCC (J per mole) given (AND EQUIDISTANT) total carbon at interphase
     *  @param GibbsBCC_AsRead the Gibbs energy in BCC (J per mole) given (AND EQUIDISTANT) total carbon at interphase
     *  @param MuSubstitutionalFCC_AsRead the chemical potential of substitutional atoms in FCC (J per mole) given the (EQUIDISTANT) carbon concentration
     *  @param MuSubstitutionalBCC_AsRead the chemical potential of substitutional atoms in BCC (J per mole) given (EQUIDISTANT) carbon concentration

 */
void ThermChemKin::MakeTablesYForGibbsVSCarbon(const double* GibbsFCC_AsRead, const double* GibbsBCC_AsRead, const double* MuSubstitutionalFCC_AsRead, const double* MuSubstitutionalBCC_AsRead) {
    mvpGibbsFCC_TableXC = new double[mvNumberOfStepsForChemTables];
    mvpGibbsBCC_TableXC = new double[mvNumberOfStepsForChemTables];
    mvpMuSubstitutionalFCC_TableXC = new double[mvNumberOfStepsForChemTables];
    mvpMuSubstitutionalBCC_TableXC = new double[mvNumberOfStepsForChemTables];
    for (int i = 0; i < mvNumberOfStepsForChemTables; i++) {
        mvpGibbsFCC_TableXC[i] = GibbsFCC_AsRead[i];
        mvpGibbsBCC_TableXC[i] = GibbsBCC_AsRead[i];
        mvpMuSubstitutionalFCC_TableXC[i] = MuSubstitutionalFCC_AsRead[i];
        mvpMuSubstitutionalBCC_TableXC[i] = MuSubstitutionalBCC_AsRead[i];
    }
}

void ThermChemKin::GetEqMuCarbonRelatedParametersABCD(double *A, double *B, double *C, double *D) {
    *A = mvMuCEqParamA;
    *B = mvMuCEqParamB;
    *C = mvMuCEqParamC;
    *D = mvMuCEqParamD;
}


/** @brief Returns the linearly interpolated value of chemical potential of substitutional atoms in BCC given the carbon concentration
     *  @param xC the carbon concentration
 */
double ThermChemKin::GetMuSubstitutionalInBCC(const double xC) const {
    // CHECK_F(mvpMuSubstitutionalBCC_TableXC != nullptr, "Gibbs energy in BCC interpolator has not been initialized, check if files are loaded first of all");
    if (mvpMuSubstitutionalBCC_TableXC == nullptr) {
        std::cerr << "Gibbs energy in BCC interpolator has not been initialized, check if files are loaded first of all" << std::endl;
        return 0.0; // or some appropriate error value
    }
    if (xC > mvXCmaxInChemTables) {
        return GetValueGibbs_Clamped(mvXCmaxInChemTables, mvpMuSubstitutionalBCC_TableXC);
    }
    if (xC < mvXCminInChemTables) {
        return GetValueGibbs_Clamped(mvXCminInChemTables, mvpMuSubstitutionalBCC_TableXC);
    }
    return GetValueGibbs_Clamped(xC, mvpMuSubstitutionalBCC_TableXC);
}

/** @brief Returns the linearly interpolated value of chemical potential of substitutional atoms in FCC given the carbon concentration
     *  @param xC the carbon concentration
 */
double ThermChemKin::GetMuSubstitutionalInFCC(const double xC) const {
    // CHECK_F(mvpMuSubstitutionalFCC_TableXC != nullptr, "Gibbs energy in FCC interpolator has not been initialized, check if files are loaded first of all");
    if (mvpMuSubstitutionalFCC_TableXC == nullptr) {
        std::cerr << "Gibbs energy in FCC interpolator has not been initialized, check if files are loaded first of all" << std::endl;
        return 0.0; // or some appropriate error value
    }
    if (xC > mvXCmaxInChemTables) {
        return GetValueGibbs_Clamped(mvXCmaxInChemTables, mvpMuSubstitutionalFCC_TableXC);
    }
    if (xC < mvXCminInChemTables) {
        return GetValueGibbs_Clamped(mvXCminInChemTables, mvpMuSubstitutionalFCC_TableXC);
    }
    return GetValueGibbs_Clamped(xC, mvpMuSubstitutionalFCC_TableXC);
}




/** @brief Returns the linearly interpolated value equilibrium (local partitioning) carbon concentration in FCC given the total interphase carbon concentration
     *  @param xCBCCNow the current (before further partitioning) carbon concentration in adjacent BCC
     *  @param xCFCCNow the current (before further partitioning) carbon concentration in adjacent FCC
     *  @return xCFCC linearly interpolated value equilibrium (local partitioning) carbon concentration in FCC
 */
double ThermChemKin::GetEqXcFCCAtThisInterface(const double xCBCCNow, const double xCFCCNow) const {
    double xCFCC;
    // CHECK_F(mvXCLocalEqFCC_TableXC != nullptr, "XC Eq. Local in FCC interpolator has not been initialized, check if files are loaded first of all");
    if (mvpXCLocalEqFCC_TableXC == nullptr) {
        std::cerr << "XC Eq. Local in FCC interpolator has not been initialized, check if files are loaded first of all" << std::endl;
        return 0.0; // or some appropriate error value
    }

    double xC = 0.5 * (xCBCCNow + xCFCCNow);

    xCFCC = GetValueLocalEqXC_Clamped(xC, mvpXCLocalEqFCC_TableXC);
    if (xC > mvXCmaxInLocalEqRelationship) {
        xCFCC = GetValueLocalEqXC_Clamped(mvXCmaxInLocalEqRelationship, mvpXCLocalEqFCC_TableXC);
    }
    if (xC < mvXCminInLocalEqRelationship) {
        xCFCC = GetValueLocalEqXC_Clamped(mvXCminInLocalEqRelationship, mvpXCLocalEqFCC_TableXC);
    }
    if (2. * xC - xCFCC < 1.e-10) {
        double xCBCC = GetValueLocalEqXC_Clamped(xC, mvpXCLocalEqBCC_TableXC);
        xCFCC = 2. * xC - xCBCC;
        if (xCBCC < 1.e-10 || xCFCC < 1.e-10) return xCFCCNow;
    }

    return xCFCC;
}



double ThermChemKin::GetValueLocalEqXC_Clamped(const double X, const double* Y) const {
    int li = static_cast<int>((X - mvXCminInLocalEqRelationship) / mvXCStepForLocalEqRelationship);
    if (li < 0) {
        return Y[0];
    }
    if (li >= mvNumberOfStepsForLocalEqRelationship) {
        return Y[mvNumberOfStepsForLocalEqRelationship - 1];
    }
    if (li == mvNumberOfStepsForLocalEqRelationship - 1) {
        return Y[li];
    } else {
        return Y[li] + (Y[li + 1] - Y[li]) * (X - mvXCminInLocalEqRelationship - li * mvXCStepForLocalEqRelationship) / mvXCStepForLocalEqRelationship;
    }
}

double ThermChemKin::GetValueGibbs_Clamped(const double X, const double* Y) const {
    int li = static_cast<int>((X - mvXCminInChemTables) / mvXCStepForChemTables);
    if (li < 0) {
        return Y[0];
    }
    if (li >= mvNumberOfStepsForChemTables) {
        return Y[mvNumberOfStepsForChemTables - 1];
    }
    if (li == mvNumberOfStepsForChemTables - 1) {
        return Y[li];
    } else {
        return Y[li] + (Y[li + 1] - Y[li]) * (X - mvXCminInChemTables - li * mvXCStepForChemTables) / mvXCStepForChemTables;
    }
}


/** @brief Returns the equilibrium concentration at liquid based on partitioning coefficient as a function of carbon
     *  @param xSolute the current (before further partitioning) solute concentration at interface
     *  @param xCarbon the current (before further partitioning) carbon concentration at interface
     *  @return k * xSolute
 */
double ThermChemKin::GetEqSoluteAtLiquidAtThisInterface(const double xSolute, const double xCarbon) const {
    double k = 1.1;
    return k * xSolute;
}