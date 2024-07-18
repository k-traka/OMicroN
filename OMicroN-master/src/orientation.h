/** @file src/orientation.h
 *
   OMicroN (optimising microstructures numerically) simulation program
 * header file containing the orientations class definitions and implementation
 */

#ifndef Orientation_H
#define Orientation_H

#include "settings.h"
#include <array>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <iostream>
#include <vector>

/// Convert rotation in the angle-axis representation to quaternion.
#define AA2Q(A, X, Y, Z) Eigen::Quaternionf(Eigen::AngleAxisf(A, Eigen::Vector3f(X, Y, Z) / sqrt(X * X + Y * Y + Z * Z))).normalized()

/// @{
const float MISORIENTATION_MAX = 64.; ///< Maximum misorientation angle. 63.8 deg for cubic symmetry
/// @}

/// @name Symmetry and orientation relationship operators
/// @{
const int N_SYM = 24; ///< Number of symmetry operators for the cubic system

/// Cubic symmetry operators (0-23)
const Eigen::Quaternionf SYM[N_SYM] = {
    AA2Q(0., 1., 0., 0.),
    AA2Q(M_PI, 1., 0., 0.),
    AA2Q(M_PI, 0., 1., 0.),
    AA2Q(M_PI, 0., 0., 1.),
    AA2Q(M_PI / 2., 1., 0., 0.),
    AA2Q(M_PI / 2., 0., 1., 0.),
    AA2Q(M_PI / 2., 0., 0., 1.),
    AA2Q(M_PI / 2., -1., 0., 0.),
    AA2Q(M_PI / 2., 0., -1., 0.),
    AA2Q(M_PI / 2., 0., 0., -1.),
    AA2Q(M_PI, 1., 1., 0.),
    AA2Q(M_PI, 1., 0., 1.),
    AA2Q(M_PI, 0., 1., 1.),
    AA2Q(M_PI, 1., -1., 0.),
    AA2Q(M_PI, -1., 0., 1.),
    AA2Q(M_PI, 0., 1., -1.),
    AA2Q(M_PI * 2. / 3., 1., 1., 1.),
    AA2Q(M_PI * 2. / 3., 1., -1., 1.),
    AA2Q(M_PI * 2. / 3., -1., 1., 1.),
    AA2Q(M_PI * 2. / 3., -1., -1., 1.),
    AA2Q(M_PI * 2. / 3., 1., 1., -1.),
    AA2Q(M_PI * 2. / 3., 1., -1., -1.),
    AA2Q(M_PI * 2. / 3., -1., 1., -1.),
    AA2Q(M_PI * 2. / 3., -1., -1., -1.)
};

const int N_CSL19a = N_SYM; ///< Number of CSL19a equivalent representations (i.e. 27 degrees around <110>) orientation relationships

/** Symmetrically equivalent operators to apply the CSL19a (Ibe and Lucke) 
 */
const Eigen::Quaternionf CSL19a[N_CSL19a] = {
    AA2Q(0.463, -0.7071, 0, 0.7071),
    CSL19a[0] * SYM[1],
    CSL19a[0] * SYM[2],
    CSL19a[0] * SYM[3],
    CSL19a[0] * SYM[4],
    CSL19a[0] * SYM[5],
    CSL19a[0] * SYM[6],
    CSL19a[0] * SYM[7],
    CSL19a[0] * SYM[8],
    CSL19a[0] * SYM[9],
    CSL19a[0] * SYM[10],
    CSL19a[0] * SYM[11],
    CSL19a[0] * SYM[12],
    CSL19a[0] * SYM[13],
    CSL19a[0] * SYM[14],
    CSL19a[0] * SYM[15],
    CSL19a[0] * SYM[16],
    CSL19a[0] * SYM[17],
    CSL19a[0] * SYM[18],
    CSL19a[0] * SYM[19],
    CSL19a[0] * SYM[20],
    CSL19a[0] * SYM[21],
    CSL19a[0] * SYM[22],
    CSL19a[0] * SYM[23]
};

/// @name Limits of Euler space for triclinic and orthorhombic sample symmetry
/// @{
const float LIM_TC[3] = { 360., 90., 90. }; ///< Triclinic symmetry {phi1, Phi, phi2}
const float LIM_OR[3] = { 90., 90., 90. }; ///< Orthorhombic symmetry {phi1, Phi, phi2}
/// @}

/**
 * @class Orientations
 * @brief Class to handle orientations and symmetry operations.
 *
 * This class manages the orientations and crystal symmetries of materials
 * in the OMicroN simulation. It provides functionalities to add and retrieve
 * orientations, convert between Euler angles and quaternions, and calculate
 * misorientations and other orientation-related metrics.
 * IN THIS CLASS ORIENTATIONS ARE STORED (PROVIDED/EXPORTED) AS EULER ANGLES BUT ALL OPERATIONS TAKE PLACE WITH QUATERNIONS
 */
class Orientations {
private:
    /// @name Geometry
    /// @{
    const float* mvpLimits; ///< Limits of (convenient zone of) Euler space
    // float mvVolumeUnit; ///< Volume unit (i.e. cell volume)
    /// @}

    /// @name Grain boundaries
    /// @{
    float mvLAGB; ///< Low angle grain boundary minimum misorientation
    float mvHAGB; ///< High angle grain boundary minimum misorientation
    /// @}

    /// @name Simulation data
    /// @{
    std::vector<Eigen::Vector3f> mvOrientation; ///< Orientations (euler angles) of each phase
    std::vector<int> mvOrientationId; ///< Orientations (euler angles) of each phase
    std::vector<int> mvCrystalSymmetry; ///< CrystalSymmetry for this ori Id (for now only cubic)
    std::vector<int> mvCellsInOrientationId; ///< How many cells have this orientation (in case euler angles are passed as grains average)
    /// @}

public:
    /** @brief Default constructor for Orientations class. */
    Orientations(void);
    //~Orientations(void); ///< Destructor

    /**
     * @brief Sets the orientation parameters from user settings.
     * 
     * @param us UserSettings object containing the parameters.
     */
    void SetOrientationParameters(const UserSettings& us)
    { 
        mvpLimits = LIM_TC;
        mvLAGB = us.mvLowerMisorientationCutOff;
        mvHAGB = us.mvHAGB;
    }

    /**
     * @brief Gets the number of orientations.
     * 
     * @return Number of orientations.
     */
    int GetOrientationCount() const {
        return mvOrientation.size();
    }

    /**
     * @brief Retrieves the Euler angles for a given orientation ID.
     * 
     * @param id Orientation ID.
     * @return Euler angles as an Eigen::Vector3f.
     */
    Eigen::Vector3f GetEulerAngles(const int id) {
        if (id < 0 || id >= mvOrientation.size()) {
            std::cerr << "Error: Orientation ID out of bounds: " << id << std::endl;
            return Eigen::Vector3f::Zero(); // Return a default value or handle the error appropriately
        }
        return mvOrientation[id];
    }

    /**
     * @brief Prints all orientations for debugging purposes.
     */
    void DebugPrintOrientations() const {
        for (size_t i = 0; i < mvOrientation.size(); ++i) {
            std::cout << "Orientation " << i << ": " << mvOrientation[i].transpose() << std::endl;
        }
    }

    /// @name Euler grid and problem size
    /// @{
    /**
     * @brief Returns the limit of Euler space for a given component.
     * 
     * @param i Component index (0:phi1, 1:Phi, 2:phi2).
     * @return Limit of Euler space for the component.
     */
    float Lim(const int i) const { return mvpLimits[i]; }

    /**
     * @brief Adds an orientation and returns its ID.
     * 
     * @param ip Crystal symmetry index.
     * @param ea Euler angles.
     * @return ID of the added orientation.
     */
    int addOriAndReturnId(const int ip, const Eigen::Vector3f& ea);


    /**
     * @brief Calculates the misorientation between two quaternions.
     * 
     * @param q0 First quaternion.
     * @param q1 Second quaternion.
     * @return Misorientation angle.
     */
    float misorientationQ(const Eigen::Quaternionf q0, const Eigen::Quaternionf q1);

    /**
     * @brief Calculates the distance of the boundary from CSL19a.
     * 
     * @param ea1 Euler angles of the first phase.
     * @param ea2 Euler angles of the second phase.
     * @return Misorientation angle.
     */
    float DistanceOfBoundaryFromCSL19a(const Eigen::Vector3f& ea1, const Eigen::Vector3f& ea2);

    /**
     * @brief Calculates the misorientation between two orientations given their IDs.
     * 
     * @param Id1 ID of the first orientation.
     * @param Id2 ID of the second orientation.
     * @return Misorientation angle.
     */
    float CalculateMisorientationBetweenTwoOriIds(const int Id1, const int Id2);

    /**
     * @brief Calculates the misorientation between two orientations given their Euler angles.
     * 
     * @param ea1 Euler angles of the first orientation.
     * @param ea2 Euler angles of the second orientation.
     * @return Misorientation angle.
     */
    float CalcMisorientationFromEulerAngles(const Eigen::Vector3f& ea1, const Eigen::Vector3f& ea2);
};

#endif // ORIENTATION_H
