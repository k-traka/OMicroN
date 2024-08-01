/** @file src/orientation.h
 *
 * OMicroN (Optimizing Microstructures Numerically) simulation program.
 * Header file containing the Orientations class definitions and implementation.
 */

#ifndef ORIENTATION_H
#define ORIENTATION_H

#include "settings.h"
#include <array>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <iostream>
#include <vector>

/// Convert rotation in the angle-axis representation to quaternion.
#define AA2Q(A, X, Y, Z) Eigen::Quaternionf(Eigen::AngleAxisf(A, Eigen::Vector3f(X, Y, Z) / sqrt(X * X + Y * Y + Z * Z))).normalized()

/// @name Constants
/// @{
const float MISORIENTATION_MAX = 63.8f; ///< Maximum misorientation angle. 63.8 degrees for cubic symmetry.
/// @}

/// @name Symmetry and Orientation Relationship Operators
/// @{
const int N_SYM = 24; ///< Number of symmetry operators for the cubic system

/// Cubic symmetry operators (0-23)
const Eigen::Quaternionf SYM[N_SYM] = {
    AA2Q(0.0f, 1.0f, 0.0f, 0.0f),
    AA2Q(M_PI, 1.0f, 0.0f, 0.0f),
    AA2Q(M_PI, 0.0f, 1.0f, 0.0f),
    AA2Q(M_PI, 0.0f, 0.0f, 1.0f),
    AA2Q(M_PI / 2.0f, 1.0f, 0.0f, 0.0f),
    AA2Q(M_PI / 2.0f, 0.0f, 1.0f, 0.0f),
    AA2Q(M_PI / 2.0f, 0.0f, 0.0f, 1.0f),
    AA2Q(M_PI / 2.0f, -1.0f, 0.0f, 0.0f),
    AA2Q(M_PI / 2.0f, 0.0f, -1.0f, 0.0f),
    AA2Q(M_PI / 2.0f, 0.0f, 0.0f, -1.0f),
    AA2Q(M_PI, 1.0f, 1.0f, 0.0f),
    AA2Q(M_PI, 1.0f, 0.0f, 1.0f),
    AA2Q(M_PI, 0.0f, 1.0f, 1.0f),
    AA2Q(M_PI, 1.0f, -1.0f, 0.0f),
    AA2Q(M_PI, -1.0f, 0.0f, 1.0f),
    AA2Q(M_PI, 0.0f, 1.0f, -1.0f),
    AA2Q(M_PI * 2.0f / 3.0f, 1.0f, 1.0f, 1.0f),
    AA2Q(M_PI * 2.0f / 3.0f, 1.0f, -1.0f, 1.0f),
    AA2Q(M_PI * 2.0f / 3.0f, -1.0f, 1.0f, 1.0f),
    AA2Q(M_PI * 2.0f / 3.0f, -1.0f, -1.0f, 1.0f),
    AA2Q(M_PI * 2.0f / 3.0f, 1.0f, 1.0f, -1.0f),
    AA2Q(M_PI * 2.0f / 3.0f, 1.0f, -1.0f, -1.0f),
    AA2Q(M_PI * 2.0f / 3.0f, -1.0f, 1.0f, -1.0f),
    AA2Q(M_PI * 2.0f / 3.0f, -1.0f, -1.0f, -1.0f)};

const int N_CSL19a = N_SYM; ///< Number of CSL19a equivalent representations (i.e., 27 degrees around <110>) orientation relationships

/** Symmetrically equivalent operators to apply the CSL19a (Ibe and Lucke) */
const Eigen::Quaternionf CSL19a[N_CSL19a] = {
    AA2Q(0.463f, -0.7071f, 0.0f, 0.7071f),
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
    CSL19a[0] * SYM[23]};

/**
 * @class Orientations
 * @brief Class to handle orientations and symmetry operations.
 *
 * This class manages the orientations and crystal symmetries of materials
 * in the OMicroN simulation. It provides functionalities to add and retrieve
 * orientations, convert between Euler angles and quaternions, and calculate
 * misorientations and other orientation-related metrics.
 *
 * Orientations are stored (provided/exported) as Euler angles, but all operations
 * are performed using quaternions.
 */
class Orientations
{
private:
    /// @name Geometry
    /// @{
    const float *mvpLimits; ///< Limits of Euler space
    /// @}

    /// @name Grain Boundaries
    /// @{
    float mvLAGB; ///< Low angle grain boundary minimum misorientation
    float mvHAGB; ///< High angle grain boundary minimum misorientation
    /// @}

    /// @name Simulation Data
    /// @{
    std::vector<Eigen::Vector3f> mvOrientation; ///< Orientations (Euler angles) of each phase
    std::vector<int> mvOrientationId;           ///< IDs of orientations (Euler angles) of each phase
    /// @}

public:
    /**
     * @brief Default constructor for the Orientations class.
     *
     * Initializes an instance of the Orientations class. This constructor does not
     * perform any specific initialization tasks beyond setting up the class instance.
     */
    Orientations(void);

    /** @brief Destructor (no need). */
    //~Orientations();

    /**
     * @brief Sets the orientation parameters from user settings.
     *
     * @param us UserSettings object containing the parameters.
     */
    void SetOrientationParameters(const UserSettings &us)
    {
        mvLAGB = us.mvLowerMisorientationCutOff;
        mvHAGB = us.mvHAGB;
    }

    /**
     * @brief Retrieves the Euler angles for a given orientation ID.
     *
     * @param id Orientation ID.
     * @return Euler angles as an Eigen::Vector3f.
     */
    Eigen::Vector3f GetEulerAngles(const int id) const
    {
        if (id < 0 || id >= mvOrientation.size())
        {
            std::cerr << "Error: Orientation ID out of bounds: " << id << std::endl;
            return Eigen::Vector3f::Zero(); // Return a default value or handle the error appropriately
        }
        return mvOrientation[id];
    }

    /**
     * @brief Prints all orientations for debugging purposes.
     */
    void DebugPrintOrientations() const
    {
        for (size_t i = 0; i < mvOrientation.size(); ++i)
        {
            std::cout << "Orientation " << i << ": " << mvOrientation[i].transpose() << std::endl;
        }
    }

    /// @name Euler Grid and Problem Size
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
     * @return ID of the added orientation.
     */
    int addOriAndReturnId(const Eigen::Vector3f &ea);

    /**
     * @brief Calculates the misorientation between two quaternions.
     *
     * @param q0 First quaternion.
     * @param q1 Second quaternion.
     * @return Misorientation angle.
     */
    float misorientationQ(const Eigen::Quaternionf q0, const Eigen::Quaternionf q1) const;

    /**
     * @brief Calculates the distance of the boundary from CSL19a.
     *
     * @param ea1 Euler angles of the first phase.
     * @param ea2 Euler angles of the second phase.
     * @return Misorientation angle.
     */
    float DistanceOfBoundaryFromCSL19a(const Eigen::Vector3f &ea1, const Eigen::Vector3f &ea2) const;

    /**
     * @brief Calculates the misorientation between two orientations given their IDs.
     *
     * @param Id1 ID of the first orientation.
     * @param Id2 ID of the second orientation.
     * @return Misorientation angle.
     */
    float CalculateMisorientationBetweenTwoOriIds(const int Id1, const int Id2) const;

    /**
     * @brief Calculates the misorientation between two orientations given their Euler angles.
     *
     * @param ea1 Euler angles of the first orientation.
     * @param ea2 Euler angles of the second orientation.
     * @return Misorientation angle.
     */
    float CalcMisorientationFromEulerAngles(const Eigen::Vector3f &ea1, const Eigen::Vector3f &ea2) const;
};

#endif // ORIENTATION_H
