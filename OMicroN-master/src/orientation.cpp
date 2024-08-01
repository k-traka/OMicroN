/* orientation.cpp
   OMicroN (optimising microstructures numerically) simulation program
   cpp-file containing orientation class implementation
*/

#include "orientation.h"
#include <cmath>
#include <cstdio>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265359
#endif

/**
 * @brief Converts degrees to radians.
 *
 * @param d Angle in degrees.
 * @return Angle in radians.
 */
static float d2r(const float d) { return d * M_PI / 180.; }

/**
 * @brief Converts radians to degrees.
 *
 * @param r Angle in radians.
 * @return Angle in degrees.
 */
static float r2d(const float r) { return r * 180. / M_PI; }

/**
 * @brief Default constructor for the Orientations class.
 *
 * This constructor initializes an instance of the Orientations class. As it is
 * a default constructor, it does not perform any specific actions or initialize
 * any member variables.
 */
Orientations::Orientations(void)
{
}

/**
 * @brief Adds an orientation and returns its ID.
 *  @param ea Euler angles.
 * @return ID of the added orientation.
 */
int Orientations::addOriAndReturnId(const Eigen::Vector3f &ea)
{
    mvOrientation.push_back(ea);
    mvOrientationId.push_back(mvOrientation.size() - 1);
    // Uncomment for debugging purposes
    // std::cout << "Added orientation ID: " << mvOrientation.size() - 1 << " with Euler angles: " << ea.transpose() << std::endl;
    return mvOrientation.size() - 1;
}

/**
 * @brief Converts Euler angles to a quaternion.
 *
 * @param ea Euler angles.
 * @return Quaternion corresponding to the Euler angles.
 */
static Eigen::Quaternionf ea2q(const Eigen::Vector3f &ea)
{
    Eigen::Quaternionf q;
    q = Eigen::Quaternionf(Eigen::AngleAxisf(d2r(ea[0]), Eigen::Vector3f::UnitZ()) * Eigen::AngleAxisf(d2r(ea[1]), Eigen::Vector3f::UnitX()) * Eigen::AngleAxisf(d2r(ea[2]), Eigen::Vector3f::UnitZ()));
    return q.normalized();
}

/**
 * @brief Calculates the misorientation angle between two orientations given their IDs.
 *
 * @param Id1 ID of the first orientation.
 * @param Id2 ID of the second orientation.
 * @return Misorientation angle in degrees.
 */
float Orientations::CalculateMisorientationBetweenTwoOriIds(const int Id1, const int Id2) const
{
    float m, w;
    Eigen::Vector3f ea1 = mvOrientation[Id1];
    Eigen::Vector3f ea2 = mvOrientation[Id2];
    Eigen::Quaternionf q0 = ea2q(ea1);
    Eigen::Quaternionf q1 = ea2q(ea2);

    Eigen::Quaternionf qs;

    m = fabs(q0.w() * q1.w() + q0.x() * q1.x() + q0.y() * q1.y() + q0.z() * q1.z());
    for (int i = 1; i < N_SYM; i++)
    {
        qs = q1 * SYM[i];
        w = fabs(q0.w() * qs.w() + q0.x() * qs.x() + q0.y() * qs.y() + q0.z() * qs.z());
        if (w > m)
            m = w;
    }
    float mis = 180.0 * 2.0 * acos(m) / M_PI;
    if (mis > MISORIENTATION_MAX)
    {
        std::cout << " check ori Id1 " << Id1 << " and Id2 " << Id2 << " which give misorientation higher than the symmetry allows " << std::endl;
        exit(-1);
    }
    return mis;
}

/**
 * @brief Calculates the misorientation angle between two orientations given their Euler angles.
 *
 * @param ea1 Euler angles of the first orientation.
 * @param ea2 Euler angles of the second orientation.
 * @return Misorientation angle in degrees.
 */
float Orientations::CalcMisorientationFromEulerAngles(const Eigen::Vector3f &ea1, const Eigen::Vector3f &ea2) const
{
    float m, w;

    Eigen::Quaternionf q0 = ea2q(ea1);
    Eigen::Quaternionf q1 = ea2q(ea2);

    Eigen::Quaternionf qs;

    m = fabs(q0.w() * q1.w() + q0.x() * q1.x() + q0.y() * q1.y() + q0.z() * q1.z());
    for (int i = 1; i < N_SYM; i++)
    {
        qs = q1 * SYM[i];
        w = fabs(q0.w() * qs.w() + q0.x() * qs.x() + q0.y() * qs.y() + q0.z() * qs.z());
        if (w > m)
            m = w;
    }
    float mis = 180.0 * 2.0 * acos(m) / M_PI;
    if (mis > MISORIENTATION_MAX)
    {
        std::cout << " check ori phi1 " << ea1[0] << " phi " << ea1[1] << " phi2 " << ea1[2] << " with phi1 " << ea2[0] << " phi " << ea2[1] << " phi2 " << ea2[2] << " which give misorientation higher than the symmetry allows " << std::endl;
        exit(-1);
    }
    return mis;
}

/**
 * @brief Calculates the misorientation angle between two quaternions.
 *
 * @param q0 First quaternion.
 * @param q1 Second quaternion.
 * @return Misorientation angle in degrees.
 */
float Orientations::misorientationQ(const Eigen::Quaternionf q0, const Eigen::Quaternionf q1) const
{
    float m, w;
    Eigen::Quaternionf qs = q0 * q1.conjugate();

    m = fabs(qs.w());
    for (int i = 1; i < N_SYM; i++)
    {
        w = fabs(SYM[i].w() * qs.w() - SYM[i].x() * qs.x() - SYM[i].y() * qs.y() - SYM[i].z() * qs.z());
        if (w > m)
            m = w;
    }
    return 180.0 * 2.0 * acos(m) / M_PI;
}

/**
 * @brief Calculates the distance of the boundary from the CSL19a for given Euler angles.
 *
 * @param ea1 Euler angles of the first orientation.
 * @param ea2 Euler angles of the second orientation.
 * @return Misorientation angle in degrees.
 */
float Orientations::DistanceOfBoundaryFromCSL19a(const Eigen::Vector3f &ea1, const Eigen::Vector3f &ea2) const
{
    float m, w;
    for (int ii = 0; ii < 2; ii++)
    { // switch symmetry must also be applied
        Eigen::Quaternionf q0 = (ii < 1) ? ea2q(ea1) : ea2q(ea2);
        Eigen::Quaternionf q1 = (ii < 1) ? ea2q(ea2) : ea2q(ea1);
        Eigen::Quaternionf RotatedQ0 = q0 * CSL19a[0].conjugate();
        Eigen::Quaternionf RotatedQ0ToQ1 = RotatedQ0 * q1.conjugate();
        if (ii < 1)
            m = fabs(RotatedQ0ToQ1.w());
        else
        {
            if (m < fabs(RotatedQ0ToQ1.w()))
                m = fabs(RotatedQ0ToQ1.w());
        }
        for (int i = 1; i < N_SYM; i++)
        {
            Eigen::Quaternionf RotatedQ0 = q0 * CSL19a[i].conjugate();
            Eigen::Quaternionf RotatedQ0ToQ1 = RotatedQ0 * q1.conjugate();
            for (int k = 1; i < N_SYM; i++)
            {
                w = fabs(SYM[k].w() * RotatedQ0ToQ1.w() - SYM[k].x() * RotatedQ0ToQ1.x() - SYM[k].y() * RotatedQ0ToQ1.y() - SYM[k].z() * RotatedQ0ToQ1.z());
                if (w > m) // larger cosine => lower misorientation angle
                    m = w;
            }
        }
    }
    float mis = 180.0 * 2.0 * acos(m) / M_PI;
    return mis;
}
