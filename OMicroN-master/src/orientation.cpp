/* orientation.cpp
   OMicroN (optimising microstructures numerically) simulation program
   cpp-file containing orientation class implementation
*/

#include "orientation.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <eigen3/Eigen/Geometry>

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
 * @brief Generates a random crystal orientation
 * @return orientation in quaternion.
 */
static Eigen::Quaternionf GenerateRandomOrientation();


/**
 * @brief Converts a quaternion to Euler angles.
 *
 * @param q Quaternion to be converted.
 * @return Euler angles (ZYX convention, in degrees).
 */
static Eigen::Vector3f q2ea(const Eigen::Quaternionf &q);

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
 * @brief Adds a random  orientation and returns its ID.
 * @return ID of the added orientation.
 */
int Orientations::addRandomOriAndReturnId()
{
     Eigen::Quaternionf qRandom = GenerateRandomOrientation();
     Eigen::Vector3f ea = q2ea(qRandom);
    mvOrientation.push_back(ea);
    mvOrientationId.push_back(mvOrientation.size() - 1);
    // Uncomment for debugging purposes
    std::cout << "Added random orientation ID: " << mvOrientation.size() - 1 << " with Euler angles: " << ea.transpose() << std::endl;
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
 * @brief Converts a quaternion to Euler angles.
 *
 * @param q Quaternion to be converted.
 * @return Euler angles (ZYX convention, in degrees).
 */
static Eigen::Vector3f q2ea(const Eigen::Quaternionf &q) {
    // Convert quaternion to rotation matrix
    Eigen::Matrix3f R = q.toRotationMatrix();

    Eigen::Vector3f ea; // Euler angles: [Z, X, Z]

    // Calculate Euler angles based on the Z-X-Z convention
    if (std::abs(R(2, 2)) < 1.0f - std::numeric_limits<float>::epsilon()) {
        // General case
        ea[1] = std::atan2(-R(2, 0), std::sqrt(R(0, 0) * R(0, 0) + R(1, 0) * R(1, 0))); // X rotation
        ea[0] = std::atan2(R(1, 0), R(0, 0));                                           // Z1 rotation
        ea[2] = std::atan2(R(2, 1), R(2, 2));                                           // Z2 rotation
    } else {
        // Handle gimbal lock (R(2,2) close to ±1)
        ea[1] = (R(2, 2) > 0) ? 0.0f : static_cast<float>(M_PI); // X rotation
        ea[0] = std::atan2(R(0, 1), R(1, 1));                   // Z1 rotation
        ea[2] = 0.0f;                                           // Z2 rotation
    }

    // Convert angles from radians to degrees
    ea = ea * (180.0f / static_cast<float>(M_PI));

    return ea;
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
 * @brief Generates a random crystal orientation
 * @return orientation in quaternion.
 */

static Eigen::Quaternionf GenerateRandomOrientation() {
    static std::random_device rd; // Seed for randomness
    static std::mt19937 gen(rd()); // Mersenne Twister generator
    static std::uniform_real_distribution<float> dist(0.0f, 1.0f); // Uniform distribution [0, 1]

    Eigen::Quaternionf q;
    float r;

    // Generate 3 random numbers for quaternion generation
    r = dist(gen);
    q.w() = q.x() = std::sqrt(1.0f - r);
    q.y() = q.z() = std::sqrt(r);
    r = 2.0f * static_cast<float>(M_PI) * dist(gen);
    q.w() *= std::sin(r);
    q.x() *= std::cos(r);
    r = 2.0f * static_cast<float>(M_PI) * dist(gen);
    q.y() *= std::sin(r);
    q.z() *= std::cos(r);

    return q.normalized();
}


// Function to apply a specified symmetry operation to a crystal direction
std::vector<float> Orientations::ApplyCubicSymmetry(const std::vector<float>& crystalDirection, int symmetryIndex)
{
    // Convert the input direction to an Eigen vector
    Eigen::Vector3f direction(crystalDirection[0], crystalDirection[1], crystalDirection[2]);

    // Apply the specified symmetry operation using the SYM array
    Eigen::Vector3f transformedDirection = SYM[symmetryIndex] * direction;

    // Return the transformed direction as a std::vector<float>
    return {transformedDirection[0], transformedDirection[1], transformedDirection[2]};
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



// Helper function to convert Euler angles (in radians) to a 3x3 rotation matrix
std::vector<std::vector<float>> Orientations::EulerToRotationMatrix(float phi1, float PHI, float phi2)
{
   
    // Precompute sines and cosines of Euler angles
    float c1 = std::cos(d2r(phi1)), s1 = std::sin(d2r(phi1));
    float c2 = std::cos(d2r(phi2)), s2 = std::sin(d2r(phi2));
    float c = std::cos(d2r(PHI)), s = std::sin(d2r(PHI));

    // Define the rotation matrix
    std::vector<std::vector<float>> rotationMatrix(3, std::vector<float>(3, 0.0f));
    rotationMatrix[0][0] = c1 * c2 - s1 * s2 * c;
    rotationMatrix[0][1] = s1 * c2 + c1 * s2 * c;
    rotationMatrix[0][2] = s2 * s;

    rotationMatrix[1][0] = -c1 * s2 - s1 * c2 * c;
    rotationMatrix[1][1] = -s1 * s2 + c1 * c2 * c;
    rotationMatrix[1][2] = c2 * s;

    rotationMatrix[2][0] = s1 * s;
    rotationMatrix[2][1] = -c1 * s;
    rotationMatrix[2][2] = c;

    return rotationMatrix;
}



// Function to convert Euler angles (phi1, Phi, phi2) to a rotation matrix
Eigen::Matrix3d Orientations::GetOriInRotationMatrix(int id) {
    // Extract the Euler angles (in radians)
    Eigen::Vector3f ea1 = mvOrientation[id];
    double phi1 = ea1[0];
        double Phi = ea1[1];
    double phi2 = ea1[2];

    // Compute trigonometric values
    double c1 = cos(phi1), s1 = sin(phi1);
    double c2 = cos(phi2), s2 = sin(phi2);
    double cPhi = cos(Phi), sPhi = sin(Phi);

    // Construct the rotation matrix using the Bunge convention
    Eigen::Matrix3d R;
    R(0, 0) = c1 * c2 - s1 * s2 * cPhi;
    R(0, 1) = -c1 * s2 - s1 * c2 * cPhi;
    R(0, 2) = s1 * sPhi;

    R(1, 0) = s1 * c2 + c1 * s2 * cPhi;
    R(1, 1) = -s1 * s2 + c1 * c2 * cPhi;
    R(1, 2) = -c1 * sPhi;

    R(2, 0) = s2 * sPhi;
    R(2, 1) = c2 * sPhi;
    R(2, 2) = cPhi;

    return R;
}