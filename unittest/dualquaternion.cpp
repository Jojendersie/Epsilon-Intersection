#include "ei/dualquaternion.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace ei;

bool test_dualquaternion()
{
    bool result = true;

    { // Test conversion Matrix <-> Quaternion
        const DualQuaternion qq0 { 0.9740f, -0.0542f, 0.1446f, 0.1660f, -0.1924f, 1.6917f, 1.7685f, 0.1407f };
        const Mat3x4 m0 { qq0 };
        const Mat3x4 m1 { 0.9523f, -0.1536f,  0.2636f, -1.0189f,
                         -0.0576f, -0.9390f, -0.3391f, -2.9236f,
                          0.2996f,  0.3077f, -0.9031f,  3.8210f};
        const Quaternion q0 { Mat3x3{m1} };
        const DualQuaternion qq1 { q0, Vec3{-1.0189f, -2.9236f, 3.8210f} };
        TEST( approx(m0, m1, 2e-4f), "Conversion DualQuaternion -> Matrix wrong!" );
        TEST( approx(qq0, qq1, 2e-4f), "Conversion Matrix -> DualQuaternion wrong!" );
        const DualQuaternion qq2 { m1 };
        TEST( qq2 == qq1, "Conversion Matrix -> DualQuaternion (direct) wrong!" );
    }

    { // Test basic math opertators
        const DualQuaternion qq0{ Quaternion{0.1f, 0.3f, 0.4f}, Vec3{0.1f, 0.2f, -0.3f} };
        const DualQuaternion qq1{ Quaternion{-0.25f, 0.5f, -0.7f}, Vec3{-0.1f, 0.0f, 0.1f} };
        const DualQuaternion qqu{ 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f };
        const Mat3x4 m0 { qq0 };
        const Mat3x4 m1 { qq1 };
        const Mat3x4 m2 = transform( m0, m1 );
        DualQuaternion res = qq0 * qq1;
        TEST( approx( Mat3x4{res}, m2 ), "Dual quaternion multiplication wrong." );
        TEST( approx( len(res), 1.0f), "Dual quaternion multiplication does not preserve unity." );
        res = (qq0 * qq1) / qq1;
        TEST( approx(qq0, res), "Dual Quaternion: multiplication and division do not cancel each other.");
        res = qq0 / qq0;
        TEST( approx(qqu, res), "qq0/qq0 is not the union quaternion" );
        res = qq1 / qq1;
        TEST( approx(qqu, res), "qq1/qq1 is not the union quaternion" );
        res = (qq0 / qq1) * qq1;
        TEST( approx(qq0, res), "Dual Quaternion: division and multiplication do not cancel each other.");
        res = qq0 * qq1;
        TEST( approx(len(res), 1.0f), "Dual Quaternion: multiplication of unit quaternions should keep the length.");
    }

    return result;
}