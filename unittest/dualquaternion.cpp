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
    }

    return result;
}