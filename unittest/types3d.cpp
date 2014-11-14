#include "ei/3dfunctions.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace ei;
using namespace std;

bool test_3dtypes()
{
    bool result = true;

    // ********************************************************************* //
    // Test volume() and surface() function
    {
        Sphere sph( Vec3(1.0f, 2.0f, 3.14159f), 0.75f );
        Box box( Vec3(1.0f, 1.0f, 1.0f), Vec3(2.0f, 2.5f, 3.0f) );
        Triangle tri( Vec3(0.0f), Vec3(1.0f, 1.0f, 0.0f), Vec3(0.0f, 0.0f, 1.0f) );
        Plane pla( Vec3(1.0f, 0.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f) );
        TEST( volume(sph) == 1.767145868f, "Volume of a sphere wrong!" );
        TEST( volume(box) == 3.0f, "Volume of a box wrong!" );
        TEST( volume(tri) == 0.0f, "Volume of a triangle wrong!" );
        TEST( volume(pla) == 0.0f, "Volume of a plane wrong!" );

        TEST( surface(sph) == 7.068583471f, "Surface of a sphere wrong!" );
        TEST( surface(box) == 13.0f, "Surface of a box wrong!" );
        TEST( surface(tri) == 0.707106781f, "Surface of a triangle wrong!" );
        TEST( surface(pla) == std::numeric_limits<float>::infinity(), "Surface of a plane wrong!" );
    }

    // Test plane construction
    {
        Plane pla0( Vec3(1.0f, 0.0f, 0.0f), 0.0f );
        Plane pla1( Vec3(1.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f) );
        Plane pla2( Vec3(0.0f, -1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, 2.0f, 4.0f) );
        TEST( all(pla0.n == pla1.n) && pla0.d == pla1.d, "The constructed planes (0,1) should be equal!" );
        TEST( all(pla0.n == pla2.n) && pla0.d == pla2.d, "The constructed planes (0,2) should be equal!" );
    }

    return result;
}