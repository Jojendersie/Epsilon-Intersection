#include "ei/3dfunctions.hpp"
#include "ei/3dintersection.hpp"
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
        Disc dis( Vec3(1.0f), normalize(Vec3(-1.0f)), 1.0f );
        Plane pla( Vec3(1.0f, 0.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f) );
        DOP dop( Vec3(1.0f, 0.0f, 0.0f), -0.5f, 1.5f );
        Ellipsoid ell( Vec3(-1.0f, -0.5f, -0.5f), Vec3(1.5f, 0.75f, 0.75f) );
        TEST( volume(sph) == 1.767145868f, "Volume of a sphere wrong!" );
        TEST( volume(box) == 3.0f, "Volume of a box wrong!" );
        TEST( volume(tri) == 0.0f, "Volume of a triangle wrong!" );
        TEST( volume(dis) == 0.0f, "Volume of a disc wrong!" );
        TEST( volume(pla) == 0.0f, "Volume of a plane wrong!" );
        TEST( volume(dop) == 0.0f, "Volume of a DOP wrong!" );
        TEST( volume(ell) == 3.53429174f, "Volume of an ellipsoid wrong!" );

        TEST( surface(sph) == 7.068583471f, "Surface of a sphere wrong!" );
        TEST( surface(box) == 13.0f, "Surface of a box wrong!" );
        TEST( surface(tri) == 0.707106781f, "Surface of a triangle wrong!" );
        TEST( surface(dis) == PI, "Surface of a disc wrong!" );
        TEST( surface(pla) == std::numeric_limits<float>::infinity(), "Surface of a plane wrong!" );
        TEST( surface(dop) == std::numeric_limits<float>::infinity(), "Surface of a DOP wrong!" );
        TEST( abs(surface(ell) / 12.0816f - 1.0f) < 0.012f, "Surface approximation of an ellipsoid too far away!" );
    }

    // Test plane construction
    {
        Plane pla0( Vec3(1.0f, 0.0f, 0.0f), 0.0f );
        Plane pla1( Vec3(1.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f) );
        Plane pla2( Vec3(0.0f, -1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, 2.0f, 4.0f) );
        TEST( all(pla0.n == pla1.n) && pla0.d == pla1.d, "The constructed planes (0,1) should be equal!" );
        TEST( all(pla0.n == pla2.n) && pla0.d == pla2.d, "The constructed planes (0,2) should be equal!" );
    }

    // Test ellipsoid construction and intersection with a point
    {
        Box box( Vec3(-2.0f, -1.0f, -1.0f), Vec3(0.0f, 0.0f, 0.0f) );
        Ellipsoid ell0( Vec3(-1.0f, -0.5f, -0.5f), Vec3(1.732050808f, 0.866025404f, 0.866025404f) );
        Ellipsoid ell1( box );
        TEST( all(ell0.center == ell1.center), "Center of ellipsoid wrong!" );
        TEST( all(ell0.radii == ell1.radii), "Radii of ellipsoid wrong!" );
        TEST( intersects( box.min, ell0 ), "The bounding box min must be contained in the ellipsoid!" );
        TEST( intersects( box.max, ell0 ), "The bounding box max must be contained in the ellipsoid!" );
        TEST( !intersects( box.max + 1e-6f, ell0 ), "The point must be outside the ellipsoid!" );
    }

    // Test sphere construction
    {
        Box box( Vec3(-2.0f, -1.0f, -1.0f), Vec3(0.0f, 0.0f, 0.0f) );
        Sphere sph0( Vec3(-1.0f, -0.5f, -0.5f), 1.224744871f );
        Sphere sph1( box );
        TEST( all(sph0.center == sph1.center), "Center of sphere wrong!" );
        TEST( sph0.radius == sph1.radius, "Radius of sphere wrong!" );
    }

    // Test box construction
    {
        Box box0( Vec3(-2.0f, -1.0f, -1.0f), Vec3(0.0f, 0.0f, 0.0f) );
        Box box1( Vec3(3.0f), Vec3(7.0f) );
        Ellipsoid ell( Vec3(-1.0f, -0.5f, -0.5f), Vec3(1.0f, 0.5f, 0.5f) );
        Triangle tri( Vec3(-2.0f, -0.5f, -0.5f), Vec3(0.0, -1.0f, 0.0f), Vec3(0.0f, 0.0f, -1.0f) );
        Sphere sph( Vec3(5.0f), 2.0f );
        Box box2( ell );
        Box box3( tri );
        Box box4( sph );
        Box box5( box0, box1 );
        TEST( all(box0.min == box2.min), "Ellipsoid bounding min is wrong!" );
        TEST( all(box0.max == box2.max), "Ellipsoid bounding min is wrong!" );
        TEST( all(box0.min == box3.min), "Triangle bounding min is wrong!" );
        TEST( all(box0.max == box3.max), "Triangle bounding max is wrong!" );
        TEST( all(box1.min == box4.min), "Sphere bounding min is wrong!" );
        TEST( all(box1.max == box4.max), "Sphere bounding max is wrong!" );
        TEST( all(box5.min == Vec3(-2.0f, -1.0f, -1.0f)), "2Box bounding min is wrong!" );
        TEST( all(box5.max == Vec3(7.0f)), "2Box bounding max is wrong!" );
    }

    return result;
}