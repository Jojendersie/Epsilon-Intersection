#include "ei/3dfunctions.hpp"
#include "ei/3dintersection.hpp"
#include "unittest.hpp"
#include "performance3d.hpp"

#include <iostream>

using namespace ei;

bool test_3dtypes()
{
    bool result = true;

    // ********************************************************************* //
    // Test volume() and surface() function
    {
        Sphere sph( Vec3(1.0f, 2.0f, 3.14159f), 0.75f );
        Box box( Vec3(1.0f, 1.0f, 1.0f), Vec3(2.0f, 2.5f, 3.0f) );
        Thetrahedron the( Vec3(1.0f, 0.0f, -1.0f/PHYTAGORAS), Vec3(-1.0f, 0.0f, -1.0f/PHYTAGORAS), Vec3(0.0f, 1.0f, 1.0f/PHYTAGORAS), Vec3(0.0f, -1.0f, 1.0f/PHYTAGORAS));
        Triangle tri( Vec3(0.0f), Vec3(1.0f, 1.0f, 0.0f), Vec3(0.0f, 0.0f, 1.0f) );
        Disc dis( Vec3(1.0f), normalize(Vec3(-1.0f)), 1.0f );
        Plane pla( Vec3(1.0f, 0.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f) );
        DOP dop( Vec3(1.0f, 0.0f, 0.0f), -0.5f, 1.5f );
        Ellipsoid ell( Vec3(-1.0f, -0.5f, -0.5f), Vec3(1.5f, 0.75f, 0.75f) );
        Ray ray( Vec3(0.0f), Vec3(1.0f, 0.0f, 0.0f) );
        Segment seg( Vec3(0.0f), Vec3(2.0f, 0.0f, 0.0f) );
        Capsule cap( Vec3(0.0f), Vec3(0.0f, 1.0f, 0.0f), 0.5f);
        Frustum fru( Vec3(0.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 1.0f, 0.0f), 1.0f, 2.0f, 0.0f, 1.0f, 0.0f, 1.0f);
        TEST( volume(sph) == 1.767145868f, "Volume of a sphere wrong!" );
        TEST( volume(box) == 3.0f, "Volume of a box wrong!" );
        TEST( volume(the) == 0.942809042f, "Volume of a tetrahedron wrong!" );
        TEST( volume(tri) == 0.0f, "Volume of a triangle wrong!" );
        TEST( volume(dis) == 0.0f, "Volume of a disc wrong!" );
        TEST( volume(pla) == 0.0f, "Volume of a plane wrong!" );
        TEST( volume(dop) == 0.0f, "Volume of a DOP wrong!" );
        TEST( volume(ell) == 3.53429174f, "Volume of an ellipsoid wrong!" );
        TEST( volume(ray) == 0.0f, "Volume of a ray wrong!" );
        TEST( volume(seg) == 0.0f, "Volume of a segment wrong!" );
        TEST( volume(cap) == 1.30899704f, "Volume of a capsule wrong!" );
        TEST( volume(fru) == 0.33333333f, "Volume of a frustum wrong!" );

        TEST( surface(sph) == 7.068583471f, "Surface of a sphere wrong!" );
        TEST( surface(box) == 13.0f, "Surface of a box wrong!" );
        TEST( surface(the) == 6.92820323f, "Surface of a tetrahedron wrong!" );
        TEST( surface(tri) == 0.707106781f, "Surface of a triangle wrong!" );
        TEST( surface(dis) == PI, "Surface of a disc wrong!" );
        TEST( surface(pla) == std::numeric_limits<float>::infinity(), "Surface of a plane wrong!" );
        TEST( surface(dop) == std::numeric_limits<float>::infinity(), "Surface of a DOP wrong!" );
        TEST( abs(surface(ell) / 12.0816f - 1.0f) < 0.012f, "Surface approximation of an ellipsoid too far away!" );
        TEST( surface(ray) == 0.0f, "Surface of a ray wrong!" );
        TEST( surface(seg) == 0.0f, "Surface of a segment wrong!" );
        TEST( surface(cap) == 6.283185307f, "Surface of a capsule wrong!" );
        Vec3 a(1.0f, 0.0f, 1.0f), b(2.0f, 0.0f, 1.0f), c(2.0f, 1.0f, 1.0f), d(1.0f, 1.0f, 1.0f);
        float refA3_A5 = 0.5f * (len(cross(a, b)) + len(cross(c, d)));
        float refA4_A6 = 0.5f * (len(cross(b, c)) + len(cross(d, a)));
        TEST( surface(fru) == refA3_A5 + refA4_A6 + 1.0f, "Surface of a frustum wrong!" );
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

    // ********************************************************************* //
    // Test distance()
    {
        Vec3 poi0(0.0f, 1.0f, 1.0f);
        Vec3 poi1(1.0f, 1.0f, 1.0f);
        Vec3 poi2(1.0f, 3.0f, 3.0f);
        Vec3 poi3(0.25f, 0.5f, 0.25f);
        Vec3 poi4(0.0f, 0.0f, -0.1f);
        Segment seg0(Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 2.0f, 2.0f));
        Segment seg1(poi1, poi2);
        Segment seg2(Vec3(-1.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f));
        Segment seg3(Vec3(-1.0f, 1.0f, 1.0f), Vec3(1.0f, 1.0f, 1.0f));
        Capsule cap0(Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 1.0f), 0.5f);
        Capsule cap1(Vec3(0.0f, 1.0f, 0.5f), Vec3(0.0f, 1.0f, 2.0f), 0.25f);
        Triangle tri0(Vec3(0.0f), Vec3(1.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 1.0f));
        Sphere sph0(poi0, 1.0f);
        Sphere sph1(poi2, 0.5f);
        TEST( distance(poi0, poi0) == 0.0f, "Distance between two equal points is 0!" );
        TEST( distance(poi0, poi1) == 1.0f, "Distance between poi0 and poi1 is 1!" );
        TEST( distance(poi0, seg0) == 0.0f, "Distance between poi0 and seg0 is 0!" );
        TEST( distance(poi1, seg0) == 1.0f, "Distance between poi1 and seg0 is 1!" );
        TEST( distance(poi2, seg0) == sqrt(3.0f), "Distance between poi2 and seg0 is sqrt(3)!" );
        TEST( distance(seg0, seg1) == 1.0f, "Distance between seg0 and seg1 is 1!" );
        TEST( distance(seg0, seg2) == 0.707106769f, "Distance between seg0 and seg2 is 1/sqrt(2)!" );
        TEST( distance(seg0, seg3) == 0.0f, "Distance between seg0 and seg3 is 0!" );
        TEST( distance(cap0, cap1) == 0.25f, "Distance between cap0 and cap1 is 0.25!" );
        TEST( distance(poi0, tri0) == 1.0f, "Distance between poi0 and tri0 is 1.0!" );
        TEST( distance(poi1, tri0) == sqrt(1.5f), "Distance between poi1 and tri0 is sqrt(1.5)!" );
        TEST( distance(poi3, tri0) == 0.5f, "Distance between poi3 and tri0 is 0.5!" );
        TEST( distance(poi0, sph0) == -1.0f, "Distance between poi0 and sph0 is -1.0!");
        TEST( distance(poi1, sph0) == 0.0f, "Distance between poi1 and sph0 is 0.0!");
        TEST( distance(poi2, sph0) == 2.0f, "Distance between poi2 and sph0 is 2.0!");
        TEST( distance(poi0, cap0) == 0.5f, "Distance between poi0 and cap0 is 0.5!");
        TEST( distance(poi4, cap0) == -0.4f, "Distance between poi4 and cap0 is -0.4!");
        TEST( distance(poi0, cap1) == -0.25f, "Distance between poi0 and cap0 is -0.25!");
        TEST( distance(sph0, cap0) == 0.0f, "Distance between sph0 and cap0 is 0!");
        TEST( distance(sph0, cap1) == 0.0f, "Distance between sph0 and cap1 is 0!");
        TEST( distance(sph1, cap0) == 2.741657387f, "Distance between sph1 and cap0 is 0!");

        performance<Vec3,Triangle,float>(distance, "distance");
    }

    return result;
}