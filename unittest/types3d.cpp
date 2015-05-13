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
        Capsule cap( Vec3(0.0f), Vec3(0.0f, 1.0f, 0.0f), 0.5f );
        Frustum fru( Vec3(0.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 1.0f, 0.0f), 1.0f, 2.0f, 0.0f, 1.0f, 0.0f, 1.0f );
        Frustum fru2( Vec3(0.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 1.0f, 0.0f), 1.0f, 2.0f, 0.0f, 1.0f, 0.5f, 1.0f );
        Frustum fru3( Vec3(0.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 1.0f, 0.0f), 0.5f, 1.0f, 0.0f, 0.5f, 0.0f, 0.5f );
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
        float vol0 = volume(fru2), vol1 = volume(fru3);
        TEST( vol0+vol1 == 0.33333333f, "Volume of frusta wrong!" );

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

        TEST( all(center(sph) == Vec3(1.0f, 2.0f, 3.14159f)), "Center of sphere wrong!" );
        TEST( all(center(box) == Vec3(1.5f, 1.75f, 2.0f)), "Center of box wrong!" );
        TEST( all(center(the) == Vec3(0.0f)), "Center of thetrahedron wrong!" );
        TEST( all(center(tri) == Vec3(1.0f/3.0f)), "Center of trinagle wrong!" );
        TEST( all(center(dis) == Vec3(1.0f)), "Center of disc wrong!" );
        TEST( all(center(ell) == Vec3(-1.0f, -0.5f, -0.5f)), "Center of ellipsoid wrong!" );
        TEST( all(center(seg) == Vec3(1.0f, 0.0f, 0.0f)), "Center of line segment wrong!" );
        TEST( all(center(cap) == Vec3(0.0f, 0.5f, 0.0f)), "Center of capsule wrong!" );
        TEST( all(center(fru) == Vec3(1.125f, 0.375f, 0.75f)), "Center of frustum wrong!" );
        Vec3 c0 = center(fru2), c1 = center(fru3);
        TEST( approx((c0*vol0 + c1*vol1) / (vol0 + vol1), Vec3(1.125f, 0.375f, 0.75f)), "Center of frustum wrong!" );
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
        OBox obo0( Vec3(0.0f), Vec3(1.0f, 2.0f, 3.0f), Quaternion(Vec3(0.0f, 0.0f, 1.0f), PI/2) );
        OBox obo1( Vec3(1.0f), Vec3(1.0f, 1.0f, 1.0f), Quaternion(Vec3(0.0f, 0.0f, 1.0f), PI/4) );
        Box box2( ell );
        Box box3( tri );
        Box box4( sph );
        Box box5( box0, box1 );
        Box box6( obo0 );
        Box box7( obo1 );
        // Box box8( OBox(box4) ); // <- This is a function? VC12 bug or standard?
        Box box8 = Box(OBox(box4));
        TEST( all(box0.min == box2.min), "Ellipsoid bounding min is wrong!" );
        TEST( all(box0.max == box2.max), "Ellipsoid bounding min is wrong!" );
        TEST( all(box0.min == box3.min), "Triangle bounding min is wrong!" );
        TEST( all(box0.max == box3.max), "Triangle bounding max is wrong!" );
        TEST( all(box1.min == box4.min), "Sphere bounding min is wrong!" );
        TEST( all(box1.max == box4.max), "Sphere bounding max is wrong!" );
        TEST( all(box5.min == Vec3(-2.0f, -1.0f, -1.0f)), "Box bounding min is wrong!" );
        TEST( all(box5.max == Vec3(7.0f)), "Box bounding max is wrong!" );
        TEST( all(box6.min == Vec3(-1.0f, -0.5f, -1.5f)), "OBox bounding min is wrong!" );
        TEST( all(box6.max == Vec3(1.0f, 0.5f, 1.5f)), "OBox bounding max is wrong!" );
        TEST( all(box7.min == Vec3(0.292893219f, 0.292893219f, 0.5f)), "OBox2 bounding min is wrong!" );
        TEST( all(box7.max == Vec3(1.707106781f, 1.707106781f, 1.5f)), "OBox2 bounding max is wrong!" );
        TEST( all(box8.min == box4.min), "Box(OBox(box4)) introduces a bias to min!" );
        TEST( all(box8.max == box4.max), "Box(OBox(box4)) introduces a bias to max!" );
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
        Sphere sph2(Vec3(-1.0f), 0.5f);
        Box box0(Vec3(0.0f), Vec3(2.0f));
        Plane pla0( Vec3(0.0f, 1.0f, 0.0f), -1.0f);
        DOP dop0( Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, 2.0f, 0.0f), Vec3(0.0f, 4.0f, 0.0f) );
        DOP dop1( Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, 4.0f, 0.0f), Vec3(0.0f, 2.0f, 0.0f) );
        DOP dop2( Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, -2.0f, 0.0f), Vec3(0.0f, 2.0f, 0.0f) );
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
        TEST( distance(sph1, cap0) == 2.741657387f, "Distance between sph1 and cap0 is 2.741657387!");
        TEST( distance(sph0, seg0) == 0.0f, "Distance between sph0 and seg0 is 0!");
        TEST( distance(sph0, cap0.seg) == 0.0f, "Distance between sph0 and cap0.seg is 0.5!");
        TEST( distance(sph1, seg2) == 3.10555124f, "Distance between sph1 and seg2 is sqrt(13)-0.5!");
        TEST( distance(poi0, box0) == 0.0f, "Distance between poi0 and box0 is 0!");
        TEST( distance(poi1, box0) == -1.0f, "Distance between poi1 and box0 is -1!");
        TEST( distance(poi2, box0) == sqrt(2.0f), "Distance between poi2 and box0 is sqrt(2)!");
        TEST( distance(sph0, box0) == 0.0f, "Distance between sph0 and box0 is 0!");
        TEST( distance(sph1, box0) == 0.914213538f, "Distance between sph1 and box0 is sqrt(2)-0.5!");
        TEST( distance(poi0, pla0) == 0.0f, "Distance between poi0 and pla0 is 0.0!");
        TEST( distance(poi2, pla0) == 2.0f, "Distance between poi2 and pla0 is 2.0!");
        TEST( distance(poi3, pla0) == -0.5f, "Distance between poi3 and pla0 is -0.5!");
        TEST( distance(sph0, pla0) == 0.0f, "Distance between sph0 and pla0 is 0.0!");
        TEST( distance(sph1, pla0) == 1.5f, "Distance between sph1 and pla0 is 1.5!");
        TEST( distance(sph2, pla0) == -1.5f, "Distance between sph2 and pla0 is -1.5!");
        TEST( distance(poi4, dop0) == -2.0f, "Distance between poi4 and dop0 is -2.0!");
        TEST( distance(poi4, dop1) == -2.0f, "Distance between poi4 and dop1 is -2.0!");
        TEST( distance(poi4, dop2) == 0.0f, "Distance between poi4 and dop2 is 0.0!");
        TEST( distance(poi0, dop2) == 0.0f, "Distance between poi0 and dop2 is -1.0!");
        TEST( distance(poi2, dop2) == 1.0f, "Distance between poi2 and dop2 is 1.0!");

        performance<Vec3,Triangle,float>(distance, "distance");
    }

    return result;
}