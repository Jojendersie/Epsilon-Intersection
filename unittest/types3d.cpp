#include "ei/3dintersection.hpp"
#include "ei/algorithm.hpp"
#include "unittest.hpp"
#include "performance3d.hpp"

#include <iostream>
#include <limits>

using namespace ei;

bool test_3dtypes()
{
    bool result = true;

    // ********************************************************************* //
    // Test volume() and surface() function
    {
        Sphere sph( Vec3(1.0f, 2.0f, 3.14159f), 0.75f );
        Box box( Vec3(1.0f, 1.0f, 1.0f), Vec3(2.0f, 2.5f, 3.0f) );
        OBox obo( Vec3(0.5f), Vec3(0.5f), Quaternion(PI, PI/2, PI/2) );
        Tetrahedron the( Vec3(1.0f, 0.0f, -1.0f/PHYTAGORAS), Vec3(-1.0f, 0.0f, -1.0f/PHYTAGORAS), Vec3(0.0f, 1.0f, 1.0f/PHYTAGORAS), Vec3(0.0f, -1.0f, 1.0f/PHYTAGORAS));
        Triangle tri( Vec3(0.0f), Vec3(1.0f, 1.0f, 0.0f), Vec3(0.0f, 0.0f, 1.0f) );
        Disc dis( Vec3(1.0f), normalize(Vec3(-1.0f)), 1.0f );
        Plane pla( Vec3(1.0f, 0.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f) );
        DOP dop( Vec3(1.0f, 0.0f, 0.0f), -0.5f, 1.5f );
        Ellipsoid ell( Vec3(-1.0f, -0.5f, -0.5f), Vec3(1.5f, 0.75f, 0.75f) );
        OEllipsoid oel( Vec3(-1.0f, -0.5f, -0.5f), Vec3(1.5f, 0.75f, 0.75f), Quaternion(0.235f, -2.352f, 1.43f));
        Ray ray( Vec3(0.0f), Vec3(1.0f, 0.0f, 0.0f) );
        Segment seg( Vec3(0.0f), Vec3(2.0f, 0.0f, 0.0f) );
        Cone con( Vec3(0.0f), Vec3(1.0f, 0.0f, 0.0f), 1.0f, 1.0f );
        Capsule cap( Vec3(0.0f), Vec3(0.0f, 1.0f, 0.0f), 0.5f );
        Frustum fru( Vec3(0.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 1.0f, 0.0f), 1.0f, 2.0f, 0.0f, 1.0f, 0.0f, 1.0f );
        Frustum fru2( Vec3(0.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 1.0f, 0.0f), 1.0f, 2.0f, 0.0f, 1.0f, 0.5f, 1.0f );
        Frustum fru3( Vec3(0.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 1.0f, 0.0f), 0.5f, 1.0f, 0.0f, 0.5f, 0.0f, 0.5f );
        TEST( volume(sph) == 1.767145868f, "Volume of a sphere wrong!" );
        TEST( volume(box) == 3.0f, "Volume of a box wrong!" );
        TEST( volume(obo) == 1.0f, "Volume of a oriented box wrong!" );
        TEST( volume(the) == 0.942809042f, "Volume of a tetrahedron wrong!" );
        TEST( volume(tri) == 0.0f, "Volume of a triangle wrong!" );
        TEST( volume(dis) == 0.0f, "Volume of a disc wrong!" );
        TEST( volume(pla) == 0.0f, "Volume of a plane wrong!" );
        TEST( volume(dop) == 0.0f, "Volume of a DOP wrong!" );
        TEST( volume(ell) == 3.53429174f, "Volume of an ellipsoid wrong!" );
        TEST( volume(oel) == 3.53429174f, "Volume of ar oriented ellipsoid wrong!" );
        TEST( volume(ray) == 0.0f, "Volume of a ray wrong!" );
        TEST( volume(seg) == 0.0f, "Volume of a segment wrong!" );
        TEST( volume(con) == 1.047197551f, "Volume of a cone wrong!" );
        TEST( volume(cap) == 1.30899704f, "Volume of a capsule wrong!" );
        TEST( volume(fru) == 0.33333333f, "Volume of a frustum wrong!" );
        float vol0 = volume(fru2), vol1 = volume(fru3);
        TEST( vol0+vol1 == 0.33333333f, "Volume of frusta wrong!" );

        TEST( surface(sph) == 7.068583471f, "Surface of a sphere wrong!" );
        TEST( surface(box) == 13.0f, "Surface of a box wrong!" );
        TEST( surface(obo) == 6.0f, "Surface of a oriented box wrong!" );
        TEST( surface(the) == 6.92820323f, "Surface of a tetrahedron wrong!" );
        TEST( surface(tri) == 0.707106781f, "Surface of a triangle wrong!" );
        TEST( surface(dis) == PI, "Surface of a disc wrong!" );
        TEST( surface(pla) == std::numeric_limits<float>::infinity(), "Surface of a plane wrong!" );
        TEST( surface(dop) == std::numeric_limits<float>::infinity(), "Surface of a DOP wrong!" );
        TEST( approx(surface(ell), 12.0816f, 0.012f), "Surface approximation of an ellipsoid too far away!" );
        TEST( approx(surface(oel), 12.0816f, 0.012f), "Surface approximation of an oriented ellipsoid too far away!" );
        TEST( surface(ray) == 0.0f, "Surface of a ray wrong!" );
        TEST( surface(seg) == 0.0f, "Surface of a segment wrong!" );
        TEST( approx(surface(con), 7.584475592f), "Surface of a cone wrong!" );
        TEST( surface(cap) == 6.283185307f, "Surface of a capsule wrong!" );
        Vec3 a(1.0f, 0.0f, 1.0f), b(2.0f, 0.0f, 1.0f), c(2.0f, 1.0f, 1.0f), d(1.0f, 1.0f, 1.0f);
        float refA3_A5 = 0.5f * (len(cross(a, b)) + len(cross(c, d)));
        float refA4_A6 = 0.5f * (len(cross(b, c)) + len(cross(d, a)));
        TEST( surface(fru) == refA3_A5 + refA4_A6 + 1.0f, "Surface of a frustum wrong!" );

        TEST( center(sph) == Vec3(1.0f, 2.0f, 3.14159f), "Center of sphere wrong!" );
        TEST( center(box) == Vec3(1.5f, 1.75f, 2.0f), "Center of box wrong!" );
        TEST( center(obo) == Vec3(0.5f), "Center of oriented box wrong!" );
        TEST( center(the) == Vec3(0.0f), "Center of thetrahedron wrong!" );
        TEST( center(tri) == Vec3(1.0f/3.0f), "Center of trinagle wrong!" );
        TEST( center(dis) == Vec3(1.0f), "Center of disc wrong!" );
        TEST( center(ell) == Vec3(-1.0f, -0.5f, -0.5f), "Center of ellipsoid wrong!" );
        TEST( center(seg) == Vec3(1.0f, 0.0f, 0.0f), "Center of line segment wrong!" );
        TEST( center(con) == Vec3(0.75f, 0.0f, 0.0f), "Center of cone wrong!" );
        TEST( center(cap) == Vec3(0.0f, 0.5f, 0.0f), "Center of capsule wrong!" );
        TEST( center(fru) == Vec3(1.125f, 0.375f, 0.75f), "Center of frustum wrong!" );
        Vec3 c0 = center(fru2), c1 = center(fru3);
        TEST( approx((c0*vol0 + c1*vol1) / (vol0 + vol1), Vec3(1.125f, 0.375f, 0.75f)), "Center of frustum wrong!" );
    }

    // Test plane construction
    {
        Plane pla0( Vec3(1.0f, 0.0f, 0.0f), 0.0f );
        Plane pla1( Vec3(1.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f) );
        Plane pla2( Vec3(0.0f, -1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, 2.0f, 4.0f) );
        TEST( (pla0.n == pla1.n) && pla0.d == pla1.d, "The constructed planes (0,1) should be equal!" );
        TEST( (pla0.n == pla2.n) && pla0.d == pla2.d, "The constructed planes (0,2) should be equal!" );
    }

    // Test ellipsoid construction and intersection with a point
    {
        Box box( Vec3(-2.0f, -1.0f, -1.0f), Vec3(0.0f, 0.0f, 0.0f) );
        Ellipsoid ell0( Vec3(-1.0f, -0.5f, -0.5f), Vec3(1.732050808f, 0.866025404f, 0.866025404f) );
        Ellipsoid ell1( box );
        TEST( (ell0.center == ell1.center), "Center of ellipsoid wrong!" );
        TEST( (ell0.radii == ell1.radii), "Radii of ellipsoid wrong!" );
        TEST( intersects( box.min, ell0 ), "The bounding box min must be contained in the ellipsoid!" );
        TEST( intersects( box.max, ell0 ), "The bounding box max must be contained in the ellipsoid!" );
        TEST( !intersects( box.max + 1e-6f, ell0 ), "The point must be outside the ellipsoid!" );
    }

    // Test sphere construction
    {
        Box box( Vec3(-2.0f, -1.0f, -1.0f), Vec3(0.0f, 0.0f, 0.0f) );
        Sphere sph0( Vec3(-1.0f, -0.5f, -0.5f), 1.224744871f );
        Sphere sph1( box );
        Vec3 points[] = {Vec3(0.0f), Vec3(1.0f), Vec3(0.1f, 0.4f, 0.3f), Vec3(-0.3f, 0.7f, 1.4f)};
        Sphere sph2(points, 4);
        for(int i = 0; i < 4; ++i)
            TEST(distance(sph2, points[i]) <= 1e-6f, "Bounding sphere does not enclose all points!");
        TEST( (sph0.center == sph1.center), "Center of sphere wrong!" );
        TEST( sph0.radius == sph1.radius, "Radius of sphere wrong!" );
    }

    // Test box construction
    {
        Box box0( Vec3(-2.0f, -1.0f, -1.0f), Vec3(0.0f, 0.0f, 0.0f) );
        Box box1( Vec3(3.0f), Vec3(7.0f) );
        Ellipsoid ell( Vec3(-1.0f, -0.5f, -0.5f), Vec3(1.0f, 0.5f, 0.5f) );
        Triangle tri( Vec3(-2.0f, -0.5f, -0.5f), Vec3(0.0f, -1.0f, 0.0f), Vec3(0.0f, 0.0f, -1.0f) );
        Sphere sph( Vec3(5.0f), 2.0f );
        OBox obo0( Vec3(0.0f), Vec3(0.5f, 1.0f, 1.5f), Quaternion(Vec3(0.0f, 0.0f, 1.0f), PI/2) );
        OBox obo1( Vec3(1.0f), Vec3(0.5f, 0.5f, 0.5f), Quaternion(Vec3(0.0f, 0.0f, 1.0f), PI/4) );
        OBox obo2( Vec3(4.0f), Vec3(0.5f, 1.0f, 1.5f), Quaternion(normalize(Vec3(0.1f, 0.5f, 1.0f)), 1.0f) );
        Vec3 points[] = {Vec3(0.0f), Vec3(1.0f), Vec3(0.1f, 0.4f, 0.3f), Vec3(-0.3f, 0.7f, 1.4f)};
        Tetrahedron tet( Vec3(-2.1f, -0.6f, -0.5f), Vec3(0.0f, -1.0f, 0.0f), Vec3(0.0f, -0.1f, -1.0f), Vec3(-0.5f, -1.1f, 0.2f) );
        Box box2( ell );
        Box box3( tri );
        Box box4( sph );
        Box box5( box0, box1 );
        Box box6( obo0 );
        Box box7( obo1 );
        // Box box8( OBox(box4) ); // <- This is a function? VC12 bug or standard?
        Box box8 = Box(OBox(box4));
        Box box9( points, 4u );
        Box box10( tet );
        Box box11( obo2 );

        TEST( box0.min == box2.min, "Ellipsoid bounding min is wrong!" );
        TEST( box0.max == box2.max, "Ellipsoid bounding min is wrong!" );
        TEST( box0.min == box3.min, "Triangle bounding min is wrong!" );
        TEST( box0.max == box3.max, "Triangle bounding max is wrong!" );
        TEST( box1.min == box4.min, "Sphere bounding min is wrong!" );
        TEST( box1.max == box4.max, "Sphere bounding max is wrong!" );
        TEST( box5.min == Vec3(-2.0f, -1.0f, -1.0f), "Box bounding min is wrong!" );
        TEST( box5.max == Vec3(7.0f), "Box bounding max is wrong!" );
        TEST( box6.min == Vec3(-1.0f, -0.5f, -1.5f), "OBox bounding min is wrong!" );
        TEST( box6.max == Vec3(1.0f, 0.5f, 1.5f), "OBox bounding max is wrong!" );
        TEST( box7.min == Vec3(0.292893219f, 0.292893219f, 0.5f), "OBox2 bounding min is wrong!" );
        TEST( box7.max == Vec3(1.707106781f, 1.707106781f, 1.5f), "OBox2 bounding max is wrong!" );
        TEST( box8.min == box4.min, "Box(OBox(box4)) introduces a bias to min!" );
        TEST( box8.max == box4.max, "Box(OBox(box4)) introduces a bias to max!" );
        TEST( box9.min == Vec3(-0.3f, 0.0f, 0.0f), "Box min of point list wrong!" );
        TEST( box9.max == Vec3(1.0f, 1.0f, 1.4f), "Box max of point list wrong!" );
        TEST( box10.min == Vec3(-2.1f, -1.1f, -1.0f), "Bounding box min of tetrahedron wrong!" );
        TEST( box10.max == Vec3(0.0f, -0.1f, 0.2f), "Bounding box max of tetrahedron wrong!" );
        TEST( approx(box11.min, Vec3(2.37966728f, 2.82336259f, 2.21573496f)), "OBox3 bounding min is wrong!" );
        TEST( approx(box11.max, Vec3(5.62033272f, 5.17663765f, 5.78426504f)), "OBox3 bounding max is wrong!" );
    }

    // Test oriented bounding box
    {
        Vec3 poi0[] = {Vec3(0.0f), Vec3(1.0f), Vec3(2.0f), Vec3(0.5f), Vec3(3.5f), Vec3(4.0f)}; // Diagonal points
        Vec3 poi1[] = {Vec3(-1.0f, 0.0f, -1.0f), Vec3(1.0f, 0.0f, -1.0f), Vec3(-1.0f, 0.0f, 1.0f), Vec3(1.0f, 0.0f, 1.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f)}; // Random points
        std::vector<Vec3> poi2;
        Mat3x3 rot = rotation(1.0f, 0.3f, -2.0f);
        Box box0(Vec3(0.0f), Vec3(0.0f)); // Create a non rotated optimal box to check later
        for( int i = 0; i < 20; ++i )
        {
            Vec3 newPoi;
            random(newPoi);
            box0.min = min(box0.min, newPoi);
            box0.max = max(box0.max, newPoi);
            poi2.push_back(rot * newPoi);
        }
        OBox obo0( poi0, 6 );
        OBox obo1( poi1, 6 );
        OBox obo2( poi2.data(), 20 );
        OBox obo3( poi0+2, 1 );
        OBox obo4( poi1, 2 );
        OBox obo5( poi0, 3 );
        TEST( approx(volume(obo0), 0.0f), "Oriented box 0 not optimal!" );
        TEST( approx(volume(obo1), 4.0f), "Oriented box 1 not optimal!" );
        TEST( volume(obo2) <= volume(box0), "Oriented box 2 not optimal!" );
        TEST( obo3.halfSides == Vec3(0.0f), "Oriented box 3 not as expected!" );
        TEST( approx(obo4.halfSides.x, 1.0f), "Oriented box 4 not as expected!" );
        TEST( approx(sum(obo5.halfSides), sqrt(3.0f)), "Oriented box 5 not as expected!" );

        OBox obo6(Disc(Vec3(2.0f), Vec3(1.0f, 0.0f, 0.0f), 0.5f));
        TEST( approx(obo6.orientation, Quaternion(Vec3(0.0f, 1.0f, 0.0f), PI/2.0f)), "Oriented box 6 has a wrong orientation!" );

        Box box1(Vec3(1.0f, 2.0f, 3.0f), Vec3(2.0f, 4.0f, 6.0f));
        OBox obo7( Quaternion(normalize(Vec3(0.1f, 0.5f, 1.0f)), 1.0f), box1 );
        OBox obo8( rotation(normalize(Vec3(0.1f, 0.5f, 1.0f)), 1.0f), box1 );
        TEST( obo7.center == Vec3(1.5f, 3.0f, 4.5f), "Oriented box 7 from box1 wrong center!" );
        TEST( obo7.halfSides == Vec3(1.54736483f, 1.383288145f, 1.670820475f), "Oriented box 7 from box1 wrong side length!" );
        TEST( obo8.center == Vec3(1.5f, 3.0f, 4.5f), "Oriented box 8 from box1 wrong center!" );
        TEST( approx(obo8.halfSides, Vec3(1.54736483f, 1.383288145f, 1.670820475f)), "Oriented box 8 from box1 wrong side length!" );
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
        Ray ray0(poi1, Vec3(0.0f, 0.707106769f, 0.707106769f));
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
        TEST( distance(poi3, box0) == -0.25f, "Distance between poi3 and box0 is -0.25!");
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
        TEST( distance(poi0, ray0) == 1.0f, "Distance between poi0 and ray0 is 1.0!");
        TEST( distance(poi1, ray0) == 0.0f, "Distance between poi1 and ray0 is 0.0!");
        TEST( approx(distance(poi2, ray0), 0.0f), "Distance between poi2 and ray0 is 0.0!");
        TEST( distance(poi4, ray0) == 1.791647287f, "Distance between poi4 and ray0 is 1.791647287f!");

        performance<Vec3,Triangle,float>(distance, "distance");
    }

    // ********************************************************************* //
    // Test convexSet algorithm
    {
        Vec3 a0[1] = {Vec3(0.0f)};
        TEST( convexSet(a0, 1) == 1, "Convex set of a0 wrong!" );
        Vec3 a1[2] = {Vec3(0.0f), Vec3(1.0f)};
        TEST( convexSet(a1, 2) == 2, "Convex set of a1 wrong!" );
        Vec3 a2[2] = {Vec3(0.0f), Vec3(0.0f)};
        TEST( convexSet(a2, 2) == 1, "Convex set of a2 wrong!" );
        Vec3 a3[3] = {Vec3(0.0f), Vec3(1.0f), Vec3(2.0f)};
        TEST( convexSet(a3, 3) == 2, "Convex set of a3 wrong!" );
        TEST( a3[1] == Vec3(2.0f), "Convex set of a3 wrong!" );
        Vec3 a4[5] = {Vec3(0.0f), Vec3(1.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(0.1f, 0.5f, 0.6f), Vec3(0.0f, 0.0f, 1.0f)};
        TEST( convexSet(a4, 5) == 4, "Convex set of a4 wrong!" );
        TEST( a4[3] == Vec3(0.0f, 0.0f, 1.0f), "Convex set of a4 wrong!" );
        Vec3 a5[5] = {Vec3(1.0f), Vec3(0.0f), Vec3(2.0f), Vec3(-1.0f), Vec3(1.5f)};
        TEST( convexSet(a5, 5) == 2, "Convex set of a5 wrong!" );
        Vec3 a6[5] = {Vec3(1.0f, 1.0f, 0.0f), Vec3(0.0f), Vec3(0.5f, 0.75f, 0.0f), Vec3(0.5f, 0.25f, 0.0f), Vec3(0.5f, 0.9f, 0.0f)};
        TEST( convexSet(a6, 5) == 4, "Convex set of a6 wrong!" );
    }

    return result;
}