#include "ei/3dintersection.hpp"
#include "unittest.hpp"
#include "performance3d.hpp"

using namespace ei;
using namespace std;


bool test_3dintersections()
{
    bool result = true;

    // Test box <-> box intersection
    {
        Box box0( Vec3(0.0f), Vec3(1.0f) );
        Box box1( Vec3(0.25f), Vec3(0.5f) );
        Box box2( Vec3(0.5f), Vec3(1.5f) );
        Box box3( Vec3(-1.0f), Vec3(0.0f) );
        Box box4( Vec3(0.4f, -1.0f, 0.4f), Vec3(0.6f, 2.0f, 0.6f) );
        Box box5( Vec3(2.0f, -1.0f, 0.0f), Vec3(3.0f, 1.0f, 1.0f) );
        TEST( intersects( box0, box1 ), "box1 is inside box0!" );
        TEST( intersects( box0, box2 ), "box2 intersects box0!" );
        TEST( intersects( box0, box3 ), "box3 touches box0!" );
        TEST( intersects( box0, box4 ), "box4 intersects box0!" );
        TEST( !intersects( box0, box5 ), "box5 does not intersect box0!" );

        performance<Box,Box>(intersects, "intersects");
    }

    // Test sphere <-> sphere intersection
    {
        Sphere sph0( Vec3(0.0f, 1.0f, 0.0f), 1.0f );
        Sphere sph1( Vec3(0.0f, 1.1f, 0.0f), 0.5f );
        Sphere sph2( Vec3(0.0f, -1.0f, 0.0f), 1.0f );
        Sphere sph3( Vec3(0.0f, -1.0f, 0.0f), 0.9f );
        TEST( intersects( sph0, sph1 ), "sph1 is inside sph0!" );
        TEST( intersects( sph0, sph2 ), "sph2 touches sph0!" );
        TEST( !intersects( sph0, sph3 ), "sph0 and sph3 do not intersect!" );

        performance<Sphere,Sphere>(intersects, "intersects");
    }

    // Test capsule <-> capsule intersection
    {
        Capsule cap0(Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 1.0f), 0.5f);
        Capsule cap1(Vec3(0.0f, 1.0f, 0.5f), Vec3(0.0f, 1.0f, 1.0f), 0.25f);
        Capsule cap2(Vec3(0.0f, 1.0f, 0.5f), Vec3(1.0f, 1.0f, 0.5f), 0.5f);
        TEST( intersects( cap1, cap2 ), "cap1 intersects cap2!" );
        TEST( intersects( cap0, cap2 ), "cap2 touches cap0!" );
        TEST( !intersects( cap0, cap1 ), "cap0 and cap1 do not intersect!" );

        performance<Capsule,Capsule>(intersects, "intersects");
    }

    // Test sphere <-> point intersection
    {
        Vec3 v0( 0.0f, 1.0f, 1.0f );
        Vec3 v1( 1e-3f, 1.0f, 1.0f );
        Sphere sph0( Vec3(0.0f, 1.0f, 0.0f), 1.0f );
        TEST( intersects( v0, sph0 ), "Point in sphere failed!" );
        TEST( !intersects( v1, sph0 ), "Point outside sphere failed!" );

        performance<Vec3,Sphere>(intersects, "intersects");
    }

    // Test point <-> box intersection
    {
        Box box0( Vec3(-0.25f), Vec3(-0.01f) );
        Vec3 poi0( 0.0f, 0.0f, 0.0f );
        Vec3 poi1( -0.01f );
        Vec3 poi2( -0.05f, -0.1f, -0.15f );
        TEST( !intersects( poi0, box0 ), "poi0 outside box0!" );
        TEST( intersects( poi1, box0 ), "poi1 touches box0!" );
        TEST( intersects( poi2, box0 ), "poi2 inside box0!" );

        performance<Vec3,Box>(intersects, "intersects");
    }

    // Test point <-> obox intersection
    {
        OBox box0( Vec3(1.0f), Vec3(0.125f, 0.5f, 0.75f), Quaternion(Vec3(1.0f, 0.0f, 0.0f), PI/2.0f) );
        Vec3 poi0( Vec3(0.85f, 1.0f, 1.0f) );
        Vec3 poi1( Vec3(0.875f, 1.25f, 0.5f) );
        Vec3 poi2( Vec3(0.9f, 1.1f, 0.9f) );
        Vec3 poi3( Vec3(1.0f, 1.25f, 1.76f) );
        TEST( !intersects( poi0, box0 ), "poi0 outside box0!" );
        TEST( intersects( poi1, box0 ), "poi1 touches box0!" );
        TEST( intersects( poi2, box0 ), "poi2 inside box0!" );
        TEST( !intersects( poi3, box0 ), "poi3 outside box0!" );

        performance<Vec3,OBox>(intersects, "intersects");
    }

    // Test point <-> dop intersection
    {
        DOP dop0( Vec3(1.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f), Vec3(1.5f, 0.0f, 0.0f) );
        Vec3 poi0( 0.0f, 1.2f, 0.0f );
        Vec3 poi1( 1.0f, 0.0f, 0.0f );
        Vec3 poi2( 1.35f, 0.0f, 0.0f );
        TEST( !intersects( poi0, dop0 ), "poi0 outside dop0!" );
        TEST( intersects( poi1, dop0 ), "poi1 touches dop0!" );
        TEST( intersects( poi2, dop0 ), "poi2 inside dop0!" );

        performance<Vec3,DOP>(intersects, "intersects");
    }

    // Test point <-> tetrahedron intersection
    {
        Tetrahedron tet0( Vec3(1.0f, 1.0f, 1.0f), Vec3(2.0f, 1.0f, 1.0f), Vec3(1.0f, 2.0f, 1.0f), Vec3(1.0f, 1.0f, 2.0f) );
        Vec3 poi0( 0.0f, 0.0f, 0.0f );
        Vec3 poi1( 1.5f, 1.0f, 1.0f );
        Vec3 poi2( 1.1f, 1.2f, 1.3f );
        Vec3 poi3( 1.34f, 1.34f, 1.34f );
        TEST( !intersects( poi0, tet0 ), "poi0 must be outside tet0!" );
        TEST( intersects( poi1, tet0 ), "poi1 must touche tet0!" );
        TEST( intersects( poi2, tet0 ), "poi2 must be inside tet0!" );
        TEST( !intersects( poi3, tet0 ), "poi3 must be outside tet0!" );

        performance<Vec3,Tetrahedron>(intersects, "intersects");
    }

    // Test sphere <-> box intersection
    {
        Sphere sph0( Vec3(-6.0f, -4.0f, 100.0f), 3.0f ); // intersects pla1
        Sphere sph1( Vec3(1.0f, 0.0f, 0.0f), 1.1f );     // intersects pla0
        Plane pla0( normalize(Vec3(1.0f, 0.5f, -0.25f)), 0.0f );
        Plane pla1( normalize(Vec3(1.0f, 1.0f, 0.0f)), 5.0f );
        TEST( !intersects( sph0, pla0 ), "sph0 and pl0 do not intersect!" );
        TEST( intersects( sph0, pla1 ), "sph0 and pla1 intersect!" );
        TEST( intersects( sph1, pla0 ), "sph1 and pl0 intersect!" );
        TEST( !intersects( sph1, pla1 ), "sph1 and pla1 do not intersect!" );

        performance<Sphere, Plane>(intersects, "intersects");
    }

    // Test sphere <-> plane intersection
    {
        Sphere sph0( Vec3(0.0f, 1.0f, 0.0f), 1.0f );
        Sphere sph1( Vec3(2.0f, 2.0f, 2.0f), 0.99f );
        Box box0( Vec3(0.0f), Vec3(1.0f) );
        TEST( intersects( sph0, box0 ), "Box sphere intersection failed!" );
        TEST( !intersects( sph1, box0 ), "sph1 is outside box0 failed!" );

        performance<Sphere, Box>(intersects, "intersects");
    }

    // Test ellipsoid <-> point intersection
    {
        Vec3 v0( 0.0f, 1.0f, 2.0f );
        Vec3 v1( 1e-8f, 1.0f, 2.0f );
        Ellipsoid ell0( Vec3(0.0f, 1.0f, 0.0f), Vec3(0.1f, 0.0f, 2.0f) );
        TEST( intersects( v0, ell0 ), "Point in degenerated ellipsoid failed!" );
        TEST( !intersects( v1, ell0 ), "Point outside degenerated ellipsoid failed!" );

        performance<Vec3,Ellipsoid>(intersects, "intersects");
    }

    // Test oellipsoid <-> point intersection
    {
        Vec3 v0( 0.0f, 2.0f, 1.0f );
        Vec3 v1( 0.0f, 1.0f, 0.5f );
        OEllipsoid oel0( Vec3(0.0f, 1.0f, 0.0f), Vec3(0.1f, 0.1f, 2.5f), Quaternion(PI/4.0f, 0.0f, 0.0f) );
        TEST( intersects( v0, oel0 ), "Point in oriented ellipsoid failed!" );
        TEST( !intersects( v1, oel0 ), "Point outside oriented ellipsoid failed!" );

        performance<Vec3,OEllipsoid>(intersects, "intersects");
    }

    // Test ray <-> plane intersection
    {
        Plane pla0( Vec3(1.0f, 0.0f, 0.0f), 0.0f );
        Plane pla1( Vec3(0.0f, 1.0f, 0.0f), 1.5f );
        Plane pla2( normalize(Vec3(1.0f, -2.0f, 43.0f)), 82935.4f );
        Ray ray0( Vec3(-1.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f) );
        Ray ray1( Vec3(0.0f, 0.0f, 0.0f), normalize(Vec3(0.0f, 1.0f, 1.0f)) );
        Ray ray2( Vec3(-1.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0f) );
        Ray ray3( Vec3(9.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0f) );
        float t;
        TEST( intersects(ray0, pla0, t) && t==1.0f, "ray0 should hit pla0!" );
        TEST( intersects(ray2, pla0, t) && t==-1.0f, "ray2 should hit pla0!" );
        TEST( intersects(ray3, pla0, t) && t==9.0f, "ray3 should hit pla0!" );
        TEST( !intersects(ray0, pla1, t), "ray0 should miss pla1!" );
        TEST( intersects(ray1, pla0, t) && t==0.0f, "ray1 should hit pla0!" );
        TEST( intersects(ray0, pla2, t) && t==-3571039.75f, "ray0 should hit pla2!" );
        TEST( intersects(ray1, pla2, t) && t==-123175.961f, "ray1 should hit pla2!" );
    }

    // Test sphere <-> ray intersection
    {
        Ray ray0( Vec3(0.0f, 1.0f, 2.0f), normalize(Vec3(0.2f, -0.9f, 0.0f)) ); // Starts inside sph1 and misses sph0
        Ray ray1( Vec3(100000.0f, 50.0f, -256.f), normalize(Vec3(-100000.0f, -50.0f, 256.0f)) ); /// Starting far away and intersecting with the origin, misses sph1
        Ray ray2( Vec3(-5.0f, 10.0f, 23.f), normalize(Vec3(5.0f, -9.0f, -22.0f)) ); // Targets center of sph1
        Ray ray3( Vec3(25.0f, -40.0f, 30.f), normalize(Vec3(5.0f, 1.0f, -2.0f)) );  // Should not hit
        Ray ray4( Vec3(0.0f, 2.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f) );     // sph0 and sph1 are behind the ray
        Sphere sph0( Vec3(0.0f, 0.0f, 0.0f), 1.0f );
        Sphere sph1( Vec3(0.0f, 1.0f, 1.0f), 1.1f );
        TEST( !intersects( ray0, sph0 ), "ray0 should miss sph0!" );
        TEST( intersects( ray0, sph1 ), "ray0 should hit sph1!" );
        TEST( intersects( ray1, sph0 ), "ray1 should hit sph0!" );
        TEST( !intersects( ray1, sph1 ), "ray1 should miss sph1!" );
        TEST( intersects( ray2, sph0 ), "ray2 should hit sph0!" );
        TEST( intersects( ray2, sph1 ), "ray2 should hit sph1!" );
        TEST( !intersects( ray3, sph0 ), "ray3 should miss sph0!" );
        TEST( !intersects( ray3, sph1 ), "ray3 should miss sph1!" );
        TEST( !intersects( ray4, sph0 ), "ray4 should miss sph0!" );
        TEST( !intersects( ray4, sph1 ), "ray4 should miss sph1!" );

        float d, d2;
        TEST( !intersects( ray0, sph0, d, d2 ), "2: ray0 should miss sph0!" );
        TEST( intersects( ray0, sph1, d, d2 ) && d==0.458257616f && d2==-0.458257616f, "ray0 should hit sph1 in a distance of 0!" );
        TEST( intersects( ray1, sph0, d, d2 ) && approx(d, 99999.3359f) && approx(d2, 100001.344f), "ray1 should hit sph0 in a distance of 99999.3359f!" );
        TEST( !intersects( ray1, sph1, d, d2 ), "2: ray1 should miss sph1!" );
        TEST( intersects( ray2, sph0, d, d2 ) && d==24.7731895f && approx(d2, 26.3591423f), "ray2 should hit sph0 in a distance of 24.7731895!" );
        TEST( intersects( ray2, sph1, d, d2 ) && d==23.1899151f && approx(d2, 25.3899155f), "ray2 should hit sph1 in a distance of 23.1899151!" );
        TEST( !intersects( ray3, sph0, d, d2 ), "2: ray3 should miss sph0!" );
        TEST( !intersects( ray3, sph1, d, d2 ), "2: ray3 should miss sph1!" );
        TEST( !intersects( ray4, sph0, d, d2 ), "2: ray4 should miss sph0!" );

        performance<Ray, Sphere, bool>(intersects, "intersects");
    }

    // Test ellipsoid <-> ray intersection
    {
        Ray ray0( Vec3(0.0f, 1.0f, 2.0f), normalize(Vec3(0.2f, -0.9f, 0.0f)) ); /// Starts inside ell1 and misses ell0
        Ray ray1( Vec3(100000.0f, 50.0f, -256.f), normalize(Vec3(-100000.0f, -50.0f, 256.0f)) ); /// Starting far away and intersecting with the origin
        Ray ray2( Vec3(-5.0f, 10.0f, 23.f), normalize(Vec3(5.0f, -9.0f, -22.0f)) ); /// Targets center of ell1
        Ray ray3( Vec3(25.0f, -40.0f, 30.f), normalize(Vec3(5.0f, 1.0f, -2.0f)) );  /// Should not hit
        Ellipsoid ell0( Vec3(0.0f, 0.0f, 0.0f), Vec3(0.5f, 0.0f, 2.0f) );
        Ellipsoid ell1( Vec3(0.0f, 1.0f, 1.0f), Vec3(1.0f, 3.0f, 2.0f) );
        TEST( !intersects( ray0, ell0 ), "ray0 should miss ell0!" );
        TEST( intersects( ray0, ell1 ), "ray0 should hit ell1!" );
        TEST( intersects( ray0.origin, ell1 ), "ray0.origin should be in ell1!" );
        TEST( intersects( Vec3(0.0f, 0.0f, 0.0f), ell1 ), "0 should be in ell1!" );
        // The next one is numerical to unstable
        TEST( intersects( ray1, ell0 ), "ray1 should hit ell0!" );
        TEST( intersects( ray1, ell1 ), "ray1 should hit ell1!" );
        TEST( !intersects( ray2, ell0 ), "ray2 should miss ell0!" );
        TEST( intersects( ray2, ell1 ), "ray2 should hit ell1!" );
        TEST( !intersects( ray3, ell0 ), "ray3 should miss ell0!" );
        TEST( !intersects( ray3, ell1 ), "ray3 should miss ell1!" );

        float d;
        TEST( intersects( ray0, ell1, d ), "ray0 should hit ell1!" );
        //TEST( intersects( ray0.origin + ray0.direction * (d-1e-6f), ell1 ), "Hit point outside ell1!" );
        //TEST( !intersects( ray0.origin + ray0.direction * (d+1e-6f), ell1 ), "Point should be outside ell1!" );
        TEST( intersects( ray2, ell1, d ), "ray2 should hit ell1!" );
        TEST( intersects( ray2.origin + ray2.direction * (d+1e-6f), ell1 ), "Hit point outside ell1!" );
        TEST( !intersects( ray2.origin + ray2.direction * (d-1e-6f), ell1 ), "Point should be outside ell1!" );

        performance<Ray,Ellipsoid,bool>(intersects, "intersects");
        performance<Ray,Ellipsoid,float,bool>(intersects, "intersects");
    }

    // Test box <-> ray intersection
    {
        Ray ray0( Vec3(-1.0f, 0.0f, 0.0f), normalize(Vec3(0.5f, 0.5f, 0.0f)) );
        Ray ray1( Vec3(3.0f, 5.0f, 8.0f), normalize(Vec3(-2.5f, -4.5f, -7.5f)) );
        Ray ray2( Vec3(1.25f, 1.25f, 1.25f), normalize(Vec3(1.5f, 0.5f, 0.5f)) );
        Ray ray3( Vec3(1.25f, 1.25f, 1.25f), normalize(Vec3(0.0f, 1.0f, 0.0f)) );
        Box box0( Vec3(0.0f, 0.0f, 0.0f), Vec3(0.5f, 1.0f, 2.0f) );
        Box box1( Vec3(1.0f, 1.0f, 1.0f), Vec3(1.5f, 1.5f, 1.5f) );
        TEST( intersects( ray0, box0 ), "ray0 should hit box0!" );
        TEST( !intersects( ray0, box1 ), "ray0 should miss box1!" );
        TEST( intersects( ray1, box0 ), "ray1 should hit box0!" );
        TEST( !intersects( ray1, box1 ), "ray1 should miss box1!" );
        TEST( !intersects( ray2, box0 ), "ray2 should miss box0!" );
        TEST( intersects( ray2, box1 ), "ray2 should hit box1!" );
        TEST( intersects( ray3, box1 ), "ray3 should hit box1!" );

        float d;
        TEST( intersects( ray3, box1, d ) && d == 0.0f, "ray3 should hit box1 with a distance of 0.0!" );
        TEST( intersects( ray0, box0, d ) && d == sqrt(2.0f), "ray0 should hit box0 with a distance of sqrt(2)!" );

        HitSide s;
        TEST( intersects( ray3, box1, d, s ) && s == HitSide::Y_POS, "ray3 should hit box1 at the top face!" );
        TEST( intersects( ray1, box0, d, s ) && s == HitSide::X_POS, "ray3 should hit box1 at the left face!" );

        performance<Ray,Box,bool>(intersects, "intersects");
        performance<Ray,Box,float,bool>(intersects, "intersects");
    }

    // Test obox <-> ray intersection
    {
        OBox obo0( Vec3(0.5f, 1.0f, 2.0f), Vec3(0.5f, 1.0f, 1.5f), Quaternion(Vec3(0.0f, 0.0f, 1.0f), PI/4.0f) );
        Ray ray0( Vec3(0.1f), Vec3(0.0f, 0.0f, 1.0f) );
        Ray ray1( Vec3(0.8f, 0.1f, 0.0f), Vec3(0.0f, 0.0f, 1.0f));
        Ray ray2( Vec3(2.5f, 0.0f, 2.0f), normalize(Vec3(-1.0f, -1.0f, 0.0f)));
        TEST( !intersects( ray0, obo0 ), "ray0 should miss obo0!" );
        TEST( intersects( ray1, obo0 ), "ray1 should hit obo0!" );
        TEST( !intersects( ray2, obo0 ), "ray2 should miss obo0!" );

        performance<OBox,Ray,bool>(intersects, "intersects");
    }

    // Test triangle <-> ray intersection
    {
        Ray ray0( Vec3(-1.0f, 0.0f, 0.0f), normalize(Vec3(0.5f, 0.5f, 0.0f)) );
        Ray ray1( Vec3(90.0f, 100.0f, -110.0f), normalize(Vec3(-88.75f, -99.5f, 111.16666f)) );
        Ray ray2( Vec3(1.0f, 1.5f, -0.9f), Vec3(1.0f, 0.0f, 0.0f));
        Triangle tri0( Vec3(0.0f, 0.0f, -1.0f), Vec3(0.0f, 2.0f, 1.0f), Vec3(0.0f, 2.0f, -1.0f) );
        Triangle tri1( Vec3(1.0f, 0.0f, 1.0f), Vec3(1.5f, 0.5f, 1.0f), Vec3(1.25f, 1.0f, 1.5f) );
        FastTriangle ftr0(tri0);
        FastTriangle ftr1(tri1);
        TEST( intersects( ray0, tri0 ), "ray0 should hit tri0!" );
        TEST( !intersects( ray0, tri1 ), "ray0 should miss tri1!" );
        TEST( !intersects( ray1, tri0 ), "ray1 should miss tri0!" );
        TEST( intersects( ray1, tri1 ), "ray1 should hit tri1!" );
        TEST( !intersects( ray2, tri0 ), "ray2 should miss tri0!" );

        float d;
        TEST( intersects( ray0, tri0, d ) && d == sqrt(2.0f), "ray0 should hit tri0 with a distance of sqrt(2)!" );
        TEST( intersects( ray1, tri1, d ) && d == 173.593903f, "ray1 should hit tri1 with a distance of 173.593903!" );
        TEST( intersects( ray0, ftr0, d ) && d == sqrt(2.0f), "ray0 should hit ftr0 with a distance of sqrt(2)!" );
        TEST( intersects( ray1, ftr1, d ) && d == 173.593887f, "ray1 should hit ftr1 with a distance of 173.593903!" );
        Vec3 bary;
        TEST( intersects( ray0, tri0, d, bary ), "Barycentric version of ray0/tri0 intersection failed!" );
        TEST( approx(ray0.origin+ray0.direction*d, tri0.v0*bary.x + tri0.v1*bary.y + tri0.v2*bary.z), "The hit point from barycentric coordinates and ray parameter should be the same (ray0, tri0)." );
        TEST( intersects( ray0, ftr0, d, bary ), "Barycentric version of ray0/ftr0 intersection failed!" );
        TEST( approx(ray0.origin+ray0.direction*d, tri0.v0*bary.x + tri0.v1*bary.y + tri0.v2*bary.z), "The hit point from barycentric coordinates and ray parameter should be the same (ray0, ftr0)." );
        TEST( intersects( ray1, tri1, d, bary ), "Barycentric version of ray1/tri1 intersection failed!" );
        TEST( approx(ray1.origin+ray1.direction*d, tri1.v0*bary.x + tri1.v1*bary.y + tri1.v2*bary.z, 1e-5f), "The hit point from barycentric coordinates and ray parameter should be the same (ray1, tri1)." );
        TEST( intersects( ray1, ftr1, d, bary ), "Barycentric version of ray1/ftr1 intersection failed!" );
        TEST( approx(ray1.origin+ray1.direction*d, tri1.v0*bary.x + tri1.v1*bary.y + tri1.v2*bary.z, 3e-5f), "The hit point from barycentric coordinates and ray parameter should be the same (ray1, ftr1)." );

        performance<Ray,Triangle,bool>(intersects, "intersects");
        performance<Ray,Triangle,float,bool>(intersects, "intersects");
        performance<Ray,FastTriangle,float,bool>(intersects, "intersects");
    }

    // Test sphere <-> triangle intersection
    {
        Triangle tri( Vec3(-1.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f) );
        Sphere sph0( Vec3(0.0f, 0.0f, 0.0f), 1.0f );
        Sphere sph1( Vec3(0.0f, 1.1f, 0.0f), 0.1f );
        Sphere sph2( Vec3(1.0f, 1.0f, 0.0f), 0.5f );
        Sphere sph3( Vec3(0.0f, 0.0f, -3.0f), 1.0f );
        TEST( intersects( tri, sph0 ), "sph0 intersects tri!" );
        TEST( intersects( tri, sph1 ), "sph1 touches tri!" );
        TEST( !intersects( tri, sph2 ), "sph2 misses tri!" );
        TEST( !intersects( tri, sph3 ), "sph3 misses tri!" );

        performance<Sphere,Triangle>(intersects, "intersects");
    }

    // Test sphere <-> capsule intersection
    {
        Capsule cap( Vec3(-1.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f), 0.5f );
        Sphere sph0( Vec3(0.0f, 0.0f, 0.0f), 1.0f );
        Sphere sph1( Vec3(0.0f, 1.1f, 0.0f), 0.6f );
        Sphere sph2( Vec3(1.0f, 1.0f, 0.0f), 0.5f );
        Sphere sph3( Vec3(0.0f, 0.0f, -3.0f), 1.0f );
        TEST( intersects( cap, sph0 ), "sph0 intersects cap!" );
        TEST( intersects( cap, sph1 ), "sph1 touches cap!" );
        TEST( intersects( cap, sph2 ), "sph2 intersects cap!" );
        TEST( !intersects( cap, sph3 ), "sph3 misses cap!" );

        performance<Sphere,Capsule>(intersects, "intersects");
    }

    // Test point <-> capsule intersection
    {
        Capsule cap( Vec3(-1.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), 0.5f );
        Vec3 poi0( 0.0f, 0.0f, 0.0f );
        Vec3 poi1( -0.5f, 0.5f, 0.0f );
        Vec3 poi2( 0.0f, 1.5f, 0.0f );
        TEST( !intersects( poi0, cap ), "poi0 outside cap!" );
        TEST( intersects( poi1, cap ), "poi1 central inside cap!" );
        TEST( intersects( poi2, cap ), "poi2 touches cap!" );

        performance<Vec3,Capsule>(intersects, "intersects");
    }

    // Test point <-> frustum intersection
    {
        FastFrustum ffr0( Vec3(0.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 1.0f, 0.0f), 1.0f, 2.0f, 0.0f, 1.0f, 0.0f, 1.0f );
        FastFrustum ffr1( Vec3(1.0f, 2.0f, 3.0f), normalize(Vec3(1.0f, 0.0f, 1.0f)), Vec3(0.0f, 1.0f, 0.0f), -1.0f, 1.0f, -0.5f, 0.5f, 0.5f, 2.0f );
        Vec3 poi0( 0.0f, 0.0f, 0.0f );
        Vec3 poi1( 1.5f, 0.5f, 0.9f );
        Vec3 poi2( -1.0f, 1.5f, 0.0f );
        Vec3 poi3( 2.0f, 2.0f, 4.0f );
        Vec3 poi4( 2.414213562f, 2.5f, 4.414213562f );
        Vec3 poi5( 2.0f, 3.0f, 5.0f );
        TEST( intersects( poi0, ffr0 ), "poi0 touches fru0!" );
        TEST( intersects( poi1, ffr0 ), "poi1 inside fru0!" );
        TEST( !intersects( poi2, ffr0 ), "poi2 outside fru0!" );
        TEST( !intersects( poi0, ffr1 ), "poi0 outside fru1!" );
        TEST( !intersects( poi1, ffr1 ), "poi1 outside fru1!" );
        TEST( !intersects( poi2, ffr1 ), "poi2 outside fru1!" );
        TEST( intersects( poi3, ffr1 ), "poi3 inside fru1!" );
        TEST( intersects( poi4, ffr1 ), "poi4 touches fru1!" );
        TEST( !intersects( poi5, ffr1 ), "poi5 outside fru1!" );

        //performance<Vec3,Capsule>(intersects, "intersects");
    }

    // Test ray <-> frustum intersection
    {
        Frustum fru0( Vec3(0.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 1.0f, 0.0f), -1.0f, 0.9f, -0.5f, 1.5f, 0.1f, 1000.0f);
        Ray ray0( Vec3(0.0f), Vec3(0.0f, 0.0f, 1.0f) ); // From apex, hitting
        Ray ray1( Vec3(0.0f), normalize(Vec3(100.0f, 0.0f, 1.0f)) ); // From apex, missing
        Ray ray2( Vec3(-0.03f, 0.07f, 50.0f), normalize(Vec3(-0.4f, 2.0f, 1.0f)) ); // From inside, hitting
        Ray ray3( Vec3(-0.999f, 1.499f, 1.0f), Vec3(0.0f, 0.0f, 1.0f) ); // Hitting a corner
        Ray ray4( ray3.origin, -ray3.direction ); // Missing
        float tmin = 0.0f, tmax = 1e38f;
        TEST( intersects( ray0, fru0, tmin, tmax ) && tmin == 0.1f && tmax == 1000.0f, "ray0 intersects fru0!" );
        TEST( !intersects( ray1, fru0, tmin, tmax ), "ray1 outside fru0!" );
        tmin = 0.0f, tmax = 1e38f;
        TEST( intersects( ray2, fru0, tmin, tmax ) && tmin == 0.0f && tmax > 0.0f, "ray2 intersects fru0!" );
        tmin = 0.0f, tmax = 1e38f;
        TEST( intersects( ray3, fru0, tmin, tmax ) && approx(tmin, 998.3333f) && approx(tmax, 999.0f), "ray3 intersects fru0!" );
        TEST( !intersects( ray4, fru0, tmin, tmax ), "ray4 outside fru0!" );
        Frustum fru1( Vec3(-1.0f, 0.0f, 50.0f), normalize(Vec3(1.0f, 1.0f, 0.0f)), Vec3(0.0f, 0.0f, 1.0f), -PHYTAGORAS, PHYTAGORAS, -PHYTAGORAS, PHYTAGORAS, PHYTAGORAS, 2 * PHYTAGORAS);
        // Missing rays
        Ray rm0( Vec3(1.0f, 0.0f, 50.0f), normalize(Vec3(-1.0f, -0.01f, 0.0f)) );
        Ray rm1( Vec3(1.0f, 0.0f, 50.0f), normalize(Vec3(1.0f, 0.95f, 0.0f)) );
        Ray rm2( Vec3(1.5f, 1.6f, 50.0f), normalize(Vec3(-1.0f, 1.0f, 0.0f)) );
        Ray rm3( Vec3(0.0f, 0.0f, 50.0f), normalize(Vec3(-1.0f, 1.0f, 0.0f)) );
        Ray rm4( Vec3(1.5f, 0.5f, 40.0f), normalize(Vec3(0.2f, 1.0f, 0.0f)) );
        Ray rm5( Vec3(2.5f, 2.0f, 50.0f), normalize(Vec3(1.0f, 0.2f, 0.0f)) );
        tmin = 0.0f, tmax = 1000.0f; TEST( !intersects( rm0, fru1, tmin, tmax ), "rm0 <=> fru1!" );
        tmin = 0.0f, tmax = 1000.0f; TEST( !intersects( rm1, fru1, tmin, tmax ), "rm1 <=> fru1!" );
        tmin = 0.0f, tmax = 1000.0f; TEST( !intersects( rm2, fru1, tmin, tmax ), "rm2 <=> fru1!" );
        tmin = 0.0f, tmax = 1000.0f; TEST( !intersects( rm3, fru1, tmin, tmax ), "rm3 <=> fru1!" );
        tmin = 0.0f, tmax = 1000.0f; TEST( !intersects( rm4, fru1, tmin, tmax ), "rm4 <=> fru1!" );
        tmin = 0.0f, tmax = 1000.0f; TEST( !intersects( rm5, fru1, tmin, tmax ), "rm5 <=> fru1!" );
        tmin = -1000.0f, tmax = 1000.0f; TEST( intersects( rm5, fru1, tmin, tmax ), "rm5.2 <=> fru1!" );
        // Hitting rays
        Ray rh0( Vec3(1.0f, 1.0f, 50.0f), normalize(Vec3(-0.3f, -0.6f, 0.4f)) );
        Ray rh1( Vec3(0.0f, 0.0f, 50.0f), normalize(Vec3(0.0f, 1.0f, 0.0f)) );
        Ray rh2( Vec3(1.1f, 0.1f, 50.0f), normalize(Vec3(-1.0f, 1.0f, 0.01f)) );
        Ray rh3( Vec3(1.0f, 1.0f, 0.0f), normalize(Vec3(0.0f, 0.0f, 1.0f)) );
        tmin = 0.0f, tmax = 0.5f; TEST( intersects( rh0, fru1, tmin, tmax ), "rh0 <=> fru1!" );
        tmin = 0.0f, tmax = 100.0f; TEST( intersects( rh1, fru1, tmin, tmax ), "rh1 <=> fru1!" );
        tmin = 0.0f, tmax = 100.0f; TEST( intersects( rh2, fru1, tmin, tmax ), "rh2 <=> fru1!" );
        tmin = 0.0f, tmax = 100.0f; TEST( intersects( rh3, fru1, tmin, tmax ), "rh3 <=> fru1!" );
    }

    // Test sphere <-> frustum intersection
    {
        Frustum fru0( Vec3(0.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 1.0f, 0.0f), 1.0f, 2.0f, 0.0f, 1.0f, 0.0f, 1.0f );
        FastFrustum ffr0( fru0 );
        Sphere sph0( Vec3(0.0f, 0.0f, 0.0f), 0.1f );    // around origin
        Sphere sph1( Vec3(1.5f, 0.5f, 0.8f), 0.1f );    // full inside
        Sphere sph2( Vec3(-1.0f, 2.5f, 0.0f), 0.5f );
        Sphere sph3( Vec3(3.0f, 1.5f, 1.25f), 0.5f );   // this one intersects three planes and should fail in the conventional test
        Sphere sph4( Vec3(-1.5f, -0.5f, -1.0f), 1.5f ); // intersects five planes and is outside
        Sphere sph5( Vec3(-1.5f, -0.5f, -1.0f), 2.0f ); // contains origin
        Sphere sph6( Vec3(0.0f, 0.25f, 0.5f), 0.5f );   // intersects left plane
        Sphere sph7( Vec3(0.0f, 1.0f, 0.5f), 0.75f );   // intersects upper left edge
        Sphere sph8( Vec3(0.0f, 1.0f, 0.5f), 0.7f );    // intersects two planes but is outside
        TEST( intersects( sph0, ffr0 ), "sph0 intersects fru0!" );
        TEST( intersects( sph1, ffr0 ), "sph1 inside fru0!" );
        TEST( !intersects( sph2, ffr0 ), "sph2 outside fru0!" );
        TEST( !intersects( sph3, ffr0 ), "sph3 outside fru0!" );
        TEST( !intersects( sph4, ffr0 ), "sph4 outside fru0!" );
        TEST( intersects( sph5, ffr0 ), "sph5 intersects fru0!" );
        TEST( intersects( sph6, ffr0 ), "sph6 intersects fru0!" );
        TEST( intersects( sph7, ffr0 ), "sph7 intersects fru0!" );
        TEST( !intersects( sph8, ffr0 ), "sph8 outside fru0!" );

        //performance<Vec3,Capsule>(intersects, "intersects");
    }

    // Test box <-> frustum intersection
    {
        FastFrustum ffr0( Vec3(1.0f, 2.0f, 3.0f), normalize(Vec3(1.0f, 0.0f, 1.0f)), Vec3(0.0f, 1.0f, 0.0f), -1.0f, 1.0f, -0.5f, 0.5f, 0.5f, 2.0f );
        Box box0( Vec3(-0.5f), Vec3(0.5f) );
        Box box1( Vec3(-10.0f, 0.0f, -10.0f), Vec3(100.0f, 1.4f, 100.0f) );
        Box box2( Vec3(0.5f), Vec3(15.5f) );
        Box box3( Vec3(2.0f, 2.0f, 4.0f), Vec3(2.1f, 2.1f, 4.1f) );
        Box box4( Vec3(2.0f, 0.0f, 2.0f), Vec3(4.0f, 4.0f, 3.5f) );
        Box box5( Vec3(3.5f, 2.0f, 3.0f), Vec3(4.0f, 2.1f, 4.0f) );
        Box box6( Vec3(1.0f, 2.0f, 4.0f), Vec3(1.5f, 3.0f, 7.0f) );
        Box box7( Vec3(3.0f, 2.25f, 3.5f), Vec3(3.5f, 2.75f, 4.0f) );
        TEST( !intersects( box0, ffr0 ), "box0 outside fru0!" );
        TEST( !intersects( box1, ffr0 ), "box1 outside fru0!" );
        TEST( intersects( box2, ffr0 ), "fru0 inside box2!" );
        TEST( intersects( box3, ffr0 ), "box3 inside fru0!" );
        TEST( intersects( box4, ffr0 ), "box4 intersects fru0!" );
        TEST( !intersects( box5, ffr0 ), "box5 outside fru0!" );
        TEST( intersects( box6, ffr0 ), "box6 intersects fru0!" );
        TEST( intersects( box7, ffr0 ), "box7 intersects fru0!" );
    }

    // Test (oriented) box <-> triangle intersection
    {
        Box box0(Vec3(-0.5f), Vec3(0.5f));
        Box box1(Vec3(1.0f), Vec3(2.0f, 3.0f, 4.0f));
        Triangle tri0(Vec3(0.0f), Vec3(0.2f), Vec3(0.1f, -0.1f, 0.2f)); // Fully inside box0, outside box1
        Triangle tri1(Vec3(0.4f), Vec3(1.5f, 10.0f, 0.4f), Vec3(1.5f, 11.0f, 2.0f)); // One vertex in box0, outside box1
        Triangle tri2(Vec3(1.5f, 2.0f, 0.0f), Vec3(1.5f, 5.0f, 2.0f), Vec3(1.5f, 4.0f, -1.0f)); // Edge overlap outside box1
        Triangle tri3(Vec3(1.5f, 2.0f, 0.0f), Vec3(1.5f, 3.5f, 5.0f), Vec3(1.5f, 4.0f, -1.0f)); // Edge overlap intersecting box1

        TEST( intersects(box0, tri0), "tri0 inside box0!" );
        TEST( !intersects(box1, tri0), "tri0 outside box1!" );
        TEST( intersects(box0, tri1), "tri1 intersects box0!" );
        TEST( !intersects(box1, tri1), "tri1 outside box1!" );
        TEST( !intersects(box0, tri2), "tri2 outside box0!" );
        TEST( !intersects(box1, tri2), "tri2 outside box1!" );
        TEST( intersects(box1, tri3), "tri3 intersects box1!" );
        performance<Triangle,Box>(intersects, "intersects");

		OBox obo0(Vec3(0.0f, 0.5f, 0.0f), Vec3(0.5f), Quaternion{0.0f, 0.0f, 0.1f});
        Triangle tri4(Vec3(0.1f, 0.0f, -0.3f), Vec3(0.1f, 0.0f, 0.3f), Vec3(2.0f, 0.0f, 0.0f));
        Triangle tri5(Vec3(0.1f, 0.0f, -0.3f), Vec3(0.1f, 0.0f, 0.3f), Vec3(-2.0f, 0.0f, 0.0f));
        TEST( !intersects(obo0, tri4), "tri4 outside obo0!" );
        TEST( intersects(obo0, tri5), "tri5 intersects obo0!" );
        performance<Triangle,OBox>(intersects, "intersects");
    }

    // Test (oriented) box <-> plane intersection
    {
        Box box0(Vec3(-0.5f), Vec3(0.5f));
        Plane pla0(normalize(Vec3(1.0f, 0.0f, 0.0f)), 0.2f); // Intersects
        Plane pla1(normalize(Vec3(1.0f, 0.0f, 0.0f)), 0.5f); // Touches
        Plane pla2(normalize(Vec3(1.0f, 0.0f, 0.0f)), 1.0f); // No intersection
        Plane pla3(normalize(Vec3(1.0f, 1.0f, 0.0f)), -0.2f); // Intersects
        Plane pla4(normalize(Vec3(1.0f, 1.0f, 0.0f)), 0.72f); // No intersection
        Plane pla5(normalize(Vec3(1.0f, 1.0f, 1.0f)), -0.6f); // Intersects
        Plane pla6(normalize(Vec3(1.0f, 1.0f, 1.0f)), 1.5f); // No intersection
        TEST( intersects(pla0, box0), "pla0 intersects box0!" );
        TEST( intersects(pla1, box0), "pla1 touches box0!" );
        TEST( !intersects(pla2, box0), "pla2 outside box0!" );
        TEST( intersects(pla3, box0), "pla3 intersects box0!" );
        TEST( !intersects(pla4, box0), "pla4 outside box0!" );
        TEST( intersects(pla5, box0), "pla5 intersects box0!" );
        TEST( !intersects(pla6, box0), "pla6 outside box0!" );
        performance<Plane,Box>(intersects, "intersects");

        OBox obo0(Vec3(1.0f, 2.0f, 1.0f), Vec3(0.5f, 0.5f, 1.0f), Quaternion(0.0f, -PI/4, 0.0f));
        Plane pla7(normalize(Vec3(1.0f, 0.0f, 1.0f)), 0.0f);
        Plane pla8(normalize(Vec3(1.0f, 0.0f, 1.0f)), -1.0f);
        Plane pla9(normalize(Vec3(1.0f, 0.0f, 1.0f)), -1.95f);
        Plane pla10(normalize(Vec3(1.0f, 0.0f, -1.0f)), 0.0f);
        Plane pla11(normalize(Vec3(1.0f, 0.0f, -1.0f)), 1.1f);
        TEST( !intersects(pla7, obo0), "pla7 outside obo0!" );
        TEST( intersects(pla8, obo0), "pla8 intersects obo0!" );
        TEST( !intersects(pla9, obo0), "pla9 outside obo0!" );
        TEST( intersects(pla10, obo0), "pla10 intersects obo0!" );
        TEST( !intersects(pla11, obo0), "pla11 outside obo0!" );
        performance<Plane,OBox>(intersects, "intersects");
    }

    // Test Cone <-> Point intersection
    {
        Cone con0(Vec3(0.0f), Vec3(1.0f, 0.0f, 0.0f), 0.5f, 1000.1f);
        Vec3 poi0(-1.0f, 0.0f, 0.0f); // Outside
        Vec3 poi1(0.0f); // Origin
        Vec3 poi2(1000.0f, 0.0f, 0.0f); // Inside far away
        Vec3 poi3(2.0f, 1.0f, 0.0f); // On boundary
        Vec3 poi4(2000.0f, 1.0f, -1.0f); // Outside far away
        Cone con1(Vec3(-1.0f), normalize(poi3 - Vec3(-1.0f)), 0.5f, 10.0f);
        TEST( !intersects(con0, poi0), "poi0 should be outside con0!" );
        TEST( intersects(con0, poi1), "poi1 should intersect con0!" );
        TEST( intersects(con0, poi2), "poi2 should intersect con0!" );
        TEST( intersects(con0, poi3), "poi3 should intersect con0!" );
        TEST( !intersects(con0, poi4), "poi4 should be outside con0!" );
        TEST( intersects(con1, poi3), "poi3 should intersect con1!" );
        TEST( !intersects(con1, poi0), "poi0 should be outside con1!" );
        TEST( !intersects(con1, poi2), "poi2 should be outside con1!" );
        performance<Cone,Vec3>(intersects, "intersects");

        FastCone fco0(con0);
        FastCone fco1(con1);
        TEST( !intersects(fco0, poi0), "poi0 should be outside fast con0!" );
        TEST( intersects(fco0, poi1), "poi1 should intersect fast con0!" );
        TEST( intersects(fco0, poi2), "poi2 should intersect fast con0!" );
        TEST( intersects(fco0, poi3), "poi3 should intersect fast con0!" );
        TEST( !intersects(fco0, poi4), "poi4 should be outside fast con0!" );
        TEST( intersects(fco1, poi3), "poi3 should intersect fast con1!" );
        TEST( !intersects(fco1, poi0), "poi0 should be outside fast con1!" );
        TEST( !intersects(fco1, poi2), "poi2 should be outside fast con1!" );
        performance<FastCone,Vec3>(intersects, "intersects");
    }

    // Test Cone <-> Triangle intersection
    {
        Triangle tri0(Vec3(0.0f), Vec3(1.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f));
        Cone con0(Vec3(0.0f, 0.0f, -1.0f), Vec3(0.0f, 0.0f, 1.0f), 0.5f, 2.0f);
        Cone con1(Vec3(0.5f, -0.01f, -1.0f), Vec3(0.0f, 0.0f, 1.0f), 0.25f, 2.0f);
        Cone con2(Vec3(0.25f, 0.25f, -1.0f), Vec3(0.0f, 0.0f, 1.0f), 0.1f, 2.0f);
        Cone con3(Vec3(-0.25f, 0.25f, -1.0f), Vec3(0.0f, 0.0f, 1.0f), 0.1f, 2.0f);
        Cone con4(Vec3(0.25f, 0.25f, -1.0f), Vec3(0.0f, 0.0f, 1.0f), 0.1f, 0.5f);
        Cone con5(Vec3(2.0f, 0.5f, 0.0f), Vec3(-1.0f, 0.0f, 0.0f), 0.1f, 1.1f);
        Cone con6(Vec3(2.0f, 0.5f, 0.0f), Vec3(-1.0f, 0.0f, 0.0f), 0.1f, 3.5f);

        TEST( intersects(tri0, con0), "tri0 should intersect con0!" ); // Vertex inside
        TEST( intersects(tri0, con1), "tri0 should intersect con1!" ); // Edge crosses
        TEST( intersects(tri0, con2), "tri0 should intersect con2!" ); // Fully occluded
        TEST( !intersects(tri0, con3), "tri0 should be outside con3!" ); // Fully outside
        TEST( !intersects(tri0, con4), "tri0 should be outside con4!" ); // Cone to short
        TEST( !intersects(tri0, con5), "tri0 should be outside con5!" ); // Cone to short
        TEST( intersects(tri0, con6), "tri0 should intersect con6!" ); // Edge crosses perpendicular
        performance<Cone, Triangle>(intersects, "intersects");
    }

    return result;
}