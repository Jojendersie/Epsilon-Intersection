#include "ei/2dintersection.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace ei;
using namespace std;

bool test_2dintersections()
{
    bool result = true;

    // ********************************************************************* //
    // Test point <-> disc
    {
        Disc2D disc( Vec2(0.0f, 0.0f), 1.0f );

        TEST( intersects(disc, Vec2(0.5f, 0.0f)), "Disc-point intersection not detected." );
        TEST( !intersects(disc, Vec2(0.75f, 0.75f)), "Wrong disc-point intersection not detected." );
    }

    // ********************************************************************* //
    // Test point <-> capsule
    {
        Capsule2D cap( Vec2(1.0f, 0.0f), Vec2(2.0f, 1.0f), 0.5f );

        TEST( intersects(cap, Vec2(0.5f, 0.0f)), "Capsule-point 0 intersection not detected." );
        TEST( intersects(cap, Vec2(1.5f, 0.5f)), "Capsule-point 1 intersection not detected." );
        TEST( intersects(cap, Vec2(1.2f, -0.2f)), "Capsule-point 2 intersection not detected." );
        TEST( !intersects(cap, Vec2(0.0f, 0.25f)), "Wrong capsule-point 3 intersection not detected." );
        TEST( !intersects(cap, Vec2(3.0f, 2.0f)), "Wrong capsule-point 4 intersection not detected." );
        TEST( !intersects(cap, Vec2(1.5f, 1.25f)), "Wrong capsule-point 5 intersection not detected." );
    }

    // ********************************************************************* //
    // Test disc <-> disc
    {
        Disc2D dis0( Vec2(0.0f, 0.0f), 1.0f );
        Disc2D dis1( Vec2(1.5f, 0.0f), 1.0f );
        Disc2D dis2( Vec2(0.0f, 1.5f), 0.5f );
        Disc2D dis3( Vec2(0.0f, -0.5f), 0.25f );

        TEST( intersects(dis0, dis1), "dis0, dis1 intersection not detected!" );
        TEST( intersects(dis0, dis2), "dis0, dis2 intersection not detected!" );
        TEST( !intersects(dis1, dis2), "dis1, dis2 false intersection detected!" );

        Vec2 location;
        TEST( intersects(dis0, dis1, location), "dis0, dis1 intersection not detected!" );
        TEST( location == Vec2(0.75f, 0.0f), "dis0, dis1 central intersection point wrong!" );
        TEST( intersects(dis0, dis2, location), "dis0, dis2 intersection not detected!" );
        TEST( location == Vec2(0.0f, 1.0f), "dis0, dis2 central intersection point wrong!" );
        TEST( intersects(dis0, dis3, location), "dis0, dis3 intersection not detected!" );
        TEST( location == Vec2(0.0f, -0.5f), "dis0, dis3 central intersection point wrong!" );
        TEST( !intersects(dis1, dis2, location), "dis1, dis2 false intersection detected!" );
    }

    // ********************************************************************* //
    // Test rect <-> rect
    {
        Rect2D rec0( Vec2(0.0f, 0.0f), Vec2(1.0f, 1.0f) );
        Rect2D rec1( Vec2(0.5f, 0.5f), Vec2(0.75f, 0.75f) );    // fully inside rec0
        Rect2D rec2( Vec2(0.25f, -1.0f), Vec2(2.0f, 2.0f) );    // left side crosses rec0
        Rect2D rec3( Vec2(0.5f, 0.5f), Vec2(1.5f, 1.5f) );      // corner overlap
        Rect2D rec4( Vec2(1.0f, 1.0f), Vec2(2.0f, 2.0f) );      // corner touches
        Rect2D rec5( Vec2(2.0f, 2.0f), Vec2(3.0f, 3.0f) );      // no intersection in any projection
        Rect2D rec6( Vec2(-1.0f, 0.25f), Vec2(-0.5f, 0.75f) );  // no intersection but in projection

        TEST( intersects(rec0, rec1), "rec0, rec1 intersection not detected!" );
        TEST( intersects(rec0, rec2), "rec0, rec2 intersection not detected!" );
        TEST( intersects(rec0, rec3), "rec0, rec3 intersection not detected!" );
        TEST( intersects(rec0, rec4), "rec0, rec4 intersection not detected!" );
        TEST( !intersects(rec0, rec5), "rec0, rec5 false intersection detected!" );
        TEST( !intersects(rec0, rec6), "rec0, rec6 false intersection detected!" );
    }

    // ********************************************************************* //
    // Test rect <-> point
    {
        Rect2D rec0( Vec2(0.0f, 0.0f), Vec2(1.0f, 1.0f) );
        Vec2 p0(0.5f, 0.5f);                                    // inside
        Vec2 p1(1.0f, 0.75f);                                   // on boundary
        Vec2 p2(0.0f, 0.0f);                                    // on corner
        Vec2 p3(-5.0f, 10.0f);                                  // far away
        Vec2 p4(0.5f, 1.1f);                                    // overlap in one dimension

        TEST( intersects(rec0, p0), "rec0, p0 intersection not detected!" );
        TEST( intersects(rec0, p1), "rec0, p1 intersection not detected!" );
        TEST( intersects(rec0, p2), "rec0, p2 intersection not detected!" );
        TEST( !intersects(rec0, p3), "rec0, p3 false intersection detected!" );
        TEST( !intersects(rec0, p4), "rec0, p4 false intersection detected!" );
    }

    // ********************************************************************* //
    // Test rect <-> disc
    {
        Rect2D rec0( Vec2(0.0f, 0.0f), Vec2(1.0f, 1.0f) );
        Disc2D dic0( Vec2(0.5f, 0.5f), 0.25f );                 // fully inside rec0
        Disc2D dic1( Vec2(1.5f, -0.5f), 10.0f );                // rec0 inside
        Disc2D dic2( Vec2(-0.25f, 1.0f), 1.0f );                // intersect
        Disc2D dic3( Vec2(-1.0f, -1.0f), 1.0f );                // touches corner
        Disc2D dic4( Vec2(-1.0f, -1.0f), 0.5f );                // no intersection in any projection
        Disc2D dic5( Vec2(-1.0f, 0.5f), 0.5f );                 // no intersection but in projection

        TEST( intersects(rec0, dic0), "rec0, dic0 intersection not detected!" );
        TEST( intersects(rec0, dic1), "rec0, dic1 intersection not detected!" );
        TEST( intersects(rec0, dic2), "rec0, dic2 intersection not detected!" );
        TEST( intersects(rec0, dic3), "rec0, dic3 intersection not detected!" );
        TEST( !intersects(rec0, dic4), "rec0, dic4 false intersection detected!" );
        TEST( !intersects(rec0, dic5), "rec0, dic5 false intersection detected!" );
    }

    // ********************************************************************* //
    // Test line <-> line
    {
        Segment2D lin0( Vec2(0.0f, 0.0f), Vec2(1.0f, 0.0f) );
        Segment2D lin1( Vec2(0.0f, 1.0f), Vec2(1.0f, 1.0f) );      // parallel
        Segment2D lin2( Vec2(0.5f, 0.0f), Vec2(1.5f, 0.0f) );      // parallel and overlapping
        Segment2D lin3( Vec2(0.5f, 0.5f), Vec2(2.0f, 2.0f) );      // somewhere else (no intersection)
        Segment2D lin4( Vec2(0.0f, 1.0f), Vec2(0.0f, 0.0f) );      // perpendicular touching in an endpoint
        Segment2D lin5( Vec2(0.5f, 0.5f), Vec2(0.75f, -0.5f) );    // intersection

        TEST( !intersects(lin0, lin1), "lin0, lin1 false intersection detected!" );
        TEST( intersects(lin0, lin2), "lin0, lin2 intersection not detected!" );
        TEST( !intersects(lin0, lin3), "lin0, lin3 false intersection detected!" );
        TEST( intersects(lin0, lin4), "lin0, lin4 intersection not detected!" );
        TEST( intersects(lin0, lin5), "lin0, lin5 intersection not detected!" );
        TEST( intersects(lin5, lin0), "lin5, lin0 intersection not detected!" );

        Vec2 location;
        TEST( !intersects(lin0, lin1, location), "lin0, lin1 false intersection w.loc. detected!" );
        TEST( !intersects(lin0, lin3, location), "lin0, lin3 false intersection w.loc. detected!" );
        TEST( intersects(lin0, lin2, location), "lin0, lin2 intersection w.loc. not detected!" );
        TEST( location == Vec2(0.75f, 0.0f), "Central point of overlap wrong!" );
        TEST( intersects(lin0, lin4, location), "lin0, lin4 intersection w.loc. not detected!" );
        TEST( location == Vec2(0.0f, 0.0f), "Intersection point of lin0, lin4 wrong!" );
        TEST( intersects(lin0, lin5, location), "lin0, lin5 intersection w.loc. not detected!" );
        TEST( location == Vec2(0.625f, 0.0f), "Intersection point of lin0, lin5 wrong!" );
        TEST( intersects(lin5, lin0, location), "lin5, lin0 intersection w.loc. not detected!" );
        TEST( location == Vec2(0.625f, 0.0f), "Intersection point of lin5, lin0 wrong!" );
    }

    return result;
}