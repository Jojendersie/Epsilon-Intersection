﻿#include "ei/2dintersection.hpp"
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
        TEST( all(location == Vec2(0.75f, 0.0f)), "dis0, dis1 central intersection point wrong!" );
        TEST( intersects(dis0, dis2, location), "dis0, dis2 intersection not detected!" );
        TEST( all(location == Vec2(0.0f, 1.0f)), "dis0, dis2 central intersection point wrong!" );
        TEST( intersects(dis0, dis3, location), "dis0, dis3 intersection not detected!" );
        TEST( all(location == Vec2(0.0f, -0.5f)), "dis0, dis3 central intersection point wrong!" );
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
        Line2D lin0( Vec2(0.0f, 0.0f), Vec2(1.0f, 0.0f) );
        Line2D lin1( Vec2(0.0f, 1.0f), Vec2(1.0f, 1.0f) );      // parallel
        Line2D lin2( Vec2(0.5f, 0.0f), Vec2(1.5f, 0.0f) );      // parallel and overlapping
        Line2D lin3( Vec2(0.5f, 0.5f), Vec2(2.0f, 2.0f) );      // somewhere else (no intersection)
        Line2D lin4( Vec2(0.0f, 0.0f), Vec2(0.0f, 1.0f) );      // perpendicular touching in an endpoint
        Line2D lin5( Vec2(0.5f, 0.5f), Vec2(0.75f, -0.5f) );    // intersection

        TEST( !intersects(lin0, lin1), "lin0, lin1 false intersection detected!" );
        TEST( intersects(lin0, lin2), "lin0, lin2 intersection not detected!" );
        TEST( !intersects(lin0, lin3), "lin0, lin3 false intersection detected!" );
        TEST( intersects(lin0, lin4), "lin0, lin4 intersection not detected!" );
        TEST( intersects(lin0, lin5), "lin0, lin5 intersection not detected!" );
        TEST( intersects(lin5, lin0), "lin5, lin0 intersection not detected!" );
    }

    return result;
}