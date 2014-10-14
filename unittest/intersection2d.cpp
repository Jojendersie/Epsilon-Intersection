#include "ei/2dintersection.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace ei;
using namespace std;

bool test_2dintersections()
{
    bool result = true;

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

    return result;
}