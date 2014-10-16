#include "ei/2dtypes.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace ei;
using namespace std;

bool test_2dtypes()
{
    bool result = true;

    // Primitives for different testes
    Disc2D dis0( Vec2(0.0f, 1.0f), 1.0f );
    Rect2D rec0( Vec2(-1.0f, -1.0f), Vec2(1.0f, 1.0f) );
    ORect2D rec1( Vec2(0.0f, 0.0f), Vec2(1.0f, 1.0f), 0.5f );
    Triangle2D tri0( Vec2(0.5f, 0.0f), Vec2(1.0f, 0.0f), Vec2(0.5f, 1.0f) );
    Ellipse2D ell0( Vec2(2.0f, 2.0f), Vec2(2.0f, 0.5f) );
    OEllipse2D ell1( Vec2(2.0f, 2.0f), Vec2(2.0f, 0.5f), -0.5f );
    Line2D lin0( Vec2(-1.0f, 0.0f), Vec2(1.0f, 1.0f) );
    Ray2D ray0( Vec2(0.0f, -5.0f), Vec2(0.0f, 1.0f) );
    Capsule2D cap0( Vec2(-1.0f, 1.0f), Vec2(1.0f, -1.0f), 0.5f );

    // ********************************************************************* //
    // Test conversion to circles
    TEST( Disc2D( Vec2(0.5f, 0.0f), Vec2(1.0f, 0.0f) ) == Disc2D( Vec2(0.75f, 0.0f), 0.25f ),
        "Circumcircle of two points not correct!"
    );

    TEST( Disc2D( Vec2(0.5f, 0.0f), Vec2(1.0f, 0.0f), Vec2(0.5f, 1.0f) ) == Disc2D( Vec2(0.75f, 0.5f), 0.5590170f ),
        "Circumcircle of three points not correct!"
    );

    TEST( Disc2D( rec0 ) == Disc2D( Vec2(0.0f, 0.0f), √2 ),
        "Circumcircle of rectangle not correct!"
    );

    TEST( Disc2D( rec1 ) == Disc2D( Vec2(0.0f, 0.0f), sqrt(0.5f) ),
        "Circumcircle of oriented rectangle not correct!"
    );

    TEST( Disc2D( tri0 ) == Disc2D( Vec2(0.75f, 0.5f), 0.5590170f ),
        "Circumcircle of triangle not correct!"
    );

    TEST( Disc2D( ell0 ) == Disc2D( Vec2(2.0f, 2.0f), 2.0f ),
        "Circumcircle of ellipse not correct!"
    );

    TEST( Disc2D( ell1 ) == Disc2D( Vec2(2.0f, 2.0f), 2.0f ),
        "Circumcircle of oriented ellipse not correct!"
    );

    TEST( Disc2D( lin0 ) == Disc2D( Vec2(0.0f, 0.5f), sqrt(1.25f) ),
        "Circumcircle of line not correct!"
    );

    TEST( Disc2D( cap0 ) == Disc2D( Vec2(0.0f, 0.0f), √2 + 0.5f ),
        "Circumcircle of triangle not correct!"
    );

    return result;
}