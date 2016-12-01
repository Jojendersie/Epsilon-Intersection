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
    Segment2D seg0( Vec2(-1.0f, 0.0f), Vec2(1.0f, 1.0f) );
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

    TEST( Disc2D( rec0 ) == Disc2D( Vec2(0.0f, 0.0f), PHYTAGORAS ),
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

    TEST( Disc2D( seg0 ) == Disc2D( Vec2(0.0f, 0.5f), sqrt(1.25f) ),
        "Circumcircle of line not correct!"
    );

    TEST( Disc2D( cap0 ) == Disc2D( Vec2(0.0f, 0.0f), PHYTAGORAS + 0.5f ),
        "Circumcircle of triangle not correct!"
    );

    TEST( area(dis0) == PI, "Area of dis0 wrong!" );
    TEST( area(rec0) == 4.0f, "Area of rec0 wrong!" );
    TEST( area(rec1) == 1.0f, "Area of rec1 wrong!" );
    TEST( area(tri0) == 0.25f, "Area of tri0 wrong!" );
    TEST( area(ell0) == PI, "Area of ell0 wrong!" );
    TEST( area(ell1) == PI, "Area of ell1 wrong!" );
    TEST( area(seg0) == 0.0f, "Area of seg0 wrong!" );
    TEST( area(ray0) == 0.0f, "Area of ray0 wrong!" );
    TEST( area(cap0) == 3.613825288f, "Area of cap0 wrong!" );

    TEST( center(dis0) == Vec2(0.0f, 1.0f), "Center of dis0 wrong!" );
    TEST( center(rec0) == Vec2(0.0f), "Center of rec0 wrong!" );
    TEST( center(rec1) == Vec2(0.0f), "Center of rec1 wrong!" );
    TEST( center(tri0) == Vec2(2/3.0f, 1/3.0f), "Center of tri0 wrong!" );
    TEST( center(ell0) == Vec2(2.0f, 2.0f), "Center of ell0 wrong!" );
    TEST( center(ell1) == Vec2(2.0f, 2.0f), "Center of ell1 wrong!" );
    TEST( center(cap0) == Vec2(0.0f), "Center of cap0 wrong!" );

    return result;
}