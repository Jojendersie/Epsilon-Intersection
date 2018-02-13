#include "ei/elementarytypes.hpp"
#include "unittest.hpp"

#include <iostream>
#include <limits>

using namespace ei;
using namespace std;

//#define UNIT_TEST_SUCCESSOR
//#define UNIT_TEST_PREDECESSOR

bool test_elementaries()
{
    bool result = true;

    // Test of type sizes
    TEST( static_cast<uint>(-1) == static_cast<unsigned int>(-1),
        "The type 'uint' is wrong defined. It should be equal to 'unsigned int'\n"
    );
    TEST( sizeof(uint8) == 1 && sizeof(int8) == 1 && sizeof(byte) == 1,
        "The 8bit integer types 'uint8', 'int8' and 'byte' are not one byte large\n"
    );
    TEST( sizeof(uint16) == 2 && sizeof(int16) == 2,
        "The 16-bit integer types 'uint16' and 'int16' are not two byte large\n"
    );
    TEST( sizeof(uint32) == 4 && sizeof(int32) == 4,
        "The 32-bit integer types 'uint32' and 'int32' are not four byte large\n"
    );
    TEST( sizeof(uint64) == 8 && sizeof(int64) == 8,
        "The 64-bit integer types 'uint64' and 'int64' are not eight byte large\n"
    );

    // Test if types are integral
    TEST( std::numeric_limits<uint8>::is_integer,
        "'uint8' is not an integer\n"
    );
    TEST( std::numeric_limits<int8>::is_integer,
        "'int8' is not an integer\n"
    );
    TEST( std::numeric_limits<byte>::is_integer,
        "'byte' is not an integer\n"
    );
    TEST( std::numeric_limits<uint16>::is_integer,
        "'uint16' is not an integer\n"
    );
    TEST( std::numeric_limits<int16>::is_integer,
        "'int16' is not an integer\n"
    );
    TEST( std::numeric_limits<uint32>::is_integer,
        "'uint32' is not an integer\n"
    );
    TEST( std::numeric_limits<int32>::is_integer,
        "'int32' is not an integer\n"
    );
    TEST( std::numeric_limits<uint64>::is_integer,
        "'uint64' is not an integer\n"
    );
    TEST( std::numeric_limits<int64>::is_integer,
        "'int64' is not an integer\n"
    );

    // frac, intfrac, floorfrac, mod
    {
        float x0 = 1.25f;
        float x1 = 3.0f;
        float x2 = -0.75f;
        float y0 = 1.0f;
        float y1 = 2.5f;
        float y2 = -1.5f;

        TEST( frac(x0) == 0.25f, "frac of x0 wrong!");
        TEST( frac(x1) == 0.0f, "frac of x1 wrong!");
        TEST( frac(x2) == -0.75f, "frac of x2 wrong!");
        int i;
        TEST( intfrac(x0, i) == 0.25f && i==1, "intfrac of x0 wrong!");
        TEST( intfrac(x1, i) == 0.0f && i==3, "intfrac of x1 wrong!");
        TEST( intfrac(x2, i) == -0.75f && i==0, "intfrac of x2 wrong!");

        TEST( floorfrac(x0, i) == 0.25f && i==1, "floorfrac of x0 wrong!");
        TEST( floorfrac(x1, i) == 0.0f && i==3, "floorfrac of x1 wrong!");
        TEST( floorfrac(x2, i) == 0.25f && i==-1, "floorfrac of x2 wrong!");

        TEST( mod(x0, y0) == 0.25f, "mod(x0, y0) wrong!");
        TEST( mod(x0, y2) == 1.25f, "mod(x0, y2) wrong!");
        TEST( mod(x1, y0) == 0.0f, "mod(x1, y0) wrong!");
        TEST( mod(x1, y1) == 0.5f, "mod(x1, y1) wrong!");
        TEST( mod(x1, y2) == 0.0f, "mod(x1, y2) wrong!");
        TEST( mod(x2, y0) == 0.25f, "mod(x2, y0) wrong!");
        TEST( mod(x2, y2) == 0.75f, "mod(x2, y2) wrong!");
        TEST( mod(7, 3) == 1, "mod(7, 3) wrong!");
        TEST( mod(7, -3) == 1, "mod(7, -3) wrong!");
        TEST( mod(-7, 3) == 2, "mod(-7, 3) wrong!");
        TEST( mod(-7, -3) == 2, "mod(-7, -3) wrong!");
        TEST( mod(4, 2) == 0, "mod(4, 2) wrong!");
        TEST( mod(-4, 2) == 0, "mod(-4, 2) wrong!");

        uint16 a = 7, b = 3;
        TEST( mod(a, b) == 1, "mod(a, b) wrong!");
    }

    // abs, sgn, sign
    {
        TEST( ei::sign(-1.0f) == -1, "sign(-1) wrong!" );
        TEST( ei::sign(1.0f) == 1, "sign(1) wrong!" );
        TEST( ei::sign(-0.0f) == 0, "sign(-0) wrong!" );
        TEST( ei::sign(0.0f) == 0, "sign(0) wrong!" );

        TEST( ei::abs(-1.0f) == 1.0f, "abs(-1.0f) wrong!" );
        TEST( ei::abs(1.0f) == 1.0f, "abs(1.0f) wrong!" );
        TEST( 1.0f/ei::abs(-0.0f) == INF, "abs(-0.0f) wrong!" );
        TEST( 1.0f/ei::abs(0.0f) == INF, "abs(0.0f) wrong!" );
        TEST( ei::abs(-1.0) == 1.0, "abs(-1.0) wrong!" );
        TEST( ei::abs(1.0) == 1.0, "abs(1.0) wrong!" );
        TEST( 1.0/ei::abs(-0.0) == INF_D, "abs(-0.0) wrong!" );
        TEST( 1.0/ei::abs(0.0) == INF_D, "abs(0.0) wrong!" );

        TEST( ei::sgn(-1.0f) == -1, "sgn(-1.0f) wrong!" );
        TEST( ei::sgn(1.0f) == 1, "sgn(1.0f) wrong!" );
        TEST( ei::sgn(-0.0f) == -1, "sgn(-0.0f) wrong!" );
        TEST( ei::sgn(0.0f) == 1, "sgn(0.0f) wrong!" );
        TEST( ei::sgn(-1.0) == -1, "sgn(-1.0) wrong!" );
        TEST( ei::sgn(1.0) == 1, "sgn(1.0) wrong!" );
        TEST( ei::sgn(-0.0) == -1, "sgn(-0.0) wrong!" );
        TEST( ei::sgn(0.0) == 1, "sgn(0.0) wrong!" );

        TEST( heaviside(-1.0f) == 0, "heaviside(-1.0f) wrong!" );
        TEST( heaviside(1.0f) == 1, "heaviside(1.0f) wrong!" );
        TEST( heaviside(-0.0f) == 0, "heaviside(-0.0f) wrong!" );
        TEST( heaviside(0.0f) == 1, "heaviside(0.0f) wrong!" );
        TEST( heaviside(-1.0) == 0, "heaviside(-1.0) wrong!" );
        TEST( heaviside(1.0) == 1, "heaviside(1.0) wrong!" );
        TEST( heaviside(-0.0) == 0, "heaviside(-0.0) wrong!" );
        TEST( heaviside(0.0) == 1, "heaviside(0.0) wrong!" );
    }

    // successor, predecessor
    {
#ifdef UNIT_TEST_SUCCESSOR
        float x = -INF;
        for(unsigned i=0; i < 2 * ((1u<<31)-(1u<<23)); ++i)
        {
            float p = successor(x);
            TEST( p > x, "successor(" << x << ") < " << p << "!" );
            x = p;
        }
#endif // UNIT_TEST_SUCCESSOR
//        TEST( successor(0.0f) == successor(-0.0f), "successor of 0 or -0 wrong!" );

#ifdef UNIT_TEST_PREDECESSOR
        float y = INF;
        for(unsigned i=0; i < 2 * ((1u<<31)-(1u<<23)); ++i)
        {
            float p = predecessor(x);
            TEST( p < y, "predecessor(" << y << ") > " << p << "!" );
            y = p;
        }
#endif // UNIT_TEST_PREDECESSOR

        TEST( successor(0.0) == successor(-0.0), "successor (double) of 0 or -0 wrong!" );
        TEST( successor(1.0f) > 1, "successor (double) of 1 wrong!" );
        TEST( successor(-1.0f) > -1, "successor (double) of -1 wrong!" );
        TEST( successor(INF_D) == INF_D, "successor (double) of INF_D wrong!" );
        TEST( successor(-INF_D) == -1.7976931348623157e+308, "successor (double) of -INF_D wrong!" );

        TEST( predecessor(0.0) == predecessor(-0.0), "predecessor (double) of 0 or -0 wrong!" );
        TEST( predecessor(1.0f) < 1, "predecessor (double) of 1 wrong!" );
        TEST( predecessor(-1.0f) < -1, "predecessor (double) of -1 wrong!" );
        TEST( predecessor(-INF_D) == -INF_D, "predecessor (double) of -INF_D wrong!" );
        TEST( predecessor(INF_D) == 1.7976931348623157e+308, "predecessor (double) of INF_D wrong!" );
    }

    return result;
}
