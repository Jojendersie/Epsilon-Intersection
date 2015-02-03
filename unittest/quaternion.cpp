#include "ei/vector.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace ei;

bool test_quaternion()
{
    bool result = true;

    Quaternion q0( normalize(Vec3(1.0f, 2.0f, -1.0f)), PI );
    Quaternion q1( 0.3f, 0.5f, -0.7f );
    Quaternion q3( normalize(Vec3(-1.0f, -2.0f, 1.0f)), PI );
    Mat3x3 m0 = rotation( normalize(Vec3(1.0f, 2.0f, -1.0f)), PI );
    Mat3x3 m1 = rotation( 0.3f, 0.5f, -0.7f );

    TEST( approx(Mat3x3(q0), m0), "Matrix(q0) does not equals m0 as expected!" );
    TEST( approx(q0, Quaternion(m0)), "q0 does not equals Quaternion(m0) as expected!" );
    TEST( approx(Mat3x3(q1), m1), "Matrix(q1) does not equals m1 as expected!" );
    TEST( approx(q1, Quaternion(m1)), "q1 does not equals Quaternion(m1) as expected!" );
    TEST( approx(axis(q0), normalize(Vec3(1.0f, 2.0f, -1.0f))), "Axis of quaternion q0 wrongly extracted!" );
    TEST( approx(axis(q3), normalize(Vec3(-1.0f, -2.0f, 1.0f))), "Axis of quaternion q3 wrongly extracted!" );
    TEST( approx(angle(q0), PI), "Angle of quaternion q0 wrongly extracted!" );
    TEST( approx(angle(q3), PI), "Angle of quaternion q3 wrongly extracted!" );

    TEST( approx(dot(q0, q0), 1.0f), "Quaternion dot product invalid" );
    TEST( approx(len(q0), 1.0f), "Quaternion len invalid" );
    q0 *= 2.0f;
    TEST( approx(len(q0), 2.0f), "Quaternion len after multiplication invalid" );
    q0 = normalize(q0);
    TEST( approx(len(q0), 1.0f), "Quaternion len after normalization invalid" );
    Quaternion q2 = conjugate(q0);
    TEST( q0.r == q2.r && q0.i == -q2.i && q0.j == -q2.j && q0.k == -q2.k, "Conjugated quaternion wrong!" );

    return result;
}
