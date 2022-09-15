#include "ei/quaternion.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace ei;

bool test_quaternion()
{
    bool result = true;

    const Quaternion q0( normalize(Vec3(1.0f, 2.0f, -1.0f)), PI );
    const Quaternion q1( 0.3f, 0.5f, -0.7f );
    const Quaternion q2( normalize(Vec3(-1.0f, -2.0f, 1.0f)), PI );
    const Mat3x3 m0 = rotation( normalize(Vec3(1.0f, 2.0f, -1.0f)), PI );
    const Mat3x3 m1 = rotation( 0.3f, 0.5f, -0.7f );

    TEST( approx(Mat3x3(q0), m0), "Matrix(q0) does not equals m0 as expected!" );
    TEST( approx(q0, Quaternion(m0)), "q0 does not equals Quaternion(m0) as expected!" );
    TEST( approx(Mat3x3(q1), m1), "Matrix(q1) does not equals m1 as expected!" );
    TEST( approx(q1, Quaternion(m1)), "q1 does not equals Quaternion(m1) as expected!" );
    TEST( approx(axis(q0), normalize(Vec3(1.0f, 2.0f, -1.0f))), "Axis of quaternion q0 wrongly extracted!" );
    TEST( approx(axis(q2), normalize(Vec3(-1.0f, -2.0f, 1.0f))), "Axis of quaternion q3 wrongly extracted!" );
    TEST( approx(angle(q0), PI), "Angle of quaternion q0 wrongly extracted!" );
    TEST( approx(angle(q2), PI), "Angle of quaternion q3 wrongly extracted!" );
    TEST( approx(angles(q1), Vec3(0.3f, 0.5f, -0.7f)), "Euler angles conversion of quaternion wrong!" );
    TEST( approx(xRow(q1), Vec3(m1(0,0), m1(0,1), m1(0,2))), "X-axis-row of quaternion q1 wrong!" );
    TEST( approx(yRow(q1), Vec3(m1(1,0), m1(1,1), m1(1,2))), "Y-axis-row of quaternion q1 wrong!" );
    TEST( approx(zRow(q1), Vec3(m1(2,0), m1(2,1), m1(2,2))), "Z-axis-row of quaternion q1 wrong!" );
    TEST( approx(xCol(q1), Vec3(m1(0,0), m1(1,0), m1(2,0))), "X-axis-col of quaternion q1 wrong!" );
    TEST( approx(yCol(q1), Vec3(m1(0,1), m1(1,1), m1(2,1))), "Y-axis-col of quaternion q1 wrong!" );
    TEST( approx(zCol(q1), Vec3(m1(0,2), m1(1,2), m1(2,2))), "Z-axis-col of quaternion q1 wrong!" );

    TEST( approx(dot(q0, q0), 1.0f), "Quaternion dot product invalid" );
    TEST( approx(len(q0), 1.0f), "Quaternion len invalid" );
    Quaternion qtmp = q0 * 2.0f;
    TEST( approx(len(qtmp), 2.0f), "Quaternion len after multiplication invalid" );
    qtmp = normalize(qtmp);
    TEST( approx(len(qtmp), 1.0f), "Quaternion len after normalization invalid" );
    Quaternion qtmp2 = conjugate(qtmp);
    TEST( qtmp.r == qtmp2.r && qtmp.i == -qtmp2.i && qtmp.j == -qtmp2.j && qtmp.k == -qtmp2.k, "Conjugated quaternion wrong!" );
    qtmp2 = ~qtmp;
    TEST( qtmp.r == qtmp2.r && qtmp.i == -qtmp2.i && qtmp.j == -qtmp2.j && qtmp.k == -qtmp2.k, "Conjugated quaternion wrong!" );
    qtmp2 = -qtmp;
    TEST( qtmp.r == -qtmp2.r && qtmp.i == -qtmp2.i && qtmp.j == -qtmp2.j && qtmp.k == -qtmp2.k, "Unary quaternion minus wrong!" );
    TEST( qtmp*qidentity() == qtmp, "Quaternion identity (rhs) wrong!" );
    TEST( qidentity()*qtmp == qtmp, "Quaternion identity (lhs) wrong!" );

    const Quaternion q4( 0.1f, 0.5f, -0.7f );
    const Quaternion q5( 0.25f, 0.5f, -0.7f );
    // Unfortunatelly the slerp is very unprecise -> large epsilon
    TEST( approx(slerp(q1, q4, 0.25f), q5), "Slerp for quaternions failed!" );
    TEST( slerp(q1, q1, 0.33f) == q1, "Slerp (equal vectors) for quaternions failed!" );

    {// Direct transformation and quaternion chaining
        const Vec3 v0(1.0f, 0.0f, 0.0f);
        const Vec3 v1(0.0f, 1.0f, 0.0f);
        const Vec3 v2 = normalize(Vec3(1.0f, 2.0f, -1.0f));
        const Quaternion q6(v0, v1);
        const Quaternion q7(v2, v1);
        const Mat3x3 m2(q6);
        TEST( approx(m2 * v0, v1), "Quaternion From-To parametrization or matrix transformation invalid!" );
        TEST( approx(transform(v0, q6), v1), "Direct quaternion transformation v0->v1 failed!" );
        TEST( approx(transform(v2, q7), v1), "Direct quaternion transformation v2->v1 failed!" );

        const Quaternion q8 = q0 * q1;
        const Mat3x3 m8 = m0 * m1;
        TEST( approx(m8, Mat3x3(q8)), "Chained rotation wrong (quaternion multiplication)." );
        Vec3 v3 = m8 * v2;
        Vec3 v4 = transform(v2, q8);
        TEST( approx(v3, v4), "Direct quaternion transformation v2->v4 failed!" );
    }

    { // Test LHS and RHS orthonormal systems
        OrthoSpace o0{m0};
        OrthoSpace o1{m1};
        Mat3x3 m2 = diag(Vec3(-1.0f, 1.0f, 1.0f)) * rotation( 0.1f, -0.5f, -0.2f );
        Mat3x3 m3 = diag(Vec3(1.0f, -1.0f, 1.0f)) * rotation( 0.0f, -0.2f, 0.0f );
        Mat3x3 m4 = diag(Vec3(1.0f, 1.0f, -1.0f)) * rotation( -0.5f, -0.1f, -0.2f );
        OrthoSpace o2{m2};
        OrthoSpace o3{m3};
        OrthoSpace o4{m4};
        OrthoSpace o5{identity3x3()};
        OrthoSpace o6{diag(Vec3(1.0f, 1.0f, -1.0f))};
        TEST(approx(Mat3x3(o0), m0), "RHS rotation m0 -> OrthoSpace -> Mat3x3 wrong.");
        TEST(approx(Mat3x3(o1), m1), "RHS rotation m1 -> OrthoSpace -> Mat3x3 wrong.");
        TEST(approx(Mat3x3(o2), m2), "LHS rotation m5 -> OrthoSpace -> Mat3x3 wrong.");
        TEST(approx(Mat3x3(o3), m3), "LHS rotation m3 -> OrthoSpace -> Mat3x3 wrong.");
        TEST(approx(Mat3x3(o4), m4), "LHS rotation m4 -> OrthoSpace -> Mat3x3 wrong.");
        TEST(o0.isRighthanded(), "o0.isRighthanded wrong.");
        TEST(o1.isRighthanded(), "o1.isRighthanded wrong.");
        TEST(o2.isLefthanded(), "o2.isLefthanded wrong.");
        TEST(o3.isLefthanded(), "o3.isLefthanded wrong.");
        TEST(o4.isLefthanded(), "o4.isLefthanded wrong.");
        TEST(o5.isRighthanded(), "o5.isRighthanded wrong.");
        TEST(o6.isLefthanded(), "o6.isLefthanded wrong.");
    }

    return result;
}
