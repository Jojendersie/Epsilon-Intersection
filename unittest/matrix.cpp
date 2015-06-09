#include "ei/vector.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace ei;
using namespace std;

template<typename T>
static void doNotOptimizeAway(T& _var)
{
    static volatile T x = _var;
}

bool test_matrix()
{
    bool result = true;

    // ********************************************************************* //
    // Test type size
    TEST( sizeof(Matrix<byte, 1, 2>) == 2,
        "sizeof(Matrix<byte, 1, 2>) is " << sizeof(Matrix<byte, 1, 2>) << " instead of 2."
    );
    TEST( sizeof(Matrix<double, 2, 1>) == 16,
        "sizeof(Matrix<double, 2, 1>) is " << sizeof(Matrix<double, 2, 1>) << " instead of 16."
    );
    TEST( sizeof(Matrix<float, 1, 3>) == 12,
        "sizeof(Matrix<float, 1, 3>) is " << sizeof(Matrix<float, 1, 3>) << " instead of 12."
    );
    TEST( sizeof(Matrix<byte, 3, 1>) == 3,
        "sizeof(Matrix<byte, 3, 1>) is " << sizeof(Matrix<byte, 3, 1>) << " instead of 3."
    );
    TEST( sizeof(Matrix<float, 4, 1>) == 16,
        "sizeof(Matrix<float, 4, 1>) is " << sizeof(Matrix<float, 4, 1>) << " instead of 16."
    );
    TEST( sizeof(Matrix<int16, 1, 4>) == 8,
        "sizeof(Matrix<int16, 1, 4>) is " << sizeof(Matrix<int16, 1, 4>) << " instead of 8."
    );
    TEST( sizeof(Matrix<uint16, 2, 2>) == 8,
        "sizeof(Matrix<uint16, 2, 2>) is " << sizeof(Matrix<uint16, 2, 2>) << " instead of 8."
    );
    TEST( sizeof(Matrix<byte, 3, 3>) == 9,
        "sizeof(Matrix<byte, 3, 3>) is " << sizeof(Matrix<byte, 3, 3>) << " instead of 9."
    );
    TEST( sizeof(Matrix<float, 4, 4>) == 64,
        "sizeof(Matrix<float, 4, 4>) is " << sizeof(Matrix<float, 4, 4>) << " instead of 64."
    );
    TEST( sizeof(Matrix<float, 10, 10>) == 400,
        "sizeof(Matrix<float, 10, 10>) is " << sizeof(Matrix<float, 10, 10>) << " instead of 400."
    );

    // ********************************************************************* //
    // When included the following line must cause a static assertion.
    // Matrix<int, 1, 1> errorMatrix;

    // ********************************************************************* //
    // Test type conversion constructors
    {
        Vec2 v0(0);     // vector with zeros
//      Vec2 v1 = 1;    // This line should not compile (no accidentally conversion)
        Vec2 v2(-1);    // Test a scalar of a different type (should yield the standard warning)
        Vec2 v3 = static_cast<Vec2>(-1); // Test a scalar of a different type without warning
        IVec2 v4(v3);   // Convert to other vector
//      v4 = v0;        // This line should not compile (no accidentally conversion)
        IVec3 v5(-1);   // Used for the truncation tests
        Vec2 v6(v5);    // Truncate and type convert
        TEST( all(v6 == v2) , "Truncation operator wrong!" );
        Mat3x3 m0(0, 1, 2, 3, 4, 5, 6, 7, 8); // Used for more truncation tests
        Mat2x2 m1(m0);  // Truncate in both dimensions
        Vec3 v7(m0);    // Get first column
        Vec3 v8(m0, 0, 1);// Get second column
        TEST( all(m1 == Mat2x2(0, 1, 3, 4)), "Truncation in both dimensions failed!" );
        TEST( all(v7 == Vec3(0, 3, 6)), "Truncation to vector failed!" );
        TEST( all(v8 == Vec3(1, 4, 7)), "Truncation to vector failed!" );
    }

    // ********************************************************************* //
    // Test vector composition constructors for column vectors
    {
        Vec2 v0(0.0f, 0.5f);
        IVec2 v1(1, 2);
        Vec3 v2(0.0f, 0.5f, 2.0f);
        Vec3 v3(0.5f, 1.0f, 2.0f);
        IVec4 v4(1, 2, 0, 0);
        Vec4 v5(0.0f, 0.5f, 1.0f, 2.0f);
        Vec4 v6(1.0f, 0.0f, 0.5f, 2.0f);
        Vec4 v7(1.0f, 2.0f, 0.0f, 0.5f);
        Vec4 v8(0.5f, 1.0f, 2.0f, 7.5f);
        TEST( all(Vec3(v0, 2.0) == v2), "vec2+scalar to vec3 failed!" );
        TEST( all(Vec3(0.5f, v1) == v3), "scalar+vec2 to vec3 failed!" );
        TEST( all(IVec4(v1, v0) == v4), "vec2+vec2 to vec4 failed!" );
        TEST( all(Vec4(v0, 1.0f, 2.0f) == v5), "vec2+2*scalar to vec4 failed!" );
        TEST( all(Vec4(1.0f, v0, 2.0f) == v6), "scalar+vec2+scalar to vec4 failed!" );
        TEST( all(Vec4(1.0f, 2.0f, v0) == v7), "2*scalar+vec2 to vec4 failed!" );
        TEST( all(Vec4(1.0f, v2) == v6), "scalar+vec3 to vec4 failed!" );
        TEST( all(Vec4(v3, 7.5f) == v8), "vec2+scalar to vec4 failed!" );
    }

    // Test vector composition constructors for row vectors
    {
        RVec2 v0(0.0f, 0.5f);
        IRVec2 v1(1, 2);
        RVec3 v2(0.0f, 0.5f, 2.0f);
        RVec3 v3(0.5f, 1.0f, 2.0f);
        IRVec4 v4(1, 2, 0, 0);
        RVec4 v5(0.0f, 0.5f, 1.0f, 2.0f);
        RVec4 v6(1.0f, 0.0f, 0.5f, 2.0f);
        RVec4 v7(1.0f, 2.0f, 0.0f, 0.5f);
        RVec4 v8(0.5f, 1.0f, 2.0f, 7.5f);
        TEST( all(RVec3(v0, 2.0) == v2), "rvec2+scalar to rvec3 failed!" );
        TEST( all(RVec3(0.5f, v1) == v3), "scalar+rvec2 to rvec3 failed!" );
        TEST( all(IRVec4(v1, v0) == v4), "rvec2+rvec2 to rvec4 failed!" );
        TEST( all(RVec4(v0, 1.0f, 2.0f) == v5), "rvec2+2*scalar to rvec4 failed!" );
        TEST( all(RVec4(1.0f, v0, 2.0f) == v6), "scalar+rvec2+scalar to rvec4 failed!" );
        TEST( all(RVec4(1.0f, 2.0f, v0) == v7), "2*scalar+rvec2 to rvec4 failed!" );
        TEST( all(RVec4(1.0f, v2) == v6), "scalar+rvec3 to rvec4 failed!" );
        TEST( all(RVec4(v3, 7.5f) == v8), "rvec2+scalar to rvec4 failed!" );
    }

    // ********************************************************************* //
    // Test of boolean operators
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(2, 2, 0, 4);
        Matrix<int, 2, 2> m2(2, 1, 0, 4);
        TEST( sum(2 == m2) == 1 , "Scalar operator == or sum() wrong!" );
        TEST( sum(m2 == 3) == 0 , "Scalar operator == or sum() wrong!" );
        TEST( sum(2 != m2) == 3 , "Scalar operator != or sum() wrong!" );
        TEST( sum(m2 != 3) == 4 , "Scalar operator != or sum() wrong!" );
        TEST( any(m0 < 2) , "Scalar operator < or any() wrong!" );
        TEST( !any(m0 < 0) , "Scalar operator < or any() wrong!" );
        TEST( sum(2 < m1) == 1 , "Scalar operator < or sum() wrong!" );
        TEST( !none(m0 <= 2) , "Scalar operator <= or none() wrong!" );
        TEST( none(m0 <= 0) , "Scalar operator <= or none() wrong!" );
        TEST( sum(2 <= m1) == 3 , "Scalar operator <= or sum() wrong!" );
        TEST( all(m0 >= 1) , "Scalar operator >= or all() wrong!" );
        TEST( !all(m1 >= 2) , "Scalar operator >= or all() wrong!" );
        TEST( sum(4 >= m1) == 4 , "Scalar operator >= or sum() wrong!" );
        TEST( any(m0 > 2) , "Scalar operator > or all() wrong!" );
        TEST( none(m0 > 4) , "Scalar operator > or all() wrong!" );
        TEST( sum(4 > m1) == 3 , "Scalar operator > or sum() wrong!" );
        TEST( sum(m0 < m1) == 1 , "Matrix operator < or sum() wrong!" );
        TEST( sum(m0 <= m1) == 3 , "Matrix operator <= or sum() wrong!" );
        TEST( sum(m0 >= m1) == 3 , "Matrix operator >= or sum() wrong!" );
        TEST( sum(m0 > m1) == 1 , "Matrix operator > or sum() wrong!" );
        TEST( sum(m0 == m1) == 2 , "Matrix operator == or sum() wrong!" );
        TEST( sum(m0 != m2) == 3 , "Matrix operator != or sum() wrong!" );
    }


    // ********************************************************************* //
    // Test of access operators and element constructors (2-4 elements)
    {
        Matrix<int, 2, 1> m0(1, 2);
        Matrix<int, 1, 3> m1(1, 2, 3);
        Matrix<int, 2, 2> m2(1, 2, 3, 4);
        m1(0,2) = 5;
        m2[2] = 5;
        TEST( m0[0] == 1 && m0[1] == 2, "Operator or constructor wrong!" );
        TEST( m1(0,0) == 1 && m1(0,1) == 2 && m1[2] == 5, "Operator or constructor wrong!" );
        TEST( m2(0,1) == 2 && m2(1,0) == 5 && m2(1,1) == 4, "Operator or constructor wrong!" );
        TEST( (const_cast<const Matrix<int, 2, 2>&>(m2))[1] == 2, "Read-only operator wrong!" );
        TEST( (const_cast<const Matrix<int, 2, 2>&>(m2))(1,0) == 5, "Read-only operator wrong!" );
    }

    // ********************************************************************* //
    // Test + and -
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(5, 6, 7, 8);
        Matrix<int, 2, 2> m2(6, 8, 10, 12);
        Matrix<float, 2, 2> m3(3.0f, 1.0f, 4.0f, 1.0f);
        Matrix<float, 2, 2> m4(4.0f, 3.0f, 7.0f, 5.0f);
        TEST( all(m0 + m1 == m2), "Component wise addition failed!" );
        TEST( all(m0 + m3 == m4), "Component wise addition failed!" );
        TEST( all(m2 - m1 == m0), "Component wise subtraction failed!" );
        TEST( all(m0 - m4 == -m3), "Component wise subtraction or unary minus failed!" );
    }

    // ********************************************************************* //
    // Test |, &, ^ and ~
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(5, 6, 7, 8);
        Matrix<int, 2, 2> m2(5, 6, 7, 12);
        Matrix<int, 2, 2> m3(1, 2, 3, 0);
        Matrix<int, 2, 2> m4(4, 4, 4, 12);
        Matrix<int, 2, 2> m5(0xfffffffe, 0xfffffffd, 0xfffffffc, 0xfffffffb);
        TEST( all((m0 | m1) == m2), "Component wise | failed!" );
        TEST( all((m0 & m1) == m3), "Component wise & failed!" );
        TEST( all((m0 ^ m1) == m4), "Component wise ^ failed!" );
        TEST( all((~m0) == m5), "Component wise ~ failed!" );
        // Generates an error message (expected)
        //Matrix<float, 2, 2> m3(3.0f, 1.0f, 4.0f, 1.0f);
        //Matrix<float, 2, 2> m4(4.0f, 3.0f, 7.0f, 5.0f);
        //Matrix<float, 2, 2> m = m3 | m4;
    }

    // ********************************************************************* //
    // Test *, / and element constructors (3,9 elements)
    {
        Matrix<int, 3, 1> v0(1, 2, 3);
        Matrix<int, 3, 1> v4(2, 1, 4);
        Matrix<int, 3, 1> v5(2, 2, 12);
        Matrix<float, 1, 3> v1(1.0f, 0.0f, -1.0f);
        Matrix<float, 1, 3> v2(0.5f, 2.0f, 0.0f);
        Matrix<float, 1, 3> v3(0.5f, 0.0f, 0.0f);
        Matrix<float, 1, 3> v6(2.0f, 3.0f, 4.0f);
        Matrix<float, 1, 3> v7(0.5f, 0.0f, -0.25f);
        Matrix<float, 3, 3> m0(1.0f, 0.0f, -1.0f, 2.0f, 0.0f, -2.0f, 3.0f, 0.0f, -3.0f);
        Matrix<float, 3, 3> m1(1.0f, 4.0f, 0.0f, 0.5f, 2.0f, 0.0f, 2.0f, 8.0f, 0.0f);
        Matrix<float, 3, 3> m2(-1.0f, -4.0f, 0.0f, -2.0f, -8.0f, 0.0f, -3.0f, -12.0f, 0.0f);
        TEST( all(v1 * v2 == v3), "Component wise multiplication failed!" );
        TEST( all(v0 * v4 == v5), "Component wise multiplication failed!" );
        TEST( all(v1 / v6 == v7), "Component wise division failed!" );
        TEST( all(v0 == v5 / v4), "Component wise division failed!" );
        TEST( v1 * v0 == -2.0f, "Row times column vector should be a scalar!" );
        TEST( all(v0 * v1 == m0), "Column times row vector should be a matrix!" );
        TEST( all(v4 * v2 == m1), "Column times row vector should be a matrix!" );
        TEST( all(m0 * m1 == m2), "Matrix multiplication wrong!" );
    }

    // ********************************************************************* //
    // Test * and element constructors (6,8,9,12,16 elements)
    {
        Matrix<int, 2, 3> m0(1, 2, 3,
                             2, 3, 4);
        Matrix<int, 3, 2> m1(1, 2,
                             2, 3,
                             3, 4);
        Matrix<int, 3, 3> m2(5, 8, 11,
                             8, 13, 18,
                             11, 18, 25);
        TEST( all(m1 * m0 == m2), "Constructor or matrix multiplication wrong!" );
        Matrix<int, 2, 4> m3(1, 2, 3, 4,
                             2, 3, 4, 5);
        Matrix<int, 4, 2> m4(1, 2,
                             2, 3,
                             3, 4,
                             4, 5);
        Matrix<int, 4, 4> m5(5, 8, 11, 14,
                             8, 13, 18, 23,
                             11, 18, 25, 32,
                             14, 23, 32, 41);
        TEST( all(m4 * m3 == m5), "Constructor or matrix multiplication wrong!" );
        Matrix<int, 3, 4> m6(1, 2, 3, 4,
                             2, 3, 4, 5,
                             0, 0, 0, 0);
        Matrix<int, 4, 3> m7(1, 2, 3,
                             2, 3, 4,
                             3, 4, 5,
                             4, 5, 6);
        TEST( all(m7 * m6 == m5), "Constructor or matrix multiplication wrong!" );
    }

    // ********************************************************************* //
    // Test scalar operators +, *, -, /
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(3, 4, 5, 6);
        Matrix<float, 1, 2> v0(1.0f, 1.0f);
        Matrix<float, 2, 1> v1(0.0f, 1.0f);
        Matrix<float, 2, 2> m2(2.0f, 2.0f, 4.0f, 4.0f);
        TEST( all(m0 + 2 == m1), "Adding a scalar failed!" );
        TEST( all(2 + m0 * 1 == m1), "Multiplying or adding a scalar failed!" );
        TEST( 2.0f * v0 * 2.0f * v1 + 0.5f == 4.5f, "Multiplying a scalar failed!" );
        TEST( all((v1 + 1.0f) * 2.0f * v0 == m2), "Mixed scalar vector operation failed!" );
        TEST( all(m0 == m1 - 2), "Scalar subtraction failed!" );
        TEST( all(2 + m0 / 1 == m1), "Scalar addition or division failed!" );
        TEST( (2.0f / v0 / 2.0f) * v1 - 0.5f == 0.5f, "Mixed scalar vector operation (/,*,-) failed!" );
        TEST( all((v1 - 0.5f) * 2.0f * v0 + 3 == m2), "Mixed scalar vector operation (-, *) failed!" );
        TEST( all((4.0f - m2) / 2.0f == (1.0f - v1) * v0), "Mixed scalar vector operation (/,*,-) failed!" );
    }

    // ********************************************************************* //
    // Test scalar |, &, ^
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(5, 6, 7, 4);
        Matrix<int, 2, 2> m2(0, 0, 0, 4);
        Matrix<int, 2, 2> m3(5, 6, 7, 0);
        TEST( all((m0 | 4) == m1), "Component wise | failed!" );
        TEST( all((4 | m0) == m1), "Component wise | failed!" );
        TEST( all((m0 & 4) == m2), "Component wise & failed!" );
        TEST( all((4 & m0) == m2), "Component wise & failed!" );
        TEST( all((m0 ^ 4) == m3), "Component wise ^ failed!" );
        TEST( all((4 ^ m0) == m3), "Component wise ^ failed!" );
        // Generates an error message (expected)
        //Matrix<float, 2, 2> m4(3.0f, 1.0f, 4.0f, 1.0f);
        //Matrix<float, 2, 2> m = m4 | 4.0f;
    }

    // ********************************************************************* //
    // Test assignment operators +=, *=, -=, /=
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(3, 4, 5, 6);
        Matrix<int, 2, 2> m2(4, 6, 8, 10);
        Matrix<float, 2, 2> m3(2.0f, 2.0f, 4.0f, 4.0f);
        Matrix<float, 2, 2> m4(5.0f, 6.0f, 9.0f, 10.0f);
        Vec2 v0(1.0f, 1.0f);
        Vec2 v1(0.0f, 1.0f);
        Vec2 v2(0.0f, 1.0f);
        Matrix<float, 1, 2> vr0(1.0f, 1.0f);
        Matrix<float, 1, 2> vr1(0.0f, 1.0f);
        Matrix<float, 1, 2> vr2(0.0f, 1.0f);
        TEST( all((m0 += m1) == m2), "Component wise += failed!" );
        m2 -= m0;
        TEST( all(m2 == Matrix<int,2,2>(0)), "Component wise -= failed!" );
        m3 += m1;
        TEST( all(m3 == m4), "Component wise += with different types failed!" );
        v0 *= v1;
        TEST( all(v0 == v1), "Component wise *= for vectors failed!" );
        (v1 += v0) /= v0 + 1.0f;
        TEST( all(v1 == v2), "Component wise /= for vectors failed!" );
        vr0 *= vr1;
        TEST( all(vr0 == vr1), "Component wise *= for row vectors failed!" );
        (vr1 += vr0) /= vr0 + 1.0f;
        TEST( all(vr1 == vr2), "Component wise /= for row vectors failed!" );

        TEST( all((vr1 += 1.0f) == vr2 + 1), "Matrix-Scalar += failed!");
        TEST( all((m3 -= 0.5f) == m4 - 0.5f), "Matrix-Scalar -= failed!");
        TEST( all((v1 *= 2.0f) == v2 / 0.5f), "Matrix-Scalar *= failed!");
        TEST( all((v1 /= 2.0f) == v2), "Matrix-Scalar /= failed!");
    }

    // ********************************************************************* //
    // Test len, lensq, dot
    {
        Vec3 v0(1.0f, 2.0f, 0.0f);
        Vec3 v1(0.0f, 0.0f, 0.5f);
        IVec3 v2(-1, 0, -4);
        Mat3x3 m0(1.0f, 1.0f, 2.0f,
                  0.5f, 0.5f, 0.5f,
                  0.5f, 1.0f, 1.0f);
        TEST( dot(v0, v1) == 0.0f, "Dot product wrong!" );
        TEST( dot(v0, v0) == 5.0f, "Dot product wrong!" );
        TEST( dot(v1, v2) == -2.0f, "Dot product of a mixed vector types wrong!" );
        TEST( lensq(v2) == 17, "Squared length of an integer vector wrong!" );
        TEST( lensq(v1) == 0.25f, "Squared length wrong!" );
        TEST( lensq(m0) == 9.0f, "Squared length of a matrix wrong!" );
        TEST( len(v2) == 4.1231056256176606, "Length of an integer vector wrong!" );
        TEST( len(v1) == 0.5f, "Length wrong!" );
        TEST( len(m0) == 3.0f, "Length of a matrix wrong!" );

        /*// Test script used to optimize dot()
        Vec3 v3(10.0f, -0.2f, 0.5f);
        uint64 a = ticks();
        for(uint i=0; i<10000000; ++i) {
            volatile float f = dot(v0, v3);
        }
        uint64 b = ticks();
        cerr << (int)(b-a) << endl;*/
    }

    // ********************************************************************* //
    // Test max, min, sum, avg, abs, sgn, sign, prod
    {
        Vec4 v0(-1.0f, -0.5f, 0.25f, 0.5f);
        Vec4 v1(0.5f, 2.5f, 1.0f, 0.0f);
        Vec4 v2(-1.0f, -0.5f, 0.25f, 0.0f);
        Vec4 v3(0.5f, 2.5f, 1.0f, 0.5f);
        Vec4 v4(1.0f, 0.5f, 0.25f, 0.5f);
        Vec4 v5(-1.0f, -0.5f, 0.25f, 0.0f);
        Vec4 v6(0.0f, 0.0f, 0.25f, 0.5f);
        Vec4 v7(0.0f, 0.0f, 1.0f, 1.0f);
        TEST( all(min(v0, v1) == v2), "Component wise minimum wrong!" );
        TEST( all(max(v0, v1) == v3), "Component wise maximum wrong!" );
        TEST( min(v0) == -1.0f, "Minimum of a vector wrong!" );
        TEST( min(v1) == 0.0f, "Minimum of a vector wrong!" );
        TEST( max(v2) == 0.25f, "Maximum of a vector wrong!" );
        TEST( max(v3) == 2.5f, "Maximum of a vector wrong!" );
        TEST( all(clamp(v0, v2, v1) == v5), "Clamp of a vector in box boundaries wrong!" );
        TEST( all(clamp(v0, 0.0f, 3.0f) == v6), "Clamp of a vector in interval wrong!" );
        TEST( all(saturate(v0 * 4.0f) == v7), "Saturate of a vector wrong!" );
        TEST( sum(v1) == 4.0f, "Sum of a vector wrong!" );
        TEST( prod(v0) == 0.0625f, "Product of a vector wrong!" );
        TEST( avg(v3) == 1.125f, "Average of a vector wrong!" );
        TEST( all(abs(v0) == v4), "Absolute value of a vector wrong!" );
        TEST( all(sign(v1) == Vec4(1.0f, 1.0f, 1.0f, 0.0f)), "Sign values of a vector wrong!" );
        TEST( all(sgn(v2) == Vec4(-1.0f, -1.0f, 1.0f, 1.0f)), "Sgn values of a vector wrong!" );
    }

    // ********************************************************************* //
    // Test round, ceil, floor, sqrt, pow
    {
        Vec4 v0(1.5f, 2.0f, -0.5f, 1.414f);
        IVec4 v1(2, 2, 0, 1);    // round(v0)
        IVec4 v2(2, 2, 0, 2);    // ceil(v0)
        IVec4 v3(1, 2, -1, 1);   // floor(v0)
        Vec4 v6(5.5f, -1.5f, -0.5f, 2.5f);
        IVec4 v7(6, -2, 0, 2);   // round(v6)
        TEST( all(round(v0) == v1), "Rounding v0 is wrong!" );
        TEST( all(ceil(v0) == v2), "Rounding v0 up is wrong!" );
        TEST( all(floor(v0) == v3), "Rounding v0 down is wrong!" );
        TEST( all(round(v6) == v7), "Rounding v0 is wrong!" );

        Vec3 v4(1.0f, 4.0f, 2.0f);
        Vec3 v5(1.0f, 2.0f, PHYTAGORAS);
        TEST( all(sqrt(v4) == v5), "sqrt of v4 is wrong!" );
        TEST( all(pow(v4,0.5f) == v5), "v4 to the power of 0.5 is wrong!" );
        TEST( approx(pow(v5,2.0f), v4), "v5 to the power of 2 is wrong!" );
    }

    // ********************************************************************* //
    // Test lerp, bilerp, slerp
    {
        Vec3 v0(0.0f);
        Vec3 v1(1.0f);
        Vec3 v2(2.0f);
        Vec3 v3(-1.0f);
        Vec3 v4(0.0f, 1.0f, 0.0f);
        Vec3 v5(1.0f, 0.0f, 0.0f);
        Vec3 v6(0.382683432f, 0.923879533f, 0.0f);
        Vec3 v7(0.25f, 0.75f, 0.0f);
        TEST( all(lerp(v0, v2, 0.5f) == v1), "Linear interpolation failed." );
        TEST( all(bilerp(v0, v1, v1, v2, 0.5f, 0.5f) == v1), "Bilinear interpolation failed." );
        TEST( lerp(1, 4, 0.5f) == 2.5f, "Linear interpolation failed." );
        TEST( all(lerp(v4, v5, 0.25f) == v7), "Linear interpolation by 0.25 failed." );
        TEST( bilerp(1, 3, 5, 7, 0.5f, 0.5f) == 4.0f, "Bilinear interpolation failed." );
        TEST( approx(slerp(v4, v5, 0.25f), v6), "Spherical linear interpolation failed." );
        TEST( approx(dot(slerp(v1, v3, 0.5f), v1), 0.0f), "Spherical linear interpolation for degenerated case failed." );
    }

    // ********************************************************************* //
    // Test crossproduct, unary minus
    {
        Vec3 v0(1.0f, 0.0f, 0.0f);
        Vec3 v1(0.0f, 1.0f, 0.0f);
        Vec3 v2(0.0f, 0.0f, 1.0f);
        Matrix<float,1,3> v5(1.0f, 0.0f, 0.0f);
        Matrix<float,1,3> v6(0.0f, 1.0f, 0.0f);
        Matrix<float,1,3> v7(0.0f, 0.0f, 1.0f);
        TEST( all(cross(v0, v1) == v2), "3D crossproduct wrong." );
        TEST( all(cross(v1, v0) == -v2), "3D crossproduct wrong." );
        TEST( all(cross(v5, v6) == v7), "3D crossproduct wrong." );
        TEST( all(cross(v6, v5) == -v7), "3D crossproduct wrong." );
        Vec2 v3(1.0f, 0.0f);
        Vec2 v4(0.0f, 1.0f);
        TEST( cross(v3, v4) == 1.0f, "2D crossproduct wrong." );
        TEST( cross(v4, v3) == -1.0f, "2D crossproduct wrong." );
    }

    // ********************************************************************* //
    // Identity and diagonal
    {
        Mat2x2 m0(1.0f, 0.0f, 0.0f, 1.0f);
        Mat3x3 m1(1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f);
        Mat4x4 m2(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
        Matrix<float,5,5> m3(0.0f); m3[0] = 1.0f; m3[6] = 1.0f; m3[12] = 1.0f; m3[18] = 1.0f; m3[24] = 1.0f;
        Mat3x3 m4(1.0f, 0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, -1.0f);
        Mat3x3 m5(diag(Vec3(1.0f, 2.0f, -1.0f)));
        TEST( all(m0 == identity2x2()), "2x2 identity matrix wrong." );
        TEST( all(m1 == identity3x3()), "3x3 identity matrix wrong." );
        TEST( all(m2 == identity4x4()), "4x4 identity matrix wrong." );
        TEST( all(m3 == ei::identity<float,5>()), "5x5 identity matrix wrong." );
        TEST( all(m4 == m5), "3x3 diagonal matrix wrong." );

        /*// Test script used to optimize identity()
        Mat3x3 ret;
        uint64 a = ticks();
        for(uint i=0; i<1000000; ++i) {
            ret = ei::identity<float,3>();
            doNotOptimizeAway(ret);
        }
        uint64 b = ticks();
        cerr << (int)(b-a) << endl;*/
    }

    // ********************************************************************* //
    // Test transpose
    {
        Mat2x2 m0(1.0f, -1.0f, 0.0f, 2.0f);
        Mat2x2 m1(1.0f, 0.0f, -1.0f, 2.0f);
        Matrix<float,1,3> v0(1.0f, 0.0f, 0.0f);
        Vec3 v1(1.0f, 0.0f, 0.0f);

        TEST( all(m0 == transpose(m1)), "m1 transposed is unequal m0!" );
        TEST( all(v0 == transpose(v1)), "v1 transposed is unequal v0!" );
    }

    // ********************************************************************* //
    // Test axis matrix
    {
        Mat3x3 m0 = axis(normalize(Vec3(1.0f, 1.0f, 0.0f)), normalize(Vec3(-1.0f, 1.0f, 0.0f)), Vec3(0.0f, 0.0f, 1.0f));
        Vec3 v0(0.0f, 0.0f, 5.0f);
        Vec3 v1(2.0f, 0.0f, 0.0f);
        Vec3 v2(1.414213562f, -1.414213562f, 0.0f);
        TEST( all(m0 * v0 == v0), "v0 should not be transformed by m0!" );
        TEST( all(m0 * v1 == v2), "v1 should be transformed into v2 by m0!" );
    }

    // ********************************************************************* //
    // Test orthonormal basis
    {
        Mat2x2 m0 = basis(Vec2(1.0f, 0.0f));
        TEST( m0(1,0) == 0.0f && m0(1,1) == 1.0f, "2D orthonormal basis of (1,0) wrong!" );
        m0 = basis(normalize(Vec2(1.0f, -2.0f)));
        TEST( dot(m0(0), m0(1)) == 0.0f, "2D orthonormal basis of (1,-2) wrong!" );
        // In 3D the spat product must be 1
        Mat3x3 m1 = basis(Vec3(0.0f, 1.0f, 0.0f));
        TEST( approx(dot(m1(0), m1(1)), 0.0f) && approx(dot(m1(0), m1(0)), 1.0f)
            && approx(dot(m1(0), m1(2)), 0.0f) && approx(dot(m1(1), m1(1)), 1.0f)
            && approx(dot(m1(2), m1(1)), 0.0f) && approx(dot(m1(2), m1(2)), 1.0f),
            "The basis of (0,1,0) is not orthonormal!" );
        m1 = basis(Vec3(0.0f, 0.0f, 1.0f));
        TEST( approx(dot(m1(0), m1(1)), 0.0f) && approx(dot(m1(0), m1(0)), 1.0f)
            && approx(dot(m1(0), m1(2)), 0.0f) && approx(dot(m1(1), m1(1)), 1.0f)
            && approx(dot(m1(2), m1(1)), 0.0f) && approx(dot(m1(2), m1(2)), 1.0f),
            "The basis of (0,0,1) is not orthonormal!" );
        m1 = basis(Vec3(-1.0f, 0.0f, 0.0f));
        TEST( approx(dot(m1(0), m1(1)), 0.0f) && approx(dot(m1(0), m1(0)), 1.0f)
            && approx(dot(m1(0), m1(2)), 0.0f) && approx(dot(m1(1), m1(1)), 1.0f)
            && approx(dot(m1(2), m1(1)), 0.0f) && approx(dot(m1(2), m1(2)), 1.0f),
            "The basis of (-1,0,0) is not orthonormal!" );
        m1 = basis(normalize(Vec3(0.5f, -0.7f, 0.1f)));
        TEST( approx(dot(m1(0), m1(1)), 0.0f) && approx(dot(m1(0), m1(0)), 1.0f)
            && approx(dot(m1(0), m1(2)), 0.0f) && approx(dot(m1(1), m1(1)), 1.0f)
            && approx(dot(m1(2), m1(1)), 0.0f) && approx(dot(m1(2), m1(2)), 1.0f),
            "The basis of normalize(0.5,-0.7,0.1) is not orthonormal!" );
    }
    
    // ********************************************************************* //
    // Test transformations
    {
        Vec3 vx(1.0f, 0.0f, 0.0f);
        Vec3 vy(0.0f, 1.0f, 0.0f);
        Vec3 vz(0.0f, 0.0f, 1.0f);
        TEST( approx(rotation(PI/2, 0.0f, 0.0f) * vy, vz), "Rotation matrix from Euler angles invalid!" );
        TEST( approx(rotation(-PI/2, 0.0f, 0.0f) * vz, vy), "Rotation matrix from Euler angles invalid!" );
        TEST( approx(rotation(0.0f, -PI/2, 0.0f) * vx, vz), "Rotation matrix from Euler angles invalid!" );
        TEST( approx(rotation(0.0f, PI/2, 0.0f) * vz, vx), "Rotation matrix from Euler angles invalid!" );
        TEST( approx(rotation(0.0f, 0.0f, PI/2) * vx, vy), "Rotation matrix from Euler angles invalid!" );
        TEST( approx(rotation(0.0f, 0.0f, -PI/2) * vy, vx), "Rotation matrix from Euler angles invalid!" );

        TEST( approx(rotation(vx, vx), identity3x3()), "Rotation matrix from vector to vector invalid!" );
        TEST( approx(rotation(vy, vz), rotation(PI/2, 0.0f, 0.0f)), "Rotation matrix from vector to vector invalid!" );

        TEST( approx(housholder(Vec3(-1.0f, 1.0f, 0.0f)) * vy, vx), "Housholder matrix invalid!" );
        TEST( approx(housholder(Vec3(1.0f, -1.0f, 0.0f)) * vy, vx), "Housholder matrix invalid!" );
    }

    // ********************************************************************* //
    // Test spherical coordinates
    {
        Vec2 v0(1.0f, 0.0f);        Vec2 s0(1.0f, 0.0f);
        Vec2 v1(1.0f, 1.0f);        Vec2 s1(PHYTAGORAS, PI/4.0f);
        Vec2 v2(-3.0f, -4.0f);      Vec2 s2(5.0f, 4.06888771f);
        Vec3 v3(1.0f, 0.0f, 0.0f);  Vec3 s3(1.0f, 0.0f, 0.0f);
        Vec3 v4(-1.0f, 2.0f, 3.0f); Vec3 s4(3.741657387f, 1.84134603f, 0.982793723f);
        Vec3 v5(0.0f, 1.0f, -1.0f); Vec3 s5(PHYTAGORAS, PI/2.0f, 5.49778748f);
        Vec<float,5> v6(0.0f, -2.0f, 2.5f, -1.0f, 0.0f);

        TEST( all(sphericalCoords(v0) == s0), "Spherical coordinates of v0 wrong!" );
        TEST( all(sphericalCoords(v1) == s1), "Spherical coordinates of v1 wrong!" );
        TEST( all(sphericalCoords(v2) == s2), "Spherical coordinates of v2 wrong!" );
        TEST( all(sphericalCoords(v3) == s3), "Spherical coordinates of v3 wrong!" );
        TEST( all(sphericalCoords(v4) == s4), "Spherical coordinates of v4 wrong!" );
        TEST( all(sphericalCoords(v5) == s5), "Spherical coordinates of v5 wrong!" );

        TEST( all(v0 == cartesianCoords(s0)), "Cartesian coordinates of s0 wrong!" );
        TEST( approx(v1, cartesianCoords(s1)), "Cartesian coordinates of s1 wrong!" );
        TEST( approx(v2, cartesianCoords(s2)), "Cartesian coordinates of s2 wrong!" );
        TEST( approx(v3, cartesianCoords(s3)), "Cartesian coordinates of s3 wrong!" );
        TEST( approx(v4, cartesianCoords(s4)), "Cartesian coordinates of s4 wrong!" );
        TEST( approx(v5, cartesianCoords(s5)), "Cartesian coordinates of s5 wrong!" );
        TEST( approx(v6, cartesianCoords(sphericalCoords(v6))), "5D vector coordinate transformation to spherical and back failed!" );
    }

    // ********************************************************************* //
    // Test LUp decomposition and inverse
    {
        Mat3x3 A0(3.0f, 2.0f, -1.0f,
                 2.0f, -2.0f, 4.0f,
                 -1.0f, 0.5f, -1.0f);
        Mat3x3 A1(0.0f, 2.0f, -4.0f,
                  4.0f, 1.0f, 0.0f,
                  8.0f, 5.0f, -6.0f);
        Mat4x4 A2(1.0f, 2.0f, 3.0f, 4.0f,
                  1.0f, 2.0f, 4.0f, 8.0f,
                  1.0f, 3.0f, 5.0f, 7.0f,
                  27.0f, 9.0f, 3.0f, 1.0f);
        Vec3 b0(1.0f, -2.0f, 0.0f);
        Vec3 x0(1.0f, -2.0f, -2.0f);
        Mat3x3 LU, X;
        UVec3 p;
        Vec3 x;
        TEST( !decomposeLUp(A1, LU, p), "Matrix A1 is singular!");
        TEST( decomposeLUp(A0, LU, p), "Matrix A0 is decomposable!" );
        x = solveLUp(LU, p, b0);
        TEST( approx(x, x0, 1e-5f), "Solution of equation system A0 x=b0 wrong!");

        // Test inverse
        X = solveLUp(LU, p, identity3x3());
        TEST( all(invert(A0) == X), "Manual inversion and invert function have different results!" );
        X = X * A0;
        TEST( approx(X, identity3x3()), "3x3 Matrix inverse bad!");
        TEST( approx(invert(A2) * A2, identity4x4(), 2e-5f), "4x4 Matrix inverse bad!");
    }

    // ********************************************************************* //
    // Test determinants
    {
        Mat2x2 m0(3.0f, 8.0f, 4.0f, 6.0f);
        Mat3x3 m1(-2.0f, 2.0f, -3.0f, -1.0f, 1.0f, 3.0f, 2.0f, 0.0f, -1.0f);
        Mat3x3 m2(6.0f, 1.0f, 1.0f, 4.0f, -2.0f, 5.0f, 2.0f, 8.0f, 7.0f);
        Mat3x3 m3(1.0f, 2.0f, 1.0f, 3.0f, 0.0f, 4.0f, 8.0f, 4.0f, 10.0f);
        Mat4x4 m4(1.0f, 2.0f, 2.0f, 1.0f, 0.0f, 2.0f, 0.0f, 0.0f, -1.0f, -4.0f, -2.0f, -1.0f, 1.0f, 4.0f, 4.0f, 1.0f);
        Mat4x4 m5(1.0f, 0.0f, 2.0f, -1.0f, 3.0f, 0.0f, 0.0f, 5.0f, 2.0f, 1.0f, 4.0f, -3.0f, 1.0f, 0.0f, 5.0f, 0.0f);
        Mat4x4 m6(3.0f, 0.0f, 0.0f, 5.0f, 1.0f, 0.0f, 2.0f, -1.0f, 2.0f, 1.0f, 4.0f, -3.0f, 1.0f, 0.0f, 5.0f, 0.0f);
        Mat4x4 m7(3.0f, 2.0f, 0.0f, 1.0f, 4.0f, 0.0f, 1.0f, 2.0f, 3.0f, 0.0f, 2.0f, 1.0f, 9.0f, 2.0f, 3.0f, 1.0f);
        TEST( determinant(m0) == -14.0f, "Determinant of m0 wrong!" );
        TEST( determinant(m1) == 18.0f, "Determinant of m1 wrong!" );
        TEST( determinant(m2) == -306.0f, "Determinant of m2 wrong!" );
        TEST( determinant(m3) == 0.0f, "Determinant of m3 wrong!" );
        TEST( determinant(m4) == 0.0f, "Determinant of m4 wrong!" );
        TEST( determinant(m5) == 30.0f, "Determinant of m5 wrong!" );
        TEST( determinant(m6) == -30.0f, "Determinant of m6 wrong!" );
        TEST( approx(determinant(m7), 24.0f, 4.0e-6f), "Determinant of m7 wrong!" );
    }

    return result;
}
