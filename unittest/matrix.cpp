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

// Try to instanciate all relevant vector types explicitly.
template class Matrix<float, 2, 1>; // Vec2;
template class Matrix<float, 3, 1>; // Vec3;
template class Matrix<float, 4, 1>; // Vec4;
template class Matrix<int32, 2, 1>; // IVec2;
template class Matrix<int32, 3, 1>; // IVec3;
template class Matrix<int32, 4, 1>; // IVec4;
template class Matrix<uint32, 2, 1>; // UVec2;
template class Matrix<uint32, 3, 1>; // UVec3;
template class Matrix<uint32, 4, 1>; // UVec4;
template class Matrix<float, 2, 2>; // Mat2x2;
template class Matrix<float, 3, 3>; // Mat3x3;
template class Matrix<float, 3, 4>; // Mat3x4;
template class Matrix<float, 4, 3>; // Mat4x3;
template class Matrix<float, 4, 4>; // Mat4x4;

bool test_matrix()
{
    bool result = true;

    // ********************************************************************* //
    // Test type size
    static_assert( sizeof(Matrix<byte, 1, 2>) == 2, "sizeof(Matrix<byte, 1, 2>) is invalid." );
    static_assert( sizeof(Matrix<double, 2, 1>) == 16, "sizeof(Matrix<double, 2, 1>) is invalid." );
    static_assert( sizeof(Matrix<float, 1, 3>) == 12, "sizeof(Matrix<float, 1, 3>) is invalid." );
    static_assert( sizeof(Matrix<byte, 3, 1>) == 3, "sizeof(Matrix<byte, 3, 1>) is invalid." );
    static_assert( sizeof(Matrix<float, 4, 1>) == 16, "sizeof(Matrix<float, 4, 1>) is invalid." );
    static_assert( sizeof(Matrix<int16, 1, 4>) == 8, "sizeof(Matrix<int16, 1, 4>) is invalid." );
    static_assert( sizeof(Matrix<uint16, 2, 2>) == 8, "sizeof(Matrix<uint16, 2, 2>) is invalid." );
    static_assert( sizeof(Matrix<byte, 3, 3>) == 9, "sizeof(Matrix<byte, 3, 3>) is invalid." );
    static_assert( sizeof(Matrix<float, 4, 4>) == 64, "sizeof(Matrix<float, 4, 4>) is invalid." );
    static_assert( sizeof(Matrix<float, 10, 10>) == 400, "sizeof(Matrix<float, 10, 10>) is invalid." );

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
        Vec3 v7(v2);
        TEST( v6 == v2 , "Truncation operator wrong!" );
        TEST( v7 == Vec3(-1.0f, -1.0f, 1.0f), "Column vector extension failed!" );

        Mat3x3 m0(0, 1, 2, 3, 4, 5, 6, 7, 8); // Used for more truncation tests
        Mat2x2 m1(m0);  // Truncate in both dimensions
        Mat4x4 m2(m0);  // Extend in both directions
        //Vec3 v7(m0);    // Get first column TODO: reanable this
        //Vec3 v8(m0, 0, 1);// Get second column
        TEST( m1 == Mat2x2(0, 1, 3, 4), "Truncation in both dimensions failed!" );
        TEST( m2 == Mat4x4(0, 1, 2, 0, 3, 4, 5, 0, 6, 7, 8, 0, 0, 0, 0, 1), "Extension in both dimensions failed!" );
//        TEST( all(v7 == Vec3(0, 3, 6)), "Truncation to vector failed!" );
 //       TEST( all(v8 == Vec3(1, 4, 7)), "Truncation to vector failed!" );

        Vec<float, 0> v9;
        Vec<float, 1> v10;

        Mat3x4 m3(m0, v7);
        Mat4x3 m4(transpose(m0), transpose(v7));
        TEST( m3 == transpose(m4), "Non-square matrix composition constructor wrong." );
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
        TEST( Vec3(v0, 2.0) == v2, "vec2+scalar to vec3 failed!" );
        TEST( Vec3(0.5f, v1) == v3, "scalar+vec2 to vec3 failed!" );
        TEST( IVec4(v1, v0) == v4, "vec2+vec2 to vec4 failed!" );
        TEST( Vec4(v0, 1.0f, 2.0f) == v5, "vec2+2*scalar to vec4 failed!" );
        TEST( Vec4(1.0f, v0, 2.0f) == v6, "scalar+vec2+scalar to vec4 failed!" );
        TEST( Vec4(1.0f, 2.0f, v0) == v7, "2*scalar+vec2 to vec4 failed!" );
        TEST( Vec4(1.0f, v2) == v6, "scalar+vec3 to vec4 failed!" );
        TEST( Vec4(v3, 7.5f) == v8, "vec2+scalar to vec4 failed!" );
    }

    // Test vector composition constructors for row vectors
    {
        Matrix<float, 1, 2> v0(0.0f, 0.5f);
        Matrix<int, 1, 2> v1(1, 2);
        Matrix<float, 1, 3> v2(0.0f, 0.5f, 2.0f);
        Matrix<float, 1, 3> v3(0.5f, 1.0f, 2.0f);
        Matrix<int, 1, 4> v4(1, 2, 0, 0);
        Matrix<float, 1, 4> v5(0.0f, 0.5f, 1.0f, 2.0f);
        Matrix<float, 1, 4> v6(1.0f, 0.0f, 0.5f, 2.0f);
        Matrix<float, 1, 4> v7(1.0f, 2.0f, 0.0f, 0.5f);
        Matrix<float, 1, 4> v8(0.5f, 1.0f, 2.0f, 7.5f);
        TEST( (Matrix<float, 1, 3>(v0, 2.0) == v2), "rvec2+scalar to rvec3 failed!" );
        TEST( (Matrix<float, 1, 3>(0.5f, v1) == v3), "scalar+rvec2 to rvec3 failed!" );
        TEST( (Matrix<int, 1, 4>(v1, v0) == v4), "rvec2+rvec2 to rvec4 failed!" );
        TEST( (Matrix<float, 1, 4>(v0, 1.0f, 2.0f) == v5), "rvec2+2*scalar to rvec4 failed!" );
        TEST( (Matrix<float, 1, 4>(1.0f, v0, 2.0f) == v6), "scalar+rvec2+scalar to rvec4 failed!" );
        TEST( (Matrix<float, 1, 4>(1.0f, 2.0f, v0) == v7), "2*scalar+rvec2 to rvec4 failed!" );
        TEST( (Matrix<float, 1, 4>(1.0f, v2) == v6), "scalar+rvec3 to rvec4 failed!" );
        TEST( (Matrix<float, 1, 4>(v3, 7.5f) == v8), "rvec2+scalar to rvec4 failed!" );
    }

    // ********************************************************************* //
    // Test of boolean function
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(2, 2, 0, 4);
        Matrix<int, 2, 2> m2(2, 1, 0, 4);
        TEST( sum(equal(2, m2)) == 1 , "Scalar operator == or sum() wrong!" );
        TEST( sum(equal(m2, 3)) == 0 , "Scalar operator == or sum() wrong!" );
        TEST( sum(neq(2, m2)) == 3 , "Scalar operator != or sum() wrong!" );
        TEST( sum(neq(m2, 3)) == 4 , "Scalar operator != or sum() wrong!" );
        TEST( any(ei::less(m0, 2)) , "Scalar operator < or any() wrong!" );
        TEST( !any(ei::less(m0, 0)) , "Scalar operator < or any() wrong!" );
        TEST( sum(ei::less(2, m1)) == 1 , "Scalar operator < or sum() wrong!" );
        TEST( !none(lesseq(m0, 2)) , "Scalar operator <= or none() wrong!" );
        TEST( none(lesseq(m0, 0)) , "Scalar operator <= or none() wrong!" );
        TEST( sum(lesseq(2, m1)) == 3 , "Scalar operator <= or sum() wrong!" );
        TEST( all(greatereq(m0, 1)) , "Scalar operator >= or all() wrong!" );
        TEST( !all(greatereq(m1, 2)) , "Scalar operator >= or all() wrong!" );
        TEST( sum(greatereq(4, m1)) == 4 , "Scalar operator >= or sum() wrong!" );
        TEST( any(ei::greater(m0, 2)) , "Scalar operator > or all() wrong!" );
        TEST( none(ei::greater(m0, 4)) , "Scalar operator > or all() wrong!" );
        TEST( sum(ei::greater(4, m1)) == 3 , "Scalar operator > or sum() wrong!" );
        TEST( sum(ei::less(m0, m1)) == 1 , "Matrix operator < or sum() wrong!" );
        TEST( sum(lesseq(m0, m1)) == 3 , "Matrix operator <= or sum() wrong!" );
        TEST( sum(greatereq(m0, m1)) == 3 , "Matrix operator >= or sum() wrong!" );
        TEST( sum(ei::greater(m0, m1)) == 1 , "Matrix operator > or sum() wrong!" );
        TEST( sum(equal(m0, m1)) == 2 , "Matrix operator == or sum() wrong!" );
        TEST( sum(neq(m0, m2)) == 3 , "Matrix operator != or sum() wrong!" );
    }


    // ********************************************************************* //
    // Test of access operators and element constructors (2-4 elements)
    {
        Matrix<int, 2, 1> m0(1, 2);
        Matrix<int, 1, 3> m1(1, 2, 3);
        Matrix<int, 2, 2> m2(1, 2, 3, 4);
        Matrix<int, 1, 4> m3(3, 1, 2, 5);
        Matrix<int, 1, 4> m4(1, 2, 5, 5);
        Matrix<int, 4, 1> m5(1, 1, 2, 2);
        m1(0,2) = 5;
        m2[2] = 5;
        TEST( m0[0] == 1 && m0[1] == 2, "Operator or constructor wrong!" );
        TEST( m1(0,0) == 1 && m1(0,1) == 2 && m1[2] == 5, "Operator or constructor wrong!" );
        TEST( m2(0,1) == 2 && m2(1,0) == 5 && m2(1,1) == 4, "Operator or constructor wrong!" );
        TEST( (const_cast<const Matrix<int, 2, 2>&>(m2))[1] == 2, "Read-only operator wrong!" );
        TEST( (const_cast<const Matrix<int, 2, 2>&>(m2))(1,0) == 5, "Read-only operator wrong!" );
        TEST( (m3.subrow<1,4>() == m1), "Range access for row vectors failed!" );
        m3.subrow<0,3>() = m1;
        TEST( m3 == m4, "Range access for row vectors failed!" );
        TEST( (m5.subcol<1,3>() == m0), "Range access for column vectors failed!" );

        Matrix<int, 3, 4> m6(1, 2, 3, 4,
                             5, 6, 7, 8,
                             9, 0, 1, 2);
        Matrix<int, 2, 2> m6s1 = m6.submat<1,3, 0,2>();
        Matrix<int, 3, 1> m6s2 = m6.submat<0,3, 3,4>();
        TEST( (m6s1 == Matrix<int, 2, 2>(5, 6, 9, 0)), "submat()->2x2 access operator failed!" );
        TEST( (m6s2 == Matrix<int, 3, 1>(4, 8, 2)), "submat()->3x1 access operator failed!" );
    }

    // ********************************************************************* //
    // Test + and -
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(5, 6, 7, 8);
        Matrix<int, 2, 2> m2(6, 8, 10, 12);
        Matrix<float, 2, 2> m3(3.0f, 1.0f, 4.0f, 1.0f);
        Matrix<float, 2, 2> m4(4.0f, 3.0f, 7.0f, 5.0f);
        TEST( m0 + m1 == m2, "Component wise addition failed!" );
        TEST( m0 + m3 == m4, "Component wise addition failed!" );
        TEST( m2 - m1 == m0, "Component wise subtraction failed!" );
        TEST( m0 - m4 == -m3, "Component wise subtraction or unary minus failed!" );
    }

    // ********************************************************************* //
    // Test |, &, ^, % and ~
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(5, 6, 7, 8);
        Matrix<int, 2, 2> m2(5, 6, 7, 12);
        Matrix<int, 2, 2> m3(1, 2, 3, 0);
        Matrix<int, 2, 2> m4(4, 4, 4, 12);
        Matrix<int, 2, 2> m5(0xfffffffe, 0xfffffffd, 0xfffffffc, 0xfffffffb);
        Matrix<int, 2, 2> m6(1, 2, 3, 8);
        Matrix<int, 2, 2> m7(4, 8, 12, 16);
        Matrix<int, 2, 2> m8(2, 3, 3, 4);
        TEST( (m0 | m1) == m2, "Component wise | failed!" );
        TEST( (m0 & m1) == m3, "Component wise & failed!" );
        TEST( (m0 ^ m1) == m4, "Component wise ^ failed!" );
        TEST( (m1 % m4) == m6, "Component wise % failed!" );
        TEST( (~m0) == m5, "Component wise ~ failed!" );
        TEST( m0 << 2 == m7, "Component wise << failed!" );
        TEST( m1 >> 1 == m8, "Component wise >> failed!" );
        Matrix<int, 2, 2> t = m0; t <<= 2;
        TEST( t == m7, "Component wise << failed!" );
        t = m1; t >>= 1;
        TEST( t == m8, "Component wise >> failed!" );
        // Generates an error message (expected)
        //Matrix<float, 2, 2> mf0(3.0f, 1.0f, 4.0f, 1.0f);
        //Matrix<float, 2, 2> mf1(4.0f, 3.0f, 7.0f, 5.0f);
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
        TEST( v1 * v2 == v3, "Component wise multiplication failed!" );
        TEST( v0 * v4 == v5, "Component wise multiplication failed!" );
        TEST( v1 / v6 == v7, "Component wise division failed!" );
        TEST( v0 == v5 / v4, "Component wise division failed!" );
        TEST( v1 * v0 == -2.0f, "Row times column vector should be a scalar!" );
        TEST( v0 * v1 == m0, "Column times row vector should be a matrix!" );
        TEST( v4 * v2 == m1, "Column times row vector should be a matrix!" );
        TEST( m0 * m1 == m2, "Matrix multiplication wrong!" );
        TEST( (transpose(v0) * m0 == Matrix<float, 1, 3>(14.0f, 0.0f, -14.0f)), "Vector * Matrix multiplication invalid!" );
        TEST( m0 * v0 == Vec3(-2.0f, -4.0f, -6.0f), "Matrix * Vector multiplication invalid!" );
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
        TEST( m1 * m0 == m2, "Constructor or matrix multiplication wrong!" );
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
        TEST( m4 * m3 == m5, "Constructor or matrix multiplication wrong!" );
        Matrix<int, 3, 4> m6(1, 2, 3, 4,
                             2, 3, 4, 5,
                             0, 0, 0, 0);
        Matrix<int, 4, 3> m7(1, 2, 3,
                             2, 3, 4,
                             3, 4, 5,
                             4, 5, 6);
        TEST( m7 * m6 == m5, "Constructor or matrix multiplication wrong!" );
    }

    // ********************************************************************* //
    // Test scalar operators +, *, -, /
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(3, 4, 5, 6);
        Matrix<float, 1, 2> v0(1.0f, 1.0f);
        Matrix<float, 2, 1> v1(0.0f, 1.0f);
        Matrix<float, 2, 2> m2(2.0f, 2.0f, 4.0f, 4.0f);
        TEST( m0 + 2 == m1, "Adding a scalar failed!" );
        TEST( 2 + m0 * 1 == m1, "Multiplying or adding a scalar failed!" );
        TEST( 2.0f * v0 * 2.0f * v1 + 0.5f == 4.5f, "Multiplying a scalar failed!" );
        TEST( (v1 + 1.0f) * 2.0f * v0 == m2, "Mixed scalar vector operation failed!" );
        TEST( m0 == m1 - 2, "Scalar subtraction failed!" );
        TEST( 2 + m0 / 1 == m1, "Scalar addition or division failed!" );
        TEST( (2.0f / v0 / 2.0f) * v1 - 0.5f == 0.5f, "Mixed scalar vector operation (/,*,-) failed!" );
        TEST( (v1 - 0.5f) * 2.0f * v0 + 3 == m2, "Mixed scalar vector operation (-, *) failed!" );
        TEST( (4.0f - m2) / 2.0f == (1.0f - v1) * v0, "Mixed scalar vector operation (/,*,-) failed!" );
    }

    // ********************************************************************* //
    // Test scalar |, &, ^, %
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(5, 6, 7, 4);
        Matrix<int, 2, 2> m2(0, 0, 0, 4);
        Matrix<int, 2, 2> m3(5, 6, 7, 0);
        Matrix<int, 2, 2> m4(1, 2, 0, 1);
        Matrix<int, 2, 2> m5(0, 1, 0, 3);
        TEST( (m0 | 4) == m1, "Component wise | failed!" );
        TEST( (4 | m0) == m1, "Component wise | failed!" );
        TEST( (m0 & 4) == m2, "Component wise & failed!" );
        TEST( (4 & m0) == m2, "Component wise & failed!" );
        TEST( (m0 ^ 4) == m3, "Component wise ^ failed!" );
        TEST( (4 ^ m0) == m3, "Component wise ^ failed!" );
        TEST( (m0 % 3) == m4, "Component wise % failed!" );
        TEST( (3 % m0) == m5, "Component wise % failed!" );
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
        Mat2x2 m3(2.0f, 2.0f, 4.0f, 4.0f);
        Mat2x2 m4(5.0f, 6.0f, 9.0f, 10.0f);
        Vec2 v0(1.0f, 1.0f);
        Vec2 v1(0.0f, 1.0f);
        Vec2 v2(0.0f, 1.0f);
        Matrix<float, 1, 2> vr0(1.0f, 1.0f);
        Matrix<float, 1, 2> vr1(0.0f, 1.0f);
        Matrix<float, 1, 2> vr2(0.0f, 1.0f);
        TEST( (m0 += m1) == m2, "Component wise += failed!" );
        m2 -= m0;
        TEST( (m2 == Matrix<int,2,2>(0)), "Component wise -= failed!" );
        m3 += m1;
        TEST( m3 == m4, "Component wise += with different types failed!" );
        v0 *= v1;
        TEST( v0 == v1, "Component wise *= for vectors failed!" );
        (v1 += v0) /= v0 + 1.0f;
        TEST( v1 == v2, "Component wise /= for vectors failed!" );
        vr0 *= vr1;
        TEST( vr0 == vr1, "Component wise *= for row vectors failed!" );
        (vr1 += vr0) /= vr0 + 1.0f;
        TEST( vr1 == vr2, "Component wise /= for row vectors failed!" );

        TEST( (vr1 += 1.0f) == vr2 + 1, "Matrix-Scalar += failed!");
        TEST( (m3 -= 0.5f) == m4 - 0.5f, "Matrix-Scalar -= failed!");
        TEST( (v1 *= 2.0f) == v2 / 0.5f, "Matrix-Scalar *= failed!");
        TEST( (v1 /= 2.0f) == v2, "Matrix-Scalar /= failed!");

        Mat2x2 m5(72.0f, 82.0f, 128.0f, 146.0f);
        m3 *= m4;
        TEST( m3 == m5, "Matrix *= Matrix failed!");
    }

    // ********************************************************************* //
    // Test assignment operators |=, &=, ^=, %=
    {
        Matrix<int, 2, 2> m00(1, 2, 3, 4);
        Matrix<int, 2, 2> m01(1, 2, 3, 4);
        Matrix<int, 2, 2> m02(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(5, 6, 7, 8);
        Matrix<int, 2, 2> m2(5, 6, 7, 12);
        Matrix<int, 2, 2> m3(1, 2, 3, 0);
        Matrix<int, 2, 2> m4(4, 4, 4, 12);
        Matrix<int, 2, 2> m5(1, 2, 3, 8);
        Matrix<int, 2, 2> m6(5, 7, 7, 13);
        Matrix<int, 2, 2> m7(1, 0, 1, 0);
        Matrix<int, 2, 2> m8(12, 12, 12, 4);
        Matrix<int, 2, 2> m9(1, 2, 0, 2);
        TEST( (m00 |= m1) == m2, "Matrix self assigning |= failed!" );
        TEST( (m01 &= m1) == m3, "Matrix self assigning &= failed!" );
        TEST( (m02 ^= m1) == m4, "Matrix self assigning ^= failed!" );
        TEST( (m1 %= m4) == m5, "Matrix self assigning %= failed!" );
        TEST( (m2 |= 1) == m6, "Matrix-scalar self assigning |= failed!" );
        TEST( (m3 &= 1) == m7, "Matrix-scalar self assigning &= failed!" );
        TEST( (m4 ^= 8) == m8, "Matrix-scalar self assigning ^= failed!" );
        TEST( (m5 %= 3) == m9, "Matrix-scalar self assigning %= failed!" );
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
        TEST( min(v0, v1) == v2, "Component wise minimum wrong!" );
        TEST( max(v0, v1) == v3, "Component wise maximum wrong!" );
        TEST( min(v0) == -1.0f, "Minimum of a vector wrong!" );
        TEST( min(v1) == 0.0f, "Minimum of a vector wrong!" );
        TEST( max(v2) == 0.25f, "Maximum of a vector wrong!" );
        TEST( max(v3) == 2.5f, "Maximum of a vector wrong!" );
        TEST( clamp(v0, v2, v1) == v5, "Clamp of a vector in box boundaries wrong!" );
        TEST( clamp(v0, 0.0f, 3.0f) == v6, "Clamp of a vector in interval wrong!" );
        TEST( saturate(v0 * 4.0f) == v7, "Saturate of a vector wrong!" );
        TEST( sum(v1) == 4.0f, "Sum of a vector wrong!" );
        TEST( prod(v0) == 0.0625f, "Product of a vector wrong!" );
        TEST( avg(v3) == 1.125f, "Average of a vector wrong!" );
        TEST( abs(v0) == v4, "Absolute value of a vector wrong!" );
        TEST( sign(v1) == Vec4(1.0f, 1.0f, 1.0f, 0.0f), "Sign values of a vector wrong!" );
        TEST( sgn(v2) == Vec4(-1.0f, -1.0f, 1.0f, 1.0f), "Sgn values of a vector wrong!" );
    }

    // ********************************************************************* //
    // Test round, ceil, floor, mod, sqrt, pow, log, log2
    {
        Vec4 v0(1.5f, 2.0f, -0.5f, 1.414f);
        IVec4 v1(2, 2, 0, 1);    // round(v0)
        IVec4 v2(2, 2, 0, 2);    // ceil(v0)
        IVec4 v3(1, 2, -1, 1);   // floor(v0)
        Vec4 v6(5.5f, -1.5f, -0.5f, 2.5f);
        IVec4 v7(6, -2, 0, 2);   // round(v6)
        Vec4 v8(0.5f, 0.0f, 0.5f, 0.414f); // mod(v0, 1.0f)
        TEST( round(v0) == v1, "Rounding v0 is wrong!" );
        TEST( ceil(v0) == v2, "Rounding v0 up is wrong!" );
        TEST( floor(v0) == v3, "Rounding v0 down is wrong!" );
        TEST( round(v6) == v7, "Rounding v6 is wrong!" );
        TEST( approx(mod(v0, 1.0f), v8), "mod(v0,1) is wrong!" );

        Vec3 v4(1.0f, 4.0f, 2.0f);
        Vec3 v5(1.0f, 2.0f, PHYTAGORAS);
        TEST( sqrt(v4) == v5, "sqrt of v4 is wrong!" );
        TEST( pow(v4,0.5f) == v5, "v4 to the power of 0.5 is wrong!" );
        TEST( approx(pow(v5,2.0f), v4), "v5 to the power of 2 is wrong!" );
        TEST( log2(v4) == Vec3(0.0f, 2.0f, 1.0f), "log2(v4) is wrong!" );
        TEST( approx(log(v4), Vec3(0.0f, 1.386294361f, 0.693147181f)), "log2(v4) is wrong!" );
    }

    // ********************************************************************* //
    // Test frac, intfrac, floorfrac
    {
        const Vec4 v0(0.6f, 1.7f, -0.7f, -2.1f);
        const Vec4 v1(0.6f, 0.7f, -0.7f, -0.1f);
        const IVec4 v2(0, 1, 0, -2);
        const IVec4 v3(0, 1, -1, -3);
        const Vec4 v4(0.6f, 0.7f, 0.3f, 0.9f);
        IVec4 tmp;
        TEST( approx(frac(v0), v1), "Vector overerload of frac() wrong!" );
        TEST( approx(intfrac(v0, tmp), v1) && tmp == v2, "Vector overerload of intfrac() wrong!" );
        TEST( approx(floorfrac(v0, tmp), v4) && tmp == v3, "Vector overerload of floorfrac() wrong!" );
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
        TEST( lerp(v0, v2, 0.5f) == v1, "Linear interpolation failed." );
        TEST( bilerp(v0, v1, v1, v2, 0.5f, 0.5f) == v1, "Bilinear interpolation failed." );
        TEST( lerp(1, 4, 0.5f) == 2.5f, "Linear interpolation failed." );
        TEST( lerp(v4, v5, 0.25f) == v7, "Linear interpolation by 0.25 failed." );
        TEST( bilerp(1, 3, 5, 7, 0.5f, 0.5f) == 4.0f, "Bilinear interpolation failed." );
        TEST( approx(slerp(v4, v5, 0.25f), v6), "Spherical linear interpolation failed." );
        TEST( approx(dot(slerp(v1, v3, 0.5f), v1), 0.0f), "Spherical linear interpolation for degenerated case failed." );
        TEST( slerp(v1, v1, 0.125f) == v1, "Spherical linear interpolation for two equal vectors failed." );
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
        TEST( cross(v0, v1) == v2, "3D crossproduct wrong." );
        TEST( cross(v1, v0) == -v2, "3D crossproduct wrong." );
        TEST( cross(v5, v6) == v7, "3D crossproduct wrong." );
        TEST( cross(v6, v5) == -v7, "3D crossproduct wrong." );
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
        TEST( m0 == identity2x2(), "2x2 identity matrix wrong." );
        TEST( m1 == identity3x3(), "3x3 identity matrix wrong." );
        TEST( m2 == identity4x4(), "4x4 identity matrix wrong." );
        TEST( (m3 == ei::identity<float,5>()), "5x5 identity matrix wrong." );
        TEST( m4 == m5, "3x3 diagonal matrix wrong." );
        auto m61 = diag(Vec<float,6>{1});
        auto m6 = ei::identity<float,6>();

        /*// Test script used to optimize identity()
        Mat3x3 ret;
        uint64 a = ticks();
        for(uint i=0; i<1000000; ++i) {
            ret = ei::identity<float,3>();
            doNotOptimizeAway(ret);
        }
        uint64 b = ticks();
        cerr << (int)(b-a) << endl;//*/
    }

    // ********************************************************************* //
    // Test transpose
    {
        Mat2x2 m0(1.0f, -1.0f, 0.0f, 2.0f);
        Mat2x2 m1(1.0f, 0.0f, -1.0f, 2.0f);
        Matrix<float,1,3> v0(1.0f, 0.0f, 0.0f);
        Vec3 v1(1.0f, 0.0f, 0.0f);

        TEST( m0 == transpose(m1), "m1 transposed is unequal m0!" );
        TEST( v0 == transpose(v1), "v1 transposed is unequal v0!" );
    }

    // ********************************************************************* //
    // Test axis matrix
    {
        Mat3x3 m0 = axis(normalize(Vec3(1.0f, 1.0f, 0.0f)), normalize(Vec3(-1.0f, 1.0f, 0.0f)), Vec3(0.0f, 0.0f, 1.0f));
        Vec3 v0(0.0f, 0.0f, 5.0f);
        Vec3 v1(2.0f, 0.0f, 0.0f);
        Vec3 v2(1.414213562f, -1.414213562f, 0.0f);
        TEST( m0 * v0 == v0, "v0 should not be transformed by m0!" );
        TEST( m0 * v1 == v2, "v1 should be transformed into v2 by m0!" );
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
        TEST( approx(rotation(PI/2, 0.0f, 0.0f) * vy, vz, 4e-6f), "Rotation matrix from Euler angles (PI/2, 0, 0) invalid!" );
        TEST( approx(rotation(-PI/2, 0.0f, 0.0f) * vz, vy, 4e-6f), "Rotation matrix from Euler angles (-PI/2, 0, 0) invalid!" );
        TEST( approx(rotation(0.0f, -PI/2, 0.0f) * vx, vz, 4e-6f), "Rotation matrix from Euler angles (0, -PI/2, 0) invalid!" );
        TEST( approx(rotation(0.0f, PI/2, 0.0f) * vz, vx, 4e-6f), "Rotation matrix from Euler angles (0, PI/2, 0) invalid!" );
        TEST( approx(rotation(0.0f, 0.0f, PI/2) * vx, vy, 4e-6f), "Rotation matrix from Euler angles (0, 0, PI/2) invalid!" );
        TEST( approx(rotation(0.0f, 0.0f, -PI/2) * vy, vx, 4e-6f), "Rotation matrix from Euler angles (0, 0, -PI/2) invalid!" );

        TEST( approx(rotation(vx, vx), identity3x3()), "Rotation matrix from vector to itself invalid!" );
        TEST( approx(rotation(vy, vz), rotation(PI/2, 0.0f, 0.0f), 4e-6f), "Rotation matrix from vector to vector invalid!" );

        TEST( approx(housholder(Vec3(-1.0f, 1.0f, 0.0f)) * vy, vx), "Housholder matrix invalid!" );
        TEST( approx(housholder(Vec3(1.0f, -1.0f, 0.0f)) * vy, vx), "Housholder matrix invalid!" );
    }

    // ********************************************************************* //
    // Test camera transformations
    {
        Mat4x4 m0 = Mat4x4(lookAt(Vec3(1.0f, 0.0f, 0.0f)));
        Mat4x4 m1 = camera(Vec3(1.0f, 1.0f, 1.0f), Vec3(0.0f));
        Mat4x4 m2 = invert(m1);
        Vec3 v0(0.0f);
        Vec3 v1(0.0f, 0.0f, 1.0f);

        TEST( approx(transform(v0, m0), v0), "lookAt() should not contain any translation!" );
        TEST( approx(transform(v1, m0), Vec3(-1.0f, 0.0f, 0.0f)), "lookAt() rotates wrongly!" );

        TEST( approx(transform(v0, m1), Vec3(0.0f, 0.0f, sqrt(3.0f))), "camera() transformation not as expected!" );
        TEST( approx(transform(v0, m2), Vec3(1.0f, 1.0f, 1.0f)), "camera() inverse transformation not as expected!" );
        TEST( approx(transform(v1 * sqrt(3.0f), m2), v0), "camera() inverse transformation not as expected!" );

        TEST( approx(transformDir(v0, m1), v0), "transformDir() should not change the length!" );
        TEST( approx(transformDir(transformDir(v1, m1), m2), v1), "transformDir() should not change the length!" );
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

        TEST( sphericalCoords(v0) == s0, "Spherical coordinates of v0 wrong!" );
        TEST( sphericalCoords(v1) == s1, "Spherical coordinates of v1 wrong!" );
        TEST( sphericalCoords(v2) == s2, "Spherical coordinates of v2 wrong!" );
        TEST( sphericalCoords(v3) == s3, "Spherical coordinates of v3 wrong!" );
        TEST( sphericalCoords(v4) == s4, "Spherical coordinates of v4 wrong!" );
        TEST( sphericalCoords(v5) == s5, "Spherical coordinates of v5 wrong!" );

        TEST( v0 == cartesianCoords(s0), "Cartesian coordinates of s0 wrong!" );
        TEST( approx(v1, cartesianCoords(s1)), "Cartesian coordinates of s1 wrong!" );
        TEST( approx(v2, cartesianCoords(s2)), "Cartesian coordinates of s2 wrong!" );
        TEST( approx(v3, cartesianCoords(s3)), "Cartesian coordinates of s3 wrong!" );
        TEST( approx(v4, cartesianCoords(s4)), "Cartesian coordinates of s4 wrong!" );
        TEST( approx(v5, cartesianCoords(s5), 4e-6f), "Cartesian coordinates of s5 wrong!" );
        TEST( approx(v6, cartesianCoords(sphericalCoords(v6)), 1e-5f), "5D vector coordinate transformation to spherical and back failed!" );
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
        TEST( invert(A0) == X, "Manual inversion and invert function have different results!" );
        X = X * A0;
        TEST( approx(X, identity3x3(), 4e-6f), "3x3 Matrix inverse bad!");
        TEST( approx(invert(A2) * A2, identity4x4(), 4e-5f), "4x4 Matrix inverse bad!");
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

    // ********************************************************************* //
    // Test orthonarmalization
    {
        Mat3x3 m0(6.0f, 1.0f, 1.0f, 4.0f, -2.0f, 5.0f, 2.0f, 8.0f, 7.0f);
        Mat3x3 m1(1.0f, 0.5f, 0.5f, 0.0f, 2.0f, 2.0f, 0.0f, 0.0f, 0.0f);
        TEST( orthonormalize(m0), "m0 is not linear dependent (orthonormalization should succeed)!" );
        m0 = transpose(m0);
        TEST( approx(lensq(m0(0)), 1.0f), "Non normalized vector x in orthonormalized base!" );
        TEST( approx(lensq(m0(1)), 1.0f), "Non normalized vector y in orthonormalized base!" );
        TEST( approx(lensq(m0(2)), 1.0f), "Non normalized vector z in orthonormalized base!" );
        TEST( approx(lensq(cross(m0(0), m0(1))), 1.0f), "x and y are not orthogonal!" );
        TEST( approx(lensq(cross(m0(0), m0(2))), 1.0f), "x and z are not orthogonal!" );
        TEST( approx(lensq(cross(m0(2), m0(1))), 1.0f), "y and z are not orthogonal!" );
        TEST( !orthonormalize(m1), "m0 is linear dependent (orthonormalization should fail)!" );

        Vec3 v0(1.0f, 0.0f, 2.0f);
        Vec3 v1(1.0f, 0.5f, 0.0f);
        Vec3 v2(0.0f, 2.0f, -2.0f);
        orthonormalize(v0, v1, v2);
        TEST( approx(lensq(v0), 1.0f), "Non normalized vector v0 after orthonormalization!" );
        TEST( approx(lensq(v1), 1.0f), "Non normalized vector v1 after orthonormalization!" );
        TEST( approx(lensq(v2), 1.0f), "Non normalized vector v2 after orthonormalization!" );
        TEST( approx(lensq(cross(v0, v1)), 1.0f), "v0 and v1 are not orthogonal!" );
        TEST( approx(lensq(cross(v0, v2)), 1.0f), "v0 and v2 are not orthogonal!" );
        TEST( approx(lensq(cross(v2, v1)), 1.0f), "v1 and v2 are not orthogonal!" );
    }


    // ********************************************************************* //
    // Test eigen vector/value decompositions
    {
        Mat2x2 Q0(0.0f, 1.0f, -1.0f, 0.0f);
        Mat2x2 Q1(1.0f, 1.0f, -1.0f, 2.0f);
        orthonormalize(Q1);
        Vec2 d0(2.0f, 3.0f);
        Vec2 d1(-1.0f, 0.0f);
        Mat2x2 A0 = transpose(Q0) * diag(d0) * Q0;
        Mat2x2 A1 = transpose(Q1) * diag(d1) * Q1;
        Mat2x2 A2(0.5f, 0.0f, 0.0f, 1.0f); // diagonal matrix
        Vec2 vtmp; Mat2x2 Qtmp;
        int itn = decomposeQl(A0, Qtmp, vtmp);
        TEST(approx(vtmp.x, 3.0f), "Eigenvalue 0 of A0 is wrong!");
        TEST(approx(vtmp.y, 2.0f), "Eigenvalue 1 of A0 is wrong!");
        TEST(approx(A0, transpose(Qtmp) * diag(vtmp) * Qtmp), "Spectral decomposition of 2x2 A0 failed!");
        itn = decomposeQl(A1, Qtmp, vtmp);
        TEST(approx(vtmp.x, 0.0f), "Eigenvalue 0 of A0 is wrong!");
        TEST(approx(vtmp.y, -1.0f), "Eigenvalue 1 of A0 is wrong!");
        TEST(approx(A1, transpose(Qtmp) * diag(vtmp) * Qtmp), "Spectral decomposition of 2x2 A1 failed!");
        decomposeQl(A2, Qtmp, vtmp);
        TEST(vtmp == Vec2(1.0f, 0.5f), "Eigenvalues of A2 are wrong!");
        TEST(Qtmp == Mat2x2(0.0f, 1.0f, 1.0f, 0.0f), "Eigenvectors of A2 are wrong!");
        decomposeQl(identity2x2(), Qtmp, vtmp);
        TEST(vtmp == Vec2(1.0f, 1.0f), "Eigenvalues of identity2x2 are wrong!");
        TEST(Qtmp == Mat2x2(1.0f, 0.0f, 0.0f, 1.0f), "Eigenvectors of identity2x2 are wrong!");
    }{
        Mat3x3 Q1(1.0f, 2.0f, 3.0f, 2.0f, 1.0f, 0.5f, 0.0f, -1.0f, 1.0f);
        orthonormalize(Q1);
        const Vec3 d1(0.7f, 1.0f, 1.1f);
        const Mat3x3 A1 = transpose(Q1) * diag(d1) * Q1;	// Arbitrary symmetric matrix
        const Mat3x3 A0(1.0f, 0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.5f);	// diagonal matrix
        Vec3 vtmp; Mat3x3 Qtmp;

        int itn = decomposeQl(A1, Qtmp, vtmp);
        TEST(approx(vtmp.x, 1.1f), "Eigenvalue 0 of A1 is wrong!");
        TEST(approx(vtmp.y, 1.0f), "Eigenvalue 1 of A1 is wrong!");
        TEST(approx(vtmp.z, 0.7f), "Eigenvalue 2 of A1 is wrong!");
        TEST(approx(A1, transpose(Qtmp) * diag(vtmp) * Qtmp), "Spectral decomposition of A1 failed!");
        itn = decomposeQlIter(A1, Qtmp, vtmp);
        TEST(approx(vtmp.x, 1.1f), "Eigenvalue 0 of A1 is wrong!");
        TEST(approx(vtmp.y, 1.0f), "Eigenvalue 1 of A1 is wrong!");
        TEST(approx(vtmp.z, 0.7f), "Eigenvalue 2 of A1 is wrong!");
        TEST(approx(A1, transpose(Qtmp) * diag(vtmp) * Qtmp), "Spectral decomposition of A1 failed!");

        decomposeQl(A0, Qtmp, vtmp);
        TEST(vtmp == Vec3(2.0f, 1.0f, 0.5f), "Eigenvalues of A0 are wrong!");
        TEST(approx(A0, transpose(Qtmp) * diag(vtmp) * Qtmp), "Spectral decomposition of A0 failed!");
        itn = decomposeQlIter(A0, Qtmp, vtmp);
        TEST(vtmp == Vec3(2.0f, 1.0f, 0.5f), "Eigenvalues of A0 are wrong!");
        TEST(approx(A0, transpose(Qtmp) * diag(vtmp) * Qtmp), "Spectral decomposition of A0 failed!");

        // More difficult case for power iteration
        const Vec3 d2(1.0f, -1.001f, 1.0001f);
        const Mat3x3 A2 = transpose(Q1) * diag(d2) * Q1;
        itn = decomposeQl(A2, Qtmp, vtmp);
        TEST(approx(vtmp.x, 1.0001f), "Eigenvalue 0 of A2 is wrong!");
        TEST(approx(vtmp.y, 1.0f), "Eigenvalue 1 of A2 is wrong!");
        TEST(approx(vtmp.z, -1.001f), "Eigenvalue 2 of A2 is wrong!");
        TEST(approx(A2, transpose(Qtmp) * diag(vtmp) * Qtmp, 5e-6f), "Spectral decomposition of A2 failed!");
        itn = decomposeQlIter(A2, Qtmp, vtmp);
        TEST(approx(vtmp.x, 1.0001f), "Eigenvalue 0 of A2 is wrong!");
        TEST(approx(vtmp.y, 1.0f), "Eigenvalue 1 of A2 is wrong!");
        TEST(approx(vtmp.z, -1.001f), "Eigenvalue 2 of A2 is wrong!");
        TEST(approx(A2, transpose(Qtmp) * diag(vtmp) * Qtmp, 2e-5f), "Spectral decomposition of A2 failed!");

        // Degenerated cases
        const Mat3x3 A3(1.0f, 0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 1.0f);
        decomposeQl(A3, Qtmp, vtmp);
        TEST(vtmp == Vec3(2.0f, 1.0f, 1.0f), "Eigenvalues of A3 are wrong!");
        TEST(approx(A3, transpose(Qtmp) * diag(vtmp) * Qtmp), "Spectral decomposition of A3 failed!");
        decomposeQlIter(A3, Qtmp, vtmp);
        TEST(vtmp == Vec3(2.0f, 1.0f, 1.0f), "Eigenvalues of A3 are wrong!");
        TEST(approx(A3, transpose(Qtmp) * diag(vtmp) * Qtmp), "Spectral decomposition of A3 failed!");

        decomposeQl(identity3x3(), Qtmp, vtmp);
        TEST(vtmp == Vec3(1.0f, 1.0f, 1.0f), "Eigenvalues of identity3x3 are wrong!");
        TEST(Qtmp == Mat3x3(1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f), "Eigenvectors of identity3x3 are wrong!");
        itn = decomposeQlIter(identity3x3(), Qtmp, vtmp);
        TEST(vtmp == Vec3(1.0f, 1.0f, 1.0f), "Eigenvalues of identity3x3 are wrong!");
        TEST(Qtmp == Mat3x3(1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f), "Eigenvectors of identity3x3 are wrong!");
    }{
        const Mat4x4 A0 {   4.0f,  -30.0f,    60.0f,   -35.0f,
                          -30.0f,  300.0f,  -675.0f,   420.0f,
                           60.0f, -675.0f,  1620.0f, -1050.0f,
                          -35.0f,  420.0f, -1050.0f,   700.0f };
        Vec4 vtmp; Mat4x4 Qtmp;
        int itn = decomposeQlIter(A0, Qtmp, vtmp);
        TEST(approx(vtmp, Vec4(2585.2538109f, 37.101491365f, 1.4780548448f, 0.16664286117f), 1e-5f), "Eigenvalues of 4x4 matrix wrong!");
        TEST(approx(Qtmp(0), Matrix<float,1,4>( 0.02919332f, -0.3287121f,  0.7914111f, -0.5145527f)), "Eigenvector 0 of 4x4 matrix wrong!");
        TEST(approx(Qtmp(1), Matrix<float,1,4>(-0.17918629f,  0.7419178f, -0.1002281f, -0.6382825f)), "Eigenvector 1 of 4x4 matrix wrong!");
        TEST(approx(Qtmp(2), Matrix<float,1,4>(-0.58207570f,  0.3705022f,  0.5095786f,  0.5140483f), 1e-5f), "Eigenvector 2 of 4x4 matrix wrong!");
        TEST(approx(Qtmp(3), Matrix<float,1,4>( 0.79260829f,  0.4519231f,  0.3224164f,  0.2521612f), 1e-5f), "Eigenvector 3 of 4x4 matrix wrong!");
    }

    // ********************************************************************* //
    // Test Cholesky decompositions
    {
        Mat2x2 L0(1.0f, 0.0f, 2.0f, 3.0f);
        Mat2x2 L1(1.0f, 0.0f, 1.0f, 500.0f);
        Mat2x2 A0 = L0 * transpose(L0);
        Mat2x2 A1 = L1 * transpose(L1);
        Mat2x2 A2(1.0f, 1.0f, 1.0f, -1.0f);
        Mat2x2 t;
        decomposeCholesky(A0, t);
        TEST(approx(t, L0), "Cholesky decomposition of A0 failed!");
        decomposeCholesky(A1, t);
        TEST(approx(t, L1), "Cholesky decomposition of A1 failed!");
        TEST(!decomposeCholesky(A2, t), "Cholesky decomposition of A2 should return false (A2 is not positive definite)!");

        Mat3x3 L3(2.0f, 0.0f, 0.0f,
                  1.0f, 0.5f, 0.0f,
                  4.0f, 8.0f, 1.0f);
        Mat3x3 L4(16.0f, 0.0f, 0.0f,
                  -1.0f, 0.125f, 0.0f,
                  400.0f, 0.0f, 1.0f);
        Mat3x3 A3 = L3 * transpose(L3);
        Mat3x3 A4 = L4 * transpose(L4);
        Mat3x3 A5(1.0f, 1.0f, 1.0f, -1.0f, -1.0f, -1.0f, 2.0f, 2.0f, 2.0f);
        Mat3x3 t2;
        decomposeCholesky(A3, t2);
        TEST(approx(t2, L3), "Cholesky decomposition of A3 failed!");
        decomposeCholesky(A4, t2);
        TEST(approx(t2, L4), "Cholesky decomposition of A4 failed!");
        TEST(!decomposeCholesky(A5, t2), "Cholesky decomposition of A5 should return false (A5 is not positive definite)!");

        Mat4x4 L6(2.0f, 0.0f, 0.0f, 0.0f,
                  1.0f, 0.5f, 0.0f, 0.0f,
                  4.0f, 8.0f, 1.0f, 0.0f,
                  2.0f, 1.0f, 8.0f, 0.25f);
        Mat4x4 L7(16.0f, 0.0f, 0.0f, 0.0f,
                  -1.0f, 0.125f, 0.0f, 0.0f,
                  400.0f, 0.0f, 1.0f, 0.0f,
                  0.0f, -128.0f, 2.0f, 1024.0f);
        Mat4x4 A6 = L6 * transpose(L6);
        Mat4x4 A7 = L7 * transpose(L7);
        Mat4x4 A8(1.0f, 1.0f, 1.0f, 1.0f, -1.0f, -1.0f, -1.0f, -1.0f, 2.0f, 2.0f, 2.0f, 2.0f, 8.0f, -8.0f, 8.0f, -8.0f);
        Mat4x4 t3;
        decomposeCholesky(A6, t3);
        TEST(approx(t3, L6), "Cholesky decomposition of A6 failed!");
        decomposeCholesky(A7, t3);
        TEST(approx(t3, L7), "Cholesky decomposition of A7 failed!");
        TEST(!decomposeCholesky(A8, t3), "Cholesky decomposition of A8 should return false (A8 is not positive definite)!");
    }

    return result;
}
