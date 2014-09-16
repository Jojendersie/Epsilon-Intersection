#include "gam/matrix.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace gam;
using namespace std;

bool test_matrix()
{
    bool result = true;

    // ********************************************************************* //
    // Test type size
    TEST( sizeof(Matrix<byte, 1, 2>) == 2,
        "sizeof(Matrix<byte, 1, 2>) is " << sizeof(Matrix<byte, 1, 2>) << " instead of 2\n"
    );
    TEST( sizeof(Matrix<double, 2, 1>) == 16,
        "sizeof(Matrix<double, 2, 1>) is " << sizeof(Matrix<double, 2, 1>) << " instead of 16\n"
    );
    TEST( sizeof(Matrix<float, 1, 3>) == 12,
        "sizeof(Matrix<float, 1, 3>) is " << sizeof(Matrix<float, 1, 3>) << " instead of 12\n"
    );
    TEST( sizeof(Matrix<byte, 3, 1>) == 3,
        "sizeof(Matrix<byte, 3, 1>) is " << sizeof(Matrix<byte, 3, 1>) << " instead of 3\n"
    );
    TEST( sizeof(Matrix<float, 4, 1>) == 16,
        "sizeof(Matrix<float, 4, 1>) is " << sizeof(Matrix<float, 4, 1>) << " instead of 16\n"
    );
    TEST( sizeof(Matrix<int16, 1, 4>) == 8,
        "sizeof(Matrix<int16, 1, 4>) is " << sizeof(Matrix<int16, 1, 4>) << " instead of 8\n"
    );
    TEST( sizeof(Matrix<uint16, 2, 2>) == 8,
        "sizeof(Matrix<uint16, 2, 2>) is " << sizeof(Matrix<uint16, 2, 2>) << " instead of 8\n"
    );
    TEST( sizeof(Matrix<byte, 3, 3>) == 9,
        "sizeof(Matrix<byte, 3, 3>) is " << sizeof(Matrix<byte, 3, 3>) << " instead of 9\n"
    );
    TEST( sizeof(Matrix<float, 4, 4>) == 64,
        "sizeof(Matrix<float, 4, 4>) is " << sizeof(Matrix<float, 4, 4>) << " instead of 64\n"
    );
    TEST( sizeof(Matrix<float, 10, 10>) == 400,
        "sizeof(Matrix<float, 10, 10>) is " << sizeof(Matrix<float, 10, 10>) << " instead of 400\n"
    );

    // ********************************************************************* //
    // When included the following line must cause a static assertion.
    // Matrix<int, 1, 1> errorMatrix;

    // ********************************************************************* //
    // Test type conversion constructors
    {
        Vec2 v0(0);  // vector with zeros
//      Vec2 v1 = 1; // This line should not compile (no accidentally conversion)
        Vec2 v2(-1); // Test a scalar of a different type (should yield the standard warning)
        Vec2 v3 = static_cast<Vec2>(-1); // Test a scalar of a different type without warning
        IVec2 v4(v3); // Convert to other vector
//      v4 = v0;      // This line should not compile (no accidentally conversion)
    }

    // ********************************************************************* //
    // Test of access operators and element constructors (2-4 elements)
    {
        Matrix<int, 2, 1> m0(1, 2);
        Matrix<int, 1, 3> m1(1, 2, 3);
        Matrix<int, 2, 2> m2(1, 2, 3, 4);
        m1(0,2) = 5;
        m2[2] = 5;
        TEST( m0[0] == 1 && m0[1] == 2, "Operator or constructor wrong!\n" );
        TEST( m1(0,0) == 1 && m1(0,1) == 2 && m1[2] == 5, "Operator or constructor wrong!\n" );
        TEST( m2(0,1) == 2 && m2(1,0) == 5 && m2(1,1) == 4, "Operator or constructor wrong!\n" );
        TEST( (const_cast<const Matrix<int, 2, 2>&>(m2))[1] == 2, "Read-only operator wrong!\n" );
        TEST( (const_cast<const Matrix<int, 2, 2>&>(m2))(1,0) == 5, "Read-only operator wrong!\n" );
    }

    // ********************************************************************* //
    // Test + and -
    {
        Matrix<int, 2, 2> m0(1, 2, 3, 4);
        Matrix<int, 2, 2> m1(5, 6, 7, 8);
        Matrix<int, 2, 2> m2(6, 8, 10, 12);
        Matrix<float, 2, 2> m3(3.0f, 1.0f, 4.0f, 1.0f);
        Matrix<float, 2, 2> m4(4.0f, 3.0f, 7.0f, 5.0f);
        TEST( m0 + m1 == m2, "Component wise addition failed!\n" );
        TEST( m0 + m3 == m4, "Component wise addition failed!\n" );
        TEST( m2 - m1 == m0, "Component wise subtraction failed!\n" );
        TEST( m0 - m4 == -m3, "Component wise subtraction or unary minus failed!\n" );
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
        TEST( v1 * v2 == v3, "Component wise multiplication failed!\n" );
        TEST( v0 * v4 == v5, "Component wise multiplication failed!\n" );
        TEST( v1 / v6 == v7, "Component wise division failed!\n" );
        TEST( v0 == v5 / v4, "Component wise division failed!\n" );
        TEST( v1 * v0 == -2.0f, "Row times column vector should be a scalar!\n" );
        TEST( v0 * v1 == m0, "Column times row vector should be a matrix!\n" );
        TEST( v4 * v2 == m1, "Column times row vector should be a matrix!\n" );
        TEST( m0 * m1 == m2, "Matrix multiplication wrong!\n" );
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
    }

    // ********************************************************************* //
    // Test max, min, sum, avg, abs, sgn, sign
    {
        Vec4 v0(-1.0f, -0.5f, 0.25f, 0.5f);
        Vec4 v1(0.5f, 2.5f, 1.0f, 0.0f);
        Vec4 v2(-1.0f, -0.5f, 0.25f, 0.0f);
        Vec4 v3(0.5f, 2.5f, 1.0f, 0.5f);
        Vec4 v4(1.0f, 0.5f, 0.25f, 0.5f);
        TEST( min(v0, v1) == v2, "Component wise minimum wrong!" );
        TEST( max(v0, v1) == v3, "Component wise maximum wrong!" );
        TEST( min(v0) == -1.0f, "Minimum of a vector wrong!" );
        TEST( min(v1) == 0.0f, "Minimum of a vector wrong!" );
        TEST( max(v2) == 0.25f, "Maximum of a vector wrong!" );
        TEST( max(v3) == 2.5f, "Maximum of a vector wrong!" );
        TEST( sum(v1) == 4.0f, "Sum of a vector wrong!" );
        TEST( avg(v3) == 1.125f, "Average of a vector wrong!" );
        TEST( abs(v0) == v4, "Absolute value of a vector wrong!" );
        TEST( sign(v1) == Vec4(1.0f, 1.0f, 1.0f, 0.0f), "Sign values of a vector wrong!" );
        TEST( sgn(v2) == Vec4(-1.0f, -1.0f, 1.0f, 1.0f), "Sgn values of a vector wrong!" );
    }

    // ********************************************************************* //
    // Test lerp, bilerp
    {
        Vec3 v0(0.0f);
        Vec3 v1(1.0f);
        Vec3 v2(2.0f);
        TEST( lerp(v0, v2, 0.5f) == v1, "Linear interpolation failed." );
        TEST( bilerp(v0, v1, v1, v2, 0.5f, 0.5f) == v1, "Bilinear interpolation failed." );
        TEST( lerp(1, 4, 0.5f) == 2.5f, "Linear interpolation failed." );
        TEST( bilerp(1, 3, 5, 7, 0.5f, 0.5f) == 4.0f, "Bilinear interpolation failed." );
    }

    return result;
}
