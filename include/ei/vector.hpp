#pragma once

#include <type_traits>

#include "elementarytypes.hpp"
#include "details/matrixcomponents.hpp"

namespace ei {

    template<typename T>
    class TQuaternion;

    /// \brief MxN row-major matrix class.
    /// \details The matrix is the basis for all matrix and vector types. It
    ///     supports all kinds of matrix <-> matrix and matrix <-> scalar
    ///     operations including adding, ... a value to all components. If not
    ///     stated else each operator works component wise.
    ///
    ///     There is a list of creating functions to build transformation
    ///     matrices: translation(), rotation(), scaling(), ...
    ///     All those functions create matrices for column vectors. Hence, a
    ///     transformation is done by multiplying the vectors from right:
    ///     rotation() * v. Thus in translation() * rotation() * v the
    ///     rotation is applied first and the translation afterwards.
    /// \tparam M Number of rows.
    /// \tparam N Number of columns.
    template<typename T, uint M, uint N>
    class Matrix: public details::Components<T, M, N>
    {
    public:
        static_assert(M * N > 1, "A matrix must have at least 1x2 or 2x1 components.");

        // There are two types of methods:
        // matrix <-> matrix: template<N0,M0,T0,N1,M1,T1> and
        // matrix <-> scalar: template<N0,M0,T0,T1>
        // Now, the problem with template type deduction is that the second
        // variant is used for matrices too (T1 = <N,M,T>). The following macro
        // forbids this wrong deduction. The enable if avoids operator-
        // instantiations with of Data2 = Matrix<X>.
        // declval is a standard conform way to find the resulting type of an
        // operation without the need of a constructor. This type deduction
        // construct inherits rules as [int + float -> float] from the
        // elementary types.
#       define RESULT_TYPE(op) typename std::enable_if<                 \
            !std::is_base_of<details::NonScalarType, T1>::value &&      \
            !std::is_base_of<details::NonScalarType, T>::value,         \
            decltype(std::declval<T>() op std::declval<T1>())           \
        >::type

        // Enable a function on a condition via template list.
        // This macro allows conditional compilation even for methods without
        // parameters and return value.
        // Therefore if must be inserted in the template list and `class` at
        // the same position in the implementation.
#       define ENABLE_IF(condition) class = typename std::enable_if<(condition), class Dummy>::type

        // The new variadic templates allow a more generic definition of all
        // those constructors but compiler support is lagging.

        /// \brief Construction without initialization. The values are undefined!
        Matrix();

        /// \brief Set all values to the same constant value.
        explicit Matrix(T _s);                                                 // TESTED

        /// \brief Convert a matrix/vector with a different elementary type.
        template<typename T1>
        explicit Matrix(const Matrix<T1,M,N>& _mat1);                          // TESTED

        /// \brief Allow explicit truncation of the dimension sizes.
        template<typename T1, uint M1, uint N1, ENABLE_IF((M < M1 && N <= N1) || (M <= M1 && N < N1))>
        explicit Matrix(const Matrix<T1,M1,N1>& _mat1, uint _rowOff = 0, uint _colOff = 0);

        /// \brief Construction from N * M scalar values (up to 16 elements).
        /// \details The template meta programming trick allows only the
        ///    compilation of the matching constructor.
        template<ENABLE_IF(N * M == 2)>
        Matrix(T _s0, T _s1);                                                  // TESTED
        template<ENABLE_IF(N * M == 3)>
        Matrix(T _s0, T _s1, T _s2);                                           // TESTED
        template<ENABLE_IF(N * M == 4)>
        Matrix(T _s0, T _s1, T _s2, T _s3);                                    // TESTED
        template<ENABLE_IF(N * M == 5)>
        Matrix(T _s0, T _s1, T _s2, T _s3, T _s4);                             // TESTED
        template<ENABLE_IF(N * M == 6)>
        Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5);                      // TESTED
        template<ENABLE_IF(N * M == 8)>
        Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7);        // TESTED
        template<ENABLE_IF(N * M == 9)>
        Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7, T _s8); // TESTED
        template<ENABLE_IF(N * M == 12)>
        Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7, T _s8, T _s9, T _s10, T _s11); // TESTED
        template<ENABLE_IF(N * M == 16)>
        Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7, T _s8, T _s9, T _s10, T _s11, T _s12, T _s13, T _s14, T _s15); // TESTED

        /// \brief Construction from mixed vectors and scalars for larger
        ///		vectors (column vectors)
        template<typename T1, typename T2, ENABLE_IF(M == 3 && N == 1)>
        Matrix(T1 _s0, Matrix<T2,2,1> _v12);                                   // TESTED
        template<typename T1, typename T2, ENABLE_IF(M == 3 && N == 1)>
        Matrix(Matrix<T1,2,1> _v01, T2 _s2);                                   // TESTED
        template<typename T1, typename T2, ENABLE_IF(M == 4 && N == 1)>
        Matrix(Matrix<T1,2,1> _v01, Matrix<T2,2,1> _v23);                      // TESTED
        template<typename T1, typename T2, typename T3, ENABLE_IF(M == 4 && N == 1)>
        Matrix(T1 _s0, T2 _s1, Matrix<T3,2,1> _v23);                           // TESTED
        template<typename T1, typename T2, typename T3, ENABLE_IF(M == 4 && N == 1)>
        Matrix(T1 _s0, Matrix<T2,2,1> _v12, T3 _s3);                           // TESTED
        template<typename T1, typename T2, typename T3, ENABLE_IF(M == 4 && N == 1)>
        Matrix(Matrix<T1,2,1> _v01, T2 _s2, T3 _s3);                           // TESTED
        template<typename T1, typename T2, ENABLE_IF(M == 4 && N == 1)>
        Matrix(T1 _s0, Matrix<T2,3,1> _v123);                                  // TESTED
        template<typename T1, typename T2, ENABLE_IF(M == 4 && N == 1)>
        Matrix(Matrix<T1,3,1> _v012, T2 _s3);                                  // TESTED
        /// \brief Construction from mixed vectors and scalars for larger
        ///		vectors (row vectors)
        template<typename T1, typename T2, ENABLE_IF(M == 1 && N == 3)>
        Matrix(T1 _s0, Matrix<T2,1,2> _v12);                                   // TESTED
        template<typename T1, typename T2, ENABLE_IF(M == 1 && N == 3)>
        Matrix(Matrix<T1,1,2> _v01, T2 _s2);                                   // TESTED
        template<typename T1, typename T2, ENABLE_IF(M == 1 && N == 4)>
        Matrix(Matrix<T1,1,2> _v01, Matrix<T2,1,2> _v23);                      // TESTED
        template<typename T1, typename T2, typename T3, ENABLE_IF(M == 1 && N == 4)>
        Matrix(T1 _s0, T2 _s1, Matrix<T3,1,2> _v23);                           // TESTED
        template<typename T1, typename T2, typename T3, ENABLE_IF(M == 1 && N == 4)>
        Matrix(T1 _s0, Matrix<T2,1,2> _v12, T3 _s3);                           // TESTED
        template<typename T1, typename T2, typename T3, ENABLE_IF(M == 1 && N == 4)>
        Matrix(Matrix<T1,1,2> _v01, T2 _s2, T3 _s3);                           // TESTED
        template<typename T1, typename T2, ENABLE_IF(M == 1 && N == 4)>
        Matrix(T1 _s0, Matrix<T2,1,3> _v123);                                  // TESTED
        template<typename T1, typename T2, ENABLE_IF(M == 1 && N == 4)>
        Matrix(Matrix<T1,1,3> _v012, T2 _s3);                                  // TESTED

        template<ENABLE_IF(M == 3 && N == 3)>
        explicit Matrix(const TQuaternion<T>& _quaternion);

        /// \brief Access a single element with two indices.
        /// \details Computes the data index _row * N + _col. Therefore
        ///    iterating over _col in inner loops is fastest.
        /// \param [in] _row Zero-based index of the row. For row vectors this
        ///    must be zero.
        /// \param [in] _col Zero-based index of the column. For column vectors
        ///    this must be zero.
        /// \returns Reference with read or write access to the element
        ///    depending on the constness of the matrix.
        T& operator() (uint _row, uint _col);                                  // TESTED
        T operator() (uint _row, uint _col) const;                             // TESTED

        /// \brief Single row access
        /// \param [in] _row Index of the row in [0,M-1].
        Matrix<T,1,N>& operator() (uint _row);                                 // TESTED
        const Matrix<T,1,N>& operator() (uint _row) const;                     // TESTED

        /// \brief Access an element by a single index treating the matrix as 1D.
        /// \param [in] _index Index in the range [0, N * M - 1].
        /// \returns Reference with read or write access to the element
        ///    depending on the constness of the matrix.
        T& operator[] (uint _index);                                           // TESTED
        T operator[] (uint _index) const;                                      // TESTED

        /// \brief Add two matrices component wise.
        /// \details Addition is commutative.
        template<typename T1>
        Matrix<RESULT_TYPE(+), M, N> operator+ (const Matrix<T1,M,N>& _mat1) const;  // TESTED
        /// \brief Subtract two matrices component wise.
        /// \details Subtraction is not commutative.
        template<typename T1>
        Matrix<RESULT_TYPE(-), M, N> operator- (const Matrix<T1,M,N>& _mat1) const;  // TESTED

        /// \brief Unary minus on all components.
        Matrix<T, M, N> operator- () const;                                    // TESTED

        /// \brief Matrix multiplication.
        /// \details Matrix multiplication is not commutative.
        /// \returns Matrix product with dimensions MxO = MxN * NxO. The result
        ///    is a scalar if M = N = 1.
        template<typename T1, uint O>
        typename std::conditional<M * O == 1, RESULT_TYPE(*), Matrix<RESULT_TYPE(*), M, O>>::type
        operator* (const Matrix<T1,N,O>& _mat1) const;                               // TESTED

        /// \brief Component wise binary or.
        /// \details Or is commutative.
        template<typename T1>
        Matrix<RESULT_TYPE(|), M, N> operator| (const Matrix<T1,M,N>& _mat1) const; // TESTED
        /// \brief Component wise binary and.
        /// \details And is commutative.
        template<typename T1>
        Matrix<RESULT_TYPE(&), M, N> operator& (const Matrix<T1,M,N>& _mat1) const; // TESTED
        /// \brief Component wise binary xor.
        /// \details Xor is commutative.
        template<typename T1>
        Matrix<RESULT_TYPE(^), M, N> operator^ (const Matrix<T1,M,N>& _mat1) const; // TESTED
        /// \brief Component wise binary not.
        Matrix<T, M, N> operator~ () const;                                         // TESTED
        /// \brief Component wise modulo (rest of integer division).
        template<typename T1>
        Matrix<RESULT_TYPE(%), M, N> operator% (const Matrix<T1,M,N>& _mat1) const; // TESTED

        /// \brief Self assigning component wise addition.
        template<typename T1>
        Matrix<T, M, N>& operator+= (const Matrix<T1,M,N>& _mat1);             // TESTED
        /// \brief Self assigning component wise subtraction.
        template<typename T1>
        Matrix<T, M, N>& operator-= (const Matrix<T1,M,N>& _mat1);             // TESTED
        /// \brief Self assigning component wise multiplication for vectors
        ///    of the same size and squared matrices.
        template<typename T1, ENABLE_IF(N == 1)>
        Matrix<T, M, 1>& operator*= (const Matrix<T1,M,1>& _mat1);             // TESTED
        template<typename T1, ENABLE_IF(M == 1)>
        Matrix<T, 1, N>& operator*= (const Matrix<T1,1,N>& _mat1);             // TESTED
        template<typename T1, ENABLE_IF(M == N)>
        Matrix<T, M, N>& operator*= (const Matrix<T1,N,M>& _mat1);
        /// \brief Self assigning component wise division for vectors
        ///    of the same size.
        template<typename T1, ENABLE_IF(N == 1)>
        Matrix<T, M, 1>& operator/= (const Matrix<T1,M,1>& _mat1);             // TESTED
        template<typename T1, ENABLE_IF(M == 1)>
        Matrix<T, 1, N>& operator/= (const Matrix<T1,1,N>& _mat1);             // TESTED
        /// \brief Self assigning component wise binary or.
        template<typename T1>
        Matrix<T, M, N>& operator|= (const Matrix<T1,M,N>& _mat1);             // TESTED
        /// \brief Self assigning component wise binary and.
        template<typename T1>
        Matrix<T, M, N>& operator&= (const Matrix<T1,M,N>& _mat1);             // TESTED
        /// \brief Self assigning component wise binary xor.
        template<typename T1>
        Matrix<T, M, N>& operator^= (const Matrix<T1,M,N>& _mat1);             // TESTED
        /// \brief Self assigning modulo (rest of integer division)
        template<typename T1>
        Matrix<T, M, N>& operator%= (const Matrix<T1,M,N>& _mat1);             // TESTED

        /// \brief Self assigning scalar addition.
        template<typename T1>
        Matrix<T, M, N>& operator+= (T1 _s);                                   // TESTED
        /// \brief Self assigning scalar subtraction.
        template<typename T1>
        Matrix<T, M, N>& operator-= (T1 _s);                                   // TESTED
        /// \brief Self assigning scalar multiplication.
        template<typename T1>
        Matrix<T, M, N>& operator*= (T1 _s);                                   // TESTED
        /// \brief Self assigning scalar division.
        template<typename T1>
        Matrix<T, M, N>& operator/= (T1 _s);                                   // TESTED
        /// \brief Self assigning scalar or.
        template<typename T1>
        Matrix<T, M, N>& operator|= (T1 _s);                                   // TESTED
        /// \brief Self assigning scalar and.
        template<typename T1>
        Matrix<T, M, N>& operator&= (T1 _s);                                   // TESTED
        /// \brief Self assigning scalar xor.
        template<typename T1>
        Matrix<T, M, N>& operator^= (T1 _s);                                   // TESTED
        /// \brief Self assigning scalar modulo (rest of integer division).
        template<typename T1>
        Matrix<T, M, N>& operator%= (T1 _s);                                   // TESTED
        /// \brief Self assigning component wise shift
        template<typename T1>
        Matrix<T, M, N>& operator>>= (T1 _s);                                  // TESTED
        template<typename T1>
        Matrix<T, M, N>& operator<<= (T1 _s);                                  // TESTED

        /// \brief Compare component wise, if two matrices are identical.
        Matrix<bool,M,N> operator== (const Matrix<T,M,N>& _mat1) const;        // TESTED
        /// \brief Compare component wise, if two matrices are distinct.
        Matrix<bool,M,N> operator!= (const Matrix<T,M,N>& _mat1) const;        // TESTED
        /// \brief Compare component wise, if elements are smaller or equal.
        Matrix<bool,M,N> operator<= (const Matrix<T,M,N>& _mat1) const;        // TESTED
        /// \brief Compare component wise, if elements are smaller.
        Matrix<bool,M,N> operator< (const Matrix<T,M,N>& _mat1) const;         // TESTED
        /// \brief Compare component wise, if elements are greater.
        Matrix<bool,M,N> operator> (const Matrix<T,M,N>& _mat1) const;         // TESTED
        /// \brief Compare component wise, if elements are greater or equal.
        Matrix<bool,M,N> operator>= (const Matrix<T,M,N>& _mat1) const;        // TESTED
    };



    // ********************************************************************* //
    // Vector operators

    /// \brief Component wise multiplication for vectors of the same size.
    template<typename T, typename T1, uint M>
    Matrix<RESULT_TYPE(*), M, 1> operator* (const Matrix<T,M,1>& _mat0, const Matrix<T1,M,1>& _mat1);  // TESTED
    template<typename T, typename T1, uint N>
    Matrix<RESULT_TYPE(*), 1, N> operator* (const Matrix<T,1,N>& _mat0, const Matrix<T1,1,N>& _mat1);  // TESTED
    /// \brief Component wise division for vectors of the same size.
    template<typename T, typename T1, uint M>
    Matrix<RESULT_TYPE(/), M, 1> operator/ (const Matrix<T,M,1>& _mat0, const Matrix<T1,M,1>& _mat1);  // TESTED
    template<typename T, typename T1, uint N>
    Matrix<RESULT_TYPE(/), 1, N> operator/ (const Matrix<T,1,N>& _mat0, const Matrix<T1,1,N>& _mat1);  // TESTED

    // ********************************************************************* //
    // Scalar operators

    /// \brief Add a scalar to all components.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(+), M, N> operator+ (const Matrix<T,M,N>& _mat, T1 _s); // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(+), M, N> operator+ (T1 _s, const Matrix<T,M,N>& _mat); // TESTED
    /// \brief Subtract a scalar from all components.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(-), M, N> operator- (const Matrix<T,M,N>& _mat, T1 _s); // TESTED
    /// \brief Subtract all components from a scalar.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(-), M, N> operator- (T1 _s, const Matrix<T,M,N>& _mat); // TESTED
    /// \brief Multiply a scalar to all components.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(*), M, N> operator* (const Matrix<T,M,N>& _mat, T1 _s); // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(*), M, N> operator* (T1 _s, const Matrix<T,M,N>& _mat); // TESTED
    /// \brief Divide all components by a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(/), M, N> operator/ (const Matrix<T,M,N>& _mat, T1 _s); // TESTED
    /// \brief Divide a scalar by each component.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(/), M, N> operator/ (T1 _s, const Matrix<T,M,N>& _mat); // TESTED

    /// \brief Binary or of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(|), M, N> operator| (const Matrix<T,M,N>& _mat, T1 _s); // TESTED
    /// \brief Binary or of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(|), M, N> operator| (T1 _s, const Matrix<T,M,N>& _mat); // TESTED
    /// \brief Binary and of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(&), M, N> operator& (const Matrix<T,M,N>& _mat, T1 _s); // TESTED
    /// \brief Binary and of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(&), M, N> operator& (T1 _s, const Matrix<T,M,N>& _mat); // TESTED
    /// \brief Binary xor of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(^), M, N> operator^ (const Matrix<T,M,N>& _mat, T1 _s); // TESTED
    /// \brief Binary xor of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(^), M, N> operator^ (T1 _s, const Matrix<T,M,N>& _mat); // TESTED
    /// \brief Modulo of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(%), M, N> operator% (const Matrix<T,M,N>& _mat, T1 _s); // TESTED
    /// \brief Modulo of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(%), M, N> operator% (T1 _s, const Matrix<T,M,N>& _mat); // TESTED
    /// \brief Component wise shift.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(>>), M, N> operator>> (const Matrix<T,M,N>& _mat, T1 _s); // TESTED
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(<<), M, N> operator<< (const Matrix<T,M,N>& _mat, T1 _s); // TESTED

    /// \brief Test all components with respect to a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> operator== (const Matrix<T,M,N>& _mat, T1 _s);          // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> operator== (T1 _s, const Matrix<T,M,N>& _mat);          // TESTED
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> operator!= (const Matrix<T,M,N>& _mat, T1 _s);          // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> operator!= (T1 _s, const Matrix<T,M,N>& _mat);          // TESTED
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> operator< (const Matrix<T,M,N>& _mat, T1 _s);           // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> operator< (T1 _s, const Matrix<T,M,N>& _mat);           // TESTED
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> operator<= (const Matrix<T,M,N>& _mat, T1 _s);          // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> operator<= (T1 _s, const Matrix<T,M,N>& _mat);          // TESTED
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> operator>= (const Matrix<T,M,N>& _mat, T1 _s);          // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> operator>= (T1 _s, const Matrix<T,M,N>& _mat);          // TESTED
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> operator> (const Matrix<T,M,N>& _mat, T1 _s);           // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> operator> (T1 _s, const Matrix<T,M,N>& _mat);           // TESTED

    // ********************************************************************* //
    // Predefined float vector and matrix types.
    /// \brief General column-vector type
    template<typename T, unsigned N> using Vec  = Matrix<T, N, 1>;
    /// \brief General row-vector type
    template<typename T, unsigned N> using RVec = Matrix<T, 1, N>;

    /// \brief 2D column-vector of type float.
    typedef Matrix<float, 2, 1> Vec2;
    /// \brief 3D column-vector of type float.
    typedef Matrix<float, 3, 1> Vec3;
    /// \brief 4D column-vector of type float.
    typedef Matrix<float, 4, 1> Vec4;

    /// \brief 2D row-vector of type float.
    typedef Matrix<float, 1, 2> RVec2;
    /// \brief 3D row-vector of type float.
    typedef Matrix<float, 1, 3> RVec3;
    /// \brief 4D row-vector of type float.
    typedef Matrix<float, 1, 4> RVec4;

    /// \brief 2x2 matrix of type float.
    typedef Matrix<float, 2, 2> Mat2x2;
    /// \brief 3x3 matrix of type float.
    typedef Matrix<float, 3, 3> Mat3x3;
    /// \brief 4x4 matrix of type float.
    typedef Matrix<float, 4, 4> Mat4x4;

    // ********************************************************************* //
    // Predefined double vector and matrix types.

    /// \brief 2D column-vector of type double.
    typedef Matrix<double, 2, 1> DVec2;
    /// \brief 3D column-vector of type double.
    typedef Matrix<double, 3, 1> DVec3;
    /// \brief 4D column-vector of type double.
    typedef Matrix<double, 4, 1> DVec4;

    /// \brief 2D row-vector of type double.
    typedef Matrix<double, 1, 2> DRVec2;
    /// \brief 3D row-vector of type double.
    typedef Matrix<double, 1, 3> DRVec3;
    /// \brief 4D row-vector of type double.
    typedef Matrix<double, 1, 4> DRVec4;

    /// \brief 2x2 matrix of type double.
    typedef Matrix<double, 2, 2> DMat2x2;
    /// \brief 3x3 matrix of type double.
    typedef Matrix<double, 3, 3> DMat3x3;
    /// \brief 4x4 matrix of type double.
    typedef Matrix<double, 4, 4> DMat4x4;

    // ********************************************************************* //
    // Predefined 32 bit integer vector and matrix types.

    /// \brief 2D column-vector of type int32.
    typedef Matrix<int32, 2, 1> IVec2;
    /// \brief 3D column-vector of type int32.
    typedef Matrix<int32, 3, 1> IVec3;
    /// \brief 4D column-vector of type int32.
    typedef Matrix<int32, 4, 1> IVec4;

    /// \brief 2D row-vector of type int32.
    typedef Matrix<int32, 1, 2> IRVec2;
    /// \brief 3D row-vector of type int32.
    typedef Matrix<int32, 1, 3> IRVec3;
    /// \brief 4D row-vector of type int32.
    typedef Matrix<int32, 1, 4> IRVec4;

    /// \brief 2x2 matrix of type int32.
    typedef Matrix<int32, 2, 2> IMat2x2;
    /// \brief 3x3 matrix of type int32.
    typedef Matrix<int32, 3, 3> IMat3x3;
    /// \brief 4x4 matrix of type int32.
    typedef Matrix<int32, 4, 4> IMat4x4;

    // ********************************************************************* //
    // Predefined 32 bit unsigned integer vector and matrix types.

    /// \brief 2D column-vector of type uint32.
    typedef Matrix<uint32, 2, 1> UVec2;
    /// \brief 3D column-vector of type uint32.
    typedef Matrix<uint32, 3, 1> UVec3;
    /// \brief 4D column-vector of type uint32.
    typedef Matrix<uint32, 4, 1> UVec4;

    /// \brief 2D row-vector of type uint32.
    typedef Matrix<uint32, 1, 2> URVec2;
    /// \brief 3D row-vector of type uint32.
    typedef Matrix<uint32, 1, 3> URVec3;
    /// \brief 4D row-vector of type uint32.
    typedef Matrix<uint32, 1, 4> URVec4;

    /// \brief 2x2 matrix of type uint32.
    typedef Matrix<uint32, 2, 2> UMat2x2;
    /// \brief 3x3 matrix of type uint32.
    typedef Matrix<uint32, 3, 3> UMat3x3;
    /// \brief 4x4 matrix of type uint32.
    typedef Matrix<uint32, 4, 4> UMat4x4;



    // ********************************************************************* //
    //                            QUATERNION TYPE                            //
    // ********************************************************************* //
    
    // ********************************************************************* //
    /// \brief 4D complex number equivalent for the representation of rotations
    /// \details The normalized form has len(q) == 1 and r>0 (RHS) / r<0 (LHS).
    ///     The second criteria makes the rotation unique because q and -q both
    ///     represent the same rotation and allows including mirroring.
    ///     This is not the standard way: Usual quaternions cannot handle
    ///     mirroring!
    template<typename T>
    class TQuaternion: public details::NonScalarType
    {
    public:
        /// \brief Construct uninitialized
        TQuaternion() {}

        /// \brief Copy construction
        TQuaternion( const TQuaternion& _other ) = default;
        /// \brief Copying assignment
        TQuaternion& operator = ( const TQuaternion& _rhs ) = default;

        /// \brief Construct from normalized axis and angle
        TQuaternion( const Vec<T,3>& _axis, T _angle );                        // TESTED

        /// \brief Create from Euler angles
        /// \details The rotations are applied in the order x, y, z:
        ///     rotationZ(_z) * rotationY(_y) * rotationX(_x)
        TQuaternion( T _x, T _y, T _z );                                       // TESTED
        TQuaternion( const Vec<T,3>& _eulerAngles ) : TQuaternion(_eulerAngles.x, _eulerAngles.y, _eulerAngles.z) {}

        /// \brief Create from rotation matrix (does a decomposition if the
        ///     matrix contains scaling).
        TQuaternion( const Matrix<T,3,3>& _matrix );                           // TESTED

        /// \brief Create from orthogonal basis vectors.
        TQuaternion( const Vec<T,3>& _xAxis, const Vec<T,3>& _yAxis, const Vec<T,3>& _zAxis );

        /// \brief Create from TQuaternion coefficients
        TQuaternion( T _i, T _j, T _k, T _r );

        /// \brief Rotate from vector to vector (rotated such that the from
        ///     vector is aligned with the to vector).
        /// \param [in] _from One certain direction vector before rotation.
        /// \param [in] _to Target direction vector. The from direction should
        ///     be aligned with the target after rotation.
        TQuaternion( const Vec<T,3>& _from, const Vec<T,3>& _to );             // TESTED

        // TODO: lookAt parametrization

        /// \brief Compare component wise, if two quaternions are identical.
        bool operator== (const TQuaternion& _q1) const;
        /// \brief Compare component wise, if two quaternions are different.
        bool operator!= (const TQuaternion& _q1) const;
        /// \brief TQuaternion multiplication is a combination of rotations.
        /// \details Non commutative (a*b != a*b)
        TQuaternion& operator*= (const TQuaternion& _q1);
        /// \brief Scale the TQuaternion
        TQuaternion& operator*= (T _s);
        /// \brief TQuaternion division   a/=b  <=>  a=a*(b^-1)=a*conjugated(b).
        TQuaternion& operator/= (const TQuaternion& _q1);
        /// \brief Scale the TQuaternion
        TQuaternion& operator/= (T _s);
        /// \brief Vector like addition
        TQuaternion& operator+= (const TQuaternion& _q1);
        /// \brief Vector like subtraction
        TQuaternion& operator-= (const TQuaternion& _q1);

        TQuaternion operator* (TQuaternion _q1) const { return _q1 *= *this; }
        TQuaternion operator* (T _s) const            { return TQuaternion(*this) *= _s; }
        TQuaternion operator/ (TQuaternion _q1) const { return _q1 /= *this; }
        TQuaternion operator/ (T _s) const            { return TQuaternion(*this) /= _s; }
        TQuaternion operator+ (TQuaternion _q1) const { return _q1 += *this; }
        TQuaternion operator- (TQuaternion _q1) const { return _q1 -= *this; }

        /// \brief Negate all components, the represented rotation is the same
        TQuaternion operator- () const;                                        // TESTED
        /// \brief Conjugate the quaternion
        TQuaternion operator~ () const;                                        // TESTED

        union {
            struct {T i, j, k, r;};         ///< Elements of 4D complex number
            T z[4];                         ///< Array access, index 3 is the real part
        };
    };

    // ********************************************************************* //
    /// \brief Returns identity element of the Hamilton-product. (Does not
    ///     rotate anything.)
    const TQuaternion<float>& qidentity();                                         // TESTED
    const TQuaternion<double>& qidentityD();

    // ********************************************************************* //
    /// \brief Scalar multiplication from left
    template<typename T>
    TQuaternion<T> operator* (T _s, TQuaternion<T> _q) {return _q *= _s; }

    // ********************************************************************* //
    /// \brief Complex conjugate: invert sign of complex components
    template<typename T>
    TQuaternion<T> conjugate(const TQuaternion<T>& _q);                        // TESTED

    /// \brief Get the rotation axis from a TQuaternion
    template<typename T>
    Vec<T,3> axis(const TQuaternion<T>& _q);                                   // TESTED

    /// \brief Get the x axis of the corresponding orthogonal system (rotation
    ///     matrix)
    template<typename T>
    Vec<T,3> xaxis(const TQuaternion<T>& _q);                                  // TESTED
    /// \brief Get the y axis of the corresponding orthogonal system (rotation
    ///     matrix)
    template<typename T>
    Vec<T,3> yaxis(const TQuaternion<T>& _q);                                  // TESTED
    /// \brief Get the z axis of the corresponding orthogonal system (rotation
    ///     matrix)
    template<typename T>
    Vec<T,3> zaxis(const TQuaternion<T>& _q);                                  // TESTED
    // TODO: row vector axis

    /// \brief Get the angle (radians) from a TQuaternion
    template<typename T>
    T angle(const TQuaternion<T>& _q);                                         // TESTED

    // ********************************************************************* //
    /// \brief Get the Euler angles (radians) from a quaternion
    template<typename T>
    Vec<T,3> angles(const TQuaternion<T>& _q);

    // ********************************************************************* //
    typedef TQuaternion<float> Quaternion;
    typedef TQuaternion<double> DQuaternion;


    // ********************************************************************* //
    //                               FUNCTIONS                               //
    // ********************************************************************* //

    // ********************************************************************* //
    /// \brief Check if the absolute difference between all elements is smaller
    ///    or equal than epsilon.
    /// \param [in] _mat0 First operand.
    /// \param [in] _mat1 Second operand.
    /// \param [in] _epsilon Maximum threshold for the difference between two
    ///    components. The default value is 1e-6.
    /// \returns true if all differences are less or equal than _epsilon.
    template<typename T, unsigned M, unsigned N>
    bool approx(const Matrix<T,M,N>& _mat0,
                const Matrix<T,M,N>& _mat1,
                T _epsilon = T(1e-6));                                        // TESTED

    // ********************************************************************* //
    /// \brief Check if the absolute difference between all elements is smaller
    ///    or equal than epsilon.
    template<typename T>
    bool approx(const TQuaternion<T>& _q0,
                const TQuaternion<T>& _q1,
                T _epsilon = T(1e-6));                                         // TESTED

    // ********************************************************************* //
    /// \brief Computes the sum of component wise products.
    /// \returns Scalar value of the sum of component products.
    template<typename T, unsigned M, unsigned N, typename T1>
    RESULT_TYPE(*) dot(const Matrix<T,M,N>& _mat0,
                       const Matrix<T1,M,N>& _mat1);                           // TESTED

    // ********************************************************************* //
    /// \brief Computes the sum of component wise products.
    /// \returns Scalar value of the sum of component products.
    float dot(const Quaternion& _q0,
              const Quaternion& _q1);

    // ********************************************************************* //
    /// \brief Computes the cross product of two 3d vectors (RHS).
    /// \returns Perpendicular vector with length |v0|·|v1|·sin(∡(v0,v1)).
    template<typename T, typename T1>
    Matrix<RESULT_TYPE(*),1,3> cross(const Matrix<T,1,3>& _v0,
                                     const Matrix<T1,1,3>& _v1);
    template<typename T, typename T1>
    Matrix<RESULT_TYPE(*),3,1> cross(const Matrix<T,3,1>& _v0,
                                     const Matrix<T1,3,1>& _v1);

    // ********************************************************************* //
    /// \brief Computes the cross product of two 2d vectors.
    /// \returns The determinant of the 2x2 matrix: v0.x·v1.y - v0.y·v1.x.
    template<typename T, typename T1>
    RESULT_TYPE(*) cross(const Matrix<T,1,2>& _v0,
                         const Matrix<T1,1,2>& _v1);
    template<typename T, typename T1>
    RESULT_TYPE(*) cross(const Matrix<T,2,1>& _v0,
                         const Matrix<T1,2,1>& _v1);

    // ********************************************************************* //
    /// \brief Computes the sum of squared components.
    /// \details This is equivalent to dot(_mat0, _mat0).
    /// \returns Squared euclidean length (scalar).
    template<typename T>
    auto lensq(const T& _elem0) -> decltype(dot(_elem0, _elem0));              // TESTED

    // ********************************************************************* //
    /// \brief Computes the root of the sum of squared components.
    /// \details This is the euclidean length for vectors and Quaternions and
    ///    the Frobenius norm for matrices.
    /// \returns Euclidean length (scalar).
    template<typename T>
    auto len(const T& _elem0) -> decltype(std::sqrt(dot(_elem0, _elem0)));     // TESTED

    // ********************************************************************* //
    /// \brief Normalizes a vector, quaternion or matrix with respect to len.
    /// \details This is equivalent to elem0 / len(_elem0).
    /// \returns Normalized vector or matrix.
    template<typename T>
    T normalize(const T& _mat0);                                               // TESTED

    // ********************************************************************* //
    /// \brief Component wise maximum.
    /// \returns A matrix with the maximum values from both inputs.
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> max(const Matrix<T,M,N>& _mat0,
                      const Matrix<T,M,N>& _mat1);                             // TESTED

    // ********************************************************************* //
    /// \brief Component wise minimum.
    /// \returns A matrix with the minimum values from both inputs.
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> min(const Matrix<T,M,N>& _mat0,
                      const Matrix<T,M,N>& _mat1);                             // TESTED

    // ********************************************************************* //
    /// \brief Maximum element from a matrix.
    /// \returns Scalar maximum value.
    template<typename T, unsigned M, unsigned N>
    T max(const Matrix<T,M,N>& _mat0);                                         // TESTED

    // ********************************************************************* //
    /// \brief Minimum element from a matrix.
    /// \returns Scalar minimum value.
    template<typename T, unsigned M, unsigned N>
    T min(const Matrix<T,M,N>& _mat0);                                         // TESTED

    // ********************************************************************* //
    /// \brief Component wise clamp to boundaries.
    /// \returns A matrix with values in the bounding box.
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> clamp(const Matrix<T,M,N>& _mat,
                        const Matrix<T,M,N>& _min,
                        const Matrix<T,M,N>& _max);                            // TESTED

    // ********************************************************************* //
    /// \brief Component wise clamp to scalar boundaries.
    /// \returns A matrix with values in the interval.
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> clamp(const Matrix<T,M,N>& _mat,
                        T _min,
                        T _max);                                               // TESTED

    /// \brief Clamp all components to [0,1]
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> saturate(const Matrix<T,M,N>& _mat)   { return clamp(_mat, static_cast<T>(0), static_cast<T>(1)); }  // TESTED

    // ********************************************************************* //
    /// \brief Round all components towards negative infinity
    template<typename T, unsigned M, unsigned N>
    Matrix<typename details::Int<sizeof(T)>::stype,M,N> floor(const Matrix<T,M,N>& _mat);   // TESTED

    // ********************************************************************* //
    /// \brief Round all components towards negative infinity
    template<typename T, unsigned M, unsigned N>
    Matrix<typename details::Int<sizeof(T)>::stype,M,N> ceil(const Matrix<T,M,N>& _mat);    // TESTED

    // ********************************************************************* //
    /// \brief Round all components towards next number (x.5 rounds up)
    template<typename T, unsigned M, unsigned N>
    Matrix<typename details::Int<sizeof(T)>::stype,M,N> round(const Matrix<T,M,N>& _mat);   // TESTED

    // ********************************************************************* //
    /// \brief Compute the square root for each component
    template<typename T, unsigned N>
    Matrix<T,1,N> sqrt(const Matrix<T,1,N>& _v0);                              // TESTED
    template<typename T, unsigned M>
    Matrix<T,M,1> sqrt(const Matrix<T,M,1>& _v0);                              // TESTED
    // The above declaration unfortunatelly hides the elementary one -> make it
    // visible again
    using std::sqrt;

    // ********************************************************************* //
    /// \brief Compute the power for each component
    template<typename T, unsigned N>
    Matrix<T,1,N> pow(const Matrix<T,1,N>& _v0, float _exponent);              // TESTED
    template<typename T, unsigned M>
    Matrix<T,M,1> pow(const Matrix<T,M,1>& _v0, float _exponent);              // TESTED
    // The above declaration unfortunatelly hides the elementary one -> make it
    // visible again
    using std::pow;

    // ********************************************************************* //
    /// \brief Element wise natural logarithm for matrices (basis e).
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> log(const Matrix<T,M,N>& _mat0);                             // TESTED
    /// \brief Element wise logarithm for matrices (basis 2).
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> log2(const Matrix<T,M,N>& _mat0);                            // TESTED
    using std::log;
    using std::log2;

    // ********************************************************************* //
    /// \brief Sum of all components.
    /// \details Can be used for boolean vectors/matrices too (number of trues).
    /// \returns Scalar sum of all values.
    template<typename T, unsigned M, unsigned N>
    decltype(std::declval<T>() + std::declval<T>()) sum(const Matrix<T,M,N>& _mat0);    // TESTED

    // ********************************************************************* //
    /// \brief Product of all components.
    /// \returns Product of all values (scalar).
    template<typename T, unsigned M, unsigned N>
    T prod(const Matrix<T,M,N>& _mat0);                                        // TESTED

    // ********************************************************************* //
    /// \brief Average of all values from a matrix.
    /// \returns Scalar average value.
    template<typename T, unsigned M, unsigned N>
    T avg(const Matrix<T,M,N>& _mat0);                                         // TESTED

    // ********************************************************************* //
    /// \brief Absolute values for all components.
    /// \returns Matrix with component wise absolute values.
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> abs(const Matrix<T,M,N>& _mat0);                             // TESTED

    // ********************************************************************* //
    /// \brief Get the component wise sign.
    /// \details There is a faster version sgn(), if you don't need to 
    ///    know about zero.
    /// \returns -1 (_x < 0), 0 (_x == 0) or 1 (_x > 0) for each component _x.
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> sign(const Matrix<T,M,N>& _mat0);                            // TESTED


    // ********************************************************************* //
    /// \brief Get the component wise sign where 0 is counted as positive.
    /// \details This function is faster than sign(). Use it if you don't need
    ///    to know about zero.
    /// \returns -1 (_x < 0) or 1 (_x >= 0) for each component _x.
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> sgn(const Matrix<T,M,N>& _mat0);                             // TESTED

    // ********************************************************************* //
    /// \brief Bilinear interpolation optimized for vectors.
    /// \param _x00 [in] Vector or matrix value. This is returned when
    ///    _t0 is zero and _t1 is zero.
    /// \param _x01 [in] Vector or matrix value. This is returned when
    ///    _t0 is one and _t1 is zero.
    /// \param _x10 [in] Vector or matrix value. This is returned when
    ///    _t0 is zero and _t1 is one.
    /// \param _x11 [in] Vector or matrix value. This is returned when
    ///    _t0 is one and _t1 is one.
    /// \param _t0 [in] Scalar interpolation parameter ("x-direction").
    /// \param _t1 [in] Scalar interpolation parameter ("y-direction").
    /// \returns lerp(lerp(_x00, _x01, _t0), lerp(_x10, _x11, _t0), _t1).
    template<typename T0, typename T1, unsigned M, unsigned N>
    Matrix<decltype(std::declval<T0>() * std::declval<T1>()),M,N>
        bilerp(Matrix<T0,M,N> _x00, Matrix<T0,M,N> _x01,
               Matrix<T0,M,N> _x10, Matrix<T0,M,N> _x11,
               T1 _t0, T1 _t1);                                                // TESTED

    // ********************************************************************* //
    /// \brief Spherical linear interpolation with constant angular speed
    /// \details Spherical interpolation is not defined, if the angle between
    ///     both vectors is pi. In that case this still method performs
    ///     an interpolation but without normalization.
    ///
    ///     Formulas from Ken Shoemake "Animating rotation with quaternion
    ///     curves" SIGGRAPH 85
    template<typename T0, typename T1, unsigned N>
    auto slerp(const Matrix<T0,1,N>& _v0, const Matrix<T0,1,N>& _v1, T1 _t) -> decltype(_v0*_t);
    template<typename T0, typename T1, unsigned M>
    auto slerp(const Matrix<T0,M,1>& _v0, const Matrix<T0,M,1>& _v1, T1 _t) -> decltype(_v0*_t);    // TESTED
    template<typename T>
    TQuaternion<T> slerp(const TQuaternion<T>& _q0, const TQuaternion<T>& _q1, T _t);  // TESTED


    // ********************************************************************* //
    /// \brief Test if at least one element of the matrix is true.
    /// \return false, if all elements off the matrix are false.
    template<unsigned M, unsigned N>
    bool any(const Matrix<bool,M,N>& _mat0);                                   // TESTED

    // ********************************************************************* //
    /// \brief Test if no element of the matrix is true.
    /// \return true, if all elements off the matrix are false.
    template<unsigned M, unsigned N>
    bool none(const Matrix<bool,M,N>& _mat0);                                  // TESTED

    // ********************************************************************* //
    /// \brief Test if all elements of the matrix are true.
    /// \return true, if all elements off the matrix are true.
    template<unsigned M, unsigned N>
    bool all(const Matrix<bool,M,N>& _mat0);                                   // TESTED


    // ********************************************************************* //
    /// \brief Transpose a matrix or vector (switch the dimensions).
    template<typename T, unsigned M, unsigned N>
    Matrix<T,N,M> transpose(const Matrix<T,M,N>& _mat0);                       // TESTED

    // ********************************************************************* //
    /// \brief Compute LUp decomposition of A such that A[p] = LU.
    /// \details L is a lower triangular matrix with zeros below the diagonal,
    ///     U is a unit upper triangular matrix with zeros above the diagonal
    ///     and all 1 on the main diagonal. Further p is a row permutation
    ///     vector.
    /// \param [in] _A The matrix to be decomposed.
    /// \param [out] _LU A matrix in the form (example for 3D):
    ///         |l00 u01 u02|           |l00   0   0|   |1 u01 u02|
    ///         |l10 l11 u12| such that |l10 l11   0| * |0   1 u12| = A[p]
    ///         |l20 l21 l22|           |l20 l21 l22|   |0   0   1|
    /// \return true if the decomposition was possible (A is non singular).
    template<typename T, unsigned N>
    bool decomposeLUp(const Matrix<T,N,N>& _A, Matrix<T,N,N>& _LU, Vec<uint,N>& _p);    // TESTED

    // ********************************************************************* //
    /// \brief Compute one or more solutions for linear equation system given
    ///     the LUp decomposition.
    /// \details Computation time: O(n²)
    /// \param [in] _LU Combined LU matrix from decomposeLUp.
    /// \param [in] _p Permutation vector from decomposeLUp.
    /// \param [in] _B Set of N column vectors which contain the target vectors.
    /// \return Solution X of LU X = B[p].
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> solveLUp(const Matrix<T,M,M>& _LU, const Matrix<uint,M,1>& _p, const Matrix<T,M,N>& _B);  // TESTED

    // ********************************************************************* //
    /// \brief Compute a spectral decomposition Q^T * diag(d) * Q of a 3x3
    ///     symmetric matrix A.
    /// \param [in] _sort Sort eigenvalues and vectors descending.
    /// \param [out] _Q The orthonormal basis where rows are the eigenvectors
    ///     corresponding to the eigenvalues _lambda of A.
    /// \param [out] _lambda Eigenvalues of A.
    /// \return Number of iterations (50 is the maximum used internally)
    template<typename T>
    int decomposeQl(const Matrix<T,3,3>& _A, Matrix<T,3,3>& _Q, Vec<T,3>& _lambda, bool _sort = false);     // TESTED

    // ********************************************************************* //
    // EXPERIMENT:
    // Rayleigh iteration is expensive per iteration and cannot guarantee that
    // the maximum λ will be found. Power iteration converges very slowly.
    // Since decomposeQl is faster then both eigenmax is discarded.
    //
    /// \brief Find the largest eigenvector and eigenvalue of a matrix.
    /// \details This uses power iteration which does not guaranteed a
    ///     convergence if eigenvalues are not distinct.
    //template<typename T, unsigned N>
    //int eigenmax(const Matrix<T,N,N>& _A, Vec<T,N>& _eigenvec, T& _eigenval);

    // ********************************************************************* //
    /// \brief Invert a quadratic matrix.
    /// \details If the matrix has no inverse the identity is returned.
    /// \return Inverse matrix or identity.
    template<typename T, unsigned N>
    Matrix<T,N,N> invert(const Matrix<T,N,N>& _mat0);                          // TESTED

    // ********************************************************************* //
    /// \brief Compute the determinant of a matrix.
    /// \details This uses fixed implementations for N=2 and N=3 and LU
    ///     decomposition for N > 3.
    template<typename T, unsigned N>
    T determinant(const Matrix<T,N,N>& _A);                                    // TESTED

    // ********************************************************************* //
    /// \brief Gram-Schmidt orthonormalization for matrices with M >= N.
    /// \details This method assumes a column space where each column is one
    ///     vector.
    ///     The method operates a stabilized Gram-Schmidt inplace.
    /// \returns false if some columns are linear dependent and not all vectors
    ///     can be orthogonalized.
    template<typename T, unsigned M, unsigned N>
    bool orthonormalize(Matrix<T,M,N>& _mat0);                                 // TESTED
    /// \brief Gram-Schmidt orthonormalization for a list of vectors.
    template<typename TVec0, typename... TVecs>
    bool orthonormalize(TVec0& _vec0, TVecs&... _vecs);                        // TESTED


    // ********************************************************************* //
    //                            TRANSFORMATIONS                            //
    // ********************************************************************* //

    // ********************************************************************* //
    /// \brief Generate the N x N identity matrix.
    template<typename T, unsigned N>
    const Matrix<T,N,N>& identity();                                           // TESTED

    // ********************************************************************* //
    /// \brief Generate the N x N diagonal matrix.
    /// \param [in] _v0 A vector with the diagonal entries.
    template<typename T, unsigned N>
    Matrix<T,N,N> diag( const Vec<T,N>& _v0 );                                 // TESTED

    // ********************************************************************* //
    /// \brief Convert a vector from Cartesian coordinates in spherical
    ///     (angular) coordinates.
    /// \details In 2d this gives (r, α) and in 3d this gives (r, θ, ϕ) for
    ///     the vectors: (r cos α, r sin α) and (r cos θ, r sin θ cos ϕ, r sin θ sin ϕ).
    ///     This is not the usual convention, but unifies dimensions!
    /// \return The N-1 spherical angles and the length of the vector.
    ///     (r, φ1, φ02, ..., φN-1) where
    ///     φ1, ..., φN-2 ∈ [0,π) and φN-1 ∈ [0,2π)
    template<typename T, unsigned N>
    Vec<T,N> sphericalCoords( const Vec<T,N>& _v0 );                           // TESTED
    template<typename T, unsigned N>
    RVec<T,N> sphericalCoords( const RVec<T,N>& _v0 );                         // TESTED

    // ********************************************************************* //
    /// \brief Convert a vector from spherical coordinates (r, φ1, φ02, ..., φN-1)
    ///     to regular Cartesian coordinates.
    /// \return The regular Cartesian vector.
    template<typename T, unsigned N>
    Vec<T,N> cartesianCoords( const Vec<T,N>& _v0 );                           // TESTED
    template<typename T, unsigned N>
    RVec<T,N> cartesianCoords( const RVec<T,N>& _v0 );                         // TESTED

    // ********************************************************************* //
    /// \brief Apply transformations in homogeneous space. This includes a
    ///     division by w after the transformation
    template<typename T, unsigned N>
    Matrix<T,N,1> transform( const Matrix<T,N,1>& _what, const Matrix<T,N+1,N+1>& _space );
    template<typename T, unsigned N>
    Matrix<T,1,N> transform( const Matrix<T,1,N>& _what, const Matrix<T,N+1,N+1>& _space );
    /// \brief Apply transformations and ignore the fourth component and the
    ///     translation.
    /// \details Use this function to transform direction vectors. Still, there
    ///     is no normalization involved and the direction might be scaled by
    ///     the matrix.
    template<typename T, unsigned N>
    Matrix<T,N,1> transformDir( const Matrix<T,N,1>& _what, const Matrix<T,N+1,N+1>& _space );
    template<typename T, unsigned N>
    Matrix<T,1,N> transformDir( const Matrix<T,1,N>& _what, const Matrix<T,N+1,N+1>& _space );
    /// \brief Apply transformations with a matrix multiplication.
    template<typename T, unsigned N>
    Matrix<T,N,1> transform( const Matrix<T,N,1>& _what, const Matrix<T,N,N>& _space );
    template<typename T, unsigned N>
    Matrix<T,1,N> transform( const Matrix<T,1,N>& _what, const Matrix<T,N,N>& _space );

    /// \brief Apply a rotation by a quaternion
    template<typename T>
    Vec<T,3> transform( const Vec<T,3>& _what, const TQuaternion<T>& _quaternion );
    template<typename T>
    RVec<T,3> transform( const RVec<T,3>& _what, const TQuaternion<T>& _quaternion );

    // ********************************************************************* //
    /// \brief Create a translation matrix in homogeneous coordinate space.
    /// \param [in] _vector Translate by/Add this vector.
    /// \details The translation matrix always has a dimension one large then
    ///    the vectors.
    ///    To transform a vector append 1 and multiply it from right:
    ///    translation() * VecX(v,1)
    template<typename T, unsigned N>
    Matrix<T,N+1,N+1> translation( const Matrix<T, N, 1>& _vector );

    // ********************************************************************* //
    /// \brief Create a scaling/diagonal matrix from vector.
    template<typename T, unsigned N>
    Matrix<T,N,N> scaling( const Matrix<T, N, 1>& _scale );

    // ********************************************************************* //
    /// \brief Create a uniform scaling/diagonal matrix from scalar.
    template<typename T, unsigned N>
    Matrix<T,N,N> scaling( T _scale );

    // ********************************************************************* //
    /// \brief Create a scaling matrix in homogeneous space.
    Matrix<float,4,4> scalingH( const Vec3& _scale );

    // ********************************************************************* //
    /// \brief Create a uniform scaling matrix in homogeneous space.
    Matrix<float,4,4> scalingH( float _scale );

    // ********************************************************************* //
    /// \brief Use vectors which span a space to build a matrix.
    /// \details The vectors become the rows of the matrix
    Mat2x2 axis( const Vec2& _x, const Vec2& _y );
    Mat3x3 axis( const Vec3& _x, const Vec3& _y, const Vec3& _z );             // TESTED
    Mat4x4 axisH( const Vec3& _x, const Vec3& _y, const Vec3& _z );
    Mat4x4 axis( const Vec4& _x, const Vec4& _y, const Vec4& _z, const Vec4& _w );

    // ********************************************************************* //
    /// \brief Create an orthonormal basis for a single direction vector
    Mat2x2 basis( const Vec2& _vector );                                       // TESTED
    Mat3x3 basis( const Vec3& _vector );                                       // TESTED

    // ********************************************************************* //
    /// \brief Rotation matrix in 2D.
    Mat2x2 rotation( float _angle );

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space around x-axis.
    Mat3x3 rotationX( float _angle );
    Mat4x4 rotationXH( float _angle );

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space around y-axis.
    Mat3x3 rotationY( float _angle );
    Mat4x4 rotationYH( float _angle );

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space around z-axis.
    Mat3x3 rotationZ( float _angle );
    Mat4x4 rotationZH( float _angle );

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space from 3 angles:
    ///     rotationZ(_z) * rotationY(_y) * rotationX(_x).
    Mat3x3 rotation( float _x, float _y, float _z );
    Mat4x4 rotationH( float _x, float _y, float _z );
    inline Mat3x3 rotation( const Vec3& _eulerAngles )         { return rotation(_eulerAngles.x, _eulerAngles.y, _eulerAngles.z); }
    inline Mat4x4 rotationH( const Vec3& _eulerAngles )        { return rotationH(_eulerAngles.x, _eulerAngles.y, _eulerAngles.z); }

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space for an arbitrary axis.
    Mat3x3 rotation( const Vec3& _axis, float _angle );
    Mat4x4 rotationH( const Vec3& _axis, float _angle );

    // ********************************************************************* //
    /// \brief Rotation matrix from one direction into another.
    Mat3x3 rotation( const Vec3& _from, const Vec3& _to );                     // TESTED
    Mat4x4 rotationH( const Vec3& _from, const Vec3& _to );

    // ********************************************************************* //
    /// \brief Rotation matrix from quaternion.
    Mat3x3 rotation( const Quaternion& _quaternion );
    Mat4x4 rotationH( const Quaternion& _quaternion );

    // ********************************************************************* //
    /// \brief Reflection matrix for a plane through the origin.
    /// \details Use this to construct mirror matrices. To simply reflect a
    ///     vector use the reflect() method (faster). Beware: _normal and _at
    ///     in the two methods are orthogonal.
    /// \param [in] _normal Normal of the reflecting plane (at the origin).
    ///     The normal must not be normalized.
    Mat3x3 housholder( const Vec3& _normal );                                  // TESTED
    Mat4x4 housholderH( const Vec3& _normal );

    // ********************************************************************* //
    /// \brief Reflect a vector at a plane.
    /// \details Reflects using I - 2 * dot(I,At) * At (householder transformation)
    /// \param [in] _incident A direction or position vector which should be
    ///     reflected.
    /// \param [in] _at The normal vector for the reflection plane (normalized!).
    template<typename T, unsigned N>
    Vec<T,N> reflect( const Vec<T,N>& _incident, const Vec<T,N>& _at );
    template<typename T, unsigned N>
    RVec<T,N> reflect( const RVec<T,N>& _incident, const RVec<T,N>& _at );

    // ********************************************************************* //
    /// \brief Create a matrix in 3D/homogeneous space where the target is on
    ///     the positive z-axis.
    /// \details This method creates an left-hand system (LHS) with positive
    ///     z-axis.
    ///
    ///     This matrix always 'looks' from the origin. Use a camera() or
    ///     combine with a translation() to add an other origin.
    /// \param [in] _target A position which should lie on the z-axis.
    /// \param [in] _up The x-axis/horizon is created perpendicular to this
    ///     vector. The up vector must not necessarily be normalized.
    Mat3x3 lookAt( const Vec3& _target, const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f));
    Mat4x4 lookAtH( const Vec3& _target, const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f));

    // ********************************************************************* //
    /// \brief Create a camera matrix in homogeneous space.
    /// \details The camera matrix must be in homogeneous space due to the
    ///    translation.
    /// \details This method creates an left-hand system (LHS) with positive
    ///    z-axis.
    Mat4x4 camera( const Vec3& _position, const Vec3& _target, const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f));

    // ********************************************************************* //
    /// \brief Create OpenGL perspective projection matrix from 6 sides.
    /// \details Assumes an LHS coordinate system.
    ///
    ///    The OpenGL frustum is defined in the [-1,-1,-1] x [1,1,1] cube.
    /// \param [in] _l Left plane x-coordinate (at near plane)
    /// \param [in] _r Right plane x-coordinate (at near plane)
    /// \param [in] _b Bottom plane y-coordinate (at near plane)
    /// \param [in] _t Bottom plane y-coordinate (at near plane)
    /// \param [in] _n Near plane
    /// \param [in] _f Far plane
    Mat4x4 perspectiveGL( float _l, float _r, float _b, float _t, float _n, float _f );

    // ********************************************************************* //
    /// \brief Create OpenGL perspective projection matrix from fovY and aspect ratio.
    /// \details Assumes an LHS coordinate system.
    ///
    ///    The OpenGL frustum is defined in the [-1,-1,-1] x [1,1,1] cube.
    /// \param [in] _fovY Field of view in the y direction, in radians.
    /// \param [in] _aspectRatio width/height of the frame buffer.
    Mat4x4 perspectiveGL( float _fovY, float _aspectRatio, float _near, float _far );

    // ********************************************************************* //
    /// \brief Create OpenGL orthographic matrix from 6 sides.
    /// \details Assumes an LHS coordinate system.
    ///
    ///    The OpenGL frustum is defined in the [-1,-1,1] x [1,1,-1] cube.
    ///    Keep in mind, that near > far because the GL-view direction is
    ///    along the negative z-axis.
    /// \param [in] _l Left plane x-coordinate (at near plane)
    /// \param [in] _r Right plane x-coordinate (at near plane)
    /// \param [in] _b Bottom plane y-coordinate (at near plane)
    /// \param [in] _t Bottom plane y-coordinate (at near plane)
    /// \param [in] _n Near plane
    /// \param [in] _f Far plane
    Mat4x4 orthographicGL( float _l, float _r, float _b, float _t, float _n, float _f );

    // ********************************************************************* //
    /// \brief Create DirectX perspective projection matrix from 6 sides.
    /// \details Assumes an LHS coordinate system.
    ///
    ///    The DirectX frustum is defined in the [-1,-1,0] x [1,1,1] cube.
    ///
    ///    This method is transposed compared to the DX documentation because
    ///    this library multiplies vectors from right.
    /// \param [in] _l Left plane x-coordinate (at near plane)
    /// \param [in] _r Right plane x-coordinate (at near plane)
    /// \param [in] _b Bottom plane y-coordinate (at near plane)
    /// \param [in] _t Bottom plane y-coordinate (at near plane)
    /// \param [in] _n Near plane
    /// \param [in] _f Far plane
    Mat4x4 perspectiveDX( float _l, float _r, float _b, float _t, float _n, float _f );

    // ********************************************************************* //
    /// \brief Create DirectX perspective projection matrix from fovY and aspect ratio.
    /// \details Assumes an LHS coordinate system.
    ///
    ///    The DirectX frustum is defined in the [-1,-1,0] x [1,1,1] cube.
    ///
    ///    This method is transposed compared to the DX documentation because
    ///    this library multiplies vectors from right.
    /// \param [in] _fovY Field of view in the y direction, in radians.
    /// \param [in] _aspectRatio width/height of the frame buffer.
    Mat4x4 perspectiveDX( float _fovY, float _aspectRatio, float _near, float _far );

    // ********************************************************************* //
    /// \brief Create DirectX orthographic projection matrix from 6 sides.
    /// \details Assumes an LHS coordinate system.
    ///
    ///    The DirectX frustum is defined in the [-1,-1,0] x [1,1,1] cube.
    ///
    ///    This method is transposed compared to the DX documentation because
    ///    this library multiplies vectors from right.
    /// \param [in] _l Left plane x-coordinate (at near plane)
    /// \param [in] _r Right plane x-coordinate (at near plane)
    /// \param [in] _b Bottom plane y-coordinate (at near plane)
    /// \param [in] _t Bottom plane y-coordinate (at near plane)
    /// \param [in] _n Near plane
    /// \param [in] _f Far plane
    Mat4x4 orthographicDX( float _l, float _r, float _b, float _t, float _n, float _f );

}

#include "details/matrixutililities.hpp"

namespace ei {
    // Include implementation.
#   include "details/matrix.inl"
#   include "details/quaternion.inl"

    /// \brief Alias for identity<float,2>().
    inline Mat2x2 identity2x2()    { return identity<float,2>(); }
    /// \brief Alias for identity<float,3>().
    inline Mat3x3 identity3x3()    { return identity<float,3>(); }             // TESTED
    /// \brief Alias for identity<float,4>().
    inline Mat4x4 identity4x4()    { return identity<float,4>(); }             // TESTED

    // Remove helper macros.
#   undef RESULT_TYPE
#   undef ENABLE_IF
}
