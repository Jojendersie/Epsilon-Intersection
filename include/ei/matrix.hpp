#pragma once

#include <type_traits>

#include "elementarytypes.hpp"
#include "details/matrixcomponents.hpp"

namespace ei {

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
#       define RESULT_TYPE(op) typename std::enable_if<					\
            !std::is_base_of<details::MatrixType, T1>::value &&			\
            !std::is_base_of<details::MatrixType, T>::value,			\
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
        explicit Matrix(const Matrix<T1,M1,N1>& _mat1);

        /// \brief Construction from N * M scalar values (up to 16 elements).
        /// \details The template meta programming trick allows only the
        ///    compilation of the matching constructor.
        template<ENABLE_IF(N * M == 2)>
        Matrix(T _s0, T _s1);                                                  // TESTED
        template<ENABLE_IF(N * M == 3)>
        Matrix(T _s0, T _s1, T _s2);                                           // TESTED
        template<ENABLE_IF(N * M == 4)>
        Matrix(T _s0, T _s1, T _s2, T _s3);                                    // TESTED
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
        /// \brief Add a scalar value to all components.

        /// \brief Matrix multiplication.
        /// \details Matrix multiplication is not commutative.
        /// \returns Matrix product with dimensions MxO = MxN * NxO. The result
        ///    is a scalar if M = N = 1.
        template<typename T1, uint O>
        typename std::conditional<M * O == 1, RESULT_TYPE(*), Matrix<RESULT_TYPE(*), M, O>>::type
        operator* (const Matrix<T1,N,O>& _mat1) const;                               // TESTED
        /// \brief Specialized version to allow simple homogeneous coordinate
        ///    transformations for column vectors. E.g. Mat4x4 * Vec3
        /// \details This transformation also divides by the last component.
      /*  template<typename T1, ENABLE_IF(N == M)>
        Matrix<RESULT_TYPE(*), N-1, 1> operator* (const Matrix<T1,N-1,1>& _mat1) const;
        /// \brief Specialized version to allow simple homogeneous coordinate
        ///    transformations for row vectors. E.g. RVec3 * Mat4x4
        /// \details This transformation also divides by the last component.
        template<typename T1, ENABLE_IF(N == M)>
        Matrix<RESULT_TYPE(*), 1, N-1> operator* (const Matrix<T1,1,N-1>& _mat1) const;*/
        /// \brief Component wise multiplication for vectors of the same size.
        template<typename T1, ENABLE_IF(N == 1)>
        Matrix<RESULT_TYPE(*), M, 1> operator* (const Matrix<T1,M,1>& _mat1) const;  // TESTED
        template<typename T1, ENABLE_IF(M == 1)>
        Matrix<RESULT_TYPE(*), 1, N> operator* (const Matrix<T1,1,N>& _mat1) const;  // TESTED
        /// \brief Component wise division for vectors of the same size.
        template<typename T1, ENABLE_IF(N == 1)>
        Matrix<RESULT_TYPE(/), M, 1> operator/ (const Matrix<T1,M,1>& _mat1) const;  // TESTED
        template<typename T1, ENABLE_IF(M == 1)>
        Matrix<RESULT_TYPE(/), 1, N> operator/ (const Matrix<T1,1,N>& _mat1) const;  // TESTED

        /// \brief Self assigning component wise addition.
        template<typename T1>
        Matrix<T, M, N>& operator+= (const Matrix<T1,M,N>& _mat1);             // TESTED
        /// \brief Self assigning component wise subtraction.
        template<typename T1>
        Matrix<T, M, N>& operator-= (const Matrix<T1,M,N>& _mat1);             // TESTED
        /// \brief Self assigning component wise multiplication for vectors
        ///    of the same size.
        template<typename T1, ENABLE_IF(N == 1)>
        Matrix<T, M, 1>& operator*= (const Matrix<T1,M,1>& _mat1);             // TESTED
        template<typename T1, ENABLE_IF(M == 1)>
        Matrix<T, 1, N>& operator*= (const Matrix<T1,1,N>& _mat1);             // TESTED
        /// \brief Self assigning component wise division for vectors
        ///    of the same size.
        template<typename T1, ENABLE_IF(N == 1)>
        Matrix<T, M, 1>& operator/= (const Matrix<T1,M,1>& _mat1);             // TESTED
        template<typename T1, ENABLE_IF(M == 1)>
        Matrix<T, 1, N>& operator/= (const Matrix<T1,1,N>& _mat1);             // TESTED

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

        /// \brief Compare component wise, if two matrices are identical.
        Matrix<bool,M,N> operator== (const Matrix<T,M,N>& _mat1) const;        // TESTED
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
   
    /// \brief Test all components with respect to a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> operator< (const Matrix<T,M,N>& _mat, T1 _s);           // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> operator< (T1 _s, const Matrix<T,M,N>& _mat);
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> operator<= (const Matrix<T,M,N>& _mat, T1 _s);          // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> operator<= (T1 _s, const Matrix<T,M,N>& _mat);
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> operator>= (const Matrix<T,M,N>& _mat, T1 _s);          // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> operator>= (T1 _s, const Matrix<T,M,N>& _mat);
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> operator> (const Matrix<T,M,N>& _mat, T1 _s);           // TESTED
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> operator> (T1 _s, const Matrix<T,M,N>& _mat);

    // ********************************************************************* //
    // Predefined float vector and matrix types.

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
    //                               FUNCTIONS                               //
    // ********************************************************************* //

    // ********************************************************************* //
    /// \brief Check if the absolute difference between all elements is smaller
    ///    or equal than epsilon.
    /// \param [in] _mat0 First operand.
    /// \param [in] _mat1 Second operand.
    /// \param [in] _epsilon Maximum threshold for the difference between two
    ///    components. The default value is 1e-6f.
    /// \returns true if all differences are less or equal than _epsilon.
    template<typename T, unsigned M, unsigned N>
    bool approx(const Matrix<T,M,N>& _mat0,
                const Matrix<T,M,N>& _mat1,
                float _epsilon = 1e-6f);

    // ********************************************************************* //
    /// \brief Computes the sum of component wise products.
    /// \returns Scalar value of the sum of component products.
    template<typename T, unsigned M, unsigned N, typename T1>
    RESULT_TYPE(*) dot(const Matrix<T,M,N>& _mat0,
                       const Matrix<T1,M,N>& _mat1);                           // TESTED

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
    template<typename T, unsigned M, unsigned N>
    T lensq(const Matrix<T,M,N>& _mat0);                                       // TESTED

    // ********************************************************************* //
    /// \brief Computes the root of the sum of squared components.
    /// \details This is the euclidean length for vectors and the Frobenius
    ///    norm for matrices.
    /// \returns Euclidean length (scalar).
    template<typename T, unsigned M, unsigned N>
    decltype(sqrt(std::declval<T>())) len(const Matrix<T,M,N>& _mat0);         // TESTED

    // ********************************************************************* //
    /// \brief Normalizes a vector or matrix with respect to len.
    /// \details This is equivalent to _mat0 / len(_mat0).
    /// \returns Normalized vector or matrix.
    template<typename T, unsigned M, unsigned N>
    Matrix<T,M,N> normalize(const Matrix<T,M,N>& _mat0);

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
    /// \brief Signed sum of all components.
    /// \returns Scalar sum of all values.
    template<typename T, unsigned M, unsigned N>
    decltype(std::declval<T>() + std::declval<T>()) sum(const Matrix<T,M,N>& _mat0);    // TESTED

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
    //                            TRANSFORMATIONS                            //
    // ********************************************************************* //

    // ********************************************************************* //
    /// \brief Generate the N x N identity matrix.
    template<typename T, unsigned N>
    Matrix<T,N,N> identity();

    // ********************************************************************* //
    /// \brief Lift vector or squared matrix to homogeneous space.
    /// \details Adds a row and a column with zeros to a matrix and sets the
    ///    new diagonal element to 1.
    ///    A     =>    A 0
    ///                0 1
    ///
    ///   Appends 1 to vectors: v    =>   (v 1)
    template<typename T, unsigned N>
    Matrix<T,N+1,N+1> homo( const Matrix<T,N,N>& _mat0 );
    template<typename T, unsigned N>
    Matrix<T,N+1,1> homo( const Matrix<T,N,1>& _v0 );
    template<typename T, unsigned N>
    Matrix<T,1,N+1> homo( const Matrix<T,1,N>& _v0 );

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
    Mat3x3 axis( const Vec3& _x, const Vec3& _y, const Vec3& _z );
    Mat4x4 axisH( const Vec3& _x, const Vec3& _y, const Vec3& _z );
    Mat4x4 axis( const Vec4& _x, const Vec4& _y, const Vec4& _z, const Vec4& _w );

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
    ///     rotationZ(_yaw) * rotationY(_pitch) * rotationX(_roll).
    // TODO: test gimbal lock
    Mat3x3 rotation( float _yaw, float _pitch, float _roll );
    Mat4x4 rotationH( float _yaw, float _pitch, float _roll );

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space for an arbitrary axis.
    Mat3x3 rotation( const Vec3& _axis, float _angle );
    Mat4x4 rotationH( const Vec3& _axis, float _angle );

    // TODO: from quaternion

    // ********************************************************************* //
    /// \brief Create a matrix in 3D/homogeneous space where the target is on
    ///    the positive z-axis.
    /// \details This method creates an left-hand system (LHS) with positive
    ///    z-axis.
    Mat3x3 lookAt( const Vec3& _target, const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f));
    Mat4x4 lookAtH( const Vec3& _target, const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f));

    // ********************************************************************* //
    /// \brief Create a camera matrix in homogeneous space.
    /// \details The camera matrix must be in homogeneous space due to the
    ///    translation.
    /// \details This method creates an left-hand system (LHS) with positive
    ///    z-axis.
    Mat4x4 camera( const Vec3& _position, const Vec3& _target, const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f));

    // TODO: orthonormalize()

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
    ///    The OpenGL frustum is defined in the [-1,-1,-1] x [1,1,1] cube.
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

	// Include implementation.
#	include "details/matrix.inl"

    /// \brief Alias for identity<float,2>().
    inline Mat2x2 identity2x2()    { return identity<float,2>(); }
    /// \brief Alias for identity<float,3>().
    inline Mat3x3 identity3x3()    { return identity<float,3>(); }
    /// \brief Alias for identity<float,4>().
    inline Mat4x4 identity4x4()    { return identity<float,4>(); }             // TESTED

	// Remove helper macros.
#	undef RESULT_TYPE
#   undef ENABLE_IF
}
