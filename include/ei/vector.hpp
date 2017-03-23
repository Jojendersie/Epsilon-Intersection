#pragma once

#include <type_traits>
#include <utility>

#include "elementarytypes.hpp"
#include "details/vectortypedef.hpp"
#include "details/vectordetailsA.hpp"

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
#       define ENABLE_IF(condition) typename = typename std::enable_if< (condition) >::type

        /// \brief Construction without initialization. The values are undefined!
        Matrix() {}

        /// \brief Convert a matrix/vector with a different elementary type.
        template<typename T1>
        explicit Matrix(const Matrix<T1,M,N>& _mat1) // TESTED
        {
            for(uint i = 0; i < N * M; ++i)
                this->m_data[i] = static_cast<T>(_mat1[i]);
        }

        /// \brief Forward to base constructors
        template<typename T1>
        explicit Matrix(T1 _a0) :
            details::Components<T,M,N>(_a0)
        {}
        template<typename T1, typename T2, typename... Args>
        Matrix(T1 _a0, T2 _a1, Args... _args) :
            details::Components<T,M,N>(_a0, _a1, std::forward<Args>(_args)...)
        {}

        /// \brief Allow explicit truncation of the dimension sizes.
        template<typename T1, uint M1, uint N1, ENABLE_IF((M < M1 && N <= N1) || (M <= M1 && N < N1))>
        explicit Matrix(const Matrix<T1,M1,N1>& _mat1, uint _rowOff = 0, uint _colOff = 0)
        {
            eiAssert( _rowOff + M <= M1, "Out of boundaries: matrix subsection wrong!" );
            eiAssert( _colOff + N <= N1, "Out of boundaries: matrix subsection wrong!" );
            // This counter avoids one index computation y*N+x in the inner loop
            uint i = 0;
            for(uint y = _rowOff; y < M+_rowOff; ++y)
                for(uint x = _colOff; x < N+_colOff; ++x)
                    this->m_data[i++] = static_cast<T>(_mat1(y, x));
        }

        /// \brief Allow explicit addition in dimension size
        /// \details Adds a row and a column with zeros to a matrix and sets the
        ///    new diagonal element to 1.
        ///    A     =>    A 0
        ///                0 1
        ///
        ///   Appends 1 to vectors: v    =>   (v 1)
        template<typename T1, uint M1, uint N1, ENABLE_IF((M > M1 && N >= N1) || (M >= M1 && N > N1))>
        explicit Matrix(const Matrix<T1,M1,N1>& _mat1) // TESTED
        {
            if(N == 1 || M == 1)
            {
                for(unsigned i = 0; i < N1 * M1; ++i)
                    this->m_data[i] = _mat1[i];
                for(unsigned i = N1 * M1; i < N * M; ++i)
                    this->m_data[i] = static_cast<T>(1);
            } else {
                // Indices for _mat1 and result
                unsigned i = 0, j = 0;
                for(unsigned y = 0; y < M1; ++y)
                {
                    // Copy MxN part
                    for(unsigned x = 0; x < N1; ++x)
                        this->m_data[j++] = _mat1[i++];
                    // New elements at the end of the row is 0
                    for(unsigned x = N1; x < N; ++x)
                        this->m_data[j++] = static_cast<T>(0);
                }
                // Fill new rows
                for(unsigned y = M1; y < M; ++y)
                    for(unsigned x = 0; x < N; ++x)
                        this->m_data[j++] = x == y ? static_cast<T>(1) : static_cast<T>(0);
            }
        }

        /// \brief Access a single element with two indices.
        /// \details Computes the data index _row * N + _col. Therefore
        ///    iterating over _col in inner loops is fastest.
        /// \param [in] _row Zero-based index of the row. For row vectors this
        ///    must be zero.
        /// \param [in] _col Zero-based index of the column. For column vectors
        ///    this must be zero.
        /// \returns Reference with read or write access to the element
        ///    depending on the constness of the matrix.
        T& operator () (uint _row, uint _col) // TESTED
        {
            eiAssertWeak(_row < M && _col < N, "Index out of bounds!");
            return this->m_data[_row * N + _col];
        }
        T operator () (uint _row, uint _col) const // TESTED
        {
            eiAssertWeak(_row < M && _col < N, "Index out of bounds!");
            return this->m_data[_row * N + _col];
        }

        /// \brief Single row access
        /// \param [in] _row Index of the row in [0,M-1].
        Matrix<T,1,N>& operator () (uint _row) // TESTED
        {
            eiAssertWeak(_row < M, "Index out of bounds!");
            return reinterpret_cast<Matrix<T,1,N>&>(this->m_data[_row * N]);
        }
        const Matrix<T,1,N>& operator () (uint _row) const // TESTED
        {
            eiAssertWeak(_row < M, "Index out of bounds!");
            return reinterpret_cast<const Matrix<T,1,N>&>(this->m_data[_row * N]);
        }

        /// \brief Access an element by a single index treating the matrix as 1D.
        /// \param [in] _index Index in the range [0, N * M - 1].
        /// \returns Reference with read or write access to the element
        ///    depending on the constness of the matrix.
        T& operator [] (uint _index) // TESTED
        {
            eiAssertWeak(_index < N * M, "Index out of bounds!");
            return this->m_data[_index];
        }
        T operator [] (uint _index) const // TESTED
        {
            eiAssertWeak(_index < N * M, "Index out of bounds!");
            return this->m_data[_index];
        }

        /// \brief Get access to a subrange [FROM, TO) in the vector.
        /// \tparam FROM First element in the output range (inclusive).
        /// \tparam TO Exclusive right boundary.
        template<uint FROM, uint TO, ENABLE_IF((N == 1) && sizeof(FROM))>
        Matrix<T, TO - FROM, 1>& subcol() // TESTED
        {
            return *reinterpret_cast<Matrix<T, TO - FROM, 1>*>(this->m_data + FROM);
        }
        template<uint FROM, uint TO, ENABLE_IF((N == 1) && sizeof(FROM))>
        const Matrix<T, TO - FROM, 1>& subcol() const
        {
            return *reinterpret_cast<const Matrix<T, TO - FROM, 1>*>(this->m_data + FROM);
        }
        template<uint FROM, uint TO, ENABLE_IF((M == 1) && sizeof(FROM))>
        Matrix<T, 1, TO - FROM>& subrow() // TESTED
        {
            return *reinterpret_cast<Matrix<T, 1, TO - FROM>*>(this->m_data + FROM);
        }
        template<uint FROM, uint TO, ENABLE_IF((M == 1) && sizeof(FROM))>
        const Matrix<T, 1, TO - FROM>& subrow() const
        {
            return *reinterpret_cast<const Matrix<T, 1, TO - FROM>*>(this->m_data + FROM);
        }

        /// \brief Add two matrices component wise.
        /// \details Addition is commutative.
        template<typename T1>
        Matrix<RESULT_TYPE(+), M, N> operator + (const Matrix<T1,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_OP(+)
        /// \brief Subtract two matrices component wise.
        /// \details Subtraction is not commutative.
        template<typename T1>
        Matrix<RESULT_TYPE(-), M, N> operator - (const Matrix<T1,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_OP(-)

        /// \brief Unary minus on all components.
        Matrix<T, M, N> operator - () const // TESTED
            EI_CODE_GEN_MAT_UNARY_OP(-)
        /// \brief Component wise binary not.
        Matrix<T, M, N> operator ~ () const // TESTED
            EI_CODE_GEN_MAT_UNARY_OP(~)

        /// \brief Matrix multiplication.
        /// \details Matrix multiplication is not commutative.
        /// \returns Matrix product with dimensions MxO = MxN * NxO. The result
        ///    is a scalar if M = N = 1.
        template<typename T1, uint O>
        typename std::conditional<M * O == 1, RESULT_TYPE(*), Matrix<RESULT_TYPE(*), M, O>>::type
        operator * (const Matrix<T1,N,O>& _mat1) const // TESTED
        {
            typename std::conditional<M * O == 1, RESULT_TYPE(*), Matrix<RESULT_TYPE(*), M, O>>::type result;
            for(uint m = 0; m < M; ++m)
            {
                for(uint o = 0; o < O; ++o)
                {
                    RESULT_TYPE(*) acc = (*this)(m,0) * _mat1(0,o);
                    for(uint n = 1; n < N; ++n)
                        acc += (*this)(m,n) * _mat1(n,o);
                        *(reinterpret_cast<RESULT_TYPE(*)*>(&result) + m * O + o) = acc;
                }
            }
            return result;
        }

        /// \brief Component wise multiplication for vectors of the same size.
        /// \details For square matrices the above matrix multiplication is used.
        template<typename T1, ENABLE_IF((M != N) && sizeof(T1))>
        Matrix<RESULT_TYPE(*), M, N> operator * (const Matrix<T1,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_OP(*)

        /// \brief Component wise division for vectors of the same size.
        template<typename T1, ENABLE_IF((M != N) && sizeof(T1))>
        Matrix<RESULT_TYPE(/), M, N> operator / (const Matrix<T1,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_OP(/)

        /// \brief Component wise binary or.
        /// \details Or is commutative.
        template<typename T1>
        Matrix<RESULT_TYPE(|), M, N> operator | (const Matrix<T1,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_OP(|)
        /// \brief Component wise binary and.
        /// \details And is commutative.
        template<typename T1>
        Matrix<RESULT_TYPE(&), M, N> operator & (const Matrix<T1,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_OP(&)
        /// \brief Component wise binary xor.
        /// \details Xor is commutative.
        template<typename T1>
        Matrix<RESULT_TYPE(^), M, N> operator ^ (const Matrix<T1,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_OP(^)
        /// \brief Component wise modulo (rest of integer division).
        template<typename T1>
        Matrix<RESULT_TYPE(%), M, N> operator % (const Matrix<T1,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_OP(%)

        /// \brief Self assigning component wise addition.
        template<typename T1>
        Matrix<T, M, N>& operator += (const Matrix<T1,M,N>& _mat1) // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(+=)
        /// \brief Self assigning component wise subtraction.
        template<typename T1>
        Matrix<T, M, N>& operator -= (const Matrix<T1,M,N>& _mat1) // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(-=)
        /// \brief Self assigning component wise multiplication for vectors
        ///    of the same size. Matrix multiplication in case of squared matrices!
        template<typename T1>
        Matrix<T, M, N>& operator *= (const Matrix<T1,M,N>& _mat1) // TESTED
        {
            if(M == N)
                *this = (*this) * _mat1;
            else
                for(uint i = 0; i < M*N; ++i)
                    (*this)[i] *= _mat1[i];
            return *this;
        }
        /// \brief Self assigning component wise division for vectors
        ///    of the same size.
        // TODO: Matrix division (mul with inverse) for squared matrices??
        template<typename T1>
        Matrix<T, M, N>& operator /= (const Matrix<T1,M,N>& _mat1) // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(/=)
        /// \brief Self assigning component wise binary or.
        template<typename T1>
        Matrix<T, M, N>& operator |= (const Matrix<T1,M,N>& _mat1) // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(|=)
        /// \brief Self assigning component wise binary and.
        template<typename T1>
        Matrix<T, M, N>& operator &= (const Matrix<T1,M,N>& _mat1) // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(&=)
        /// \brief Self assigning component wise binary xor.
        template<typename T1>
        Matrix<T, M, N>& operator ^= (const Matrix<T1,M,N>& _mat1) // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(^=)
        /// \brief Self assigning modulo (rest of integer division)
        template<typename T1>
        Matrix<T, M, N>& operator %= (const Matrix<T1,M,N>& _mat1) // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(%=)

        /// \brief Self assigning scalar addition.
        template<typename T1>
        Matrix<T, M, N>& operator += (T1 _s) // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(+=)
        /// \brief Self assigning scalar subtraction.
        template<typename T1>
        Matrix<T, M, N>& operator -= (T1 _s) // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(-=)
        /// \brief Self assigning scalar multiplication.
        template<typename T1>
        Matrix<T, M, N>& operator *= (T1 _s) // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(*=)
        /// \brief Self assigning scalar division.
        template<typename T1>
        Matrix<T, M, N>& operator /= (T1 _s) // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(/=)
        /// \brief Self assigning scalar or.
        template<typename T1>
        Matrix<T, M, N>& operator |= (T1 _s) // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(|=)
        /// \brief Self assigning scalar and.
        template<typename T1>
        Matrix<T, M, N>& operator &= (T1 _s) // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(&=)
        /// \brief Self assigning scalar xor.
        template<typename T1>
        Matrix<T, M, N>& operator ^= (T1 _s) // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(^=)
        /// \brief Self assigning scalar modulo (rest of integer division).
        template<typename T1>
        Matrix<T, M, N>& operator %= (T1 _s) // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(%=)
        /// \brief Self assigning component wise shift
        template<typename T1>
        Matrix<T, M, N>& operator >>= (T1 _s) // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(>>=)
        template<typename T1>
        Matrix<T, M, N>& operator <<= (T1 _s) // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(<<=)

        /// \brief Compare component wise, if two matrices are identical.
        bool operator == (const Matrix<T,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_BOOL_ALL_OP(==)
        /// \brief Compare component wise, if two matrices are distinct.
        bool operator != (const Matrix<T,M,N>& _mat1) const // TESTED
        {
            for(uint i = 0; i < N * M; ++i)
                if((*this)[i] != _mat1[i]) return true;
                    return false;
        }
        /// \brief Compare component wise, if elements are smaller or equal.
        bool operator <= (const Matrix<T,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_BOOL_ALL_OP(<=)
        /// \brief Compare component wise, if elements are smaller.
        bool operator < (const Matrix<T,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_BOOL_ALL_OP(<)
        /// \brief Compare component wise, if elements are greater.
        bool operator > (const Matrix<T,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_BOOL_ALL_OP(>)
        /// \brief Compare component wise, if elements are greater or equal.
        bool operator >= (const Matrix<T,M,N>& _mat1) const // TESTED
            EI_CODE_GEN_MAT_MAT_BOOL_ALL_OP(>=)
    };


    // ********************************************************************* //
    // Scalar operators

    /// \brief Add a scalar to all components.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(+), M, N> operator + (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(+)
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(+), M, N> operator + (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(+)
    /// \brief Subtract a scalar from all components.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(-), M, N> operator - (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(-)
    /// \brief Subtract all components from a scalar.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(-), M, N> operator - (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(-)
    /// \brief Multiply a scalar to all components.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(*), M, N> operator * (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(*)
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(*), M, N> operator * (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(*)
    /// \brief Divide all components by a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(/), M, N> operator / (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(/)
    /// \brief Divide a scalar by each component.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(/), M, N> operator / (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(/)

    /// \brief Binary or of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(|), M, N> operator | (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(|)
    /// \brief Binary or of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(|), M, N> operator | (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(|)
    /// \brief Binary and of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(&), M, N> operator & (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(&)
    /// \brief Binary and of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(&), M, N> operator & (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(&)
    /// \brief Binary xor of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(^), M, N> operator ^ (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(^)
    /// \brief Binary xor of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(^), M, N> operator ^ (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(^)
    /// \brief Modulo of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(%), M, N> operator % (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(%)
    /// \brief Modulo of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    Matrix<RESULT_TYPE(%), M, N> operator % (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(%)
    /// \brief Component wise shift.
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(>>), M, N> operator >> (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(>>)
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(>>), M, N> operator >> (T1 _s, const Matrix<T,M,N>& _mat)
        EI_CODE_GEN_SCALAR_MAT_OP(>>)
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(<<), M, N> operator << (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(<<)
    template<typename T, uint M, uint N, typename T1>
    Matrix<RESULT_TYPE(<<), M, N> operator << (T1 _s, const Matrix<T,M,N>& _mat)
        EI_CODE_GEN_SCALAR_MAT_OP(<<)

    /// \brief Test all components with respect to a scalar.
    template<typename T, uint M, uint N, typename T1>
    bool operator == (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_ALL_OP(==)
    template<typename T1, typename T, uint M, uint N>
    bool operator == (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_ALL_OP(==)
    template<typename T, uint M, uint N, typename T1>
    bool operator != (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            if(_mat[i] != _s) return true;
        return false;
    }
    template<typename T1, typename T, uint M, uint N>
    bool operator != (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            if(_s != _mat[i]) return true;
        return false;
    }
    template<typename T, uint M, uint N, typename T1>
    bool operator < (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_ALL_OP(<)
    template<typename T1, typename T, uint M, uint N>
    bool operator < (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_ALL_OP(<)
    template<typename T, uint M, uint N, typename T1>
    bool operator <= (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_ALL_OP(<=)
    template<typename T1, typename T, uint M, uint N>
    bool operator <= (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_ALL_OP(<=)
    template<typename T, uint M, uint N, typename T1>
    bool operator >= (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_ALL_OP(>=)
    template<typename T1, typename T, uint M, uint N>
    bool operator >= (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_ALL_OP(>=)
    template<typename T, uint M, uint N, typename T1>
    bool operator > (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_ALL_OP(>)
    template<typename T1, typename T, uint M, uint N>
    bool operator > (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_ALL_OP(>)


    /// \brief Compare component wise, if two matrices are identical.
    template<typename T, uint M, uint N>
    Matrix<bool,M,N> equal (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(==)
    /// \brief Compare component wise, if two matrices are distinct.
    template<typename T, uint M, uint N>
    Matrix<bool,M,N> neq (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(!=)
    /// \brief Compare component wise, if elements are smaller or equal.
    template<typename T, uint M, uint N>
    Matrix<bool,M,N> lesseq (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(<=)
    /// \brief Compare component wise, if elements are smaller.
    template<typename T, uint M, uint N>
    Matrix<bool,M,N> less (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(<)
    /// \brief Compare component wise, if elements are greater.
    template<typename T, uint M, uint N>
    Matrix<bool,M,N> greater (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(>)
    /// \brief Compare component wise, if elements are greater or equal.
    template<typename T, uint M, uint N>
    Matrix<bool,M,N> greatereq (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(>=)

    /// \brief Test if all components compare equal to a scalar.
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> equal (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(==)
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> equal (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(==)
    template<typename T, uint M, uint N, typename T1>
    /// \brief Test if any component is non equal to a scalar.
    Matrix<bool, M, N> neq (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(!=)
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> neq (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(!=)
    template<typename T, uint M, uint N, typename T1>
    /// \brief Test if all components compare to a scalar.
    Matrix<bool, M, N> less (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(<)
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> less (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(<)
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> lesseq (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(<=)
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> lesseq (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(<=)
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> greatereq (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(>=)
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> greatereq (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(>=)
    template<typename T, uint M, uint N, typename T1>
    Matrix<bool, M, N> greater (const Matrix<T,M,N>& _mat, T1 _s) // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(>)
    template<typename T1, typename T, uint M, uint N>
    Matrix<bool, M, N> greater (T1 _s, const Matrix<T,M,N>& _mat) // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(>)



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
        TQuaternion( const Vec<T,3>& _axis, T _angle ) // TESTED
        {
            eiAssert( approx(lensq(_axis), 1.0f), "Expected a normalized axis vector!" );
            _angle *= 0.5f;
            T sinA = sin(_angle);
            r = cos(_angle);
            // Assert normalization condition
            if( r < static_cast<T>(0) ) {r = -r; sinA = -sinA;}
            i = sinA * _axis.x;
            j = sinA * _axis.y;
            k = sinA * _axis.z;
        }

        /// \brief Create from Euler angles
        /// \details The rotations are applied in the order x, y, z:
        ///     rotationZ(_z) * rotationY(_y) * rotationX(_x)
        TQuaternion( T _x, T _y, T _z ) // TESTED
        {
            double halfAngle;

            halfAngle = _x * 0.5;
            double sinX = sin(halfAngle);
            double cosX = cos(halfAngle);

            halfAngle = _y * 0.5;
            double sinY = sin(halfAngle);
            double cosY = cos(halfAngle);

            halfAngle = _z * 0.5;
            double sinZ = sin(halfAngle);
            double cosZ = cos(halfAngle);

            double cZcY = cosZ * cosY;
            double cZsY = cosZ * sinY;
            double sZcY = sinZ * cosY;
            double sZsY = sinZ * sinY;

            i = T(sinX * cZcY - cosX * sZsY);
            j = T(cosX * cZsY + sinX * sZcY);
            k = T(cosX * sZcY - sinX * cZsY);
            r = T(cosX * cZcY + sinX * sZsY);

            // Assert normalization condition
            if( r < static_cast<T>(0) )
            {
                r = -r;
                i = -i;
                j = -j;
                k = -k;
            }

            //*this = normalize(*this);
        }
        TQuaternion( const Vec<T,3>& _eulerAngles ) :
            TQuaternion(_eulerAngles.x, _eulerAngles.y, _eulerAngles.z)
        {}

        /// \brief Create from rotation matrix (does a decomposition if the
        ///     matrix contains scaling).
        TQuaternion( const Matrix<T,3,3>& _matrix ) : // TESTED
            TQuaternion<T>(transpose(_matrix(0)), transpose(_matrix(1)), transpose(_matrix(2)))
        {}

        /// \brief Create from orthogonal basis vectors.
        TQuaternion( const Vec<T,3>& _xAxis, const Vec<T,3>& _yAxis, const Vec<T,3>& _zAxis )
        {
            // Check handness
            //eiAssert(dot(cross(_xAxis, _yAxis), _m(2)) > 0.0f, "Quaternions cannot handle reflections. The matrix must be RHS.");
            T handness = dot(cross(_xAxis, _yAxis), _zAxis);
            Vec<T,3> zAxis;
            if(handness < T(0)) zAxis = -_zAxis;
            else zAxis = _zAxis;
            eiAssert(approx(abs(handness), T(1), 1e-4f), "System is not orthonormal!");

            // Build TQuaternion<T> from rotation matrix
            // Src: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
            T trace = _xAxis.x + _yAxis.y + zAxis.z;
            if( trace > 0 )
            {
                float s = T(0.5) / sqrt( trace + T(1) );
                i = (  zAxis.y - _yAxis.z ) * s;
                j = ( _xAxis.z -  zAxis.x ) * s;
                k = ( _yAxis.x - _xAxis.y ) * s;
                r = T(0.25) / s;
            } else {
                if( _xAxis.x > _yAxis.y && _xAxis.x > zAxis.z )
                {
                    float s = T(2) * sqrt( T(1) + _xAxis.x - _yAxis.y - zAxis.z );
                    i = T(0.25) * s;
                    j = ( _xAxis.y + _yAxis.x ) / s;
                    k = ( _xAxis.z +  zAxis.x ) / s;
                    r = (  zAxis.y - _yAxis.z ) / s;
                } else if( _yAxis.y > zAxis.z )
                {
                    float s = T(2) * sqrt( T(1) + _yAxis.y - _xAxis.x - zAxis.z );
                    i = ( _xAxis.y + _yAxis.x ) / s;
                    j = T(0.25) * s;
                    k = ( _yAxis.z +  zAxis.y ) / s;
                    r = ( _xAxis.z -  zAxis.x ) / s;
                } else {
                    float s = T(2) * sqrt( T(1) + zAxis.z - _xAxis.x - _yAxis.y );
                    i = ( _xAxis.z +  zAxis.x ) / s;
                    j = ( _yAxis.z +  zAxis.y ) / s;
                    k = T(0.25) * s;
                    r = ( _yAxis.x - _xAxis.y ) / s;
                }
            }//*/

             /*r = sqrt( max( T(0), T(1) + _m.m00 + _m.m11 + _m.m22 ) ) * T(0.5);
             i = sqrt( max( T(0), T(1) + _m.m00 - _m.m11 - _m.m22 ) ) * T(0.5) * sgn(_m.m21 - _m.m12);
             j = sqrt( max( T(0), T(1) - _m.m00 + _m.m11 - _m.m22 ) ) * T(0.5) * sgn(_m.m02 - _m.m20);
             k = sqrt( max( T(0), T(1) - _m.m00 - _m.m11 + _m.m22 ) ) * T(0.5) * sgn(_m.m10 - _m.m01);//*/

            *this = normalize(*this);

            // Assert additional normalization condition
            if( r * handness < static_cast<T>(0) ) {i = -i; j = -j; k = -k; r = -r;}
        }

        /// \brief Create from TQuaternion coefficients
        TQuaternion( T _i, T _j, T _k, T _r ) :
            i(_i), j(_j), k(_k), r(_r)
        {}

        /// \brief Rotate from vector to vector (rotated such that the from
        ///     vector is aligned with the to vector).
        /// \param [in] _from One certain direction vector before rotation.
        /// \param [in] _to Target direction vector. The from direction should
        ///     be aligned with the target after rotation.
        TQuaternion( const Vec<T,3>& _from, const Vec<T,3>& _to ) // TESTED
        {
            Vec<T,3> from = normalize(_from);
            Vec<T,3> to = normalize(_to);
            // half angle trick from http://physicsforgames.blogspot.de/2010/03/quaternion-tricks.html
            Vec<T,3> half = normalize(from + to);
            // Opposite vectors or one vector 0.0 -> 180° rotation
            if(half.x != half.x)
            {
                if(approx(abs(from.y), 1.0f))
                    half = Vec3(0.0f, 0.0f, 1.0f);
                else half = normalize(cross(from, Vec3(0.0f, 1.0f, 0.0f)));
            }

            // cos(theta) = dot product since both vectors are normalized
            r = dot(from, half);
            eiAssert(r >= T(0), "Normalization condition violated!");
            // Axis from cross product -> already multiplied with sin(theta)
            i = from.y*half.z - from.z*half.y;
            j = from.z*half.x - from.x*half.z;
            k = from.x*half.y - from.y*half.x;
        }

        // TODO: lookAt parametrization

        /// \brief Compare component wise, if two quaternions are identical.
        bool operator == (const TQuaternion& _q1) const
        {
            return r==_q1.r && i==_q1.i && j==_q1.j && k==_q1.k;
        }
        /// \brief Compare component wise, if two quaternions are different.
        bool operator!= (const TQuaternion& _q1) const
        {
            return r!=_q1.r || i!=_q1.i || j!=_q1.j || k!=_q1.k;
        }

        /// \brief TQuaternion multiplication is a combination of rotations.
        /// \details Non commutative (a*b != a*b)
        TQuaternion& operator *= (const TQuaternion& _q1)
        {
            // Preserve the sign for handness: the result can contain a mirroring, if
            // and only if one of the two arguments has a mirroring part.
            T handness = r*_q1.r;
            T nr = r*_q1.r - i*_q1.i - j*_q1.j - k*_q1.k;
            T ni = r*_q1.i + i*_q1.r + j*_q1.k - k*_q1.j;
            T nj = r*_q1.j + j*_q1.r + k*_q1.i - i*_q1.k;
            k = r*_q1.k + k*_q1.r + i*_q1.j - j*_q1.i;
            r = nr;
            i = ni;
            j = nj;
            if( r * handness < static_cast<T>(0) ) {i = -i; j = -j; k = -k; r = -r;}
            return *this;
        }

        /// \brief Scale the TQuaternion
        TQuaternion& operator *= (T _s)
        {
            eiAssert(_s >= T(0), "Using a negative scalar changes handness!");
            i*=_s; j*=_s; k*=_s; r*=_s;
            return *this;
        }

        /// \brief TQuaternion division   a/=b  <=>  a=a*(b^-1)=a*conjugated(b).
        TQuaternion& operator /= (const TQuaternion& _q1)
        {
            T handness = r*_q1.r;
            T nr =   r*_q1.r + i*_q1.i + j*_q1.j + k*_q1.k;
            T ni = - r*_q1.i + i*_q1.r - j*_q1.k + k*_q1.j;
            T nj = - r*_q1.j + j*_q1.r - k*_q1.i + i*_q1.k;
            k = - r*_q1.k + k*_q1.r - i*_q1.j + j*_q1.i;
            r = nr;
            i = ni;
            j = nj;
            if( r * handness < static_cast<T>(0) ) {i = -i; j = -j; k = -k; r = -r;}
            return *this;
        }

        /// \brief Scale the TQuaternion
        TQuaternion& operator /= (T _s)
        {
            eiAssert(_s >= T(0), "Using a negative scalar changes handness!");
            i/=_s; j/=_s; k/=_s; r/=_s;
            return *this;
        }

        /// \brief Vector like addition
        TQuaternion& operator += (const TQuaternion& _q1)
        {
            i+=_q1.i; j+=_q1.j; k+=_q1.k; r+=_q1.r;
            return *this;
        }

        /// \brief Vector like subtraction
        TQuaternion& operator -= (const TQuaternion& _q1)
        {
            i-=_q1.i; j-=_q1.j; k-=_q1.k; r-=_q1.r;
            return *this;
        }

        TQuaternion operator * (const TQuaternion& _q1) const
        {
           // TQuaternion q0 = *this;
            return TQuaternion(*this) *= _q1;
        }
        TQuaternion operator * (T _s) const
        {
            return TQuaternion(*this) *= _s;
        }
        TQuaternion operator / (TQuaternion _q1) const
        {
            return _q1 /= *this;
        }
        TQuaternion operator / (T _s) const
        {
            return TQuaternion(*this) /= _s;
        }
        TQuaternion operator + (TQuaternion _q1) const
        {
            return _q1 += *this;
        }
        TQuaternion operator - (TQuaternion _q1) const
        {
            return _q1 -= *this;
        }

        /// \brief Negate all components, the represented rotation is the same
        TQuaternion operator - () const // TESTED
        {
            return TQuaternion<T>(-i, -j, -k, -r);
        }

        /// \brief Conjugate the quaternion
        TQuaternion operator ~ () const // TESTED
        {
            return TQuaternion<T>(-i, -j, -k, r);
        }

        union {
            struct {T i, j, k, r;};         ///< Elements of 4D complex number
            T z[4];                         ///< Array access, index 3 is the real part
        };
    };

}



#include "details/vectordetailsB.hpp"



namespace ei {

    // ********************************************************************* //
    //                               FUNCTIONS                               //
    // ********************************************************************* //

    // ********************************************************************* //
    /// \brief Returns identity element of the Hamilton-product. (Does not
    ///     rotate anything.)
    inline const TQuaternion<float>& qidentity() // TESTED
    {
        return details::QUATERNION_IDENTITY;
    }

    inline const TQuaternion<double>& qidentityD()
    {
        return details::QUATERNIOND_IDENTITY;
    }

    // ********************************************************************* //
    /// \brief Scalar multiplication from left
    template<typename T>
    inline TQuaternion<T> operator* (T _s, TQuaternion<T> _q)
    {
        return _q *= _s;
    }

    // ********************************************************************* //
    /// \brief Complex conjugate: invert sign of complex components
    template<typename T>
    inline TQuaternion<T> conjugate(const TQuaternion<T>& _q) // TESTED
    {
        return TQuaternion<T>(-_q.i, -_q.j, -_q.k, _q.r);
    }

    /// \brief Get the rotation axis from a TQuaternion
    template<typename T>
    inline Vec<T,3> axis(const TQuaternion<T>& _q) // TESTED
    {
        return Vec<T,3>(-_q.i, -_q.j, -_q.k) / max(T(EPSILON), std::sqrt(T(1)-_q.r*_q.r));
    }

    /// \brief Get the x axis of the corresponding orthogonal system (rotation
    ///     matrix)
    template<typename T>
    inline Vec<T,3> xaxis(const TQuaternion<T>& _q) // TESTED
    {
        return Vec<T,3>( T(1)-T(2)*(_q.j*_q.j+_q.k*_q.k), T(2)*(_q.i*_q.j-_q.k*_q.r), T(2)*(_q.i*_q.k+_q.j*_q.r) );
    }

    /// \brief Get the y axis of the corresponding orthogonal system (rotation
    ///     matrix)
    template<typename T>
    inline Vec<T,3> yaxis(const TQuaternion<T>& _q) // TESTED
    {
        return Vec<T,3>( T(2)*(_q.i*_q.j+_q.k*_q.r), T(1)-T(2)*(_q.i*_q.i+_q.k*_q.k), T(2)*(_q.j*_q.k-_q.i*_q.r) );
    }

    /// \brief Get the z axis of the corresponding orthogonal system (rotation
    ///     matrix)
    template<typename T>
    inline Vec<T,3> zaxis(const TQuaternion<T>& _q) // TESTED
    {
        T h = _q.r < T(0) ? T(-1) : T(1);
        T h2 = h * 2;
        return Vec<T,3>( h2*(_q.i*_q.k-_q.j*_q.r), h2*(_q.j*_q.k+_q.i*_q.r), h-h2*(_q.i*_q.i+_q.j*_q.j) );
    }
    // TODO: row vector axis

    /// \brief Get the angle (radians) from a TQuaternion
    template<typename T>
    inline T angle(const TQuaternion<T>& _q) // TESTED
    {
        return acos(_q.r) * T(2);
    }

    // ********************************************************************* //
    /// \brief Get the Euler angles (radians) from a quaternion
    template<typename T>
    inline Vec<T,3> angles(const TQuaternion<T>& _q)
    {
        // TODO: handness?
        Vec<T,3> angles;
        // Derivation from http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToEuler/index.htm
        // but changed angles because of else convention
        const double m20half = _q.j * _q.r - _q.i * _q.k;

        if(approx(m20half, 0.5))
        {
            angles.x = 0.0f;
            angles.y = PI/2.0f;
            angles.z = static_cast<T>(-2.0 * atan2(_q.i, _q.r));
        }
        else if(approx(m20half, -0.5))
        {
            angles.x = 0.0f;
            angles.y = -PI/2.0f;
            angles.z = static_cast<T>(2.0 * atan2(_q.i, _q.r));
        }
        else
        {
            const double sqr = _q.r * _q.r;
            const double sqi = _q.i * _q.i;
            const double sqj = _q.j * _q.j;
            const double sqk = _q.k * _q.k;
            angles.x = static_cast<T>(atan2(2.0 * (_q.j * _q.k + _q.i * _q.r), -sqi - sqj + sqk + sqr));
            angles.y = static_cast<T>(asin( clamp(m20half * 2.0, -1.0, 1.0) ));
            angles.z = static_cast<T>(atan2(2.0 * (_q.i * _q.j + _q.k * _q.r),  sqi - sqj - sqk + sqr));
        }
        return angles;
    }

    // ********************************************************************* //
    /// \brief Check if the absolute difference between all elements is smaller
    ///    or equal than epsilon.
    /// \param [in] _mat0 First operand.
    /// \param [in] _mat1 Second operand.
    /// \param [in] _epsilon Maximum threshold for the difference between two
    ///    components. The default value is 1e-6.
    /// \returns true if all differences are less or equal than _epsilon.
    template<typename T, unsigned M, unsigned N>
    inline bool approx(const Matrix<T,M,N>& _mat0,
                       const Matrix<T,M,N>& _mat1,
                       T _epsilon = T(1e-6))  // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            // if(!(abs(_mat1[i] - _mat0[i]) <= _epsilon)) return false;
            if(!approx(_mat1[i], _mat0[i], _epsilon))
                return false;
                return true;
    }

    // ********************************************************************* //
    /// \brief Check if the absolute difference between all elements is smaller
    ///    or equal than epsilon.
    template<typename T>
    bool approx(const TQuaternion<T>& _q0,
                const TQuaternion<T>& _q1,
                T _epsilon = T(1e-6)) // TESTED
    {
        return abs(_q0.r - _q1.r) <= _epsilon
            && abs(_q0.i - _q1.i) <= _epsilon
            && abs(_q0.j - _q1.j) <= _epsilon
            && abs(_q0.k - _q1.k) <= _epsilon;
    }

    // ********************************************************************* //
    /// \brief Computes the sum of component wise products.
    /// \returns Scalar value of the sum of component products.
    template<typename T, unsigned M, unsigned N, typename T1>
    inline RESULT_TYPE(*) dot(const Matrix<T,M,N>& _mat0,
                              const Matrix<T1,M,N>& _mat1) // TESTED
    {
        RESULT_TYPE(*) sum = _mat0[0] * _mat1[0];
        for(uint i = 1; i < N * M; ++i)
            sum += _mat0[i] * _mat1[i];
        return sum;
    }

    // ********************************************************************* //
    /// \brief Computes the sum of component wise products.
    /// \returns Scalar value of the sum of component products.
    inline float dot(const Quaternion& _q0,
                     const Quaternion& _q1)
    {
        return _q0.r*_q1.r + _q0.i*_q1.i + _q0.j*_q1.j + _q0.k*_q1.k;
    }

    // ********************************************************************* //
    /// \brief Computes the cross product of two 3d vectors (RHS).
    /// \returns Perpendicular vector with length |v0|·|v1|·sin(∡(v0,v1)).
    template<typename T, typename T1, unsigned M, unsigned N, ENABLE_IF((N==1 && M==3) || (N==3 && M==1))>
    inline Matrix<RESULT_TYPE(*),M,N> cross(const Matrix<T,M,N>& _v0,
                                            const Matrix<T1,M,N>& _v1)
    {
        return Matrix<RESULT_TYPE(*),M,N>(_v0.y * _v1.z - _v0.z * _v1.y,
            _v0.z * _v1.x - _v0.x * _v1.z,
            _v0.x * _v1.y - _v0.y * _v1.x);
    }

    // ********************************************************************* //
    /// \brief Computes the cross product of two 2d vectors.
    /// \returns The determinant of the 2x2 matrix: v0.x·v1.y - v0.y·v1.x.
    template<typename T, typename T1, unsigned M, unsigned N, ENABLE_IF((N==1 && M==2) || (N==2 && M==1))>
    inline RESULT_TYPE(*) cross(const Matrix<T,M,N>& _v0,
                                const Matrix<T1,M,N>& _v1)
    {
        return _v0.x * _v1.y - _v0.y * _v1.x;
    }

    // ********************************************************************* //
    /// \brief Computes the sum of squared components.
    /// \details This is equivalent to dot(_mat0, _mat0).
    /// \returns Squared euclidean length (scalar).
    template<typename T>
    inline auto lensq(const T& _elem0) -> decltype(dot(_elem0, _elem0)) // TESTED
    {
        return dot(_elem0, _elem0);
    }

    // ********************************************************************* //
    /// \brief Computes the root of the sum of squared components.
    /// \details This is the euclidean length for vectors and Quaternions and
    ///    the Frobenius norm for matrices.
    /// \returns Euclidean length (scalar).
    template<typename T>
    inline auto len(const T& _elem0) -> decltype(std::sqrt(dot(_elem0, _elem0))) // TESTED
    {
        return sqrt(dot(_elem0, _elem0));
    }

    // ********************************************************************* //
    /// \brief Normalizes a vector, quaternion or matrix with respect to len.
    /// \details This is equivalent to elem0 / len(_elem0).
    /// \returns Normalized vector or matrix.
    template<typename T>
    inline T normalize(const T& _mat0) // TESTED
    {
        return _mat0 / len(_mat0);
    }

    // ********************************************************************* //
    /// \brief Component wise maximum.
    /// \returns A matrix with the maximum values from both inputs.
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> max(const Matrix<T,M,N>& _mat0,
                             const Matrix<T,M,N>& _mat1) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = max(_mat0[i], _mat1[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Component wise minimum.
    /// \returns A matrix with the minimum values from both inputs.
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> min(const Matrix<T,M,N>& _mat0,
                             const Matrix<T,M,N>& _mat1) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = min(_mat0[i], _mat1[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Maximum element from a matrix.
    /// \returns Scalar maximum value.
    template<typename T, unsigned M, unsigned N>
    inline T max(const Matrix<T,M,N>& _mat0) // TESTED
    {
        T result = _mat0[0];
        for(uint i = 1; i < N * M; ++i)
            result = max(_mat0[i], result);
        return result;
    }

    // ********************************************************************* //
    /// \brief Minimum element from a matrix.
    /// \returns Scalar minimum value.
    template<typename T, unsigned M, unsigned N>
    inline T min(const Matrix<T,M,N>& _mat0) // TESTED
    {
        T result = _mat0[0];
        for(uint i = 1; i < N * M; ++i)
            result = min(_mat0[i], result);
        return result;
    }

    // ********************************************************************* //
    /// \brief Component wise clamp to boundaries.
    /// \returns A matrix with values in the bounding box.
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> clamp(const Matrix<T,M,N>& _mat,
                               const Matrix<T,M,N>& _min,
                               const Matrix<T,M,N>& _max) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = clamp(_mat[i], _min[i], _max[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Component wise clamp to scalar boundaries.
    /// \returns A matrix with values in the interval.
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> clamp(const Matrix<T,M,N>& _mat,
                               T _min,
                               T _max) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = clamp(_mat[i], _min, _max);
        return result;
    }

    /// \brief Clamp all components to [0,1]
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> saturate(const Matrix<T,M,N>& _mat) // TESTED
    {
        return clamp(_mat, static_cast<T>(0), static_cast<T>(1));
    }

    // ********************************************************************* //
    /// \brief Round all components towards negative infinity
    template<typename T, unsigned M, unsigned N>
    inline Matrix<typename details::Int<sizeof(T)>::stype,M,N> floor(const Matrix<T,M,N>& _mat) // TESTED
    {
        Matrix<typename details::Int<sizeof(T)>::stype,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = floor(_mat[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Round all components towards negative infinity
    template<typename T, unsigned M, unsigned N>
    inline Matrix<typename details::Int<sizeof(T)>::stype,M,N> ceil(const Matrix<T,M,N>& _mat) // TESTED
    {
        Matrix<typename details::Int<sizeof(T)>::stype,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = ceil(_mat[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Round all components towards next number (x.5 rounds up)
    template<typename T, unsigned M, unsigned N>
    inline Matrix<typename details::Int<sizeof(T)>::stype,M,N> round(const Matrix<T,M,N>& _mat) // TESTED
    {
        Matrix<typename details::Int<sizeof(T)>::stype,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = round(_mat[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Get the smallest positive number m, such that x=y*c+m with c
    ///     in Z, for each component.
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> mod(const Matrix<T,M,N>& _x, T _y) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = mod(_x[i], _y);
        return result;
    }

    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> mod(const Matrix<T,M,N>& _x, const Matrix<T,M,N>& _y)
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = mod(_x[i], _y[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Compute the square root for each component
    template<typename T, unsigned M, unsigned N, ENABLE_IF((N==1) || (M==1))>
    inline Matrix<T,M,N> sqrt(const Matrix<T,M,N>& _v0) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < M * N; ++i)
            result[i] = std::sqrt(_v0[i]);
        return result;
    }

    // Unfortunately, the above declaration hides the elementary one -> make it
    // visible again
    using std::sqrt;

    // ********************************************************************* //
    /// \brief Compute the power for each component
    template<typename T, unsigned M, unsigned N, ENABLE_IF((N==1) || (M==1))>
    inline Matrix<T,M,N> pow(const Matrix<T,M,N>& _v0, float _exponent) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < M * N; ++i)
            result[i] = std::pow(_v0[i], _exponent);
        return result;
    }

    // Unfortunately, the above declaration hides the elementary one -> make it
    // visible again
    using std::pow;

    // ********************************************************************* //
    /// \brief Element wise natural logarithm for matrices (basis e).
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> log(const Matrix<T,M,N>& _mat0) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < M * N; ++i)
            result[i] = std::log(_mat0[i]);
        return result;
    }

    /// \brief Element wise logarithm for matrices (basis 2).
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> log2(const Matrix<T,M,N>& _mat0) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < M * N; ++i)
            result[i] = std::log2(_mat0[i]);
        return result;
    }

    using std::log;
    using std::log2;

    // ********************************************************************* //
    /// \brief Sum of all components.
    /// \details Can be used for boolean vectors/matrices too (number of trues).
    /// \returns Scalar sum of all values.
    template<typename T, unsigned M, unsigned N>
    inline decltype(std::declval<T>() + std::declval<T>()) sum(const Matrix<T,M,N>& _mat0) // TESTED
    {
        decltype(std::declval<T>() + std::declval<T>()) result = _mat0[0];
        for(uint i = 1; i < N * M; ++i)
            result += _mat0[i];
        return result;
    }

    // ********************************************************************* //
    /// \brief Product of all components.
    /// \returns Product of all values (scalar).
    template<typename T, unsigned M, unsigned N>
    inline T prod(const Matrix<T,M,N>& _mat0) // TESTED
    {
        T result = _mat0[0];
        for(uint i = 1; i < N * M; ++i)
            result *= _mat0[i];
        return result;
    }

    // ********************************************************************* //
    /// \brief Average of all values from a matrix.
    /// \returns Scalar average value.
    template<typename T, unsigned M, unsigned N>
    inline T avg(const Matrix<T,M,N>& _mat0) // TESTED
    {
        return sum(_mat0) / T(M * N);
    }

    // ********************************************************************* //
    /// \brief Absolute values for all components.
    /// \returns Matrix with component wise absolute values.
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> abs(const Matrix<T,M,N>& _mat0) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = abs(_mat0[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Get the component wise sign.
    /// \details There is a faster version sgn(), if you don't need to 
    ///    know about zero.
    /// \returns -1 (_x < 0), 0 (_x == 0) or 1 (_x > 0) for each component _x.
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> sign(const Matrix<T,M,N>& _mat0) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = sign(_mat0[i]);
        return result;
    }


    // ********************************************************************* //
    /// \brief Get the component wise sign where 0 is counted as positive.
    /// \details This function is faster than sign(). Use it if you don't need
    ///    to know about zero.
    /// \returns -1 (_x < 0) or 1 (_x >= 0) for each component _x.
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> sgn(const Matrix<T,M,N>& _mat0) // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = sgn(_mat0[i]);
        return result;
    }

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
    inline Matrix<decltype(std::declval<T0>() * std::declval<T1>()),M,N>
        bilerp(Matrix<T0,M,N> _x00, Matrix<T0,M,N> _x01,
               Matrix<T0,M,N> _x10, Matrix<T0,M,N> _x11,
               T1 _t0, T1 _t1) // TESTED
    {
        // For vectors it is faster to compute the four factors first
        T1 one = static_cast<T1>(1);
        T1 t0i = (one - _t0);
        T1 t1i = (one - _t1);
        T1 t00 = t0i * t1i;
        T1 t01 = _t0 * t1i;
        T1 t10 = t0i * _t1;
        T1 t11 = _t0 * _t1;
        return _x00 * t00 + _x01 * t01
            + _x10 * t10 + _x11 * t11;
    }

    // ********************************************************************* //
    /// \brief Spherical linear interpolation with constant angular speed
    /// \details Spherical interpolation is not defined, if the angle between
    ///     both vectors is pi. In that case this still method performs
    ///     an interpolation but without normalization.
    ///
    ///     Formulas from Ken Shoemake "Animating rotation with quaternion
    ///     curves" SIGGRAPH 85
    template<typename T0, typename T1, unsigned M, unsigned N, ENABLE_IF((N==1) || (M==1))>
    inline auto slerp(const Matrix<T0,M,N>& _v0, const Matrix<T0,M,N>& _v1, T1 _t) -> decltype(_v0*_t)
    {
        T1 theta = acos( clamp(dot(_v0,_v1), static_cast<T1>(-1.0), static_cast<T1>(1.0)) );
        T1 so = sin( theta );
        // Special cases for so->0 reduce to linear interpolation
        if(so == static_cast<T1>(0)) return _v0 + (_v1 - _v0) * _t;
        T1 f0 = sin( theta * (static_cast<T1>(1.0)-_t) ) / so;
        T1 f1 = sin( theta * _t ) / so;
        decltype(_v0*_t) result;
        for(uint i = 0; i < M * N; ++i)
            result[i] = _v0[i] * f0 + _v1[i] * f1;
        return result;
    }

    template<typename T>
    TQuaternion<T> slerp(const TQuaternion<T>& _q0, const TQuaternion<T>& _q1, T _t) // TESTED
    {
        // TODO: handness?
        // http://en.wikipedia.org/wiki/Slerp
        T theta = acos( clamp(dot(_q0,_q1), T(-1), T(1)) );
        T so = sin( theta );
        if(approx(so, T(0)))
        {
            // Converges towards linear interpolation for small so
            return TQuaternion<T>(_q0.i + (_q1.i - _q0.i) * _t,
                                  _q0.j + (_q1.j - _q0.j) * _t,
                                  _q0.k + (_q1.k - _q0.k) * _t,
                                  _q0.r + (_q1.r - _q0.r) * _t);
        }
        T f0 = sin( theta * (1.0f-_t) ) / so;
        T f1 = sin( theta * _t ) / so;
        return TQuaternion<T>(_q0.i * f0 + _q1.i * f1,
                              _q0.j * f0 + _q1.j * f1,
                              _q0.k * f0 + _q1.k * f1,
                              _q0.r * f0 + _q1.r * f1);
    }


    // ********************************************************************* //
    /// \brief Test if at least one element of the matrix is true.
    /// \return false, if all elements off the matrix are false.
    template<unsigned M, unsigned N>
    inline bool any(const Matrix<bool,M,N>& _mat0) // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            if(_mat0[i]) return true;
        return false;
    }

    // ********************************************************************* //
    /// \brief Test if no element of the matrix is true.
    /// \return true, if all elements off the matrix are false.
    template<unsigned M, unsigned N>
    inline bool none(const Matrix<bool,M,N>& _mat0) // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            if(_mat0[i]) return false;
        return true;
    }

    // ********************************************************************* //
    /// \brief Test if all elements of the matrix are true.
    /// \return true, if all elements off the matrix are true.
    template<unsigned M, unsigned N>
    bool all(const Matrix<bool,M,N>& _mat0) // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            if(!_mat0[i]) return false;
        return true;
    }


    // ********************************************************************* //
    /// \brief Transpose a matrix or vector (switch the dimensions).
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,N,M> transpose(const Matrix<T,M,N>& _mat0) // TESTED
    {
        Matrix<T,N,M> result;
        // This counter avoids one index computation y*N+x in the inner loop
        uint i = 0;
        for(uint x = 0; x < N; ++x)
            for(uint y = 0; y < M; ++y)
                result[i++] = _mat0(y,x);
        return result;
    }

    // ********************************************************************* //
    /// \brief Gram-Schmidt orthonormalization for matrices with M >= N.
    /// \details This method assumes a column space where each column is one
    ///     vector.
    ///     The method operates a stabilized Gram-Schmidt inplace.
    /// \returns false if some columns are linear dependent and not all vectors
    ///     can be orthogonalized.
    template<typename T, unsigned M, unsigned N>
    inline bool orthonormalize(Matrix<T,M,N>& _mat0) // TESTED
    {
        static_assert( M >= N, "Number of vectors N must be smaller than their dimension to be orthogonal." );

        // For each column
        for(uint x = 0; x < N; ++x)
        {
            // Normalize column
            float norm = sq(_mat0(0,x));
            for(uint y = 1; y < M; ++y)
                norm += sq(_mat0(y,x));
            if(norm <= 1e-30f) return false;
            norm = sqrt(norm);
            for(uint y = 0; y < M; ++y)
                _mat0(y,x) /= norm;

            // Remove current direction from all following vectors.
            for(uint j = x+1; j < N; ++j)
            {
                // Dot product for projection
                float projScale = _mat0(0,j) * _mat0(0,x);
                for(uint y = 1; y < M; ++y)
                    projScale += _mat0(y,j) * _mat0(y,x);
                // Projection and subtraction in one step.
                for(uint y = 0; y < M; ++y)
                    _mat0(y,j) -= projScale * _mat0(y,x);
            }
        }

        return true;
    }

    /// \brief Gram-Schmidt orthonormalization for a list of vectors.
    template<typename TVec0, typename... TVecs>
    inline bool orthonormalize(TVec0& _vec0, TVecs&... _vecs) // TESTED
    {
        float norm = len(_vec0);
        if(norm <= 1e-30f) return false;
        _vec0 /= norm;

        // Remove current vector from all following
        details::removeProjectedPart(_vec0, _vecs...);
        // continue with the next vector as reference
        return orthonormalize(_vecs...);
    }


    // ********************************************************************* //
    //                            TRANSFORMATIONS                            //
    // ********************************************************************* //

    // ********************************************************************* //
    /// \brief Generate the N x N identity matrix.
    template<typename T, unsigned N>
    inline const Matrix<T,N,N>& identity() // TESTED
    {
        static Matrix<T,N,N> result(diag(Vec<T,N>(1)));
        return result;
    }

    // Faster implementations for known sizes (does not branch due to static)
    template<>
    inline const Matrix<float,2,2>& identity<float,2>()
    {
        return details::MAT2X2_IDENTITY;
    }
    template<>
    inline const Matrix<float,3,3>& identity<float,3>()
    {
        return details::MAT3X3_IDENTITY;
    }
    template<>
    inline const Matrix<float,4,4>& identity<float,4>()
    {
        return details::MAT4X4_IDENTITY;
    }

    /// \brief Alias for identity<float,2>().
    inline Mat2x2 identity2x2()    { return identity<float,2>(); }
    /// \brief Alias for identity<float,3>().
    inline Mat3x3 identity3x3()    { return identity<float,3>(); } // TESTED
    /// \brief Alias for identity<float,4>().
    inline Mat4x4 identity4x4()    { return identity<float,4>(); } // TESTED

    // ********************************************************************* //
    /// \brief Generate the N x N diagonal matrix.
    /// \param [in] _v0 A vector with the diagonal entries.
    template<typename T, unsigned N>
    inline Matrix<T,N,N> diag( const Vec<T,N>& _v0 ) // TESTED
    {
        Matrix<T,N,N> result(T(0));
        for(uint n = 0; n < N; ++n)
            result[n * N + n] = _v0[n];
        return result;
    }

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
    inline Vec<T,N> sphericalCoords( const Vec<T,N>& _v0 ) // TESTED
    {
        static_assert(N >= 2, "In 1D cartesian and spherical coordinates are the same!");
        Vec<T,N> result;
        // Accumulate the squared length over the iterations
        result[0] = sq(_v0[N-1]) + sq(_v0[N-2]);
        result[N-1] = atan2(_v0[N-1], _v0[N-2]);
        if(result[N-1] < 0.0f) result[N-1] += 2.0f * PI;
        // if( _v0[N-1] < 0.0f ) result[N-1] = 2.0f * PI - result[N-1];
        for(uint i = 2; i < N; ++i)
        {
            result[0] += sq(_v0[N-i-1]);
            result[N-i] = acos(clamp(_v0[N-i-1]/sqrt(result[0]), -1.0f, 1.0f));
        }
        result[0] = sqrt(result[0]);
        return result;
    }

    template<typename T, unsigned N>
    inline RVec<T,N> sphericalCoords( const RVec<T,N>& _v0 ) // TESTED
    {
        return *reinterpret_cast<RVec<T,N>*>(&sphericalCoords(*reinterpret_cast<Vec<T,N>*>(&_v0)));
    }

    // ********************************************************************* //
    /// \brief Convert a vector from spherical coordinates (r, φ1, φ02, ..., φN-1)
    ///     to regular Cartesian coordinates.
    /// \return The regular Cartesian vector.
    template<typename T, unsigned N>
    inline Vec<T,N> cartesianCoords( const Vec<T,N>& _v0 ) // TESTED
    {
        eiAssertWeak(_v0[0] > 0.0f, "Expected the length to be greater 0!");
        Vec<T,N> result;
        float tmp = _v0[0];
        for(uint i = 0; i < N-1; ++i)
        {
            result[i] = tmp * cos(_v0[i+1]);
            tmp *= sin(_v0[i+1]);
        }
        result[N-1] = tmp;
        return result;
    }

    template<typename T, unsigned N>
    inline RVec<T,N> cartesianCoords( const RVec<T,N>& _v0 ) // TESTED
    {
        return *reinterpret_cast<RVec<T,N>*>(&cartesianCoords(*reinterpret_cast<Vec<T,N>*>(&_v0)));
    }

    // ********************************************************************* //
    /// \brief Apply transformations in homogeneous space. This includes a
    ///     division by w after the transformation
    template<typename T, unsigned N>
    inline Matrix<T,N,1> transformDiv( const Matrix<T, N, 1>& _what,
                                       const Matrix<T, N+1, N+1>& _space )
    {
        T t[N+1];
        // Multiply Matrix * Vector(_what,1)
        for(uint y = 0; y <= N; ++y)
        {
            // Initialize with the last component * 1
            t[y] = _space(y, N);
            // Add the other N factors
            for(uint x = 0; x < N; ++x)
                t[y] += _space(y,x) * _what[x];
        }
        // Create a reduced vector with divided components
        Matrix<T,N,1> result;
        for(uint i = 0; i < N; ++i)
            result[i] = t[i] / t[N];
        return result;
    }
    template<typename T, unsigned N>
    inline Matrix<T,1,N> transformDiv( const Matrix<T, 1, N>& _what,
                                       const Matrix<T, N+1, N+1>& _space )
    {
        T t[N+1];
        // Multiply Vector(_what,1) * Matrix
        for(uint x = 0; x <= N; ++x)
        {
            // Initialize with the last component * 1
            t[x] = _space(N, x);
            // Add the other N factors
            for(uint y = 0; y < N; ++y)
                t[x] += _what[y] * _space(y,x);
        }
        // Create a reduced vector with divided components
        Matrix<T,1,N> result;
        for(uint i = 0; i < N; ++i)
            result[i] = t[i] / t[N];
        return result;
    }

    /// \brief Apply transformations in 3x4/4x3 space (rotation + translation).
    ///     This does NOT include a division by w.
    template<typename T, unsigned N>
    inline Matrix<T,N,1> transform( const Matrix<T, N, 1>& _what,
                                    const Matrix<T, N+1, N+1>& _space ) // TESTED
    {
        Matrix<T,N,1> result;
        // Multiply Matrix * Vector(_what,1)
        for(uint y = 0; y < N; ++y)
        {
            // Initialize with the last component * 1
            result[y] = _space(y, N);
            // Add the other N factors
            for(uint x = 0; x < N; ++x)
                result[y] += _space(y,x) * _what[x];
        }
        return result;
    }
    template<typename T, unsigned N>
    inline Matrix<T,1,N> transform( const Matrix<T, 1, N>& _what,
                                    const Matrix<T, N+1, N+1>& _space )
    {
        Matrix<T,1,N> result;
        // Multiply Vector(_what,1) * Matrix
        for(uint x = 0; x < N; ++x)
        {
            // Initialize with the last component * 1
            result[x] = _space(N, x);
            // Add the other N factors
            for(uint y = 0; y < N; ++y)
                result[x] += _what[y] * _space(y,x);
        }
        return result;
    }

    /// \brief Apply transformations and ignore the fourth component and the
    ///     translation.
    /// \details Use this function to transform direction vectors. Still, there
    ///     is no normalization involved and the direction might be scaled by
    ///     the matrix.
    template<typename T, unsigned N>
    inline Matrix<T,N,1> transformDir( const Matrix<T, N, 1>& _what,
                                       const Matrix<T, N+1, N+1>& _space ) // TESTED
    {
        Matrix<T,N,1> result;
        // Multiply Matrix * Vector(_what,0)
        for(uint y = 0; y < N; ++y)
        {
            // Initialize with the first component
            result[y] = _space(y,0) * _what[0];
            // Add the other N-1 factors
            for(uint x = 1; x < N; ++x)
                result[y] += _space(y,x) * _what[x];
        }
        return result;
    }
    template<typename T, unsigned N>
    inline Matrix<T,1,N> transformDir( const Matrix<T, 1, N>& _what,
                                       const Matrix<T, N+1, N+1>& _space )
    {
        Matrix<T,1,N> result;
        // Multiply Vector(_what,0) * Matrix
        for(uint x = 0; x < N; ++x)
        {
            // Initialize with the first component
            result[x] = _space(0,x) * _what[0];
            // Add the other N-1 factors
            for(uint y = 1; y < N; ++y)
                result[x] += _what[y] * _space(y,x);
        }
        return result;
    }

    /// \brief Apply transformations with a matrix multiplication.
    template<typename T, unsigned M>
    inline Matrix<T,M,1> transform( const Matrix<T,M,1>& _what,
                                    const Matrix<T,M,M>& _space )
    {
        return _space * _what;
    }
    template<typename T, unsigned N>
    inline Matrix<T,1,N> transform( const Matrix<T,1,N>& _what,
                                    const Matrix<T,N,N>& _space )
    {
        return _what * _space;
    }

    /// \brief Apply a rotation by a quaternion (q v q-1 with v=(0, _v.x, _v.y, _v.z)).
    template<typename T, unsigned M, unsigned N, ENABLE_IF((M==1) || (N==1))>
    inline Matrix<T,M,N> transform( const Matrix<T,M,N>& _what, const TQuaternion<T>& _quaternion )
    {
        T handness = _quaternion.r < T(0) ? T(-1) : T(1);
        // http://physicsforgames.blogspot.de/2010/03/quaternion-tricks.html
        T x1 = _quaternion.j*_what.z - _quaternion.k*_what.y;
        T y1 = _quaternion.k*_what.x - _quaternion.i*_what.z;
        T z1 = _quaternion.i*_what.y - _quaternion.j*_what.x;

        return Matrix<T,M,N>(
             _what.x + 2.0f * (_quaternion.r*x1 + _quaternion.j*z1 - _quaternion.k*y1),
             _what.y + 2.0f * (_quaternion.r*y1 + _quaternion.k*x1 - _quaternion.i*z1),
            (_what.z + 2.0f * (_quaternion.r*z1 + _quaternion.i*y1 - _quaternion.j*x1)) * handness
            );
    }

    // ********************************************************************* //
    /// \brief Create a translation matrix in homogeneous coordinate space.
    /// \param [in] _vector Translate by/Add this vector.
    /// \details The translation matrix always has a dimension one large then
    ///    the vectors.
    ///    To transform a vector append 1 and multiply it from right:
    ///    translation() * VecX(v,1)
    template<typename T, unsigned N>
    inline Matrix<T,N+1,N+1> translation( const Matrix<T, N, 1>& _vector )
    {
        Matrix<T,N+1,N+1> result = identity<T,N+1>();
        for(uint i = 0; i < N; ++i)
            result[i * (N+1) + N] = _vector[i];
        return result;
    }

    // ********************************************************************* //
    /// \brief Create a scaling/diagonal matrix from vector.
    template<typename T, unsigned N>
    inline Matrix<T,N,N> scaling( const Matrix<T, N, 1>& _scale )
    {
        Matrix<T,N,N> result(T(0));
        for(uint n = 0; n < N; ++n)
            result[n * N + n] = _scale[n];
        return result;
    }

    // ********************************************************************* //
    /// \brief Create a uniform scaling/diagonal matrix from scalar.
    template<typename T, unsigned N>
    inline Matrix<T,N,N> scaling( T _scale )
    {
        Matrix<T,N,N> result(T(0));
        for(uint n = 0; n < N; ++n)
            result[n * N + n] = _scale;
        return result;
    }

    // ********************************************************************* //
    /// \brief Use vectors which span a space to build a matrix.
    /// \details The vectors become the rows of the matrix
    // TODO: variadic template variant
    inline Mat2x2 axis( const Vec2& _x, const Vec2& _y )
    {
        return Mat2x2(_x.x, _x.y,
                      _y.x, _y.y);
    }
    inline Mat3x3 axis( const Vec3& _x, const Vec3& _y, const Vec3& _z ) // TESTED
    {
        return Mat3x3(_x.x, _x.y, _x.z,
                      _y.x, _y.y, _y.z,
                      _z.x, _z.y, _z.z);
    }
    inline Mat4x4 axis( const Vec4& _x, const Vec4& _y, const Vec4& _z, const Vec4& _w )
    {
        return Mat4x4(_x.x, _x.y, _x.z, _x.w,
                      _y.x, _y.y, _y.z, _y.w,
                      _z.x, _z.y, _z.z, _z.w,
                      _w.x, _w.y, _w.z, _w.w);
    }

    // ********************************************************************* //
    /// \brief Create an orthonormal basis for a single direction vector
    inline Mat2x2 basis( const Vec2& _vector ) // TESTED
    {
        return Mat2x2(_vector.x, _vector.y,
                     -_vector.y, _vector.x);
    }

    inline Mat3x3 basis( const Vec3& _vector ) // TESTED
    {
        eiAssert(approx(len(_vector), 1.0f), "Expected normalized direction vector!");
        Vec3 y;
        if(abs(_vector.x) >= 1.0f) y = Vec3(0.0f, 1.0f, 0.0f);
        else y = normalize(Vec3(0.0f, -_vector.z, _vector.y));
        return axis(_vector, y, cross(_vector, y));
    }

    // ********************************************************************* //
    /// \brief Rotation matrix in 2D.
    inline Mat2x2 rotation( float _angle )
    {
        float sinA = sin(_angle);
        float cosA = cos(_angle);
        return Mat2x2(cosA, -sinA,
                      sinA,  cosA);
    }

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D space around x-axis.
    inline Mat3x3 rotationX( float _angle )
    {
        float sinA = sin(_angle);
        float cosA = cos(_angle);
        return Mat3x3(1.0f, 0.0f,  0.0f,
                      0.0f, cosA, -sinA,
                      0.0f, sinA,  cosA);
    }


    // ********************************************************************* //
    /// \brief Rotation matrix in 3D space around y-axis.
    inline Mat3x3 rotationY( float _angle )
    {
        float sinA = sin(_angle);
        float cosA = cos(_angle);
        return Mat3x3( cosA, 0.0f, sinA,
                       0.0f, 1.0f, 0.0f,
                      -sinA, 0.0f, cosA);
    }

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space around z-axis.
    inline Mat3x3 rotationZ( float _angle )
    {
        float sinA = sin(_angle);
        float cosA = cos(_angle);
        return Mat3x3(cosA, -sinA, 0.0f,
                      sinA,  cosA, 0.0f,
                      0.0f,  0.0f, 1.0f);
    }

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space from 3 angles:
    ///     rotationZ(_z) * rotationY(_y) * rotationX(_x).
    inline Mat3x3 rotation( float _x, float _y, float _z )
    {
        float sinA = sin(_z);
        float cosA = cos(_z);
        float sinB = sin(_y);
        float cosB = cos(_y);
        float sinC = sin(_x);
        float cosC = cos(_x);
        return Mat3x3(cosA * cosB, cosA * sinB * sinC - sinA * cosC, cosA * sinB * cosC + sinA * sinC,
                      sinA * cosB, sinA * sinB * sinC + cosA * cosC, sinA * sinB * cosC - cosA * sinC,
                     -sinB,        cosB * sinC,                      cosB * cosC                     );
    }

    inline Mat3x3 rotation( const Vec3& _eulerAngles )         { return rotation(_eulerAngles.x, _eulerAngles.y, _eulerAngles.z); }

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space for an arbitrary axis.
    inline Mat3x3 rotation( const Vec3& _axis, float _angle )
    {
        float sinA = sin(_angle);
        float cosA = cos(_angle);
        float iCosA = 1.0f - cosA;
        return Mat3x3(_axis.x * _axis.x * iCosA + cosA,           _axis.x * _axis.y * iCosA - _axis.z * sinA, _axis.x * _axis.z * iCosA + _axis.y * sinA,
                      _axis.x * _axis.y * iCosA + _axis.z * sinA, _axis.y * _axis.y * iCosA + cosA,           _axis.y * _axis.z * iCosA - _axis.x * sinA,
                      _axis.x * _axis.z * iCosA - _axis.y * sinA, _axis.y * _axis.z * iCosA + _axis.x * sinA, _axis.z * _axis.z * iCosA + cosA       );
    }

    // ********************************************************************* //
    /// \brief Rotation matrix from one direction into another.
    inline Mat3x3 rotation( const Vec3& _from, const Vec3& _to ) // TESTED
    {
        // Get lengths for normalization
        eiAssert(approx(len(_from), 1.0f), "Expected a normalized direction vector '_from'.");
        eiAssert(approx(len(_to), 1.0f), "Expected a normalized direction vector '_from'.");
        //float lf = len(_from);
        //float lt = len(_to);
        Vec3 axis = cross(_from, _to);
        // Compute sin(alpha) from cross product lf * lt * sin(alpha) and normalize
        float sinA = len(axis);
        if(sinA) axis /= sinA;
        //sinA /= lf * lt;
        float cosA = sqrt((1.0f - sinA) * (1.0f + sinA));
        // Create axis-angle matrix
        float iCosA = 1.0f - cosA;
        return Mat3x3(axis.x * axis.x * iCosA + cosA,          axis.x * axis.y * iCosA - axis.z * sinA, axis.x * axis.z * iCosA + axis.y * sinA,
                      axis.x * axis.y * iCosA + axis.z * sinA, axis.y * axis.y * iCosA + cosA,          axis.y * axis.z * iCosA - axis.x * sinA,
                      axis.x * axis.z * iCosA - axis.y * sinA, axis.y * axis.z * iCosA + axis.x * sinA, axis.z * axis.z * iCosA + cosA         );
    }

    // ********************************************************************* //
    /// \brief Rotation matrix from quaternion.
    inline Mat3x3 rotation( const Quaternion& _quaternion )
    {
        return Mat3x3(_quaternion);
    }

    // ********************************************************************* //
    /// \brief Reflection matrix for a plane through the origin.
    /// \details Use this to construct mirror matrices. To simply reflect a
    ///     vector use the reflect() method (faster). Beware: _normal and _at
    ///     in the two methods are orthogonal.
    /// \param [in] _normal Normal of the reflecting plane (at the origin).
    ///     The normal must not be normalized.
    inline Mat3x3 housholder( const Vec3& _normal ) // TESTED
    {
        float norm = 2.0f / dot(_normal, _normal);
        float nx = norm * _normal.x;
        float ny = norm * _normal.y;
        float nz = norm * _normal.z;
        float nxy = nx * _normal.y;
        float nxz = nx * _normal.z;
        float nyz = ny * _normal.z;
        return Mat3x3(1.0f - nx * _normal.x,      - nxy,                 - nxz,
                           - nxy,            1.0f - ny * _normal.y,      - nyz,
                           - nxz,                 - nyz,            1.0f - nz * _normal.z);
    }

    // ********************************************************************* //
    /// \brief Reflect a vector at a plane.
    /// \details Reflects using I - 2 * dot(I,At) * At (householder transformation)
    /// \param [in] _incident A direction or position vector which should be
    ///     reflected.
    /// \param [in] _at The normal vector for the reflection plane (normalized!).
    template<typename T, unsigned N>
    inline Vec<T,N> reflect( const Vec<T,N>& _incident, const Vec<T,N>& _at )
    {
        eiAssertWeak(approx(lensq(_at), 1.0f), "The reflection normal must be normalized!");
        return _incident - (static_cast<T>(2) * dot(_incident, _at)) * _at;
    }

    template<typename T, unsigned N>
    inline RVec<T,N> reflect( const RVec<T,N>& _incident, const RVec<T,N>& _at )
    {
        eiAssertWeak(approx(lensq(_at), 1.0f), "The reflection normal must be normalized!");
        return _incident - (static_cast<T>(2) * dot(_incident, _at)) * _at;
    }

    // ********************************************************************* //
    /// \brief Create a matrix in 3D space where the target is on
    ///     the positive z-axis.
    /// \details This method creates an left-hand system (LHS) with positive
    ///     z-axis.
    ///
    ///     This matrix always 'looks' from the origin. Use a camera() or
    ///     combine with a translation() to add an other origin.
    /// \param [in] _target A position which should lie on the z-axis.
    /// \param [in] _up The x-axis/horizon is created perpendicular to this
    ///     vector. The up vector must not necessarily be normalized.
    inline Mat3x3 lookAt( const Vec3& _target, const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f))
    {
        Vec3 zAxis = normalize(_target);
        Vec3 xAxis = normalize(cross(_up, zAxis));
        Vec3 yAxis = cross(zAxis, xAxis);
        return axis( xAxis, yAxis, zAxis );
    }

    // ********************************************************************* //
    /// \brief Create a camera matrix in homogeneous space.
    /// \details The camera matrix must be in homogeneous space due to the
    ///    translation.
    /// \details This method creates an left-hand system (LHS) with positive
    ///    z-axis.
    inline Mat4x4 camera( const Vec3& _position,
                          const Vec3& _target,
                          const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f) )
    {
        return Mat4x4(lookAt( _target - _position, _up )) * translation( -_position );
    }

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
    inline Mat4x4 perspectiveGL( float _l, float _r, float _b, float _t, float _n, float _f )
    {
        return Mat4x4(2.0f*_n / (_r-_l), 0.0f,              (_l+_r) / (_l-_r),  0.0f,
                      0.0f,              2.0f*_n / (_t-_b), (_b+_t) / (_b-_t),  0.0f,
                      0.0f,              0.0f,              (_f+_n) / (_f-_n), -2.0f*_n*_f / (_f-_n),
                      0.0f,              0.0f,              1.0f,               0.0f);
    }

    // ********************************************************************* //
    /// \brief Create OpenGL perspective projection matrix from fovY and aspect ratio.
    /// \details Assumes an LHS coordinate system.
    ///
    ///    The OpenGL frustum is defined in the [-1,-1,-1] x [1,1,1] cube.
    /// \param [in] _fovY Field of view in the y direction, in radians.
    /// \param [in] _aspectRatio width/height of the frame buffer.
    inline Mat4x4 perspectiveGL( float _fovY, float _aspectRatio, float _near, float _far )
    {
        // cot(x) == tan(π/2 - x)
        float h = tan(PI * 0.5f -_fovY / 2.0f);
        float w = h / _aspectRatio;
        return Mat4x4(w,    0.0f, 0.0f,              0.0f,
                      0.0f, h,    0.0f,              0.0f,
                      0.0f, 0.0f, (_far+_near) / (_far-_near), -2.0f*_near*_far / (_far-_near),
                      0.0f, 0.0f, 1.0f,              0.0f);
    }

    // ********************************************************************* //
    /// \brief Create inverse OpenGL perspective projection matrix from fovY and aspect ratio.
    /// \details This is the inverse matrix of perspectiveGL computed analytical.
    /// \param [in] _fovY Field of view in the y direction, in radians.
    /// \param [in] _aspectRatio width/height of the frame buffer.
    inline Mat4x4 inversePerspectiveGL( float _fovY, float _aspectRatio, float _near, float _far )
    {
        float h = tan(PI * 0.5f -_fovY / 2.0f);
        float w = h / _aspectRatio;
        return Mat4x4(1.0f / w, 0.0f, 0.0f, 0.0f,
                      0.0f, 1.0f / h, 0.0f, 0.0f,
                      0.0f, 0.0f, 0.0f, 1.0f,
                      0.0f, 0.0f, (_far-_near) / (-2.0f*_near*_far), (_far+_near) / (2.0f*_near*_far) );
    }

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
    inline Mat4x4 orthographicGL( float _l, float _r, float _b, float _t, float _n, float _f )
    {
        return Mat4x4(2.0f / (_r-_l), 0.0f, 0.0f, -(_r+_l) / (_r-_l),
                      0.0f, 2.0f / (_t-_b), 0.0f, -(_t+_b) / (_t-_b),
                      0.0f, 0.0f, 2.0f / (_f-_n), -(_f+_n) / (_f-_n),
                      0.0f, 0.0f, 0.0f,           1.0f);
    }

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
    inline Mat4x4 perspectiveDX( float _l, float _r, float _b, float _t, float _n, float _f )
    {
        return Mat4x4(2.0f*_n / (_r-_l), 0.0f,              (_l+_r) / (_l-_r), 0.0f,
                      0.0f,              2.0f*_n / (_t-_b), (_b+_t) / (_b-_t), 0.0f,
                      0.0f,              0.0f,              _f / (_f-_n),      -_n*_f / (_f-_n),
                      0.0f,              0.0f,              1.0f,              0.0f);
    }

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
    inline Mat4x4 perspectiveDX( float _fovY, float _aspectRatio, float _near, float _far )
    {
        // cot(x) == tan(π/2 - x)
        float h = tan(PI * 0.5f -_fovY / 2.0f);
        float w = h / _aspectRatio;
        return Mat4x4(w,    0.0f, 0.0f,         0.0f,
                      0.0f, h,    0.0f,         0.0f,
                      0.0f, 0.0f, _far / (_far-_near), -_near*_far/(_far-_near),
                      0.0f, 0.0f, 1.0f,         0.0f);
    }

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
    inline Mat4x4 orthographicDX( float _l, float _r, float _b, float _t, float _n, float _f )
    {
        return Mat4x4(2.0f / (_r-_l), 0.0f, 0.0f, -(_r+_l) / (_r-_l),
                      0.0f, 2.0f / (_t-_b), 0.0f, -(_t+_b) / (_t-_b),
                      0.0f, 0.0f, 1.0f / (_f-_n), -(_f+_n) / (_f-_n),
                      0.0f, 0.0f, 0.0f,           1.0f);
    }


    // ********************************************************************* //
    // Algebra                                                               //
    // ********************************************************************* //

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
    inline bool decomposeLUp(const Matrix<T,N,N>& _A,
                             Matrix<T,N,N>& _LU,
                             Vec<uint,N>& _p) // TESTED
    {
        // LUP decomposition algorithm by Cormen et al. with both matrices L and U
        // combined to one.

        // Copy and work inplace
        _LU = _A;

        // Init permutation to none
        for(uint i = 0; i < N; ++i) _p[i] = i;
        // Loop over rows in A
        for(uint r = 0; r < N-1; ++r)
        {
            // Search a non zero pivot element in current column
            T pivot = static_cast<T>(0);
            uint p = 0;
            for(uint i = r; i < N; ++i) if(abs(_LU(i,r)) > pivot) {
                pivot = abs(_LU(i,r));
                p = i;
            }
            if(pivot == static_cast<T>(0)) return false;
            // Swap the lines
            if(p != r)
            {
                // Permutation vector
                uint tmpUI = _p[r]; _p[r] = _p[p]; _p[p] = tmpUI;
                // Complete row: the first r-1 elements are already part of U
                for(uint c = 0; c < N; ++c) {
                    T tmpT = _LU(r,c); _LU(r,c) = _LU(p,c); _LU(p,c) = tmpT;
                }
            }
            // Gauß elimination but keep the diagonal element
            for(uint i = r+1; i < N; ++i)
            {
                _LU(i,r) /= _LU(r,r);
                for(uint j = r+1; j < N; ++j)
                    _LU(i,j) -= _LU(i,r) * _LU(r,j);
            }
        }
        return _LU(N-1,N-1) != 0.0f;
    }

    // ********************************************************************* //
    /// \brief Compute one or more solutions for linear equation system given
    ///     the LUp decomposition.
    /// \details Computation time: O(n²)
    /// \param [in] _LU Combined LU matrix from decomposeLUp.
    /// \param [in] _p Permutation vector from decomposeLUp.
    /// \param [in] _B Set of N column vectors which contain the target vectors.
    /// \return Solution X of LU X = B[p].
    template<typename T, unsigned M, unsigned N>
    inline Matrix<T,M,N> solveLUp(const Matrix<T,M,M>& _LU,
                                  const Matrix<uint,M,1>& _p,
                                  const Matrix<T,M,N>& _B) // TESTED
    {
        Matrix<T,M,N> X;
        for(uint n = 0; n < N; ++n)
        {
            // Compute L Y = P B
            for(uint i = 0; i < M; ++i)
            {
                T sum = static_cast<T>(0);
                for(uint j = 0; j < i; ++j)
                    sum += _LU(i,j) * X(j,n);
                X(i,n) = _B(_p[i],n) - sum;
            }
            // Compute U X = Y
            for(int i = M-1; i >= 0; --i)
            {
                T sum = static_cast<T>(0);
                for(uint j = i+1; j < M; ++j)
                    sum += _LU(i,j) * X(j,n);
                X(i,n) = (X(i,n) - sum) / _LU(i,i);
            }
        }
        return X;
    }

    // ********************************************************************* //
    /// \brief Compute a spectral decomposition Q^T * diag(d) * Q of 2x2 and 3x3
    ///     symmetric matrices A.
    /// \param [in] _sort Sort eigenvalues and vectors descending.
    ///
    ///     The results for 2x2 matrices are always sorted. For 3x3 matices
    ///     sorting takes extra time and can be disabled if not needed.
    /// \param [out] _Q The orthonormal basis where rows are the eigenvectors
    ///     corresponding to the eigenvalues _lambda of A.
    /// \param [out] _lambda Eigenvalues of A.
    /// \return Number of iterations (50 is the maximum used internally) or -1
    ///     if no solution can be found (complex eigenvalues).
    template<typename T>
    inline int decomposeQl(const Matrix<T,2,2>& _A, Matrix<T,2,2>& _Q, Vec<T,2>& _lambda)
    {
        T p = -_A[0] - _A[3];
        T q = _A[0] * _A[3] - _A[1] * _A[2];
        T discriminant = p*p - static_cast<T>(4)*q;
        if(discriminant < static_cast<T>(0)) return -1;
        T dsqrt = sqrt(discriminant);
        // Numerically stable solution for both roots
        // The roots are guaranteed to be sorted descending
        if(p > 0) {
            _lambda.x = -static_cast<T>(2)*q/(p + dsqrt);
            _lambda.y = (-p - dsqrt)/static_cast<T>(2);
        } else {
            _lambda.x = (-p + dsqrt)/static_cast<T>(2);
            _lambda.y = -static_cast<T>(2)*q/(p - dsqrt);
        }
        // Compute normalized eigenvectors
        if(_A[1] != static_cast<T>(0)) {
            _Q(0) = normalize(RVec<T,2>(static_cast<T>(1), (_lambda.x - _A[0]) / _A[1]));
            _Q(1) = normalize(RVec<T,2>(static_cast<T>(1), (_lambda.y - _A[0]) / _A[1]));
        } else if(_A[2] != static_cast<T>(0)) {
            _Q(0) = normalize(RVec<T,2>(static_cast<T>(1), (_lambda.x - _A[3]) / _A[2]));
            _Q(1) = normalize(RVec<T,2>(static_cast<T>(1), (_lambda.y - _A[3]) / _A[2]));
        } else {
            // Fallback identity if input was a diagonal matrix
            _Q[0] = static_cast<T>(1); _Q[1] = static_cast<T>(0);
            _Q[2] = static_cast<T>(0); _Q[3] = static_cast<T>(1);
        }
        return 1;
    }

    /// \brief Iterative spectral decomposition for 3x3 matrices.
    // Implementation from http://www.melax.com/diag.html
    // Other can be found on http://stackoverflow.com/questions/4372224/fast-method-for-computing-3x3-symmetric-matrix-spectral-decomposition
    template<typename T>
    int decomposeQl(const Matrix<T,3,3>& _A, Matrix<T,3,3>& _Q, Vec<T,3>& _lambda, bool _sort = true) // TESTED
    {
        int i = 0;
        Quaternion q = qidentity();// TODO type T
        while(i < 50)
        {
            ++i;
            _Q = Mat3x3(q);
            Vec<T,3> offdiag; // elements (1,2), (0,2) and (0,1) of the diagonal matrix D = Q * A * Q^T
            Matrix<T,3,3> D = _Q * _A * transpose(_Q);
            _lambda.x = D(0,0); _lambda.y = D(1,1); _lambda.z = D(2,2);
            offdiag.x = D(1,2); offdiag.y = D(0,2); offdiag.z = D(0,1);
            // Find index of largest element of offdiag
            Vec<T,3> absod = abs(offdiag);
            int k0 = (absod.x > absod.y && absod.x > absod.z) ? 0 : (absod.y > absod.z) ? 1 : 2;
            int k1 = (k0+1)%3;
            int k2 = (k0+2)%3;
            if(offdiag[k0] == static_cast<T>(0)) break;  // Diagonal matrix converged
            float t = (_lambda[k2] - _lambda[k1]) / (static_cast<T>(2) * offdiag[k0]);
            float sgnt = (t > static_cast<T>(0)) ? static_cast<T>(1) : static_cast<T>(-1);
            t = sgnt / (t*sgnt + sqrt(t*t + static_cast<T>(1)));
            float c = sqrt(t*t + static_cast<T>(1));
            if(c == static_cast<T>(1)) break;           // reached numeric limit
            c = static_cast<T>(1) / c;
            // Create a jacobi rotation for this iteration
            Quaternion jr;
            jr.z[k0] = sgnt * sqrt((1.0f - c) / 2.0f);
            jr.z[k1] = jr.z[k2] = 0.0f;
            jr.r = sqrt((1.0f - jr.z[k0]) * (1.0f + jr.z[k0]));
            if(jr.r == 1.0f) break;                     // another numeric limit
            q = normalize(jr * q);
        }
        _Q = Mat3x3(q);
        _lambda.x = _Q(0,0)*dot(_A(0), _Q(0)) + _Q(0,1)*dot(_A(1), _Q(0)) + _Q(0,2)*dot(_A(2), _Q(0));
        _lambda.y = _Q(1,0)*dot(_A(0), _Q(1)) + _Q(1,1)*dot(_A(1), _Q(1)) + _Q(1,2)*dot(_A(2), _Q(1));
        _lambda.z = _Q(2,0)*dot(_A(0), _Q(2)) + _Q(2,1)*dot(_A(1), _Q(2)) + _Q(2,2)*dot(_A(2), _Q(2));

        // Sort (Network 0,2 0,1 1,2)
        if(_sort)
        {
            if(_lambda.x < _lambda.z)
            {
                float ts = _lambda.x; _lambda.x = _lambda.z; _lambda.z = ts;
                Matrix<T,1,3> tv = _Q(0); _Q(0) = _Q(2); _Q(2) = tv;
            }
            if(_lambda.x < _lambda.y)
            {
                float ts = _lambda.x; _lambda.x = _lambda.y; _lambda.y = ts;
                Matrix<T,1,3> tv = _Q(0); _Q(0) = _Q(1); _Q(1) = tv;
            }
            if(_lambda.y < _lambda.z)
            {
                float ts = _lambda.y; _lambda.y = _lambda.z; _lambda.z = ts;
                Matrix<T,1,3> tv = _Q(1); _Q(1) = _Q(2); _Q(2) = tv;
            }
        }

        return i;
    }

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
    /*{
        // "Random" initialization
        for(int i = 0; i < N; ++i)
            _eigenvec[i] = 1.0f / (i+1);

        // Rayleigh quotient iteration
        // Problem: does not guarantee to converge to the largest eigenvalue
        /*_eigenvec = normalize(_eigenvec);
        _eigenval = dot(_eigenvec, _A * _eigenvec);
        int i = 0;
        for(; i < 30; ++i)     // 30 is the maximum iteration number
        {
            Matrix<T,N,N> AshiftLU;
            Vec<uint,N> Ashiftp;
            if(!decomposeLUp(_A - diag(Vec<T,N>(_eigenval)), AshiftLU, Ashiftp))
                break;
            Vec<T, N> b1 = solveLUp(AshiftLU, Ashiftp, _eigenvec);
            b1 = normalize(b1);
            _eigenval = dot(_eigenvec, _A * _eigenvec);
            bool converged = approx(b1, _eigenvec);
            _eigenvec = b1;
            if(converged) break;
        }
        _eigenval = dot(_eigenvec, _A * _eigenvec);
        return i+1;//* /

        // Power iteration
        int i;
        for(i = 0; i < 300; ++i)     // 300 is the maximum iteration number
        {
            Vec<T,N> b1 = normalize(_A * _eigenvec);
            // Custom modification: 1.5x the step-width made into the correct
            // direction.
            //b1 = normalize(b1 - reflect(_eigenvec, b1));
            //b1 = -reflect(_eigenvec, b1);
            bool converged = approx(b1, _eigenvec);
            _eigenvec = b1;
            if(converged) break;
        }
        // Compute the eigenvalue with the Rayleigh quotient (x^T A x) / (x^T x).
        // Since the vector is normalized the denominator is not required.
        _eigenval = dot(_eigenvec, _A * _eigenvec);//* /
        return i;
    }*/

    // ********************************************************************* //
    /// \brief Compute the Cholesky decomposition L*L^T of symmetric positive
    ///     definite matrices A.
    /// \param [out] _L A lower triangular matrix such that L*L^T = A.
    /// \return true if the decomposition was possible (A is symmetric positive
    ///     definite).
    ///
    ///     The code will not explicitly test if the matrix is symmetric. It
    ///     simply assumes it. I.e. you might get a return value of true
    ///     even for non symmetric matrices.
    template<typename T, unsigned N>
    inline bool decomposeCholesky(const Matrix<T,N,N>& _A, Matrix<T,N,N>& _L) // TESTED
    {
        for(uint y = 0; y < N; ++y)
        {
            for(uint x = 0; x <= y; ++x)
            {
                float sum = _A(y,x);
                for(uint i = 0; i < x; ++i)
                    sum -= _L(y,i) * _L(x,i);
                // Handle diagonal elements else than the lower ones
                if(y > x)
                    _L(y,x) = sum / _L(x,x);
                else if(sum > 0.0f)
                    _L(x,x) = sqrt(sum);
                // Not positive definite!
                else return false;
            }
            // Zero out the remaining elements
            for(uint x = y + 1; x < N; ++x)
                _L(y,x) = static_cast<T>(0);
        }
        // Fill the upper triangular matrix with 0
        return true;
    }

    // Some toy specializations: It is unsure if the compiler may reach the same
    // code with the loop implementation only. Surly, the specializations are faster
    // in debug mode.
    template<typename T>
    inline bool decomposeCholesky(const Matrix<T,2,2>& _A, Matrix<T,2,2>& _L) // TESTED
    {
        if(_A[0] <= 0.0f) return false;
        _L[0] = sqrt(_A[0]);
        _L[1] = 0.0f;
        _L[2] = _A[1] / _L[0];
        float l10sq = _A[3] - _L[2] * _L[2];
        _L[3] = sqrt(l10sq);
        return l10sq > 0.0f;
    }

    template<typename T>
    inline bool decomposeCholesky(const Matrix<T,3,3>& _A, Matrix<T,3,3>& _L) // TESTED
    {
        if(_A[0] <= 0.0f) return false;
        _L[0] = sqrt(_A[0]);
        _L[1] = 0.0f;
        _L[2] = 0.0f;
        _L[3] = _A[3] / _L[0];
        float tmp = _A[4] - _L[3] * _L[3];
        if(tmp <= 0.0f) return false;
        _L[4] = sqrt(tmp);
        _L[5] = 0.0f;
        _L[6] = _A[6] / _L[0];
        _L[7] = (_A[7] - _L[6] * _L[3]) / _L[4];
        tmp = _A[8] - _L[6] * _L[6] - _L[7] * _L[7];
        if(tmp <= 0.0f) return false;
        _L[8] = sqrt(tmp);
        return true;
    }

    // ********************************************************************* //
    /// \brief Invert a quadratic matrix.
    /// \details If the matrix has no inverse the identity is returned.
    /// \return Inverse matrix or identity.
    template<typename T, unsigned N>
    Matrix<T,N,N> invert(const Matrix<T,N,N>& _mat0) // TESTED
    {
        Matrix<T,N,N> LU;
        Vec<uint32,N> p;
        if( decomposeLUp( _mat0, LU, p ) )
            return solveLUp( LU, p, identity<T,N>() );
        else return identity<T,N>();
    }

    // ********************************************************************* //
    /// \brief Compute the determinant of a matrix.
    /// \details This uses fixed implementations for N=2 and N=3 and LU
    ///     decomposition for N > 3.
    template<typename T, unsigned N>
    inline T determinant(const Matrix<T,N,N>& _A) // TESTED
    {
        Matrix<T,N,N> LU;
        Vec<uint32,N> p;
        if( decomposeLUp( _A, LU, p ) )
        {
            // The diagonal product of L gives the value
            T det = -LU[0];
            // The number of permutations the sign (#p even -> positive)
            if( p[0] != 0 ) det = -det;
            for(int i = 1; i < N; ++i) {
                det *= LU[i + N * i];
                if( p[i] != i ) det = -det;
            }
            return det;
        }
        return static_cast<T>(0);
    }

    template<typename T>
    inline T determinant(const Matrix<T,2,2>& _A)
    {
        return _A[0]*_A[3] - _A[1]*_A[2];
    }

    template<typename T>
    inline T determinant(const Matrix<T,3,3>& _A)
    {
        return _A[0]*_A[4]*_A[8] + _A[1]*_A[5]*_A[6] + _A[2]*_A[3]*_A[7]
              -_A[2]*_A[4]*_A[6] - _A[1]*_A[3]*_A[8] - _A[0]*_A[5]*_A[7];
    }

}

// Remove helper macros.
#undef RESULT_TYPE
#undef ENABLE_IF
