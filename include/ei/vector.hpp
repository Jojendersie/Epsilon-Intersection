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
#       define RESULT_TYPE(op) std::enable_if_t<                   \
            !std::is_base_of<details::NonScalarType, T1>::value && \
            !std::is_base_of<details::NonScalarType, T>::value,    \
            decltype(std::declval<T>() op std::declval<T1>())      \
        >

        // Enable a function on a condition via template list.
        // This macro allows conditional compilation even for methods without
        // parameters and return value.
        // Therefore if must be inserted in the template list and `class` at
        // the same position in the implementation.
#       define ENABLE_IF(condition) typename = std::enable_if_t< (condition) >

        /// \brief Construction without initialization. The values are undefined!
        EIAPI Matrix() noexcept = default;

        /// \brief Convert a matrix/vector with a different elementary type.
        template<typename T1>
        EIAPI explicit Matrix(const Matrix<T1,M,N>& _mat1) noexcept // TESTED
        {
            for(uint i = 0; i < N * M; ++i)
                this->m_data[i] = static_cast<T>(_mat1[i]);
        }

        /// \brief Forward to base constructors
        template<typename T1, ENABLE_IF((!std::is_base_of<details::NonScalarType, T1>::value))>
        constexpr EIAPI explicit Matrix(T1 _a0) noexcept :
            details::Components<T,M,N>(_a0)
        {}
        template<typename T1, typename T2, typename... Args>
        constexpr EIAPI Matrix(T1 _a0, T2 _a1, Args... _args) noexcept :
            details::Components<T,M,N>(_a0, _a1, std::forward<Args>(_args)...)
        {}

        /// \brief Allow explicit truncation of the dimension sizes.
        template<typename T1, uint M1, uint N1, ENABLE_IF((M < M1 && N <= N1) || (M <= M1 && N < N1))>
        EIAPI explicit Matrix(const Matrix<T1,M1,N1>& _mat1, uint _rowOff = 0, uint _colOff = 0) noexcept
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
        template<typename T1, uint N1, ENABLE_IF((M == 1 && N > N1))>
        EIAPI explicit Matrix(const Matrix<T1,1,N1>& _mat1) noexcept // TESTED
        {
            for(unsigned i = 0; i < N1; ++i)
                this->m_data[i] = static_cast<T>(_mat1[i]);
            for(unsigned i = N1; i < N; ++i)
                this->m_data[i] = static_cast<T>(1);
        }
        template<typename T1, uint M1, ENABLE_IF((N == 1 && M > M1))>
        EIAPI explicit Matrix(const Matrix<T1,M1,1>& _mat1) noexcept // TESTED
        {
            for(unsigned i = 0; i < M1; ++i)
                this->m_data[i] = static_cast<T>(_mat1[i]);
            for(unsigned i = M1; i < M; ++i)
                this->m_data[i] = static_cast<T>(1);
        }
        template<typename T1, uint M1, uint N1, ENABLE_IF((M > M1 && N >= N1 && N > 1) || (M >= M1 && N > N1 && M > 1))>
        EIAPI explicit Matrix(const Matrix<T1,M1,N1>& _mat1) noexcept // TESTED
        {
                // Indices for _mat1 and result
                unsigned i = 0, j = 0;
                for(unsigned y = 0; y < M1; ++y)
                {
                    // Copy MxN part
                    for(unsigned x = 0; x < N1; ++x)
                        this->m_data[j++] = static_cast<T>(_mat1[i++]);
                    // New elements at the end of the row is 0
                    for(unsigned x = N1; x < N; ++x)
                        this->m_data[j++] = static_cast<T>(0);
                }
                // Fill new rows
                for(unsigned y = M1; y < M; ++y)
                    for(unsigned x = 0; x < N; ++x)
                        this->m_data[j++] = x == y ? static_cast<T>(1) : static_cast<T>(0);
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
        EIAPI T& operator () (uint _row, uint _col) noexcept // TESTED
        {
            eiAssertWeak(_row < M && _col < N, "Index out of bounds!");
            return this->m_data[_row * N + _col];
        }
        EIAPI T operator () (uint _row, uint _col) const noexcept // TESTED
        {
            eiAssertWeak(_row < M && _col < N, "Index out of bounds!");
            return this->m_data[_row * N + _col];
        }

        /// \brief Single row access
        /// \param [in] _row Index of the row in [0,M-1].
        EIAPI Matrix<T,1,N>& operator () (uint _row) noexcept // TESTED
        {
            eiAssertWeak(_row < M, "Index out of bounds!");
            return reinterpret_cast<Matrix<T,1,N>&>(this->m_data[_row * N]);
        }
        EIAPI const Matrix<T,1,N>& operator () (uint _row) const noexcept // TESTED
        {
            eiAssertWeak(_row < M, "Index out of bounds!");
            return reinterpret_cast<const Matrix<T,1,N>&>(this->m_data[_row * N]);
        }

        /// \brief Access an element by a single index treating the matrix as 1D.
        /// \param [in] _index Index in the range [0, N * M - 1].
        /// \returns Reference with read or write access to the element
        ///    depending on the constness of the matrix.
        EIAPI T& operator [] (uint _index) noexcept // TESTED
        {
            eiAssertWeak(_index < N * M, "Index out of bounds!");
            return this->m_data[_index];
        }
        EIAPI T operator [] (uint _index) const noexcept // TESTED
        {
            eiAssertWeak(_index < N * M, "Index out of bounds!");
            return this->m_data[_index];
        }

        /// \brief Get access to a subrange [FROM, TO) in the vector.
        /// \tparam FROM First element in the output range (inclusive).
        /// \tparam TO Exclusive right boundary.
        template<uint FROM, uint TO>//, ENABLE_IF((N == 1) && (FROM < TO) && (TO <= M))>
        EIAPI Matrix<T, TO - FROM, 1>& subcol() noexcept // TESTED
        {
            static_assert(N == 1, "THIS must be a column vector.");
            static_assert((FROM >= 0) && (FROM < TO) && (TO <= M), "Invalid parameter range.");
            return *reinterpret_cast<Matrix<T, TO - FROM, 1>*>(this->m_data + FROM);
        }
        template<uint FROM, uint TO>//, ENABLE_IF((N == 1) && (FROM < TO) && (TO <= M))>
        EIAPI const Matrix<T, TO - FROM, 1>& subcol() const noexcept
        {
            static_assert(N == 1, "THIS must be a column vector.");
            static_assert((FROM >= 0) && (FROM < TO) && (TO <= M), "Invalid parameter range.");
            return *reinterpret_cast<const Matrix<T, TO - FROM, 1>*>(this->m_data + FROM);
        }
        template<uint FROM, uint TO>//, ENABLE_IF((M == 1) && (FROM < TO) && (TO <= N))>
        EIAPI Matrix<T, 1, TO - FROM>& subrow() noexcept // TESTED
        {
            static_assert(M == 1, "THIS must be a row vector.");
            static_assert((FROM >= 0) && (FROM < TO) && (TO <= N), "Invalid parameter range.");
            return *reinterpret_cast<Matrix<T, 1, TO - FROM>*>(this->m_data + FROM);
        }
        template<uint FROM, uint TO>//, ENABLE_IF((M == 1) && (FROM < TO) && (TO <= N))>
        EIAPI const Matrix<T, 1, TO - FROM>& subrow() const noexcept
        {
            static_assert(M == 1, "THIS must be a row vector.");
            static_assert((FROM >= 0) && (FROM < TO) && (TO <= N), "Invalid parameter range.");
            return *reinterpret_cast<const Matrix<T, 1, TO - FROM>*>(this->m_data + FROM);
        }
        template<uint R_FROM, uint R_TO, uint C_FROM, uint C_TO>
        EIAPI Matrix<T, R_TO - R_FROM, C_TO - C_FROM> submat() const noexcept // TESTED
        {
            static_assert((R_FROM >= 0) && (R_FROM < R_TO) && (R_TO <= M), "Invalid parameter range for rows.");
            static_assert((C_FROM >= 0) && (C_FROM < C_TO) && (C_TO <= N), "Invalid parameter range for columns.");
            Matrix<T, R_TO - R_FROM, C_TO - C_FROM> result;
            uint readIdx = R_FROM * N + C_FROM;
            uint writeIdx = 0;
            for(uint r = R_FROM; r < R_TO; ++r)
            {
                for(uint c = C_FROM; c < C_TO; ++c)
                {
                    result.m_data[writeIdx] = this->m_data[readIdx];
                    ++writeIdx;
                    ++readIdx;
                }
                readIdx += N - (C_TO - C_FROM);
            }
            return result;
        }

        /// \brief Add two matrices component wise.
        /// \details Addition is commutative.
        template<typename T1>
        EIAPI Matrix<RESULT_TYPE(+), M, N> operator + (const Matrix<T1,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_OP(+)
        /// \brief Subtract two matrices component wise.
        /// \details Subtraction is not commutative.
        template<typename T1>
        EIAPI Matrix<RESULT_TYPE(-), M, N> operator - (const Matrix<T1,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_OP(-)

        /// \brief Unary minus on all components.
        template<typename T1 = T, ENABLE_IF(std::is_signed<T1>::value)>
        EIAPI Matrix<T, M, N> operator - () const noexcept // TESTED
            EI_CODE_GEN_MAT_UNARY_OP(-)
        /// \brief Component wise binary not.
        template<typename T1 = T, ENABLE_IF(std::is_integral<T1>::value)>
        EIAPI Matrix<T, M, N> operator ~ () const noexcept // TESTED
            EI_CODE_GEN_MAT_UNARY_OP(~)

        /// \brief Matrix multiplication.
        /// \details Matrix multiplication is not commutative.
        /// \returns Matrix product with dimensions MxO = MxN * NxO. The result
        ///    is a scalar if M = N = 1.
        template<typename T1, uint O>
        EIAPI std::conditional_t<M * O == 1, RESULT_TYPE(*), Matrix<RESULT_TYPE(*), M, O>>
        operator * (const Matrix<T1,N,O>& _mat1) const noexcept // TESTED
        {
			std::conditional_t<M * O == 1, RESULT_TYPE(*), Matrix<RESULT_TYPE(*), M, O>> result{};
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
        EIAPI Matrix<RESULT_TYPE(*), M, N> operator * (const Matrix<T1,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_OP(*)

        /// \brief Component wise division for vectors of the same size.
        template<typename T1, ENABLE_IF((M != N) && sizeof(T1))>
        EIAPI Matrix<RESULT_TYPE(/), M, N> operator / (const Matrix<T1,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_OP(/)

        /// \brief Component wise binary or.
        /// \details Or is commutative.
        template<typename T1>
        EIAPI Matrix<RESULT_TYPE(|), M, N> operator | (const Matrix<T1,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_OP(|)
        /// \brief Component wise binary and.
        /// \details And is commutative.
        template<typename T1>
        EIAPI Matrix<RESULT_TYPE(&), M, N> operator & (const Matrix<T1,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_OP(&)
        /// \brief Component wise binary xor.
        /// \details Xor is commutative.
        template<typename T1>
        EIAPI Matrix<RESULT_TYPE(^), M, N> operator ^ (const Matrix<T1,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_OP(^)
        /// \brief Component wise modulo (rest of integer division).
        template<typename T1>
        EIAPI Matrix<RESULT_TYPE(%), M, N> operator % (const Matrix<T1,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_OP(%)

        /// \brief Self assigning component wise addition.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator += (const Matrix<T1,M,N>& _mat1) noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(+=)
        /// \brief Self assigning component wise subtraction.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator -= (const Matrix<T1,M,N>& _mat1) noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(-=)
        /// \brief Self assigning component wise multiplication for vectors
        ///    of the same size. Matrix multiplication in case of squared matrices!
        template<typename T1, ENABLE_IF((M != N) && sizeof(T1))>
        EIAPI Matrix<T, M, N>& operator *= (const Matrix<T1,M,N>& _mat1) noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(*=)
        template<typename T1>
        EIAPI Matrix<T, N, N>& operator *= (const Matrix<T1,M,M>& _mat1) noexcept // TESTED
        {
            // Use matrix multiplication
            *this = (*this) * _mat1;
            return *this;
        }
        /// \brief Self assigning component wise division for vectors
        ///    of the same size.
        // TODO: Matrix division (mul with inverse) for squared matrices??
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator /= (const Matrix<T1,M,N>& _mat1) noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(/=)
        /// \brief Self assigning component wise binary or.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator |= (const Matrix<T1,M,N>& _mat1) noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(|=)
        /// \brief Self assigning component wise binary and.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator &= (const Matrix<T1,M,N>& _mat1) noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(&=)
        /// \brief Self assigning component wise binary xor.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator ^= (const Matrix<T1,M,N>& _mat1) noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(^=)
        /// \brief Self assigning modulo (rest of integer division)
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator %= (const Matrix<T1,M,N>& _mat1) noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_SEFL_OP(%=)

        /// \brief Self assigning scalar addition.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator += (T1 _s) noexcept // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(+=)
        /// \brief Self assigning scalar subtraction.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator -= (T1 _s) noexcept // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(-=)
        /// \brief Self assigning scalar multiplication.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator *= (T1 _s) noexcept // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(*=)
        /// \brief Self assigning scalar division.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator /= (T1 _s) noexcept // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(/=)
        /// \brief Self assigning scalar or.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator |= (T1 _s) noexcept // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(|=)
        /// \brief Self assigning scalar and.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator &= (T1 _s) noexcept // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(&=)
        /// \brief Self assigning scalar xor.
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator ^= (T1 _s) noexcept // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(^=)
        /// \brief Self assigning scalar modulo (rest of integer division).
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator %= (T1 _s) noexcept // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(%=)
        /// \brief Self assigning component wise shift
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator >>= (T1 _s) noexcept // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(>>=)
        template<typename T1>
        EIAPI Matrix<T, M, N>& operator <<= (T1 _s) noexcept // TESTED
            EI_CODE_GEN_MAT_SCALAR_SEFL_OP(<<=)

        /// \brief Compare component wise, if two matrices are identical.
        EIAPI bool operator == (const Matrix<T,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_BOOL_ALL_OP(==)
        /// \brief Compare component wise, if two matrices are distinct.
        EIAPI bool operator != (const Matrix<T,M,N>& _mat1) const noexcept // TESTED
        {
            for(uint i = 0; i < N * M; ++i)
                if((*this)[i] != _mat1[i]) return true;
                    return false;
        }
        /// \brief Compare component wise, if elements are smaller or equal.
        EIAPI bool operator <= (const Matrix<T,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_BOOL_ALL_OP(<=)
        /// \brief Compare component wise, if elements are smaller.
        EIAPI bool operator < (const Matrix<T,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_BOOL_ALL_OP(<)
        /// \brief Compare component wise, if elements are greater.
        EIAPI bool operator > (const Matrix<T,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_BOOL_ALL_OP(>)
        /// \brief Compare component wise, if elements are greater or equal.
        EIAPI bool operator >= (const Matrix<T,M,N>& _mat1) const noexcept // TESTED
            EI_CODE_GEN_MAT_MAT_BOOL_ALL_OP(>=)
    };


    // ********************************************************************* //
    // Scalar operators

    /// \brief Add a scalar to all components.
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(+), M, N> operator + (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(+)
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<RESULT_TYPE(+), M, N> operator + (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(+)
    /// \brief Subtract a scalar from all components.
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(-), M, N> operator - (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(-)
    /// \brief Subtract all components from a scalar.
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<RESULT_TYPE(-), M, N> operator - (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(-)
    /// \brief Multiply a scalar to all components.
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(*), M, N> operator * (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(*)
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<RESULT_TYPE(*), M, N> operator * (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(*)
    /// \brief Divide all components by a scalar.
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(/), M, N> operator / (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(/)
    /// \brief Divide a scalar by each component.
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<RESULT_TYPE(/), M, N> operator / (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(/)

    /// \brief Binary or of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(|), M, N> operator | (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(|)
    /// \brief Binary or of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<RESULT_TYPE(|), M, N> operator | (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(|)
    /// \brief Binary and of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(&), M, N> operator & (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(&)
    /// \brief Binary and of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<RESULT_TYPE(&), M, N> operator & (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(&)
    /// \brief Binary xor of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(^), M, N> operator ^ (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(^)
    /// \brief Binary xor of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<RESULT_TYPE(^), M, N> operator ^ (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(^)
    /// \brief Modulo of all components and a scalar.
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(%), M, N> operator % (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(%)
    /// \brief Modulo of a scalar and all components.
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<RESULT_TYPE(%), M, N> operator % (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_OP(%)
    /// \brief Component wise shift.
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(>>), M, N> operator >> (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(>>)
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(>>), M, N> operator >> (T1 _s, const Matrix<T,M,N>& _mat) noexcept
        EI_CODE_GEN_SCALAR_MAT_OP(>>)
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(<<), M, N> operator << (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_OP(<<)
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<RESULT_TYPE(<<), M, N> operator << (T1 _s, const Matrix<T,M,N>& _mat) noexcept
        EI_CODE_GEN_SCALAR_MAT_OP(<<)

    /// \brief Test all components with respect to a scalar.
    template<typename T, uint M, uint N, typename T1>
    EIAPI bool operator == (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_ALL_OP(==)
    template<typename T1, typename T, uint M, uint N>
    EIAPI bool operator == (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_ALL_OP(==)
    template<typename T, uint M, uint N, typename T1>
    EIAPI bool operator != (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            if(_mat[i] != _s) return true;
        return false;
    }
    template<typename T1, typename T, uint M, uint N>
    EIAPI bool operator != (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            if(_s != _mat[i]) return true;
        return false;
    }
    template<typename T, uint M, uint N, typename T1>
    EIAPI bool operator < (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_ALL_OP(<)
    template<typename T1, typename T, uint M, uint N>
    EIAPI bool operator < (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_ALL_OP(<)
    template<typename T, uint M, uint N, typename T1>
    EIAPI bool operator <= (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_ALL_OP(<=)
    template<typename T1, typename T, uint M, uint N>
    EIAPI bool operator <= (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_ALL_OP(<=)
    template<typename T, uint M, uint N, typename T1>
    EIAPI bool operator >= (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_ALL_OP(>=)
    template<typename T1, typename T, uint M, uint N>
    EIAPI bool operator >= (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_ALL_OP(>=)
    template<typename T, uint M, uint N, typename T1>
    EIAPI bool operator > (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_ALL_OP(>)
    template<typename T1, typename T, uint M, uint N>
    EIAPI bool operator > (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_ALL_OP(>)


    /// \brief Compare component wise, if two matrices are identical.
    template<typename T, uint M, uint N>
    EIAPI Matrix<bool,M,N> equal (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) noexcept // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(==)
    /// \brief Compare component wise, if two matrices are distinct.
    template<typename T, uint M, uint N>
    EIAPI Matrix<bool,M,N> neq (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) noexcept // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(!=)
    /// \brief Compare component wise, if elements are smaller or equal.
    template<typename T, uint M, uint N>
    EIAPI Matrix<bool,M,N> lesseq (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) noexcept // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(<=)
    /// \brief Compare component wise, if elements are smaller.
    template<typename T, uint M, uint N>
    EIAPI Matrix<bool,M,N> less (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) noexcept // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(<)
    /// \brief Compare component wise, if elements are greater.
    template<typename T, uint M, uint N>
    EIAPI Matrix<bool,M,N> greater (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) noexcept // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(>)
    /// \brief Compare component wise, if elements are greater or equal.
    template<typename T, uint M, uint N>
    EIAPI Matrix<bool,M,N> greatereq (const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1) noexcept // TESTED
        EI_CODE_GEN_MAT_MAT_BOOL_OP(>=)

    /// \brief Test if all components compare equal to a scalar.
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<bool, M, N> equal (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(==)
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<bool, M, N> equal (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(==)
    /// \brief Test if any component is non equal to a scalar.
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<bool, M, N> neq (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(!=)
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<bool, M, N> neq (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(!=)
    template<typename T, uint M, uint N, typename T1>
    /// \brief Test if all components compare to a scalar.
    EIAPI Matrix<bool, M, N> less (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(<)
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<bool, M, N> less (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(<)
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<bool, M, N> lesseq (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(<=)
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<bool, M, N> lesseq (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(<=)
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<bool, M, N> greatereq (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(>=)
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<bool, M, N> greatereq (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(>=)
    template<typename T, uint M, uint N, typename T1>
    EIAPI Matrix<bool, M, N> greater (const Matrix<T,M,N>& _mat, T1 _s) noexcept // TESTED
        EI_CODE_GEN_MAT_SCALAR_BOOL_OP(>)
    template<typename T1, typename T, uint M, uint N>
    EIAPI Matrix<bool, M, N> greater (T1 _s, const Matrix<T,M,N>& _mat) noexcept // TESTED
        EI_CODE_GEN_SCALAR_MAT_BOOL_OP(>)

// Remove helper macros from vectordetailsA.hpp.
#undef EI_CODE_GEN_MAT_MAT_OP
#undef EI_CODE_GEN_MAT_UNARY_OP
#undef EI_CODE_GEN_MAT_MAT_SEFL_OP
#undef EI_CODE_GEN_MAT_SCALAR_SEFL_OP
#undef EI_CODE_GEN_MAT_MAT_BOOL_OP
#undef EI_CODE_GEN_MAT_SCALAR_OP
#undef EI_CODE_GEN_SCALAR_MAT_OP
#undef EI_CODE_GEN_MAT_SCALAR_BOOL_OP
#undef EI_CODE_GEN_SCALAR_MAT_BOOL_OP


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
    EIAPI bool approx(const Matrix<T,M,N>& _mat0,
        const Matrix<T,M,N>& _mat1,
        T _epsilon = T(1e-6)) noexcept  // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            // if(!(abs(_mat1[i] - _mat0[i]) <= _epsilon)) return false;
            if(!approx(_mat1[i], _mat0[i], _epsilon))
                return false;
        return true;
    }

    // ********************************************************************* //
    /// \brief Computes the sum of component wise products.
    /// \returns Scalar value of the sum of component products.
    template<typename T, unsigned M, unsigned N, typename T1>
    constexpr EIAPI RESULT_TYPE(*) dot(const Matrix<T,M,N>& _mat0,
        const Matrix<T1,M,N>& _mat1) noexcept // TESTED
    {
        RESULT_TYPE(*) sum = _mat0[0] * _mat1[0];
        for(uint i = 1; i < N * M; ++i)
            sum += _mat0[i] * _mat1[i];
        return sum;
    }

    // ********************************************************************* //
    /// \brief Computes the cross product of two 3d vectors (RHS).
    /// \returns Perpendicular vector with length |v0|·|v1|·sin(∡(v0,v1)).
    template<typename T, typename T1, unsigned M, unsigned N, ENABLE_IF((N==1 && M==3) || (N==3 && M==1))>
    constexpr EIAPI Matrix<RESULT_TYPE(*),M,N> cross(const Matrix<T,M,N>& _v0,
        const Matrix<T1,M,N>& _v1) noexcept
    {
        return Matrix<RESULT_TYPE(*),M,N>{_v0.y * _v1.z - _v0.z * _v1.y,
            _v0.z * _v1.x - _v0.x * _v1.z,
            _v0.x * _v1.y - _v0.y * _v1.x};
    }

    // ********************************************************************* //
    /// \brief Computes the cross product of two 2d vectors.
    /// \returns The determinant of the 2x2 matrix: v0.x·v1.y - v0.y·v1.x.
    template<typename T, typename T1, unsigned M, unsigned N, ENABLE_IF((N==1 && M==2) || (N==2 && M==1))>
    constexpr EIAPI RESULT_TYPE(*) cross(const Matrix<T,M,N>& _v0,
        const Matrix<T1,M,N>& _v1) noexcept
    {
        return _v0.x * _v1.y - _v0.y * _v1.x;
    }

    // ********************************************************************* //
    /// \brief Computes the sum of squared components.
    /// \details This is equivalent to dot(_mat0, _mat0).
    /// \returns Squared euclidean length (scalar).
    template<typename T>
    constexpr EIAPI auto lensq(const T& _elem0) noexcept -> decltype(dot(_elem0, _elem0)) // TESTED
    {
        return dot(_elem0, _elem0);
    }

    // ********************************************************************* //
    /// \brief Computes the root of the sum of squared components.
    /// \details This is the euclidean length for vectors and Quaternions and
    ///    the Frobenius norm for matrices.
    /// \returns Euclidean length (scalar).
    template<typename T>
    constexpr EIAPI auto len(const T& _elem0) noexcept -> decltype(std::sqrt(dot(_elem0, _elem0))) // TESTED
    {
        return sqrt(dot(_elem0, _elem0));
    }

    // ********************************************************************* //
    /// \brief Normalizes a vector, quaternion or matrix with respect to len.
    /// \details This is equivalent to elem0 / len(_elem0).
    /// \returns Normalized vector or matrix.
    template<typename T>
    EIAPI T normalize(const T& _mat0) noexcept // TESTED
    {
        return _mat0 / len(_mat0);
    }

    // ********************************************************************* //
    /// \brief Component wise maximum.
    /// \returns A matrix with the maximum values from both inputs.
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,M,N> max(const Matrix<T,M,N>& _mat0,
        const Matrix<T,M,N>& _mat1) noexcept // TESTED
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
    EIAPI Matrix<T,M,N> min(const Matrix<T,M,N>& _mat0,
                            const Matrix<T,M,N>& _mat1) noexcept // TESTED
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
    EIAPI T max(const Matrix<T,M,N>& _mat0) noexcept // TESTED
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
    EIAPI T min(const Matrix<T,M,N>& _mat0) noexcept // TESTED
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
    EIAPI Matrix<T,M,N> clamp(const Matrix<T,M,N>& _mat,
        const Matrix<T,M,N>& _min,
        const Matrix<T,M,N>& _max) noexcept // TESTED
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
    EIAPI Matrix<T,M,N> clamp(const Matrix<T,M,N>& _mat,
        T _min,
        T _max) noexcept // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = clamp(_mat[i], _min, _max);
        return result;
    }

    /// \brief Clamp all components to [0,1]
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,M,N> saturate(const Matrix<T,M,N>& _mat) noexcept // TESTED
    {
        return clamp(_mat, static_cast<T>(0), static_cast<T>(1));
    }

    // ********************************************************************* //
    /// \brief Round all components towards negative infinity
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<typename details::Int<sizeof(T)>::stype,M,N> floor(const Matrix<T,M,N>& _mat) noexcept // TESTED
    {
        Matrix<typename details::Int<sizeof(T)>::stype,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = floor(_mat[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Round all components towards negative infinity
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<typename details::Int<sizeof(T)>::stype,M,N> ceil(const Matrix<T,M,N>& _mat) noexcept // TESTED
    {
        Matrix<typename details::Int<sizeof(T)>::stype,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = ceil(_mat[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Round all components towards next number (x.5 rounds up)
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<typename details::Int<sizeof(T)>::stype,M,N> round(const Matrix<T,M,N>& _mat) noexcept // TESTED
    {
        Matrix<typename details::Int<sizeof(T)>::stype,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = round(_mat[i]);
        return result;
    }

     // ********************************************************************* //
    /// \brief Get the fraction in (-1,1) using elementary frac().
    /// \param _x [in] The number to be splitted.
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,M,N> frac(const Matrix<T,M,N>& _x) noexcept
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = frac(_x[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Divide a number into the integer and fractional part.
    /// \param _x [in] The number to be splitted.
    /// \param _int [out] The integer part of the number (rounded to zero).
    /// \returns The fraction of the number in (-1,1).
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,M,N> intfrac(const Matrix<T,M,N>& _x, Matrix<Sint<sizeof(T)>,M,N>& _int) noexcept
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = intfrac(_x[i], _int[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Divide a number into an integer and a positive fractional part.
    /// \details This method uses f-floor(f) instead of f-int(f) and therefore
    ///      has a continous behavior around zero
    /// \param _x [in] The number to be splitted.
    /// \param _int [out] The integer part of the number.
    /// \returns The fraction of the number in [0,1).
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,M,N> floorfrac(const Matrix<T,M,N>& _x, Matrix<Sint<sizeof(T)>,M,N>& _int) noexcept
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = floorfrac(_x[i], _int[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Get the smallest positive number m, such that x=y*c+m with c
    ///     in Z, for each component.
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,M,N> mod(const Matrix<T,M,N>& _x, T _y) noexcept // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = mod(_x[i], _y);
        return result;
    }

    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,M,N> mod(const Matrix<T,M,N>& _x, const Matrix<T,M,N>& _y) noexcept
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < N * M; ++i)
            result[i] = mod(_x[i], _y[i]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Compute the square root for each component
    template<typename T, unsigned M, unsigned N, ENABLE_IF((N==1) || (M==1))>
    EIAPI Matrix<T,M,N> sqrt(const Matrix<T,M,N>& _v0) noexcept // TESTED
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
    EIAPI Matrix<T,M,N> pow(const Matrix<T,M,N>& _v0, float _exponent) noexcept // TESTED
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
    /// \brief Compute the power for each component
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,M,N> exp(const Matrix<T,M,N>& _v0) noexcept
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < M * N; ++i)
            result[i] = std::exp(_v0[i]);
        return result;
    }

    // Unfortunately, the above declaration hides the elementary one -> make it
    // visible again
    using std::exp;

    // ********************************************************************* //
    /// \brief Element wise natural logarithm for matrices (basis e).
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,M,N> log(const Matrix<T,M,N>& _mat0) noexcept // TESTED
    {
        Matrix<T,M,N> result;
        for(uint i = 0; i < M * N; ++i)
            result[i] = std::log(_mat0[i]);
        return result;
    }

    /// \brief Element wise logarithm for matrices (basis 2).
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,M,N> log2(const Matrix<T,M,N>& _mat0) noexcept // TESTED
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
    EIAPI decltype(std::declval<T>() + std::declval<T>()) sum(const Matrix<T,M,N>& _mat0) noexcept // TESTED
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
    EIAPI T prod(const Matrix<T,M,N>& _mat0) noexcept // TESTED
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
    EIAPI T avg(const Matrix<T,M,N>& _mat0) noexcept // TESTED
    {
        return sum(_mat0) / T(M * N);
    }

    // ********************************************************************* //
    /// \brief Absolute values for all components.
    /// \returns Matrix with component wise absolute values.
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,M,N> abs(const Matrix<T,M,N>& _mat0) noexcept // TESTED
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
    EIAPI Matrix<T,M,N> sign(const Matrix<T,M,N>& _mat0) noexcept // TESTED
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
    EIAPI Matrix<T,M,N> sgn(const Matrix<T,M,N>& _mat0) noexcept // TESTED
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
    EIAPI Matrix<decltype(std::declval<T0>() * std::declval<T1>()),M,N>
        bilerp(Matrix<T0,M,N> _x00, Matrix<T0,M,N> _x01,
               Matrix<T0,M,N> _x10, Matrix<T0,M,N> _x11,
               T1 _t0, T1 _t1) noexcept // TESTED
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
    EIAPI auto slerp(const Matrix<T0,M,N>& _v0, const Matrix<T0,M,N>& _v1, T1 _t) noexcept -> decltype(_v0*_t)
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

    // ********************************************************************* //
    /// \brief Test if at least one element of the matrix is true.
    /// \return false, if all elements off the matrix are false.
    template<unsigned M, unsigned N>
    EIAPI bool any(const Matrix<bool,M,N>& _mat0) noexcept // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            if(_mat0[i]) return true;
        return false;
    }

    // ********************************************************************* //
    /// \brief Test if no element of the matrix is true.
    /// \return true, if all elements off the matrix are false.
    template<unsigned M, unsigned N>
    EIAPI bool none(const Matrix<bool,M,N>& _mat0) noexcept // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            if(_mat0[i]) return false;
        return true;
    }

    // ********************************************************************* //
    /// \brief Test if all elements of the matrix are true.
    /// \return true, if all elements off the matrix are true.
    template<unsigned M, unsigned N>
    EIAPI bool all(const Matrix<bool,M,N>& _mat0) noexcept // TESTED
    {
        for(uint i = 0; i < N * M; ++i)
            if(!_mat0[i]) return false;
        return true;
    }


    // ********************************************************************* //
    /// \brief Transpose a matrix or vector (switch the dimensions).
    template<typename T, unsigned M, unsigned N>
    EIAPI Matrix<T,N,M> transpose(const Matrix<T,M,N>& _mat0) noexcept // TESTED
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
    EIAPI bool orthonormalize(Matrix<T,M,N>& _mat0) noexcept // TESTED
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

    namespace details {
        // Recursive helper for orthonormalization
        // Recursion end
        template<typename TVec0>
        EIAPI void removeProjectedPart(const TVec0&)
        {
        }
        template<typename TVec0, typename TVec1, typename... TVecs>
        EIAPI void removeProjectedPart(const TVec0& _vec0, TVec1& _vec1, TVecs&... _vecs)
        {
            _vec1 -= dot(_vec0, _vec1) * _vec0;
            removeProjectedPart(_vec0, _vecs...);
        }
    }

    /// \brief Gram-Schmidt orthonormalization for a list of vectors.
    /// \details The first vector will only be normalized, the others orthogonalized on top.
    template<typename TVec0, typename... TVecs>
    EIAPI bool orthonormalize(TVec0& _vec0, TVecs&... _vecs) noexcept // TESTED
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
    /// \brief Generate the N x N diagonal matrix.
    /// \param [in] _v0 A vector with the diagonal entries.
    template<typename T, unsigned N>
    EIAPI Matrix<T,N,N> diag( const Vec<T,N>& _v0 ) noexcept // TESTED
    {
        Matrix<T,N,N> result { T(0) };
        for(uint n = 0; n < N; ++n)
            result[n * N + n] = _v0[n];
        return result;
    }

    // ********************************************************************* //
    /// \brief Generate the N x N identity matrix.
    template<typename T, unsigned N>
    EIAPI Matrix<T,N,N> identity() noexcept // TESTED
    {
        return diag(Vec<T,N>{1});
    }

    /// \brief Alias for identity<float,2>().
    constexpr EIAPI Mat2x2 identity2x2() noexcept    { return Mat2x2{1.0f, 0.0f, 0.0f, 1.0f}; }
    /// \brief Alias for identity<float,3>().
    constexpr EIAPI Mat3x3 identity3x3() noexcept    { return Mat3x3{1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f}; } // TESTED
    /// \brief Alias for identity<float,4>().
    constexpr EIAPI Mat4x4 identity4x4() noexcept    { return Mat4x4{1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f}; } // TESTED

    // ********************************************************************* //
    /// \brief Convert a vector from Cartesian coordinates in spherical
    ///     (angular) coordinates.
    /// \details In 2d this gives (r, α) and in 3d this gives (r, θ, ϕ) for
    ///     the vectors: (r cos α, r sin α) and (r cos θ, r sin θ cos ϕ, r sin θ sin ϕ).
    ///     This is not the usual convention, but unifies dimensions!
    /// \return The N-1 spherical angles and the length of the vector.
    ///     (r, φ1, φ02, ..., φN-1) where
    ///     φ1, ..., φN-2 ∈ [0,π) and φN-1 ∈ [0,2π)
    template<typename T, unsigned M, unsigned N, ENABLE_IF((M==1) || (N==1))>
    EIAPI Matrix<T,M,N> sphericalCoords( const Matrix<T,M,N>& _v0 ) noexcept // TESTED
    {
        static_assert(N*M >= 2, "In 1D cartesian and spherical coordinates are the same!");
        Matrix<T,M,N> result;
        // Accumulate the squared length over the iterations
        result[0] = sq(_v0[N*M-1]) + sq(_v0[N*M-2]);
        result[N*M-1] = atan2(_v0[N*M-1], _v0[N*M-2]);
        if(result[N*M-1] < 0.0f) result[N*M-1] += 2.0f * PI;
        // if( _v0[N-1] < 0.0f ) result[N-1] = 2.0f * PI - result[N-1];
        for(uint i = 2; i < N*M; ++i)
        {
            result[0] += sq(_v0[N*M-i-1]);
            result[N*M-i] = acos(clamp(_v0[N*M-i-1]/sqrt(result[0]), -1.0f, 1.0f));
        }
        result[0] = sqrt(result[0]);
        return result;
    }

    // ********************************************************************* //
    /// \brief Convert a vector from spherical coordinates (r, φ1, φ02, ..., φN-1)
    ///     to regular Cartesian coordinates.
    /// \return The regular Cartesian vector.
    template<typename T, unsigned M, unsigned N, ENABLE_IF((M==1) || (N==1))>
    EIAPI Matrix<T,M,N> cartesianCoords( const Matrix<T,M,N>& _v0 ) noexcept // TESTED
    {
        eiAssertWeak(_v0[0] > 0.0f, "Expected the length to be greater 0!");
        Matrix<T,M,N> result;
        float tmp = _v0[0];
        for(uint i = 0; i < M*N-1; ++i)
        {
            result[i] = tmp * cos(_v0[i+1]);
            tmp *= sin(_v0[i+1]);
        }
        result[M*N-1] = tmp;
        return result;
    }

    // ********************************************************************* //
    /// \brief Apply transformations in homogeneous space. This includes a
    ///     division by w after the transformation
    template<typename T, unsigned N>
    EIAPI Vec<T,N> transformDiv( const Vec<T, N>& _what,
        const Matrix<T, N+1, N+1>& _space ) noexcept
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
        Vec<T,N> result;
        for(uint i = 0; i < N; ++i)
            result[i] = t[i] / t[N];
        return result;
    }

    /// \brief Apply transformations in 3x4/4x4 space (rotation + translation).
    ///     This does NOT include a division by w.
    template<typename T, unsigned N>
    EIAPI Vec<T,N> transform( const Vec<T, N>& _what,
        const Matrix<T, N+1, N+1>& _space ) noexcept // TESTED
    {
        Vec<T,N> result;
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
    EIAPI Vec<T,N> transform( const Vec<T, N>& _what,
        const Matrix<T, N, N+1>& _space ) noexcept // TESTED
    {
        Vec<T,N> result;
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

    /// \brief Apply transformations and ignore the fourth component and the
    ///     translation.
    /// \details Use this function to transform direction vectors. Still, there
    ///     is no normalization involved and the direction might be scaled by
    ///     the matrix.
    template<typename T, unsigned N>
    EIAPI Vec<T,N> transformDir( const Vec<T, N>& _what,
        const Matrix<T, N+1, N+1>& _space ) noexcept // TESTED
    {
        Vec<T,N> result;
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
    EIAPI Vec<T,N> transformDir( const Vec<T, N>& _what,
        const Matrix<T, N, N+1>& _space ) noexcept // TESTED
    {
        Vec<T,N> result;
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

    /// \brief Apply transformations with a matrix multiplication.
    template<typename T, unsigned M>
    EIAPI Vec<T,M> transform( const Vec<T,M>& _what,
        const Matrix<T,M,M>& _space ) noexcept
    {
        return _space * _what;
    }

    // ********************************************************************* //
    /// \brief Create a translation matrix in homogeneous coordinate space.
    /// \param [in] _vector Translate by/Add this vector.
    /// \details The translation matrix always has a dimension one large then
    ///    the vectors.
    ///    To transform a vector append 1 and multiply it from right:
    ///    translation() * VecX(v,1)
    template<typename T, unsigned N>
    constexpr EIAPI Matrix<T,N+1,N+1> translation( const Vec<T, N>& _vector ) noexcept
    {
        Matrix<T,N+1,N+1> result = identity<T,N+1>();
        for(uint i = 0; i < N; ++i)
            result[i * (N+1) + N] = _vector[i];
        return result;
    }

    // ********************************************************************* //
    /// \brief Create a scaling/diagonal matrix from vector.
    template<typename T, unsigned N>
    constexpr EIAPI Matrix<T,N,N> scaling( const Vec<T, N>& _scale ) noexcept
    {
        Matrix<T,N,N> result { T(0) };
        for(uint n = 0; n < N; ++n)
            result[n * N + n] = _scale[n];
        return result;
    }

    // ********************************************************************* //
    /// \brief Create a uniform scaling/diagonal matrix from scalar.
    template<typename T, unsigned N>
    constexpr EIAPI Matrix<T,N,N> scaling( T _scale ) noexcept
    {
        Matrix<T,N,N> result { T(0) };
        for(uint n = 0; n < N; ++n)
            result[n * N + n] = _scale;
        return result;
    }

    // ********************************************************************* //
    /// \brief Use vectors which span a space to build a matrix.
    /// \details The vectors become the rows of the matrix
    // TODO: variadic template variant
    constexpr EIAPI Mat2x2 axis( const Vec2& _x, const Vec2& _y ) noexcept
    {
        return Mat2x2(_x.x, _x.y,
                      _y.x, _y.y);
    }
    constexpr EIAPI Mat3x3 axis( const Vec3& _x, const Vec3& _y, const Vec3& _z ) noexcept // TESTED
    {
        return Mat3x3(_x.x, _x.y, _x.z,
                      _y.x, _y.y, _y.z,
                      _z.x, _z.y, _z.z);
    }
    constexpr EIAPI Mat4x4 axis( const Vec4& _x, const Vec4& _y, const Vec4& _z, const Vec4& _w ) noexcept
    {
        return Mat4x4(_x.x, _x.y, _x.z, _x.w,
                      _y.x, _y.y, _y.z, _y.w,
                      _z.x, _z.y, _z.z, _z.w,
                      _w.x, _w.y, _w.z, _w.w);
    }

    // ********************************************************************* //
    /// \brief Create a vector which is perpendicular to the input one
    template<typename T>
    constexpr EIAPI Vec<T,2u> perpendicular( const Vec<T,2u>& _vector ) noexcept
    {
        return Vec<T,2u> { -_vector.y, _vector.x };
    }
    template<typename T>
    constexpr EIAPI Matrix<T, 1, 2> perpendicular( const Matrix<T, 1, 2>& _vector ) noexcept
    {
        return Matrix<T, 1, 2>{ -_vector.y, _vector.x };
    }

    template<typename T>
    constexpr EIAPI Vec<T, 3> perpendicular( const Vec<T, 3>& _vector ) noexcept
    {
        return abs(_vector.z) < abs(_vector.x) ?
            Vec<T, 3>{-_vector.y, _vector.x, 0} :
            Vec<T, 3>{0, -_vector.z, _vector.y};
    }

    template<typename T>
    constexpr EIAPI Matrix<T, 1, 3> perpendicular( const Matrix<T, 1, 3>& _vector ) noexcept
    {
        return abs(_vector.z) < abs(_vector.x) ?
            Matrix<T, 1, 3>{-_vector.y, _vector.x, 0.0f} :
            Matrix<T, 1, 3>{0.0f, -_vector.z, _vector.y};
    }

    // ********************************************************************* //
    /// \brief Create an orthonormal basis for a single direction vector
    constexpr EIAPI Mat2x2 basis( const Vec2& _vector ) noexcept // TESTED
    {
        return Mat2x2{_vector.x, _vector.y,
                     -_vector.y, _vector.x};
    }

    EIAPI Mat3x3 basis( const Vec3& _vector ) noexcept // TESTED
    {
        eiAssert(approx(len(_vector), 1.0f), "Expected normalized direction vector!");
        Vec3 y;
        if(abs(_vector.z) < abs(_vector.x))
            y = Vec3(-_vector.y, _vector.x, 0.0f) / sqrt(_vector.y*_vector.y + _vector.x*_vector.x);
        else y = Vec3(0.0f, -_vector.z, _vector.y) / sqrt(_vector.y*_vector.y + _vector.z*_vector.z);
        return axis(_vector, y, cross(_vector, y));
    }

    // ********************************************************************* //
    /// \brief Rotation matrix in 2D.
    EIAPI Mat2x2 rotation( float _angle ) noexcept
    {
        float sinA = sin(_angle);
        float cosA = cos(_angle);
        return Mat2x2(cosA, -sinA,
                      sinA,  cosA);
    }

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D space around x-axis.
    EIAPI Mat3x3 rotationX( float _angle ) noexcept
    {
        float sinA = sin(_angle);
        float cosA = cos(_angle);
        return Mat3x3(1.0f, 0.0f,  0.0f,
                      0.0f, cosA, -sinA,
                      0.0f, sinA,  cosA);
    }


    // ********************************************************************* //
    /// \brief Rotation matrix in 3D space around y-axis.
    EIAPI Mat3x3 rotationY( float _angle ) noexcept
    {
        float sinA = sin(_angle);
        float cosA = cos(_angle);
        return Mat3x3( cosA, 0.0f, sinA,
                       0.0f, 1.0f, 0.0f,
                      -sinA, 0.0f, cosA);
    }

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space around z-axis.
    EIAPI Mat3x3 rotationZ( float _angle ) noexcept
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
    EIAPI Mat3x3 rotation( float _x, float _y, float _z ) noexcept
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

    EIAPI Mat3x3 rotation( const Vec3& _eulerAngles )         { return rotation(_eulerAngles.x, _eulerAngles.y, _eulerAngles.z); }

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space for an arbitrary axis.
    EIAPI Mat3x3 rotation( const Vec3& _axis, float _angle ) noexcept
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
    EIAPI Mat3x3 rotation( const Vec3& _from, const Vec3& _to ) noexcept // TESTED
    {
        // Get lengths for normalization
        eiAssert(approx(len(_from), 1.0f), "Expected a normalized direction vector '_from'.");
        eiAssert(approx(len(_to), 1.0f), "Expected a normalized direction vector '_from'.");
        const Vec3 axis = cross(_from, _to);
        // Compute sin(alpha)² from cross product len(_from) * len(_to) * sin(alpha)
        const float lsq = lensq(axis);
        const float cosA = dot(_from, _to);
        // Create axis-angle matrix
        const float iCosA = (1.0f - cosA) / (lsq != 0.0f ? lsq : 1.0f);
        return Mat3x3(axis.x * axis.x * iCosA + cosA,   axis.x * axis.y * iCosA - axis.z, axis.x * axis.z * iCosA + axis.y,
                      axis.x * axis.y * iCosA + axis.z, axis.y * axis.y * iCosA + cosA,   axis.y * axis.z * iCosA - axis.x,
                      axis.x * axis.z * iCosA - axis.y, axis.y * axis.z * iCosA + axis.x, axis.z * axis.z * iCosA + cosA   );
    }

    // ********************************************************************* //
    /// \brief Rotation matrix in 3D/homogeneous space around a ray (arbitrary
    //      position and normalized axis).
    EIAPI Mat3x4 translatedRotation( const Vec3& _origin, const Vec3& _axis, float _angle ) noexcept
    {
        const Mat3x3 rot { rotation(_axis, _angle) };
        // To rotate around an arbitrary position we have to compute
        // translate(o) * rotation() * translate(-o). I.e. first move to origin
        // and later move back. The 3x3 rotation part stays untouched, but there
        // is a rotated translation in the last column.
        const Vec3 translation = _origin + rot * -_origin;
        return Mat3x4 { rot, translation };
    }

    // ********************************************************************* //
    /// \brief Reflection matrix for a plane through the origin.
    /// \details Use this to construct mirror matrices. To simply reflect a
    ///     vector use the reflect() method (faster). Beware: _normal and _at
    ///     in the two methods are orthogonal.
    /// \param [in] _normal Normal of the reflecting plane (at the origin).
    ///     The normal must not be normalized.
    constexpr EIAPI Mat3x3 housholder( const Vec3& _normal ) noexcept // TESTED
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
    template<typename T, unsigned M, unsigned N, ENABLE_IF((M==1) || (N==1))>
    EIAPI Matrix<T,M,N> reflect( const Matrix<T,M,N>& _incident, const Matrix<T,M,N>& _at ) noexcept
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
    EIAPI Mat3x3 lookAt( const Vec3& _target, const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f)) noexcept
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
    EIAPI Mat4x4 camera( const Vec3& _position,
        const Vec3& _target,
        const Vec3& _up = Vec3(0.0f, 1.0f, 0.0f) ) noexcept
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
    constexpr EIAPI Mat4x4 perspectiveGL( float _l, float _r, float _b, float _t, float _n, float _f ) noexcept
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
    EIAPI Mat4x4 perspectiveGL( float _fovY, float _aspectRatio, float _near, float _far ) noexcept
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
    EIAPI Mat4x4 inversePerspectiveGL( float _fovY, float _aspectRatio, float _near, float _far ) noexcept
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
    constexpr EIAPI Mat4x4 orthographicGL( float _l, float _r, float _b, float _t, float _n, float _f ) noexcept
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
    constexpr EIAPI Mat4x4 perspectiveDX( float _l, float _r, float _b, float _t, float _n, float _f ) noexcept
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
    EIAPI Mat4x4 perspectiveDX( float _fovY, float _aspectRatio, float _near, float _far ) noexcept
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
    constexpr EIAPI Mat4x4 orthographicDX( float _l, float _r, float _b, float _t, float _n, float _f ) noexcept
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
    EIAPI bool decomposeLUp(const Matrix<T,N,N>& _A,
                             Matrix<T,N,N>& _LU,
                             Vec<uint,N>& _p) noexcept // TESTED
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
    EIAPI Matrix<T,M,N> solveLUp(const Matrix<T,M,M>& _LU,
                                  const Matrix<uint,M,1>& _p,
                                  const Matrix<T,M,N>& _B) noexcept // TESTED
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
    ///
    ///     The results are always sorted descending for the eigen values.
    /// \param [out] _Q The orthonormal basis where rows are the eigenvectors
    ///     corresponding to the eigenvalues _lambda of A.
    /// \param [out] _lambda Eigenvalues of A.
    /// \return Number of iterations (50 is the maximum used internally) or -1
    ///     if no solution can be found (complex eigenvalues).
    template<typename T>
    EIAPI int decomposeQl(const Matrix<T,2u,2u>& _A, Matrix<T,2u,2u>& _Q, Vec<T,2u>& _lambda) noexcept
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
            // Rows might be swapped if eigen values are sorted inversely.
            if(_A[0] >= _A[3])
            {
                _Q[0] = static_cast<T>(1); _Q[1] = static_cast<T>(0);
                _Q[2] = static_cast<T>(0); _Q[3] = static_cast<T>(1);
            } else
            {
                _Q[0] = static_cast<T>(0); _Q[1] = static_cast<T>(1);
                _Q[2] = static_cast<T>(1); _Q[3] = static_cast<T>(0);
            }
        }
        return 1;
    }

    template<typename T>
    EIAPI int decomposeQl(const Matrix<T,3u,3u>& _A, Matrix<T,3u,3u>& _Q, Vec<T,3u>& _lambda) noexcept // TESTED
    {
        // It follows some substitution magic from https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices.
        // Another useful source is the paper Efficient numerical diagonalization of hermitian 3x3 matrices
        // (see https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/ and https://arxiv.org/pdf/physics/0610206.pdf)
        /*float q = (_A(0,0) + _A(1,1) + _A(2,2)) / 3.0f;
        float p1 = sq(_A(0,1)) + sq(_A(0,2)) + sq(_A(1,2));
        float p2 = sq(_A(0,0) - q) + sq(_A(1,1) - q) + sq(_A(2,2) - q) + 2.0f * p1;
        float p = sqrt(p2 / 6.0f);
        Mat3x3 B = (_A - q * identity3x3()) / p;
        float det = determinant(B) / 2.0f;
        float phi = acos(clamp(det, -1.0f, 1.0f)) / 3.0f;
        // Eigen values satisfy _lambda.x >= _lambda.y >= _lambda.z
        _lambda.x = q + 2 * p * cos(phi);
        _lambda.z = q + 2 * p * cos(phi + 2*PI/3.0f);
        _lambda.y = 3 * q - _lambda.x - _lambda.z;*/
        // The above code suffers from numerical issues. Using doubles helps a bit:
        double q = (double(_A[0]) + double(_A[4]) + double(_A[8])) / 3.0;
        double p1 = sq(double(_A[1])) + sq(double(_A[2])) + sq(double(_A[5]));
        // Check for diagonal matrix
        if(p1 <= q * 1.0e-015) // Something between 4 and 5 ULP error threshold
        {
            // Find the order of the eigenvalues to set lambda and Q in sorted order.
            // Sort (Network 0,2 0,1 1,2)
            int i0 = 0, i1 = 1, i2 = 2;
            if(_A[0] < _A[8]) {i0 = 2; i2 = 0;}
            if(_A[i0*4] < _A[4]) {i1 = i0; i0 = 1;}
            if(_A[i1*4] < _A[i2*4]) {int i = i2; i2 = i1; i1 = i;}
            _lambda[i0] = _A[0];
            _lambda[i1] = _A[4];
            _lambda[i2] = _A[8];
            _Q(i0) = RVec<T,3>(1.0, 0.0, 0.0);
            _Q(i1) = RVec<T,3>(0.0, 1.0, 0.0);
            _Q(i2) = RVec<T,3>(0.0, 0.0, 1.0);
            return 1;
        }
        double p2 = sq(_A[0] - q) + sq(_A[4] - q) + sq(_A[8] - q) + 2.0 * p1;
        double p = sqrt(p2 / 6.0);
        DMat3x3 Bd = (_A - q * identity3x3()) / p;
        //double det = determinant(Bd) / 2.0;
        double det2 = (Bd[4]*(Bd[0]*Bd[8] - Bd[2]*Bd[2]) + Bd[1]*Bd[5]*Bd[2] * 2.0
                     - Bd[1]*Bd[1]*Bd[8] - Bd[0]*Bd[5]*Bd[5]) / 2.0;
        double phi = acos(clamp(det2, -1.0, 1.0)) / 3.0;
        // Eigen values satisfy _lambda.x >= _lambda.y >= _lambda.z
        _lambda.x = float(q + p * 2.0 * cos(phi));
        _lambda.z = float(q + p * 2.0 * cos(phi + 2*3.14159265358979323846/3.0));
        _lambda.y = float(3 * q - _lambda.x - _lambda.z);

        // Compute eigenvectors for the eigenvalues.
        // The two independent columns of A-lambda*I are perpendicular to the eigenvector.
        // Since we do not know which of the columns are independent we compute 2 cross
        // products and use the one with the better condition.
        // There is one problem: if two eigenvalues are equal the vectors are not defined
        // uniquely. This may happen for the first xor the second block below. If it
        // happens for the first we don't now the plane in which to find the vector,
        // which is then defined by the second block -> defer to later.
        bool deferQ0 = false;
        Mat3x3 B = _A - _lambda.x * identity3x3();
        Matrix<float, 1, 3> canditate1 = cross(B(0), B(1));
        Matrix<float, 1, 3> canditate2 = cross(B(0), B(2));
        float l1 = lensq(canditate1);
        float l2 = lensq(canditate2);
        if(l1 + l2 == 0.0f) deferQ0 = true;
        else if(l1 > l2) _Q(0) = canditate1 / sqrt(l1);
        else _Q(0) = canditate2 / sqrt(l2);
        // Repeat for the most different eigenvector (much more stable and can handle
        // two identically eigenvalues)
        B = _A - _lambda.z * identity3x3();
        canditate1 = cross(B(1), B(0));
        canditate2 = cross(B(2), B(0));
        l1 = lensq(canditate1);
        l2 = lensq(canditate2);
        if(l1 + l2 == 0.0f) {
            eiAssert(!deferQ0, "Cannot stably compute eigenvectors.");
            _Q(2) = normalize(perpendicular(_Q(0)));
        }
        else if(l1 > l2) _Q(2) = canditate1 / sqrt(l1);
        else _Q(2) = canditate2 / sqrt(l2);
        if(deferQ0) _Q(0) = normalize(perpendicular(_Q(2)));
        // The third eigenvector is simply the cross product of the other two.
        _Q(1) = cross(_Q(0), _Q(2));

        return 1;
    }

    /// \brief Iterative spectral decomposition for general symmetric matrices.
    // Implementation based on using https://en.wikipedia.org/wiki/Jacobi_rotation
    // Other can be found on http://stackoverflow.com/questions/4372224/fast-method-for-computing-3x3-symmetric-matrix-spectral-decomposition
    // and http://www.melax.com/diag.html
    template<typename T, uint N>
    EIAPI int decomposeQlIter(const Matrix<T,N,N>& _A, Matrix<T,N,N>& _Q, Vec<T,N>& _lambda, bool _sort = true) noexcept // TESTED
    {
        int i = 0;
        _Q = ei::identity<T,N>();
        Matrix<T,N,N> D = _A;	// This matrix is going to be a diagonal matrix over time
        int maxIter = (N * (N - 1)) / 2		// Number of off-diagonal entries
            * (ei::ilog2(N)+1);				// Number of likely sweeps (heuristic)
        for(;i < maxIter; ++i)
        {
            // Find index of largest element of off-side the diagonal
            int k = 0, l = 1;
            for(int row = 0; row < N-1; ++row)
                for(int col = row+1; col < N; ++col)
                    if(ei::abs(D(row,col)) > ei::abs(D(k,l))) { k=row; l=col; }
            // Converged to diagonal matrix?
            if(D(k,l) == static_cast<T>(0)) break;
            // Compute tan(theta) (see wikipedia)
            T beta = (D(l,l) - D(k,k)) / (static_cast<T>(2) * D(k,l));
            T sgnBeta = (beta > static_cast<T>(0)) ? static_cast<T>(1) : static_cast<T>(-1);
            T t = sgnBeta / (beta*sgnBeta + sqrt(beta*beta + static_cast<T>(1)));
            T c = sqrt(t*t + static_cast<T>(1));
            c = static_cast<T>(1) / c;
            T s = c * t;
            // Update matrix entries of Q
            for(int h = 0; h < N; ++h) {
                T qkh = _Q(k,h);
                T qlh = _Q(l,h);
                _Q(k,h) = qkh * c - qlh * s;
                _Q(l,h) = qlh * c + qkh * s;
            }
            // Update matrix entries of D
            T r = s / (static_cast<T>(1) + c);
            for(int h = 0; h < N; ++h) if(h!=k && h!=l) {
                T dhk = D(h,k);
                T dhl = D(h,l);
                D(h,k) = D(k,h) = dhk - s * (dhl + r * dhk);
                D(h,l) = D(l,h) = dhl + s * (dhk - r * dhl);
            }
            D(k,k) = D(k,k) - t * D(k,l);
            D(l,l) = D(l,l) + t * D(k,l);
            D(k,l) = D(l,k) = static_cast<T>(0);
        }
        for(int k = 0; k < N; ++k)
            _lambda[k] = D(k,k);

        if(_sort)
        {
            // Selection sort; TODO use better method.
            for(int k = 0; k < N-1; ++k)
            {
                int m = k;
                for(int l = k+1; l < N; ++l)
                    if(_lambda[l] > _lambda[m]) m = l;
                if(k != m) { // Swap
                    float ts = _lambda[k]; _lambda[k] = _lambda[m]; _lambda[m] = ts;
                    Matrix<T,1,N> tv = _Q(k); _Q(k) = _Q(m); _Q(m) = tv;
                }
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
    EIAPI bool decomposeCholesky(const Matrix<T,N,N>& _A, Matrix<T,N,N>& _L) noexcept // TESTED
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
    EIAPI bool decomposeCholesky(const Matrix<T,2u,2u>& _A, Matrix<T,2u,2u>& _L) noexcept // TESTED
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
    EIAPI bool decomposeCholesky(const Matrix<T,3u,3u>& _A, Matrix<T,3u,3u>& _L) noexcept // TESTED
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
    EIAPI Matrix<T,N,N> invert(const Matrix<T,N,N>& _mat0) noexcept // TESTED
    {
        Matrix<T,N,N> LU;
        Vec<uint32,N> p;
        if( decomposeLUp( _mat0, LU, p ) )
            return solveLUp( LU, p, identity<T,N>() );
        else return identity<T,N>();
    }

    /// \brief Invert a 2x2 matrix.
    /// \details Direct inversion of small matrices is faster than decomposition.
    template<typename T>
    EIAPI Matrix<T,2,2> invert(const Matrix<T,2,2>& _mat) noexcept // TESTED
    {
        float det = _mat[0] * _mat[3] - _mat[1] * _mat[2];
        if(det == 0.0f) return identity<T,2>();
        return Matrix<T,2,2>{
             _mat[3] / det, -_mat[1] / det,
            -_mat[2] / det,  _mat[0] / det
        };
    }

    // ********************************************************************* //
    /// \brief Compute the determinant of a matrix.
    /// \details This uses fixed implementations for N=2 and N=3 and LU
    ///     decomposition for N > 3.
    template<typename T, unsigned N>
    EIAPI T determinant(const Matrix<T,N,N>& _A) noexcept // TESTED
    {
        Matrix<T,N,N> LU;
        Vec<uint32,N> p;
        if( decomposeLUp( _A, LU, p ) )
        {
            // The diagonal product of L gives the value
            T det = -LU[0];
            // The number of permutations the sign (#p even -> positive)
            if( p[0] != 0 ) det = -det;
            for(uint i = 1; i < N; ++i) {
                det *= LU[i + N * i];
                if( p[i] != i ) det = -det;
            }
            return det;
        }
        return static_cast<T>(0);
    }

    template<typename T>
    constexpr EIAPI T determinant(const Matrix<T,2,2>& _A) noexcept
    {
        return _A[0]*_A[3] - _A[1]*_A[2];
    }

    template<typename T>
    constexpr EIAPI T determinant(const Matrix<T,3,3>& _A) noexcept
    {
        return _A[0]*_A[4]*_A[8] + _A[1]*_A[5]*_A[6] + _A[2]*_A[3]*_A[7]
              -_A[2]*_A[4]*_A[6] - _A[1]*_A[3]*_A[8] - _A[0]*_A[5]*_A[7];
    }

}

// Remove helper macros.
#undef RESULT_TYPE
#undef ENABLE_IF
