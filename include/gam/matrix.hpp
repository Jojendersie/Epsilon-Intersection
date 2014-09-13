#pragma once

#include <type_traits>

#include "elementarytypes.hpp"
#include "details/matrixcomponents.hpp"

namespace gam {

	/// \brief MxN row-major matrix class.
	/// \details The matrix is the basis for all matrix and vector types. It
	///     supports all kinds of matrix <-> matrix and matrix <-> scalar
	///     operations including adding, ... a value to all components. If not
	///     stated else each operator works component wise.
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
#		define RESULT_TYPE(op) typename std::enable_if<					\
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
		Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5);
		template<ENABLE_IF(N * M == 8)>
		Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7);
		template<ENABLE_IF(N * M == 9)>
		Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7, T _s8);
		template<ENABLE_IF(N * M == 12)>
		Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7, T _s8, T _s9, T _s10, T _s11);
		template<ENABLE_IF(N * M == 16)>
		Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7, T _s8, T _s9, T _s10, T _s11, T _s12, T _s13, T _s14, T _s15);

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
		Matrix<RESULT_TYPE(+), M, N> operator+ (const Matrix<T1,M,N>& _mat1);  // TESTED
		/// \brief Subtract two matrices component wise.
		/// \details Subtraction is not commutative.
		template<typename T1>
		Matrix<RESULT_TYPE(-), M, N> operator- (const Matrix<T1,M,N>& _mat1);  // TESTED
		/// \brief Unary minus on all components.
		Matrix<T, M, N> operator- ();                                          // TESTED
		/// \brief Add a scalar value to all components.

		/// \brief Matrix multiplication.
		/// \details Matrix multiplication is not commutative.
		/// \returns Matrix product with dimensions MxO = MxN * NxO. The result
		///    is a scalar if M = N = 1.
		template<typename T1, uint O>
		typename std::conditional<M * O == 1, RESULT_TYPE(*), Matrix<RESULT_TYPE(*), M, O>>::type
		operator* (const Matrix<T1,N,O>& _mat1);                               // TESTED
		/// \brief Component wise multiplication for vectors of the same size.
		template<typename T1, ENABLE_IF(N == 1)>
		Matrix<RESULT_TYPE(*), M, 1> operator* (const Matrix<T1,M,1>& _mat1);  // TESTED
		template<typename T1, ENABLE_IF(M == 1)>
		Matrix<RESULT_TYPE(*), 1, N> operator* (const Matrix<T1,1,N>& _mat1);  // TESTED
		/// \brief Component wise division for vectors of the same size.
		template<typename T1, ENABLE_IF(N == 1)>
		Matrix<RESULT_TYPE(/), M, 1> operator/ (const Matrix<T1,M,1>& _mat1);  // TESTED
		template<typename T1, ENABLE_IF(M == 1)>
		Matrix<RESULT_TYPE(/), 1, N> operator/ (const Matrix<T1,1,N>& _mat1);  // TESTED

		/// \brief Compare if two matrices are identical (using elementary !=).
		bool operator== (const Matrix<T,M,N>& _mat1);                          // TESTED
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
	//								 FUNCTIONS								 //
	// ********************************************************************* //

	// ********************************************************************* //
	/// \brief Check if the absolute difference between all elements is smaller
	///    than epsilon.
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





	// Include implementation.
#	include "details/matrix.inl"

	// Remove helper macros.
#	undef RESULT_TYPE
#   undef ENABLE_IF
}
