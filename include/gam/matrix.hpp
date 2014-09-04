#pragma once

#include <type_traits>

#include "gam/elementarytypes.hpp"
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

		// The new variadic templates allow a more generic definition of all
		// those constructors but compiler support is lagging.

		/// \brief Construction without initialization. The values are undefined!
		Matrix();

		/// \brief Construction from N * M scalar values (up to 16 elements).
		/// \details The template meta programming trick allows only the
		///    compilation of the matching constructor.
		template<class = typename std::enable_if<N * M == 2, class Dummy>::type>
		Matrix(T _v0, T _v1);                                                  // TESTED
		template<class = typename std::enable_if<N * M == 3, class Dummy>::type>
		Matrix(T _v0, T _v1, T _v2);                                           // TESTED
		template<class = typename std::enable_if<N * M == 4, class Dummy>::type>
		Matrix(T _v0, T _v1, T _v2, T _v3);                                    // TESTED

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
#		define RESULT_TYPE(Op) typename std::enable_if<					\
			!std::is_base_of<details::MatrixType, T1>::value &&			\
			!std::is_base_of<details::MatrixType, T>::value,			\
			decltype(std::declval<T>() Op std::declval<T1>())           \
		>::type


		/// \brief Add two matrices component wise.
		template<typename T1>
		Matrix<RESULT_TYPE(+), M, N> operator+ (const Matrix<T1,M,N>& _mat1);  // TESTED

		/// \brief Compare if two matrices are identical (using elementary !=).
		bool operator== (const Matrix<T,M,N>& _mat1);                          // TESTED
	};

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





	// Include implementation.
#	include "details/matrix.inl"

	// Remove helper macro.
#	undef RESULT_TYPE
}
