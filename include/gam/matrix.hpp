#pragma once

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
	template<typename T, unsigned M, unsigned N>
	class Matrix: public details::Components<T, M, N>
	{
		static_assert(M * N > 1, "A matrix must have at least 1x2 or 2x1 components.");


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


	// Include implementation.
#	include "details/matrix.inl"
}
