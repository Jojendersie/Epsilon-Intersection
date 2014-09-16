#pragma once

#include "config.hpp"
#include "details/inttemplate.hpp"

#ifndef USE_ELEMENTARIES_WITHOUT_NAMESPACE
namespace gam {
#endif

	/// \brief Short name for unsigned / unsigned int.
	/// \details The number of bit might vary on different systems.
    typedef unsigned int uint;												   // TESTED

	

	// Declaration of fixed sized int8, uint8 and byte = uint8
	typedef details::Int<1>::utype uint8;                                      // TESTED
	typedef details::Int<1>::utype byte;                                       // TESTED
	typedef details::Int<1>::stype int8;                                       // TESTED

	// Declaration of fixed sized int16, uint16
	typedef details::Int<2>::utype uint16;                                     // TESTED
	typedef details::Int<2>::stype int16;                                      // TESTED

	// Declaration of fixed sized int32, uint32
	typedef details::Int<4>::utype uint32;                                     // TESTED
	typedef details::Int<4>::stype int32;                                      // TESTED

	// Declaration of fixed sized int64, uint64
	typedef details::Int<8>::utype uint64;                                     // TESTED
	typedef details::Int<8>::stype int64;                                      // TESTED


#ifndef USE_ELEMENTARIES_WITHOUT_NAMESPACE
}
#endif

namespace gam {
	// ********************************************************************* //
	//								 FUNCTIONS								 //
	// ********************************************************************* //

	// ********************************************************************* //
	/// \brief Compute the square x*x.
	template<typename T>
	T sq(T _x);

	// ********************************************************************* //
	/// \brief Get the maximum from x and y.
	/// \details In case of x and y are equal, x is returned. This only makes
	///    a difference if you are sorting object types with more than the
	///    compared value.
	template<typename T>
	T max(T _x, T _y);

	// ********************************************************************* //
	/// \brief Get the minimum from x and y.
	/// \details In case of x and y are equal, x is returned. This only makes
	///    a difference if you are sorting object types with more than the
	///    compared value.
	template<typename T>
	T min(T _x, T _y);

	// ********************************************************************* //
	/// \brief Get the absolute value.
	template<typename T>
	T abs(T _x);

	// ********************************************************************* //
	/// \brief Get the sign of a value.
	/// \details There is a faster version sgn(), if you don't need to 
	///    know about zero.
	/// \returns -1 (_x < 0), 0 (_x == 0) or 1 (_x > 0)
	template<typename T>
	T sign(T _x);

	// ********************************************************************* //
	/// \brief Get the sign of a value where 0 is counted as positive.
	/// \details This function is faster than sign(). Use it if you don't need
	///    to know about zero.
	/// \returns -1 (_x < 0) or 1 (_x >= 0)
	template<typename T>
	T sgn(T _x);

	// Include implementation.
#	include "details/elementary.inl"
}