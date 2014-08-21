#pragma once

#include "config.hpp"

namespace details {
	/// \brief Template construct to create integer types by a size value.
	/// \details The specializations define signed and unsigned types derived
	///		from a standard type. Using Int<2>::stype... will chose a
	///		specialization with the required size automatically.
	///
	///		The DUMMY parameter avoids the creation of two equal
	///		specializations in case one or multiple types have the same size.
	template <int BYTES, int DUMMY = 0> struct Int;
	template <> struct Int<sizeof(char)>
	{
		typedef signed char stype;
		typedef unsigned char utype;
	};

	template <> struct Int<sizeof(short), sizeof(short) == sizeof(char) ? 1 : 0>
	{
		typedef signed short stype;
		typedef unsigned short utype;
	};

	template <> struct Int<sizeof(int), sizeof(int) == sizeof(short) ? 2 : 0>
	{
		typedef signed int stype;
		typedef unsigned int utype;
	};

	template <> struct Int<sizeof(long), sizeof(long) == sizeof(int) ? 3 : 0>
	{
		typedef signed long stype;
		typedef unsigned long utype;
	};

	template <> struct Int<sizeof(long long), sizeof(long long) == sizeof(long) ? 4 : 0>
	{
		typedef signed long long stype;
		typedef unsigned long long utype;
	};
}


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
