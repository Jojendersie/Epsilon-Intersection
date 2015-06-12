#pragma once

#include "elementarytypes.hpp"

namespace ei
{
    /// Signed fixed point numbers. You can freely chose a size
    /// (8, 16 or 32*n) and the number of binary digits after the binary point.
    template<unsigned BitSize, unsigned FractionalBits>
    struct Fix
    {
        enum {NUM_INTS = BitSize/32};
        static_assert(FractionalBits < BitSize, "Number of fractional bits larger than entire type! Remember: there is a sign bit.");
        static_assert(NUM_INTS*32 == BitSize, "BitSize is not a multiple of 32!");

        enum {
            NUM_INT_DIGITS2  = BitSize - FractionalBits,
            NUM_INT_DIGITS10 = int((BitSize - FractionalBits) / 3.321928095 + 1),
            NUM_FRAC_DIGITS2 = FractionalBits,
            NUM_FRAC_DIGITS10 = FractionalBits,
        };

        int32 intrep[NUM_INTS];      ///< Internal integer representation
    private:
        typedef Fix<BitSize, FractionalBits> T;
    };

    // Specialization for 32 bit type
    template<unsigned FractionalBits>
    struct Fix<32, FractionalBits>
    {
        static_assert(FractionalBits < 32, "Number of fractional bits larger than entire type! Remember: there is a sign bit.");

        enum {
            NUM_INT_DIGITS2  = 32 - FractionalBits - 1,
            NUM_INT_DIGITS10 = int(NUM_INT_DIGITS2 / 3.321928095 + 1),
            NUM_FRAC_DIGITS2 = FractionalBits,
            NUM_FRAC_DIGITS10 = int(FractionalBits / 3.321928095 + 1),
        };

        int32 intrep;               ///< Internal integer representation

        /// Empty default construction
        Fix() {}

        /// Construct from float
        explicit Fix(float _value);
        /// Construct from double
        explicit Fix(double _value);
        /// Construct from int
        explicit Fix(int _value);

        /// Cast explicit to float
        explicit operator float () const;
        /// Cast explicit to double
        explicit operator double () const;

        /// Convert into a string. The buffer must be allocated and large enough.
        /// \param [in] _buf A buffer for the result with at least TODO bytes.
        //const char* Fix<32, FractionalBits>::toString(char* _buf) const;
    };

    // ********************************************************************* //
    // IMPLEMENTATION Fix<32>                                                //
    // ********************************************************************* //

    // ********************************************************************* //
    // Construction
    template<unsigned FractionalBits>
    Fix<32, FractionalBits>::Fix(float _value)
    {
        intrep = int32(_value * double(1 << FractionalBits));
    }

    template<unsigned FractionalBits>
    Fix<32, FractionalBits>::Fix(double _value)
    {
        intrep = int32(_value * double(1 << FractionalBits));
    }

    template<unsigned FractionalBits>
    Fix<32, FractionalBits>::Fix(int _value)
    {
        intrep = int32(_value << FractionalBits);
    }

    // ********************************************************************* //
    // Casts
    template<unsigned FractionalBits>
    Fix<32, FractionalBits>::operator float () const
    {
        // Mantissa of double > 32 bit -> save to compute anything in floats
        return float(intrep / double(1 << FractionalBits));
    }

    template<unsigned FractionalBits>
    Fix<32, FractionalBits>::operator double () const
    {
        // Mantissa of double > 32 bit -> save to compute anything in floats
        return intrep / double(1 << FractionalBits);
    }

    // ********************************************************************* //
    // Operators
    template<unsigned FractionalBits>
    Fix<32, FractionalBits> operator + (Fix<32, FractionalBits> _lhs, Fix<32, FractionalBits> _rhs)
    {
        Fix<32, FractionalBits> result;
        result.intrep = _lhs.intrep + _rhs.intrep;
        return result;
    }

    template<unsigned FractionalBits>
    Fix<32, FractionalBits> operator - (Fix<32, FractionalBits> _lhs, Fix<32, FractionalBits> _rhs)
    {
        Fix<32, FractionalBits> result;
        result.intrep = _lhs.intrep - _rhs.intrep;
        return result;
    }

    /// Unary minus
    template<unsigned FractionalBits>
    Fix<32, FractionalBits> operator - (Fix<32, FractionalBits> _x)
    {
        // Use ~x + 1 = -x from two's complement properties.
        _x.intrep = ~_x.intrep + 1;
        return _x;
    }

    template<unsigned FractionalBits>
    Fix<32, FractionalBits> operator * (Fix<32, FractionalBits> _lhs, Fix<32, FractionalBits> _rhs)
    {
        // Do computation with larger number to avoid intermediate overflow
        Fix<32, FractionalBits> result;
        result.intrep = int32(int64(_lhs.intrep) * int64(_rhs.intrep) >> FractionalBits);
        return result;
    }

    template<unsigned FractionalBits>
    Fix<32, FractionalBits> operator / (Fix<32, FractionalBits> _lhs, Fix<32, FractionalBits> _rhs)
    {
        // Do computation with larger number to avoid intermediate overflow
        Fix<32, FractionalBits> result;
        result.intrep = int32((int64(_lhs.intrep) << FractionalBits) / int64(_rhs.intrep));
        return result;
    }

    template<unsigned FractionalBits>
    bool operator == (Fix<32,FractionalBits> _lhs, Fix<32,FractionalBits> _rhs)
    {
        return _lhs.intrep == _rhs.intrep;
    }
}