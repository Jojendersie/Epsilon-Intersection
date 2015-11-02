#pragma once

#include "elementarytypes.hpp"

namespace details {
    // Helper shift for signed right operands. It corresponds to a right shift
    // for positive values and to a left shift for negative ones.
    ei::uint32 xshift(ei::uint32 _value, int _count)
    {
        return (_count < 0) ? _value << -_count : _value >> _count;
    }
    ei::uint64 xshift(ei::uint64 _value, int _count)
    {
        return (_count < 0) ? _value << -_count : _value >> _count;
    }

    // Chech if there are bits which overflow when shifted left by a given amount
    bool overflows(ei::uint32 _value, int _leftShift)
    {
        /// Create a bit mask with 111..0000 where the number of ones is _leftShift.
        return (_value & (0xffffffffu << (32-_leftShift))) != 0;
    }
    bool overflows(ei::uint64 _value, int _leftShift)
    {
        /// Create a bit mask with 111..0000 where the number of ones is _leftShift.
        return (_value & (0xffffffffffffffffu << (64-_leftShift))) != 0;
    }

    // Compute negated two complement value inplace
    template<int N>
    void negate(ei::int32* intrep)
    {
        int carry = 1;
        for(int i = N-1; i >= 0; --i)
        {
            intrep[i] = ~intrep[i] + carry;
            if(intrep[i]) carry = 0;
        }
    }
}

namespace ei
{
    /// Signed fixed point numbers. You can freely choose a size
    /// (8, 16 or 32*n) and the number of binary digits after the binary point.
    ///
    /// The internal representation is based on two complements and inherites
    /// some properties. The negative maximum is always one step large in its
    /// absolute value than its positive: |nmax| > max.
    /// E.g.: Fix<32,31> has the minimum -1 and the maximum (2^31-1)/(2^31)
    /// (that is 0.9999999995343387126922607421875).
    /// This also includes that abs(nmax) = nmax, because it immediatelly
    /// overflows on conversion (note that this is the same for common integer
    /// types. E.g. for int32: abs(-2147483648) == -2147483648).
    ///
    /// Fixed point numbers can be converted from and two float/double. Thereby
    /// too large and too small values are clamped to the positive and negative
    /// maxima.
    /// NaNs are converted to nmax.
    template<unsigned BitSize, unsigned FractionalBits>
    struct Fix
    {
        enum {NUM_INTS = BitSize/32};
        static_assert(FractionalBits < BitSize, "Number of fractional bits larger than entire type! Remember: there is a sign bit.");
        static_assert(NUM_INTS*32 == BitSize, "BitSize is not a multiple of 32!");
        static_assert(NUM_INTS >= 2, "There is a 32 bit specialization. This class should only be compiled for larger types!");

        enum {
            NUM_INT_DIGITS2  = BitSize - FractionalBits - 1,
            NUM_INT_DIGITS10 = int((BitSize - FractionalBits) / 3.321928095 + 1),
            NUM_FRAC_DIGITS2 = FractionalBits,
            NUM_FRAC_DIGITS10 = FractionalBits,
        };

        /// Construct from float
        explicit Fix(float _value);
        /// Construct from double
        explicit Fix(double _value);

        /// Cast explicit (not automatic) to float
        explicit operator float () const;
        /// Cast explicit (not automatic) to double
        explicit operator double () const;

        /// TODO: there is no? guarantee what happens with negative numbers on shift
        /// maybe it is necessary to store all but the first as uint32.
        int32 intrep[NUM_INTS];      ///< Internal integer representation with [0] as the most significant word
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

        /// Cast explicit (not automatic) to float
        explicit operator float () const;
        /// Cast explicit (not automatic) to double
        explicit operator double () const;

        /// Convert into a string. The buffer must be allocated and large enough.
        /// \param [in] _buf A buffer for the result with at least TODO bytes.
        //const char* Fix<32, FractionalBits>::toString(char* _buf) const;
    };

    // ********************************************************************* //
    // IMPLEMENTATION Fix<32 * N>                                            //
    // ********************************************************************* //
    template<unsigned BitSize, unsigned FractionalBits>
    Fix<BitSize, FractionalBits>::Fix(float _value)
    {
        // Zero out most words (not all might be influenced)
        for(int i = 0; i < NUM_INTS; ++i) intrep[i] = 0;
        // Split float into its parts
        int32 ix = *reinterpret_cast<int32*>(&_value);
        uint32 mantissa = (ix & 0x007fffff);
        int32 exponent = ((ix >> 23) & 0x000000ff) - 127;
        int32 sign = ix & 0x80000000;
        // Not null or denormalized -> add implicit leading bit
        if(exponent > -127)
            mantissa |= 0x00800000;
        // Put it in up to two words shifted
        int32 rightShift = NUM_INT_DIGITS2 - 8 - exponent;
        // Shift left -> only first word influenced. Overflow possible.
        if(rightShift <= 0)
        {
            if(details::overflows(mantissa, -rightShift + 1))
                if(sign) intrep[0] = 0x80000000;
                else {intrep[0] = 0x7fffffff; for(int i=1; i<NUM_INTS;++i) intrep[i] = 0xffffffff;}
            else intrep[0] = mantissa << -rightShift;
        } else {
            // Otherwise there are at most two words influenced. The range of
            // influence is in the bits 8+rightshift till 31+rightshift (0 indexed).
            int lw = (8+rightShift) / 32;
            int rw = (31+rightShift) / 32;
            if(lw < NUM_INTS)
                intrep[lw] = details::xshift(mantissa, rightShift - lw * 32);
            if(lw != rw && rw < NUM_INTS)
                intrep[rw] = details::xshift(mantissa, rightShift - rw * 32);
            if(sign) details::negate<NUM_INTS>(intrep);
        }
    }

    template<unsigned BitSize, unsigned FractionalBits>
    Fix<BitSize, FractionalBits>::Fix(double _value)
    {
        // Zero out most words (not all might be influenced)
        for(int i = 0; i < NUM_INTS; ++i) intrep[i] = 0;
        // Split double into its parts
        int64 ix = *reinterpret_cast<int64*>(&_value);
        uint64 mantissa = (ix & 0x000fffffffffffff);
        int32 exponent = ((ix >> 52) & 0x000007ff) - 1023;
        int64 sign = ix & 0x8000000000000000;
        // Not null or denormalized -> add implicit leading bit
        if(exponent > -1023)
            mantissa |= 0x0010000000000000;
        // Put it in up to three words shifted
        int32 rightShift = NUM_INT_DIGITS2 - 11 - exponent;
        // Shift left -> only first two word influenced. Overflow possible.
        if(rightShift <= 0)
        {
            if(details::overflows(mantissa, -rightShift + 1))
                if(sign) intrep[0] = 0x80000000;
                else {intrep[0] = 0x7fffffff; for(int i=1; i<NUM_INTS;++i) intrep[i] = 0xffffffff;}
            else {
                mantissa <<= -rightShift;
                intrep[0] = mantissa >> 32;
                intrep[1] = mantissa & 0xffffffffu;
            }
        } else {
            // Otherwise there are at most three words influenced. The range of
            // influence is in the bits 11+rightshift till 63+rightshift (0 indexed).
            int lw = (11+rightShift) / 32;
            int rw = (63+rightShift) / 32;
            for(int w = lw; w<=rw && w<NUM_INTS; ++w)
                intrep[w] = details::xshift(mantissa, rightShift - w * 32 + 32) & 0xffffffffu;
            if(sign) details::negate<NUM_INTS>(intrep);
        }
    }

    // ********************************************************************* //
    // Casts
    template<unsigned BitSize, unsigned FractionalBits>
    Fix<BitSize, FractionalBits>::operator float () const
    {
        // Find the two most significant words and put them into a 64 bit int.
        int64 mantissa = 0;
        int significantWord = 0; // Index of the first word != 0
        for(; significantWord < NUM_INTS-1; ++significantWord)
        {
            if(intrep[significantWord])
            {
                mantissa = (uint64(intrep[significantWord]) << 32) | uint32(intrep[significantWord+1]);
                break;
            }
        }
        if(significantWord == NUM_INTS-1)
            mantissa = int64(intrep[significantWord]) << 32;
        else ++significantWord;

        // Convert into float and shift by changing the exponent
        return float(mantissa) * pow(2.0f, int((NUM_INTS-significantWord-1) * 32 - FractionalBits));
    }

    template<unsigned BitSize, unsigned FractionalBits>
    Fix<BitSize, FractionalBits>::operator double () const
    {
        // Make number positive
        T x = *this;
        if(intrep[0] & 0x80000000) details::negate<NUM_INTS>(x.intrep);
        // Find the three most significant words and put them into a 64 bit int.
        int64 mantissa = 0;
        int msb = (NUM_INTS-2) * 32; // Index of the first bit != 0
        for(int w = 0; w < NUM_INTS; ++w)
        {
            if(x.intrep[w])
            {
                // Get position of the most significant bit -1
                double tmp = (double)(x.intrep[w]|1);
                msb = 1054 - (((*reinterpret_cast<uint64*>(&tmp)) >> 52) & 0x7ff) - 1;
                // Pack mantissa from up to three words
                mantissa = (uint64(x.intrep[w]) << (32+msb));
                if(w+1 < NUM_INTS) mantissa |= (uint64(uint32(x.intrep[w+1])) << msb);
                if(w+2 < NUM_INTS) mantissa |= (uint32(x.intrep[w+2]) >> (32-msb));
                msb += w*32;
                break;
            }
        }

        // Convert into double and shift by changing the exponent
        if(intrep[0] & 0x80000000)
            return -double(mantissa) * pow(2.0, int((NUM_INTS-2) * 32 - msb - FractionalBits));
        else return double(mantissa) * pow(2.0, int((NUM_INTS-2) * 32 - msb - FractionalBits));
    }

    // ********************************************************************* //
    // IMPLEMENTATION Fix<32>                                                //
    // ********************************************************************* //

    // ********************************************************************* //
    // Construction
    template<unsigned FractionalBits>
    Fix<32, FractionalBits>::Fix(float _value)
    {
        // Convert to a large int by shifting (i.e. change exponent).
        double val = _value * double(1u << FractionalBits);
        // Check for under and overflow
        if(val >= double(1u << 31)) intrep = 0x7fffffff;
        else if(val <= -double(1u << 31)) intrep = 0x80000000;
        else intrep = int32(val);
    }

    template<unsigned FractionalBits>
    Fix<32, FractionalBits>::Fix(double _value)
    {
        // Convert to a large int by shifting (i.e. change exponent).
        double val = _value * double(1u << FractionalBits);
        // Check for under and overflow
        if(val >= double(1u << 31)) intrep = 0x7fffffff;
        else if(val <= -double(1u << 31)) intrep = 0x80000000;
        else intrep = int32(val);
    }

    template<unsigned FractionalBits>
    Fix<32, FractionalBits>::Fix(int _value)
    {
        if(details::overflows(unsigned(_value & 0x7fffffff), FractionalBits+1))
            intrep = _value < 0 ? 0x80000000 : 0x7fffffff;
        else intrep = int32(_value << FractionalBits);
    }

    // ********************************************************************* //
    // Casts
    template<unsigned FractionalBits>
    Fix<32, FractionalBits>::operator float () const
    {
        // Mantissa of double > 32 bit -> save to compute anything in floats
        return float(intrep / double(1u << FractionalBits));
    }

    template<unsigned FractionalBits>
    Fix<32, FractionalBits>::operator double () const
    {
        // Mantissa of double > 32 bit -> save to compute anything in floats
        return intrep / double(1u << FractionalBits);
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

    template<unsigned FractionalBits>
    bool operator != (Fix<32,FractionalBits> _lhs, Fix<32,FractionalBits> _rhs)
    {
        return _lhs.intrep != _rhs.intrep;
    }

    template<unsigned FractionalBits>
    bool operator <= (Fix<32,FractionalBits> _lhs, Fix<32,FractionalBits> _rhs)
    {
        return _lhs.intrep <= _rhs.intrep;
    }

    template<unsigned FractionalBits>
    bool operator >= (Fix<32,FractionalBits> _lhs, Fix<32,FractionalBits> _rhs)
    {
        return _lhs.intrep >= _rhs.intrep;
    }

    template<unsigned FractionalBits>
    bool operator < (Fix<32,FractionalBits> _lhs, Fix<32,FractionalBits> _rhs)
    {
        return _lhs.intrep < _rhs.intrep;
    }

    template<unsigned FractionalBits>
    bool operator > (Fix<32,FractionalBits> _lhs, Fix<32,FractionalBits> _rhs)
    {
        return _lhs.intrep > _rhs.intrep;
    }
}