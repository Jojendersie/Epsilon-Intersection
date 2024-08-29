#pragma once

#include "configure.hpp"
#include "details/inttemplate.hpp"

#include <type_traits>
#include <limits>
#include <cmath>
#include <math.h>
#include <cstring>

namespace eitypes {
    /// \brief Short name for unsigned / unsigned int.
    /// \details The number of bit might vary on different systems.
    typedef unsigned int uint;                                                 // TESTED

    // Declaration of fixed sized int8, uint8 and byte = uint8
    typedef ei::details::Int<1>::utype uint8;                                  // TESTED
    typedef ei::details::Int<1>::utype byte;                                   // TESTED
    typedef ei::details::Int<1>::stype int8;                                   // TESTED

    // Declaration of fixed sized int16, uint16
    typedef ei::details::Int<2>::utype uint16;                                 // TESTED
    typedef ei::details::Int<2>::stype int16;                                  // TESTED

    // Declaration of fixed sized int32, uint32
    typedef ei::details::Int<4>::utype uint32;                                 // TESTED
    typedef ei::details::Int<4>::stype int32;                                  // TESTED

    // Declaration of fixed sized int64, uint64
    typedef ei::details::Int<8>::utype uint64;                                 // TESTED
    typedef ei::details::Int<8>::stype int64;                                  // TESTED
}

#ifdef EI_GLOBAL_ELEMENTARIES
    using namespace eitypes;
#endif
namespace ei {
    using namespace eitypes;
}

namespace ei { namespace details {
    /// \brief Dummy class to detect correct types for matrix <-> matrix
    ///     and matrix <-> scalar operations (and the same for quaternions).
    /// \details The overloading mechanism fails when both types of operations
    ///     are templated, because the matrix <-> scalar is chosen even if
    ///     both operants are matrices.
    class NonScalarType    {};

    template<typename T, typename F>
    EIAPI T hard_cast(F _from)
    {
        static_assert(sizeof(T) == sizeof(F), "Cannot cast types of different sizes");
        T to;
        std::memcpy(&to, &_from, sizeof(T));
        return to;
    }

}} // namespace ei::details

namespace ei {
    // ********************************************************************* //
    //                             MATH CONSTANTS                            //
    // ********************************************************************* //
    constexpr float PI = 3.141592654f;
    constexpr float E = 2.718281828f;
    constexpr float GOLDEN_RATIO = 1.61803398875f;
    constexpr float PHYTAGORAS = 1.4142135623f;
    constexpr float FMIN = std::numeric_limits<float>::min();
    // The cmath header has an ugly macro with name INFINITY -> name conflict
    constexpr float INF { std::numeric_limits<float>::infinity() };
    constexpr double INF_D { std::numeric_limits<double>::infinity() };
    // Unicode names for the above constants
#ifdef EI_USE_UNICODE_NAMES
    constexpr float π = PI;
    constexpr float Φ = GOLDEN_RATIO;
    constexpr float ℇ = E;
#endif

    // ********************************************************************* //
    //                               FUNCTIONS                               //
    // ********************************************************************* //

    // ********************************************************************* //
    /// \brief Compute the square x*x.
    template<typename T>
    constexpr EIAPI T sq(T _x) noexcept
    {
        return _x * _x;
    }

    // ********************************************************************* //
    /// \brief Get the maximum from x and y.
    /// \details In case of x and y are equal, x is returned. This only makes
    ///    a difference if you are sorting object types with more than the
    ///    compared value.
    template<typename T>
    constexpr EIAPI T max(T _x, T _y) noexcept // TESTED
    {
        return _x < _y ? _y : _x;
    }

    /// \brief Get the maximum of any number of arguments
    /// \details In case of equal arguments the left most one is returned
    template<typename T, typename... Ttail>
    constexpr EIAPI T max(T _first, Ttail... _tail) noexcept // TESTED
    {
        return max(_first, max(_tail...));
    }

    // ********************************************************************* //
    /// \brief Get the minimum from x and y.
    /// \details In case of x and y are equal, x is returned. This only makes
    ///    a difference if you are sorting object types with more than the
    ///    compared value.
    template<typename T>
    constexpr EIAPI T min(T _x, T _y) noexcept // TESTED
    {
        return _x > _y ? _y : _x;
    }

    /// \brief Get the maximum of any number of arguments
    /// \details In case of equal arguments the left most one is returned
    template<typename T, typename... Ttail>
    constexpr EIAPI T min(T _first, Ttail... _tail) noexcept // TESTED
    {
        return min(_first, min(_tail...));
    }

    // ********************************************************************* //
    /// \brief Clamp a value to the boundaries.
    template<typename T>
    constexpr EIAPI T clamp(T _x, T _min, T _max) noexcept // TESTED
    {
        return _x > _max ? _max : (_x < _min ? _min : _x);
    }

    // ********************************************************************* //
    /// \brief Clamp a value to [0,1] interval.
    template<typename T>
    constexpr EIAPI T saturate(T _x) noexcept
    {
        return _x > static_cast<T>(1) ? static_cast<T>(1) : (_x < static_cast<T>(0) ? static_cast<T>(0) : _x);
    }


    // ********************************************************************* //
    /// \brief Safe division x / y, prevents division by zero.
    /// \details This always produces a finite result except _x is already
    ///     infinite or _x or _y is nan.
    template<typename T>
    constexpr EIAPI T sdiv(T _x, T _y) noexcept
    {
        if (_y == static_cast<T>(0))
            return _x < 0.0f ? -std::numeric_limits<T>::max() : (_x > 0.0f ? std::numeric_limits<T>::max() : 0.0f);
        return _x / _y;
    }

    // ********************************************************************* //
    /// \brief Get the absolute value.
    template<typename T>
    constexpr EIAPI T abs(T _x) noexcept // TESTED
    {
        return _x < static_cast<T>(0) ? -_x : _x;
    }

    constexpr EIAPI float abs(float _x) noexcept // TESTED
    {
        // In float there is a negative -0 -> simple ?: does not work
        // The following bit manipulation does not work in constexpr due to undefined behavior.
        //using details::hard_cast;
        //return hard_cast<float>(hard_cast<uint32>(_x) & 0x7fffffff);
        // The following are valid. TODO: benchmark
        return _x < 0.0f ? -_x : (_x == 0.0f ? 0.0f : _x);
        //return _x < 0.0f ? -_x : _x + 0.0f; // IEEE float assumption
    }
    constexpr EIAPI double abs(double _x) noexcept // TESTED
    {
        return _x < 0.0 ? -_x : (_x == 0.0 ? 0.0 : _x);
    }

    // ********************************************************************* //
    /// \brief Get the sign of a value.
    /// \details There is a faster version sgn(), if you don't need to
    ///    know about zero.
    /// \returns -1 (_x < 0), 0 (_x == 0) or 1 (_x > 0)
    template<typename T>
    constexpr EIAPI T sign(T _x) noexcept // TESTED
    {
        return _x < static_cast<T>(0) ? static_cast<T>(-1)
            : (_x > static_cast<T>(0) ? static_cast<T>(1) : static_cast<T>(0));
    }

    // ********************************************************************* //
    /// \brief Get the sign of a value where the sign of 0 is counted too.
    /// \details This function should be faster than sign().
    /// \returns -1 (_x <= -0) or 1 (_x >= 0)
    template<typename T>
    constexpr EIAPI T sgn(T _x) noexcept // TESTED
    {
        return _x < static_cast<T>(0) ? static_cast<T>(-1) : static_cast<T>(1);
    }
    constexpr EIAPI float sgn(float _x) noexcept // TESTED
    {
        if(_x == 0.0f) return std::signbit(_x) ? -1.0f : 1.0f;
        return _x < 0.0f ? -1.0f : 1.0f;
    }
    constexpr EIAPI double sgn(double _x) noexcept // TESTED
    {
        if(_x == 0.0) return std::signbit(_x) ? -1.0 : 1.0;
        return _x < 0.0 ? -1.0 : 1.0;
    }

    // ********************************************************************* //
    /// \brief Get 0 for (_x <= -0) or 1 for (_x >= 0).
    /// \returns 0 (_x <= -0) or 1 (_x >= 0)
    template<typename T>
    constexpr EIAPI int heaviside(T _x) noexcept
    {
        return _x < static_cast<T>(0) ? 0 : 1;
    }
    EIAPI int heaviside(float _x) noexcept // TESTED
    {
        return details::hard_cast<uint32>(_x) & 0x80000000 ? 0 : 1;
    }
    EIAPI int heaviside(double _x) noexcept // TESTED
    {
        return details::hard_cast<uint64>(_x) & 0x8000000000000000ull ? 0 : 1;
    }

    // ********************************************************************* //
    /// \brief Check if the relative difference between two scalars is less
    ///    or equal than epsilon.
    /// \details Computes 2*|x0-x1|/max(|x0|+|x1|, offset), i.e. the ratio
    ///     between the absolute difference and the clamped average absolute
    ///     value.
    ///
    ///     For numbers towards 0 the clamping (offset = 0.03125) produces an
    ///     absolute error again
    /// \param [in] _epsilon A percentage threshold for the difference
    ///    between two elements. The default value is 1e-6.
    /// \returns true if the difference is less or equal than _epsilon.
    template<typename T, class = std::enable_if_t<!std::is_base_of<details::NonScalarType, T>::value>>
    constexpr EIAPI bool approx(T _x0, T _x1, T _epsilon = T(1e-6)) noexcept // TESTED
    {
        // Use an offset of 1.0 for comparisons to zero.
        T sum = max(abs(_x0) + abs(_x1), static_cast<T>(1.0));
        return (static_cast<T>(2) * abs(_x1 - _x0) / sum) <= _epsilon;
    }

    // ********************************************************************* //
    /// \brief Round value towards negative infinity.
    template<typename T, class = std::enable_if_t<!std::is_base_of<details::NonScalarType, T>::value>>
    constexpr EIAPI Sint<sizeof(T)> floor(T _x) noexcept
    {
        Sint<sizeof(T)> r = static_cast<Sint<sizeof(T)>>(_x);
        return r - static_cast<Sint<sizeof(T)>>((_x<static_cast<T>(0)) && (_x-r!=static_cast<T>(0)));
    }

    // ********************************************************************* //
    /// \brief Round value towards positive infinity.
    template<typename T, class = std::enable_if_t<!std::is_base_of<details::NonScalarType, T>::value>>
    constexpr EIAPI Sint<sizeof(T)> ceil(T _x) noexcept
    {
        Sint<sizeof(T)> r = static_cast<Sint<sizeof(T)>>(_x);
        return r + static_cast<Sint<sizeof(T)>>((_x>static_cast<T>(0)) && (_x-r!=static_cast<T>(0)));
    }


    // ********************************************************************* //
    /// \brief Round value towards next integral number (0.5 rounds to even).
    template<typename T, class = std::enable_if_t<!std::is_base_of<details::NonScalarType, T>::value>>
    constexpr EIAPI Sint<sizeof(T)> round(T _x) noexcept
    {
        // Round up
        //return floor(_x + static_cast<T>(0.5));
        // Round to even
        Sint<sizeof(T)> r = static_cast<Sint<sizeof(T)>>(_x);
        r -= static_cast<Sint<sizeof(T)>>((_x<static_cast<T>(0)) && (_x!=static_cast<T>(r)));   // Subtract 1 if _x < 0 -> r is floor(_x)
        T f = _x - r;    // Fractional part (positive)
        if(f < static_cast<T>(0.5)) return r;
        if(f > static_cast<T>(0.5)) return r+1;
        // f is 0.5
        if(r & 1) return r + 1;
        return r;
    }

    // ********************************************************************* //
    /// \brief Get the fraction in (-1,1) with f-int(f).
    /// \param _x [in] The number to be splitted.
    template<typename T>
    constexpr EIAPI T frac(T _x) noexcept // TESTED
    {
        return _x - static_cast<Sint<sizeof(T)>>(_x);
    }

    // ********************************************************************* //
    /// \brief Divide a number into the integer and fractional part.
    /// \param _x [in] The number to be splitted.
    /// \param _int [out] The integer part of the number.
    /// \returns The fraction of the number in (-1,1).
    template<typename T>
    constexpr EIAPI T intfrac(T _x, Sint<sizeof(T)>& _int) noexcept // TESTED
    {
        _int = static_cast<Sint<sizeof(T)>>(_x);
        return _x - _int;
    }

    // ********************************************************************* //
    /// \brief Divide a number into an integer and a positive fractional part.
    /// \details This method uses f-floor(f) instead of f-int(f) and therefore
    ///      has a continous behavior around zero
    /// \param _x [in] The number to be splitted.
    /// \param _int [out] The integer part of the number (rounded to zero).
    /// \returns The fraction of the number in [0,1).
    template<typename T>
    constexpr EIAPI T floorfrac(T _x, Sint<sizeof(T)>& _int) noexcept // TESTED
    {
        _int = floor(_x);
        return _x - _int;
    }

    // ********************************************************************* //
    /// \brief Get the smallest positive number m such that x=y*c+m with c in Z.
    /// \returns The mathematically defined positive modulus.
    template<typename T>
    constexpr EIAPI T mod(T _x, T _y) noexcept // TESTED
    {
        eiAssert(_y != 0.0f, "Modulo 0 is not defined!");
        T m = fmod(_x, _y);
        return m < 0 ? m+abs(_y) : m;
    }

    // Pure integer specializations
    template<>
    constexpr EIAPI int8 mod<int8>(int8 _x, int8 _y) noexcept
    {
        eiAssert(_y != 0, "Modulo 0 is not defined!");
        int8 m = _x % _y;
        return m < 0 ? m+abs(_y) : m;
    }
    template<>
    constexpr EIAPI int16 mod<int16>(int16 _x, int16 _y) noexcept
    {
        eiAssert(_y != 0, "Modulo 0 is not defined!");
        int16 m = _x % _y;
        return m < 0 ? m+abs(_y) : m;
    }
    template<>
    constexpr EIAPI int32 mod<int32>(int32 _x, int32 _y) noexcept
    {
        eiAssert(_y != 0, "Modulo 0 is not defined!");
        int32 m = _x % _y;
        return m < 0 ? m+abs(_y) : m;
    }
    template<>
    constexpr EIAPI int64 mod<int64>(int64 _x, int64 _y) noexcept
    {
        eiAssert(_y != 0, "Modulo 0 is not defined!");
        int64 m = _x % _y;
        return m < 0 ? m+abs(_y) : m;
    }

    template<>
    constexpr EIAPI uint8 mod<uint8>(uint8 _x, uint8 _y) noexcept { return _x % _y; }
    template<>
    constexpr EIAPI uint16 mod<uint16>(uint16 _x, uint16 _y) noexcept { return _x % _y; }
    template<>
    constexpr EIAPI uint32 mod<uint32>(uint32 _x, uint32 _y) noexcept { return _x % _y; }
    template<>
    constexpr EIAPI uint64 mod<uint64>(uint64 _x, uint64 _y) noexcept { return _x % _y; }

    // ********************************************************************* //
    /// \brief Compute floor(log2), i.e. the position of the most significant bit.
    /// \details Expects IEEE and litte-endian (standard PC)
    template<typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
    int ilog2(T _x) noexcept                    // TESTED
    {
        if(_x <= 0) return -2147483647-1; // Actually NaN
        double f = static_cast<double>(_x);
        return (details::hard_cast<uint64>(f) >> 52) - 1023;
    }


    // ********************************************************************* //
    /// \brief Linear interpolation.
    /// \details There are two formulations:
    ///    * x * (1 - t) + y * t
    ///    * x + (y - x) * t
    ///    The second one does not need a constant 1 (might be the type does
    ///    not support this) and is faster by one scalar operation.
    /// \param _x0 [in] Scalar, vector or matrix value. This is returned when
    ///    _t is zero.
    /// \param _x1 [in] Scalar, vector or matrix value. This is returned when
    ///    _t is one.
    /// \param _t [in] Interpolation parameter. Can be scalar or vector.
    /// \returns x + (y - x) * t where the type is derived from the operands.
    template<typename T0, typename T1>
    constexpr EIAPI auto lerp(T0 _x0, T0 _x1, T1 _t) noexcept -> decltype(_x0*_t) // TESTED
    {
        return _x0 + (_x1 - _x0) * _t;
    }

    // ********************************************************************* //
    /// \brief Bilinear interpolation optimized for scalars.
    /// \param _x00 [in] Scalar value. This is returned when
    ///    _t0 is zero and _t1 is zero.
    /// \param _x01 [in] Scalar value. This is returned when
    ///    _t0 is one and _t1 is zero.
    /// \param _x10 [in] Scalar value. This is returned when
    ///    _t0 is zero and _t1 is one.
    /// \param _x11 [in] Scalar value. This is returned when
    ///    _t0 is one and _t1 is one.
    /// \param _t0 [in] Scalar interpolation parameter ("x-direction").
    /// \param _t1 [in] Scalar interpolation parameter ("y-direction").
    /// \returns lerp(lerp(_x00, _x01, _t0), lerp(_x10, _x11, _t0), _t1).
    template<typename T0, typename T1>
    constexpr EIAPI auto bilerp(T0 _x00, T0 _x01,
                                 T0 _x10, T0 _x11,
                                 T1 _t0, T1 _t1) noexcept -> decltype(_x00*_t0) // TESTED
    {
        return lerp(lerp(_x00, _x01, _t0),
                    lerp(_x10, _x11, _t0), _t1);
    }

    // ********************************************************************* //
    /// \brief Evalutate the smoothstep polynomial 3t^2 - 2t^3.
    /// \param _t [in] The value to be inserted into the polynomial. The useful
    ///    definition interval is in [0,1].
    template<typename T>
    constexpr EIAPI T smoothstep(T _t) noexcept
    {
        return _t * _t * (T(3) - T(2) * _t);
    }

    // ********************************************************************* //
    /// \brief Evalutate the smootherstep polynomial 6t^5 - 15t^4 + 10t^3.
    /// \param _t [in] The value to be inserted into the polynomial. The useful
    ///    definition interval is in [0,1].
    template<typename T>
    constexpr EIAPI T smootherstep(T _t) noexcept
    {
        return _t * _t * _t * (_t * (_t * T(6) - T(15)) + T(10));
    }

    // ********************************************************************* //
    /// \brief Compute the next machine representable number in positive direction.
    /// \details Iterates overall floating point numbers, even through
    ///    denormalized range. Thereby the successor of -0 and 0 are both the
    ///    same smallest positive float. Further -INF is transformed into
    ///    -FLOAT_MAX and +INF remains +INF.
    EIAPI float successor(float _number) noexcept // TESTED
    {
        // Handling NaN
        if(_number != _number) return _number;
        // Infinity
        if(_number == INF) return _number;
        if(_number == -INF) return details::hard_cast<float>(0xff7fffffu);

        uint32 bits = details::hard_cast<uint32>(_number);
        uint32 sign = bits & 0x80000000u;
        uint32 mantissa = bits & 0x7fffffffu;
        if( sign && mantissa )
            return details::hard_cast<float>( sign | (mantissa - 1) );
        else return details::hard_cast<float>( mantissa + 1 );
    }

    EIAPI double successor(double _number) noexcept
    {
        // Handling NaN
        if(_number != _number) return _number;
        // Infinity
        if(_number == INF_D) return _number;
        if(_number == -INF_D) return details::hard_cast<double>(0xffeffffffffffffful);

        uint64 bits = details::hard_cast<uint64>(_number);
        uint64 sign = bits & 0x8000000000000000ul;
        uint64 mantissa = bits & 0x7ffffffffffffffful;
        if( sign && mantissa )
            return details::hard_cast<double>( sign | (mantissa - 1) );
        else return details::hard_cast<double>( mantissa + 1 );
    }

    /// \brief Compute the next machine representable number in negative direction.
    /// \details Iterates overall floating point numbers, even through
    ///    denormalized range. Thereby the predecessor of -0 and 0 are both the
    ///    same smallest negative float. Further -INF remains -INF and + INF
    ///    is transformed into FLOAT_MAX.
    EIAPI float predecessor(float _number) noexcept // TESTED
    {
        // Handling NaN
        if(_number != _number) return _number;
        // Infinity
        if(_number == -INF) return _number;
        if(_number == INF) return details::hard_cast<float>(0x7f7fffffu);

        uint32 bits = details::hard_cast<uint32>(_number);
        uint32 sign = bits & 0x80000000u;
        uint32 mantissa = bits & 0x7fffffffu;
        if( sign || !mantissa )
            return details::hard_cast<float>( 0x80000000u | (mantissa + 1) );
        else return details::hard_cast<float>( mantissa - 1 );
    }

    EIAPI double predecessor(double _number) noexcept
    {
        // Handling NaN
        if(_number != _number) return _number;
        // Infinity
        if(_number == -INF) return _number;
        if(_number == INF) return details::hard_cast<double>(0x7feffffffffffffful);

        uint64 bits = details::hard_cast<uint64>(_number);
        uint64 sign = bits & 0x8000000000000000ul;
        uint64 mantissa = bits & 0x7ffffffffffffffful;
        if( sign || !mantissa )
            return details::hard_cast<double>( 0x8000000000000000ul | (mantissa + 1) );
        else return details::hard_cast<double>( mantissa - 1 );
    }

    /// \brief Helper method to solve ax^2 + bx + c (numerically more stable than naive method).
    /// \returns Solutions x1, x2 (can be the same if a~0).
    template<typename T>
    constexpr EIAPI bool solveSquarePoly(T a, T b, T c, T& x1, T& x2)
    {
        // Special case: linear equation
        if(a == 0) {
            if(b == 0) return false;
            x1 = -c / b;
            x2 = x1;
            return true;
        }
        T discriminant = b*b - static_cast<T>(4)*a*c;
        if(discriminant < 0) return false;
        T dsqrt = sqrt(discriminant);
        if(b > 0) {
            x1 = -(b + dsqrt)/(static_cast<T>(2)*a);
            x2 = -static_cast<T>(2)*c/(b + dsqrt);
        } else {
            x1 = -static_cast<T>(2)*c/(b - dsqrt);
            x2 = -(b - dsqrt)/(static_cast<T>(2)*a);
        }
        return true;
    }


    // ********************************************************************* //
    //                 NORMALIZED INTEGERS (FIXED POINT)                     //
    // ********************************************************************* //

    /// \brief Storage type to discretize [-1,1] into a signed integer or [0,1]
    ///      into an unsigned integer.
    /// \details The encoding correctly represents -1, 0 and 1.
    template<typename T, uint Bits = sizeof(T)*8>
    struct NormalizedInt
    {
    private:
        // Template helpers to get some compiletime contants
        template<typename T1, typename Dummy = void> struct MinInterval {
            static constexpr float value = -1.0f;
        };
        template<typename Dummy> struct MinInterval<typename details::Int<sizeof(T)>::utype, Dummy> {
            static constexpr float value = 0.0f;
        };
        template<typename T1, typename Dummy = void> struct MaxValue {
            static constexpr uint64 value = 0xffffffffffffffffull >> (64u - Bits + 1u);
        };
        template<typename Dummy> struct MaxValue<typename details::Int<sizeof(T)>::utype, Dummy> {
            static constexpr uint64 value = 0xffffffffffffffffull >> (64u - Bits);
        };
    public:
        static_assert(std::is_integral<T>::value, "NormalizedInt can only use integer types as basis type.");
        static_assert(Bits <= sizeof(T)*8, "Base type as not enough bits.");
        //static constexpr float INTERVAL_MIN = []{ if(std::is_signed<T>::value) return -1.0f; else return 0.0f; }();
        static constexpr float INTERVAL_MIN = MinInterval<T>::value;
        static constexpr float INTERVAL_MAX = 1.0f;
        static constexpr uint BITS = Bits;
        //static constexpr uint64 MAX_POSITIVE_VALUE = []{ if(std::is_signed<T>::value) return 1ull << (BITS-1); else return 1ull << BITS; }();
        static constexpr uint64 MAX_POSITIVE_VALUE = MaxValue<T>::value;
        static constexpr uint64 MASK = (1ull << BITS) - 1;

        NormalizedInt() = default;

        template<typename T1, typename = std::enable_if_t<std::is_integral<T>::value && BITS <= sizeof(T1)*8>>
        constexpr EIAPI explicit NormalizedInt(T1 v) noexcept : value(v & MASK)
        {
            // Signed types in 2-complement must have all significant bits set.
            // This is violated if BITS < 8*sizeof(T).
            if(value > static_cast<T>(MAX_POSITIVE_VALUE))
                value |= 0xffffffffffffffffull << BITS;
        }

        constexpr EIAPI explicit NormalizedInt(float v)
        {
            eiAssert(v >= INTERVAL_MIN && v <= INTERVAL_MAX, "Cannot map an integer outside the unit interval to nint");
            value = floor(v * MAX_POSITIVE_VALUE + 0.5f) & MASK;
        }

        constexpr EIAPI  explicit NormalizedInt(double v)
        {
            eiAssert(v >= INTERVAL_MIN && v <= INTERVAL_MAX, "Cannot map an integer outside the unit interval to nint");
            value = floor(v * MAX_POSITIVE_VALUE + 0.5) & MASK;
        }

        constexpr EIAPI explicit operator float () const noexcept
        {
            return value / float(MAX_POSITIVE_VALUE);
        }

        constexpr EIAPI explicit operator double () const noexcept
        {
            return value / double(MAX_POSITIVE_VALUE);
        }

        template<typename T1, typename = std::enable_if_t<std::is_integral<T>::value && BITS <= sizeof(T1)*8>>
        constexpr EIAPI explicit operator T1 () const noexcept { return T1(value); }
    private:
        T value;
    };

    using nint8 = NormalizedInt<int8>;
    using nint16 = NormalizedInt<int16>;
    using nint32 = NormalizedInt<int32>;
    using nint64 = NormalizedInt<int64>;
    using nuint8 = NormalizedInt<uint8>;
    using nuint16 = NormalizedInt<uint16>;
    using nuint32 = NormalizedInt<uint32>;
    using nuint64 = NormalizedInt<uint64>;
}