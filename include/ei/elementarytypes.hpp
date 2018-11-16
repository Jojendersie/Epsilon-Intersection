#pragma once

#include "configure.hpp"
#include "details/inttemplate.hpp"

#include <type_traits>
// Ugly to include this here but you will miss sqrt,sin,... otherwise
#include <cmath>

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
    T hard_cast(F _from)
    {
        static_assert(sizeof(T) == sizeof(F), "Cannot cast types of different sizes");
        return *reinterpret_cast<T*>(&_from);
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
    // The cmath header has an ugly macro with name INFINITY -> name conflict
#pragma warning(push)
#pragma warning(disable:4056) // overflow in fp constant arithmetic
    const float INF = 1e30f * 1e30f;
    const double INF_D = 1e300 * 1e300;
#pragma warning(pop)
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
    inline T sq(T _x) noexcept
    {
        return _x * _x;
    }

    // ********************************************************************* //
    /// \brief Get the maximum from x and y.
    /// \details In case of x and y are equal, x is returned. This only makes
    ///    a difference if you are sorting object types with more than the
    ///    compared value.
    template<typename T>
    inline T max(T _x, T _y) noexcept // TESTED
    {
        return _x < _y ? _y : _x;
    }

    /// \brief Get the maximum of any number of arguments
    /// \details In case of equal arguments the left most one is returned
    template<typename T, typename... Ttail>
    inline T max(T _first, Ttail... _tail) noexcept // TESTED
    {
        return max(_first, max(_tail...));
    }

    // ********************************************************************* //
    /// \brief Get the minimum from x and y.
    /// \details In case of x and y are equal, x is returned. This only makes
    ///    a difference if you are sorting object types with more than the
    ///    compared value.
    template<typename T>
    inline T min(T _x, T _y) noexcept // TESTED
    {
        return _x > _y ? _y : _x;
    }

    /// \brief Get the maximum of any number of arguments
    /// \details In case of equal arguments the left most one is returned
    template<typename T, typename... Ttail>
    inline T min(T _first, Ttail... _tail) noexcept // TESTED
    {
        return min(_first, min(_tail...));
    }

    // ********************************************************************* //
    /// \brief Clamp a value to the boundaries.
    template<typename T>
    inline T clamp(T _x, T _min, T _max) noexcept // TESTED
    {
        return _x > _max ? _max : (_x < _min ? _min : _x);
    }

    // ********************************************************************* //
    /// \brief Clamp a value to [0,1] interval.
    template<typename T>
    inline T saturate(T _x) noexcept
    {
        return _x > static_cast<T>(1) ? static_cast<T>(1) : (_x < static_cast<T>(0) ? static_cast<T>(0) : _x);
    }

    // ********************************************************************* //
    /// \brief Get the absolute value.
    template<typename T>
    inline T abs(T _x) noexcept // TESTED
    {
        return _x < static_cast<T>(0) ? -_x : _x;
    }

    inline float abs(float _x) noexcept // TESTED
    {
        using details::hard_cast;
        return hard_cast<float>(hard_cast<uint32>(_x) & 0x7fffffff);
    }
    inline double abs(double _x) noexcept // TESTED
    {
        using details::hard_cast;
        return hard_cast<double>(hard_cast<uint64>(_x) & 0x7fffffffffffffffull);
    }

    // ********************************************************************* //
    /// \brief Get the sign of a value.
    /// \details There is a faster version sgn(), if you don't need to 
    ///    know about zero.
    /// \returns -1 (_x < 0), 0 (_x == 0) or 1 (_x > 0)
    template<typename T>
    inline T sign(T _x) noexcept // TESTED
    {
        return _x < static_cast<T>(0) ? static_cast<T>(-1)
            : (_x > static_cast<T>(0) ? static_cast<T>(1) : static_cast<T>(0));
    }

    // ********************************************************************* //
    /// \brief Get the sign of a value where the sign of 0 is counted too.
    /// \details This function should be faster than sign().
    /// \returns -1 (_x <= -0) or 1 (_x >= 0)
    template<typename T>
    inline T sgn(T _x) noexcept // TESTED
    {
        return _x < static_cast<T>(0) ? static_cast<T>(-1) : static_cast<T>(1);
    }
    inline float sgn(float _x) noexcept // TESTED
    {
        return details::hard_cast<uint32>(_x) & 0x80000000 ? -1.0f : 1.0f;
    }
    inline double sgn(double _x) noexcept // TESTED
    {
        return details::hard_cast<uint64>(_x) & 0x8000000000000000ull ? -1.0 : 1.0;
    }

    // ********************************************************************* //
    /// \brief Get 0 for (_x <= -0) or 1 for (_x >= 0).
    /// \returns 0 (_x <= -0) or 1 (_x >= 0)
    template<typename T>
    inline int heaviside(T _x) noexcept
    {
        return _x < static_cast<T>(0) ? 0 : 1;
    }
    inline int heaviside(float _x) noexcept // TESTED
    {
        return details::hard_cast<uint32>(_x) & 0x80000000 ? 0 : 1;
    }
    inline int heaviside(double _x) noexcept // TESTED
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
    inline bool approx(T _x0, T _x1, T _epsilon = T(1e-6)) noexcept // TESTED
    {
        // Use an offset of 1.0 for comparisons to zero.
        T sum = max(abs(_x0) + abs(_x1), static_cast<T>(1.0));
        return (static_cast<T>(2) * abs(_x1 - _x0) / sum) <= _epsilon;
    }

    // ********************************************************************* //
    /// \brief Round value towards negative infinity.
    template<typename T, class = std::enable_if_t<!std::is_base_of<details::NonScalarType, T>::value>>
    inline Sint<sizeof(T)> floor(T _x) noexcept
    {
        Sint<sizeof(T)> r = static_cast<Sint<sizeof(T)>>(_x);
        return r - static_cast<Sint<sizeof(T)>>((_x<static_cast<T>(0)) && (_x-r!=static_cast<T>(0)));
    }

    // ********************************************************************* //
    /// \brief Round value towards positive infinity.
    template<typename T, class = std::enable_if_t<!std::is_base_of<details::NonScalarType, T>::value>>
    inline Sint<sizeof(T)> ceil(T _x) noexcept
    {
        Sint<sizeof(T)> r = static_cast<Sint<sizeof(T)>>(_x);
        return r + static_cast<Sint<sizeof(T)>>((_x>static_cast<T>(0)) && (_x-r!=static_cast<T>(0)));
    }


    // ********************************************************************* //
    /// \brief Round value towards next integral number (0.5 rounds to even).
    template<typename T, class = std::enable_if_t<!std::is_base_of<details::NonScalarType, T>::value>>
    inline Sint<sizeof(T)> round(T _x) noexcept
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
    inline T frac(T _x) noexcept // TESTED
    {
        return _x - static_cast<Sint<sizeof(T)>>(_x);
    }

    // ********************************************************************* //
    /// \brief Divide a number into the integer and fractional part.
    /// \param _x [in] The number to be splitted.
    /// \param _int [out] The integer part of the number.
    /// \returns The fraction of the number in (-1,1).
    template<typename T>
    inline T intfrac(T _x, Sint<sizeof(T)>& _int) noexcept // TESTED
    {
        _int = static_cast<Sint<sizeof(T)>>(_x);
        return _x - _int;
    }

    // ********************************************************************* //
    /// \brief Divide a number into an integer and a positive fractional part.
    /// \details This method uses f-floor(f) instead of f-int(f) and therefore
    ///      has a continous behavior around zero
    /// \param _x [in] The number to be splitted.
    /// \param _int [out] The integer part of the number.
    /// \returns The fraction of the number in [0,1).
    template<typename T>
    inline T floorfrac(T _x, Sint<sizeof(T)>& _int) noexcept // TESTED
    {
        _int = floor(_x);
        return _x - _int;
    }

    // ********************************************************************* //
    /// \brief Get the smallest positive number m such that x=y*c+m with c in Z.
    /// \returns The mathematically defined positive modulus.
    template<typename T>
    inline T mod(T _x, T _y) noexcept // TESTED
    {
        eiAssert(_y != 0.0f, "Modulo 0 is not defined!");
        T m = fmod(_x, _y);
        return m < 0 ? m+abs(_y) : m;
    }

    // Pure integer specializations
    template<>
    inline int8 mod<int8>(int8 _x, int8 _y) noexcept
    {
        eiAssert(_y != 0, "Modulo 0 is not defined!");
        int8 m = _x % _y;
        return m < 0 ? m+abs(_y) : m;
    }
    template<>
    inline int16 mod<int16>(int16 _x, int16 _y) noexcept
    {
        eiAssert(_y != 0, "Modulo 0 is not defined!");
        int16 m = _x % _y;
        return m < 0 ? m+abs(_y) : m;
    }
    template<>
    inline int32 mod<int32>(int32 _x, int32 _y) noexcept
    {
        eiAssert(_y != 0, "Modulo 0 is not defined!");
        int32 m = _x % _y;
        return m < 0 ? m+abs(_y) : m;
    }
    template<>
    inline int64 mod<int64>(int64 _x, int64 _y) noexcept
    {
        eiAssert(_y != 0, "Modulo 0 is not defined!");
        int64 m = _x % _y;
        return m < 0 ? m+abs(_y) : m;
    }

    template<>
    inline uint8 mod<uint8>(uint8 _x, uint8 _y) noexcept { return _x % _y; }
    template<>
    inline uint16 mod<uint16>(uint16 _x, uint16 _y) noexcept { return _x % _y; }
    template<>
    inline uint32 mod<uint32>(uint32 _x, uint32 _y) noexcept { return _x % _y; }
    template<>
    inline uint64 mod<uint64>(uint64 _x, uint64 _y) noexcept { return _x % _y; }

    // ********************************************************************* //
    /// \brief Compute floor(log2), i.e. the position of the most significant bit.
    /// \details Expects IEEE and litte-endian (standard PC)
    template<typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
    constexpr int ilog2(T _x) noexcept                    // TESTED
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
    inline auto lerp(T0 _x0, T0 _x1, T1 _t) noexcept -> decltype(_x0*_t) // TESTED
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
    inline auto bilerp(T0 _x00, T0 _x01,
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
    inline T smoothstep(T _t) noexcept
    {
        return _t * _t * (T(3) - T(2) * _t);
    }

    // ********************************************************************* //
    /// \brief Evalutate the smootherstep polynomial 6t^5 - 15t^4 + 10t^3.
    /// \param _t [in] The value to be inserted into the polynomial. The useful
    ///    definition interval is in [0,1].
    template<typename T>
    inline T smootherstep(T _t) noexcept
    {
        return _t * _t * _t * (_t * (_t * T(6) - T(15)) + T(10));
    }

    // ********************************************************************* //
    /// \brief Compute the next machine representable number in positive direction.
    /// \details Iterates overall floating point numbers, even through
    ///    denormalized range. Thereby the successor of -0 and 0 are both the
    ///    same smallest positive float. Further -INF is transformed into
    ///    -FLOAT_MAX and +INF remains +INF.
    constexpr inline float successor(float _number) noexcept // TESTED
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

    constexpr inline double successor(double _number) noexcept
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
    constexpr inline float predecessor(float _number) noexcept // TESTED
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

    constexpr inline double predecessor(double _number) noexcept
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

    /// \brif Helper method to solve ax^2 + bx + c (numerically more stable than naive method).
    /// \returns Solutions x1 >= x2 (x2 is always the greater of the two results).
    inline bool solveSquarePoly(float a, float b, float c, float& x1, float& x2)
    {
        float discriminant = b*b - 4.0f*a*c;
        if(discriminant < 0.0f) return false;
        float dsqrt = sqrt(discriminant);
        if(b > 0) {
            x1 = (-b - dsqrt)/(2.0f*a);
            x2 = -2.0f*c/(b + dsqrt);
        } else {
            x1 = -2.0f*c/(b - dsqrt);
            x2 = (-b + dsqrt)/(2.0f*a);
        }
        return true;
    }
}