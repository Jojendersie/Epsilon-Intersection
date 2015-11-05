#pragma once

#include "configure.hpp"
#include "details/inttemplate.hpp"

// Ugly to include this here but you will miss sqrt,sin,... otherwise
#include <cmath>

namespace eitypes {
    /// \brief Short name for unsigned / unsigned int.
    /// \details The number of bit might vary on different systems.
    typedef unsigned int uint;                                                 // TESTED

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
}

#ifdef EI_GLOBAL_ELEMENTARIES
    using namespace eitypes;
#endif
namespace ei {
    using namespace eitypes;
}

namespace details {
    /// \brief Dummy class to detect correct types for matrix <-> matrix
    ///     and matrix <-> scalar operations (and the same for quaternions).
    /// \details The overloading mechanism fails when both types of operations
    ///     are templated, because the matrix <-> scalar is chosen even if
    ///     both operants are matrices.
    class NonScalarType    {};

    // Avoid including <limits> by defining infinity itself.
    union ReinterpretFloat {
        float f;
        Int<4>::utype i;
        ReinterpretFloat(Int<4>::utype _i) : i(_i) {}
        ReinterpretFloat(float _f) : f(_f) {}
    };
    const ReinterpretFloat F_INF = 0x7f800000u;

    union ReinterpretDouble {
        double f;
        Int<8>::utype i;
        ReinterpretDouble(Int<8>::utype _i) : i(_i) {}
        ReinterpretDouble(double _f) : f(_f) {}
    };
    const ReinterpretDouble D_INF = 0x7ff0000000000000ul;
}

namespace ei {
    // ********************************************************************* //
    //                             MATH CONSTANTS                            //
    // ********************************************************************* //
    const float PI = 3.141592654f;
    const float E = 2.718281828f;
    const float GOLDEN_RATIO = 1.61803398875f;
    const float PHYTAGORAS = 1.4142135623f;
    // The cmath header has an ugly macro with name INFINITY -> name conflict
    const float INF = details::F_INF.f;
    const double INF_D = details::D_INF.f;
    // Unicode names for the above constants
#ifdef EI_USE_UNICODE_NAMES
    const float π = PI;
    const float Φ = GOLDEN_RATIO;
    const float √2 = PHYTAGORAS;
    const float ℇ = E;
    const float ∞ = INF;
#endif

    // ********************************************************************* //
    //                               FUNCTIONS                               //
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
    T max(T _x, T _y);                                                         // TESTED
    /// \brief Get the maximum of any number of arguments
    /// \details In case of equal arguments the left most one is returned
    template<typename T, typename... Ttail>
    T max(T _first, Ttail... _tail);                                           // TESTED

    // ********************************************************************* //
    /// \brief Get the minimum from x and y.
    /// \details In case of x and y are equal, x is returned. This only makes
    ///    a difference if you are sorting object types with more than the
    ///    compared value.
    template<typename T>
    T min(T _x, T _y);                                                         // TESTED
    /// \brief Get the maximum of any number of arguments
    /// \details In case of equal arguments the left most one is returned
    template<typename T, typename... Ttail>
    T min(T _first, Ttail... _tail);                                           // TESTED

    // ********************************************************************* //
    /// \brief Clamp a value to the boundaries.
    template<typename T>
    T clamp(T _x, T _min, T _max);                                             // TESTED

    // ********************************************************************* //
    /// \brief Clamp a value to [0,1] interval.
    template<typename T>
    T saturate(T _x);

    // ********************************************************************* //
    /// \brief Get the absolute value.
    template<typename T>
    T abs(T _x);                                                               // TESTED

    // ********************************************************************* //
    /// \brief Get the sign of a value.
    /// \details There is a faster version sgn(), if you don't need to 
    ///    know about zero.
    /// \returns -1 (_x < 0), 0 (_x == 0) or 1 (_x > 0)
    template<typename T>
    T sign(T _x);                                                              // TESTED

    // ********************************************************************* //
    /// \brief Get the sign of a value where 0 is counted as positive.
    /// \details This function is faster than sign(). Use it if you don't need
    ///    to know about zero.
    /// \returns -1 (_x < 0) or 1 (_x >= 0)
    template<typename T>
    T sgn(T _x);                                                               // TESTED

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
    template<typename T, class = typename std::enable_if<!std::is_base_of<details::NonScalarType, T>::value, class Dummy>::type>
    bool approx(T _x0, T _x1, T _epsilon = T(1e-6));                           // TESTED

    // ********************************************************************* //
    /// \brief Round value towards negative infinity.
    template<typename T, class = typename std::enable_if<!std::is_base_of<details::NonScalarType, T>::value, class Dummy>::type>
    typename details::Int<sizeof(T)>::stype floor(T _x);

    // ********************************************************************* //
    /// \brief Round value towards positive infinity.
    template<typename T, class = typename std::enable_if<!std::is_base_of<details::NonScalarType, T>::value, class Dummy>::type>
    typename details::Int<sizeof(T)>::stype ceil(T _x);

    // ********************************************************************* //
    /// \brief Round value towards next integral number (0.5 rounds up).
    template<typename T, class = typename std::enable_if<!std::is_base_of<details::NonScalarType, T>::value, class Dummy>::type>
    typename details::Int<sizeof(T)>::stype round(T _x);

    // ********************************************************************* //
    /// \brief Get the fraction in (-1,1) with f-int(f).
    /// \param _x [in] The number to be splitted.
    template<typename T>
    T frac(T _x);                                                              // TESTED

    // ********************************************************************* //
    /// \brief Divide a number into the integer and fractional part.
    /// \param _x [in] The number to be splitted.
    /// \param _int [out] The integer part of the number.
    /// \returns The fraction of the number in (-1,1).
    template<typename T>
    T intfrac(T _x, typename details::Int<sizeof(T)>::stype& _int);            // TESTED

    // ********************************************************************* //
    /// \brief Divide a number into an integer and a positive fractional part.
    /// \details This method uses f-floor(f) instead of f-int(f) and therefore
    ///      has a continous behavior around zero
    /// \param _x [in] The number to be splitted.
    /// \param _int [out] The integer part of the number.
    /// \returns The fraction of the number in [0,1).
    template<typename T>
    T floorfrac(T _x, typename details::Int<sizeof(T)>::stype& _int);          // TESTED

    // ********************************************************************* //
    /// \brief Get the smallest positive number m such that x=y*c+m with c in Z.
    /// \returns The mathematically defined positive modulus.
    template<typename T>
    T mod(T _x, T _y);                                                         // TESTED

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
    auto lerp(T0 _x0, T0 _x1, T1 _t) -> decltype(_x0*_t);                      // TESTED

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
    auto bilerp(T0 _x00, T0 _x01,
                T0 _x10, T0 _x11,
                T1 _t0, T1 _t1) -> decltype(_x00*_t0);                         // TESTED

    // ********************************************************************* //
    /// \brief Evalutate the smoothstep polynomial 3t^2 - 2t^3.
    /// \param _t [in] The value to be inserted into the polynomial. The useful
    ///    definition interval is in [0,1].
    template<typename T>
    T smoothstep(T _t);

    // ********************************************************************* //
    /// \brief Evalutate the smootherstep polynomial 6t^5 - 15t^4 + 10t^3.
    /// \param _t [in] The value to be inserted into the polynomial. The useful
    ///    definition interval is in [0,1].
    template<typename T>
    T smootherstep(T _t);

    // ********************************************************************* //
    /// \brief Compute the next machine representable number in positive direction.
    /// \details Iterates overall floating point numbers, even through
    ///    denormalized range. Thereby the successor of -0 and 0 are both the
    ///    same smallest positive float. Further -INF is transformed into
    ///    -FLOAT_MAX and +INF remains +INF.
    float successor(float _number);                                            // TESTED
    double successor(double _number);

    /// \brief Compute the next machine representable number in negative direction.
    /// \details Iterates overall floating point numbers, even through
    ///    denormalized range. Thereby the predecessor of -0 and 0 are both the
    ///    same smallest negative float. Further -INF remains -INF and + INF
    ///    is transformed into FLOAT_MAX.
    float predecessor(float _number);                                          // TESTED
    double predecessor(double _number);

    // Include implementation.
#   include "details/elementary.inl"
}