// ************************************************************************* //
template<typename T>
inline T sq(T _x)
{
    return _x * _x;
}

// ************************************************************************* //
template<typename T>
inline T max(T _x, T _y)
{
    return _x < _y ? _y : _x;
}

// ************************************************************************* //
template<typename T>
inline T min(T _x, T _y)
{
    return _x > _y ? _y : _x;
}

// ************************************************************************* //
template<typename T>
inline T clamp(T _x, T _min, T _max)
{
    return _x > _max ? _max : (_x < _min ? _min : _x);
}

// ************************************************************************* //
template<typename T>
inline T saturate(T _x)
{
    return _x > static_cast<T>(1) ? static_cast<T>(1) : (_x < static_cast<T>(0) ? static_cast<T>(0) : _x);
}

// ************************************************************************* //
template<typename T>
inline T abs(T _x)
{
    return _x < static_cast<T>(0) ? -_x : _x;
}

// ************************************************************************* //
template<typename T>
inline T sign(T _x)
{
    return _x < static_cast<T>(0) ? static_cast<T>(-1)
        : (_x > static_cast<T>(0) ? static_cast<T>(1) : static_cast<T>(0));
}

// ************************************************************************* //
template<typename T>
inline T sgn(T _x)
{
    // Use explicit construction T() or casting
    return _x < static_cast<T>(0) ? static_cast<T>(-1) : static_cast<T>(1);
}

// ************************************************************************* //
template<typename T>
inline bool approx(T _x0, T _x1, float _epsilon)
{
    return abs(_x1 - _x0) <= _epsilon;
}

// ************************************************************************* //
template<typename T, class>
inline typename details::Int<sizeof(T)>::stype floor(T _x)
{
    typename details::Int<sizeof(T)>::stype r = static_cast<typename details::Int<sizeof(T)>::stype>(_x);
    return r - static_cast<typename details::Int<sizeof(T)>::stype>((_x<static_cast<T>(0)) && (_x-r!=static_cast<T>(0)));
}

// ************************************************************************* //
template<typename T, class>
inline typename details::Int<sizeof(T)>::stype ceil(T _x)
{
    typename details::Int<sizeof(T)>::stype r = static_cast<typename details::Int<sizeof(T)>::stype>(_x);
    return r + static_cast<typename details::Int<sizeof(T)>::stype>((_x>static_cast<T>(0)) && (_x-r!=static_cast<T>(0)));
}

// ************************************************************************* //
template<typename T, class>
inline typename details::Int<sizeof(T)>::stype round(T _x)
{
    // Round up
    //return floor(_x + static_cast<T>(0.5));
    // Round to even
    typename details::Int<sizeof(T)>::stype r = static_cast<typename details::Int<sizeof(T)>::stype>(_x);
    r -= static_cast<typename details::Int<sizeof(T)>::stype>((_x<static_cast<T>(0)) && (_x!=static_cast<T>(r)));   // Subtract 1 if _x < 0 -> r is floor(_x)
    T f = _x - r;    // Fractional part (positive)
    if(f < static_cast<T>(0.5)) return r;
    if(f > static_cast<T>(0.5)) return r+1;
    // f is 0.5
    if(r & 1) return r + 1;
    return r;
}

// ********************************************************************* //
template<typename T>
T frac(T _x)
{
    return _x - static_cast<typename details::Int<sizeof(T)>::stype>(_x);
}

// ********************************************************************* //
template<typename T>
T intfrac(T _x, typename details::Int<sizeof(T)>::stype& _int)
{
    _int = static_cast<typename details::Int<sizeof(T)>::stype>(_x);
    return _x - _int;
}

// ********************************************************************* //
template<typename T>
T floorfrac(T _x, typename details::Int<sizeof(T)>::stype& _int)
{
    _int = floor(_x);
    return _x - _int;
}

// ********************************************************************* //
template<typename T>
T mod(T _x, T _y)
{
    eiAssert(_y != 0.0f, "Modulu 0 is not defined!");
    T m = fmod(_x, _y);
    return m < 0 ? m+abs(_y) : m;
}

// Pure integer specializations
template<>
inline int8 mod<int8>(int8 _x, int8 _y)
{
    eiAssert(_y != 0, "Modulu 0 is not defined!");
    int8 m = _x % _y;
    return m < 0 ? m+abs(_y) : m;
}
template<>
inline int16 mod<int16>(int16 _x, int16 _y)
{
    eiAssert(_y != 0, "Modulu 0 is not defined!");
    int16 m = _x % _y;
    return m < 0 ? m+abs(_y) : m;
}
template<>
inline int32 mod<int32>(int32 _x, int32 _y)
{
    eiAssert(_y != 0, "Modulu 0 is not defined!");
    int32 m = _x % _y;
    return m < 0 ? m+abs(_y) : m;
}
template<>
inline int64 mod<int64>(int64 _x, int64 _y)
{
    eiAssert(_y != 0, "Modulu 0 is not defined!");
    int64 m = _x % _y;
    return m < 0 ? m+abs(_y) : m;
}

template<>
inline uint8 mod<uint8>(uint8 _x, uint8 _y) { return _x % _y; }
template<>
inline uint16 mod<uint16>(uint16 _x, uint16 _y) { return _x % _y; }
template<>
inline uint32 mod<uint32>(uint32 _x, uint32 _y) { return _x % _y; }
template<>
inline uint64 mod<uint64>(uint64 _x, uint64 _y) { return _x % _y; }

// ************************************************************************* //
template<typename T0, typename T1>
auto lerp(T0 _x0, T0 _x1, T1 _t) -> decltype(_x0*_t)
{
    return _x0 + (_x1 - _x0) * _t;
}

// ************************************************************************* //
template<typename T0, typename T1>
auto bilerp(T0 _x00, T0 _x01,
            T0 _x10, T0 _x11,
            T1 _t0, T1 _t1) -> decltype(_x00*_t0)
{
    return lerp(lerp(_x00, _x01, _t0),
                lerp(_x10, _x11, _t0), _t1);
}
