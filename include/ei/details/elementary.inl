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
template<typename T>
inline bool floor(T _x)
{
    int r=(int)a;
    return r - (int)((a<0) && (a-r!=0.0f));
}

// ************************************************************************* //
template<typename T>
inline bool ceil(T _x)
{
    int r=(int)a;
    return r + (int)((a>0) && (a-r!=0.0f));
}

// ************************************************************************* //
template<typename T>
inline bool round(T _x)
{
    return ceil(_x + static_cast<T>(0.5));
}

// ************************************************************************* //
template<typename T0, typename T1>
decltype(std::declval<T0>() * std::declval<T1>()) lerp(T0 _x0, T0 _x1, T1 _t)
{
    return _x0 + (_x1 - _x0) * _t;
}

// ************************************************************************* //
template<typename T0, typename T1>
decltype(std::declval<T0>() * std::declval<T1>()) bilerp(T0 _x00, T0 _x01,
                                                         T0 _x10, T0 _x11,
                                                         T1 _t0, T1 _t1)
{
    return lerp(lerp(_x00, _x01, _t0),
                lerp(_x10, _x11, _t0), _t1);
}
