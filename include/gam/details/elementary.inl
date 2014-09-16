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
