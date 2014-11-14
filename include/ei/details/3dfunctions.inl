// ************************************************************************* //
inline float volume( const Sphere& _sphere)
{
    return 4.0f / 3.0f * PI * _sphere.radius * _sphere.radius * _sphere.radius;
}

inline float volume( const Box& _box)
{
    Vec3 size = _box.max - _box.min;
    return size.x * size.y * size.z;
}

inline float volume( const Triangle& _triangle)
{
    return 0.0f;
}

inline float volume( const Plane& _plane)
{
    return 0.0f;
}

// ************************************************************************* //
inline float surface( const Sphere& _sphere)
{
    return 4.0f * PI * sq(_sphere.radius);
}

inline float surface( const Box& _box)
{
    Vec3 size = _box.max - _box.min;
    return 2.0f * (size.x * size.y + size.x * size.z + size.y * size.z);
}

inline float surface( const Triangle& _triangle)
{
    // Heron's formula is much more expensive than cross product because
    // the 3 side lengths must be computed first.
    return len( cross(_triangle.v1 - _triangle.v0, _triangle.v2 - _triangle.v0) ) * 0.5f;
}

inline float surface( const Plane& _plane)
{
    // Avoid including <limits> by defining infinity itself.
    union {
        float f;
        uint32 i;
    } inf;
    inf.i = 0x7f800000;
    return inf.f;
}