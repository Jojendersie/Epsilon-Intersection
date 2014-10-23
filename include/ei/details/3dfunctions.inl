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