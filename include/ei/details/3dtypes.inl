// ************************************************************************* //
inline Sphere::Sphere( const Vec3& _center, float _radius ) :
    center(_center),
    radius(_radius)
{
}

// ************************************************************************* //
inline Box::Box( const Vec3& _min, const Vec3& _max ) :
    min(_min),
    max(_max)
{
    assertlvl1( all(_max >= _min),
        "Minimum coordinates must be smaller or equal the maximum." );
}

// ************************************************************************* //
inline Box::Box( const Box& _box0, const Box& _box1 ) :
    min(ei::min(_box0.min, _box1.min)),
    max(ei::max(_box0.max, _box1.max))
{
    assertlvl1( all(max >= min),
        "Minimum coordinates must be smaller or equal the maximum." );
}