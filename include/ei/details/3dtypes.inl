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