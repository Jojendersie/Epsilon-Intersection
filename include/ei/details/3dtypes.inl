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
    eiAssert( all(_max >= _min),
        "Minimum coordinates must be smaller or equal the maximum." );
}

// ************************************************************************* //
inline Box::Box( const Box& _box0, const Box& _box1 ) :
    min(ei::min(_box0.min, _box1.min)),
    max(ei::max(_box0.max, _box1.max))
{
    eiAssert( all(max >= min),
        "Minimum coordinates must be smaller or equal the maximum." );
}

// ************************************************************************* //
inline Box::Box( const Sphere& _sphere ) :
    min(_sphere.center - _sphere.radius),
    max(_sphere.center + _sphere.radius)
{
    eiAssertWeak( all(max >= min),
        "Subtraction or addition of a scalar failed or sphere had negative radius!" );
}

// ************************************************************************* //
inline Box::Box( const Triangle& _triangle ) :
    min(ei::min(_triangle.v0, ei::min(_triangle.v1, _triangle.v2))),
    max(ei::max(_triangle.v0, ei::max(_triangle.v1, _triangle.v2)))
{
    eiAssertWeak( all(max >= min),
        "min() or max() failed for a vector!" );
}

// ************************************************************************* //
inline Triangle::Triangle( const Vec3& _v0, const Vec3& _v1, const Vec3& _v2 ) :
    v0(_v0),
    v1(_v1),
    v2(_v2)
{
}

// ************************************************************************* //
inline Vec3& Triangle::v(int _index)
{
    eiAssertWeak(_index >= 0 && _index < 3, "A triangle only has 3 vertices!");
    return reinterpret_cast<Vec3*>(this)[_index];
}

inline const Vec3& Triangle::v(int _index) const
{
    eiAssertWeak(_index >= 0 && _index < 3, "A triangle only has 3 vertices!");
    return reinterpret_cast<const Vec3*>(this)[_index];
}

// ************************************************************************* //
inline Plane::Plane(const Vec3& _normal, float _d) :
    n(_normal),
    d(_d)
{
    eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
}

// ************************************************************************* //
inline Plane::Plane(const Vec3& _normal, const Vec3& _support) :
    n(_normal),
    d(-dot(_normal, _support))
{
    eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
}

// ************************************************************************* //
inline Plane::Plane(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2)
{
    n = normalize(cross(_v1 - _v0, _v2 - _v0));
    d = -dot(n, _v0);
}


// ************************************************************************* //
inline Ellipsoid::Ellipsoid(const Vec3& _center, const Vec3& _radii) :
    center(_center),
    radii(max(_radii, Vec3(1e-16f)))
{
}

// ************************************************************************* //
inline Ellipsoid::Ellipsoid(const Box& _box)
{
    eiAssert( all(_box.max >= _box.min), "Invalid bounding box." );
    center = (_box.max + _box.min) * 0.5f;
    /// sqrt(n) * side length / 2, where n is the number of dimensions with
    /// an extension (side length 0 allows to generate ellipses or rays)
    Vec3 sideLen = _box.max - _box.min;
    radii = (sqrt((float)sum(sideLen != 0.0f)) * 0.5f) * sideLen;
    radii = max(radii, Vec3(1e-16f));
}

// ************************************************************************* //
inline Ray::Ray(const Vec3& _origin, const Vec3& _direction) :
    origin(_origin),
    direction(_direction)
{
    eiAssert(approx(lensq(_direction), 1.0f), "Insert a normalized normal!");
}