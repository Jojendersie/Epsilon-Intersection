// ************************************************************************* //
inline Sphere::Sphere( const Vec3& _center, float _radius ) :
    center(_center),
    radius(_radius)
{
}

inline Sphere::Sphere( const Box& _box ) :
    center((_box.min + _box.max) * 0.5f),
    radius(len(_box.max - _box.min) * 0.5f)
{
    eiAssert( all(_box.max >= _box.min), "Invalid bounding box." );
}


// ************************************************************************* //
inline Box::Box( const Vec3& _min, const Vec3& _max ) :
    min(_min),
    max(_max)
{
    eiAssert( all(_max >= _min),
        "Minimum coordinates must be smaller or equal the maximum." );
}

inline Box::Box( const Box& _box0, const Box& _box1 ) :
    min(ei::min(_box0.min, _box1.min)),
    max(ei::max(_box0.max, _box1.max))
{
    eiAssert( all(max >= min),
        "Minimum coordinates must be smaller or equal the maximum." );
}

inline Box::Box( const Sphere& _sphere ) :
    min(_sphere.center - _sphere.radius),
    max(_sphere.center + _sphere.radius)
{
    eiAssertWeak( all(max >= min),
        "Subtraction or addition of a scalar failed or sphere had negative radius!" );
}

inline Box::Box( const Triangle& _triangle ) :
    min(ei::min(_triangle.v0, ei::min(_triangle.v1, _triangle.v2))),
    max(ei::max(_triangle.v0, ei::max(_triangle.v1, _triangle.v2)))
{
    eiAssertWeak( all(max >= min),
        "min() or max() failed for a vector!" );
}

inline Box::Box( const Ellipsoid& _ellipsoid ) :
    min(_ellipsoid.center - _ellipsoid.radii),
    max(_ellipsoid.center + _ellipsoid.radii)
{
    eiAssertWeak(all(_ellipsoid.radii >= 0.0f), "Invalid ellipsoid!");
}


// ************************************************************************* //
inline Thetrahedron::Thetrahedron( const Vec3& _v0, const Vec3& _v1, const Vec3& _v2, const Vec3& _v3 ) :
    v0(_v0),
    v1(_v1),
    v2(_v2),
    v3(_v3)
{
}

inline Vec3& Thetrahedron::v(int _index)
{
    eiAssertWeak(_index >= 0 && _index < 4, "A thetrahedron only has 4 vertices!");
    return reinterpret_cast<Vec3*>(this)[_index];
}

inline const Vec3& Thetrahedron::v(int _index) const
{
    eiAssertWeak(_index >= 0 && _index < 4, "A thetrahedron only has 3 vertices!");
    return reinterpret_cast<const Vec3*>(this)[_index];
}


// ************************************************************************* //
inline Triangle::Triangle( const Vec3& _v0, const Vec3& _v1, const Vec3& _v2 ) :
    v0(_v0),
    v1(_v1),
    v2(_v2)
{
}

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
inline Disc::Disc(const Vec3& _center, const Vec3& _normal, float _radius) :
    center(_center),
    normal(_normal),
    radius(_radius)
{
    eiAssert(_radius >= 0.0f, "Expected a positive radius!");
}


// ************************************************************************* //
inline Plane::Plane(const Vec3& _normal, float _d) :
    n(_normal),
    d(_d)
{
    eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
}

inline Plane::Plane(const Vec3& _normal, const Vec3& _support) :
    n(_normal),
    d(-dot(_normal, _support))
{
    eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
}

inline Plane::Plane(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2)
{
    n = normalize(cross(_v1 - _v0, _v2 - _v0));
    d = -dot(n, _v0);
}


// ************************************************************************* //
inline DOP::DOP(const Vec3& _normal, float _d0, float _d1) :
    n(_normal),
    d0(_d0),
    d1(_d1)
{
    eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
}

inline DOP::DOP(const Vec3& _normal, const Vec3& _support0, const Vec3& _support1) :
    n(_normal),
    d0(-dot(_normal, _support0)),
    d1(-dot(_normal, _support1))
{
    eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
}


// ************************************************************************* //
inline Ellipsoid::Ellipsoid(const Vec3& _center, const Vec3& _radii) :
    center(_center),
    radii(max(_radii, Vec3(1e-16f)))
{
}

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


// ************************************************************************* //
inline Segment::Segment(const Vec3& _a, const Vec3& _b) :
    a(_a),
    b(_b)
{
}

inline Segment::Segment(const Ray& _ray, float _distance) :
    a(_ray.origin),
    b(_ray.origin + _ray.direction * _distance)
{
    eiAssertWeak(approx(lensq(_ray.direction), 1.0f), "The input ray is not normalized!");
}


// ************************************************************************* //
inline Capsule::Capsule(const Vec3& _a, const Vec3& _b, float _radius) :
    seg(_a, _b),
    radius(_radius)
{
    eiAssertWeak(_radius >= 0.0f, "Radius must be positive!");
}

inline Capsule::Capsule(const Segment& _line, float _radius) :
    seg(_line),
    radius(_radius)
{
    eiAssertWeak(_radius >= 0.0f, "Radius must be positive!");
}


// ************************************************************************* //
inline Frustum::Frustum(const Vec3& _apex, const Vec3& _direction, const Vec3& _up, float _l, float _r, float _b, float _t, float _n, float _f) :
    l(_l), r(_r), b(_b), t(_t), n(_n), f(_f),
    apex(_apex),
    up(_up),
    direction(_direction)
{
    eiAssert(approx(lensq(_direction), 1.0f), "Insert a normalized direction!");
    eiAssert(approx(lensq(_up), 1.0f), "Insert a normalized up vector!");
    eiAssert(_n < _f && 0 <= _n, "Near and far frustum planes are sorted wrongly.");
    eiAssert(_l < _r, "Left and right frustum planes are sorted wrongly.");
    eiAssert(_b < _t, "Top and bottom frustum planes are sorted wrongly.");
}


// ************************************************************************* //
inline FastFrustum::FastFrustum(const Vec3& _apex, const Vec3& _direction, const Vec3& _up, float _l, float _r, float _b, float _t, float _n, float _f) :
    FastFrustum(Frustum(_apex, _direction, _up, _l, _r, _b, _t, _n, _f))
{
}

inline Frustum& FastFrustum::operator = (const Frustum& _frustum)
{
    const_cast<DOP&>(nf) = nf;
    const_cast<Plane&>(l) = l;
    const_cast<Plane&>(r) = r;
    const_cast<Plane&>(b) = b;
    const_cast<Plane&>(t) = t;
}