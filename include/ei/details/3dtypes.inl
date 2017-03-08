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
    eiAssert( _box.max >= _box.min, "Invalid bounding box." );
}

inline Sphere::Sphere( const Vec3& _p0, const Vec3& _p1 ) :
    center((_p0 + _p1) * 0.5f),
    radius(len(_p0 - _p1) * 0.5f)
{
}


// ************************************************************************* //
inline Box::Box( const Vec3& _point ) :
    min(_point),
    max(_point)
{
}

template<typename... Args>
inline Box::Box( const Vec3& _point, Args... _morePoints ) :
    Box(_morePoints...)
{
    min = ei::min(min, _point);
    max = ei::max(max, _point);
}

inline Box::Box( const Box& _box0, const Box& _box1 ) :
    min(ei::min(_box0.min, _box1.min)),
    max(ei::max(_box0.max, _box1.max))
{
    eiAssert( max >= min,
        "Minimum coordinates must be smaller or equal the maximum." );
}

inline Box::Box( const Sphere& _sphere ) :
    min(_sphere.center - _sphere.radius),
    max(_sphere.center + _sphere.radius)
{
    eiAssertWeak( max >= min,
        "Subtraction or addition of a scalar failed or sphere had negative radius!" );
}

inline Box::Box( const Triangle& _triangle ) :
    min(ei::min(_triangle.v0, _triangle.v1, _triangle.v2)),
    max(ei::max(_triangle.v0, _triangle.v1, _triangle.v2))
{
    eiAssertWeak( max >= min,
        "min() or max() failed for a vector!" );
}

inline Box::Box( const Tetrahedron& _tetrahedron ) :
    min(ei::min(_tetrahedron.v0, _tetrahedron.v1, _tetrahedron.v2, _tetrahedron.v3)),
    max(ei::max(_tetrahedron.v0, _tetrahedron.v1, _tetrahedron.v2, _tetrahedron.v3))
{
    eiAssertWeak( max >= min,
        "min() or max() failed for a vector!" );
}

inline Box::Box( const Ellipsoid& _ellipsoid ) :
    min(_ellipsoid.center - _ellipsoid.radii),
    max(_ellipsoid.center + _ellipsoid.radii)
{
    eiAssertWeak(_ellipsoid.radii >= 0.0f, "Invalid ellipsoid!");
}


// ************************************************************************* //
inline OBox::OBox( const Vec3& _center, const Vec3& _halfSides, const Quaternion& _orientation ) :
    center(_center),
    halfSides(_halfSides),
    orientation(_orientation)
{
}

inline OBox::OBox( const Box& _box ) :
    center((_box.max + _box.min) * 0.5f),
    halfSides((_box.max - _box.min) * 0.5f),
    orientation(qidentity())
{
}

inline OBox::OBox( const Disc& _disc ) :
    center(_disc.center),
    halfSides(_disc.radius, _disc.radius, 0.0f),
    orientation( Vec3(0.0f, 0.0f, 1.0f), _disc.normal)
{
}


// ************************************************************************* //
inline Tetrahedron::Tetrahedron( const Vec3& _v0, const Vec3& _v1, const Vec3& _v2, const Vec3& _v3 ) :
    v0(_v0),
    v1(_v1),
    v2(_v2),
    v3(_v3)
{
}

inline Vec3& Tetrahedron::v(int _index)
{
    eiAssertWeak(_index >= 0 && _index < 4, "A thetrahedron only has 4 vertices!");
    return reinterpret_cast<Vec3*>(this)[_index];
}

inline const Vec3& Tetrahedron::v(int _index) const
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
    d0(max(_d0, _d1)),
    d1(min(_d0, _d1))
{
    eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
}

inline DOP::DOP(const Vec3& _normal, const Vec3& _support0, const Vec3& _support1) :
    DOP(_normal, -dot(_normal, _support0), -dot(_normal, _support1))
{
}


// ************************************************************************* //
inline Ellipsoid::Ellipsoid(const Vec3& _center, const Vec3& _radii) :
    center(_center),
    radii(max(_radii, Vec3(1e-16f)))
{
}

inline Ellipsoid::Ellipsoid(const Box& _box)
{
    eiAssert( _box.max >= _box.min, "Invalid bounding box." );
    center = (_box.max + _box.min) * 0.5f;
    /// sqrt(n) * side length / 2, where n is the number of dimensions with
    /// an extension (side length 0 allows to generate ellipses or rays)
    Vec3 sideLen = _box.max - _box.min;
    radii = (sqrt((float)sum(neq(sideLen, 0.0f))) * 0.5f) * sideLen;
    radii = max(radii, Vec3(1e-16f));
}

// ************************************************************************* //
inline OEllipsoid::OEllipsoid(const Vec3& _center, const Vec3& _radii, const Quaternion& _orientation) :
    center(_center),
    radii(_radii),
    orientation(_orientation)
{
}

inline OEllipsoid::OEllipsoid(const Box& _box)
{
    eiAssert( _box.max >= _box.min, "Invalid box." );
    center = (_box.max + _box.min) * 0.5f;
    /// sqrt(n) * side length / 2, where n is the number of dimensions with
    /// an extension (side length 0 allows to generate ellipses or rays)
    Vec3 sideLen = _box.max - _box.min;
    radii = (sqrt((float)sum(neq(sideLen, 0.0f))) * 0.5f) * sideLen;
    radii = max(radii, Vec3(1e-16f));
    orientation = qidentity();
}

inline OEllipsoid::OEllipsoid(const OBox& _box)
{
    eiAssert( _box.halfSides >= 0.0f, "Invalid box." );
    center = _box.center;
    /// sqrt(n) * side length / 2, where n is the number of dimensions with
    /// an extension (side length 0 allows to generate ellipses or rays)
    radii = sqrt((float)sum(neq(_box.halfSides, 0.0f))) * _box.halfSides;
    radii = max(radii, Vec3(1e-16f));
    orientation = _box.orientation;
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

inline FastFrustum& FastFrustum::operator = (const FastFrustum& _frustum)
{
    const_cast<DOP&>(nf) = _frustum.nf;
    const_cast<Plane&>(l) = _frustum.l;
    const_cast<Plane&>(r) = _frustum.r;
    const_cast<Plane&>(b) = _frustum.b;
    const_cast<Plane&>(t) = _frustum.t;
}


// ************************************************************************* //
// VOLUME AND SURFACE METHODS
// ************************************************************************* //

inline float volume(const Sphere& _sphere)
{
    return 4.0f / 3.0f * PI * _sphere.radius * _sphere.radius * _sphere.radius;
}

inline float volume(const Box& _box)
{
    Vec3 size = _box.max - _box.min;
    return size.x * size.y * size.z;
}

inline float volume(const OBox& _obox)
{
    return 8.0f * _obox.halfSides.x * _obox.halfSides.y * _obox.halfSides.z;
}

inline float volume(const Tetrahedron& _thetrahedron)
{
    return dot(_thetrahedron.v3-_thetrahedron.v0, cross(_thetrahedron.v2-_thetrahedron.v0, _thetrahedron.v1-_thetrahedron.v0)) / 6.0f;
}

inline float volume(const Triangle&)
{
    return 0.0f;
}

inline float volume(const Disc&)
{
    return 0.0f;
}

inline float volume(const Plane&)
{
    return 0.0f;
}

inline float volume(const DOP&)
{
    return 0.0f;
}

inline float volume(const Ellipsoid& _ellipsoid)
{
    return 4.0f / 3.0f * PI * _ellipsoid.radii.x * _ellipsoid.radii.y * _ellipsoid.radii.z;
}

inline float volume(const OEllipsoid& _oellipsoid)
{
    return 4.0f / 3.0f * PI * _oellipsoid.radii.x * _oellipsoid.radii.y * _oellipsoid.radii.z;
}

inline float volume(const Ray&)
{
    return 0.0f;
}

inline float volume(const Segment&)
{
    return 0.0f;
}

inline float volume(const Capsule& _capsule)
{
    return PI * sq(_capsule.radius) * (_capsule.radius*4.0f/3.0f + len(_capsule.seg.b-_capsule.seg.a));
}

// ************************************************************************* //
inline float surface(const Sphere& _sphere)
{
    return 4.0f * PI * sq(_sphere.radius);
}

inline float surface(const Box& _box)
{
    Vec3 size = _box.max - _box.min;
    return 2.0f * (size.x * size.y + size.x * size.z + size.y * size.z);
}

inline float surface(const OBox& _obox)
{
    return 8.0f * (_obox.halfSides.x * _obox.halfSides.y + _obox.halfSides.x * _obox.halfSides.z + _obox.halfSides.y * _obox.halfSides.z);
}

inline float surface(const Tetrahedron& _thetra)
{
    // Analogous to a triangle (repeated four times)
    Vec3 a = _thetra.v1 - _thetra.v0;
    Vec3 b = _thetra.v2 - _thetra.v0;
    Vec3 c = _thetra.v3 - _thetra.v0;
    Vec3 d = _thetra.v2 - _thetra.v1;
    Vec3 e = _thetra.v3 - _thetra.v1;
    return 0.5f * (len( cross(a, b) )
                 + len( cross(a, c) )
                 + len( cross(b, c) )
                 + len( cross(d, e) ));
}

inline float surface(const Triangle& _triangle)
{
    // Heron's formula is much more expensive than cross product because
    // the 3 side lengths must be computed first.
    return len( cross(_triangle.v1 - _triangle.v0, _triangle.v2 - _triangle.v0) ) * 0.5f;
}

inline float surface(const Disc& _disc)
{
    return PI * _disc.radius;
}

inline float surface(const Plane&)
{
    return INF;
}

inline float surface(const DOP&)
{
    return INF;
}

inline float surface(const Ellipsoid& _ellipsoid)
{
    // Use approximation (Knud Thomsen's formula) only! Everything else is a
    // lot larger.
    Vec3 pr( pow(_ellipsoid.radii.x, 1.6075f),
             pow(_ellipsoid.radii.y, 1.6075f),
             pow(_ellipsoid.radii.z, 1.6075f) );
    return 4.0f * PI * pow((pr.x * pr.y + pr.x * pr.z + pr.y * pr.z) / 3.0f, 0.622083981f);
}

inline float surface(const OEllipsoid& _oellipsoid)
{
    // Use approximation (Knud Thomsen's formula) only! Everything else is a
    // lot larger.
    Vec3 pr( pow(_oellipsoid.radii.x, 1.6075f),
             pow(_oellipsoid.radii.y, 1.6075f),
             pow(_oellipsoid.radii.z, 1.6075f) );
    return 4.0f * PI * pow((pr.x * pr.y + pr.x * pr.z + pr.y * pr.z) / 3.0f, 0.622083981f);
}

inline float surface(const Ray&)
{
    return 0.0f;
}

inline float surface(const Segment&)
{
    return 0.0f;
}

inline float surface(const Capsule& _capsule)
{
    return 2 * PI * _capsule.radius * (2 * _capsule.radius + len(_capsule.seg.b-_capsule.seg.a));
}


// ************************************************************************* //
// CENTER METHODS
// ************************************************************************* //
inline Vec3 center(const Sphere& _sphere)
{
    return _sphere.center;
}

inline Vec3 center(const Box& _box)
{
    return (_box.min + _box.max) * 0.5f;
}

inline Vec3 center(const OBox& _obox)
{
    return _obox.center;
}

inline Vec3 center(const Tetrahedron& _thetrahedron)
{
    return (_thetrahedron.v0 + _thetrahedron.v1 + _thetrahedron.v2 + _thetrahedron.v3) / 4.0f;
}

inline Vec3 center(const Triangle& _triangle)
{
    return (_triangle.v0 + _triangle.v1 + _triangle.v2) / 3.0f;
}

inline Vec3 center(const Disc& _disc)
{
    return _disc.center;
}

inline Vec3 center(const Ellipsoid& _ellipsoid)
{
    return _ellipsoid.center;
}

inline Vec3 center(const Segment& _line)
{
    return (_line.a + _line.b) * 0.5f;
}

inline Vec3 center(const Capsule& _capsule)
{
    return (_capsule.seg.a + _capsule.seg.b) * 0.5f;
}