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
    return _obox.sides.x * _obox.sides.y * _obox.sides.z;
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
    return 2.0f * (_obox.sides.x * _obox.sides.y + _obox.sides.x * _obox.sides.z + _obox.sides.y * _obox.sides.z);
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
inline float distance(const Vec3& _point0, const Vec3& _point1)
{
    return len(_point1 - _point0);
}

// ************************************************************************* //
inline float distance(const Vec3& _point, const Sphere& _sphere)
{
    return distance(_sphere.center, _point) - _sphere.radius;
}

// ************************************************************************* //
inline float distance(const Vec3& _point, const Capsule& _capsule)
{
    return distance(_point, _capsule.seg) - _capsule.radius;
}

// ************************************************************************* //
inline float distance(const Vec3& _point, const Plane& _plane)
{
    eiAssert(approx(1.0f, len(_plane.n)), "The plane is not normalized!");
    return dot(_plane.n, _point) + _plane.d;
}

// ************************************************************************* //
inline float distance(const Vec3& _point, const DOP& _dop)
{
    eiAssert(approx(1.0f, len(_dop.n)), "The plane is not normalized!");
    float d = dot(_dop.n, _point);
    if(d < -_dop.d0) return d + _dop.d0;
    if(d > -_dop.d1) return d + _dop.d1;
    return 0.0f;
    // There are three cases which all result in the same expression.
    // if(d < -_dop.d0) return -d - _dop.d0;
    // if(d > -_dop.d1) return d + _dop.d1;
    // Point is between both planes, return the shorter (negative) distance
    //return max(-d - _dop.d0, d + _dop.d1);
}

// ************************************************************************* //
inline float distance(const Sphere& _sphere, const Plane& _plane)
{
    eiAssert(approx(1.0f, len(_plane.n)), "The plane is not normalized!");
    eiAssert(_sphere.radius >= 0.0f, "Invalid sphere!");
    float d = dot(_plane.n, _sphere.center) + _plane.d;
    // Keep sign for the plane distances
    if(d < 0.0f) return min(d + _sphere.radius, 0.0f);
    else return max(d - _sphere.radius, 0.0f);
}

// ************************************************************************* //
inline float distance(const Sphere& _sphere, const Segment& _segment)
{
    return max(0.0f, distance(_sphere.center, _segment) - _sphere.radius);
}

// ************************************************************************* //
inline float distance(const Sphere& _sphere, const Capsule& _capsule)
{
    return max(0.0f, distance(_sphere.center, _capsule.seg) - _capsule.radius - _sphere.radius);
}

// ************************************************************************* //
inline float distance(const Capsule& _capsule0, const Capsule& _capsule1)
{
    return max(0.0f, distance(_capsule0.seg, _capsule1.seg) - _capsule0.radius - _capsule1.radius);
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