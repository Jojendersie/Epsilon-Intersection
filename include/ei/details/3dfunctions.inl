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

inline float volume(const Thetrahedron& _thetrahedron)
{
    return dot(_thetrahedron.v3-_thetrahedron.v0, cross(_thetrahedron.v2-_thetrahedron.v0, _thetrahedron.v1-_thetrahedron.v0)) / 6.0f;
}

inline float volume(const Triangle& _triangle)
{
    return 0.0f;
}

inline float volume(const Disc& _disc)
{
    return 0.0f;
}

inline float volume(const Plane& _plane)
{
    return 0.0f;
}

inline float volume(const DOP& _dop)
{
    return 0.0f;
}

inline float volume(const Ellipsoid& _ellipsoid)
{
    return 4.0f / 3.0f * PI * _ellipsoid.radii.x * _ellipsoid.radii.y * _ellipsoid.radii.z;
}

inline float volume(const Ray& _ray)
{
    return 0.0f;
}

inline float volume(const Line& _line)
{
    return 0.0f;
}

inline float volume(const Capsule& _capsule)
{
    return PI * sq(_capsule.radius) * (_capsule.radius*4.0f/3.0f + len(_capsule.b-_capsule.a));
}

inline float volume(const Frustum& _frustum)
{
    eiAssert( _frustum.r > _frustum.l, "Bad parametrization!" );
    eiAssert( _frustum.t > _frustum.b, "Bad parametrization!" );
    eiAssert( _frustum.f > _frustum.n, "Bad parametrization!" );
    // Use formular: V = (f-n)/3 * (a1 + sqrt(a1*a2) + a2) where a1 is the area on the
    // far plane and a2 that on the near plane.
    // Since a2 = sq(n/f) * a1 the formular simplifies to:
    // (f-n)/3 * (a1 + a1 * (n/f) + a1 * sq(n/f))
    float a1 = (_frustum.r - _frustum.l) * (_frustum.t - _frustum.b);
    float ratio = (_frustum.n / _frustum.f);
    return (_frustum.f - _frustum.n) / 3.0f * (a1 + (a1 + a1 * ratio) * ratio);
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

inline float surface(const Thetrahedron& _thetra)
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

inline float surface(const Plane& _plane)
{
    return INF;
}

inline float surface(const DOP& _dop)
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

inline float surface(const Ray& _ray)
{
    return 0.0f;
}

inline float surface(const Line& _line)
{
    return 0.0f;
}

inline float surface(const Capsule& _capsule)
{
    return 2 * PI * _capsule.radius * (2 * _capsule.radius + len(_capsule.b-_capsule.a));
}

inline float surface(const Frustum& _frustum)
{
    // Calculate area of all 6 sides and add them
    float hratio = (_frustum.f - _frustum.n) / _frustum.f;
    float ratio = _frustum.n / _frustum.f;
    // Far side
    float a1 = (_frustum.r - _frustum.l) * (_frustum.t - _frustum.b);
    // Near side
    float a2 = sq(ratio) * a1;
    // Top and bottom side
    float a3_a5 = (sqrt(sq(_frustum.f)+sq(_frustum.b)) + sqrt(sq(_frustum.f)+sq(_frustum.t))) * 0.5f * ((_frustum.r - _frustum.l) * (1.0f + ratio));
    // Left and right side
    float a4_a6 = (sqrt(sq(_frustum.f)+sq(_frustum.l)) + sqrt(sq(_frustum.f)+sq(_frustum.r))) * 0.5f * ((_frustum.t - _frustum.b) * (1.0f + ratio));
    return a1 + a2 + a3_a5 + a4_a6;
}
