#include "ei/3dfunctions.hpp"

namespace ei {

// ************************************************************************* //
float volume(const Frustum& _frustum)
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
float surface(const Frustum& _frustum)
{
    // Calculate area of all 6 sides and add them
    float hratio = (_frustum.f - _frustum.n) / _frustum.f;
    float ratio = _frustum.n / _frustum.f;
    // Far side
    float a1 = (_frustum.r - _frustum.l) * (_frustum.t - _frustum.b);
    // Near side
    float a2 = sq(ratio) * a1;
    // Top and bottom side
    float a3_a5 = (sqrt(sq(_frustum.f)+sq(_frustum.b)) + sqrt(sq(_frustum.f)+sq(_frustum.t))) * hratio * 0.5f * ((_frustum.r - _frustum.l) * (1.0f + ratio));
    // Left and right side
    float a4_a6 = (sqrt(sq(_frustum.f)+sq(_frustum.l)) + sqrt(sq(_frustum.f)+sq(_frustum.r))) * hratio * 0.5f * ((_frustum.t - _frustum.b) * (1.0f + ratio));
    return a1 + a2 + a3_a5 + a4_a6;
}

// ************************************************************************* //
float distance(const Vec3& _point, const Segment& _line)
{
    Vec3 u = _line.b - _line.a;
    Vec3 w = _point - _line.a;

    float s = saturate(dot(u, w) / dot(u, u));
    // _line.a + u * s - _point
    return len(u * s - w);
}

// ************************************************************************* //
float distance(const Vec3& _point, const Triangle& _triangle)
{
    Vec3 a = _triangle.v1-_triangle.v0;
    Vec3 b = _triangle.v2-_triangle.v1;
    Vec3 c = _triangle.v2-_triangle.v0;
    Vec3 n = cross(c, a);

    // The point is inside if it is not on the right hand of any edge
    Vec3 p_v0 = _point - _triangle.v0;
    Vec3 p_v1 = _point - _triangle.v1;
    bool outSideA = dot(cross(p_v0, a), n) < 0.0f;
    bool outSideB = dot(cross(p_v1, b), n) < 0.0f;
    bool outSideC = dot(cross(c, p_v0), n) < 0.0f;   // Inverse sign because of use from p-v0

    if(!(outSideA || outSideB || outSideC))
        // Distance to the plane
        return dot(normalize(n), _point-_triangle.v0);

    // Minimum is somewhere on the three edges.
    float dist;
    dist =           lensq(p_v0) - sq(dot(p_v0, a)) / lensq(a);
    dist = min(dist, lensq(p_v1) - sq(dot(p_v1, b)) / lensq(b));
    dist = min(dist, lensq(p_v0) - sq(dot(p_v0, c)) / lensq(c));
    return sqrt(dist);
}

// ************************************************************************* //
float distance(const Segment& _line0, const Segment& _line1)
{
    // http://geomalgorithms.com/a07-_distance.html
    Vec3 u = _line0.b - _line0.a;
    Vec3 v = _line1.b - _line1.a;
    Vec3 w = _line1.a - _line0.a;
    float a = dot(u, u);
    float b = dot(u, v);
    float c = dot(v, v);
    float p = dot(u, w);
    float q = dot(v, w);
    float n = a * c - b * b;
    if(n == 0.0f)
    {
        // Parallel lines
        // If the two lines "overlap" the distance is equal over the whole
        // range. Otherwise 2 of the endpoints are closest. So in any case the
        // distance of one of the end points to the other lines must be the
        // minimum distance.
        float s0 = saturate(p / a);
        float d = len(s0 * u - w);
        Vec3 x = _line1.b - _line0.a;
        float s1 = saturate(dot(u, x) / a);
        d = min(d, len(s1 * u - x));
        float t0 = saturate(-q / c);
        d = min(d, len(t0 * v + w));
        Vec3 y = _line0.b - _line1.a;
        float t1 = saturate(dot(u, y) / c);
        d = min(d, len(t1 * v - y));
        return d;
    }
    float s = (b * q - c * p) / -n;
    float t = (a * q - b * p) / -n;
    // len(_line0.a + saturate(s) * u - (_line1.a + saturate(t) * v))
    return len(saturate(s) * u - w - saturate(t) * v);
}

}
