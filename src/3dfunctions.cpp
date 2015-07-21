#include "ei/3dtypes.hpp"

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
Vec3 center(const Frustum& _frustum)
{
    // For a full pyramid the centroid is located on the line segment from apex
    // to base centroid at 3/4 of that segment.
    Vec3 baseCentroidDir = (_frustum.l + _frustum.r) * 0.5f * cross(_frustum.up, _frustum.direction)
        + (_frustum.t + _frustum.b) * 0.5f * _frustum.up
        + _frustum.f * _frustum.direction;
    // Since this is a trunctad pyramid the height must be determined:
    // a1 = (r-l) * (t-b)    ground side area
    // a2 = a1 * sq(n/f)     intercept theorem
    // Relative distance of centroid to base centroid
    // d = (a1 + 2*sqrt(a1*a2) + 3*a2) / 4 * (a1 + sqrt(a1*a2) + a2)  http://mathworld.wolfram.com/PyramidalFrustum.html
    //   = (3*n*n + 2*f*n + f*f) / (n*n + f*n + f*f)
    // Distance from apex da = 1-d
    float da = 1 - (3*sq(_frustum.n) + 2*_frustum.n*_frustum.f + sq(_frustum.f))
                   / (4 * (sq(_frustum.n) + _frustum.n*_frustum.f + sq(_frustum.f)));
    // baseCentroidDir is the vector from apex to base and ´da´ the relative
    // distance from near to centroid plane. We need relative distance from
    // apex to centroid, which is ´(n+da*(f-n))/f´;
    return _frustum.apex + ((_frustum.n + da*(_frustum.f-_frustum.n))/_frustum.f) * baseCentroidDir;
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

// ************************************************************************* //
float distance(const Vec3& _point, const Box& _box)
{
    eiAssert(all(_box.min <= _box.max), "Invalid box!");
    // Connection vector to box corner
    Vec3 d = _box.min - _point;
    if(d.x < 0.0f) d.x = _point.x - _box.max.x;
    if(d.y < 0.0f) d.y = _point.y - _box.max.y;
    if(d.z < 0.0f) d.z = _point.z - _box.max.z;
    // Inner distance if all single distances are negative
    float innerDist = min(0.0f, max(d.x, max(d.y, d.z)));
    // Point is outside if len is greater than 0
    return len(max(d,Vec3(0.0f))) + innerDist;
}

// ************************************************************************* //
float distance(const Sphere& _sphere, const Box& _box)
{
    eiAssert(all(_box.min <= _box.max), "Invalid box!");
    eiAssert(_sphere.radius >= 0.0f, "Invalid sphere!");
    Vec3 d = max(_box.min - _sphere.center, _sphere.center - _box.max);
    return max(len(max(d, Vec3(0.0f))) - _sphere.radius, 0.0f);
}




// ************************************************************************* //
// OBJECT TRANSFORMATIONS
// ************************************************************************* //

OBox transform(const Box& _box, const Quaternion& _rotation)
{
    OBox box(_box);
    box.orientation *= _rotation;
    box.center = transform(box.center, _rotation);
    return box;
}

OBox transform(const OBox& _box, const Quaternion& _rotation)
{
    return OBox(transform(_box.center, _rotation), _box.sides, _box.orientation * _rotation);
}

Box transform(const Box& _box, const Vec3& _translation)
{
    return Box(_box.min + _translation, _box.max + _translation);
}

OBox transform(const OBox& _box, const Vec3& _translation)
{
    return OBox(_box.center + _translation, _box.sides, _box.orientation);
}

OBox transform(const Box& _box, const Quaternion& _rotation, const Vec3& _translation)
{
    OBox box(_box);
    box.orientation *= _rotation;
    box.center = transform(box.center, _rotation) + _translation;
    return box;
}

OBox transform(const OBox& _box, const Quaternion& _rotation, const Vec3& _translation)
{
    return OBox(transform(_box.center, _rotation) + _translation,
        _box.sides, _box.orientation * _rotation);
}

}
