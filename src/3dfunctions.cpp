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
