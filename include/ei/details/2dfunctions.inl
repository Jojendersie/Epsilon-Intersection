// ************************************************************************* //
// AREA METHODS
// ************************************************************************* //
inline float area(const Disc2D& _disc)
{
    return _disc.radius * _disc.radius * PI;
}

inline float area(const Rect2D& _rect)
{
    eiAssert(all(_rect.max >= _rect.min), "Rect max-boundary must be larger than its min-boundary!");
    return (_rect.max.x - _rect.min.x) * (_rect.max.y - _rect.min.y);
}

inline float area(const ORect2D& _orect)
{
    eiAssert(all(_orect.size >= 0.0f), "Rect must have positive side lengths!");
    return _orect.size.x * _orect.size.y;
}

inline float area(const Triangle2D& _triangle)
{
    // Use determinant rule for the cross product of two sides
    return 0.5f * ((_triangle.v1.x-_triangle.v0.x) * (_triangle.v2.y-_triangle.v0.y)
                 - (_triangle.v2.x-_triangle.v0.x) * (_triangle.v1.y-_triangle.v0.y));
}

inline float area(const Ellipse2D& _ellipse)
{
    return PI * _ellipse.radii.x * _ellipse.radii.y;
}

inline float area(const OEllipse2D& _oellipse)
{
    return PI * _oellipse.radii.x * _oellipse.radii.y;
}

inline float area(const Segment2D& _segment)
{
    return 0.0f;
}

inline float area(const Ray2D& _ray)
{
    return 0.0f;
}

inline float area(const Capsule2D& _capsule)
{
    return PI * _capsule.radius * _capsule.radius + 2.0f * _capsule.radius * len(_capsule.p1-_capsule.p0);
}

// ************************************************************************* //
// CENTER METHODS
// ************************************************************************* //
inline Vec2 center(const Disc2D& _disc)
{
    return _disc.center;
}

inline Vec2 center(const Rect2D& _rect)
{
    return 0.5f * (_rect.min + _rect.max);
}

inline Vec2 center(const ORect2D& _orect)
{
    return _orect.center;
}

inline Vec2 center(const Triangle2D& _triangle)
{
    return (_triangle.v0 + _triangle.v1 + _triangle.v2) / 3.0f;
}

inline Vec2 center(const Ellipse2D& _ellipse)
{
    return _ellipse.center;
}

inline Vec2 center(const OEllipse2D& _oellipse)
{
    return _oellipse.center;
}

inline Vec2 center(const Segment2D& _segment)
{
    return (_segment.p0 + _segment.p1) * 0.5f;
}

inline Vec2 center(const Capsule2D& _capsule)
{
    return (_capsule.p0 + _capsule.p1) * 0.5f;
}
