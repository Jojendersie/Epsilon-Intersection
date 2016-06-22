// ************************************************************************* //
inline Disc2D::Disc2D( Vec2 _center, float _radius ) :
    center(_center),
    radius(_radius)
{
}

// ************************************************************************* //
inline Disc2D::Disc2D( Vec2 _p0, Vec2 _p1 )
{
    center = (_p0 + _p1) * 0.5f;
    radius = len(_p0 - _p1) * 0.5f;
}

// ********************************************************************* //
inline Disc2D::Disc2D( const Triangle2D& _triangle ) :
    Disc2D( _triangle.v0, _triangle.v1, _triangle.v2 )
{
}

// ************************************************************************* //
inline Disc2D::Disc2D( const Rect2D& _rect )
{
    eiAssert( all(_rect.max >= _rect.min),
        "The input rect is degenerated! All components of max should be larger than in min." );
    // For a regular n-gone the center is simply the center of all vertices.
    // The radius is then equal to all vertices.
    center = (_rect.max + _rect.min) * 0.5f;
    radius = len(_rect.max - center);
}

// ************************************************************************* //
inline Disc2D::Disc2D( const ORect2D& _rect )
{
    eiAssert( all(_rect.size >= 0),
        "Side lengths of a rectangle should never be negative." );
    // The center is already given and the distance to all corners is equal.
    center = _rect.center;
    radius = len(_rect.size) * 0.5f;
}

// ************************************************************************* //
inline Disc2D::Disc2D( const Ellipse2D& _ellipse )
{
    eiAssert( all(_ellipse.radii >= 0),
        "Radii of an ellipse should never be negative." );
    center = _ellipse.center;
    radius = max( _ellipse.radii );
}

// ************************************************************************* //
inline Disc2D::Disc2D( const OEllipse2D& _ellipse )
{
    eiAssert( all(_ellipse.radii >= 0),
        "Radii of an ellipse should never be negative." );
    // Rotation does not change anything
    center = _ellipse.center;
    radius = max( _ellipse.radii );
}

// ************************************************************************* //
inline Disc2D::Disc2D( const Capsule2D& _capsule )
{
    center = (_capsule.seg.a + _capsule.seg.b) * 0.5f;
    radius = len(_capsule.seg.a - _capsule.seg.b) * 0.5f + _capsule.radius;
}

// ************************************************************************* //
inline Disc2D::Disc2D( const Segment2D& _line )
{
    center = (_line.a + _line.b) * 0.5f;
    radius = len(_line.a - _line.b) * 0.5f;
}

// ************************************************************************* //
inline bool Disc2D::operator== ( const Disc2D& _circle ) const
{
    return all( center == _circle.center ) && radius == _circle.radius;
}

// ************************************************************************* //
inline Rect2D::Rect2D( const Vec2& _min, const Vec2& _max ) :
    min(_min),
    max(_max)
{
    eiAssert( all(_min <= _max),
        "Minimum coordinates must be smaller or equal the maximum." );
}

// ************************************************************************* //
inline ORect2D::ORect2D( const Vec2& _center, const Vec2& _size, float _angle ) :
    center(_center),
    size(_size),
    angle(_angle)
{
    eiAssert( all(_size >= 0.0f), "Side lengths must be positive." );
}

// ************************************************************************* //
inline Triangle2D::Triangle2D( const Vec2& _v0, const Vec2& _v1, const Vec2& _v2 ) :
    v0(_v0),
    v1(_v1),
    v2(_v2)
{
}

// ************************************************************************* //
inline Ellipse2D::Ellipse2D( const Vec2& _center, const Vec2& _radii ) :
    center(_center),
    radii(_radii)
{
    eiAssert( all(_radii >= 0),
        "Radii of an ellipse should never be negative." );
}

// ************************************************************************* //
inline OEllipse2D::OEllipse2D( const Vec2& _center, const Vec2& _radii, float _angle ) :
    center(_center),
    radii(_radii),
    angle(_angle)
{
    eiAssert( all(_radii >= 0),
        "Radii of an ellipse should never be negative." );
}

// ************************************************************************* //
inline Segment2D::Segment2D( const Vec2& _a, const Vec2& _b ) :
    a(_a),
    b(_b)
{
}

// ************************************************************************* //
inline Ray2D::Ray2D( const Vec2& _origin, const Vec2& _direction ) :
    origin(_origin),
    direction(_direction)
{
    eiAssert( approx(len(_direction), 1.0f),
        "The direction vector must be normalized." );
}

// ************************************************************************* //
inline Capsule2D::Capsule2D( const Vec2& _a, const Vec2& _b, float _radius ) :
    seg(_a, _b),
    radius(_radius)
{
}

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
    return PI * _capsule.radius * _capsule.radius + 2.0f * _capsule.radius * len(_capsule.seg.b-_capsule.seg.a);
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
    return (_segment.a + _segment.b) * 0.5f;
}

inline Vec2 center(const Capsule2D& _capsule)
{
    return (_capsule.seg.a + _capsule.seg.b) * 0.5f;
}
