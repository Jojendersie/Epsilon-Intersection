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
// INTERSECTION METHODS
// ************************************************************************* //
inline bool intersects( const Vec3& _point, const Box& _box )
{
    if(_point.x < _box.min.x) return false;
    if(_point.y < _box.min.y) return false;
    if(_point.z < _box.min.z) return false;
    if(_point.x > _box.max.x) return false;
    if(_point.y > _box.max.y) return false;
    return _point.z <= _box.max.z;
    //return all(_point >= _box.min) && all(_point <= _box.max);
}

// ************************************************************************* //
inline bool intersects( const Vec3& _point, const DOP& _dop )
{
    float p = -dot(_point, _dop.n);
    if( p > _dop.d0 ) return false;
    return p > _dop.d1;
}