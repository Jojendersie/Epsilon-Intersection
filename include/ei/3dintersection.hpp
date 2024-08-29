#pragma once

#include "3dtypes.hpp"

namespace ei {

    /// \brief Descriptor for hit results with boxes.
    /// \details Sides are relative to the local system of the box.
    ///     For OBoxs this is not the absolute coordinate. E.g. if a box is
    ///     rotated by π/2 two coordinates seem to be exchanged.
    enum struct HitSide
    {
        X_NEG = 0x01,         ///< Left.
        X_POS = 0x02,         ///< Right.
        X     = 0x03,         ///< Any side in x-direction
        Y_NEG = 0x04,         ///< Bottom.
        Y_POS = 0x08,         ///< Top.
        Y     = 0x0c,         ///< Any side in y-direction
        Z_NEG = 0x10,         ///< Front.
        Z_POS = 0x20,         ///< Back.
        Z     = 0x30,         ///< Any side in z-direction
    };

    // ********************************************************************* //
    // DISTANCE FUNCTIONS                                                    //
    // ********************************************************************* //

    /// \brief Get the euclidean distance between two objects.
    /// \details The distance for point-solid queries can be negative. All
    ///     other geometries return 0 if they intersect.
    EIAPI float distance(const Vec3& _point0, const Vec3& _point1)            // TESTED
    {
        return len(_point1 - _point0);
    }

    EIAPI float distance(const Vec3& _point, const Segment& _line)            // TESTED
    {
        Vec3 u = _line.b - _line.a;
        Vec3 w = _point - _line.a;

        float s = saturate(dot(u, w) / dot(u, u));
        // _line.a + u * s - _point
        return len(u * s - w);
    }

    EIAPI float distance(const Vec3& _point, const Triangle& _triangle)       // TESTED
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

    EIAPI float distance(const Vec3& _point, const Sphere& _sphere)           // TESTED
    {
        return distance(_sphere.center, _point) - _sphere.radius;
    }

    EIAPI float distance(const Vec3& _point, const Capsule& _capsule)         // TESTED
    {
        return distance(_point, _capsule.seg) - _capsule.radius;
    }

    EIAPI float distance(const Vec3& _point, const Ray& _ray)
    {
        // Project point to the ray
        Vec3 o = _point - _ray.origin;
        float odotd = max(0.0f, dot(o, _ray.direction));
        // Compute difference between projection and original
        o -= _ray.direction * odotd;
        return len(o);
    }

    EIAPI float distance(const Vec3& _point, const Box& _box)                 // TESTED
    {
        eiAssert(_box.min <= _box.max, "Invalid box!");
        // Connection vector to box corner
        Vec3 d = _box.min - _point;
        if(d.x < 0.0f) d.x = max(d.x, _point.x - _box.max.x);
        if(d.y < 0.0f) d.y = max(d.y, _point.y - _box.max.y);
        if(d.z < 0.0f) d.z = max(d.z, _point.z - _box.max.z);
        // Inner distance if all single distances are negative
        float innerDist = min(0.0f, max(d.x, d.y, d.z));
        // Point is outside if len is greater than 0
        return len(max(d,Vec3(0.0f))) + innerDist;
    }

    EIAPI float distance(const Vec3& _point, const Plane& _plane)             // TESTED
    {
        eiAssert(approx(1.0f, len(_plane.n)), "The plane is not normalized!");
        return dot(_plane.n, _point) + _plane.d;
    }


    /// \returns -x if dot(n,_point) larger then both, x if smaller and 0 if it
    ///     is between both planes
    EIAPI float distance(const Vec3& _point, const DOP& _dop)
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

    EIAPI float distance(const Sphere& _sphere, const Segment& _segment)      // TESTED
    {
        return max(0.0f, distance(_sphere.center, _segment) - _sphere.radius);
    }

    EIAPI float distance(const Sphere& _sphere, const Capsule& _capsule)      // TESTED
    {
        return max(0.0f, distance(_sphere.center, _capsule.seg) - _capsule.radius - _sphere.radius);
    }

    EIAPI float distance(const Sphere& _sphere, const Box& _box)              // TESTED
    {
        eiAssert(_box.min <= _box.max, "Invalid box!");
        eiAssert(_sphere.radius >= 0.0f, "Invalid sphere!");
        Vec3 d = max(_box.min - _sphere.center, _sphere.center - _box.max);
        return max(len(max(d, Vec3(0.0f))) - _sphere.radius, 0.0f);
    }

    EIAPI float distance(const Sphere& _sphere, const Plane& _plane)          // TESTED
    {
        eiAssert(approx(1.0f, len(_plane.n)), "The plane is not normalized!");
        eiAssert(_sphere.radius >= 0.0f, "Invalid sphere!");
        float d = dot(_plane.n, _sphere.center) + _plane.d;
        // Keep sign for the plane distances
        if(d < 0.0f) return min(d + _sphere.radius, 0.0f);
        else return max(d - _sphere.radius, 0.0f);
    }
    
    EIAPI float distance(const Segment& _line0, const Segment& _line1)        // TESTED
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

    EIAPI float distance(const Capsule& _capsule0, const Capsule& _capsule1)  // TESTED
    {
        return max(0.0f, distance(_capsule0.seg, _capsule1.seg) - _capsule0.radius - _capsule1.radius);
    }

    EIAPI float distance(const Segment& _line, const Vec3& _point)            { return distance(_point, _line); }
    EIAPI float distance(const Triangle& _triangle, const Vec3& _point)       { return distance(_point, _triangle); }
    EIAPI float distance(const Sphere& _sphere, const Vec3& _point)           { return distance(_point, _sphere); }
    EIAPI float distance(const Capsule& _capsule, const Vec3& _point)         { return distance(_point, _capsule); }
    EIAPI float distance(const Capsule& _capsule, const Sphere& _sphere)      { return distance(_sphere, _capsule); }
    EIAPI float distance(const Box& _box, const Vec3& _point)                 { return distance(_point, _box); }
    EIAPI float distance(const Box& _box, const Sphere& _sphere)              { return distance(_sphere, _box); }
    EIAPI float distance(const Plane& _plane, const Vec3& _point)             { return distance(_point, _plane); }
    EIAPI float distance(const Plane& _plane, const Sphere& _sphere)          { return distance(_sphere, _plane); }
    EIAPI float distance(const DOP& _dop, const Vec3& _point)                 { return distance(_point, _dop); }
    EIAPI float distance(const Ray& _ray, const Vec3& _point)                 { return distance(_point, _ray); }


    // ********************************************************************* //
    // INTERSECTION TESTS                                                    //
    // ********************************************************************* //

    /// \brief Do two spheres intersect, touch or is one inside the other?
    /// \return true if both spheres have at least one point in common.
    EIAPI bool intersects( const Sphere& _sphere0, const Sphere& _sphere1 )   // TESTED
    {
        return lensq(_sphere1.center-_sphere0.center) <= sq(_sphere0.radius + _sphere1.radius);
    }

    /// \brief Does a point lie inside a sphere/on the boundary?
    /// \return true if the point is inside or on the boundary.
    EIAPI bool intersects( const Vec3& _point, const Sphere& _sphere )        // TESTED
    {
        return lensq(_point-_sphere.center) <= sq(_sphere.radius);
    }

    EIAPI bool intersects( const Sphere& _sphere, const Vec3& _point )  { return intersects( _point, _sphere ); }

    /// \brief Do a sphere and a box intersect?
    /// \return true if the sphere and the box have at least one point in common.
    EIAPI bool intersects( const Sphere& _sphere, const Box& _box )               // TESTED
    {
        float distSq = sq(_sphere.radius);
        if (_sphere.center.x < _box.min.x) distSq -= sq(_sphere.center.x - _box.min.x);
        else if (_sphere.center.x > _box.max.x) distSq -= sq(_sphere.center.x - _box.max.x);
        if(distSq < 0.0f) return false;
        if (_sphere.center.y < _box.min.y) distSq -= sq(_sphere.center.y - _box.min.y);
        else if (_sphere.center.y > _box.max.y) distSq -= sq(_sphere.center.y - _box.max.y);
        if(distSq < 0.0f) return false;
        if (_sphere.center.z < _box.min.z) distSq -= sq(_sphere.center.z - _box.min.z);
        else if (_sphere.center.z > _box.max.z) distSq -= sq(_sphere.center.z - _box.max.z);
        // Vectorized alternative (much slower)
        //distSq -= lensq(max(Vec3(0.0f), _box.min - _sphere.center));
        //distSq -= lensq(max(Vec3(0.0f), _sphere.center - _box.max));
        return distSq > 0.0f;
    }

    EIAPI bool intersects( const Box& _box, const Sphere& _sphere )  { return intersects( _sphere, _box ); }

    /// \brief Does a point lie inside a box/on the boundary?
    /// \return true if the point is inside or on the boundary.
    EIAPI bool intersects( const Vec3& _point, const Box& _box )              // TESTED
    {
        if(_point.x < _box.min.x) return false;
        if(_point.y < _box.min.y) return false;
        if(_point.z < _box.min.z) return false;
        if(_point.x > _box.max.x) return false;
        if(_point.y > _box.max.y) return false;
        return _point.z <= _box.max.z;
        //return all(_point >= _box.min) && all(_point <= _box.max);
    }

    EIAPI bool intersects( const Box& _box, const Vec3& _point )  { return intersects( _point, _box ); }

    /// \brief Does a point lie inside a box/on the boundary?
    /// \return true if both boxes have at least one point in common.
    EIAPI bool intersects( const Box& _box0, const Box& _box1 )               // TESTED
    {
        eiAssert(_box0.min <= _box0.max, "Box0 is invalid.");
        eiAssert(_box1.min <= _box1.max, "Box1 is invalid.");
        // There must be an intersection if the sum of side length is larger
        // than that of the bounding box.
//        return all((max(_box0.max, _box1.max) - min(_box0.min, _box1.min))
//            <= ((_box0.max - _box0.min) + (_box1.max - _box1.min)));
        // Non-vector variant is faster
        return (max(_box0.max.x, _box1.max.x) - min(_box0.min.x, _box1.min.x)) <= ((_box0.max.x - _box0.min.x) + (_box1.max.x - _box1.min.x))
            && (max(_box0.max.y, _box1.max.y) - min(_box0.min.y, _box1.min.y)) <= ((_box0.max.y - _box0.min.y) + (_box1.max.y - _box1.min.y))
            && (max(_box0.max.z, _box1.max.z) - min(_box0.min.z, _box1.min.z)) <= ((_box0.max.z - _box0.min.z) + (_box1.max.z - _box1.min.z));
    }

    /// \brief Do two capsules intersect, touch or is one inside the other?
    /// \return true if both spheres have at least one point in common.
    EIAPI bool intersects( const Capsule& _capsule0, const Capsule& _capsule1 ) // TESTED
    {
        return distance(_capsule0.seg, _capsule1.seg) <= (_capsule0.radius + _capsule1.radius);
    }

    /// \brief Does a point lie inside an ellipsoid/on the boundary?
    /// \return true if the point is inside or on the boundary.
    EIAPI bool intersects( const Vec3& _point, const Ellipsoid& _ellipsoid )  // TESTED
    {
        // Use ellipsoid equation.
        //return lensq(Vec<double,3>((_point - _ellipsoid.center) / _ellipsoid.radii)) <= 1.0;
        //return lensq((_point - _ellipsoid.center) / _ellipsoid.radii) <= 1.0;
        Vec3 o = (_point - _ellipsoid.center);
        if( _ellipsoid.radii.x != 0.0f ) o.x /= _ellipsoid.radii.x; else o.x *= 3.402823466e+38f;
        if( _ellipsoid.radii.y != 0.0f ) o.y /= _ellipsoid.radii.y; else o.y *= 3.402823466e+38f;
        if( _ellipsoid.radii.z != 0.0f ) o.z /= _ellipsoid.radii.z; else o.z *= 3.402823466e+38f;
        return double(o.x*o.x) + double(o.y*o.y) + double(o.z*o.z) <= 1.0;
        // The following is a try to get a more stable solution without
        // addition, but it is less stable in the end.
        //Vec3 t = (_point - _ellipsoid.center) / _ellipsoid.radii;
        //return exp(t.x*t.x)*exp(t.y*t.y)*exp(t.z*t.z) <= 2.718281828f;
    }

    EIAPI bool intersects( const Ellipsoid& _ellipsoid, const Vec3& _point )  { return intersects( _point, _ellipsoid ); }

    /// \brief Does a point lie inside an oriented ellipsoid/on the boundary?
    /// \return true if the point is inside or on the boundary.
    EIAPI bool intersects( const Vec3& _point, const OEllipsoid& _oellipsoid ) // TESTED
    {
        Vec3 o = transform(_point - _oellipsoid.center, _oellipsoid.orientation);
        // Use ellipsoid equation.
        if( _oellipsoid.radii.x != 0.0f ) o.x /= _oellipsoid.radii.x; else o.x *= 3.402823466e+38f;
        if( _oellipsoid.radii.y != 0.0f ) o.y /= _oellipsoid.radii.y; else o.y *= 3.402823466e+38f;
        if( _oellipsoid.radii.z != 0.0f ) o.z /= _oellipsoid.radii.z; else o.z *= 3.402823466e+38f;
        return double(o.x*o.x) + double(o.y*o.y) + double(o.z*o.z) <= 1.0;
    }

    EIAPI bool intersects( const OEllipsoid& _oellipsoid, const Vec3& _point )  { return intersects( _point, _oellipsoid ); }

    /// \brief Does a point lie inside a DOP or on the boundary?
    /// \return true if the point is inside or on the boundary.
    EIAPI bool intersects( const Vec3& _point, const DOP& _dop )              // TESTED
    {
        float p = -dot(_point, _dop.n);
        if( p > _dop.d0 ) return false;
        return p > _dop.d1;
    }

    EIAPI bool intersects( const DOP& _dop, const Vec3& _point )  { return intersects( _point, _dop ); }

    /// \brief Intersection test between point and oriented box.
    /// \return true if the point and the oriented box have at least on point in common.
    EIAPI bool intersects( const Vec3& _point, const OBox& _obox )            // TESTED
    {
        Vec3 alignedPoint = _point - _obox.center;
        alignedPoint = transform( alignedPoint, conjugate(_obox.orientation) );
        if(alignedPoint.x < - _obox.halfSides.x) return false;
        if(alignedPoint.y < - _obox.halfSides.y) return false;
        if(alignedPoint.z < - _obox.halfSides.z) return false;
        if(alignedPoint.x >   _obox.halfSides.x) return false;
        if(alignedPoint.y >   _obox.halfSides.y) return false;
        return alignedPoint.z <= _obox.halfSides.z;
    }

    EIAPI bool intersects( const OBox& _obox, const Vec3& _point )          { return intersects(_point, _obox); }


    /// \brief Get the intersection distance of a ray and a plane.
    /// \param [out,opt] _distance The distance to the intersection
    /// \return false if the ray is (near) parallel to the plane.
    ///     The plane counts as parallel when the intersection is at +-inf.
    ///     A ray inside the plane counts as intersecting at distance=0.
    EIAPI bool intersects( const Ray& _ray, const Plane& _plane, float& _distance ) // TESTED
    {
        const float dist_perpendicular = dot(_ray.origin, _plane.n) + _plane.d;
        // Origin in plane? Need to capture the special case to avoid NaN
        // if the ray is parallel as well.
        if (dist_perpendicular == 0.0f)
        {
            _distance = 0.0f;
            return true;
        }
        _distance = -dist_perpendicular / dot(_ray.direction, _plane.n);
        return isfinite(_distance);
    }

    EIAPI bool intersects( const Plane& _plane, const Ray& _ray, float& _distance ) { return intersects(_ray, _plane, _distance); }


    /// \brief Do a ray and a sphere intersect or touch?
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    ///
    ///     If the ray starts inside the distance to the exit point is returned.
    /// \param[out,opt] _distance2 The distance to the other intersection point.
    ///     This is smaller zero and _distance > _distance2, if the ray started inside.
    ///     Otherwise, it is the distance to the exit point with _distance2 >= _distance
    /// \return true if the ray has at least one point in common with the sphere
    EIAPI bool intersects( const Ray& _ray, const Sphere& _sphere )
    {
        // Go towards closest point and compare its distance to the radius
        Vec3 o = _sphere.center - _ray.origin;
        float odotd = max(0.0f, dot(o, _ray.direction)); // The max0 ensures we do not travel along negative ray direction (there might be an intersection behind the ray)
        o -= _ray.direction * odotd;
        return lensq(o) <= _sphere.radius * _sphere.radius;
    }

    EIAPI bool intersects( const Ray& _ray, const Sphere& _sphere, float& _distance, float& _distance2)
    {
        float rSq = _sphere.radius * _sphere.radius;
        // Go towards closest point and compare its distance to the radius
        Vec3 o = _sphere.center - _ray.origin;
        float odoto = dot(o, o);				// Distance ray.origin <-> sphere.center
        float odotd = dot(o, _ray.direction);	// Distance ray.origin <-> closest point to sphere center
        // If the ray starts outside and the closest point is behind us there is no intersection
        bool outside = odoto > rSq;
        if(odotd < 0.0f && outside) return false;
        // Jump to the closest point for numerical reasons
        o -= _ray.direction * odotd;
        float dClosestSq = dot(o, o);
        // If the closest point is outside the sphere there is no intersection
        if(dClosestSq > rSq) return false;
        // Compute the t for the closest intersection
        float innerSegLenHalf = sqrt(rSq - dClosestSq);
        // If the ray started inside, the distance for the first point must
        // be added instead substracted.
        if(outside) {
            _distance  = odotd - innerSegLenHalf;
            _distance2 = odotd + innerSegLenHalf;
        } else {
            _distance  = odotd + innerSegLenHalf;
            _distance2 = odotd - innerSegLenHalf;
        }
        return true;
    }

    EIAPI bool intersects( const Sphere& _sphere, const Ray& _ray ) { return intersects(_ray, _sphere); }
    EIAPI bool intersects( const Sphere& _sphere, const Ray& _ray, float& _distance, float& _distance2 ) { return intersects(_ray, _sphere, _distance, _distance2); }

    /// \brief Do a ray and an ellipsoid intersect or touch?
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    ///
    ///     If the ray starts inside 0 is returned.
    /// \return true if there is at least one common point between ray and
    ///     ellipsoid.
    EIAPI bool intersects( const Ray& _ray, const Ellipsoid& _ellipsoid )     // TESTED
    {
        // Translate to origin
        Vec3 o = _ray.origin - _ellipsoid.center;

        // Go to sphere for numerical stability
        float r = max(_ellipsoid.radii);
        float t = dot(o, _ray.direction);
        t = - t - r;
        if(t > 0.0f)
            o += _ray.direction * t;

        // Scale system
        o /= _ellipsoid.radii;
        float odoto = dot(o, o);

        // Test if quadratic equation for the hit point has a solution
        Vec3 d = _ray.direction / _ellipsoid.radii;
        float odotd = dot(o, d);
        float ddotd = dot(d, d);
        float phalf = odotd;
        float q = (odoto - 1.0f);
        return phalf * (phalf / ddotd) - q >= 0.0f;
    }

    EIAPI bool intersects( const Ray& _ray, const Ellipsoid& _ellipsoid, float& _distance )
    {
        // Translate to origin
        Vec3 o = _ray.origin - _ellipsoid.center;

        // Go to sphere for numerical stability
        float r = max(_ellipsoid.radii);
        float t = dot(o, _ray.direction);
        t = max(0.0f, - t - r);
        if(t > 0.0f)
            o += _ray.direction * t;

        // Scale system
        o /= _ellipsoid.radii;
        float odoto = dot(o, o);

        // Test if quadratic equation for the hit point has a solution
        Vec3 d = _ray.direction / _ellipsoid.radii;//max(_ellipsoid.radii, Vec3(1.0f));
        float odotd = dot(o, d);
        float ddotd = dot(d, d);
        float phalf = odotd / ddotd;
        float q = (odoto - 1.0f) / ddotd;
        float rad = phalf * phalf - q;
        if( rad < 0.0f ) return false;
        rad = sqrt(rad);
        phalf -= t;
        _distance = max(-phalf - rad, 0.0f);//(-phalf - rad < 0.0f) ? (-phalf + rad) : (-phalf - rad);
        return -phalf + rad >= 0.0f;
    }

    EIAPI bool intersects( const Ellipsoid& _ellipsoid, const Ray& _ray )  { return intersects( _ray, _ellipsoid ); }
    EIAPI bool intersects( const Ellipsoid& _ellipsoid, const Ray& _ray, float& _distance )  { return intersects( _ray, _ellipsoid, _distance ); }

    /// \brief Do a ray and a box intersect or touch?
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    ///
    ///     If the ray starts inside 0 is returned.
    /// \param [out,opt] _distanceExit The ray parameter (distance) for the
    ///     second intersection point in positive direction.
    ///
    ///     This is always the exit point. If the ray starts on the boundary
    ///     and shows away _distance and _distanceExit are the same (0).
    /// \param [out,opt] _side Returns which side of the box was hit.
    /// \return true if there is at least one common point between ray and box
    EIAPI bool intersects( const Ray& _ray, const Box& _box )                 // TESTED
    {
        float t0 = (_box.min.x - _ray.origin.x) / _ray.direction.x;
        float t1 = (_box.max.x - _ray.origin.x) / _ray.direction.x;
        float tmin = min(t0, t1);
        float tmax = max(t0, t1);
        if(tmin != tmin) tmin = -INF;
        if(tmax != tmax) tmax = INF;
        if(tmax < 0.0f) return false;
        t0 = (_box.min.y - _ray.origin.y) / _ray.direction.y;
        t1 = (_box.max.y - _ray.origin.y) / _ray.direction.y;
        float min2 = min(t0, t1);
        float max2 = max(t0, t1);
        if(min2 != min2) min2 = -INF;
        if(max2 != max2) max2 = INF;
        tmin = max(tmin, min2);
        tmax = min(tmax, max2);
        if(tmax < 0.0f || tmin > tmax) return false;
        t0 = (_box.min.z - _ray.origin.z) / _ray.direction.z;
        t1 = (_box.max.z - _ray.origin.z) / _ray.direction.z;
        min2 = min(t0, t1);
        max2 = max(t0, t1);
        if(min2 != min2) min2 = -INF;
        if(max2 != max2) max2 = INF;
        tmin = max(tmin, min2);
        tmax = min(tmax, max2);
        return (tmax >= 0.0f) && (tmin <= tmax);
        /*Vec3 tbot = (_box.min - _ray.origin) / _ray.direction;
        Vec3 ttop = (_box.max - _ray.origin) / _ray.direction;
        Vec3 tmin = min(ttop, tbot);
        Vec3 tmax = max(ttop, tbot);
        float firstHit = max(0.0f, max(tmin));
        float lastHit = min(tmax);
        return firstHit <= lastHit;*/
    }

    EIAPI bool intersects( const Ray& _ray, const Box& _box, float& _distance )
    {
        float tmin, tmax, tymin, tymax, tzmin, tzmax;
        const Vec3 * bbounds = (const Vec3 *)&_box;

        tmin  = (bbounds[1-heaviside(_ray.direction.x)].x - _ray.origin.x) / _ray.direction.x;  if(isnan(tmin)) tmin = 0.0f;
        tmax  = (bbounds[  heaviside(_ray.direction.x)].x - _ray.origin.x) / _ray.direction.x;  if(isnan(tmax)) tmax = 0.0f;
        tymin = (bbounds[1-heaviside(_ray.direction.y)].y - _ray.origin.y) / _ray.direction.y;  if(isnan(tymin)) tymin = 0.0f;
        tymax = (bbounds[  heaviside(_ray.direction.y)].y - _ray.origin.y) / _ray.direction.y;  if(isnan(tymax)) tymax = 0.0f;

        if((tmin > tymax) || (tymin > tmax))
            return false;
        if(tymin > tmin) tmin = tymin;
        if(tymax < tmax) tmax = tymax;
        if(tmax < 0.0f)
            return false;

        tzmin = (bbounds[1-heaviside(_ray.direction.z)].z - _ray.origin.z) / _ray.direction.z;  if(isnan(tzmin)) tzmin = 0.0f;
        tzmax = (bbounds[  heaviside(_ray.direction.z)].z - _ray.origin.z) / _ray.direction.z;  if(isnan(tzmax)) tzmax = 0.0f;

        if ((tmin > tzmax) || (tzmin > tmax) || (tzmax < 0.0f))
            return false;

        _distance = max(tzmin, tmin, 0.0f);
        return true;
    }

    EIAPI bool intersects( const Ray& _ray, const Box& _box, float& _distance, HitSide& _side )
    {
        // Box intersection algorithm:
        // http://people.csail.mit.edu/amy/papers/box-jgt.pdf
        // Projection to all planes
        float min0, tmin, max0, tmax;
        HitSide tside;
        if( _ray.direction.x >= 0 ) {
            min0 = (_box.min.x - _ray.origin.x) / _ray.direction.x;
            max0 = (_box.max.x - _ray.origin.x) / _ray.direction.x;
            _side = min0 < 0.0f ? HitSide::X_POS : HitSide::X_NEG;
        } else {
            max0 = (_box.min.x - _ray.origin.x) / _ray.direction.x;
            min0 = (_box.max.x - _ray.origin.x) / _ray.direction.x;
            _side = min0 < 0.0f ? HitSide::X_NEG : HitSide::X_POS;
        }

        if( _ray.direction.y >= 0 ) {
            tmin = (_box.min.y - _ray.origin.y) / _ray.direction.y;
            tmax = (_box.max.y - _ray.origin.y) / _ray.direction.y;
            tside = tmin < 0.0f ? HitSide::Y_POS : HitSide::Y_NEG;
        } else {
            tmax = (_box.min.y - _ray.origin.y) / _ray.direction.y;
            tmin = (_box.max.y - _ray.origin.y) / _ray.direction.y;
            tside = tmin < 0.0f ? HitSide::Y_NEG : HitSide::Y_POS;
        }

        if( (min0 > tmax) || (tmin > max0) )
            return false;
        if( tmin > min0 ) { min0 = tmin; _side = tside; }
        if( tmax < max0 ) max0 = tmax;

        if( _ray.direction.z >= 0 ) {
            tmin = (_box.min.z - _ray.origin.z) / _ray.direction.z;
            tmax = (_box.max.z - _ray.origin.z) / _ray.direction.z;
            tside = tmin < 0.0f ? HitSide::Z_POS : HitSide::Z_NEG;
        } else {
            tmax = (_box.min.z - _ray.origin.z) / _ray.direction.z;
            tmin = (_box.max.z - _ray.origin.z) / _ray.direction.z;
            tside = tmin < 0.0f ? HitSide::Z_NEG : HitSide::Z_POS;
        }

        if( (min0 > tmax) || (tmin > max0) )
            return false;
        if( 0.0f > max0 || 0.0f > tmax ) return false;
        if( tmin > min0 ) _side = tside;

        _distance = max(min0, tmin, 0.0f);
        return true;
    }

    EIAPI bool intersects( const Ray& _ray, const Box& _box, float& _distance, float& _distanceExit )
    {
        float tmin, tmax, tymin, tymax, tzmin, tzmax;
        const Vec3 * bbounds = (const Vec3 *)&_box;

        tmax  = (bbounds[_ray.direction.x < 0.0f ? 0 : 1].x - _ray.origin.x) / _ray.direction.x;
        tmin  = (bbounds[_ray.direction.x < 0.0f ? 1 : 0].x - _ray.origin.x) / _ray.direction.x;
        tymin = (bbounds[_ray.direction.y < 0.0f ? 1 : 0].y - _ray.origin.y) / _ray.direction.y;
        tymax = (bbounds[_ray.direction.y < 0.0f ? 0 : 1].y - _ray.origin.y) / _ray.direction.y;

        if((tmin > tymax) || (tymin > tmax))
            return false;
        if(tymin > tmin) tmin = tymin;
        if(tymax < tmax) tmax = tymax;
        if(tmax < 0.0f)
            return false;

        tzmin = (bbounds[_ray.direction.z < 0.0f ? 1 : 0].z - _ray.origin.z) / _ray.direction.z;
        tzmax = (bbounds[_ray.direction.z < 0.0f ? 0 : 1].z - _ray.origin.z) / _ray.direction.z;

        if ((tmin > tzmax) || (tzmin > tmax) || (tzmax < 0.0f))
            return false;

        _distance = max(tzmin, tmin, 0.0f);
        _distanceExit = min(tzmax, tmax);
        return true;
    }

    EIAPI bool intersects( const FastRay& _ray, const Box& _box, float& _distance )
    {
        float tmin, tmax, tymin, tymax, tzmin, tzmax;
        const Vec3 * bbounds = (const Vec3 *)&_box;

        /*tmax  = (bbounds[_ray.invDirection.x < 0.0f ? 0 : 1].x - _ray.origin.x) * _ray.invDirection.x;
        tmin  = (bbounds[_ray.invDirection.x < 0.0f ? 1 : 0].x - _ray.origin.x) * _ray.invDirection.x;
        tymin = (bbounds[_ray.invDirection.y < 0.0f ? 1 : 0].y - _ray.origin.y) * _ray.invDirection.y;
        tymax = (bbounds[_ray.invDirection.y < 0.0f ? 0 : 1].y - _ray.origin.y) * _ray.invDirection.y;*/
        tmax  = bbounds[_ray.invDirection.x < 0.0f ? 0 : 1].x * _ray.invDirection.x - _ray.oDivDir.x;
        tmin  = bbounds[_ray.invDirection.x < 0.0f ? 1 : 0].x * _ray.invDirection.x - _ray.oDivDir.x;
        tymin = bbounds[_ray.invDirection.y < 0.0f ? 1 : 0].y * _ray.invDirection.y - _ray.oDivDir.y;
        tymax = bbounds[_ray.invDirection.y < 0.0f ? 0 : 1].y * _ray.invDirection.y - _ray.oDivDir.y;

        if((tmin > tymax) || (tymin > tmax))
            return false;
        if(tymin > tmin) tmin = tymin;
        if(tymax < tmax) tmax = tymax;
        if(tmax < 0.0f)
            return false;

        //tzmin = (bbounds[_ray.invDirection.z < 0.0f ? 1 : 0].z - _ray.origin.z) * _ray.invDirection.z;
        //tzmax = (bbounds[_ray.invDirection.z < 0.0f ? 0 : 1].z - _ray.origin.z) * _ray.invDirection.z;
        tzmin = bbounds[_ray.invDirection.z < 0.0f ? 1 : 0].z * _ray.invDirection.z - _ray.oDivDir.z;
        tzmax = bbounds[_ray.invDirection.z < 0.0f ? 0 : 1].z * _ray.invDirection.z - _ray.oDivDir.z;

        if((tmin > tzmax) || (tzmin > tmax) || (tzmax < 0.0f))
            return false;

        _distance = max(tzmin, tmin, 0.0f);
        return true;
    }

    EIAPI bool intersects( const Box& _box, const Ray& _ray )  { return intersects( _ray, _box ); }
    EIAPI bool intersects( const Box& _box, const Ray& _ray, float& _distance )  { return intersects( _ray, _box, _distance ); }
    EIAPI bool intersects( const Box& _box, const Ray& _ray, float& _distance, HitSide& _side )  { return intersects( _ray, _box, _distance, _side ); }
    EIAPI bool intersects( const Box& _box, const Ray& _ray, float& _distance, float& _distanceExit )  { return intersects( _ray, _box, _distance, _distanceExit ); }
    EIAPI bool intersects( const Box& _box, const FastRay& _ray, float& _distance )  { return intersects( _ray, _box, _distance ); }

    /// \brief Do an oriented box and a ray intersect or touch?
    EIAPI bool intersects( const Ray& _ray, const OBox& _obox )               // TESTED
    {
        // Transform the ray such that the box is centered at the origin
        Ray ray;
        ray.origin = transform( _ray.origin - _obox.center, conjugate(_obox.orientation) );
        ray.direction = transform( _ray.direction, conjugate(_obox.orientation) );
        return intersects( ray, Box(-_obox.halfSides, _obox.halfSides) );
    }

    EIAPI bool intersects( const Ray& _ray, const OBox& _obox, float& _distance )
    {
        // Transform the ray such that the box is centered at the origin
        Ray ray;
        ray.origin = transform( _ray.origin - _obox.center, conjugate(_obox.orientation) );
        ray.direction = transform( _ray.direction, conjugate(_obox.orientation) );
        return intersects( ray, Box(-0.5f*_obox.halfSides, 0.5f*_obox.halfSides), _distance );
    }

    EIAPI bool intersects( const OBox& _obox, const Ray& _ray ) { return intersects( _ray, _obox ); }
    EIAPI bool intersects( const OBox& _obox, const Ray& _ray, float& _distance ) { return intersects( _ray, _obox, _distance ); }

    /// \brief Do a ray and a triangle intersect or touch?
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    /// \param [out,opt] _barycentric The barycentric coordinates of the hit point
    ///     on the triangle.
    /// \return true if there is at least one common point between ray and triangle.
    ///     This point has a ray parameter >= 0 (no negative direction).
    EIAPI bool intersects( const Ray& _ray, const Triangle& _triangle )       // TESTED
    {
        // Möller and Trumbore
        Vec3 e0 = _triangle.v1 - _triangle.v0;
        Vec3 e1 = _triangle.v0 - _triangle.v2;
        // Normal scaled by 2A
        Vec3 normal = cross( e1, e0 );

        float dist2A = dot( normal, _ray.direction );
        Vec3 o = (_triangle.v0 - _ray.origin);
        Vec3 d = cross( _ray.direction, o );

        float barycentricCoord1 = dot( d, e1 ) / dist2A;
        if(barycentricCoord1 < 0.0f) return false;
        float barycentricCoord2 = dot( d, e0 ) / dist2A;
        if(barycentricCoord2 < 0.0f) return false;
        return barycentricCoord1 + barycentricCoord2 <= 1.0f && dot(normal, o) * dist2A >= 0.0f;
    }

    EIAPI bool intersects( const Ray& _ray, const Triangle& _triangle, float& _distance )   // TESTED
    {
        // Möller and Trumbore
        Vec3 e0 = _triangle.v1 - _triangle.v0;
        Vec3 e1 = _triangle.v0 - _triangle.v2;
        // Normal scaled by 2A
        Vec3 normal = cross( e1, e0 );

        float cosT2A = dot( normal, _ray.direction );
        Vec3 o = (_triangle.v0 - _ray.origin);
        Vec3 d = cross( _ray.direction, o );

        float barycentricCoord1 = dot( d, e1 ) / cosT2A;
        if(barycentricCoord1 < -EPSILON || barycentricCoord1 != barycentricCoord1) return false;
        float barycentricCoord2 = dot( d, e0 ) / cosT2A;
        if(barycentricCoord2 < -EPSILON || barycentricCoord2 != barycentricCoord2) return false;
        if(barycentricCoord1 + barycentricCoord2 > 1.0f+EPSILON) return false;

        // Projection to plane. The 2A from normal is canceled out
        _distance = dot( normal, o ) / cosT2A;
        return _distance >= 0.0f;
    }

    EIAPI bool intersects( const Ray& _ray, const FastTriangle& _triangle, float& _distance )   // TESTED
    {
        float cosT = dot( _triangle.normal, _ray.direction );
        float cosT2AInv = 0.5f / (cosT * _triangle.area);
        Vec3 o = (_triangle.v0 - _ray.origin);
        Vec3 d = cross( _ray.direction, o );

        float barycentricCoord1 = -dot( d, _triangle.e02 ) * cosT2AInv;
        if(barycentricCoord1 < -EPSILON || barycentricCoord1 != barycentricCoord1) return false;
        float barycentricCoord2 = dot( d, _triangle.e01 ) * cosT2AInv;
        if(barycentricCoord2 < -EPSILON || barycentricCoord2 != barycentricCoord2) return false;
        if(barycentricCoord1 + barycentricCoord2 > 1.0f+EPSILON) return false;

        // Projection to plane. The 2A from normal is canceled out
        _distance = dot( _triangle.normal, o ) / cosT;
        return _distance >= 0.0f;
    }

    EIAPI bool intersects( const Ray& _ray, const Triangle& _triangle, float& _distance, Vec3& _barycentric )    // TESTED
    {
        // Möller and Trumbore
        Vec3 e0 = _triangle.v1 - _triangle.v0;
        Vec3 e1 = _triangle.v0 - _triangle.v2;
        // Normal scaled by 2A
        Vec3 normal = cross( e1, e0 );

        float cosT2A = dot( normal, _ray.direction );
        Vec3 o = (_triangle.v0 - _ray.origin);
        Vec3 d = cross( _ray.direction, o );

        _barycentric.y = dot( d, e1 ) / cosT2A;
        if(_barycentric.y < -EPSILON) return false;
        _barycentric.z = dot( d, e0 ) / cosT2A;
        if(_barycentric.z < -EPSILON) return false;
        _barycentric.x = 1.0f - (_barycentric.y + _barycentric.z);
        // Do one check on NaN - if any other coordinate is NaN x will be NaN too
        if(_barycentric.x < -EPSILON || _barycentric.x != _barycentric.x) return false;

        // Projection to plane. The 2A from normal is canceled out
        _distance = dot( normal, o ) / cosT2A;
        return _distance >= 0.0f;
    }

    EIAPI bool intersects( const Ray& _ray, const FastTriangle& _triangle, float& _distance, Vec3& _barycentric )  // TESTED
    {
        float cosT = dot( _triangle.normal, _ray.direction );
        float cosT2AInv = 0.5f / (cosT * _triangle.area);
        Vec3 o = (_triangle.v0 - _ray.origin);
        Vec3 d = cross( _ray.direction, o );

        _barycentric.y = -dot( d, _triangle.e02 ) * cosT2AInv;
        if(_barycentric.y < -EPSILON) return false;
        _barycentric.z = dot( d, _triangle.e01 ) * cosT2AInv;
        if(_barycentric.z < -EPSILON) return false;
        _barycentric.x = 1.0f - (_barycentric.y + _barycentric.z);
        // Do one check on NaN - if any other coordinate is NaN x will be NaN too
        if(_barycentric.x < -EPSILON || _barycentric.x != _barycentric.x) return false;

        // Projection to plane. The 2A from normal is canceled out
        _distance = dot( _triangle.normal, o ) / cosT;
        return _distance >= 0.0f;
    }

    EIAPI bool intersects( const Triangle& _triangle, const Ray& _ray )  { return intersects( _ray, _triangle ); }
    EIAPI bool intersects( const Triangle& _triangle, const Ray& _ray, float& _distance )  { return intersects( _ray, _triangle, _distance ); }
    EIAPI bool intersects( const FastTriangle& _triangle, const Ray& _ray, float& _distance )  { return intersects( _ray, _triangle, _distance ); }
    EIAPI bool intersects( const Triangle& _triangle, const Ray& _ray, float& _distance, Vec3& _barycentric )  { return intersects( _ray, _triangle, _distance, _barycentric ); }
    EIAPI bool intersects( const FastTriangle& _triangle, const Ray& _ray, float& _distance, Vec3& _barycentric )  { return intersects( _ray, _triangle, _distance, _barycentric ); }

    /// \brief Do a sphere and a plane intersect or touch?
    /// \return true if there is at least one point in common.
    EIAPI bool intersects( const Sphere& _sphere, const Plane& _plane )       // TESTED
    {
        return abs(dot(_plane.n, _sphere.center) + _plane.d) <= _sphere.radius;
    }

    EIAPI bool intersects( const Plane& _plane, const Sphere& _sphere )       { return intersects(_sphere, _plane); }

    /// \brief Does the sphere touches the triangle?
    /// \return true if the sphere and the triangle have at least one point in common.
    EIAPI bool intersects( const Sphere& _sphere, const Triangle& _triangle ) // TESTED
    {
        Vec3 a = _triangle.v1-_triangle.v0;
        Vec3 b = _triangle.v2-_triangle.v1;
        Vec3 c = _triangle.v2-_triangle.v0;
        Vec3 n = normalize(cross(c, a));

        // Distance to the plane
        float dist = dot(n, _sphere.center-_triangle.v0);
        if(dist > _sphere.radius) return false;

        // The point is inside if it is not on the right hand of any edge
        Vec3 p_v0 = _sphere.center - _triangle.v0;
        Vec3 p_v1 = _sphere.center - _triangle.v1;
        bool outSideA = dot(cross(p_v0, a), n) < 0.0f;
        bool outSideB = dot(cross(p_v1, b), n) < 0.0f;
        bool outSideC = dot(cross(c, p_v0), n) < 0.0f;   // Inverse sign because of use from p-v0

        if(!(outSideA || outSideB || outSideC))
            return true;

        // Minimum is somewhere on the three edges.
        dist =           lensq(p_v0) - sq(dot(p_v0, a)) / lensq(a);
        dist = min(dist, lensq(p_v1) - sq(dot(p_v1, b)) / lensq(b));
        dist = min(dist, lensq(p_v0) - sq(dot(p_v0, c)) / lensq(c));
        return dist <= sq(_sphere.radius);
    }

    EIAPI bool intersects( const Triangle& _triangle, const Sphere& _sphere ) { return intersects(_sphere, _triangle); }

    /// \brief Intersection test between sphere and capsule.
    /// \return true if the sphere and the capsule have at least one point in common.
    EIAPI bool intersects( const Sphere& _sphere, const Capsule& _capsule )   // TESTED
    {
        return distance(_sphere.center, _capsule.seg) <= _capsule.radius + _sphere.radius;
    }

    EIAPI bool intersects( const Capsule& _capsule, const Sphere& _sphere )   { return intersects(_sphere, _capsule); }

    /// \brief Intersection test between point and capsule.
    /// \return true if the point and the capsule have at least one point in common.
    EIAPI bool intersects( const Vec3& _point, const Capsule& _capsule )      // TESTED
    {
        return distance(_point, _capsule.seg) <= _capsule.radius;
    }

    EIAPI bool intersects( const Capsule& _capsule, const Vec3& _point )      { return intersects(_point, _capsule); }

    /// \brief Intersection test between point and frustum.
    /// \return true if the point and the frustum have at least one point in common.
    EIAPI bool intersects( const Vec3& _point, const FastFrustum& _frustum )  // TESTED
    {
        if(distance(_point, _frustum.nf) > 0.0f) return false;
        if(distance(_point, _frustum.l) < 0.0f) return false;
        if(distance(_point, _frustum.r) < 0.0f) return false;
        if(distance(_point, _frustum.b) < 0.0f) return false;
        if(distance(_point, _frustum.t) < 0.0f) return false;
        return true;
    }

    EIAPI bool intersects( const FastFrustum& _frustum, const Vec3& _point )  { return intersects(_point, _frustum); }


    /// \brief Intersection test between ray and frustum.
    /// \param _tmin [in/out] Start of the ray at input. First intersection between
    ///     the frustum and the ray or the ray's input tmin if the ray starts inside.
    /// \param _tmax [in/out] End of the ray at input. Exit point from the frustum
    ///     or the ray's input tmax if the ray interval ends inside.
    /// \return true if the ray is inside the frustum or intersects the boundary.
    EIAPI bool intersects( const Ray& _ray, const Frustum& _frustum, float& _tmin, float& _tmax ) // TESTED
    {
        // Test the near and far plane in 3D (the normal of both is "direction")
        // Perpendicular distances:
        const Vec3 local_origin = _ray.origin - _frustum.apex;
        const float oz = dot(local_origin, _frustum.direction);
        const float dz = dot(_ray.direction, _frustum.direction);
        const float dn = oz - _frustum.n;
        const float df = oz - _frustum.f;
        // Projected distance on the ray.
        const float tn = -sdiv(dn, dz);
        const float tf = -sdiv(df, dz);
        _tmin = max(_tmin, min(tn, tf));
        _tmax = min(_tmax, max(tn, tf));
        if(_tmin > _tmax)
            return false;

        // The segment is now between near and far plane.
        // When we project the ray into the frustum space, we can do 2D line
        // segment tests to get the distances to all other four planes.
        // Note that oz and dz are origin and direction coordinates in local space along z.
        const Vec3 right = cross(_frustum.up, _frustum.direction); // TODO: store with frustum?
        const float ox = dot(local_origin, right);
        const float dx = dot(_ray.direction, right);
        const float oy = dot(local_origin, _frustum.up);
        const float dy = dot(_ray.direction, _frustum.up);

        // Looking at the yz-plane from the side, we have the two lines
        // yzray = (oy,oz) + (dy,dz) * t and yzfrust = (0,0) + (b,f) * s.
        // From yzray == yzfrust we get an equation system with 2 variables which we
        // want to solve for t. The same goes for all other planes
        const float detl = _frustum.f * dx - _frustum.l * dz;
        const float detr = _frustum.f * dx - _frustum.r * dz;
        const float detb = _frustum.f * dy - _frustum.b * dz;
        const float dett = _frustum.f * dy - _frustum.t * dz;
        const float tl = sdiv(_frustum.l * oz - _frustum.f * ox, detl);
        const float tr = sdiv(_frustum.r * oz - _frustum.f * ox, detr);
        const float tb = sdiv(_frustum.b * oz - _frustum.f * oy, detb);
        const float tt = sdiv(_frustum.t * oz - _frustum.f * oy, dett);

        // While it looks like we could, we cannot simply do an interval overlap test like
        // in classic AABox szenarios. Because 4 planes go through one point and produce a
        // negative frustum behind the apex again, the interval [tl,tr] or [tb,tt] could be
        // outside, while the ray is still intersecting.
        // 
        // For each plane we can determine if we hit it from outside or inside.
        // If we hit from the outside it is an entry point. Note that both intersections in
        // each dimension can be entry or exit points, but those we can simply order and take
        // the latest like in common AABox tests.
        // To check for inside/outside we can compute the sign of the determinant for the
        // two direction vectors (dx,dz) and (l,f) (example for left plane).
        // Note that we already computed those determinants for the equations above!
        if (detl >= 0.0f) _tmin = max(_tmin, tl);
        else _tmax = min(_tmax, tl);
        if (detr < 0.0f) _tmin = max(_tmin, tr);
        else _tmax = min(_tmax, tr);
        if (detb >= 0.0f) _tmin = max(_tmin, tb);
        else _tmax = min(_tmax, tb);
        if (dett < 0.0f) _tmin = max(_tmin, tt);
        else _tmax = min(_tmax, tt);
        return _tmin <= _tmax;
    }


    /// \brief Intersection test between sphere and frustum.
    /// \return true if the sphere and the frustum have at least one point in common.
    EIAPI bool intersects( const Sphere& _sphere, const FastFrustum& _frustum )    // TESTED
    {
        // The usual test by comparing all distances to the sphere radius are
        // wrong because they yield false positives when the sphere is outside
        // but intersects all 3 planes in a corner.

        // The following computes the distance of the sphere and early outs
        // if the distance is larger than the radius.
        /*float dsq = sq(_sphere.radius);
        dsq -= sq(max(0.0f, distance(_sphere.center, _frustum.nf)));
        if(dsq < 0.0f) return false;
        dsq -= sq(min(0.0f, distance(_sphere.center, _frustum.l)));
        if(dsq < 0.0f) return false;
        dsq -= sq(min(0.0f, distance(_sphere.center, _frustum.r)));
        if(dsq < 0.0f) return false;
        dsq -= sq(min(0.0f, distance(_sphere.center, _frustum.b)));
        if(dsq < 0.0f) return false;
        dsq -= sq(min(0.0f, distance(_sphere.center, _frustum.t)));
        if(dsq < 0.0f) return false;
        return true;*/

        // Projection approach: track the closest point.
        // Whenever the center is outside and the sphere intersects a plane
        // project the current closest point along the current possible space.
        Vec3 refPoint = _sphere.center;
        Vec3 refData; // Semantic depends on state 0:unused, 1:plane normal, 2:ray direction
        int state = 0;
        float d = distance(refPoint, _frustum.nf);
        if(abs(d) > 0.0f)
        {
            if(abs(d) > _sphere.radius) return false;
            // Intersects but is outside -> project.
            refPoint -= d * _frustum.nf.n;
            refData = _frustum.nf.n;
            state = 1;
        }
        d = distance(refPoint, _frustum.l);
        if(d < 0.0f)
        {
            if(d < -_sphere.radius) return false;
            // Intersects but is outside -> react different on state.
            if(state == 0) {
                refPoint -= d * _frustum.l.n;
                refData = _frustum.l.n;
            } else {
                // Move inside last plane as projection
                Vec3 dir = cross(_frustum.l.n, refData);
                Vec3 projDir = cross(dir, refData);
                refPoint -= projDir * (d / -dot(projDir, _frustum.l.n));
                refData = dir;
            }
            ++state;
        }
        d = distance(refPoint, _frustum.b);
        if(d < 0.0f)
        {
            if(d < -_sphere.radius) return false;
            // Intersects but is outside -> react different on state.
            if(state == 0) {
                refPoint -= d * _frustum.b.n;
                refData = _frustum.b.n;
            } else if(state == 1) {
                // Move inside last plane as projection
                Vec3 dir = cross(_frustum.b.n, refData);
                Vec3 projDir = cross(dir, refData);
                refPoint -= projDir * (d / -dot(projDir, _frustum.b.n));
                refData = dir;
            } else {
                // Find ray plane intersection
                refPoint += refData * (d / -dot(refData, _frustum.b.n));
                return lensq(refPoint - _sphere.center) <= _sphere.radius * _sphere.radius;
            }
            ++state;
        }
        d = distance(refPoint, _frustum.r);
        if(d < 0.0f)
        {
            if(d < -_sphere.radius) return false;
            // Intersects but is outside -> react different on state.
            if(state == 0) {
                refPoint -= d * _frustum.r.n;
                refData = _frustum.r.n;
            } else if(state == 1) {
                // Move inside last plane as projection
                Vec3 dir = cross(_frustum.r.n, refData);
                Vec3 projDir = cross(dir, refData);
                refPoint -= projDir * (d / -dot(projDir, _frustum.r.n));
                refData = dir;
            } else {
                // Find ray plane intersection
                refPoint += refData * (d / -dot(refData, _frustum.r.n));
                return lensq(refPoint - _sphere.center) <= _sphere.radius * _sphere.radius;
            }
            ++state;
        }
        d = distance(refPoint, _frustum.t);
        if(d < 0.0f)
        {
            if(d < -_sphere.radius) return false;
            // Intersects but is outside -> react different on state.
            if(state == 0) {
                refPoint -= d * _frustum.t.n;
                refData = _frustum.t.n;
            } else if(state == 1) {
                // Move inside last plane as projection
                Vec3 dir = cross(_frustum.t.n, refData);
                Vec3 projDir = cross(dir, refData);
                refPoint -= projDir * (d / -dot(projDir, _frustum.t.n));
                refData = dir;
            } else {
                // Find ray plane intersection
                refPoint += refData * (d / -dot(refData, _frustum.t.n));
                return lensq(refPoint - _sphere.center) <= _sphere.radius * _sphere.radius;
            }
            ++state;
        }

        if(state == 2)
            return lensq(refPoint - _sphere.center) <= _sphere.radius * _sphere.radius;
        return true;
    }

    EIAPI bool intersects( const FastFrustum& _frustum, const Sphere& _sphere ) { return intersects(_sphere, _frustum); }

    namespace details {
        EIAPI bool separates( const Vec3& _dir, const Vec3* _box, const Vec3* _fru )
        {
            // Min and max values for the projected vertices of box and frustum
            float bmin = dot(_dir, _box[0]), fmin = dot(_dir, _fru[0]);
            float bmax = bmin, fmax = fmin;
            for(int i=1; i<8; ++i)
            {
                float d = dot(_dir, _box[i]);
                if(d < bmin) bmin = d;
                else if(d > bmax) bmax = d;
                d = dot(_dir, _fru[i]);
                if(d < fmin) fmin = d;
                else if(d > fmax) fmax = d;
            }
            return bmin > fmax || bmax < fmin;
        }
    }

    /// \brief Intersection test between box and frustum.
    /// \return true if the box and the frustum have at least one point in common.
    EIAPI bool intersects( const Box& _box, const FastFrustum& _frustum )       // TESTED
    {
        // The usual box test (first part) results in false positives if a box
        // is between 0 or 1 plane pairs and intersects the others (3 or 2)
        // outside the frustum.
        // One failsafe algorithm is to find a separating plane. If none can be
        // found the two convex objects intersect. For polyhedron the number
        // of possible plane normals is finite (a side of one of the objects or
        // a cross product between two edges of the different objects.
        // http://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf

        // Test against box sides first (cheaper because of axis alignment)
        int out[6] = {0}; // track point status for early out
        for(int i=0; i<8; ++i)
        {
            int in = 0;
            if( _frustum.vertices[i].x < _box.min.x ) ++out[0]; else ++in;
            if( _frustum.vertices[i].x > _box.max.x ) ++out[1]; else ++in;
            if( _frustum.vertices[i].y < _box.min.y ) ++out[2]; else ++in;
            if( _frustum.vertices[i].y > _box.max.y ) ++out[3]; else ++in;
            if( _frustum.vertices[i].z < _box.min.z ) ++out[4]; else ++in;
            if( _frustum.vertices[i].z > _box.max.z ) ++out[5]; else ++in;
            // One vertex inside the box?
            if(in == 6) return true;
        }
        // Fully separated by one of the sides?
        for(int i=0; i<6; ++i) {
            if(out[i] == 8) return false;
            out[i] = 0;
        }

        // Test against frustum sides
        Vec3 boxv[8] = {_box.min, Vec3(_box.min.x, _box.min.y, _box.max.z),
                        Vec3(_box.min.x, _box.max.y, _box.min.z), Vec3(_box.min.x, _box.max.y, _box.max.z),
                        Vec3(_box.max.x, _box.min.y, _box.min.z), Vec3(_box.max.x, _box.min.y, _box.max.z),
                        Vec3(_box.max.x, _box.max.y, _box.min.z), _box.max};
        for(int i=0; i<8; ++i)
        {
            int in = 0;
            float d = distance(_frustum.nf, boxv[i]);
            if( d < 0.0f ) ++out[0];
            else if( d > 0.0f ) ++out[1];
            else in+=2;
            if( distance(_frustum.l, boxv[i]) < 0.0f ) ++out[2]; else ++in;
            if( distance(_frustum.r, boxv[i]) < 0.0f ) ++out[3]; else ++in;
            if( distance(_frustum.b, boxv[i]) < 0.0f ) ++out[4]; else ++in;
            if( distance(_frustum.t, boxv[i]) < 0.0f ) ++out[5]; else ++in;
            if(in == 6) return true;
        }
        for(int i=0; i<6; ++i)
            if(out[i] == 8) return false;

        // Test against the other possible planes
        Vec3 n; // current plane normal
        n = cross(Vec3(1.0f, 0.0f, 0.0f), _frustum.vertices[1] - _frustum.vertices[0]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(1.0f, 0.0f, 0.0f), _frustum.vertices[2] - _frustum.vertices[0]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(1.0f, 0.0f, 0.0f), _frustum.vertices[4] - _frustum.vertices[0]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(1.0f, 0.0f, 0.0f), _frustum.vertices[5] - _frustum.vertices[1]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(1.0f, 0.0f, 0.0f), _frustum.vertices[6] - _frustum.vertices[2]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(1.0f, 0.0f, 0.0f), _frustum.vertices[7] - _frustum.vertices[3]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;

        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[1] - _frustum.vertices[0]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[2] - _frustum.vertices[0]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[4] - _frustum.vertices[0]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[5] - _frustum.vertices[1]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[6] - _frustum.vertices[2]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[7] - _frustum.vertices[3]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;

        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[1] - _frustum.vertices[0]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[2] - _frustum.vertices[0]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[4] - _frustum.vertices[0]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[5] - _frustum.vertices[1]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[6] - _frustum.vertices[2]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[7] - _frustum.vertices[3]);
        if(details::separates(n, boxv, _frustum.vertices)) return false;

        return true;
    }

    EIAPI bool intersects( const FastFrustum& _frustum, const Box& _box )     { return intersects(_box, _frustum); }

    /// \brief Intersection test between point and tetrahedron.
    /// \return true if the point and the tetrahedron have at least one point in common.
    EIAPI bool intersects( const Vec3& _point, const Tetrahedron& _tetrahedron )    // TESTED
    {
        // Check if the point is on the inner side of each plane.
        // First determine edges for the cross products to find normals
        Vec3 e01 = _tetrahedron.v1 - _tetrahedron.v0;
        Vec3 e02 = _tetrahedron.v2 - _tetrahedron.v0;
        Vec3 e03 = _tetrahedron.v3 - _tetrahedron.v0;
        Vec3 n = cross(e01, e02);
        // Check if the fourth point of the tetrahedron and _point are on the same side.
        // Reuse as much as possible operations.
        Vec3 p = _point - _tetrahedron.v0;
        float dp = dot(n, p);
        float dt = dot(n, e03);
        if(dp * dt < 0.0f) return false;
        Vec3 e13 = _tetrahedron.v3 - _tetrahedron.v1;
        n = cross(e13, e01);
        dp = dot(n, p);
        dt = dot(n, e02);
        if(dp * dt < 0.0f) return false;
        Vec3 e23 = _tetrahedron.v3 - _tetrahedron.v2;
        n = cross(e23, e02);
        dp = dot(n, p);
        dt = dot(n, e01);
        if(dp * dt < 0.0f) return false;
        n = cross(e23, e13);
        dp = dot(n, _tetrahedron.v3 - _point);
        dt = dot(n, e03);
        return dp * dt >= 0.0f;

        // Much slower alternative (could be faster for batch operations when
        // J^-1 or JLUp is provided):
        // Find a transformation matrix that projects the tetrahedron into a
        // normalized one: http://www.iue.tuwien.ac.at/phd/hollauer/node29.html
        /*Mat3x3 J = axis(_tetrahedron.v1 - _tetrahedron.v0, _tetrahedron.v2 - _tetrahedron.v1, _tetrahedron.v3 - _tetrahedron.v1);
        // Now p - v0 = J * x where x is p in normalized coordinates
        Mat3x3 JLU; Vec<uint, 3> permutation;
        decomposeLUp(J, JLU, permutation);
        Vec3 normp = solveLUp(JLU, permutation, _point - _tetrahedron.v0);
        if(any(normp < 0.0f)) return false;
        return dot(normp, Vec3(1.0f)) <= 0.577350269f;*/
    }

    EIAPI bool intersects( const Tetrahedron& _tetrahedron, const Vec3& _point ) { return intersects(_point, _tetrahedron); }

    /// \brief Intersection test between triangle and box (based on SAT).
    /// \return true if the triangle and the box have at least one point in common.
    EIAPI bool intersects( const Triangle& _triangle, const Box& _box )          // TESTED
    {
        // Implementation based on SAT.
        // Like "Fast 3D Triangle-Box Overlap Testing" from Tomas Akenine-Möller.

        // *** Case: A box face separates the triangle. Since the box is axis
        // aligned, a projection is simply selecting one of the coordinates.
        if(any(less(max(_triangle.v0, _triangle.v1, _triangle.v2), _box.min)))
            return false;
        if(any(greater(min(_triangle.v0, _triangle.v1, _triangle.v2), _box.max)))
            return false;

        // Some precomputations for upcoming cases.
        Vec3 e0 = _triangle.v1 - _triangle.v0;
        Vec3 e1 = _triangle.v2 - _triangle.v0;
        Vec3 triNormal = cross(e0, e1); // No need for normalization!
        Vec3 boxCenter = center(_box);
        Vec3 boxSideHalf = (_box.max - _box.min) * 0.5f;
        Vec3 v0 = _triangle.v0 - boxCenter;

        // *** Case: The triangle plane separates the box.
        // Instead looping over all eight corners we only compute the minimum and
        // the maximum projected coordinate. If they have different signs the triangle
        // plane (not the triangle) intersects the box.
        // The minimum and maximum depend on the plane normal direction.
        // The combination which maximizes the dot product is the larger coord (max)
        // if the coord in the normal is positive and the smaller one if it is negative.
        // The same goes for minimizing.
        float triPlaneOffset = -dot(triNormal, v0);
        float projMax = dot(abs(triNormal), boxSideHalf);
        if(abs(triPlaneOffset) > projMax) // Triangle plane is farther away then the maximum possible box coordinate -> separates
            return false;

        // *** Case: Separating plane is maybe defined by one of the 9
        // box<->triangle edge combinations.
        // We can do some optimizations because of axis alignment (cross products with
        // unit axis vectors).
        Vec3 e2 = _triangle.v2 - _triangle.v1;
        Vec3 v1 = _triangle.v1 - boxCenter;
        Vec3 v2 = _triangle.v2 - boxCenter;
        const Vec3 planeDirs[9] = {
            Vec3(0.0f, e0.z, -e0.y),
            Vec3(0.0f, e1.z, -e1.y),
            Vec3(0.0f, e2.z, -e2.y),
            Vec3(-e0.z, 0.0f, e0.x),
            Vec3(-e1.z, 0.0f, e1.x),
            Vec3(-e2.z, 0.0f, e2.x),
            Vec3(e0.y, -e0.x, 0.0f),
            Vec3(e1.y, -e1.x, 0.0f),
            Vec3(e2.y, -e2.x, 0.0f)
        };
        for(int i = 0; i < 9; ++i)
        {
            projMax = dot(abs(planeDirs[i]), boxSideHalf);
            float projTri0 = dot(planeDirs[i], v0);
            float projTri1 = dot(planeDirs[i], v1);
            float projTri2 = dot(planeDirs[i], v2);
            if(min(projTri0, projTri1, projTri2) > projMax || max(projTri0, projTri1, projTri2) < -projMax)
                return false;
        }

        // No separating axis found
        return true;
    }

    EIAPI bool intersects( const Box& _box, const Triangle& _triangle ) { return intersects(_triangle, _box); }

    /// \brief Intersection test between triangle and oriented box (based on SAT).
    /// \return true if the triangle and the box have at least one point in common.
    EIAPI bool intersects( const Triangle& _triangle, const OBox& _obox )     // TESTED
    {
        // Transform triangle into local box space
        Triangle alignedTriangle;
        Mat3x3 rotation(conjugate(_obox.orientation));
        alignedTriangle.v0 = rotation * (_triangle.v0 - _obox.center);
        alignedTriangle.v1 = rotation * (_triangle.v1 - _obox.center);
        alignedTriangle.v2 = rotation * (_triangle.v2 - _obox.center);
       // return intersects(alignedTriangle, Box(-_obox.halfSides, _obox.halfSides));

        // For details see Triangle <-> Box intersection. The following is
        // copied and simplified with _obox.min = -_obox.max = -_obox.halfSides.
        if(any(less(max(alignedTriangle.v0, alignedTriangle.v1, alignedTriangle.v2), -_obox.halfSides)))
            return false;
        if(any(greater(min(alignedTriangle.v0, alignedTriangle.v1, alignedTriangle.v2), _obox.halfSides)))
            return false;

        Vec3 e0 = alignedTriangle.v1 - alignedTriangle.v0;
        Vec3 e1 = alignedTriangle.v2 - alignedTriangle.v0;
        Vec3 triNormal = cross(e0, e1); // No need for normalization!
        float triPlaneOffset = dot(triNormal, alignedTriangle.v0);

        float projMax = dot(abs(triNormal), _obox.halfSides);
        if(abs(triPlaneOffset) > projMax) // Triangle plane is farther away then the maximum possible box coordinate -> separates
            return false;

        Vec3 e2 = alignedTriangle.v2 - alignedTriangle.v1;
        const Vec3 planeDirs[9] = {
            Vec3(0.0f, e0.z, -e0.y),
            Vec3(0.0f, e1.z, -e1.y),
            Vec3(0.0f, e2.z, -e2.y),
            Vec3(-e0.z, 0.0f, e0.x),
            Vec3(-e1.z, 0.0f, e1.x),
            Vec3(-e2.z, 0.0f, e2.x),
            Vec3(e0.y, -e0.x, 0.0f),
            Vec3(e1.y, -e1.x, 0.0f),
            Vec3(e2.y, -e2.x, 0.0f)
        };
        for(int i = 0; i < 9; ++i)
        {
            projMax = dot(abs(planeDirs[i]), _obox.halfSides);
            float projTri0 = dot(planeDirs[i], alignedTriangle.v0);
            float projTri1 = dot(planeDirs[i], alignedTriangle.v1);
            float projTri2 = dot(planeDirs[i], alignedTriangle.v2);
            if(min(projTri0, projTri1, projTri2) > projMax || max(projTri0, projTri1, projTri2) < -projMax)
                return false;
        }

        // No separating axis found
        return true;
    }

    EIAPI bool intersects( const OBox& _obox, const Triangle& _triangle ) { return intersects(_triangle, _obox); }

    /// \brief Intersection test between plane and box.
    /// \return true if the plane and the box have at least one point in common.
    EIAPI bool intersects( const Plane& _plane, const Box& _box )
    {
        // Instead looping over all eight corners we only compute the minimum and
        // the maximum projected coordinate. If they have different signs the
        // plane intersects the box.
        // The minimum and maximum depend on the plane normal direction.
        // The combination which maximizes the dot product is the larger coord (max)
        // if the normal is positive and the smaller one if it is negative in
        // the respective component.
        // The same goes for minimizing.
        /*float projMin = dot(_plane.n, Vec3(_plane.n.x > 0 ? _box.min.x : _box.max.x,
            _plane.n.y > 0 ? _box.min.y : _box.max.y,
            _plane.n.z > 0 ? _box.min.z : _box.max.z)) + _plane.d;
        float projMax = dot(_plane.n, Vec3(_plane.n.x < 0 ? _box.min.x : _box.max.x,
            _plane.n.y < 0 ? _box.min.y : _box.max.y,
            _plane.n.z < 0 ? _box.min.z : _box.max.z)) + _plane.d;
        return projMin * projMax <= 0.0f; // Same sign -> separates*/

        // Weird: the following one is faster.
        float boxLocalPlaneOffset = _plane.d + dot(_plane.n, center(_box));
        float projMax = dot(abs(_plane.n), _box.max - _box.min) * 0.5f;
        return abs(boxLocalPlaneOffset) <= projMax;
    }

    EIAPI bool intersects( const Box& _box, const Plane& _plane ) { return intersects(_plane, _box); }

    /// \brief Intersection test between plane and oriented box.
    /// \return true if the plane and the oriented box have at least one point in common.
    EIAPI bool intersects( const Plane& _plane, const OBox& _obox )
    {
        float boxLocalPlaneOffset = _plane.d + dot(_plane.n, _obox.center);
        Vec3 boxLocalPlaneNormal = transform(_plane.n, conjugate(_obox.orientation));
        float projMax = dot(abs(boxLocalPlaneNormal), _obox.halfSides);
        return abs(boxLocalPlaneOffset) <= projMax; // Plane is farther away then the maximum possible box coordinate -> separates
    }

    EIAPI bool intersects( const OBox& _obox, const Plane& _plane ) { return intersects(_plane, _obox); }

    /// \brief Intersection test between cone and point.
    EIAPI bool intersects( const Vec3& _point, const Cone& _cone )
    {
        // The projection of the point to the central ray gives a distance dp.
        // With dp the radius of the cone can be computed and compared to the
        // distance d2 between point and projected point.
        Vec3 oToP = _point - _cone.centralRay.origin;
        float dp = dot(oToP, _cone.centralRay.direction);
        if(dp < 0.0f || dp > _cone.height) return false; // Early out: behind cone origin or base
        // Use Pythagoras to avoid computation of d2 and the projected point.
        float dOPSq = dot(oToP, oToP);
        //float coneRadius = _cone.tanTheta * dp;
        // Maybe use an alternative for numerical more stable forms?
        //     dOPSq - dp * dp <= sq(_cone.tanTheta * dp)
        // <=> dOPSq - dp * dp <= sq(_cone.tanTheta) * dp * dp
        // <=> dOPSq / (dp * dp) - 1 <= sq(_cone.tanTheta)
        // <=> (sqrt(dOPSq) / dp - 1) * (sqrt(dOPSq) / dp + 1) <= sq(_cone.tanTheta)
        return dOPSq - dp * dp <= sq(_cone.tanTheta * dp);
    }

    EIAPI bool intersects( const Cone& _cone, const Vec3& _point ) { return intersects(_point, _cone); }
    
    EIAPI bool intersects( const Vec3& _point, const FastCone& _cone )
    {
        // Given the cos(theta)^2 the point can be compared more
        // easily by computing its cosine over the dot product.
        // This would require a square root in the length computation.
        // Therefore it is cheaper to compute and compare cosine squares.
        Vec3 oToP = _point - _cone.centralRay.origin;
        float dp = dot(oToP, _cone.centralRay.direction);
        if(dp < 0.0f || dp > _cone.height) return false; // Early out: behind cone origin or base
        return dp * dp >= _cone.cosThetaSq * dot(oToP, oToP);
    }

    EIAPI bool intersects( const FastCone& _cone, const Vec3& _point ) { return intersects(_point, _cone); }

    /// \brief Intersection test between cone and triangle.
    EIAPI bool intersects( const Triangle& _triangle, const Cone& _cone )
    {
        // Idea: first get the three closest points between each of the edges
        // and the cone central segment. If at least one of them is inside we
        // can see a part of the triangle. If all are outside the central ray
        // must hit the triangle or the triangle is fully outside.

        float cosThetaSq = 1.0f / (1.0f + _cone.tanTheta * _cone.tanTheta);

        /*// Segment - segment distance: http://geomalgorithms.com/a07-_distance.html
        Vec3 e0 = _triangle.v1 - _triangle.v0;
        Vec3 vToO = _cone.centralRay.origin - _triangle.v0;
        float a = dot(e0, e0);
        float b = dot(e0, _cone.centralRay.direction);
        float p = dot(vToO, e0);
        float q = dot(vToO, _cone.centralRay.direction);
        float n = a - b * b;
        if(n != 0.0f) // Ignore the parallel case - another edge should be non parallel
        {
            float s = (b * q -     p) / -n; // Ray parameter on edge
            float t = (a * q - b * p) / -n; // Ray parameter on cone
            // Clamp s, such that the point is not behind the 
            Vec3 closestOnEdge = _triangle.v0 + saturate(s) * e0;
            // return dOPSq - dp * dp <= sq(_cone.tanTheta * dp);
            if(lensq(saturate(s) * e0 - vToO - clamp(t, 0.0f, _cone.height) * _cone.centralRay.direction))
        }*/


        // [1] https://www.geometrictools.com/Documentation/IntersectionTriangleCone.pdf
        // First check if any of the vertices is inside. If not store on which side of
        // the cone half space they are.
        Vec3 p[3];
        float dp[3];
        float pp[3];
        for(int i = 0; i < 3; ++i)
        {
            p[i] = _triangle.v(i) - _cone.centralRay.origin;
            dp[i] = dot(p[i], _cone.centralRay.direction);
            pp[i] = dot(p[i], p[i]);
            if(dp[i] >= 0 && dp[i] <= _cone.height && dp[i] * dp[i] >= cosThetaSq * pp[i])
                return true;
        }

        // If all points are on the back side, the triangle cannot be hit
        if(dp[0] < 0 && dp[1] < 0 && dp[2] < 0) return false;
        // If all points are farther away than the cone height, the triangle cannot be hit
        if(dp[0] > _cone.height && dp[1] > _cone.height && dp[2] > _cone.height) return false;

        // Now, test if any of the sides intersects the cone.
        for(int i0 = 0; i0 < 3; ++i0)
        {
            int i1 = (i0 + 1) % 3;
            if((dp[i0] >= 0 || dp[i1] >= 0)
                && (dp[i0] <= _cone.height || dp[i1] <= _cone.height)) // If both vertices are on the wrong side there is no intersection.
            {
                // Compute the intersection of the edge with the cone.
                // Cone relative parametric edge: Q(t) = (v0 - o) + t * (v1 - v0)
                //                                     = p + t * e
                // Parametric double cone: dot(d, p/|p|) = cos θ
                // Setting p = Q(t), multiplying with |Q(t)|, squaring
                // and putting all to the same side yields:
                //    dot(d, Q(t))^2 - cos^2 θ * |Q(t)|^2 = 0
                // which must be solved (c2 t^2 + 2 c1 t + c0 = 0).
                Vec3 e = _triangle.v(i1) - _triangle.v(i0);
                float de = dot(e, _cone.centralRay.direction);
                float ee = dot(e, e);
                float c2 = de * de - cosThetaSq * ee;
                if(c2 < 0.0f) // If c2 isn't negative there cannot be an intersection (see [1]).
                {
                    float pe = dot(p[i0], e);
                    float c0 = dp[i0] * dp[i0] - cosThetaSq * pp[i0];         // Regularization as in the geometrictools doc?
                    float c1 = dp[i0] * de - cosThetaSq * pe;
                    float t0, t1;
                    if(!solveSquarePoly(c2, c1 * 2.0f, c0, t0, t1)
                        || t0 < 0.0f || t1 > 1.0f) continue;
                    eiAssert(t0 <= 1.0f || t1 >= 0.0f, "At least one solution must be in [0,1].");
                    // The edge intersects the infinite double cone. Make sure that
                    // one of the intersection point is really in the cone.
                    // Since we know the point is in the cone we do not need the full test
                    // and the projected distance comparison is sufficient.
                    float dx = dot(p[i0] + e * saturate(t0), _cone.centralRay.direction);
                    if(dx >= 0.0f && dx <= _cone.height) return true;
                    dx = dot(p[i0] + e * saturate(t1), _cone.centralRay.direction);
                    if(dx >= 0.0f && dx <= _cone.height) return true;
                }
                // Original conditions from [1] without explicit intersection test.
                // Does not include the boundary conditions of our finite cone.
                /*if(c2 < 0.0f && c1*c1 >= c0*c2)
                {
                    if(onConeSide[i0] && onConeSide[i1] && 0 <= c1 && c1 <= -c2) return true;
                    bool opCon = c2 * dp[i0] <= c1 * de;
                    if(onConeSide[i0] && opCon && 0 <= c1) return true;
                }*/
            }
        }

        // Test if the center ray hits the triangle after Möller and Trumbore
        Vec3 e0 = _triangle.v1 - _triangle.v0;
        Vec3 e1 = _triangle.v0 - _triangle.v2;
        // Normal scaled by 2A
        Vec3 normal = cross( e1, e0 );

        float dist2A = dot( normal, _cone.centralRay.direction );
        Vec3 d = cross( _cone.centralRay.direction, p[0] );

        float barycentricCoord1 = dot( d, e1 ) / dist2A;
        if(barycentricCoord1 < 0.0f) return false;
        float barycentricCoord2 = dot( d, e0 ) / dist2A;
        if(barycentricCoord2 < 0.0f) return false;
        if(barycentricCoord1 + barycentricCoord2 > 1.0f) return false;

        // Projection to plane. The 2A from normal is canceled out
        float distance = dot( normal, p[0] ) / dist2A;
        return distance >= 0.0f && distance <= _cone.height;
    }

    EIAPI bool intersects( const Cone& _cone, const Triangle& _triangle ) { return intersects(_triangle, _cone); }
}
