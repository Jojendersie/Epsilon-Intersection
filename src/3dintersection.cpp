#include "ei/3dintersection.hpp"

namespace ei {

    // ************************************************************************* //
    // DISTANCE METHOD
    // ************************************************************************* //
    float distance(const Vec3& _point0, const Vec3& _point1)
    {
        return len(_point1 - _point0);
    }

    // ************************************************************************* //
    float distance(const Vec3& _point, const Sphere& _sphere)
    {
        return distance(_sphere.center, _point) - _sphere.radius;
    }

    // ************************************************************************* //
    float distance(const Vec3& _point, const Capsule& _capsule)
    {
        return distance(_point, _capsule.seg) - _capsule.radius;
    }

    // ************************************************************************* //
    float distance(const Vec3& _point, const Plane& _plane)
    {
        eiAssert(approx(1.0f, len(_plane.n)), "The plane is not normalized!");
        return dot(_plane.n, _point) + _plane.d;
    }

    // ************************************************************************* //
    float distance(const Vec3& _point, const DOP& _dop)
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
    float distance(const Sphere& _sphere, const Plane& _plane)
    {
        eiAssert(approx(1.0f, len(_plane.n)), "The plane is not normalized!");
        eiAssert(_sphere.radius >= 0.0f, "Invalid sphere!");
        float d = dot(_plane.n, _sphere.center) + _plane.d;
        // Keep sign for the plane distances
        if(d < 0.0f) return min(d + _sphere.radius, 0.0f);
        else return max(d - _sphere.radius, 0.0f);
    }

    // ************************************************************************* //
    float distance(const Sphere& _sphere, const Segment& _segment)
    {
        return max(0.0f, distance(_sphere.center, _segment) - _sphere.radius);
    }

    // ************************************************************************* //
    float distance(const Sphere& _sphere, const Capsule& _capsule)
    {
        return max(0.0f, distance(_sphere.center, _capsule.seg) - _capsule.radius - _sphere.radius);
    }

    // ************************************************************************* //
    float distance(const Capsule& _capsule0, const Capsule& _capsule1)
    {
        return max(0.0f, distance(_capsule0.seg, _capsule1.seg) - _capsule0.radius - _capsule1.radius);
    }

    // ************************************************************************* //
    // INTERSECTION METHODS
    // ********************************************************************* //
    bool intersects( const Sphere& _sphere0, const Sphere& _sphere1 )
    {
        return lensq(_sphere1.center-_sphere0.center) <= sq(_sphere0.radius + _sphere1.radius);
    }

    // ********************************************************************* //
    bool intersects( const Vec3& _point, const Sphere& _sphere )
    {
        return lensq(_point-_sphere.center) <= sq(_sphere.radius);
    }

    // ********************************************************************* //
    bool intersects( const Sphere& _sphere, const Box& _box )
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

    // ********************************************************************* //
    bool intersects( const Box& _box0, const Box& _box1 )
    {
        eiAssert(all(_box0.min <= _box0.max), "Box0 is invalid.");
        eiAssert(all(_box1.min <= _box1.max), "Box1 is invalid.");
        // There must be an intersection if the sum of side length is larger
        // than that of the bounding box.
//        return all((max(_box0.max, _box1.max) - min(_box0.min, _box1.min))
//            <= ((_box0.max - _box0.min) + (_box1.max - _box1.min)));
        // Non-vector variant is faster
        return (max(_box0.max.x, _box1.max.x) - min(_box0.min.x, _box1.min.x)) <= ((_box0.max.x - _box0.min.x) + (_box1.max.x - _box1.min.x))
            && (max(_box0.max.y, _box1.max.y) - min(_box0.min.y, _box1.min.y)) <= ((_box0.max.y - _box0.min.y) + (_box1.max.y - _box1.min.y))
            && (max(_box0.max.z, _box1.max.z) - min(_box0.min.z, _box1.min.z)) <= ((_box0.max.z - _box0.min.z) + (_box1.max.z - _box1.min.z));
    }

    // ********************************************************************* //
    bool intersects( const Capsule& _capsule0, const Capsule& _capsule1 )
    {
        return distance(_capsule0.seg, _capsule1.seg) <= (_capsule0.radius + _capsule1.radius);
    }

    // ********************************************************************* //
    bool intersects( const Vec3& _point, const Ellipsoid& _ellipsoid )
    {
        // Use ellipsoid equation.
        //return lensq(Vec<double,3>((_point - _ellipsoid.center) / _ellipsoid.radii)) <= 1.0;
        //return lensq((_point - _ellipsoid.center) / _ellipsoid.radii) <= 1.0;
        Vec3 o = (_point - _ellipsoid.center);
        if( _ellipsoid.radii.x != 0.0f ) o.x /= _ellipsoid.radii.x; else o.x *= 3.402823466e+38f;
        if( _ellipsoid.radii.y != 0.0f ) o.y /= _ellipsoid.radii.y; else o.y *= 3.402823466e+38f;
        if( _ellipsoid.radii.z != 0.0f ) o.z /= _ellipsoid.radii.z; else o.z *= 3.402823466e+38f;
        return double(o.x*o.x) + double(o.y*o.y) + double(o.z*o.z) <= 1.0;
        // The following is a try to get a more staple solution without
        // addition, but it is less stable in the end.
        //Vec3 t = (_point - _ellipsoid.center) / _ellipsoid.radii;
        //return exp(t.x*t.x)*exp(t.y*t.y)*exp(t.z*t.z) <= 2.718281828f;
    }

    // ********************************************************************* //
    bool intersects( const Ray& _ray, const Ellipsoid& _ellipsoid )
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

    // ********************************************************************* //
    bool intersects( const Ray& _ray, const Ellipsoid& _ellipsoid, float& _distance )
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

    // ********************************************************************* //
    bool intersects( const Ray& _ray, const Box& _box )
    {
        float t0 = (_box.min.x - _ray.origin.x) / _ray.direction.x;
        float t1 = (_box.max.x - _ray.origin.x) / _ray.direction.x;
        float tmin = min(t0, t1);
        float tmax = max(t0, t1);
        if(tmax < 0.0f) return false;
        t0 = (_box.min.y - _ray.origin.y) / _ray.direction.y;
        t1 = (_box.max.y - _ray.origin.y) / _ray.direction.y;
        tmin = max(tmin, min(t0, t1));
        tmax = min(tmax, max(t0, t1));
        if(tmax < 0.0f || tmin > tmax) return false;
        t0 = (_box.min.z - _ray.origin.z) / _ray.direction.z;
        t1 = (_box.max.z - _ray.origin.z) / _ray.direction.z;
        tmin = max(tmin, min(t0, t1));
        tmax = min(tmax, max(t0, t1));
        return (tmax >= 0.0f) && (tmin <= tmax);
        /*Vec3 tbot = (_box.min - _ray.origin) / _ray.direction;
        Vec3 ttop = (_box.max - _ray.origin) / _ray.direction;
        Vec3 tmin = min(ttop, tbot);
        Vec3 tmax = max(ttop, tbot);
        float firstHit = max(0.0f, max(tmin));
        float lastHit = min(tmax);
        return firstHit <= lastHit;*/
    }

    // ********************************************************************* //
    bool intersects( const Ray& _ray, const Box& _box, float& _distance )
    {
        float t0 = (_box.min.x - _ray.origin.x) / _ray.direction.x;
        float t1 = (_box.max.x - _ray.origin.x) / _ray.direction.x;
        float tmin = min(t0, t1);
        float tmax = max(t0, t1);
        if(tmax < 0.0f) return false;
        t0 = (_box.min.y - _ray.origin.y) / _ray.direction.y;
        t1 = (_box.max.y - _ray.origin.y) / _ray.direction.y;
        tmin = max(tmin, min(t0, t1));
        tmax = min(tmax, max(t0, t1));
        if(tmax < 0.0f || tmin > tmax) return false;
        t0 = (_box.min.z - _ray.origin.z) / _ray.direction.z;
        t1 = (_box.max.z - _ray.origin.z) / _ray.direction.z;
        tmin = max(tmin, min(t0, t1));
        tmax = min(tmax, max(t0, t1));
        _distance = max(tmin, 0.0f);//tmin < 0.0f ? tmax : tmin;
        return tmin <= tmax;
    }

    // ********************************************************************* //
    bool intersects( const Ray& _ray, const Box& _box, float& _distance, HitSide& _side )
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
            tside = min0 < 0.0f ? HitSide::Y_POS : HitSide::Y_NEG;
        } else {
            tmax = (_box.min.y - _ray.origin.y) / _ray.direction.y;
            tmin = (_box.max.y - _ray.origin.y) / _ray.direction.y;
            tside = min0 < 0.0f ? HitSide::Y_NEG : HitSide::Y_POS;
        }

        if( (min0 > tmax) || (tmin > max0) )
            return false;
        if( tmin > min0 ) { min0 = tmin; _side = tside; }
        if( tmax < max0 ) max0 = tmax;

        if( _ray.direction.z >= 0 ) {
            tmin = (_box.min.z - _ray.origin.z) / _ray.direction.z;
            tmax = (_box.max.z - _ray.origin.z) / _ray.direction.z;
            tside = min0 < 0.0f ? HitSide::Z_POS : HitSide::Z_NEG;
        } else {
            tmax = (_box.min.z - _ray.origin.z) / _ray.direction.z;
            tmin = (_box.max.z - _ray.origin.z) / _ray.direction.z;
            tside = min0 < 0.0f ? HitSide::Z_NEG : HitSide::Z_POS;
        }

        if( (min0 > tmax) || (tmin > max0) )
            return false;
        if( 0.0f > max0 || 0.0f > tmax ) return false;
        if( tmin > min0 ) _side = tside;

        _distance = max(min0, tmin, 0.0f);
        return true;
    }

    // ********************************************************************* //
    bool intersects( const Ray& _ray, const Box& _box, float& _distance, float& _distanceExit )
    {
        float t0 = (_box.min.x - _ray.origin.x) / _ray.direction.x;
        float t1 = (_box.max.x - _ray.origin.x) / _ray.direction.x;
        float tmin = min(t0, t1);
        float tmax = max(t0, t1);
        if(tmax < 0.0f) return false;
        t0 = (_box.min.y - _ray.origin.y) / _ray.direction.y;
        t1 = (_box.max.y - _ray.origin.y) / _ray.direction.y;
        tmin = max(tmin, min(t0, t1));
        tmax = min(tmax, max(t0, t1));
        if(tmax < 0.0f || tmin > tmax) return false;
        t0 = (_box.min.z - _ray.origin.z) / _ray.direction.z;
        t1 = (_box.max.z - _ray.origin.z) / _ray.direction.z;
        tmin = max(tmin, min(t0, t1));
        tmax = min(tmax, max(t0, t1));
        _distance = max(tmin, 0.0f);
        _distanceExit = tmax;
        return tmin <= tmax;
    }

    // ********************************************************************* //
    bool intersects( const Ray& _ray, const Triangle& _triangle )
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
        return barycentricCoord1 + barycentricCoord2 <= 1.0f;
    }

    // ********************************************************************* //
    bool intersects( const Ray& _ray, const Triangle& _triangle, float& _distance )
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
        if(barycentricCoord1 < -EPSILON || barycentricCoord1 != barycentricCoord1) return false;
        float barycentricCoord2 = dot( d, e0 ) / dist2A;
        if(barycentricCoord2 < -EPSILON || barycentricCoord2 != barycentricCoord2) return false;
        if(barycentricCoord1 + barycentricCoord2 > 1.0f) return false;
        
        // Projection to plane. The 2A from normal is canceled out
        _distance = dot( normal, o ) / dist2A;
        return _distance >= 0.0f;
    }

    // ********************************************************************* //
    bool intersects( const Ray& _ray, const Triangle& _triangle, float& _distance, Vec3& _barycentric )
    {
        // Möller and Trumbore
        Vec3 e0 = _triangle.v1 - _triangle.v0;
        Vec3 e1 = _triangle.v0 - _triangle.v2;
        // Normal scaled by 2A
        Vec3 normal = cross( e1, e0 );

        float dist2A = dot( normal, _ray.direction );
        Vec3 o = (_triangle.v0 - _ray.origin);
        Vec3 d = cross( _ray.direction, o );

        _barycentric.y = dot( d, e1 ) / dist2A;
        if(_barycentric.y < -EPSILON) return false;
        _barycentric.z = dot( d, e0 ) / dist2A;
        if(_barycentric.z < -EPSILON) return false;
        _barycentric.x = 1.0f - (_barycentric.y + _barycentric.z);
        // Do one check on NaN - if any other coordinate is NaN x will be NaN too
        if(_barycentric.x < -EPSILON || _barycentric.x != _barycentric.x) return false;
        
        // Projection to plane. The 2A from normal is canceled out
        _distance = dot( normal, o ) / dist2A;
        return _distance >= 0.0f;
    }

    // ********************************************************************* //
    bool intersects( const Sphere& _sphere, const Plane& _plane )
    {
        return abs(dot(_plane.n, _sphere.center) + _plane.d) <= _sphere.radius;
    }

    // ********************************************************************* //
    bool intersects( const Sphere& _sphere, const Triangle& _triangle )
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

    // ********************************************************************* //
    bool intersects( const Sphere& _sphere, const Capsule& _capsule )
    {
        return distance(_sphere.center, _capsule.seg) <= _capsule.radius + _sphere.radius;
    }

    // ********************************************************************* //
    bool intersects( const Vec3& _point, const Capsule& _capsule )
    {
        return distance(_point, _capsule.seg) <= _capsule.radius;
    }

    // ********************************************************************* //
    bool intersects( const Vec3& _point, const FastFrustum& _frustum )
    {
        if(distance(_point, _frustum.nf) > 0.0f) return false;
        if(distance(_point, _frustum.l) < 0.0f) return false;
        if(distance(_point, _frustum.r) < 0.0f) return false;
        if(distance(_point, _frustum.b) < 0.0f) return false;
        if(distance(_point, _frustum.t) < 0.0f) return false;
        return true;
    }

    // ************************************************************************* //
    bool intersects( const Vec3& _point, const Box& _box )
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
    bool intersects( const Vec3& _point, const OBox& _obox )
    {
        Vec3 alignedPoint = _point - _obox.center;
        alignedPoint = transform( alignedPoint, _obox.orientation );
        if(alignedPoint.x < - _obox.sides.x / 2.0f) return false;
        if(alignedPoint.y < - _obox.sides.y / 2.0f) return false;
        if(alignedPoint.z < - _obox.sides.z / 2.0f) return false;
        if(alignedPoint.x >   _obox.sides.x / 2.0f) return false;
        if(alignedPoint.y >   _obox.sides.y / 2.0f) return false;
        return alignedPoint.z <= _obox.sides.z / 2.0f;
    }

    // ************************************************************************* //
    bool intersects( const Ray& _ray, const Sphere& _sphere )
    {
        // Go towards closest point and compare its distance to the radius
        Vec3 o = _ray.origin - _sphere.center;
        float odotd = dot(o, _ray.direction);
        o += _ray.direction * -odotd;
        return lensq(o) <= _sphere.radius * _sphere.radius;
    }

    bool intersects( const Ray& _ray, const Sphere& _sphere, float& _distance )
    {
        // Go towards closest point and compare its distance to the radius
        Vec3 o = _ray.origin - _sphere.center;
        float le = len(o) - 1.0f;
        float odotd = -dot(o, _ray.direction);
        o += _ray.direction * odotd;
        float distSq = lensq(o);
        float rSq = _sphere.radius * _sphere.radius;
        if(distSq > rSq) return false;
        // Compute the correct distance via phytagoras
        _distance = max(0.0f, odotd - sqrt(rSq - distSq));
        return true;
    }

    // ************************************************************************* //
    bool intersects( const Vec3& _point, const DOP& _dop )
    {
        float p = -dot(_point, _dop.n);
        if( p > _dop.d0 ) return false;
        return p > _dop.d1;
    }

    // ********************************************************************* //
    bool intersects( const Sphere& _sphere, const FastFrustum& _frustum )
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

    // ********************************************************************* //
    inline bool separates( const Vec3& _dir, const Vec3* _box, const Vec3* _fru )
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

    bool intersects( const Box& _box, const FastFrustum& _frustum )
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
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(1.0f, 0.0f, 0.0f), _frustum.vertices[2] - _frustum.vertices[0]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(1.0f, 0.0f, 0.0f), _frustum.vertices[4] - _frustum.vertices[0]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(1.0f, 0.0f, 0.0f), _frustum.vertices[5] - _frustum.vertices[1]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(1.0f, 0.0f, 0.0f), _frustum.vertices[6] - _frustum.vertices[2]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(1.0f, 0.0f, 0.0f), _frustum.vertices[7] - _frustum.vertices[3]);
        if(separates(n, boxv, _frustum.vertices)) return false;

        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[1] - _frustum.vertices[0]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[2] - _frustum.vertices[0]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[4] - _frustum.vertices[0]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[5] - _frustum.vertices[1]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[6] - _frustum.vertices[2]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 1.0f, 0.0f), _frustum.vertices[7] - _frustum.vertices[3]);
        if(separates(n, boxv, _frustum.vertices)) return false;

        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[1] - _frustum.vertices[0]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[2] - _frustum.vertices[0]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[4] - _frustum.vertices[0]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[5] - _frustum.vertices[1]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[6] - _frustum.vertices[2]);
        if(separates(n, boxv, _frustum.vertices)) return false;
        n = cross(Vec3(0.0f, 0.0f, 1.0f), _frustum.vertices[7] - _frustum.vertices[3]);
        if(separates(n, boxv, _frustum.vertices)) return false;

        return true;
    }
}