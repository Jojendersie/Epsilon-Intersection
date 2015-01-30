#include "ei/3dintersection.hpp"

namespace ei {

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
        t = max(0.0f, - t - r);
        if(t > 0.0f)
            o += _ray.direction * t;//*/

        /*Vec3 scale = max(Vec3(0.0f), 1.0f - _ellipsoid.radii) * _ray.direction;
        scale = Vec3(scale.y + scale.z, scale.x + scale.z, scale.x + scale.y);
        o += scale * _ray.direction;*/

        // Scale system
        o /= _ellipsoid.radii;
        float odoto = dot(o, o);

        // Test if quadratic equation for the hit point has a solution
        Vec3 d = _ray.direction / _ellipsoid.radii;
        float odotd = dot(o, d);
        float ddotd = dot(d, d);
        float phalf = odotd / ddotd;
        float q = (odoto - 1.0f) / ddotd;
        return phalf * phalf - q >= 0.0f;
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
}