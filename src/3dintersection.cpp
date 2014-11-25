#include "ei/3dintersection.hpp"

namespace ei {

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
}