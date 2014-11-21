#include "ei/3dintersection.hpp"

namespace ei {

    // ********************************************************************* //
    bool intersects( const Vec3& _point, const Ellipsoid& _ellipsoid )
    {
        // Use ellipsoid equation.
        return lensq(Vec<double,3>((_point - _ellipsoid.center) / _ellipsoid.radii)) <= 1.0;
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
        /*float r = max(_ellipsoid.radii);
        float t = dot(o, _ray.direction);
        t = max(0.0f, - t - r);
        o += _ray.direction * t;//*/

        /*Vec3 scale = max(Vec3(0.0f), 1.0f - _ellipsoid.radii) * _ray.direction;
        scale = Vec3(scale.y + scale.z, scale.x + scale.z, scale.x + scale.y);
        o += scale * _ray.direction;*/

        // Scale system
        o /= _ellipsoid.radii;
        float odoto = dot(o, o);
        // Origin inside?
        if(odoto <= 1.0f) return true;

        // Test if quadratic equation for the hit point has a solution
        Vec3 d = _ray.direction / _ellipsoid.radii;//max(_ellipsoid.radii, Vec3(1.0f));
        float odotd = dot(o, d);
        float ddotd = dot(d, d);
        float phalf = odotd / ddotd;
        float q = odoto / ddotd;
        return phalf * phalf - q >= 0.0f;
    }
}