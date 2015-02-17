#include "ei/3dtypes.hpp"

namespace ei {

    // ********************************************************************* //
    FastFrustum::FastFrustum(const Frustum& _frustum)
    {
        // Initialization of planes is difficult in the list, so the const-cast
        // only defers the initialization a bit.

        Vec3 far = _frustum.f * _frustum.direction + _frustum.apex;
        const_cast<DOP&>(nf) = DOP(_frustum.direction, _frustum.n * _frustum.direction + _frustum.apex, far);
        // Get third axis and central off point
        Vec3 right = cross(_frustum.up, _frustum.direction);
        // Use two vectors in the planes to derive the normal.
        Vec3 onPlane = far + _frustum.l*right;
        const_cast<Plane&>(l) = Plane(normalize(cross(_frustum.up, onPlane)), onPlane);
        onPlane = far + _frustum.r*right;
        const_cast<Plane&>(r) = Plane(normalize(cross(onPlane, _frustum.up)), onPlane);
        onPlane = far + _frustum.b*_frustum.up;
        const_cast<Plane&>(b) = Plane(normalize(cross(onPlane, right)), onPlane);
        onPlane = far + _frustum.t*_frustum.up;
        const_cast<Plane&>(t) = Plane(normalize(cross(right, onPlane)), onPlane);
    }

}