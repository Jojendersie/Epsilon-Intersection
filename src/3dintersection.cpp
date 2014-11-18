#include "ei/3dintersection.hpp"

namespace ei {

    // ********************************************************************* //
    bool intersects( const Vec3& _point, const Ellipsoid& _ellipsoid )
    {
        // Use ellipsoid equation.
        return lensq((_point - _ellipsoid.center) / _ellipsoid.radii) <= 1.0f;
    }
}