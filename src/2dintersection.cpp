#include "ei/2dintersection.hpp"

namespace ei {

    // ********************************************************************* //
    bool intersects( const Disc2D& _circle0, const Disc2D& _circle1 )
    {
        return lensq(_circle0.center - _circle1.center)
            <= sq(_circle0.radius + _circle1.radius);
    }

    bool intersects( const Disc2D& _circle0, const Disc2D& _circle1, Vec2& _outInfo )
    {
        // Compute intermediate results to determine _outInfo later.
        Vec2 dir = _circle1.center - _circle0.center;
        float d = len(dir);
        float overlap = _circle0.radius + _circle1.radius - d;
        if( overlap >= 0.0f )
        {
            // Compute the half size of the overlap.
            overlap *= 0.5f;
            if( _circle0.radius < _circle1.radius )
            {
                // When one is inside the other one clamp.
                overlap = min( overlap, _circle0.radius );
                // Go the overlap distance relative to the boundary of circle0 back
                // into the direction of the center.
                _outInfo = _circle0.center + dir * ((_circle0.radius - overlap) / d);
            } else {
                // The same with again if circle1 is the smaller one.
                overlap = min( overlap, _circle1.radius );
                _outInfo = _circle1.center - dir * ((_circle1.radius - overlap) / d);
            }
            return true;
        } else return false;
    }
}
