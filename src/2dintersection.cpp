#include "ei/2dintersection.hpp"

namespace ei {

    // ********************************************************************* //
    bool intersects( const Vec2& _point, const Disc2D& _disc )
    {
        return lensq(_point - _disc.center) <= sq(_disc.radius);
    }

    // ********************************************************************* //
    bool intersects( const Disc2D& _disc0, const Disc2D& _disc1 )
    {
        return lensq(_disc0.center - _disc1.center)
            <= sq(_disc0.radius + _disc1.radius);
    }

    bool intersects( const Disc2D& _disc0, const Disc2D& _disc1, Vec2& _outInfo )
    {
        // Compute intermediate results to determine _outInfo later.
        Vec2 dir = _disc1.center - _disc0.center;
        float d = len(dir);
        float overlap = _disc0.radius + _disc1.radius - d;
        if( overlap >= 0.0f )
        {
            if( d == 0.0f ) { _outInfo = _disc0.center; return true; }
            // Compute the half size of the overlap.
            overlap *= 0.5f;
            if( _disc0.radius < _disc1.radius )
            {
                // When one is inside the other one clamp.
                overlap = min( overlap, _disc0.radius );
                // Go the overlap distance relative to the boundary of disc0 back
                // into the direction of the center.
                _outInfo = _disc0.center + dir * ((_disc0.radius - overlap) / d);
            } else {
                // The same with again if circle1 is the smaller one.
                overlap = min( overlap, _disc1.radius );
                _outInfo = _disc1.center - dir * ((_disc1.radius - overlap) / d);
            }
            return true;
        } else return false;
    }

    // ********************************************************************* //
    bool intersects( const Rect2D& _rect0, const Rect2D& _rect1 )
    {
        // There must be an intersection if the sum of side length is larger
        // than that of the bounding rectangle.
        return all((max(_rect0.max, _rect1.max) - min(_rect0.min, _rect1.min))
            <= ((_rect0.max - _rect0.min) + (_rect1.max - _rect1.min)));
    }

    // ********************************************************************* //
    bool intersects( const Disc2D& _disc, const Rect2D& _rect )
    {
        return (_disc.center.x + _disc.radius >= _rect.min.x)
            && (_disc.center.x - _disc.radius <= _rect.max.x)
            && (_disc.center.y + _disc.radius >= _rect.min.y)
            && (_disc.center.y - _disc.radius <= _rect.max.y);
    }
}
