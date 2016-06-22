#include "ei/2dintersection.hpp"

namespace ei {

    // ************************************************************************* //
    float distanceSq(const Vec2& _point, const Segment2D& _line)
    {
        Vec2 u = _line.b - _line.a;
        Vec2 w = _point - _line.a;

        float s = saturate(dot(u, w) / dot(u, u));
        // _line.a + u * s - _point
        return lensq(u * s - w);
    }



    // ********************************************************************* //
    bool intersects( const Vec2& _point, const Capsule2D& _capsule )
    {
        return distanceSq(_point, _capsule.seg) <= _capsule.radius * _capsule.radius;
    }

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
    bool intersects( const Vec2& _point, const Rect2D& _rect )
    {
        return all(_point >= _rect.min) && all(_point <= _rect.max);
    }

    // ********************************************************************* //
    bool intersects( const Disc2D& _disc, const Rect2D& _rect )
    {
        return all( _disc.center + _disc.radius >= _rect.min )
            && all( _disc.center - _disc.radius <= _rect.max );
    }

    // ********************************************************************* //
    bool intersects( const Segment2D& _line0, const Segment2D& _line1 )
    {
        Vec2 dir0 = _line0.b - _line0.a;
        Vec2 dir1 = _line1.b - _line1.a;
        // Solve _line0.p0 + s * dir0 == _line1.p0 + t * dir1
        float den = dot( dir0, Vec2(-dir1.y, dir1.x) );
        if( den == 0.0f ) // Parallel or identical?
            return dot(_line0.a - _line1.a, Vec2(-dir1.y, dir1.x)) == 0.0f;

        // Compute and check first ray parameter
        float s = dot(_line1.a - _line0.a, Vec2(-dir1.y, dir1.x)) / den;
        if( s < 0.0f || s > 1.0f ) return false;

        // Compute and check second ray parameter
        float t = dot(_line1.a - _line0.a, Vec2(-dir0.y, dir0.x)) / den;
        return t >= 0.0f && t <= 1.0f;
    }

    bool intersects( const Segment2D& _line0, const Segment2D& _line1, Vec2& _outInfo )
    {
        Vec2 dir0 = _line0.b - _line0.a;
        Vec2 dir1 = _line1.b - _line1.a;
        // Solve _line0.p0 + s * dir0 == _line1.p0 + t * dir1
        float den = dot( dir0, Vec2(-dir1.y, dir1.x) );
        if( den == 0.0f ) // Parallel or identical?
        {
            if( dot(_line0.a - _line1.a, Vec2(-dir1.y, dir1.x)) != 0.0f )
                return false;
            // Find the center of overlap
            float d0 = dot(dir0, dir0);                   // squared distance from _line0.a to _line0.b
            float d1 = dot(dir0, _line1.a - _line0.a);    // signed distance from _line0.a to _line1.a
            float d2 = dot(dir0, _line1.b - _line0.a);    // signed distance from _line0.a to _line1.b
            // "sort" the 4 relative distance and use the two in the middle
            if( d1 >= 0.0f )
            {
                if( d2 >= 0.0f )
                {
                    if( d0 > d1 && d0 > d2 ) _outInfo = (_line1.b + _line1.a) * 0.5f;
                    else if( d1 > d2 )       _outInfo = (_line1.b + _line0.b) * 0.5f;
                    else                     _outInfo = (_line1.a + _line0.b) * 0.5f;
                } else {
                    if( d1 < d0 )            _outInfo = (_line1.a + _line0.a) * 0.5f;
                    else                     _outInfo = (_line0.b + _line0.a) * 0.5f;
                }
            } else {
                if( d2 < d0 )                _outInfo = (_line1.b + _line0.a) * 0.5f;
                else                         _outInfo = (_line0.b + _line0.a) * 0.5f;
            }
            return true;
        }

        // Compute and check first ray parameter
        float s = dot(_line1.a - _line0.a, Vec2(-dir1.y, dir1.x)) / den;
        if( s < 0.0f || s > 1.0f ) return false;

        // Compute and check second ray parameter
        float t = dot(_line1.a - _line0.a, Vec2(-dir0.y, dir0.x)) / den;
        if( t < 0.0f || t > 1.0f ) return false;

        // Compute output point
        _outInfo = _line0.a + dir0 * s;
        return true;
    }
}
