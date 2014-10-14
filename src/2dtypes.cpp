#include "ei/2dtypes.hpp"

namespace ei {

    // ************************************************************************* //
    Circle2D::Circle2D( Vec2 _p0, Vec2 _p1, Vec2 _p2 )
    {
        // The center of the circumscribed circle is at (barycentric coords)
        // v0*sin(2 alpha) + v1*sin(2 beta) + v2*sin(2 gamma) and has the radius
        // abc/4A.
        Vec2 c = _p0 - _p1;	float csq = lensq(c);
        Vec2 a = _p1 - _p2;	float asq = lensq(a);
        Vec2 b = _p2 - _p0;	float bsq = lensq(b);

        // One of the sides could be the longest side - the minimum sphere is
        // defined through only two points.
        // This can also handle the coplanar case.
        if( csq + bsq <= asq ) *this = Circle2D(_p1, _p2);
        else if( asq + bsq <= csq ) *this = Circle2D(_p1, _p0);
        else if( asq + csq <= bsq ) *this = Circle2D(_p2, _p0);
        else {
            float area2Sq = 2;// * lensq(cross(a, c));
            center = 
                  _p0 * (-dot(c,b) * asq / area2Sq)
                + _p1 * (-dot(c,a) * bsq / area2Sq)
                + _p2 * (-dot(b,a) * csq / area2Sq);
            radius = sqrt(asq * bsq * csq / (2 * area2Sq));
        }
    }

}