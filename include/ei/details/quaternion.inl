
// ************************************************************************* //
inline Quaternion::Quaternion( const Vec3& _axis, float _angle )
{
    eiAssert( approx(lensq(_axis), 1.0f), "Expected a normalized axis vector!" );
    _angle *= 0.5f;
    float sinA = sin(_angle);
    r = cos(_angle);
    // Assert normalization condition
    if( r < 0.0f ) {r = -r; sinA = -sinA;}
    i = sinA * _axis.x;
    j = sinA * _axis.y;
    k = sinA * _axis.z;
}

// ************************************************************************* //
inline Quaternion::Quaternion( float _x, float _y, float _z )
{
    double halfAngle;

    halfAngle = _x * 0.5;
    double sinX = sin(halfAngle);
    double cosX = cos(halfAngle);

    halfAngle = _y * 0.5;
    double sinY = sin(halfAngle);
    double cosY = cos(halfAngle);

    halfAngle = _z * 0.5;
    double sinZ = sin(halfAngle);
    double cosZ = cos(halfAngle);

    double cZcY = cosZ * cosY;
    double cZsY = cosZ * sinY;
    double sZcY = sinZ * cosY;
    double sZsY = sinZ * sinY;

    i = float(sinX * cZcY - cosX * sZsY);
    j = float(cosX * cZsY + sinX * sZcY);
    k = float(cosX * sZcY - sinX * cZsY);
    r = float(cosX * cZcY + sinX * sZsY);

    // Assert normalization condition
    if( r < 0.0f )
    {
        r = -r;
        i = -i;
        j = -j;
        k = -k;
    }

    //*this = normalize(*this);
}


// ************************************************************************* //
inline Quaternion::Quaternion( const Mat3x3& _m )
{
    // Ignore de-orthogonalization

    // Build Quaternion from rotation matrix
    // Src: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    /*float trace = M.m11 + M.m22 + M.m33;
    if( trace > 0 )
    {
        float s = 0.5f / sqrt( trace + 1.0f );
        i = ( M.m32 - M.m23 ) * s;
        j = ( M.m13 - M.m31 ) * s;
        k = ( M.m21 - M.m12 ) * s;
        r = 0.25f / s;
    } else {
        if( M.m11 > M.m22 && M.m11 > M.m33 )
        {
            float s = 2.0f * sqrtf( 1.0f + M.m11 - M.m22 - M.m33 );
            i = 0.25f * s;
            j = ( M.m12 + M.m21 ) / s;
            k = ( M.m13 + M.m31 ) / s;
            r = ( M.m32 - M.m23 ) / s;
        } else if( M.m22 > M.m33 )
        {
            float s = 2.0f * sqrtf( 1.0f + M.m22 - M.m11 - M.m33 );
            i = ( M.m12 + M.m21 ) / s;
            j = 0.25f * s;
            k = ( M.m23 + M.m32 ) / s;
            r = ( M.m13 - M.m31 ) / s;
        } else {
            float s = 2.0f * sqrtf( 1.0f + M.m33 - M.m11 - M.m22 );
            i = ( M.m13 + M.m31 ) / s;
            j = ( M.m23 + M.m32 ) / s;
            k = 0.25f * s;
            r = ( M.m21 - M.m12 ) / s;
        }
    }*/

    r = sqrt( max( 0.0f, 1.0f + _m.m00 + _m.m11 + _m.m22 ) ) * 0.5f;
    i = sqrt( max( 0.0f, 1.0f + _m.m00 - _m.m11 - _m.m22 ) ) * 0.5f * sign(_m.m21 - _m.m12);
    j = sqrt( max( 0.0f, 1.0f - _m.m00 + _m.m11 - _m.m22 ) ) * 0.5f * sign(_m.m02 - _m.m20);
    k = sqrt( max( 0.0f, 1.0f - _m.m00 - _m.m11 + _m.m22 ) ) * 0.5f * sign(_m.m10 - _m.m01);
}

inline Quaternion::Quaternion( float _r, float _i, float _j, float _k ) :
    r(_r), i(_i), j(_j), k(_k)
{
}

// ************************************************************************* //
inline Quaternion::Quaternion( const Vec3& _from, const Vec3& _to )
{
    // Get lengths for normalization
    float lf = len(_from);
    float lt = len(_to);
    Vec3 axis = cross(_from, _to);
    // Compute sin(alpha)^2 from cross product lf * lt * sin(alpha) and normalize
    float sa = len(axis);
    axis /= sa;
    sa /= lf * lt;
    // sin(alpha) = sin(2*theta)
    float theta = asin(sa) * 0.5f;
    float st = sin(theta);
    r = sqrt((1.0f - st) * (1.0f + st));      // cos(theta)
    i = st * axis.x;
    j = st * axis.y;
    k = st * axis.z;
}

// ************************************************************************* //
inline bool Quaternion::operator== (const Quaternion& _q1) const
{
    return r==_q1.r && i==_q1.i && j==_q1.j && k==_q1.k;
}

// ************************************************************************* //
inline bool Quaternion::operator!= (const Quaternion& _q1) const
{
    return r!=_q1.r || i!=_q1.i || j!=_q1.j || k!=_q1.k;
}

// ************************************************************************* //
inline Quaternion& Quaternion::operator*= (const Quaternion& _q1)
{
    float nr = r*_q1.r - i*_q1.i - j*_q1.j - k*_q1.k;
    float ni = r*_q1.i + i*_q1.r + j*_q1.k - k*_q1.j;
    float nj = r*_q1.j + j*_q1.r + k*_q1.i - i*_q1.k;
           k = r*_q1.k + k*_q1.r + i*_q1.j - j*_q1.i;
    r = nr;
    i = ni;
    j = nj;
    return *this;
}

// ************************************************************************* //
inline Quaternion& Quaternion::operator*= (float _s)
{
    i*=_s; j*=_s; k*=_s; r*=_s;
    return *this;
}

// ************************************************************************* //
inline Quaternion& Quaternion::operator/= (const Quaternion& _q1)
{
    float nr =   r*_q1.r + i*_q1.i + j*_q1.j + k*_q1.k;
    float ni = - r*_q1.i + i*_q1.r - j*_q1.k + k*_q1.j;
    float nj = - r*_q1.j + j*_q1.r - k*_q1.i + i*_q1.k;
           k = - r*_q1.k + k*_q1.r - i*_q1.j + j*_q1.i;
    r = nr;
    i = ni;
    j = nj;
    return *this;
}

// ************************************************************************* //
inline Quaternion& Quaternion::operator/= (float _s)
{
    i/=_s; j/=_s; k/=_s; r/=_s;
    return *this;
}

// ********************************************************************* //
inline Quaternion& Quaternion::operator+= (const Quaternion& _q1)
{
    i+=_q1.i; j+=_q1.j; k+=_q1.k; r+=_q1.r;
    return *this;
}

// ********************************************************************* //
inline Quaternion& Quaternion::operator-= (const Quaternion& _q1)
{
    i-=_q1.i; j-=_q1.j; k-=_q1.k; r-=_q1.r;
    return *this;
}



// ********************************************************************* //
inline Quaternion conjugate(const Quaternion& _q)
{
    return Quaternion(_q.r, -_q.i, -_q.j, -_q.k);
}

inline Vec3 axis(const Quaternion& _q)
{
    return Vec3(-_q.i, -_q.j, -_q.k) / max(EPSILON, sqrt(1.0f-_q.r*_q.r));
}

inline float angle(const Quaternion& _q)
{
    return acos(_q.r) * 2.0f;
}

// ************************************************************************* //
inline Vec3 angles(const Quaternion& _q)
{
    Vec3 angles;
    // Derivation from http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToEuler/index.htm
    // but changed angles because of else convention
    const double m20half = _q.j * _q.r - _q.i * _q.k;

    if(approx(m20half, 0.5))
    {
        angles.x = 0.0f;
        angles.y = PI/2.0f;
        angles.z = float(-2.0 * atan2(_q.i, _q.r));
    }
    else if(approx(m20half, -0.5))
    {
        angles.x = 0.0f;
        angles.y = -PI/2.0f;
        angles.z = float(2.0 * atan2(_q.i, _q.r));
    }
    else
    {
        const double sqr = _q.r * _q.r;
        const double sqi = _q.i * _q.i;
        const double sqj = _q.j * _q.j;
        const double sqk = _q.k * _q.k;
        angles.x = (float)atan2(2.0 * (_q.j * _q.k + _q.i * _q.r), -sqi - sqj + sqk + sqr);
        angles.y = (float)asin( clamp(m20half * 2.0, -1.0, 1.0) );
        angles.z = (float)atan2(2.0 * (_q.i * _q.j + _q.k * _q.r),  sqi - sqj - sqk + sqr);
    }
    return angles;
}

// ************************************************************************* //
inline bool approx(const Quaternion& _q0,
                   const Quaternion& _q1,
                   float _epsilon)
{
    return abs(_q0.r - _q1.r) <= _epsilon
        && abs(_q0.i - _q1.i) <= _epsilon
        && abs(_q0.j - _q1.j) <= _epsilon
        && abs(_q0.k - _q1.k) <= _epsilon;
}

// ************************************************************************* //
inline Quaternion slerp(const Quaternion& _q0, const Quaternion& _q1, float _t)
{
    // http://en.wikipedia.org/wiki/Slerp
    float theta = acos( clamp(dot(_q0,_q1), -1.0f, 1.0f) );
    float so = sin( theta );
    if(approx(so, 0.0f))
    {
        // Converges towards linear interpolation for small so
        return Quaternion(_q0.r + (_q1.r - _q0.r) * _t,
                          _q0.i + (_q1.i - _q0.i) * _t,
                          _q0.j + (_q1.j - _q0.j) * _t,
                          _q0.k + (_q1.k - _q0.k) * _t);
    }
    float f0 = sin( theta * (1.0f-_t) ) / so;
    float f1 = sin( theta * _t ) / so;
    return Quaternion(_q0.r * f0 + _q1.r * f1,
                      _q0.i * f0 + _q1.i * f1,
                      _q0.j * f0 + _q1.j * f1,
                      _q0.k * f0 + _q1.k * f1);
}

// ************************************************************************* //
inline Vec3 transform( const Vec3& _v, const Quaternion& _q )
{
    // q v q-1 with v=(0, _v.x, _v.y, _v.z) expanded with Maxima
    return Vec3(
        _v.x*_q.r*_q.r-2*_v.y*_q.k*_q.r+2*_q.j*_v.z*_q.r-_v.x*_q.k*_q.k+2*_q.i*_v.z*_q.k-_v.x*_q.j*_q.j+2*_q.i*_v.y*_q.j+_v.x*_q.i*_q.i,
        _v.y*_q.r*_q.r+2*_v.x*_q.k*_q.r-2*_q.i*_v.z*_q.r-_v.y*_q.k*_q.k+2*_q.j*_v.z*_q.k+_v.y*_q.j*_q.j+2*_v.x*_q.i*_q.j-_q.i*_q.i*_v.y,
        _v.z*_q.r*_q.r-2*_v.x*_q.j*_q.r+2*_q.i*_v.y*_q.r+_v.z*_q.k*_q.k+2*_v.y*_q.j*_q.k+2*_v.x*_q.i*_q.k-_q.j*_q.j*_v.z-_q.i*_q.i*_v.z
    );
}

inline Vec3 transform( const RVec3& _v, const Quaternion& _q )
{
    return Vec3(
        _v.x*_q.r*_q.r-2*_v.y*_q.k*_q.r+2*_q.j*_v.z*_q.r-_v.x*_q.k*_q.k+2*_q.i*_v.z*_q.k-_v.x*_q.j*_q.j+2*_q.i*_v.y*_q.j+_v.x*_q.i*_q.i,
        _v.y*_q.r*_q.r+2*_v.x*_q.k*_q.r-2*_q.i*_v.z*_q.r-_v.y*_q.k*_q.k+2*_q.j*_v.z*_q.k+_v.y*_q.j*_q.j+2*_v.x*_q.i*_q.j-_q.i*_q.i*_v.y,
        _v.z*_q.r*_q.r-2*_v.x*_q.j*_q.r+2*_q.i*_v.y*_q.r+_v.z*_q.k*_q.k+2*_v.y*_q.j*_q.k+2*_v.x*_q.i*_q.k-_q.j*_q.j*_v.z-_q.i*_q.i*_v.z
    );
}
