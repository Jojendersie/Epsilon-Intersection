
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

    double cXcY = cosX * cosY;
    double cXsY = cosX * sinY;
    double sXcY = sinX * cosY;
    double sXsY = sinX * sinY;

    i = float(sinZ * cXcY - cosZ * sXsY);
    j = float(cosZ * cXsY + sinZ * sXcY);
    k = float(cosZ * sXcY - sinZ * cXsY);
    r = float(cosZ * cXcY + sinZ * sXsY);

    // Assert normalization condition
    if( r < 0.0f )
    {
        r = -r;
        i = -i;
        j = -j;
        k = -k;
    }
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
inline bool approx(const Quaternion& _q0,
                   const Quaternion& _q1,
                   float _epsilon)
{
    return abs(_q0.r - _q1.r) <= _epsilon
        && abs(_q0.i - _q1.i) <= _epsilon
        && abs(_q0.j - _q1.j) <= _epsilon
        && abs(_q0.k - _q1.k) <= _epsilon;
}