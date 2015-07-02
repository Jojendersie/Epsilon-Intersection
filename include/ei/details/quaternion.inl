
// ************************************************************************* //
template<typename T>
TQuaternion<T>::TQuaternion( const Vec<T,3>& _axis, T _angle )
{
    eiAssert( approx(lensq(_axis), 1.0f), "Expected a normalized axis vector!" );
    _angle *= 0.5f;
    T sinA = sin(_angle);
    r = cos(_angle);
    // Assert normalization condition
    if( r < static_cast<T>(0) ) {r = -r; sinA = -sinA;}
    i = sinA * _axis.x;
    j = sinA * _axis.y;
    k = sinA * _axis.z;
}

// ************************************************************************* //
template<typename T>
inline TQuaternion<T>::TQuaternion( T _x, T _y, T _z )
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

    i = T(sinX * cZcY - cosX * sZsY);
    j = T(cosX * cZsY + sinX * sZcY);
    k = T(cosX * sZcY - sinX * cZsY);
    r = T(cosX * cZcY + sinX * sZsY);

    // Assert normalization condition
    if( r < static_cast<T>(0) )
    {
        r = -r;
        i = -i;
        j = -j;
        k = -k;
    }

    //*this = normalize(*this);
}


// ************************************************************************* //
template<typename T>
TQuaternion<T>::TQuaternion( const Matrix<T,3,3>& _m )
{
    // Ignore de-orthogonalization

    // Build TQuaternion<T> from rotation matrix
    // Src: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    T trace = _m.m00 + _m.m11 + _m.m22;
    if( trace > 0 )
    {
        float s = T(0.5) / sqrt( trace + T(1) );
        i = ( _m.m21 - _m.m12 ) * s;
        j = ( _m.m02 - _m.m20 ) * s;
        k = ( _m.m10 - _m.m01 ) * s;
        r = T(0.25) / s;
    } else {
        if( _m.m00 > _m.m11 && _m.m00 > _m.m22 )
        {
            float s = T(2) * sqrt( T(1) + _m.m00 - _m.m11 - _m.m22 );
            i = T(0.25) * s;
            j = ( _m.m01 + _m.m10 ) / s;
            k = ( _m.m02 + _m.m20 ) / s;
            r = ( _m.m21 - _m.m12 ) / s;
        } else if( _m.m11 > _m.m22 )
        {
            float s = T(2) * sqrt( T(1) + _m.m11 - _m.m00 - _m.m22 );
            i = ( _m.m01 + _m.m10 ) / s;
            j = T(0.25) * s;
            k = ( _m.m12 + _m.m21 ) / s;
            r = ( _m.m02 - _m.m20 ) / s;
        } else {
            float s = T(2) * sqrt( T(1) + _m.m22 - _m.m00 - _m.m11 );
            i = ( _m.m02 + _m.m20 ) / s;
            j = ( _m.m12 + _m.m21 ) / s;
            k = T(0.25) * s;
            r = ( _m.m10 - _m.m01 ) / s;
        }
    }//*/

    /*r = sqrt( max( T(0), T(1) + _m.m00 + _m.m11 + _m.m22 ) ) * T(0.5);
    i = sqrt( max( T(0), T(1) + _m.m00 - _m.m11 - _m.m22 ) ) * T(0.5) * sgn(_m.m21 - _m.m12);
    j = sqrt( max( T(0), T(1) - _m.m00 + _m.m11 - _m.m22 ) ) * T(0.5) * sgn(_m.m02 - _m.m20);
    k = sqrt( max( T(0), T(1) - _m.m00 - _m.m11 + _m.m22 ) ) * T(0.5) * sgn(_m.m10 - _m.m01);//*/

    *this = normalize(*this);

    // Assert additional normalization condition
    if( r < static_cast<T>(0) ) {i = -i; j = -j; k = -k; r = -r;}
}

template<typename T>
TQuaternion<T>::TQuaternion( const Vec<T,3>& _xAxis, const Vec<T,3>& _yAxis, const Vec<T,3>& _zAxis ) :
    TQuaternion<T>(axis(_xAxis, _yAxis, _zAxis))
{
}

template<typename T>
TQuaternion<T>::TQuaternion( T _i, T _j, T _k, T _r ) :
    i(_i), j(_j), k(_k), r(_r)
{
    // Assert additional normalization condition
    if( r < static_cast<T>(0) )
    {
        i = -i;
        j = -j;
        k = -k;
        r = -r;
    }
}

// ************************************************************************* //
template<typename T>
TQuaternion<T>::TQuaternion( const Vec<T,3>& _from, const Vec<T,3>& _to )
{
    Vec<T,3> from = normalize(_from);
    Vec<T,3> to = normalize(_to);
    // half angle trick from http://physicsforgames.blogspot.de/2010/03/quaternion-tricks.html
    Vec<T,3> half = normalize(from + to);

    // cos(theta) = dot product since both vectors are normalized
    r = dot(from, half);
    // Axis from cross product -> already multiplied with sin(theta)
    i = from.y*half.z - from.z*half.y;
    j = from.z*half.x - from.x*half.z;
    k = from.x*half.y - from.y*half.x;
}

// ************************************************************************* //
inline const TQuaternion<float>& qidentity()
{
    return details::QUATERNION_IDENTITY;
}

inline const TQuaternion<double>& qidentityD()
{
    return details::QUATERNIOND_IDENTITY;
}

// ************************************************************************* //
template<typename T>
bool TQuaternion<T>::operator== (const TQuaternion<T>& _q1) const
{
    return r==_q1.r && i==_q1.i && j==_q1.j && k==_q1.k;
}

// ************************************************************************* //
template<typename T>
bool TQuaternion<T>::operator!= (const TQuaternion<T>& _q1) const
{
    return r!=_q1.r || i!=_q1.i || j!=_q1.j || k!=_q1.k;
}

// ************************************************************************* //
template<typename T>
TQuaternion<T>& TQuaternion<T>::operator*= (const TQuaternion<T>& _q1)
{
    T nr = r*_q1.r - i*_q1.i - j*_q1.j - k*_q1.k;
    T ni = r*_q1.i + i*_q1.r + j*_q1.k - k*_q1.j;
    T nj = r*_q1.j + j*_q1.r + k*_q1.i - i*_q1.k;
       k = r*_q1.k + k*_q1.r + i*_q1.j - j*_q1.i;
    r = nr;
    i = ni;
    j = nj;
    return *this;
}

// ************************************************************************* //
template<typename T>
TQuaternion<T>& TQuaternion<T>::operator*= (T _s)
{
    i*=_s; j*=_s; k*=_s; r*=_s;
    return *this;
}

// ************************************************************************* //
template<typename T>
TQuaternion<T>& TQuaternion<T>::operator/= (const TQuaternion<T>& _q1)
{
    T nr =   r*_q1.r + i*_q1.i + j*_q1.j + k*_q1.k;
    T ni = - r*_q1.i + i*_q1.r - j*_q1.k + k*_q1.j;
    T nj = - r*_q1.j + j*_q1.r - k*_q1.i + i*_q1.k;
       k = - r*_q1.k + k*_q1.r - i*_q1.j + j*_q1.i;
    r = nr;
    i = ni;
    j = nj;
    return *this;
}

// ************************************************************************* //
template<typename T>
TQuaternion<T>& TQuaternion<T>::operator/= (T _s)
{
    i/=_s; j/=_s; k/=_s; r/=_s;
    return *this;
}

// ************************************************************************* //
template<typename T>
TQuaternion<T>& TQuaternion<T>::operator+= (const TQuaternion<T>& _q1)
{
    i+=_q1.i; j+=_q1.j; k+=_q1.k; r+=_q1.r;
    return *this;
}

// ************************************************************************* //
template<typename T>
TQuaternion<T>& TQuaternion<T>::operator-= (const TQuaternion<T>& _q1)
{
    i-=_q1.i; j-=_q1.j; k-=_q1.k; r-=_q1.r;
    return *this;
}



// ************************************************************************* //
template<typename T>
TQuaternion<T> conjugate(const TQuaternion<T>& _q)
{
    return TQuaternion<T>(-_q.i, -_q.j, -_q.k, _q.r);
}

template<typename T>
Vec<T,3> axis(const TQuaternion<T>& _q)
{
    return Vec<T,3>(-_q.i, -_q.j, -_q.k) / max(EPSILON, sqrt(T(1)-_q.r*_q.r));
}

template<typename T>
Vec<T,3> xaxis(const TQuaternion<T>& _q)
{
    return Vec<T,3>( T(1)-T(2)*(_q.j*_q.j+_q.k*_q.k), T(2)*(_q.i*_q.j-_q.k*_q.r), T(2)*(_q.i*_q.k+_q.j*_q.r) );
}
template<typename T>
Vec<T,3> yaxis(const TQuaternion<T>& _q)
{
    return Vec<T,3>( T(2)*(_q.i*_q.j+_q.k*_q.r), T(1)-T(2)*(_q.i*_q.i+_q.k*_q.k), T(2)*(_q.j*_q.k-_q.i*_q.r) );
}
template<typename T>
Vec<T,3> zaxis(const TQuaternion<T>& _q)
{
    return Vec<T,3>( T(2)*(_q.i*_q.k-_q.j*_q.r), T(2)*(_q.j*_q.k+_q.i*_q.r), T(1)-T(2)*(_q.i*_q.i+_q.j*_q.j) );
}

template<typename T>
T angle(const TQuaternion<T>& _q)
{
    return acos(_q.r) * T(2);
}

// ************************************************************************* //
template<typename T>
Vec<T,3> angles(const TQuaternion<T>& _q)
{
    Vec<T,3> angles;
    // Derivation from http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToEuler/index.htm
    // but changed angles because of else convention
    const double m20half = _q.j * _q.r - _q.i * _q.k;

    if(approx(m20half, 0.5))
    {
        angles.x = 0.0f;
        angles.y = PI/2.0f;
        angles.z = static_cast<T>(-2.0 * atan2(_q.i, _q.r));
    }
    else if(approx(m20half, -0.5))
    {
        angles.x = 0.0f;
        angles.y = -PI/2.0f;
        angles.z = static_cast<T>(2.0 * atan2(_q.i, _q.r));
    }
    else
    {
        const double sqr = _q.r * _q.r;
        const double sqi = _q.i * _q.i;
        const double sqj = _q.j * _q.j;
        const double sqk = _q.k * _q.k;
        angles.x = static_cast<T>(atan2(2.0 * (_q.j * _q.k + _q.i * _q.r), -sqi - sqj + sqk + sqr));
        angles.y = static_cast<T>(asin( clamp(m20half * 2.0, -1.0, 1.0) ));
        angles.z = static_cast<T>(atan2(2.0 * (_q.i * _q.j + _q.k * _q.r),  sqi - sqj - sqk + sqr));
    }
    return angles;
}

// ************************************************************************* //
template<typename T>
bool approx(const TQuaternion<T>& _q0,
            const TQuaternion<T>& _q1,
            T _epsilon)
{
    return abs(_q0.r - _q1.r) <= _epsilon
        && abs(_q0.i - _q1.i) <= _epsilon
        && abs(_q0.j - _q1.j) <= _epsilon
        && abs(_q0.k - _q1.k) <= _epsilon;
}

// ************************************************************************* //
template<typename T>
TQuaternion<T> slerp(const TQuaternion<T>& _q0, const TQuaternion<T>& _q1, T _t)
{
    // http://en.wikipedia.org/wiki/Slerp
    T theta = acos( clamp(dot(_q0,_q1), T(-1), T(1)) );
    T so = sin( theta );
    if(approx(so, T(0)))
    {
        // Converges towards linear interpolation for small so
        return TQuaternion<T>(_q0.i + (_q1.i - _q0.i) * _t,
                              _q0.j + (_q1.j - _q0.j) * _t,
                              _q0.k + (_q1.k - _q0.k) * _t,
                              _q0.r + (_q1.r - _q0.r) * _t);
    }
    T f0 = sin( theta * (1.0f-_t) ) / so;
    T f1 = sin( theta * _t ) / so;
    return TQuaternion<T>(_q0.i * f0 + _q1.i * f1,
                          _q0.j * f0 + _q1.j * f1,
                          _q0.k * f0 + _q1.k * f1,
                          _q0.r * f0 + _q1.r * f1);
}

// ************************************************************************* //
template<typename T>
Vec<T,3> transform( const Vec<T,3>& _v, const TQuaternion<T>& _q )
{
    // http://physicsforgames.blogspot.de/2010/03/quaternion-tricks.html
    float x1 = _q.j*_v.z - _q.k*_v.y;
    float y1 = _q.k*_v.x - _q.i*_v.z;
    float z1 = _q.i*_v.y - _q.j*_v.x;

    return Vec<T,3>(
        _v.x + 2.0f * (_q.r*x1 + _q.j*z1 - _q.k*y1),
        _v.y + 2.0f * (_q.r*y1 + _q.k*x1 - _q.i*z1),
        _v.z + 2.0f * (_q.r*z1 + _q.i*y1 - _q.j*x1)
    );

    // q v q-1 with v=(0, _v.x, _v.y, _v.z) expanded with Maxima
/*    return Vec<T,3>(
        _v.x*_q.r*_q.r-2*_v.y*_q.k*_q.r+2*_q.j*_v.z*_q.r-_v.x*_q.k*_q.k+2*_q.i*_v.z*_q.k-_v.x*_q.j*_q.j+2*_q.i*_v.y*_q.j+_v.x*_q.i*_q.i,
        _v.y*_q.r*_q.r+2*_v.x*_q.k*_q.r-2*_q.i*_v.z*_q.r-_v.y*_q.k*_q.k+2*_q.j*_v.z*_q.k+_v.y*_q.j*_q.j+2*_v.x*_q.i*_q.j-_q.i*_q.i*_v.y,
        _v.z*_q.r*_q.r-2*_v.x*_q.j*_q.r+2*_q.i*_v.y*_q.r+_v.z*_q.k*_q.k+2*_v.y*_q.j*_q.k+2*_v.x*_q.i*_q.k-_q.j*_q.j*_v.z-_q.i*_q.i*_v.z
    );*/
}

template<typename T>
RVec<T,3> transform( const RVec<T,3>& _v, const TQuaternion<T>& _q )
{
    float x1 = _q.j*_v.z - _q.k*_v.y;
    float y1 = _q.k*_v.x - _q.i*_v.z;
    float z1 = _q.i*_v.y - _q.j*_v.x;

    return Vec<T,3>(
        _v.x + 2.0f * (_q.r*x1 + _q.j*z1 - _q.k*y1),
        _v.y + 2.0f * (_q.r*y1 + _q.k*x1 - _q.i*z1),
        _v.z + 2.0f * (_q.r*z1 + _q.i*y1 - _q.j*x1)
    );

/*    return RVec<T,3>(
        _v.x*_q.r*_q.r-2*_v.y*_q.k*_q.r+2*_q.j*_v.z*_q.r-_v.x*_q.k*_q.k+2*_q.i*_v.z*_q.k-_v.x*_q.j*_q.j+2*_q.i*_v.y*_q.j+_v.x*_q.i*_q.i,
        _v.y*_q.r*_q.r+2*_v.x*_q.k*_q.r-2*_q.i*_v.z*_q.r-_v.y*_q.k*_q.k+2*_q.j*_v.z*_q.k+_v.y*_q.j*_q.j+2*_v.x*_q.i*_q.j-_q.i*_q.i*_v.y,
        _v.z*_q.r*_q.r-2*_v.x*_q.j*_q.r+2*_q.i*_v.y*_q.r+_v.z*_q.k*_q.k+2*_v.y*_q.j*_q.k+2*_v.x*_q.i*_q.k-_q.j*_q.j*_v.z-_q.i*_q.i*_v.z
    );*/
}
