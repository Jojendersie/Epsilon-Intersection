
// ************************************************************************* //
Quaternion::Quaternion( const Vec3& _axis, float _angle )
{
    eiAssert( approx(lensq(_axis), 1.0f), "Expected a normalized axis vector!" );
    _angle *= 0.5f;
    float sinA = sin(_angle);
    r = cos(_angle);
    i = sinA * _axis.x;
    j = sinA * _axis.y;
    k = sinA * _axis.z;
}

// ************************************************************************* //
Quaternion::Quaternion( float _x, float _y, float _z )
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

    double cYcZ = cosY * cosZ;
    double sYcZ = sinY * cosZ;
    double cYsZ = cosY * sinZ;
    double sYsZ = sinY * sinZ;

    i = float(sinX * cYcZ - cosX * sYsZ);
    j = float(cosX * sYcZ + sinX * cYsZ);
    k = float(cosX * cYsZ - sinX * sYcZ);
    r = float(cosX * cYcZ + sinX * sYsZ);

    return normalize();
}


// ************************************************************************* //
Quaternion::Quaternion( const Mat3x3& _m )
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