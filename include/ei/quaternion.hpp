#pragma once

#include "vector.hpp"

namespace ei {

    // ********************************************************************* //
    //                            QUATERNION TYPE                            //
    // ********************************************************************* //
    
    // ********************************************************************* //
    /// \brief 4D complex number equivalent for the representation of rotations
    /// \details The normalized form has len(q) == 1 and r>0 (RHS) / r<0 (LHS).
    ///     The second criteria makes the rotation unique because q and -q both
    ///     represent the same rotation and allows including mirroring.
    ///     This is not the standard way: Usual quaternions cannot handle
    ///     mirroring!
    template<typename T>
    class TQuaternion: public details::NonScalarType
    {
    public:
        /// \brief Construct uninitialized
        constexpr TQuaternion() noexcept = default;

        /// \brief Copy construction
        constexpr TQuaternion( const TQuaternion& _other ) noexcept = default;
        /// \brief Copying assignment
        constexpr TQuaternion& operator = ( const TQuaternion& _rhs ) noexcept = default;

        /// \brief Construct from normalized axis and angle
        constexpr TQuaternion( const Vec<T,3>& _axis, T _angle ) noexcept // TESTED
        {
            eiAssert( approx(lensq(_axis), static_cast<T>(1)), "Expected a normalized axis vector!" );
            _angle *= 0.5f;
            T sinA = sin(_angle);
            r = cos(_angle);
            // Assert normalization condition
            if( r < static_cast<T>(0) ) {r = -r; sinA = -sinA;}
            i = sinA * _axis.x;
            j = sinA * _axis.y;
            k = sinA * _axis.z;
        }

        /// \brief Create from Euler angles
        /// \details The rotations are applied in the order x, y, z:
        ///     rotationZ(_z) * rotationY(_y) * rotationX(_x)
        constexpr TQuaternion( T _x, T _y, T _z ) noexcept // TESTED
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
        constexpr TQuaternion( const Vec<T,3>& _eulerAngles ) noexcept :
            TQuaternion(_eulerAngles.x, _eulerAngles.y, _eulerAngles.z)
        {}

        /// \brief Create from rotation matrix (does a decomposition if the
        ///     matrix contains scaling).
        constexpr TQuaternion( const Matrix<T,3,3>& _matrix ) noexcept : // TESTED
            TQuaternion<T>(transpose(_matrix(0)), transpose(_matrix(1)), transpose(_matrix(2)))
        {}

        /// \brief Create from orthogonal basis vectors.
        constexpr TQuaternion( const Vec<T,3>& _xAxis, const Vec<T,3>& _yAxis, const Vec<T,3>& _zAxis ) noexcept 
        {
            // Check handness
            //eiAssert(dot(cross(_xAxis, _yAxis), _m(2)) > 0.0f, "Quaternions cannot handle reflections. The matrix must be RHS.");
            T handness = dot(cross(_xAxis, _yAxis), _zAxis);
            Vec<T,3> zAxis;
            if(handness < T(0)) zAxis = -_zAxis;
            else zAxis = _zAxis;
            eiAssert(approx(abs(handness), T(1), 1e-4f), "System is not orthonormal!");

            // Build TQuaternion<T> from rotation matrix
            // Src: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
            T trace = _xAxis.x + _yAxis.y + zAxis.z;
            if( trace > 0 )
            {
                float s = T(0.5) / sqrt( trace + T(1) );
                i = (  zAxis.y - _yAxis.z ) * s;
                j = ( _xAxis.z -  zAxis.x ) * s;
                k = ( _yAxis.x - _xAxis.y ) * s;
                r = T(0.25) / s;
            } else {
                if( _xAxis.x > _yAxis.y && _xAxis.x > zAxis.z )
                {
                    float s = T(2) * sqrt( T(1) + _xAxis.x - _yAxis.y - zAxis.z );
                    i = T(0.25) * s;
                    j = ( _xAxis.y + _yAxis.x ) / s;
                    k = ( _xAxis.z +  zAxis.x ) / s;
                    r = (  zAxis.y - _yAxis.z ) / s;
                } else if( _yAxis.y > zAxis.z )
                {
                    float s = T(2) * sqrt( T(1) + _yAxis.y - _xAxis.x - zAxis.z );
                    i = ( _xAxis.y + _yAxis.x ) / s;
                    j = T(0.25) * s;
                    k = ( _yAxis.z +  zAxis.y ) / s;
                    r = ( _xAxis.z -  zAxis.x ) / s;
                } else {
                    float s = T(2) * sqrt( T(1) + zAxis.z - _xAxis.x - _yAxis.y );
                    i = ( _xAxis.z +  zAxis.x ) / s;
                    j = ( _yAxis.z +  zAxis.y ) / s;
                    k = T(0.25) * s;
                    r = ( _yAxis.x - _xAxis.y ) / s;
                }
            }//*/

             /*r = sqrt( max( T(0), T(1) + _m.m00 + _m.m11 + _m.m22 ) ) * T(0.5);
             i = sqrt( max( T(0), T(1) + _m.m00 - _m.m11 - _m.m22 ) ) * T(0.5) * sgn(_m.m21 - _m.m12);
             j = sqrt( max( T(0), T(1) - _m.m00 + _m.m11 - _m.m22 ) ) * T(0.5) * sgn(_m.m02 - _m.m20);
             k = sqrt( max( T(0), T(1) - _m.m00 - _m.m11 + _m.m22 ) ) * T(0.5) * sgn(_m.m10 - _m.m01);//*/

            *this = normalize(*this);

            // Assert additional normalization condition
            if( r * handness < static_cast<T>(0) ) {i = -i; j = -j; k = -k; r = -r;}
        }

        /// \brief Create from TQuaternion coefficients
        constexpr TQuaternion( T _i, T _j, T _k, T _r ) noexcept :
            i(_i), j(_j), k(_k), r(_r)
        {}

        /// \brief Rotate from vector to vector (rotated such that the from
        ///     vector is aligned with the to vector).
        /// \param [in] _from One certain direction vector before rotation.
        /// \param [in] _to Target direction vector. The from direction should
        ///     be aligned with the target after rotation.
        constexpr TQuaternion( const Vec<T,3>& _from, const Vec<T,3>& _to ) noexcept // TESTED
        {
            Vec<T,3> from = normalize(_from);
            Vec<T,3> to = normalize(_to);
            // half angle trick from http://physicsforgames.blogspot.de/2010/03/quaternion-tricks.html
            Vec<T,3> half = normalize(from + to);
            // Opposite vectors or one vector 0.0 -> 180° rotation
            if(half.x != half.x)
            {
                if(approx(ei::abs(from.y), static_cast<T>(1)))
                    half = Vec<T,3>(0, 0, 1);
                else half = normalize(cross(from, Vec<T,3>(0, 1, 0)));
            }

            // cos(theta) = dot product since both vectors are normalized
            r = dot(from, half);
            eiAssert(r >= T(0), "Normalization condition violated!");
            // Axis from cross product -> already multiplied with sin(theta)
            i = from.y*half.z - from.z*half.y;
            j = from.z*half.x - from.x*half.z;
            k = from.x*half.y - from.y*half.x;
        }

        // TODO: lookAt parametrization

        /// \brief Compare component wise, if two quaternions are identical.
        constexpr bool operator == (const TQuaternion& _q1) const noexcept
        {
            return r==_q1.r && i==_q1.i && j==_q1.j && k==_q1.k;
        }
        /// \brief Compare component wise, if two quaternions are different.
        constexpr bool operator!= (const TQuaternion& _q1) const noexcept
        {
            return r!=_q1.r || i!=_q1.i || j!=_q1.j || k!=_q1.k;
        }

        /// \brief TQuaternion multiplication is a combination of rotations.
        /// \details Non commutative (a*b != a*b)
        constexpr TQuaternion& operator *= (const TQuaternion& _q1) noexcept
        {
            // Preserve the sign for handness: the result can contain a mirroring, if
            // and only if one of the two arguments has a mirroring part.
            T handness = r*_q1.r;
            T nr = r*_q1.r - i*_q1.i - j*_q1.j - k*_q1.k;
            T ni = r*_q1.i + i*_q1.r + j*_q1.k - k*_q1.j;
            T nj = r*_q1.j + j*_q1.r + k*_q1.i - i*_q1.k;
            k = r*_q1.k + k*_q1.r + i*_q1.j - j*_q1.i;
            r = nr;
            i = ni;
            j = nj;
            if( r * handness < static_cast<T>(0) ) {i = -i; j = -j; k = -k; r = -r;}
            return *this;
        }

        /// \brief Scale the TQuaternion
        constexpr TQuaternion& operator *= (T _s) noexcept
        {
            eiAssert(_s >= T(0), "Using a negative scalar changes handness!");
            i*=_s; j*=_s; k*=_s; r*=_s;
            return *this;
        }

        /// \brief TQuaternion division   a/=b  <=>  a=a*(b^-1)=a*conjugated(b).
        constexpr TQuaternion& operator /= (const TQuaternion& _q1) noexcept
        {
            T handness = r*_q1.r;
            T nr =   r*_q1.r + i*_q1.i + j*_q1.j + k*_q1.k;
            T ni = - r*_q1.i + i*_q1.r - j*_q1.k + k*_q1.j;
            T nj = - r*_q1.j + j*_q1.r - k*_q1.i + i*_q1.k;
            k = - r*_q1.k + k*_q1.r - i*_q1.j + j*_q1.i;
            r = nr;
            i = ni;
            j = nj;
            if( r * handness < static_cast<T>(0) ) {i = -i; j = -j; k = -k; r = -r;}
            return *this;
        }

        /// \brief Scale the TQuaternion
        constexpr TQuaternion& operator /= (T _s) noexcept
        {
            eiAssert(_s >= T(0), "Using a negative scalar changes handness!");
            i/=_s; j/=_s; k/=_s; r/=_s;
            return *this;
        }

        /// \brief Vector like addition
        constexpr TQuaternion& operator += (const TQuaternion& _q1) noexcept
        {
            i+=_q1.i; j+=_q1.j; k+=_q1.k; r+=_q1.r;
            return *this;
        }

        /// \brief Vector like subtraction
        constexpr TQuaternion& operator -= (const TQuaternion& _q1) noexcept
        {
            i-=_q1.i; j-=_q1.j; k-=_q1.k; r-=_q1.r;
            return *this;
        }

        constexpr TQuaternion operator * (const TQuaternion& _q1) const noexcept
        {
           // TQuaternion q0 = *this;
            return TQuaternion(*this) *= _q1;
        }
        constexpr TQuaternion operator * (T _s) const noexcept
        {
            return TQuaternion(*this) *= _s;
        }
        constexpr TQuaternion operator / (TQuaternion _q1) const noexcept
        {
            return _q1 /= *this;
        }
        constexpr TQuaternion operator / (T _s) const noexcept
        {
            return TQuaternion(*this) /= _s;
        }
        constexpr TQuaternion operator + (TQuaternion _q1) const noexcept
        {
            return _q1 += *this;
        }
        constexpr TQuaternion operator - (TQuaternion _q1) const noexcept
        {
            return _q1 -= *this;
        }

        /// \brief Negate all components, the represented rotation is the same
        constexpr TQuaternion operator - () const noexcept // TESTED
        {
            return TQuaternion<T>(-i, -j, -k, -r);
        }

        /// \brief Conjugate the quaternion
        constexpr TQuaternion operator ~ () const noexcept // TESTED
        {
            return TQuaternion<T>(-i, -j, -k, r);
        }

        union {
            struct {T i, j, k, r;};         ///< Elements of 4D complex number
            T z[4];                         ///< Array access, index 3 is the real part
        };
    };

    // ********************************************************************* //
    /// \brief Returns identity element of the Hamilton-product. (Does not
    ///     rotate anything.)
    constexpr inline const TQuaternion<float> qidentity() noexcept // TESTED
    {
        return ei::TQuaternion<float>(0.0f, 0.0f, 0.0f, 1.0f);
    }

    constexpr inline const TQuaternion<double> qidentityD() noexcept
    {
        return ei::TQuaternion<double>(0.0, 0.0, 0.0, 1.0);
    }

    // ********************************************************************* //
    /// \brief Scalar multiplication from left
    template<typename T>
    constexpr inline TQuaternion<T> operator* (T _s, TQuaternion<T> _q) noexcept
    {
        return _q *= _s;
    }

    // ********************************************************************* //
    /// \brief Complex conjugate: invert sign of complex components
    template<typename T>
    constexpr inline TQuaternion<T> conjugate(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return TQuaternion<T>(-_q.i, -_q.j, -_q.k, _q.r);
    }

    /// \brief Get the rotation axis from a TQuaternion
    template<typename T>
    constexpr inline Vec<T,3> axis(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return Vec<T,3>(-_q.i, -_q.j, -_q.k) / max(T(EPSILON), std::sqrt(T(1)-_q.r*_q.r));
    }

    /// \brief Get the x axis of the corresponding orthogonal system (rotation
    ///     matrix)
    template<typename T>
    constexpr inline Vec<T,3> xaxis(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return Vec<T,3>( T(1)-T(2)*(_q.j*_q.j+_q.k*_q.k), T(2)*(_q.i*_q.j-_q.k*_q.r), T(2)*(_q.i*_q.k+_q.j*_q.r) );
    }

    /// \brief Get the y axis of the corresponding orthogonal system (rotation
    ///     matrix)
    template<typename T>
    constexpr inline Vec<T,3> yaxis(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return Vec<T,3>( T(2)*(_q.i*_q.j+_q.k*_q.r), T(1)-T(2)*(_q.i*_q.i+_q.k*_q.k), T(2)*(_q.j*_q.k-_q.i*_q.r) );
    }

    /// \brief Get the z axis of the corresponding orthogonal system (rotation
    ///     matrix)
    template<typename T>
    constexpr inline Vec<T,3> zaxis(const TQuaternion<T>& _q) noexcept // TESTED
    {
        T h = _q.r < T(0) ? T(-1) : T(1);
        T h2 = h * 2;
        return Vec<T,3>( h2*(_q.i*_q.k-_q.j*_q.r), h2*(_q.j*_q.k+_q.i*_q.r), h-h2*(_q.i*_q.i+_q.j*_q.j) );
    }
    // TODO: row vector axis

    /// \brief Get the angle (radians) from a TQuaternion
    template<typename T>
    constexpr inline T angle(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return acos(_q.r) * T(2);
    }

    // ********************************************************************* //
    /// \brief Get the Euler angles (radians) from a quaternion
    template<typename T>
    constexpr inline Vec<T,3> angles(const TQuaternion<T>& _q) noexcept
    {
        // TODO: handness?
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

    // ********************************************************************* //
    /// \brief Check if the absolute difference between all elements is smaller
    ///    or equal than epsilon.
    template<typename T>
    constexpr bool approx(const TQuaternion<T>& _q0,
                          const TQuaternion<T>& _q1,
                          T _epsilon = T(1e-6)) noexcept // TESTED
    {
        return abs(_q0.r - _q1.r) <= _epsilon
            && abs(_q0.i - _q1.i) <= _epsilon
            && abs(_q0.j - _q1.j) <= _epsilon
            && abs(_q0.k - _q1.k) <= _epsilon;
    }

    // ********************************************************************* //
    /// \brief Computes the sum of component wise products.
    /// \returns Scalar value of the sum of component products.
    constexpr inline float dot(const Quaternion& _q0,
                               const Quaternion& _q1) noexcept
    {
        return _q0.r*_q1.r + _q0.i*_q1.i + _q0.j*_q1.j + _q0.k*_q1.k;
    }

    // ********************************************************************* //
    /// \brief Spherical linear interpolation with constant angular speed
    template<typename T>
    constexpr TQuaternion<T> slerp(const TQuaternion<T>& _q0, const TQuaternion<T>& _q1, T _t) noexcept // TESTED
    {
        // TODO: handness?
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

    // ********************************************************************* //
    /// \brief Rotation matrix from quaternion.
    constexpr inline Mat3x3 rotation( const Quaternion& _quaternion ) noexcept
    {
        return Mat3x3(_quaternion);
    }

    // ********************************************************************* //
    /// \brief Apply a rotation by a quaternion (q v q-1 with v=(0, _v.x, _v.y, _v.z)).
    template<typename T, unsigned M, unsigned N, typename = std::enable_if_t<(M==1) || (N==1)>>
    constexpr inline Matrix<T,M,N> transform( const Matrix<T,M,N>& _what, const TQuaternion<T>& _quaternion ) noexcept
    {
        T handness = _quaternion.r < T(0) ? T(-1) : T(1);
        // http://physicsforgames.blogspot.de/2010/03/quaternion-tricks.html
        T x1 = _quaternion.j*_what.z - _quaternion.k*_what.y;
        T y1 = _quaternion.k*_what.x - _quaternion.i*_what.z;
        T z1 = _quaternion.i*_what.y - _quaternion.j*_what.x;

        return Matrix<T,M,N>(
             _what.x + 2.0f * (_quaternion.r*x1 + _quaternion.j*z1 - _quaternion.k*y1),
             _what.y + 2.0f * (_quaternion.r*y1 + _quaternion.k*x1 - _quaternion.i*z1),
            (_what.z + 2.0f * (_quaternion.r*z1 + _quaternion.i*y1 - _quaternion.j*x1)) * handness
            );
    }

} // namespace ei
