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
        constexpr EIAPI TQuaternion( const Vec<T,3>& _axis, T _angle ) noexcept // TESTED
        {
            eiAssert( approx(lensq(_axis), static_cast<T>(1)), "Expected a normalized axis vector!" );
            _angle *= 0.5f;
            T sinA = sin(_angle);
            r = cos(_angle);
            i = sinA * _axis.x;
            j = sinA * _axis.y;
            k = sinA * _axis.z;
        }

        /// \brief Create from Euler angles
        /// \details The rotations are applied in the order x, y, z:
        ///     rotationZ(_z) * rotationY(_y) * rotationX(_x)
        constexpr EIAPI TQuaternion( T _x, T _y, T _z ) noexcept // TESTED
        {
            double halfAngle = _x * 0.5;
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

            //*this = normalize(*this);
        }
        constexpr EIAPI TQuaternion( const Vec<T,3>& _eulerAngles ) noexcept :
            TQuaternion(_eulerAngles.x, _eulerAngles.y, _eulerAngles.z)
        {}

        /// \brief Create from rotation matrix (does a decomposition if the
        ///     matrix contains scaling).
        constexpr EIAPI TQuaternion( const Matrix<T,3,3>& _matrix ) noexcept : // TESTED
            TQuaternion<T>(transpose(_matrix(0)), transpose(_matrix(1)), transpose(_matrix(2)))
        {}

        /// \brief Create from orthogonal basis vectors.
        constexpr EIAPI TQuaternion( const Vec<T,3>& _xAxis, const Vec<T,3>& _yAxis, const Vec<T,3>& _zAxis ) noexcept 
        {
            // Check handness
            //eiAssert(dot(cross(_xAxis, _yAxis), _m(2)) > 0.0f, "Quaternions cannot handle reflections. The matrix must be RHS.");
#if defined(DEBUG) || defined(_DEBUG)
            T handness = dot(cross(_xAxis, _yAxis), _zAxis);
#endif
            eiAssert(handness > 0.0f, "Quaternions cannot handle reflections. The matrix must be RHS.");
            eiAssert(approx(handness, T(1), 1e-4f), "System is not orthonormal!");

            // Build TQuaternion<T> from rotation matrix
            // Src: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
            T trace = _xAxis.x + _yAxis.y + _zAxis.z;
            if( trace > 0 )
            {
                float s = T(0.5) / sqrt( trace + T(1) );
                i = ( _zAxis.y - _yAxis.z ) * s;
                j = ( _xAxis.z - _zAxis.x ) * s;
                k = ( _yAxis.x - _xAxis.y ) * s;
                r = T(0.25) / s;
            } else {
                if( _xAxis.x > _yAxis.y && _xAxis.x > _zAxis.z )
                {
                    float s = T(2) * sqrt( T(1) + _xAxis.x - _yAxis.y - _zAxis.z );
                    i = T(0.25) * s;
                    j = ( _xAxis.y + _yAxis.x ) / s;
                    k = ( _xAxis.z + _zAxis.x ) / s;
                    r = ( _zAxis.y - _yAxis.z ) / s;
                } else if( _yAxis.y > _zAxis.z )
                {
                    float s = T(2) * sqrt( T(1) + _yAxis.y - _xAxis.x - _zAxis.z );
                    i = ( _xAxis.y + _yAxis.x ) / s;
                    j = T(0.25) * s;
                    k = ( _yAxis.z + _zAxis.y ) / s;
                    r = ( _xAxis.z - _zAxis.x ) / s;
                } else {
                    float s = T(2) * sqrt( T(1) + _zAxis.z - _xAxis.x - _yAxis.y );
                    i = ( _xAxis.z + _zAxis.x ) / s;
                    j = ( _yAxis.z + _zAxis.y ) / s;
                    k = T(0.25) * s;
                    r = ( _yAxis.x - _xAxis.y ) / s;
                }
            }//*/

             /*r = sqrt( max( T(0), T(1) + _m.m00 + _m.m11 + _m.m22 ) ) * T(0.5);
             i = sqrt( max( T(0), T(1) + _m.m00 - _m.m11 - _m.m22 ) ) * T(0.5) * sgn(_m.m21 - _m.m12);
             j = sqrt( max( T(0), T(1) - _m.m00 + _m.m11 - _m.m22 ) ) * T(0.5) * sgn(_m.m02 - _m.m20);
             k = sqrt( max( T(0), T(1) - _m.m00 - _m.m11 + _m.m22 ) ) * T(0.5) * sgn(_m.m10 - _m.m01);//*/

            *this = normalize(*this);
        }

        /// \brief Create from TQuaternion coefficients
        constexpr EIAPI TQuaternion( T _i, T _j, T _k, T _r ) noexcept :
            i(_i), j(_j), k(_k), r(_r)
        {}

        /// \brief Rotate from vector to vector (rotated such that the from
        ///     vector is aligned with the to vector).
        /// \param [in] _from One certain direction vector before rotation.
        /// \param [in] _to Target direction vector. The from direction should
        ///     be aligned with the target after rotation.
        constexpr EIAPI TQuaternion( const Vec<T,3>& _from, const Vec<T,3>& _to ) noexcept // TESTED
        {
            eiAssert(approx(len(_from),1.0f), "Input (_from) must be normalized direction vector.");
            eiAssert(approx(len(_to),1.0f), "Input (_to) must be normalized direction vector.");
            // half angle trick from http://physicsforgames.blogspot.de/2010/03/quaternion-tricks.html
            Vec<T,3> half = (_from + _to);
            T hlensq = lensq(half);
            // Opposite vectors or one vector 0.0 -> 180° rotation
            if(hlensq < 1e-10f)
            {
                half = perpendicular(_from);
                hlensq = lensq(half);
            }
            half /= sqrt(hlensq);

            // cos(theta) = dot product since both vectors are normalized
            r = dot(_from, half);
            // Axis from cross product -> already multiplied with sin(theta)
            i = _from.y*half.z - _from.z*half.y;
            j = _from.z*half.x - _from.x*half.z;
            k = _from.x*half.y - _from.y*half.x;
        }

        /// \brief Compare component wise, if two quaternions are identical.
        constexpr EIAPI bool operator == (const TQuaternion& _q1) const noexcept
        {
            return (r== _q1.r && i== _q1.i && j== _q1.j && k== _q1.k)
                || (r==-_q1.r && i==-_q1.i && j==-_q1.j && k==-_q1.k);
        }
        /// \brief Compare component wise, if two quaternions are different.
        constexpr EIAPI bool operator!= (const TQuaternion& _q1) const noexcept
        {
            return (r!= _q1.r || i!= _q1.i || j!= _q1.j || k!= _q1.k)
                && (r!=-_q1.r || i!=-_q1.i || j!=-_q1.j || k!=-_q1.k);
        }

        template<typename T1>
        constexpr EIAPI explicit operator Matrix<T1,3,3>() const noexcept
        {
            // Rotation composition from quaternion (remaining rest direct in matrix)
            // See http://de.wikipedia.org/wiki/Quaternion#Bezug_zu_orthogonalen_Matrizen for
            // details.
            T f2i  = 2.0f * i;
            T f2j  = 2.0f * j;
            T f2k  = 2.0f * k;
            T f2ri = f2i  * r;
            T f2rj = f2j  * r;
            T f2rk = f2k  * r;
            T f2ii = f2i  * i;
            T f2ij = f2j  * i;
            T f2ik = f2k  * i;
            T f2jj = f2j  * j;
            T f2jk = f2k  * j;
            T f2kk = f2k  * k;

            return Matrix<T1,3,3>{
                1.0f - ( f2jj + f2kk ), f2ij - f2rk,            f2ik + f2rj,
                f2ij + f2rk,            1.0f - ( f2ii + f2kk ), f2jk - f2ri,
                f2ik - f2rj,            f2jk + f2ri,            1.0f - ( f2ii + f2jj )
            };
        }

        /// \brief TQuaternion multiplication is a combination of rotations.
        /// \details Non commutative (a*b != b*a)
        constexpr EIAPI TQuaternion& operator *= (const TQuaternion& _q1) noexcept // TESTED
        {
            T nr = r*_q1.r - i*_q1.i - j*_q1.j - k*_q1.k;
            T ni = r*_q1.i + i*_q1.r + j*_q1.k - k*_q1.j;
            T nj = r*_q1.j + j*_q1.r + k*_q1.i - i*_q1.k;
            k    = r*_q1.k + k*_q1.r + i*_q1.j - j*_q1.i;
            r = nr;
            i = ni;
            j = nj;
            return *this;
        }

        /// \brief Scale the TQuaternion
        constexpr EIAPI TQuaternion& operator *= (T _s) noexcept
        {
            i*=_s; j*=_s; k*=_s; r*=_s;
            return *this;
        }

        /// \brief TQuaternion division   a/=b  <=>  a=a*(b^-1)=a*conjugated(b).
        constexpr EIAPI TQuaternion& operator /= (const TQuaternion& _q1) noexcept
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

        /// \brief Scale the TQuaternion
        constexpr EIAPI TQuaternion& operator /= (T _s) noexcept
        {
            i/=_s; j/=_s; k/=_s; r/=_s;
            return *this;
        }

        /// \brief Vector like addition
        constexpr EIAPI TQuaternion& operator += (const TQuaternion& _q1) noexcept
        {
            i+=_q1.i; j+=_q1.j; k+=_q1.k; r+=_q1.r;
            return *this;
        }

        /// \brief Vector like subtraction
        constexpr EIAPI TQuaternion& operator -= (const TQuaternion& _q1) noexcept
        {
            i-=_q1.i; j-=_q1.j; k-=_q1.k; r-=_q1.r;
            return *this;
        }

        constexpr EIAPI TQuaternion operator * (const TQuaternion& _q1) const noexcept // TESTED
        {
            return TQuaternion(*this) *= _q1;
        }
        constexpr EIAPI TQuaternion operator * (T _s) const noexcept
        {
            return TQuaternion(*this) *= _s;
        }
        constexpr TQuaternion operator / (TQuaternion _q1) const noexcept
        {
            return _q1 /= *this;
        }
        constexpr EIAPI TQuaternion operator / (T _s) const noexcept
        {
            return TQuaternion(*this) /= _s;
        }
        constexpr EIAPI TQuaternion operator + (TQuaternion _q1) const noexcept
        {
            return _q1 += *this;
        }
        constexpr TQuaternion operator - (TQuaternion _q1) const noexcept
        {
            return _q1 -= *this;
        }

        /// \brief Negate all components, the represented rotation is the same
        constexpr EIAPI TQuaternion operator - () const noexcept // TESTED
        {
            return TQuaternion<T>(-i, -j, -k, -r);
        }

        /// \brief Conjugate the quaternion
        constexpr EIAPI TQuaternion operator ~ () const noexcept // TESTED
        {
            return TQuaternion<T>(-i, -j, -k, r);
        }

        union {
            struct {T i, j, k, r;};         ///< Elements of 4D complex number
            T z[4];                         ///< Array access, index 3 is the real part
            struct { Vec<T, 3> complex; float real; };   ///< Access to the complex part as a vector
        };
    };

    // ********************************************************************* //
    /// \brief Returns identity element of the Hamilton-product. (Does not
    ///     rotate anything.)
    constexpr EIAPI const TQuaternion<float> qidentity() noexcept // TESTED
    {
        return ei::TQuaternion<float>(0.0f, 0.0f, 0.0f, 1.0f);
    }

    constexpr EIAPI const TQuaternion<double> qidentityD() noexcept
    {
        return ei::TQuaternion<double>(0.0, 0.0, 0.0, 1.0);
    }

    // ********************************************************************* //
    /// \brief Scalar multiplication from left
    template<typename T>
    constexpr EIAPI TQuaternion<T> operator* (T _s, TQuaternion<T> _q) noexcept
    {
        return _q *= _s;
    }

    // ********************************************************************* //
    /// \brief Complex conjugate: invert sign of complex components
    template<typename T>
    constexpr EIAPI TQuaternion<T> conjugate(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return TQuaternion<T>(-_q.i, -_q.j, -_q.k, _q.r);
    }

    /// \brief Get the rotation axis from a TQuaternion
    template<typename T>
    constexpr EIAPI Vec<T,3> axis(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return Vec<T,3>(_q.i, _q.j, _q.k) / std::sqrt(max(T(EPSILON), T(1)-_q.r*_q.r));
    }

    /// \brief Get the first row of the corresponding rotation matrix
    template<typename T>
    constexpr EIAPI Vec<T,3> xRow(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return Vec<T,3>( T(1)-T(2)*(_q.j*_q.j+_q.k*_q.k), T(2)*(_q.i*_q.j-_q.k*_q.r), T(2)*(_q.i*_q.k+_q.j*_q.r) );
    }

    /// \brief Get the second row of the corresponding rotation matrix
    template<typename T>
    constexpr EIAPI Vec<T,3> yRow(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return Vec<T,3>( T(2)*(_q.i*_q.j+_q.k*_q.r), T(1)-T(2)*(_q.i*_q.i+_q.k*_q.k), T(2)*(_q.j*_q.k-_q.i*_q.r) );
    }

    /// \brief Get the third row of the corresponding rotation matrix
    template<typename T>
    constexpr EIAPI Vec<T,3> zRow(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return Vec<T,3>( T(2)*(_q.i*_q.k-_q.j*_q.r), T(2)*(_q.j*_q.k+_q.i*_q.r), T(1)-T(2)*(_q.i*_q.i+_q.j*_q.j) );
    }

    /// \brief Get the first column of the corresponding rotation matrix.
    /// \details This is equivalent the x vector (1,0,0) transformed into the new space.
    template<typename T>
    constexpr EIAPI Vec<T,3> xCol(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return Vec<T,3>( T(1)-T(2)*(_q.j*_q.j+_q.k*_q.k), T(2)*(_q.i*_q.j+_q.k*_q.r), T(2)*(_q.i*_q.k-_q.j*_q.r) );
    }

    /// \brief Get the second column of the corresponding rotation matrix.
    /// \details This is equivalent the y vector (0,1,0) transformed into the new space.
    template<typename T>
    constexpr EIAPI Vec<T,3> yCol(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return Vec<T,3>( T(2)*(_q.i*_q.j-_q.k*_q.r), T(1)-T(2)*(_q.i*_q.i+_q.k*_q.k), T(2)*(_q.j*_q.k+_q.i*_q.r) );
    }

    /// \brief Get the third column of the corresponding rotation matrix.
    /// \details This is equivalent the z vector (0,0,1) transformed into the new space.
    template<typename T>
    constexpr EIAPI Vec<T,3> zCol(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return Vec<T,3>( T(2)*(_q.i*_q.k+_q.j*_q.r), T(2)*(_q.j*_q.k-_q.i*_q.r), T(1)-T(2)*(_q.i*_q.i+_q.j*_q.j) );
    }

    /// \brief Get the angle (radians) from a TQuaternion
    template<typename T>
    constexpr EIAPI T angle(const TQuaternion<T>& _q) noexcept // TESTED
    {
        return acos(_q.r) * T(2);
    }

    // ********************************************************************* //
    /// \brief Get the Euler angles (radians) from a quaternion
    template<typename T>
    constexpr EIAPI Vec<T,3> angles(const TQuaternion<T>& _q) noexcept
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
    constexpr EIAPI bool approx(const TQuaternion<T>& _q0,
                                const TQuaternion<T>& _q1,
                                T _epsilon = T(1e-6)) noexcept // TESTED
    {
        TQuaternion<T> qt = _q0.r * _q1.r < 0.0f ? -_q1 : _q1;
        return abs(_q0.r - qt.r) <= _epsilon
            && abs(_q0.i - qt.i) <= _epsilon
            && abs(_q0.j - qt.j) <= _epsilon
            && abs(_q0.k - qt.k) <= _epsilon;
    }

    // ********************************************************************* //
    /// \brief Computes the sum of component wise products.
    /// \returns Scalar value of the sum of component products.
    constexpr EIAPI float dot(const Quaternion& _q0,
                               const Quaternion& _q1) noexcept
    {
        return _q0.r*_q1.r + _q0.i*_q1.i + _q0.j*_q1.j + _q0.k*_q1.k;
    }

    // ********************************************************************* //
    /// \brief Spherical linear interpolation with constant angular speed
    template<typename T>
    constexpr EIAPI TQuaternion<T> slerp(const TQuaternion<T>& _q0, const TQuaternion<T>& _q1, T _t) noexcept // TESTED
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

    // ********************************************************************* //
    /// \brief Rotation matrix from quaternion.
    constexpr EIAPI Mat3x3 rotation( const Quaternion& _quaternion ) noexcept
    {
        return Mat3x3(_quaternion);
    }

    // ********************************************************************* //
    /// \brief Apply a rotation by a quaternion (q v q-1 with v=(0, _v.x, _v.y, _v.z)).
    template<typename T, unsigned M, unsigned N, typename = std::enable_if_t<(M==1) || (N==1)>>
    constexpr EIAPI Matrix<T,M,N> transform( const Matrix<T,M,N>& _what, const TQuaternion<T>& _quaternion ) noexcept
    {
        // http://physicsforgames.blogspot.de/2010/03/quaternion-tricks.html
        T x1 = _quaternion.j*_what.z - _quaternion.k*_what.y;
        T y1 = _quaternion.k*_what.x - _quaternion.i*_what.z;
        T z1 = _quaternion.i*_what.y - _quaternion.j*_what.x;

        return Matrix<T,M,N>(
             _what.x + 2.0f * (_quaternion.r*x1 + _quaternion.j*z1 - _quaternion.k*y1),
             _what.y + 2.0f * (_quaternion.r*y1 + _quaternion.k*x1 - _quaternion.i*z1),
             _what.z + 2.0f * (_quaternion.r*z1 + _quaternion.i*y1 - _quaternion.j*x1)
            );
    }



    // ********************************************************************* //
    //                       ORTHORNORMAL SPACE TYPE                         //
    // ********************************************************************* //
    /// \brief Storage class to store any orthonormal space, including those
    ///     with det() = -1 (including a reflection).
    /// \details There is no concatenation or other operation on the compressed
    ///     spaces because they would be non trivial. The compression always
    ///     assumes the mirroring to take place in the z-axis. However,
    ///     rotating another orthonormal space changes this direction and it
    ///     would be necessary to remember the reflection direction OR to
    ///     compute a fully new rotation with z-axis reflection.
    ///     The first option would render the compression useless and the
    ///     second is essentially OrthoSpace(Mat3x3(ortho1) * Mat3x3(ortho2)).
    ///     Since this is expensive it is better not to hide its costs.
    ///     If you are sure that there are only rotations use the faster
    ///     OrthoSpace(Quaternion(ortho1) * Quaternion(ortho2)).
    template<typename T>
    class TOrthoSpace : public details::NonScalarType
    {
    public:
        TOrthoSpace() = default;

        /// \brief Initialize from normalized quaternion
        constexpr EIAPI explicit TOrthoSpace(const TQuaternion<T>& _q) noexcept : // TESTED
            m_quaternion(_q * sgn(_q.r))
        {
            eiAssert(approx(len(_q), T(1)), "Quaternion must be normalized.");
        }

        /// \brief Initialize from any 3x3 orthonorml 
        constexpr EIAPI explicit TOrthoSpace(const Matrix<T,3,3>& _m) noexcept  // TESTED
        {
            T handness = dot(cross(_m(0), _m(1)), _m(2));
            m_quaternion = TQuaternion<T>(transpose(_m(0)), transpose(_m(1)), transpose(handness * _m(2)));
            // Assert additional normalization condition for the handness
            if(sgn(m_quaternion.r) < static_cast<T>(0))
                m_quaternion = -m_quaternion;
            // Store the handness
            if(handness < static_cast<T>(0))
                m_quaternion.r = -m_quaternion.r;
        }

        /// \brief Reconstruct the full orthonorml system.
        template<typename T1>
        constexpr EIAPI explicit operator Matrix<T1,3,3>() const noexcept {  // TESTED
            // Remove the handness sign and get the rotation matrix
            Matrix<T1,3,3> rot{TQuaternion<T1>(m_quaternion.i, m_quaternion.j, m_quaternion.k, ei::abs(m_quaternion.r))};
            if(isLefthanded())
                rot(2) = -rot(2);
            return rot;
        }

        /// \brief Get the rotation. If the stored system contains a reflection it is lost.
        template<typename T1>
        constexpr EIAPI explicit operator TQuaternion<T1>() const noexcept {  // TESTED
            return TQuaternion<T1>(m_quaternion.i, m_quaternion.j, m_quaternion.k, ei::abs(m_quaternion.r));
        }

        constexpr EIAPI bool isRighthanded() const { return sgn(m_quaternion.r) == 1.0f; }  // TESTED
        constexpr EIAPI bool isLefthanded() const { return sgn(m_quaternion.r) == -1.0f; }  // TESTED

        constexpr EIAPI bool operator == (const TOrthoSpace& _other) const {
            // Unlike quaternions where q = -q we have a unique representation.
            return m_quaternion.r == _other.m_quaternion.r
                && m_quaternion.i == _other.m_quaternion.i
                && m_quaternion.j == _other.m_quaternion.j
                && m_quaternion.k == _other.m_quaternion.k;
        }

        constexpr EIAPI bool operator != (const TOrthoSpace& _other) const {
            // Unlike quaternions where q = -q we have a unique representation.
            return m_quaternion.r != _other.m_quaternion.r
                || m_quaternion.i != _other.m_quaternion.i
                || m_quaternion.j != _other.m_quaternion.j
                || m_quaternion.k != _other.m_quaternion.k;
        }

        // Get access to the internal quaternion.
        // WARNING: the quternion has additional information encoded and cannot be
        // used as a quaternion.
        constexpr EIAPI const TQuaternion<T>& data() const noexcept { return m_quaternion; }
        constexpr EIAPI TQuaternion<T>& data() noexcept { return m_quaternion; }
    private:
        // Use a standard quaternion.
        // The sign of the determinant can be encoded in the r component.
        // Since +q=-q we can make sure to use the one with q.r>0. Then the
        // sign of the determinant is simply put into the guarenteed positive r.
        // Important: we need -0 and +0 from float to be sure the sign is always
        // encoded.
        TQuaternion<T> m_quaternion;
    };

    using OrthoSpace = TOrthoSpace<float>;

    // ********************************************************************* //
    /// \brief Check if the absolute difference between all elements is smaller
    ///    or equal than epsilon.
    template<typename T>
    constexpr EIAPI bool approx(const TOrthoSpace<T>& _o0,
                          const TOrthoSpace<T>& _o1,
                          T _epsilon = T(1e-6)) noexcept // TESTED
    {
        return abs(_o0.data().r - _o1.data().r) <= _epsilon
            && abs(_o0.data().i - _o1.data().i) <= _epsilon
            && abs(_o0.data().j - _o1.data().j) <= _epsilon
            && abs(_o0.data().k - _o1.data().k) <= _epsilon;
    }

} // namespace ei
