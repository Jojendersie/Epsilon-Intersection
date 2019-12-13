#pragma once

#include "quaternion.hpp"

namespace ei {

    /// \brief A dual quaternion consits of two quaternions q0 + qε where qε
    ///     has the unit ε with the property ε² = 0.
    /// \details Dual quaternions can describe a rotation and a translation
    ///     together. Thereby q0 = (cos(θ/2) + (x i,y j,z k) sin(θ/2)) describes
    ///     the (half) rotation around axis (x,y,z).
    ///     Analogously qε = 1 + ε/2(tx i,ty j,tz k) describes a (half) translation
    ///     by (tx, ty, tz).
    ///     To represent some vector v in dual quaternion space one can write
    ///     q = 1 + ε(vx i, vy j, vz k)
    ///
    ///     Also see "Skinning with Dual Quaternions"
    ///     https://www.cs.utah.edu/~ladislav/kavan07skinning/kavan07skinning.pdf
    template<typename T>
    class TDualQuaternion: public details::NonScalarType
    {
    public:
        TQuaternion<T> q0, qe;

        /// \brief Construct uninitialized
        EIAPI TDualQuaternion() noexcept {
#if defined(DEBUG) || defined(_DEBUG)
            q0.r = q0.i = q0.j = q0.k = qe.r = qe.i = qe.j = qe.k = details::IniVal<T>::value;
#endif
            // Empty on release
        }

        /// \brief Copy construction
        constexpr TDualQuaternion( const TDualQuaternion& _other ) noexcept = default;
        /// \brief Copying assignment
        constexpr TDualQuaternion& operator = ( const TDualQuaternion& _rhs ) noexcept = default;

        /// \brief Create from coefficients
        constexpr EIAPI TDualQuaternion( T _i, T _j, T _k, T _r, T _ie, T _je, T _ke, T _re ) noexcept :
            q0{_i, _j, _k, _r},
            qe{_ie,_je,_ke,_re}
        {}

        /// \brief Create from quaternion and translation
        constexpr TDualQuaternion( const TQuaternion<T>& _q0, const Vec<T,3>& _translation ) : // TESTED
            q0{ _q0 },
            qe{ TQuaternion<T>{_translation.x * 0.5f, _translation.y * 0.5f, _translation.z * 0.5f, 0.0f} * _q0 }
        {}

        /// \brief Create from 3x4 matrix
        template<typename T1>
        constexpr TDualQuaternion( const Matrix<T1,3,4>& _mat ) noexcept : // TESTED
            TDualQuaternion( TQuaternion<T>{ Matrix<T1,3,3>{_mat} }, Vec<T,3>{_mat.m03, _mat.m13, _mat.m23} )
        {}

        template<typename T1>
        constexpr EIAPI explicit operator Matrix<T1,3,4>() const noexcept // TESTED
        {
            eiAssertWeak( approx(static_cast<T>(1), dot(q0, q0), 1e-4f), "Not a union dual quaternion!" );
            constexpr T1 ONE = static_cast<T1>(1);
            constexpr T1 TWO = static_cast<T1>(2);
            T1 f2i = TWO * q0.i;
            T1 f2j = TWO * q0.j;
            T1 f2k = TWO * q0.k;
            T1 f2r = TWO * q0.r;
            T1 t0 = -qe.r*f2i + qe.i*f2r - qe.j*f2k + qe.k*f2j;
            T1 t1 = -qe.r*f2j + qe.i*f2k + qe.j*f2r - qe.k*f2i;
            T1 t2 = -qe.r*f2k - qe.i*f2j + qe.j*f2i + qe.k*f2r;
            T1 f2ri = f2i * q0.r;
            T1 f2rj = f2j * q0.r;
            T1 f2rk = f2k * q0.r;
            T1 f2ii = f2i * q0.i;
            T1 f2ij = f2j * q0.i;
            T1 f2ik = f2k * q0.i;
            T1 f2jj = f2j * q0.j;
            T1 f2jk = f2k * q0.j;
            T1 f2kk = f2k * q0.k;

            return Matrix<T1,3,4>{
                ONE - ( f2jj + f2kk ), f2ij - f2rk,           f2ik + f2rj,           t0,
                f2ij + f2rk,           ONE - ( f2ii + f2kk ), f2jk - f2ri,           t1,
                f2ik - f2rj,           f2jk + f2ri,           ONE - ( f2ii + f2jj ), t2
            };
        }

        /// \brief Compare component wise, if two quaternions are identical.
        constexpr EIAPI bool operator == (const TDualQuaternion& _q1) const noexcept // TESTED
        {
            return _q1.q0 == q0 && _q1.qe == qe;
        }
        /// \brief Compare component wise, if two quaternions are different.
        constexpr EIAPI bool operator!= (const TDualQuaternion& _q1) const noexcept
        {
            return _q1.q0 != q0 || _q1.qe != qe;
        }

        /// \brief TDualQuaternion multiplication is a combination of the transformations.
        /// \details Non commutative (a*b != b*a)
        constexpr EIAPI TDualQuaternion& operator *= (const TDualQuaternion& _q1) noexcept // TESTED
        {
            qe = q0 * _q1.qe + qe * _q1.q0;
            q0 = q0 * _q1.q0;
            return *this;
        }

        /// \brief Scale the TQuaternion
        constexpr EIAPI TDualQuaternion& operator *= (T _s) noexcept
        {
            q0 *= _s; qe *= _s;
            return *this;
        }

        /// \brief TDualQuaternion division (expects unit quaternions).
        /// \details Non commutative (a*b != b*a)
        ///     a/=b <=> a=a*(b⁻¹)=a*-(conjugate(b0)/len(b0) - ε conjugate(b0) bε conjugate(b0) / lensq(b0))
        ///      unit => a*-(conjugate(b0) - ε bε * conjugate(b0)²)
        constexpr EIAPI TDualQuaternion& operator /= (const TDualQuaternion& _q1) noexcept // TESTED
        {
            Quaternion q10inv = conjugate(_q1.q0);
            qe = (q0 * q10inv * _q1.qe - qe) * q10inv;
            q0 = q0 * q10inv;
            return *this;
        }

        /// \brief Divied the TQuaternion by scalar
        constexpr EIAPI TDualQuaternion& operator /= (T _s) noexcept
        {
            q0 /= _s; qe /= _s;
            return *this;
        }

        /// \brief Vector like addition
        constexpr EIAPI TDualQuaternion& operator += (const TDualQuaternion& _q1) noexcept
        {
            q0 += _q1.q0;
            qe += _q1.qe;
            return *this;
        }

        /// \brief Vector like subtraction
        constexpr EIAPI TDualQuaternion& operator -= (const TDualQuaternion& _q1) noexcept
        {
            q0 -= _q1.q0;
            qe -= _q1.qe;
            return *this;
        }

        constexpr EIAPI TDualQuaternion operator * (const TDualQuaternion& _q1) const noexcept // TESTED
        {
            return TDualQuaternion(*this) *= _q1;
        }
        constexpr EIAPI TDualQuaternion operator * (T _s) const noexcept
        {
            return TDualQuaternion(*this) *= _s;
        }
        constexpr TDualQuaternion operator / (TDualQuaternion _q1) const noexcept // TESTED
        {
            return TDualQuaternion(*this) /= _q1;
        }
        constexpr EIAPI TDualQuaternion operator / (T _s) const noexcept
        {
            return TDualQuaternion(*this) /= _s;
        }
        constexpr EIAPI TDualQuaternion operator + (TDualQuaternion _q1) const noexcept
        {
            return _q1 += *this;
        }
        constexpr TDualQuaternion operator - (TDualQuaternion _q1) const noexcept
        {
            return _q1 -= *this;
        }


        /// \brief Negate all components, the represented rotation is the same
        constexpr EIAPI TDualQuaternion operator - () const noexcept
        {
            return TDualQuaternion<T>(-q0.i, -q0.j, -q0.k, -q0.r, -qe.i, -qe.j, -qe.k, -qe.r);
        }

        /// \brief Complex conjugate the quaternions
        constexpr EIAPI TDualQuaternion operator ~ () const noexcept
        {
            return TDualQuaternion<T>(-q0.i, -q0.j, -q0.k, q0.r, -qe.i, -qe.j, -qe.k, qe.r);
        }
    };

    // ********************************************************************* //
    /// \brief Returns identity element of the Hamilton-product. (Does not
    ///     rotate or translate anything.)
    template<typename T = float>
    constexpr EIAPI const TDualQuaternion<T> qqidentity() noexcept
    {
        return ei::TDualQuaternion<T>(static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(1),
                                      static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
    }

    // ********************************************************************* //
    /// \brief Scalar multiplication from left
    template<typename T>
    constexpr EIAPI TDualQuaternion<T> operator* (T _s, TDualQuaternion<T> _q) noexcept
    {
        return _q *= _s;
    }


    // ********************************************************************* //
    /// \brief Complex conjugate: invert sign of complex components
    template<typename T>
    constexpr EIAPI TDualQuaternion<T> conjugate(const TDualQuaternion<T>& _q) noexcept // TESTED
    {
        return TDualQuaternion<T>(-_q.q0.i, -_q.q0.j, -_q.q0.k, _q.q0.r, -_q.qe.i, -_q.qe.j, -_q.qe.k, _q.qe.r);
    }

    // ********************************************************************* //
    /// \brief Dual conjugate: invert sign of dual components
    template<typename T>
    constexpr EIAPI TDualQuaternion<T> dualconjugate(const TDualQuaternion<T>& _q) noexcept // TESTED
    {
        return TDualQuaternion<T>(_q.q0.i, _q.q0.j, _q.q0.k, _q.q0.r, -_q.qe.i, -_q.qe.j, -_q.qe.k, -_q.qe.r);
    }

    // ********************************************************************* //
    /// \brief Check if the absolute difference between all elements is smaller
    ///    or equal than epsilon.
    template<typename T>
    constexpr EIAPI bool approx(const TDualQuaternion<T>& _q0,
                                const TDualQuaternion<T>& _q1,
                                T _epsilon = T(1e-6)) noexcept // TESTED
    {
        return approx(_q0.q0, _q1.q0, _epsilon) && approx(_q0.qe, _q1.qe, _epsilon);
    }

    // ********************************************************************* //
    /// \brief Computes the magnitude squared of the dual quatenrion.
    /// \returns Scalar value of the sum of component products.
    template<typename T>
    constexpr EIAPI float lensq(const TDualQuaternion<T>& _q) noexcept
    {
        // q q' where q' is the complex conjugate.
        // The full norm is: q0 q0' + ε (q0 qε' + qε q0')
        T lq0 = dot(_q.q0, _q.q0);
        eiAssertWeak( approx(static_cast<T>(0), dot(_q.q0, _q.qe)), "Not a union dual quaternion!" );
        return lq0;	// Return only real part
    }


    // ********************************************************************* //
    /// \brief Computes the magnitude of the dual quatenrion.
    /// \returns Scalar value of the sum of component products.
    template<typename T>
    constexpr EIAPI float len(const TDualQuaternion<T>& _q) noexcept // TESTED
    {
        // The full norm is len(q0) + ε dot(q0, qε) / len(q0).
        // Only the real part seems to be of interest (not sure about this)
        T lq0 = len(_q.q0);
        eiAssertWeak( approx(static_cast<T>(0), dot(_q.q0, _q.qe)), "Not a union dual quaternion!" );
        return lq0;
    }

    // ********************************************************************* //
    /// \brief Apply a rotation by a quaternion (ignore translation)
    template<typename T, unsigned M, unsigned N, typename = std::enable_if_t<(M==1) || (N==1)>>
    constexpr EIAPI Matrix<T,M,N> transformDir( const Matrix<T,M,N>& _what, const TDualQuaternion<T>& _qq ) noexcept // TESTED
    {
        return transform(_what, _qq.q0);
    }

    // ********************************************************************* //
    /// \brief Apply full transformation (rotation + translation).
    template<typename T, unsigned M, unsigned N, typename = std::enable_if_t<(M==1) || (N==1)>>
    constexpr EIAPI Matrix<T,M,N> transform( const Matrix<T,M,N>& _what, const TDualQuaternion<T>& _qq ) noexcept // TESTED
    {
        constexpr T TWO = static_cast<T>(2);
        const T f2i = TWO * _qq.q0.i;
        const T f2j = TWO * _qq.q0.j;
        const T f2k = TWO * _qq.q0.k;
        const T f2r = TWO * _qq.q0.r;
        //T t0 = -_qq.qe.r*f2i + _qq.qe.i*f2r - _qq.qe.j*f2k + _qq.qe.k*f2j;
        //T t1 = -_qq.qe.r*f2j + _qq.qe.i*f2k + _qq.qe.j*f2r - _qq.qe.k*f2i;
        //T t2 = -_qq.qe.r*f2k - _qq.qe.i*f2j + _qq.qe.j*f2i + _qq.qe.k*f2r;
        // http://physicsforgames.blogspot.de/2010/03/quaternion-tricks.html
        const T x1 = _qq.q0.j*_what.z - _qq.q0.k*_what.y;
        const T y1 = _qq.q0.k*_what.x - _qq.q0.i*_what.z;
        const T z1 = _qq.q0.i*_what.y - _qq.q0.j*_what.x;

        //return Matrix<T,M,N>(
        //     _what.x + f2r*x1 + f2j*z1 - f2k*y1 + t0,
        //     _what.y + f2r*y1 + f2k*x1 - f2i*z1 + t1,
        //     _what.z + f2r*z1 + f2i*y1 - f2j*x1 + t2
        //);

        // Optimized further by rearanging floats
        const T x1i = x1 + _qq.qe.i;
        const T y1j = y1 + _qq.qe.j;
        const T z1k = z1 + _qq.qe.k;
        return Matrix<T,M,N>(
            _what.x + f2r*x1i + f2j*z1k - f2k*y1j - f2i*_qq.qe.r,
            _what.y + f2r*y1j + f2k*x1i - f2i*z1k - f2j*_qq.qe.r,
            _what.z + f2r*z1k + f2i*y1j - f2j*x1i - f2k*_qq.qe.r
        );
    }

} // namespace ei