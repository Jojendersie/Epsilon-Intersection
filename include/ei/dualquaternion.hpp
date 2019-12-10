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
            r = i = j = k = re = ie = je = ke = details::IniVal<T>::value;
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
            q0{_q0},
            qe{TQuaternion<T>{_translation.x * 0.5f, _translation.y * 0.5f, _translation.z * 0.5f, 0.0f} * _q0}
        {}

        template<typename T1>
        constexpr EIAPI explicit operator Matrix<T1,3,4>() const noexcept // TESTED
        {
            constexpr T1 ONE = static_cast<T1>(1);
            constexpr T1 TWO = static_cast<T1>(2);
            T1 t0 = TWO * (-qe.r*q0.i + qe.i*q0.r - qe.j*q0.k + qe.k*q0.j);
            T1 t1 = TWO * (-qe.r*q0.j + qe.i*q0.k + qe.j*q0.r - qe.k*q0.i);
            T1 t2 = TWO * (-qe.r*q0.k - qe.i*q0.j + qe.j*q0.i + qe.k*q0.r);
            T1 f2i  = TWO * q0.i;
            T1 f2j  = TWO * q0.j;
            T1 f2k  = TWO * q0.k;
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
    };


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

} // namespace ei