namespace details {

    /// \brief Specialized component access for small vectors and matrices.
    /// \details This is the fallback for larger vectors without single component
    ///     access.
    template<typename T, unsigned M, unsigned N> struct Components: public NonScalarType
    {
    protected:
        // For vectors with M==0 or N==0 take one dummy element
        T m_data[M * N < 1 ? 1 : N * M];

    public:
        Components() {}
        /// \brief Construct from exactly N*M arguments.
        /// \details The int argument is a dummy to prevent some compilers (vc120) from generating
        ///     two constructors with 0 arguments.
        template<typename T1, typename... Args>
        Components(T1 _a0, Args... _args)// : m_data{ _a0, T(_args)... }
        {
            static_assert(sizeof...(Args)+1 == M*N, "Wrong number of arguments!");
            init<0>(_a0, _args...);
        }
        /// \brief Initialize all members from single scalar
        template<typename T1>
        explicit Components(T1 _s)
        {
            for(uint i = 0; i < N * M; ++i)
                this->m_data[i] = static_cast<T>(_s);
        }

    private:
        // Helper to unroll the M*N argument construction. The Syntax with
        // : m_data{ _a0, T(_args)... } is much prettier, but does not work in vc120.
        template<int Index>
        void init() {}
        template<int Index, typename T1, typename... Args>
        void init(T1 _a0, Args... _args) { m_data[Index] = _a0; init<Index+1>(_args...); }
    };


    /// \brief Specialized version for 1 component row or column vector.
    template<typename T> struct Components<T, 1, 1>: public NonScalarType
    {
        union {
            T x;
            T r;
            T m_data[1];
        };

        Components() {}
        template<typename T1> explicit Components(T1 _s) : x(static_cast<T>(_s)) {}
    };

    /// \brief Specialized version for 2 component row and column vectors.
    template<typename T> struct Components<T, 2, 1>: public NonScalarType
    {
        union {
            struct { T x, y; };
            struct { T r, g; };
            T m_data[2];
        };

        Components() {}
        template<typename T1> explicit Components(T1 _s) : x(static_cast<T>(_s)), y(static_cast<T>(_s)) {}
        template<typename T1, typename T2>
        Components(T1 _s0, T2 _s1) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1)) {}
        /*template<typename T1>
        Components(T _s0, Components<T1, 1, 1> _s1) : x(_s0), y(_s1.x) {}
        template<typename T1>
        Components(Components<T1, 1, 1> _s0, T _s1) : x(_s0.x), y(_s1) {}*/
    };
    template<typename T> struct Components<T, 1, 2>: public NonScalarType
    {
        union {
            struct { T x, y; };
            struct { T r, g; };
            T m_data[2];
        };

        Components() {}
        template<typename T1> explicit Components(T1 _s) : x(static_cast<T>(_s)), y(static_cast<T>(_s)) {}
        template<typename T1, typename T2>
        Components(T1 _s0, T2 _s1) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1)) {}
        /*template<typename T1>
        Components(T _s0, Components<T1, 1, 1> _s1) : x(_s0), y(_s1.x) {}
        template<typename T1>
        Components(Components<T1, 1, 1> _s0, T _s1) : x(_s0.x), y(_s1) {}*/
    };

    /// \brief Specialized version for 3 component row and column vectors.
    template<typename T> struct Components<T, 3, 1>: public NonScalarType
    {
        union {
            struct { T x, y, z; };
            struct { T r, g, b; };
            T m_data[3];
        };

        Components() {}
        template<typename T1> explicit Components(T1 _s) : x(static_cast<T>(_s)), y(static_cast<T>(_s)), z(static_cast<T>(_s)) {}
        template<typename T1, typename T2, typename T3>
        Components(T1 _s0, T2 _s1, T3 _s2) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1)), z(static_cast<T>(_s2)) {}
        template<typename T1, typename T2>
        Components(T1 _s0, const Components<T2, 2, 1>& _s1) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1.x)), z(static_cast<T>(_s1.y)) {}
        template<typename T1, typename T2>
        Components(const Components<T1, 2, 1>& _s0, T2 _s1) : x(static_cast<T>(_s0.x)), y(static_cast<T>(_s0.y)), z(static_cast<T>(_s1)) {}
    };
    template<typename T> struct Components<T, 1, 3>: public NonScalarType
    {
        union {
            struct { T x, y, z; };
            struct { T r, g, b; };
            T m_data[3];
        };

        Components() {}
        template<typename T1> explicit Components(T1 _s) : x(static_cast<T>(_s)), y(static_cast<T>(_s)), z(static_cast<T>(_s)) {}
        template<typename T1, typename T2, typename T3>
        Components(T1 _s0, T2 _s1, T3 _s2) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1)), z(static_cast<T>(_s2)) {}
        template<typename T1, typename T2>
        Components(T1 _s0, const Components<T2, 1, 2>& _s1) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1.x)), z(static_cast<T>(_s1.y)) {}
        template<typename T1, typename T2>
        Components(const Components<T1, 1, 2>& _s0, T2 _s1) : x(static_cast<T>(_s0.x)), y(static_cast<T>(_s0.y)), z(static_cast<T>(_s1)) {}
    };

    /// \brief Specialized version for 4 component row and column vectors and
    ///      2x2 matrices.
    template <typename T> struct Components<T, 4, 1>: public NonScalarType
    {
        union {
            struct { T x, y, z, w; };
            struct { T r, g, b, a; };
            T m_data[4];
        };

        Components() {}
        template<typename T1> explicit Components(T1 _s) : x(static_cast<T>(_s)), y(static_cast<T>(_s)), z(static_cast<T>(_s)), w(static_cast<T>(_s)) {}
        template<typename T1, typename T2, typename T3, typename T4>
        Components(T1 _s0, T2 _s1, T3 _s2, T4 _s3) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1)), z(static_cast<T>(_s2)), w(static_cast<T>(_s3)) {}
        template<typename T1>
        Components(T _s0, const Components<T1, 3, 1>& _s1) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1.x)), z(static_cast<T>(_s1.y)), w(static_cast<T>(_s1.z)) {}
        template<typename T1>
        Components(const Components<T1, 3, 1>& _s0, T _s1) : x(static_cast<T>(_s0.x)), y(static_cast<T>(_s0.y)), z(static_cast<T>(_s0.z)), w(static_cast<T>(_s1)) {}
        template<typename T1>
        Components(const Components<T1, 2, 1>& _s0, T _s1, T _s2) : x(static_cast<T>(_s0.x)), y(static_cast<T>(_s0.y)), z(static_cast<T>(_s1)), w(static_cast<T>(_s2)) {}
        template<typename T1>
        Components(T _s0, const Components<T1, 2, 1>& _s1, T _s2) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1.x)), z(static_cast<T>(_s1.y)), w(static_cast<T>(_s2)) {}
        template<typename T1>
        Components(T _s0, T _s1, const Components<T1, 2, 1>& _s2) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1)), z(static_cast<T>(_s2.x)), w(static_cast<T>(_s2.y)) {}
        template<typename T1, typename T2>
        Components(const Components<T1, 2, 1>& _s0, const Components<T2, 2, 1>& _s1) : x(static_cast<T>(_s0.x)), y(static_cast<T>(_s0.y)), z(static_cast<T>(_s1.x)), w(static_cast<T>(_s1.y)) {}
    };
    template <typename T> struct Components<T, 1, 4>: public NonScalarType
    {
        union {
            struct { T x, y, z, w; };
            struct { T r, g, b, a; };
            T m_data[4];
        };

        Components() {}
        template<typename T1> explicit Components(T1 _s) : x(static_cast<T>(_s)), y(static_cast<T>(_s)), z(static_cast<T>(_s)), w(static_cast<T>(_s)) {}
        template<typename T1, typename T2, typename T3, typename T4>
        Components(T1 _s0, T2 _s1, T3 _s2, T4 _s3) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1)), z(static_cast<T>(_s2)), w(static_cast<T>(_s3)) {}
        template<typename T1>
        Components(T _s0, const Components<T1, 1, 3>& _s1) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1.x)), z(static_cast<T>(_s1.y)), w(static_cast<T>(_s1.z)) {}
        template<typename T1>
        Components(const Components<T1, 1, 3>& _s0, T _s1) : x(static_cast<T>(_s0.x)), y(static_cast<T>(_s0.y)), z(static_cast<T>(_s0.z)), w(static_cast<T>(_s1)) {}
        template<typename T1>
        Components(const Components<T1, 1, 2>& _s0, T _s1, T _s2) : x(static_cast<T>(_s0.x)), y(static_cast<T>(_s0.y)), z(static_cast<T>(_s1)), w(static_cast<T>(_s2)) {}
        template<typename T1>
        Components(T _s0, const Components<T1, 1, 2>& _s1, T _s2) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1.x)), z(static_cast<T>(_s1.y)), w(static_cast<T>(_s2)) {}
        template<typename T1>
        Components(T _s0, T _s1, const Components<T1, 1, 2>& _s2) : x(static_cast<T>(_s0)), y(static_cast<T>(_s1)), z(static_cast<T>(_s2.x)), w(static_cast<T>(_s2.y)) {}
        template<typename T1, typename T2>
        Components(const Components<T1, 1, 2>& _s0, const Components<T2, 1, 2>& _s1) : x(static_cast<T>(_s0.x)), y(static_cast<T>(_s0.y)), z(static_cast<T>(_s1.x)), w(static_cast<T>(_s1.y)) {}
    };
    template <typename T> struct Components<T, 2, 2>: public NonScalarType
    {
        union {
            struct { T m00, m01,
                       m10, m11;
            };
            T m_data[4];
        };

        Components() {}
        template<typename T1> explicit Components(T1 _s) :
            m00(static_cast<T>(_s)), m01(static_cast<T>(_s)),
            m10(static_cast<T>(_s)), m11(static_cast<T>(_s))
        {}
        template<typename T1, typename T2, typename T3, typename T4>
        Components(T1 _s0, T2 _s1, T3 _s2, T4 _s3) :
            m00(static_cast<T>(_s0)), m01(static_cast<T>(_s1)),
            m10(static_cast<T>(_s2)), m11(static_cast<T>(_s3))
        {}
        template<typename T1>
        Components(const Components<T1, 2, 1>& _s0, const Components<T1, 2, 1>& _s1) :
            m00(static_cast<T>(_s0.x)), m01(static_cast<T>(_s1.x)),
            m10(static_cast<T>(_s0.y)), m11(static_cast<T>(_s1.y))
        {}
        template<typename T1>
        Components(const Components<T1, 1, 2>& _s0, const Components<T1, 1, 2>& _s1) :
            m00(static_cast<T>(_s0.x)), m01(static_cast<T>(_s0.y)),
            m10(static_cast<T>(_s1.x)), m11(static_cast<T>(_s1.y))
        {}
    };

    /// \brief Specialized version for 3x3 matrices.
    template <typename T> struct Components<T, 3, 3>: public NonScalarType
    {
        union {
            struct { T m00, m01, m02,
                       m10, m11, m12,
                       m20, m21, m22;
            };
            T m_data[9];
        };

        Components() {}
        template<typename T1> explicit Components(T1 _s) :
            m00(static_cast<T>(_s)), m01(static_cast<T>(_s)), m02(static_cast<T>(_s)),
            m10(static_cast<T>(_s)), m11(static_cast<T>(_s)), m12(static_cast<T>(_s)),
            m20(static_cast<T>(_s)), m21(static_cast<T>(_s)), m22(static_cast<T>(_s))
        {}
        template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
        Components(T1 _s0, T2 _s1, T3 _s2, T4 _s3, T5 _s4, T6 _s5, T7 _s6, T8 _s7, T9 _s8) :
            m00(static_cast<T>(_s0)), m01(static_cast<T>(_s1)), m02(static_cast<T>(_s2)),
            m10(static_cast<T>(_s3)), m11(static_cast<T>(_s4)), m12(static_cast<T>(_s5)),
            m20(static_cast<T>(_s6)), m21(static_cast<T>(_s7)), m22(static_cast<T>(_s8))
        {}
        template<typename T1>
        Components(const Components<T1, 3, 1>& _s0, const Components<T1, 3, 1>& _s1, const Components<T1, 3, 1>& _s2) :
            m00(static_cast<T>(_s0.x)), m01(static_cast<T>(_s1.x)), m02(static_cast<T>(_s2.x)),
            m10(static_cast<T>(_s0.y)), m11(static_cast<T>(_s1.y)), m12(static_cast<T>(_s2.y)),
            m20(static_cast<T>(_s0.z)), m21(static_cast<T>(_s1.z)), m22(static_cast<T>(_s2.z))
        {}
        template<typename T1>
        Components(const Components<T1, 1, 3>& _s0, const Components<T1, 1, 3>& _s1, const Components<T1, 1, 3>& _s2) :
            m00(static_cast<T>(_s0.x)), m01(static_cast<T>(_s0.y)), m02(static_cast<T>(_s0.z)),
            m10(static_cast<T>(_s1.x)), m11(static_cast<T>(_s1.y)), m12(static_cast<T>(_s1.z)),
            m20(static_cast<T>(_s2.x)), m21(static_cast<T>(_s2.y)), m22(static_cast<T>(_s2.z))
        {}
        explicit Components(const ei::TQuaternion<T>& _quaternion)
        {
            // Rotation composition from quaternion (remaining rest direct in matrix)
            // See http://de.wikipedia.org/wiki/Quaternion#Bezug_zu_orthogonalen_Matrizen for
            // details.
            T f2i  = 2.0f * _quaternion.i;
            T f2j  = 2.0f * _quaternion.j;
            T f2k  = 2.0f * _quaternion.k;
            T f2ri = f2i  * _quaternion.r;
            T f2rj = f2j  * _quaternion.r;
            T f2rk = f2k  * _quaternion.r;
            T f2ii = f2i  * _quaternion.i;
            T f2ij = f2j  * _quaternion.i;
            T f2ik = f2k  * _quaternion.i;
            T f2jj = f2j  * _quaternion.j;
            T f2jk = f2k  * _quaternion.j;
            T f2kk = f2k  * _quaternion.k;

            this->m00 = 1.0f - ( f2jj + f2kk ); this->m01 = f2ij - f2rk;            this->m02 = f2ik + f2rj;
            this->m10 = f2ij + f2rk;            this->m11 = 1.0f - ( f2ii + f2kk ); this->m12 = f2jk - f2ri;
            this->m20 = f2ik - f2rj;            this->m21 = f2jk + f2ri;            this->m22 = 1.0f - ( f2ii + f2jj );
        }
    };

    /// \brief Specialized version for 4x4 matrices.
    template <typename T> struct Components<T, 4, 4>: public NonScalarType
    {
        union {
            struct { T m00, m01, m02, m03,
                       m10, m11, m12, m13,
                       m20, m21, m22, m23,
                       m30, m31, m32, m33;
            };
            T m_data[16];
        };

        Components() {}
        template<typename T1> explicit Components(T1 _s) :
            m00(static_cast<T>(_s)), m01(static_cast<T>(_s)), m02(static_cast<T>(_s)), m03(static_cast<T>(_s)),
            m10(static_cast<T>(_s)), m11(static_cast<T>(_s)), m12(static_cast<T>(_s)), m13(static_cast<T>(_s)),
            m20(static_cast<T>(_s)), m21(static_cast<T>(_s)), m22(static_cast<T>(_s)), m23(static_cast<T>(_s)),
            m30(static_cast<T>(_s)), m31(static_cast<T>(_s)), m32(static_cast<T>(_s)), m33(static_cast<T>(_s))
        {}
        Components(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7, T _s8, T _s9, T _s10, T _s11, T _s12, T _s13, T _s14, T _s15) :
            m00(static_cast<T>(_s0)), m01(static_cast<T>(_s1)), m02(static_cast<T>(_s2)), m03(static_cast<T>(_s3)),
            m10(static_cast<T>(_s4)), m11(static_cast<T>(_s5)), m12(static_cast<T>(_s6)), m13(static_cast<T>(_s7)),
            m20(static_cast<T>(_s8)), m21(static_cast<T>(_s9)), m22(static_cast<T>(_s10)), m23(static_cast<T>(_s11)),
            m30(static_cast<T>(_s12)), m31(static_cast<T>(_s13)), m32(static_cast<T>(_s14)), m33(static_cast<T>(_s15))
        {}
        template<typename T1>
        Components(const Components<T1, 4, 1>& _s0, const Components<T1, 4, 1>& _s1, const Components<T1, 4, 1>& _s2, const Components<T1, 4, 1>& _s3) :
            m00(static_cast<T>(_s0.x)), m01(static_cast<T>(_s1.x)), m02(static_cast<T>(_s2.x)), m03(static_cast<T>(_s3.x)),
            m10(static_cast<T>(_s0.y)), m11(static_cast<T>(_s1.y)), m12(static_cast<T>(_s2.y)), m13(static_cast<T>(_s3.y)),
            m20(static_cast<T>(_s0.z)), m21(static_cast<T>(_s1.z)), m22(static_cast<T>(_s2.z)), m23(static_cast<T>(_s3.z)),
            m30(static_cast<T>(_s0.w)), m31(static_cast<T>(_s1.w)), m32(static_cast<T>(_s2.w)), m33(static_cast<T>(_s3.w))
        {}
        template<typename T1>
        Components(const Components<T1, 1, 4>& _s0, const Components<T1, 1, 4>& _s1, const Components<T1, 1, 4>& _s2, const Components<T1, 1, 4>& _s3) :
            m00(static_cast<T>(_s0.x)), m01(static_cast<T>(_s0.y)), m02(static_cast<T>(_s0.z)), m03(static_cast<T>(_s0.w)),
            m10(static_cast<T>(_s1.x)), m11(static_cast<T>(_s1.y)), m12(static_cast<T>(_s1.z)), m13(static_cast<T>(_s1.w)),
            m20(static_cast<T>(_s2.x)), m21(static_cast<T>(_s2.y)), m22(static_cast<T>(_s2.z)), m23(static_cast<T>(_s2.w)),
            m30(static_cast<T>(_s3.x)), m31(static_cast<T>(_s3.y)), m32(static_cast<T>(_s3.z)), m33(static_cast<T>(_s3.w))
        {}
    };

}


// ************************************************************************** //
// Code generators for highly redundant operator implementations.             //
// ************************************************************************** //
#define EI_CODE_GEN_MAT_MAT_OP(op) \
{ \
    Matrix<RESULT_TYPE(op), M, N> result; \
    for(uint i = 0; i < N * M; ++i) \
        result[i] = (*this)[i] op _mat1[i]; \
    return result; \
}

#define EI_CODE_GEN_MAT_UNARY_OP(op) \
{ \
    Matrix<T, M, N> result; \
    for(uint i = 0; i < N * M; ++i) \
        result[i] = op (*this)[i]; \
    return result; \
}

#define EI_CODE_GEN_MAT_MAT_SEFL_OP(op) \
{ \
    for(uint i = 0; i < N * M; ++i) \
        (*this)[i] op _mat1[i]; \
    return *this; \
}

#define EI_CODE_GEN_MAT_SCALAR_SEFL_OP(op) \
{ \
    for(uint i = 0; i < N * M; ++i) \
        (*this)[i] op _s; \
    return *this; \
}

#define EI_CODE_GEN_MAT_MAT_BOOL_OP(op) \
{ \
    Matrix<bool, M, N> result; \
    for(uint i = 0; i < N * M; ++i) \
        result[i] = (*this)[i] op _mat1[i]; \
    return result; \
}

// ************************************************************************** //
#define EI_CODE_GEN_MAT_SCALAR_OP(op) \
{ \
    Matrix<RESULT_TYPE(op), M, N> result; \
    for(uint i = 0; i < N * M; ++i) \
        result[i] = _mat[i] op _s; \
    return result; \
}

#define EI_CODE_GEN_SCALAR_MAT_OP(op) \
{ \
    Matrix<RESULT_TYPE(op), M, N> result; \
    for(uint i = 0; i < N * M; ++i) \
        result[i] = _s op _mat[i]; \
    return result; \
}

// ************************************************************************** //
#define EI_CODE_GEN_MAT_SCALAR_BOOL_OP(op) \
{ \
    Matrix<bool, M, N> result; \
    for(uint i = 0; i < N * M; ++i) \
        result[i] = _mat[i] op _s; \
    return result; \
}

#define EI_CODE_GEN_SCALAR_MAT_BOOL_OP(op) \
{ \
    Matrix<bool, M, N> result; \
    for(uint i = 0; i < N * M; ++i) \
        result[i] = _s op _mat[i]; \
    return result; \
}
