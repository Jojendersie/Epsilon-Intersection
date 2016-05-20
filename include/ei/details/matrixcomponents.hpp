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
        template<typename... Args>
        Components(Args... _args) : m_data{ T(_args)... }
        {
            static_assert(sizeof...(Args) == M*N, "Wrong number of arguments!");
        }
    };


    /// \brief Specialized version for 1 component row or column vector.
    template<typename T> struct Components<T, 1, 1>: public NonScalarType
    {
        union {
            T x;
            T r;
            T m_data[1];
        };
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
        template<typename T1, typename T2>
        Components(T1 _s0, T2 _s1) : x(_s0), y(_s1) {}
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
        template<typename T1, typename T2>
        Components(T1 _s0, T2 _s1) : x(_s0), y(_s1) {}
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
        template<typename T1, typename T2, typename T3>
        Components(T1 _s0, T2 _s1, T3 _s2) : x(_s0), y(_s1), z(_s2) {}
        template<typename T1, typename T2>
        Components(T1 _s0, const Components<T2, 2, 1>& _s1) : x(_s0), y(_s1.x), z(_s1.y) {}
        template<typename T1, typename T2>
        Components(const Components<T1, 2, 1>& _s0, T2 _s1) : x(_s0.x), y(_s0.y), z(_s1) {}
    };
    template<typename T> struct Components<T, 1, 3>: public NonScalarType
    {
        union {
            struct { T x, y, z; };
            struct { T r, g, b; };
            T m_data[3];
        };

        Components() {}
        template<typename T1, typename T2, typename T3>
        Components(T1 _s0, T2 _s1, T3 _s2) : x(_s0), y(_s1), z(_s2) {}
        template<typename T1>
        Components(T _s0, const Components<T1, 1, 2>& _s1) : x(_s0), y(_s1.x), z(_s1.y) {}
        template<typename T1>
        Components(const Components<T1, 1, 2>& _s0, T _s1) : x(_s0.x), y(_s0.y), z(_s1) {}
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
        template<typename T1, typename T2, typename T3, typename T4>
        Components(T1 _s0, T2 _s1, T3 _s2, T4 _s3) : x(_s0), y(_s1), z(_s2), w(_s3) {}
        template<typename T1>
        Components(T _s0, const Components<T1, 3, 1>& _s1) : x(_s0), y(_s1.x), z(_s1.y), w(_s1.z) {}
        template<typename T1>
        Components(const Components<T1, 3, 1>& _s0, T _s1) : x(_s0.x), y(_s0.y), z(_s0.z), w(_s1) {}
        template<typename T1>
        Components(const Components<T1, 2, 1>& _s0, T _s1, T _s2) : x(_s0.x), y(_s0.y), z(_s1), w(_s2) {}
        template<typename T1>
        Components(T _s0, const Components<T1, 2, 1>& _s1, T _s2) : x(_s0), y(_s1.x), z(_s1.y), w(_s2) {}
        template<typename T1>
        Components(T _s0, T _s1, const Components<T1, 2, 1>& _s2) : x(_s0), y(_s1), z(_s2.x), w(_s2.y) {}
        template<typename T1, typename T2>
        Components(const Components<T1, 2, 1>& _s0, const Components<T2, 2, 1>& _s1) : x(_s0.x), y(_s0.y), z(_s1.x), w(_s1.y) {}
    };
    template <typename T> struct Components<T, 1, 4>: public NonScalarType
    {
        union {
            struct { T x, y, z, w; };
            struct { T r, g, b, a; };
            T m_data[4];
        };

        Components() {}
        template<typename T1, typename T2, typename T3, typename T4>
        Components(T1 _s0, T2 _s1, T3 _s2, T4 _s3) : x(_s0), y(_s1), z(_s2), w(_s3) {}
        template<typename T1>
        Components(T _s0, const Components<T1, 1, 3>& _s1) : x(_s0), y(_s1.x), z(_s1.y), w(_s1.z) {}
        template<typename T1>
        Components(const Components<T1, 1, 3>& _s0, T _s1) : x(_s0.x), y(_s0.y), z(_s0.z), w(_s1) {}
        template<typename T1>
        Components(const Components<T1, 1, 2>& _s0, T _s1, T _s2) : x(_s0.x), y(_s0.y), z(_s1), w(_s2) {}
        template<typename T1>
        Components(T _s0, const Components<T1, 1, 2>& _s1, T _s2) : x(_s0), y(_s1.x), z(_s1.y), w(_s2) {}
        template<typename T1>
        Components(T _s0, T _s1, const Components<T1, 1, 2>& _s2) : x(_s0), y(_s1), z(_s2.x), w(_s2.y) {}
        template<typename T1, typename T2>
        Components(const Components<T1, 1, 2>& _s0, const Components<T2, 1, 2>& _s1) : x(_s0.x), y(_s0.y), z(_s1.x), w(_s1.y) {}
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
        template<typename T1, typename T2, typename T3, typename T4>
        Components(T1 _s0, T2 _s1, T3 _s2, T4 _s3) :
            m00(_s0), m01(_s1),
            m10(_s2), m11(_s3)
        {}
        template<typename T1>
        Components(const Components<T1, 2, 1>& _s0, const Components<T1, 2, 1>& _s1) :
            m00(_s0.x), m01(_s1.x),
            m10(_s0.y), m11(_s1.y)
        {}
        template<typename T1>
        Components(const Components<T1, 1, 2>& _s0, const Components<T1, 1, 2>& _s1) :
            m00(_s0.x), m01(_s0.y),
            m10(_s1.x), m11(_s1.y)
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
        template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
        Components(T1 _s0, T2 _s1, T3 _s2, T4 _s3, T5 _s4, T6 _s5, T7 _s6, T8 _s7, T9 _s8) :
            m00(_s0), m01(_s1), m02(_s2),
            m10(_s3), m11(_s4), m12(_s5),
            m20(_s6), m21(_s7), m22(_s8)
        {}
        template<typename T1>
        Components(const Components<T1, 3, 1>& _s0, const Components<T1, 3, 1>& _s1, const Components<T1, 3, 1>& _s2) :
            m00(_s0.x), m01(_s1.x), m02(_s2.x),
            m10(_s0.y), m11(_s1.y), m12(_s2.y),
            m20(_s0.z), m21(_s1.z), m22(_s2.z)
        {}
        template<typename T1>
        Components(const Components<T1, 1, 3>& _s0, const Components<T1, 1, 3>& _s1, const Components<T1, 1, 3>& _s2) :
            m00(_s0.x), m01(_s0.y), m02(_s0.z),
            m10(_s1.x), m11(_s1.y), m12(_s1.z),
            m20(_s2.x), m21(_s2.y), m22(_s2.z)
        {}
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
        Components(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7, T _s8, T _s9, T _s10, T _s11, T _s12, T _s13, T _s14, T _s15) :
            m00(_s0), m01(_s1), m02(_s2), m03(_s3),
            m10(_s4), m11(_s5), m12(_s6), m13(_s7),
            m20(_s8), m21(_s9), m22(_s10), m23(_s11),
            m30(_s12), m31(_s13), m32(_s14), m33(_s15)
        {}
        template<typename T1>
        Components(const Components<T1, 4, 1>& _s0, const Components<T1, 4, 1>& _s1, const Components<T1, 4, 1>& _s2, const Components<T1, 4, 1>& _s3) :
            m00(_s0.x), m01(_s1.x), m02(_s2.x), m03(_s3.x),
            m10(_s0.y), m11(_s1.y), m12(_s2.y), m13(_s3.y),
            m20(_s0.z), m21(_s1.z), m22(_s2.z), m23(_s3.z),
            m30(_s0.w), m31(_s1.w), m32(_s2.w), m33(_s3.w)
        {}
        template<typename T1>
        Components(const Components<T1, 1, 4>& _s0, const Components<T1, 1, 4>& _s1, const Components<T1, 1, 4>& _s2, const Components<T1, 1, 4>& _s3) :
            m00(_s0.x), m01(_s0.y), m02(_s0.z), m03(_s0.w),
            m10(_s1.x), m11(_s1.y), m12(_s1.z), m13(_s1.w),
            m20(_s2.x), m21(_s2.y), m22(_s2.z), m23(_s2.w),
            m30(_s3.x), m31(_s3.y), m32(_s3.z), m33(_s3.w)
        {}
    };

}
