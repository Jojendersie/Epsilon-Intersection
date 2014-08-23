namespace details {

	/// \brief Dummy class to detect correct types for matrix <-> matrix
	///     and matrix <-> scalar operations.
	/// \details The overloading mechanism fails when both types of operations
	///     are templated, because the matrix <-> scalar is chosen even if
	///     both operants are matrices.
	class MatrixType	{};

	/// \brief Specialized component access for small vectors and matrices.
	/// \details This is the fallback for larger vectors without single component
	///		access.
	template<typename T, unsigned M, unsigned N> struct Components: public MatrixType
	{
	protected:
		T m_data[M * N];
	};

	/// \brief Specialized version for 2 component row and column vectors.
	template<typename T> struct Components<T, 2, 1>: public MatrixType
	{
		union {
			struct { T x, y; };
			struct { T r, g; };
			T m_data[2];
		};
	};
	template<typename T> struct Components<T, 1, 2>: public MatrixType
	{
		union {
			struct { T x, y; };
			struct { T r, g; };
			T m_data[2];
		};
	};

	/// \brief Specialized version for 3 component row and column vectors.
	template<typename T> struct Components<T, 3, 1>: public MatrixType
	{
		union {
			struct { T x, y, z; };
			struct { T r, g, b; };
			T m_data[3];
		};
	};
	template<typename T> struct Components<T, 1, 3>: public MatrixType
	{
		union {
			struct { T x, y, z; };
			struct { T r, g, b; };
			T m_data[3];
		};
	};

	/// \brief Specialized version for 4 component row and column vectors and
	///      2x2 matrices.
	template <typename T> struct Components<T, 4, 1>: public MatrixType
	{
		union {
			struct { T x, y, z, w; };
			struct { T r, g, b, a; };
			T m_data[4];
		};
	};
	template <typename T> struct Components<T, 1, 4>: public MatrixType
	{
		union {
			struct { T x, y, z, w; };
			struct { T r, g, b, a; };
			T m_data[4];
		};
	};
	template <typename T> struct Components<T, 2, 2>: public MatrixType
	{
		union {
			struct { T m00, m01,
				       m10, m11;
			};
			T m_data[4];
		};
	};

	/// \brief Specialized version for 3x3 matrices.
	template <typename T> struct Components<T, 3, 3>: public MatrixType
	{
		union {
			struct { T m00, m01, m02,
				       m10, m11, m12,
				       m20, m21, m22;
			};
			T m_data[9];
		};
	};

	/// \brief Specialized version for 4x4 matrices.
	template <typename T> struct Components<T, 4, 4>: public MatrixType
	{
		union {
			struct { T m00, m01, m02, m03,
			           m10, m11, m12, m13,
				       m20, m21, m22, m23,
				       m30, m31, m32, m33;
			};
			T m_data[16];
		};
	};
}
