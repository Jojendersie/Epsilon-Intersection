namespace ei {

    // ************************************************************************** //
    // Basic template classes
    template<typename T>
    class TQuaternion;

    template<typename T, uint M, uint N>
    class Matrix;

    // ************************************************************************** //
    // Predefined float vector and matrix types.
    /// \brief General column-vector type
    template<typename T, uint N> using Vec  = Matrix<T, N, 1u>;
    /// \brief General row-vector type
    template<typename T, uint N> using RVec = Matrix<T, 1u, N>;

    /// \brief 2D column-vector of type float.
    typedef Matrix<float, 2u, 1u> Vec2;
    /// \brief 3D column-vector of type float.
    typedef Matrix<float, 3u, 1u> Vec3;
    /// \brief 4D column-vector of type float.
    typedef Matrix<float, 4u, 1u> Vec4;

    /// \brief 2x2 matrix of type float.
    typedef Matrix<float, 2u, 2u> Mat2x2;
    /// \brief 3x3 matrix of type float.
    typedef Matrix<float, 3u, 3u> Mat3x3;
    /// \brief 3x4 (3 rows) matrix of type float.
    typedef Matrix<float, 3u, 4u> Mat3x4;
    /// \brief 4x3 (3 columns) matrix of type float.
    typedef Matrix<float, 4u, 3u> Mat4x3;
    /// \brief 4x4 matrix of type float.
    typedef Matrix<float, 4u, 4u> Mat4x4;

    // ************************************************************************** //
    // Predefined double vector and matrix types.

    /// \brief 2D column-vector of type double.
    typedef Matrix<double, 2u, 1u> DVec2;
    /// \brief 3D column-vector of type double.
    typedef Matrix<double, 3u, 1u> DVec3;
    /// \brief 4D column-vector of type double.
    typedef Matrix<double, 4u, 1u> DVec4;

    /// \brief 2x2 matrix of type double.
    typedef Matrix<double, 2u, 2u> DMat2x2;
    /// \brief 3x3 matrix of type double.
    typedef Matrix<double, 3u, 3u> DMat3x3;
    /// \brief 3x4 (3 rows) matrix of type float.
    typedef Matrix<double, 3u, 4u> DMat3x4;
    /// \brief 4x3 (3 columns) matrix of type float.
    typedef Matrix<double, 4u, 3u> DMat4x3;
    /// \brief 4x4 matrix of type double.
    typedef Matrix<double, 4u, 4u> DMat4x4;

    // ************************************************************************** //
    // Predefined 32 bit integer vector and matrix types.

    /// \brief 2D column-vector of type int32.
    typedef Matrix<int32, 2u, 1u> IVec2;
    /// \brief 3D column-vector of type int32.
    typedef Matrix<int32, 3u, 1u> IVec3;
    /// \brief 4D column-vector of type int32.
    typedef Matrix<int32, 4u, 1u> IVec4;

    /// \brief 2x2 matrix of type int32.
    typedef Matrix<int32, 2u, 2u> IMat2x2;
    /// \brief 3x3 matrix of type int32.
    typedef Matrix<int32, 3u, 3u> IMat3x3;
    /// \brief 4x4 matrix of type int32.
    typedef Matrix<int32, 4u, 4u> IMat4x4;

    // ************************************************************************** //
    // Predefined 32 bit unsigned integer vector and matrix types.

    /// \brief 2D column-vector of type uint32.
    typedef Matrix<uint32, 2u, 1u> UVec2;
    /// \brief 3D column-vector of type uint32.
    typedef Matrix<uint32, 3u, 1u> UVec3;
    /// \brief 4D column-vector of type uint32.
    typedef Matrix<uint32, 4u, 1u> UVec4;

    /// \brief 2x2 matrix of type uint32.
    typedef Matrix<uint32, 2u, 2u> UMat2x2;
    /// \brief 3x3 matrix of type uint32.
    typedef Matrix<uint32, 3u, 3u> UMat3x3;
    /// \brief 4x4 matrix of type uint32.
    typedef Matrix<uint32, 4u, 4u> UMat4x4;

    // ********************************************************************* //
    typedef TQuaternion<float> Quaternion;
    typedef TQuaternion<double> DQuaternion;
}
