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
    template<typename T, unsigned N> using Vec  = Matrix<T, N, 1>;
    /// \brief General row-vector type
    template<typename T, unsigned N> using RVec = Matrix<T, 1, N>;

    /// \brief 2D column-vector of type float.
    typedef Matrix<float, 2, 1> Vec2;
    /// \brief 3D column-vector of type float.
    typedef Matrix<float, 3, 1> Vec3;
    /// \brief 4D column-vector of type float.
    typedef Matrix<float, 4, 1> Vec4;

    /// \brief 2D row-vector of type float.
    typedef Matrix<float, 1, 2> RVec2;
    /// \brief 3D row-vector of type float.
    typedef Matrix<float, 1, 3> RVec3;
    /// \brief 4D row-vector of type float.
    typedef Matrix<float, 1, 4> RVec4;

    /// \brief 2x2 matrix of type float.
    typedef Matrix<float, 2, 2> Mat2x2;
    /// \brief 3x3 matrix of type float.
    typedef Matrix<float, 3, 3> Mat3x3;
    /// \brief 4x4 matrix of type float.
    typedef Matrix<float, 4, 4> Mat4x4;

    // ************************************************************************** //
    // Predefined double vector and matrix types.

    /// \brief 2D column-vector of type double.
    typedef Matrix<double, 2, 1> DVec2;
    /// \brief 3D column-vector of type double.
    typedef Matrix<double, 3, 1> DVec3;
    /// \brief 4D column-vector of type double.
    typedef Matrix<double, 4, 1> DVec4;

    /// \brief 2D row-vector of type double.
    typedef Matrix<double, 1, 2> DRVec2;
    /// \brief 3D row-vector of type double.
    typedef Matrix<double, 1, 3> DRVec3;
    /// \brief 4D row-vector of type double.
    typedef Matrix<double, 1, 4> DRVec4;

    /// \brief 2x2 matrix of type double.
    typedef Matrix<double, 2, 2> DMat2x2;
    /// \brief 3x3 matrix of type double.
    typedef Matrix<double, 3, 3> DMat3x3;
    /// \brief 4x4 matrix of type double.
    typedef Matrix<double, 4, 4> DMat4x4;

    // ************************************************************************** //
    // Predefined 32 bit integer vector and matrix types.

    /// \brief 2D column-vector of type int32.
    typedef Matrix<int32, 2, 1> IVec2;
    /// \brief 3D column-vector of type int32.
    typedef Matrix<int32, 3, 1> IVec3;
    /// \brief 4D column-vector of type int32.
    typedef Matrix<int32, 4, 1> IVec4;

    /// \brief 2D row-vector of type int32.
    typedef Matrix<int32, 1, 2> IRVec2;
    /// \brief 3D row-vector of type int32.
    typedef Matrix<int32, 1, 3> IRVec3;
    /// \brief 4D row-vector of type int32.
    typedef Matrix<int32, 1, 4> IRVec4;

    /// \brief 2x2 matrix of type int32.
    typedef Matrix<int32, 2, 2> IMat2x2;
    /// \brief 3x3 matrix of type int32.
    typedef Matrix<int32, 3, 3> IMat3x3;
    /// \brief 4x4 matrix of type int32.
    typedef Matrix<int32, 4, 4> IMat4x4;

    // ************************************************************************** //
    // Predefined 32 bit unsigned integer vector and matrix types.

    /// \brief 2D column-vector of type uint32.
    typedef Matrix<uint32, 2, 1> UVec2;
    /// \brief 3D column-vector of type uint32.
    typedef Matrix<uint32, 3, 1> UVec3;
    /// \brief 4D column-vector of type uint32.
    typedef Matrix<uint32, 4, 1> UVec4;

    /// \brief 2D row-vector of type uint32.
    typedef Matrix<uint32, 1, 2> URVec2;
    /// \brief 3D row-vector of type uint32.
    typedef Matrix<uint32, 1, 3> URVec3;
    /// \brief 4D row-vector of type uint32.
    typedef Matrix<uint32, 1, 4> URVec4;

    /// \brief 2x2 matrix of type uint32.
    typedef Matrix<uint32, 2, 2> UMat2x2;
    /// \brief 3x3 matrix of type uint32.
    typedef Matrix<uint32, 3, 3> UMat3x3;
    /// \brief 4x4 matrix of type uint32.
    typedef Matrix<uint32, 4, 4> UMat4x4;

    // ********************************************************************* //
    typedef TQuaternion<float> Quaternion;
    typedef TQuaternion<double> DQuaternion;
}
