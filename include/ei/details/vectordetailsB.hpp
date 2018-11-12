namespace ei { namespace details {
    // ********************************************************************* //
    // Matrix constants
    static const ei::Matrix<float,2,2> MAT2X2_IDENTITY =
        ei::Mat2x2(1.0f, 0.0f,
            0.0f, 1.0f);

    static const ei::Matrix<float,3,3> MAT3X3_IDENTITY =
        ei::Mat3x3(1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f);

    static const ei::Matrix<float,4,4> MAT4X4_IDENTITY =
        ei::Mat4x4(1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f);

    // Quaternion constants
    static const ei::TQuaternion<float> QUATERNION_IDENTITY =
        ei::TQuaternion<float>(0.0f, 0.0f, 0.0f, 1.0f);
    static const ei::TQuaternion<double> QUATERNIOND_IDENTITY =
        ei::TQuaternion<double>(0.0, 0.0, 0.0, 1.0);

    // Recursive helper for orthonormalization
    // Recursion end
    template<typename TVec0>
    inline void removeProjectedPart(const TVec0&)
    {
    }
    template<typename TVec0, typename TVec1, typename... TVecs>
    inline void removeProjectedPart(const TVec0& _vec0, TVec1& _vec1, TVecs&... _vecs)
    {
        _vec1 -= dot(_vec0, _vec1) * _vec0;
        removeProjectedPart(_vec0, _vecs...);
    }
}} // namespace ei::details

// Remove helper macros from vectordetailsA.hpp.
#undef EI_CODE_GEN_MAT_MAT_OP
#undef EI_CODE_GEN_MAT_UNARY_OP
#undef EI_CODE_GEN_MAT_MAT_SEFL_OP
#undef EI_CODE_GEN_MAT_SCALAR_SEFL_OP
#undef EI_CODE_GEN_MAT_MAT_BOOL_OP
#undef EI_CODE_GEN_MAT_SCALAR_OP
#undef EI_CODE_GEN_SCALAR_MAT_OP
#undef EI_CODE_GEN_MAT_SCALAR_BOOL_OP
#undef EI_CODE_GEN_SCALAR_MAT_BOOL_OP
