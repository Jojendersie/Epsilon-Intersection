namespace ei { namespace details {
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
