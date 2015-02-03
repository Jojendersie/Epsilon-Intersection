// ************************************************************************* //
//                               CONSTRUCTORS                                //
// ************************************************************************* //

// ************************************************************************* //
template<typename T, uint M, uint N>
Matrix<T, M, N>::Matrix()
{
}

// ************************************************************************* //
template<typename T, uint M, uint N>
Matrix<T, M, N>::Matrix(T _s)
{
    for(uint i = 0; i < N * M; ++i)
        m_data[i] = _s;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1>
Matrix<T, M, N>::Matrix(const Matrix<T1,M,N>& _mat1)
{
    for(uint i = 0; i < N * M; ++i)
        m_data[i] = static_cast<T>(_mat1[i]);
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, uint M1, uint N1, class>
Matrix<T, M, N>::Matrix(const Matrix<T1,M1,N1>& _mat1, uint _rowOff = 0, uint _colOff = 0)
{
    eiAssert( _rowOff + M <= M1, "Out of boundaries: matrix subsection wrong!" );
    eiAssert( _colOff + N <= N1, "Out of boundaries: matrix subsection wrong!" );
    // This counter avoids one index computation y*N+x in the inner loop
    uint i = 0;
    for(uint y = _rowOff; y < M+_rowOff; ++y)
        for(uint x = _colOff; x < N+_colOff; ++x)
            m_data[i++] = static_cast<T>(_mat1(y, x));
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _s0, T _s1)
{
    static_assert(M * N == 2, "Constructor should not exist!");
    m_data[0] = _s0;
    m_data[1] = _s1;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _s0, T _s1, T _s2)
{
    static_assert(M * N == 3, "Constructor should not exist!");
    m_data[0] = _s0;
    m_data[1] = _s1;
    m_data[2] = _s2;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _s0, T _s1, T _s2, T _s3)
{
    static_assert(M * N == 4, "Constructor should not exist!");
    m_data[0] = _s0;  m_data[1] = _s1;
    m_data[2] = _s2;  m_data[3] = _s3;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _s0, T _s1, T _s2, T _s3, T _s4)
{
    static_assert(M * N == 5, "Constructor should not exist!");
    m_data[0] = _s0;  m_data[1] = _s1;
    m_data[2] = _s2;  m_data[3] = _s3;
    m_data[4] = _s4;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5)
{
    static_assert(M * N == 6, "Constructor should not exist!");
    m_data[0] = _s0;  m_data[1] = _s1;
    m_data[2] = _s2;  m_data[3] = _s3;
    m_data[4] = _s4;  m_data[5] = _s5;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7)
{
    static_assert(M * N == 8, "Constructor should not exist!");
    m_data[0] = _s0;  m_data[1] = _s1;
    m_data[2] = _s2;  m_data[3] = _s3;
    m_data[4] = _s4;  m_data[5] = _s5;
    m_data[6] = _s6;  m_data[7] = _s7;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7, T _s8)
{
    static_assert(M * N == 9, "Constructor should not exist!");
    m_data[0] = _s0;  m_data[1] = _s1;  m_data[2] = _s2;
    m_data[3] = _s3;  m_data[4] = _s4;  m_data[5] = _s5;
    m_data[6] = _s6;  m_data[7] = _s7;  m_data[8] = _s8;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7, T _s8, T _s9, T _s10, T _s11)
{
    static_assert(M * N == 12, "Constructor should not exist!");
    m_data[0] = _s0;  m_data[1] = _s1;    m_data[2] = _s2;
    m_data[3] = _s3;  m_data[4] = _s4;    m_data[5] = _s5;
    m_data[6] = _s6;  m_data[7] = _s7;    m_data[8] = _s8;
    m_data[9] = _s9;  m_data[10] = _s10;  m_data[11] = _s11;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _s0, T _s1, T _s2, T _s3, T _s4, T _s5, T _s6, T _s7, T _s8, T _s9, T _s10, T _s11, T _s12, T _s13, T _s14, T _s15)
{
    static_assert(M * N == 16, "Constructor should not exist!");
    m_data[0] = _s0;    m_data[1] = _s1;    m_data[2] = _s2;    m_data[3] = _s3;
    m_data[4] = _s4;    m_data[5] = _s5;    m_data[6] = _s6;    m_data[7] = _s7;
    m_data[8] = _s8;    m_data[9] = _s9;    m_data[10] = _s10;  m_data[11] = _s11;
    m_data[12] = _s12; 	m_data[13] = _s13;  m_data[14] = _s14;  m_data[15] = _s15;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, typename T2, class>
Matrix<T, M, N>::Matrix(T1 _s0, Matrix<T2,2,1> _v12)
{
    m_data[0] = static_cast<T>(_s0);
    m_data[1] = static_cast<T>(_v12.x);
    m_data[2] = static_cast<T>(_v12.y);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, class>
Matrix<T, M, N>::Matrix(Matrix<T1,2,1> _v01, T2 _s2)
{
    m_data[0] = static_cast<T>(_v01.x);
    m_data[1] = static_cast<T>(_v01.y);
    m_data[2] = static_cast<T>(_s2);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, class>
Matrix<T, M, N>::Matrix(Matrix<T1,2,1> _v01, Matrix<T2,2,1> _v23)
{
    m_data[0] = static_cast<T>(_v01.x);
    m_data[1] = static_cast<T>(_v01.y);
    m_data[2] = static_cast<T>(_v23.x);
    m_data[3] = static_cast<T>(_v23.y);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, typename T3, class>
Matrix<T, M, N>::Matrix(T1 _s0, T2 _s1, Matrix<T3,2,1> _v23)
{
    m_data[0] = static_cast<T>(_s0);
    m_data[1] = static_cast<T>(_s1);
    m_data[2] = static_cast<T>(_v23.x);
    m_data[3] = static_cast<T>(_v23.y);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, typename T3, class>
Matrix<T, M, N>::Matrix(T1 _s0, Matrix<T2,2,1> _v12, T3 _s3)
{
    m_data[0] = static_cast<T>(_s0);
    m_data[1] = static_cast<T>(_v12.x);
    m_data[2] = static_cast<T>(_v12.y);
    m_data[3] = static_cast<T>(_s3);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, typename T3, class>
Matrix<T, M, N>::Matrix(Matrix<T1,2,1> _v01, T2 _s2, T3 _s3)
{
    m_data[0] = static_cast<T>(_v01.x);
    m_data[1] = static_cast<T>(_v01.y);
    m_data[2] = static_cast<T>(_s2);
    m_data[3] = static_cast<T>(_s3);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, class>
Matrix<T, M, N>::Matrix(T1 _s0, Matrix<T2,3,1> _v123)
{
    m_data[0] = static_cast<T>(_s0);
    m_data[1] = static_cast<T>(_v123.x);
    m_data[2] = static_cast<T>(_v123.y);
    m_data[3] = static_cast<T>(_v123.z);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, class>
Matrix<T, M, N>::Matrix(Matrix<T1,3,1> _v012, T2 _s3)
{
    m_data[0] = static_cast<T>(_v012.x);
    m_data[1] = static_cast<T>(_v012.y);
    m_data[2] = static_cast<T>(_v012.z);
    m_data[3] = static_cast<T>(_s3);
}


template<typename T, uint M, uint N>
template<typename T1, typename T2, class>
Matrix<T, M, N>::Matrix(T1 _s0, Matrix<T2,1,2> _v12)
{
    m_data[0] = static_cast<T>(_s0);
    m_data[1] = static_cast<T>(_v12.x);
    m_data[2] = static_cast<T>(_v12.y);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, class>
Matrix<T, M, N>::Matrix(Matrix<T1,1,2> _v01, T2 _s2)
{
    m_data[0] = static_cast<T>(_v01.x);
    m_data[1] = static_cast<T>(_v01.y);
    m_data[2] = static_cast<T>(_s2);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, class>
Matrix<T, M, N>::Matrix(Matrix<T1,1,2> _v01, Matrix<T2,1,2> _v23)
{
    m_data[0] = static_cast<T>(_v01.x);
    m_data[1] = static_cast<T>(_v01.y);
    m_data[2] = static_cast<T>(_v23.x);
    m_data[3] = static_cast<T>(_v23.y);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, typename T3, class>
Matrix<T, M, N>::Matrix(T1 _s0, T2 _s1, Matrix<T3,1,2> _v23)
{
    m_data[0] = static_cast<T>(_s0);
    m_data[1] = static_cast<T>(_s1);
    m_data[2] = static_cast<T>(_v23.x);
    m_data[3] = static_cast<T>(_v23.y);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, typename T3, class>
Matrix<T, M, N>::Matrix(T1 _s0, Matrix<T2,1,2> _v12, T3 _s3)
{
    m_data[0] = static_cast<T>(_s0);
    m_data[1] = static_cast<T>(_v12.x);
    m_data[2] = static_cast<T>(_v12.y);
    m_data[3] = static_cast<T>(_s3);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, typename T3, class>
Matrix<T, M, N>::Matrix(Matrix<T1,1,2> _v01, T2 _s2, T3 _s3)
{
    m_data[0] = static_cast<T>(_v01.x);
    m_data[1] = static_cast<T>(_v01.y);
    m_data[2] = static_cast<T>(_s2);
    m_data[3] = static_cast<T>(_s3);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, class>
Matrix<T, M, N>::Matrix(T1 _s0, Matrix<T2,1,3> _v123)
{
    m_data[0] = static_cast<T>(_s0);
    m_data[1] = static_cast<T>(_v123.x);
    m_data[2] = static_cast<T>(_v123.y);
    m_data[3] = static_cast<T>(_v123.z);
}

template<typename T, uint M, uint N>
template<typename T1, typename T2, class>
Matrix<T, M, N>::Matrix(Matrix<T1,1,3> _v012, T2 _s3)
{
    m_data[0] = static_cast<T>(_v012.x);
    m_data[1] = static_cast<T>(_v012.y);
    m_data[2] = static_cast<T>(_v012.z);
    m_data[3] = static_cast<T>(_s3);
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(const Quaternion& _quaternion)
{
    // Rotation composition from quaternion (remaining rest direct in matrix)
    // See http://de.wikipedia.org/wiki/Quaternion#Bezug_zu_orthogonalen_Matrizen for
    // details.
    float f2i  = 2.0f * _quaternion.i;
    float f2j  = 2.0f * _quaternion.j;
    float f2k  = 2.0f * _quaternion.k;
    float f2ri = f2i  * _quaternion.r;
    float f2rj = f2j  * _quaternion.r;
    float f2rk = f2k  * _quaternion.r;
    float f2ii = f2i  * _quaternion.i;
    float f2ij = f2j  * _quaternion.i;
    float f2ik = f2k  * _quaternion.i;
    float f2jj = f2j  * _quaternion.j;
    float f2jk = f2k  * _quaternion.j;
    float f2kk = f2k  * _quaternion.k;

    m00 = 1.0f - ( f2jj + f2kk ); m01 = f2ij - f2rk;            m02 = f2ik + f2rj;
    m10 = f2ij + f2rk;            m11 = 1.0f - ( f2ii + f2kk ); m12 = f2jk - f2ri;
    m20 = f2ik - f2rj;            m21 = f2jk + f2ri;            m22 = 1.0f - ( f2ii + f2jj );
}

// ************************************************************************* //
//                               OPERATORS                                   //
// ************************************************************************* //

// ************************************************************************* //
template<typename T, uint M, uint N>
T& Matrix<T, M, N>::operator() (uint _row, uint _col)
{
    eiAssertWeak(_row < M && _col < N, "Index out of bounds!");
    return m_data[_row * N + _col];
}

template<typename T, uint M, uint N>
T Matrix<T, M, N>::operator() (uint _row, uint _col) const
{
    eiAssertWeak(_row < M && _col < N, "Index out of bounds!");
    return m_data[_row * N + _col];
}

template<typename T, uint M, uint N>
Matrix<T,1,N>& Matrix<T, M, N>::operator() (uint _row)
{
    eiAssertWeak(_row < M, "Index out of bounds!");
    return reinterpret_cast<Matrix<T,1,N>&>(m_data[_row * N]);
}

template<typename T, uint M, uint N>
const Matrix<T,1,N>& Matrix<T, M, N>::operator() (uint _row) const
{
    eiAssertWeak(_row < M, "Index out of bounds!");
    return reinterpret_cast<const Matrix<T,1,N>&>(m_data[_row * N]);
}


// ************************************************************************* //
template<typename T, uint M, uint N>
T& Matrix<T, M, N>::operator[] (uint _index)
{
    eiAssertWeak(_index < N * M, "Index out of bounds!");
    return m_data[_index];
}

template<typename T, uint M, uint N>
T Matrix<T, M, N>::operator[] (uint _index) const
{
    eiAssertWeak(_index < N * M, "Index out of bounds!");
    return m_data[_index];
}

// ************************************************************************* //
#define CODE_GEN_MAT_MAT_OP(op)                         \
template<typename T, uint M, uint N>                    \
template<typename T1>                                   \
Matrix<RESULT_TYPE(op), M, N> Matrix<T, M, N>::operator op (const Matrix<T1,M,N>& _mat1) const\
{                                                       \
    Matrix<RESULT_TYPE(op), M, N> result;               \
    for(uint i = 0; i < N * M; ++i)                     \
        result[i] = (*this)[i] op _mat1[i];             \
    return result;                                      \
}

// ************************************************************************* //
CODE_GEN_MAT_MAT_OP(+)
CODE_GEN_MAT_MAT_OP(-)

// ************************************************************************* //
template<typename T, uint M, uint N>
Matrix<T, M, N> Matrix<T, M, N>::operator- () const
{
    Matrix<T, M, N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = -(*this)[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, uint O>
typename std::conditional<M * O == 1, RESULT_TYPE(*), Matrix<RESULT_TYPE(*), M, O>>::type
Matrix<T, M, N>::operator* (const Matrix<T1,N,O>& _mat1) const
{
    typename std::conditional<M * O == 1, RESULT_TYPE(*), Matrix<RESULT_TYPE(*), M, O>>::type result;
    for(uint m = 0; m < M; ++m)
    {
        for(uint o = 0; o < O; ++o)
        {
            RESULT_TYPE(*) acc = (*this)(m,0) * _mat1(0,o);
            for(uint n = 1; n < N; ++n)
                acc += (*this)(m,n) * _mat1(n,o);
            *(reinterpret_cast<RESULT_TYPE(*)*>(&result) + m * O + o) = acc;
        }
    }
    return result;
}

// ************************************************************************* //
template<typename T, typename T1, uint M>
Matrix<RESULT_TYPE(*), M, 1> operator* (const Matrix<T,M,1>& _mat0, const Matrix<T1,M,1>& _mat1)
{
    Matrix<RESULT_TYPE(*), M, 1> result;
    for(uint i = 0; i < M; ++i)
        result[i] = _mat0[i] * _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, typename T1, uint N>
Matrix<RESULT_TYPE(*), 1, N> operator* (const Matrix<T,1,N>& _mat0, const Matrix<T1,1,N>& _mat1)
{
    Matrix<RESULT_TYPE(*), 1, N> result;
    for(uint i = 0; i < N; ++i)
        result[i] = _mat0[i] * _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, typename T1, uint M>
Matrix<RESULT_TYPE(/), M, 1> operator/ (const Matrix<T,M,1>& _mat0, const Matrix<T1,M,1>& _mat1)
{
    Matrix<RESULT_TYPE(*), M, 1> result;
    for(uint i = 0; i < M; ++i)
        result[i] = _mat0[i] / _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, typename T1, uint N>
Matrix<RESULT_TYPE(/), 1, N> operator/ (const Matrix<T,1,N>& _mat0, const Matrix<T1,1,N>& _mat1)
{
    Matrix<RESULT_TYPE(*), 1, N> result;
    for(uint i = 0; i < N; ++i)
        result[i] = _mat0[i] / _mat1[i];
    return result;
}

// ************************************************************************* //
CODE_GEN_MAT_MAT_OP(|)
CODE_GEN_MAT_MAT_OP(&)
CODE_GEN_MAT_MAT_OP(^)

// ************************************************************************* //
template<typename T, uint M, uint N>
Matrix<T, M, N> Matrix<T, M, N>::operator~ () const
{
    Matrix<T, M, N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = ~(*this)[i];
    return result;
}

// ************************************************************************* //
// putting the implementation here causes error C2244 in visual studio.
template<typename T, uint M, uint N>
template<typename T1>
Matrix<T, M, N>& Matrix<T, M, N>::operator+= (const Matrix<T1,M,N>& _mat1)
{
    for(uint i = 0; i < N * M; ++i)
        (*this)[i] += _mat1[i];
    return *this;
}

template<typename T, uint M, uint N>
template<typename T1>
Matrix<T, M, N>& Matrix<T, M, N>::operator-= (const Matrix<T1,M,N>& _mat1)
{
    for(uint i = 0; i < N * M; ++i)
        (*this)[i] -= _mat1[i];
    return *this;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<T, M, 1>& Matrix<T, M, N>::operator*= (const Matrix<T1,M,1>& _mat1)
{
    for(uint i = 0; i < M; ++i)
        (*this)[i] *= _mat1[i];
    return *this;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<T, 1, N>& Matrix<T, M, N>::operator*= (const Matrix<T1,1,N>& _mat1)
{
    for(uint i = 0; i < N; ++i)
        (*this)[i] *= _mat1[i];
    return *this;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<T, M, N>& Matrix<T, M, N>::operator*= (const Matrix<T1,N,M>& _mat1)
{
    *this = (*this) * _mat1;
    return *this;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<T, M, 1>& Matrix<T, M, N>::operator/= (const Matrix<T1,M,1>& _mat1)
{
    for(uint i = 0; i < M; ++i)
        (*this)[i] /= _mat1[i];
    return *this;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<T, 1, N>& Matrix<T, M, N>::operator/= (const Matrix<T1,1,N>& _mat1)
{
    for(uint i = 0; i < N; ++i)
        (*this)[i] /= _mat1[i];
    return *this;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1>
Matrix<T, M, N>& Matrix<T, M, N>::operator+= (T1 _s)
{
    for(uint i = 0; i < N * M; ++i)
        (*this)[i] += _s;
    return *this;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1>
Matrix<T, M, N>& Matrix<T, M, N>::operator-= (T1 _s)
{
    for(uint i = 0; i < N * M; ++i)
        (*this)[i] -= _s;
    return *this;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1>
Matrix<T, M, N>& Matrix<T, M, N>::operator*= (T1 _s)
{
    for(uint i = 0; i < N * M; ++i)
        (*this)[i] *= _s;
    return *this;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1>
Matrix<T, M, N>& Matrix<T, M, N>::operator/= (T1 _s)
{
    for(uint i = 0; i < N * M; ++i)
        (*this)[i] /= _s;
    return *this;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline Matrix<bool,M,N> Matrix<T, M, N>::operator== (const Matrix<T,M,N>& _mat1) const
{
    Matrix<bool,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = (*this)[i] == _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline Matrix<bool,M,N> Matrix<T, M, N>::operator!= (const Matrix<T,M,N>& _mat1) const
{
    Matrix<bool,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = (*this)[i] != _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline Matrix<bool,M,N> Matrix<T, M, N>::operator<= (const Matrix<T,M,N>& _mat1) const
{
    Matrix<bool,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = (*this)[i] <= _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline Matrix<bool,M,N> Matrix<T, M, N>::operator< (const Matrix<T,M,N>& _mat1) const
{
    Matrix<bool,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = (*this)[i] < _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline Matrix<bool,M,N> Matrix<T, M, N>::operator> (const Matrix<T,M,N>& _mat1) const
{
    Matrix<bool,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = (*this)[i] > _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline Matrix<bool,M,N> Matrix<T, M, N>::operator>= (const Matrix<T,M,N>& _mat1) const
{
    Matrix<bool,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = (*this)[i] >= _mat1[i];
    return result;
}

// ********************************************************************* //
#define CODE_GEN_MAT_SCALAR_OP(op)                                                  \
template<typename T, uint M, uint N, typename T1>                                   \
inline Matrix<RESULT_TYPE(op), M, N> operator op (const Matrix<T,M,N>& _mat, T1 _s) \
{                                                                                   \
    Matrix<RESULT_TYPE(op), M, N> result;                                           \
    for(uint i = 0; i < N * M; ++i)                                                 \
        result[i] = _mat[i] op _s;                                                  \
    return result;                                                                  \
}                                                                                   \
                                                                                    \
template<typename T1, typename T, uint M, uint N>                                   \
inline Matrix<RESULT_TYPE(op), M, N> operator op (T1 _s, const Matrix<T,M,N>& _mat) \
{                                                                                   \
    Matrix<RESULT_TYPE(op), M, N> result;                                           \
    for(uint i = 0; i < N * M; ++i)                                                 \
        result[i] = _s op _mat[i];                                                  \
    return result;                                                                  \
}

// ********************************************************************* //
CODE_GEN_MAT_SCALAR_OP(+)
CODE_GEN_MAT_SCALAR_OP(-)
CODE_GEN_MAT_SCALAR_OP(*)
CODE_GEN_MAT_SCALAR_OP(/)

// ********************************************************************* //
CODE_GEN_MAT_SCALAR_OP(|)
CODE_GEN_MAT_SCALAR_OP(&)
CODE_GEN_MAT_SCALAR_OP(^)

// ********************************************************************* //
#define CODE_GEN_MAT_SCALAR_RELATION(op)                                            \
template<typename T, uint M, uint N, typename T1>                                   \
inline Matrix<bool, M, N> operator op (const Matrix<T,M,N>& _mat, T1 _s)            \
{                                                                                   \
    Matrix<bool, M, N> result;                                                      \
    for(uint i = 0; i < N * M; ++i)                                                 \
        result[i] = _mat[i] op _s;                                                  \
    return result;                                                                  \
}                                                                                   \
                                                                                    \
template<typename T1, typename T, uint M, uint N>                                   \
inline Matrix<bool, M, N> operator op (T1 _s, const Matrix<T,M,N>& _mat)            \
{                                                                                   \
    Matrix<bool, M, N> result;                                                      \
    for(uint i = 0; i < N * M; ++i)                                                 \
        result[i] = _s op _mat[i];                                                  \
    return result;                                                                  \
}

// ********************************************************************* //
CODE_GEN_MAT_SCALAR_RELATION(==)
CODE_GEN_MAT_SCALAR_RELATION(!=)
CODE_GEN_MAT_SCALAR_RELATION(<)
CODE_GEN_MAT_SCALAR_RELATION(<=)
CODE_GEN_MAT_SCALAR_RELATION(>=)
CODE_GEN_MAT_SCALAR_RELATION(>)


#undef CODE_GEN_MAT_SCALAR_OP
#undef CODE_GEN_MAT_SCALAR_RELATION


// ************************************************************************* //
//                               FUNCTIONS                                   //
// ************************************************************************* //

// ************************************************************************* //
template<typename T, uint M, uint N>
inline bool approx(const Matrix<T,M,N>& _mat0,
                   const Matrix<T,M,N>& _mat1,
                   float _epsilon)
{
    for(uint i = 0; i < N * M; ++i)
        if(abs(_mat1[i] - _mat0[i]) > _epsilon) return false;
    return true;
}

// ************************************************************************* //
inline bool approx(const Quaternion& _q0,
                   const Quaternion& _q1,
                   float _epsilon)
{
    return abs(_q0.r - _q1.r) <= _epsilon
        && abs(_q0.i - _q1.i) <= _epsilon
        && abs(_q0.j - _q1.j) <= _epsilon
        && abs(_q0.k - _q1.k) <= _epsilon;
}

// ********************************************************************* //
template<typename T, unsigned M, unsigned N, typename T1>
inline RESULT_TYPE(*) dot(const Matrix<T,M,N>& _mat0,
                          const Matrix<T1,M,N>& _mat1)
{
    RESULT_TYPE(*) sum = _mat0[0] * _mat1[0];
    for(uint i = 1; i < N * M; ++i)
        sum += _mat0[i] * _mat1[i];
    return sum;
}

// ********************************************************************* //
inline float dot(const Quaternion& _q0,
          const Quaternion& _q1)
{
    return _q0.r*_q1.r + _q0.i*_q1.i + _q0.j*_q1.j + _q0.k*_q1.k;
}

// ********************************************************************* //
template<typename T, typename T1>
inline Matrix<RESULT_TYPE(*),1,3> cross(const Matrix<T,1,3>& _v0,
                                        const Matrix<T1,1,3>& _v1)
{
    return Matrix<RESULT_TYPE(*),1,3>(_v0.y * _v1.z - _v0.z * _v1.y,
                                      _v0.z * _v1.x - _v0.x * _v1.z,
                                      _v0.x * _v1.y - _v0.y * _v1.x);
}

// ********************************************************************* //
template<typename T, typename T1>
inline Matrix<RESULT_TYPE(*),3,1> cross(const Matrix<T,3,1>& _v0,
                                        const Matrix<T1,3,1>& _v1)
{
    return Matrix<RESULT_TYPE(*),3,1>(_v0.y * _v1.z - _v0.z * _v1.y,
                                      _v0.z * _v1.x - _v0.x * _v1.z,
                                      _v0.x * _v1.y - _v0.y * _v1.x);
}

// ********************************************************************* //
template<typename T, typename T1>
inline RESULT_TYPE(*) cross(const Matrix<T,1,2>& _v0,
                            const Matrix<T1,1,2>& _v1)
{
    return _v0.x * _v1.y - _v0.y * _v1.x;
}

// ********************************************************************* //
template<typename T, typename T1>
inline RESULT_TYPE(*) cross(const Matrix<T,2,1>& _v0,
                            const Matrix<T1,2,1>& _v1)
{
    return _v0.x * _v1.y - _v0.y * _v1.x;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline T lensq(const Matrix<T,M,N>& _mat0)
{
    return dot(_mat0, _mat0);
}

// ************************************************************************* //
template<typename T>
inline auto len(const T& _elem0) -> decltype(std::sqrt(dot(_elem0, _elem0)))
{
    return sqrt(dot(_elem0, _elem0));
}

// ************************************************************************* //
template<typename T>
inline T normalize(const T& _elem0)
{
    return _elem0 / len(_elem0);
    // TODO: Test if this one is faster
    //return _mat0 * (1.0f / len(_mat0));
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<T,M,N> max(const Matrix<T,M,N>& _mat0,
                         const Matrix<T,M,N>& _mat1)
{
    Matrix<T,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = max(_mat0[i], _mat1[i]);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<T,M,N> min(const Matrix<T,M,N>& _mat0,
                         const Matrix<T,M,N>& _mat1)
{
    Matrix<T,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = min(_mat0[i], _mat1[i]);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline T max(const Matrix<T,M,N>& _mat0)
{
    T result = _mat0[0];
    for(uint i = 1; i < N * M; ++i)
        result = max(_mat0[i], result);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline T min(const Matrix<T,M,N>& _mat0)
{
    T result = _mat0[0];
    for(uint i = 1; i < N * M; ++i)
        result = min(_mat0[i], result);
    return result;
}

// ********************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<T,M,N> clamp(const Matrix<T,M,N>& _mat,
                    const Matrix<T,M,N>& _min,
                    const Matrix<T,M,N>& _max)
{
    Matrix<T,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = clamp(_mat[i], _min[i], _max[i]);
    return result;
}

// ********************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<T,M,N> clamp(const Matrix<T,M,N>& _mat,
                    T _min,
                    T _max)
{
    Matrix<T,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = clamp(_mat[i], _min, _max);
    return result;
}

// ********************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<typename details::Int<sizeof(T)>::stype,M,N> floor(const Matrix<T,M,N>& _mat)
{
    Matrix<typename details::Int<sizeof(T)>::stype,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = floor(_mat[i]);
    return result;
}

// ********************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<typename details::Int<sizeof(T)>::stype,M,N> ceil(const Matrix<T,M,N>& _mat)
{
    Matrix<typename details::Int<sizeof(T)>::stype,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = ceil(_mat[i]);
    return result;
}

// ********************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<typename details::Int<sizeof(T)>::stype,M,N> round(const Matrix<T,M,N>& _mat)
{
    Matrix<typename details::Int<sizeof(T)>::stype,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = round(_mat[i]);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,1,N> sqrt(const Matrix<T,1,N>& _v0)
{
    Matrix<T,1,N> result;
    for(uint i = 0; i < N; ++i)
        result[i] = sqrt(_v0[i]);
    return result;
}

template<typename T, unsigned M>
Matrix<T,M,1> sqrt(const Matrix<T,M,1>& _v0)
{
    Matrix<T,M,1> result;
    for(uint i = 0; i < M; ++i)
        result[i] = sqrt(_v0[i]);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,1,N> pow(const Matrix<T,1,N>& _v0, float _exponent)
{
    Matrix<T,1,N> result;
    for(uint i = 0; i < N; ++i)
        result[i] = pow(_v0[i], _exponent);
    return result;
}

template<typename T, unsigned M>
Matrix<T,M,1> pow(const Matrix<T,M,1>& _v0, float _exponent)
{
    Matrix<T,M,1> result;
    for(uint i = 0; i < M; ++i)
        result[i] = pow(_v0[i], _exponent);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline decltype(std::declval<T>() + std::declval<T>()) sum(const Matrix<T,M,N>& _mat0)
{
    decltype(std::declval<T>() + std::declval<T>()) result = _mat0[0];
    for(uint i = 1; i < N * M; ++i)
        result += _mat0[i];
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline T prod(const Matrix<T,M,N>& _mat0)
{
    T result = _mat0[0];
    for(uint i = 1; i < N * M; ++i)
        result *= _mat0[i];
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline T avg(const Matrix<T,M,N>& _mat0)
{
    return sum(_mat0) / T(M * N);
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<T,M,N> abs(const Matrix<T,M,N>& _mat0)
{
    Matrix<T,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = abs(_mat0[i]);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<T,M,N> sign(const Matrix<T,M,N>& _mat0)
{
    Matrix<T,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = sign(_mat0[i]);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<T,M,N> sgn(const Matrix<T,M,N>& _mat0)
{
    Matrix<T,M,N> result;
    for(uint i = 0; i < N * M; ++i)
        result[i] = sgn(_mat0[i]);
    return result;
}


// ************************************************************************* //
template<typename T0, typename T1, unsigned M, unsigned N>
Matrix<decltype(std::declval<T0>() * std::declval<T1>()),M,N>
    bilerp(Matrix<T0,M,N> _x00, Matrix<T0,M,N> _x01,
           Matrix<T0,M,N> _x10, Matrix<T0,M,N> _x11,
           T1 _t0, T1 _t1)
{
    // For vectors it is faster to compute the four factors first
    T1 one = static_cast<T1>(1);
    T1 t0i = (one - _t0);
    T1 t1i = (one - _t1);
    T1 t00 = t0i * t1i;
    T1 t01 = _t0 * t1i;
    T1 t10 = t0i * _t1;
    T1 t11 = _t0 * _t1;
    return _x00 * t00 + _x01 * t01
         + _x10 * t10 + _x11 * t11;
}

// ************************************************************************* //
template<unsigned M, unsigned N>
bool any(const Matrix<bool,M,N>& _mat0)
{
    for(uint i = 0; i < N * M; ++i)
        if(_mat0[i]) return true;
    return false;
}

// ************************************************************************* //
template<unsigned M, unsigned N>
bool none(const Matrix<bool,M,N>& _mat0)
{
    for(uint i = 0; i < N * M; ++i)
        if(_mat0[i]) return false;
    return true;
}

// ************************************************************************* //
template<unsigned M, unsigned N>
bool all(const Matrix<bool,M,N>& _mat0)
{
    for(uint i = 0; i < N * M; ++i)
        if(!_mat0[i]) return false;
    return true;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
Matrix<T,N,M> transpose(const Matrix<T,M,N>& _mat0)
{
    Matrix<T,N,M> result;
    // This counter avoids one index computation y*N+x in the inner loop
    uint i = 0;
    for(uint x = 0; x < N; ++x)
        for(uint y = 0; y < M; ++y)
            result[i++] = _mat0(y,x);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned N>
bool decomposeLUp(const Matrix<T,N,N>& _A, Matrix<T,N,N>& _LU, Vec<uint,N>& _p)
{
    // LUP decomposition algorithm by Cormen et al. with both matrices L and U
    // combined to one.

    // Copy and work inplace
    _LU = _A;

    // Init permutation to none
    for(uint i = 0; i < N; ++i) _p[i] = i;
    // Loop over rows in A
    for(uint r = 0; r < N-1; ++r)
    {
        // Search a non zero pivot element in current column
        T pivot = static_cast<T>(0);
        uint p = 0;
        for(uint i = r; i < N; ++i) if(abs(_LU(i,r)) > pivot) {
            pivot = abs(_LU(i,r));
            p = i;
        }
        if(pivot == static_cast<T>(0)) return false;
        // Swap the lines
        if(p != r)
        {
            // Permutation vector
            uint tmpUI = _p[r]; _p[r] = _p[p]; _p[p] = tmpUI;
            // Complete row: the first r-1 elements are already part of U
            for(uint c = 0; c < N; ++c) {
                T tmpT = _LU(r,c); _LU(r,c) = _LU(p,c); _LU(p,c) = tmpT;
            }
        }
        // Gauß elimination but keep the diagonal element
        for(uint i = r+1; i < N; ++i)
        {
            _LU(i,r) /= _LU(r,r);
            for(uint j = r+1; j < N; ++j)
                _LU(i,j) -= _LU(i,r) * _LU(r,j);
        }
    }
    return _LU(N-1,N-1) != 0.0f;
}

// ********************************************************************* //
template<typename T, unsigned M, unsigned N>
Matrix<T,M,N> solveLUp(const Matrix<T,M,M>& _LU, const Matrix<uint,M,1>& _p, const Matrix<T,M,N>& _B)
{
    Matrix<T,M,N> X;
    for(uint n = 0; n < N; ++n)
    {
        // Compute L Y = P B
        for(uint i = 0; i < M; ++i)
        {
            T sum = static_cast<T>(0);
            for(uint j = 0; j < i; ++j)
                sum += _LU(i,j) * X(j,n);
            X(i,n) = _B(_p[i],n) - sum;
        }
        // Compute U X = Y
        for(int i = M-1; i >= 0; --i)
        {
            T sum = static_cast<T>(0);
            for(uint j = i+1; j < M; ++j)
                sum += _LU(i,j) * X(j,n);
            X(i,n) = (X(i,n) - sum) / _LU(i,i);
        }
    }
    return X;
}


// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N,N> invert(const Matrix<T,N,N>& _mat0)
{
    Matrix<T,N,N> LU;
    Vec<uint32,N> p;
    if( decomposeLUp( _mat0, LU, p ) )
        return solveLUp( LU, p, identity<T,N>() );
    else return identity<T,N>();
}


// ************************************************************************* //
//                              TRANSFORMATIONS                              //
// ************************************************************************* //

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N,N> identity()
{
    Matrix<T,N,N> result(T(0));
    for(uint n = 0; n < N; ++n)
        result[n*n] = T(1);
    return result;
}

// TODO: test if this specialization is faster of if the compiler optimizes
// identity() enough
template<>
inline Matrix<float,2,2> identity<float,2>()
{
    return Mat2x2(1.0f, 0.0f,
                  0.0f, 1.0f);
}
template<>
inline Matrix<float,3,3> identity<float,3>()
{
    return Mat3x3(1.0f, 0.0f, 0.0f,
                  0.0f, 1.0f, 0.0f,
                  0.0f, 0.0f, 1.0f);
}
template<>
inline Matrix<float,4,4> identity<float,4>()
{
    return Mat4x4(1.0f, 0.0f, 0.0f, 0.0f,
                  0.0f, 1.0f, 0.0f, 0.0f,
                  0.0f, 0.0f, 1.0f, 0.0f,
                  0.0f, 0.0f, 0.0f, 1.0f);
}

// ********************************************************************* //
template<typename T, unsigned N>
Vec<T,N> sphericalCoords( const Vec<T,N>& _v0 )
{
    static_assert(N >= 2, "In 1D cartesian and spherical coordinates are the same!");
    Vec<T,N> result;
    // Accumulate the squared length over the iterations
    result[0] = sq(_v0[N-1]) + sq(_v0[N-2]);
    result[N-1] = atan2(_v0[N-1], _v0[N-2]);
    if(result[N-1] < 0.0f) result[N-1] += 2.0f * PI;
   // if( _v0[N-1] < 0.0f ) result[N-1] = 2.0f * PI - result[N-1];
    for(uint i = 2; i < N; ++i)
    {
        result[0] += sq(_v0[N-i-1]);
        result[N-i] = acos(clamp(_v0[N-i-1]/sqrt(result[0]), -1.0f, 1.0f));
    }
    result[0] = sqrt(result[0]);
    return result;
}

template<typename T, unsigned N>
RVec<T,N> sphericalCoords( const RVec<T,N>& _v0 )
{
    return *reinterpret_cast<RVec<T,N>(&sphericalCoords(*reinterpret_cast<Vec<T,N>(&_v0)));
}

// ********************************************************************* //
template<typename T, unsigned N>
Vec<T,N> cartesianCoords( const Vec<T,N>& _v0 )
{
	eiAssertWeak(_v0[0] > 0.0f, "Expected the length to be greater 0!");
    Vec<T,N> result;
    float tmp = _v0[0];
    for(uint i = 0; i < N-1; ++i)
    {
        result[i] = tmp * cos(_v0[i+1]);
        tmp *= sin(_v0[i+1]);
    }
    result[N-1] = tmp;
    return result;
}

template<typename T, unsigned N>
RVec<T,N> cartesianCoords( const RVec<T,N>& _v0 )
{
    return *reinterpret_cast<RVec<T,N>(&cartesianCoords(*reinterpret_cast<Vec<T,N>(&_v0)));
}

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N,1> transform( const Matrix<T,N,1>& _what, const Matrix<T,N+1,N+1>& _space )
{
    Matrix<T,N+1,1> t;
    // Multiply Matrix * Vector(_what,1)
    for(uint y = 0; y <= N; ++y)
    {
        // Initialize with the last component * 1
        t[y] = _space(y,N);
        // Add the other N factors
        for(uint x = 0; x < N; ++x)
            t[y] += _space(y,x) * _what[x];
    }
    // Create a reduced vector with divided components
    Matrix<T,N,1> result;
    for(uint i = 0; i < N; ++i)
        result[i] = t[i] / t[N];
    return result;
}

template<typename T, unsigned N>
Matrix<T,1,N> transform( const Matrix<T,1,N>& _what, const Matrix<T,N+1,N+1>& _space )
{
    Matrix<T,1,N+1> t;
    // Multiply Vector(_what,1) * Matrix
    for(uint x = 0; x <= N; ++x)
    {
        // Initialize with the last component * 1
        t[x] = _space(N,x);
        // Add the other N factors
        for(uint y = 0; y < N; ++y)
            t[x] += _what[y] * _space(y,x);
    }
    // Create a reduced vector with divided components
    Matrix<T,1,N> result;
    for(uint i = 0; i < N; ++i)
        result[i] = t[i] / t[N];
    return result;
}

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N,1> transformDir( const Matrix<T,N,1>& _what, const Matrix<T,N+1,N+1>& _space )
{
    Matrix<T,N,1> result;
    // Multiply Matrix * Vector(_what,0)
    for(uint y = 0; y < N; ++y)
    {
        // Initialize with the first component
        t[y] = _space(y,0) * _what[0];
        // Add the other N-1 factors
        for(uint x = 1; x < N; ++x)
            t[y] += _space(y,x) * _what[x];
    }
    return result;
}

template<typename T, unsigned N>
Matrix<T,1,N> transformDir( const Matrix<T,1,N>& _what, const Matrix<T,N+1,N+1>& _space )
{
    Matrix<T,1,N> result;
    // Multiply Vector(_what,0) * Matrix
    for(uint x = 0; x < N; ++x)
    {
        // Initialize with the first component
        t[x] = _what[0] * _space(0,x);
        // Add the other N-1 factors
        for(uint y = 1; y < N; ++y)
            t[x] += _what[y] * _space(y,x);
    }
    return result;
}

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N,1> transform( const Matrix<T,N,1>& _what, const Matrix<T,N,N>& _space )
{
    return _space * _what;
}

template<typename T, unsigned N>
Matrix<T,1,N> transform( const Matrix<T,1,N>& _what, const Matrix<T,N,N>& _space )
{
    return _what * _space;
}

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N+1,N+1> translation( const Matrix<T, N, 1>& _vector )
{
    Matrix<T,N+1,N+1> result = identity<T,N+1>();
    for(uint i = 0; i < N; ++i)
        result[i * (N+1) + N] = _vector[i];
    return result;
}

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N,N> scaling( const Matrix<T, N, 1>& _scale )
{
    Matrix<T,N,N> result(T(0));
    for(uint n = 0; n < N; ++n)
        result[n*n] = _scale[n];
    return result;
}

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N,N> scaling( T _scale )
{
    Matrix<T,N,N> result(T(0));
    for(uint n = 0; n < N; ++n)
        result[n*n] = _scale;
    return result;
}

// ************************************************************************* //
inline Matrix<float,4,4> scalingH( const Vec3& _scale )
{
    return Mat4x4(_scale.x, 0.0f,     0.0f,     0.0f,
                  0.0f,     _scale.y, 0.0f,     0.0f,
                  0.0f,     0.0f,     _scale.z, 0.0f,
                  0.0f,     0.0f,     0.0f,     1.0f);
}

// ************************************************************************* //
inline Matrix<float,4,4> scalingH( float _scale )
{
    return Mat4x4(_scale, 0.0f,   0.0f,   0.0f,
                  0.0f,   _scale, 0.0f,   0.0f,
                  0.0f,   0.0f,   _scale, 0.0f,
                  0.0f,   0.0f,   0.0f,   1.0f);
}

// ************************************************************************* //
inline Mat2x2 axis( const Vec2& _x, const Vec2& _y )
{
    return Mat2x2(_x.x, _x.y,
                  _y.x, _y.y);
}

inline Mat3x3 axis( const Vec3& _x, const Vec3& _y, const Vec3& _z )
{
    return Mat3x3(_x.x, _x.y, _x.z,
                  _y.x, _y.y, _y.z,
                  _z.x, _z.y, _z.z);
}

inline Mat4x4 axisH( const Vec3& _x, const Vec3& _y, const Vec3& _z )
{
    return Mat4x4(_x.x, _x.y, _x.z, 0.0f,
                  _y.x, _y.y, _y.z, 0.0f,
                  _z.x, _z.y, _z.z, 0.0f,
                  0.0f, 0.0f, 0.0f, 1.0f);
}

inline Mat4x4 axis( const Vec4& _x, const Vec4& _y, const Vec4& _z, const Vec4& _w )
{
    return Mat4x4(_x.x, _x.y, _x.z, _x.w,
                  _y.x, _y.y, _y.z, _y.w,
                  _z.x, _z.y, _z.z, _z.w,
                  _w.x, _w.y, _w.z, _w.w);
}

// ************************************************************************* //
inline Mat2x2 basis( const Vec2& _vector )
{
    return Mat2x2(_vector.x, _vector.y,
                  -_vector.y, _vector.x);
}

inline Mat3x3 basis( const Vec3& _vector )
{
    eiAssert(approx(len(_vector), 1.0f), "Expected normalized direction vector!");
    Vec3 y;
    if(abs(_vector.x) >= 1.0f) y = Vec3(0.0f, 1.0f, 0.0f);
    else y = normalize(Vec3(0.0f, -_vector.z, _vector.y));
    return axis(_vector, y, cross(_vector, y));
}

// ************************************************************************* //
inline Mat2x2 rotation( float _angle )
{
    float sinA = sin(_angle);
    float cosA = cos(_angle);
    return Mat2x2(cosA, -sinA,
                  sinA,  cosA);
}

// ************************************************************************* //
inline Mat3x3 rotationX( float _angle )
{
    float sinA = sin(_angle);
    float cosA = cos(_angle);
    return Mat3x3(1.0f, 0.0f,  0.0f,
                  0.0f, cosA, -sinA,
                  0.0f, sinA,  cosA);
}

inline Mat4x4 rotationXH( float _angle )
{
    return details::incrementDims(rotationX( _angle ));
}

// ************************************************************************* //
inline Mat3x3 rotationY( float _angle )
{
    float sinA = sin(_angle);
    float cosA = cos(_angle);
    return Mat3x3( cosA, 0.0f, sinA,
                   0.0f, 1.0f, 0.0f,
                  -sinA, 0.0f, cosA);
}

inline Mat4x4 rotationYH( float _angle )
{
    return details::incrementDims(rotationY( _angle ));
}

// ************************************************************************* //
inline Mat3x3 rotationZ( float _angle )
{
    float sinA = sin(_angle);
    float cosA = cos(_angle);
    return Mat3x3(cosA, -sinA, 0.0f,
                  sinA,  cosA, 0.0f,
                  0.0f,  0.0f, 1.0f);
}

inline Mat4x4 rotationZH( float _angle )
{
    return details::incrementDims(rotationZ( _angle ));
}

// ************************************************************************* //
inline Mat3x3 rotation( float _yaw, float _pitch, float _roll )
{
    float sinA = sin(_yaw);
    float cosA = cos(_yaw);
    float sinB = sin(_pitch);
    float cosB = cos(_pitch);
    float sinC = sin(_roll);
    float cosC = cos(_roll);
    return Mat3x3(cosA * cosB, cosA * sinB * sinC - sinA * cosC, cosA * sinB * cosC + sinA * sinC,
                  sinA * cosB, sinA * sinB * sinC + cosA * cosC, sinA * sinB * cosC - cosA * sinC,
                  -sinB,       cosB * sinC,                      cosB * cosC                     );
}

inline Mat4x4 rotationH( float _yaw, float _pitch, float _roll )
{
    return details::incrementDims(rotation( _yaw, _pitch, _roll ));
}

// ************************************************************************* //
inline Mat3x3 rotation( const Vec3& _v, float _angle )
{
    float sinA = sin(_angle);
    float cosA = cos(_angle);
    float iCosA = 1.0f - cosA;
    return Mat3x3(_v.x * _v.x * iCosA + cosA,        _v.x * _v.y * iCosA - _v.z * sinA, _v.x * _v.z * iCosA + _v.y * sinA,
                  _v.x * _v.y * iCosA + _v.z * sinA, _v.y * _v.y * iCosA + cosA,        _v.y * _v.z * iCosA - _v.x * sinA,
                  _v.x * _v.z * iCosA - _v.y * sinA, _v.y * _v.z * iCosA + _v.x * sinA, _v.z * _v.z * iCosA + cosA       );
}

inline Mat4x4 rotationH( const Vec3& _v, float _angle )
{
    return details::incrementDims(rotation( _v, _angle ));
}

// ************************************************************************* //
inline Mat3x3 rotation( const Quaternion& _quaternion )
{
    return Mat3x3(_quaternion);
}

inline Mat4x4 rotationH( const Quaternion& _quaternion )
{
    return details::incrementDims(Mat3x3(_quaternion));
}

// ************************************************************************* //
template<typename T, unsigned N>
inline Vec<T,N> reflect( const Vec<T,N>& _incident, const Vec<T,N>& _at )
{
    eiAssertWeak(approx(lensq(_at), 1.0f), "The reflection normal must be normalized!");
    return _incident - (static_cast<T>(2) * dot(_incident, _at)) * _at;
}

template<typename T, unsigned N>
inline RVec<T,N> reflect( const RVec<T,N>& _incident, const RVec<T,N>& _at )
{
    eiAssertWeak(approx(lensq(_at), 1.0f), "The reflection normal must be normalized!");
    return _incident - (static_cast<T>(2) * dot(_incident, _at)) * _at;
}

// ************************************************************************* //
inline Mat3x3 lookAt( const Vec3& _target, const Vec3& _up )
{
    Vec3 zAxis = normalize(_target);
    Vec3 xAxis = normalize(cross(zAxis, _up));
    Vec3 yAxis = cross(xAxis, zAxis);
    return axis( xAxis, yAxis, zAxis );
}

// ************************************************************************* //
inline Mat4x4 lookAtH( const Vec3& _target, const Vec3& _up )
{
    return details::incrementDims(lookAt( _target, _up ));
}

// ************************************************************************* //
inline Mat4x4 camera( const Vec3& _position, const Vec3& _target, const Vec3& _up )
{
    return lookAtH( _target - _position, _up ) * translation( -_position );
}

// ************************************************************************* //
inline Mat4x4 perspectiveGL( float l, float r, float b, float t, float n, float f )
{
    return Mat4x4(2.0f*n / (r-l), 0.0f,           (l+r) / (l-r),  0.0f,
                  0.0f,           2.0f*n / (t-b), (b+t) / (b-t),  0.0f,
                  0.0f,           0.0f,           (f+n) / (f-n), -2.0f*n*f / (f-n),
                  0.0f,           0.0f,           1.0f,           0.0f);
}

// ************************************************************************* //
inline Mat4x4 perspectiveGL( float _fovY, float _aspectRatio, float _n, float _f )
{
    // cot(x) == tan(π/2 - x)
    float h = tan(π * 0.5f -_fovY / 2.0f);
    float w = h / _aspectRatio;
    return Mat4x4(w,    0.0f, 0.0f,              0.0f,
                  0.0f, h,    0.0f,              0.0f,
                  0.0f, 0.0f, (_f+_n) / (_f-_n), -2.0f*_n*_f / (_f-_n),
                  0.0f, 0.0f, 1.0f,              0.0f);
}

// ************************************************************************* //
inline Mat4x4 orthographicGL( float _l, float _r, float _b, float _t, float _n, float _f )
{
    return Mat4x4(2.0f / (_r-_l), 0.0f, 0.0f, -(_r+_l) / (_r-_l),
                  0.0f, 2.0f / (_t-_b), 0.0f, -(_t+_b) / (_t-_b),
                  0.0f, 0.0f, 2.0f / (_f-_n), -(_f+_n) / (_f-_n),
                  0.0f, 0.0f, 0.0f,           1.0f);
}

// ************************************************************************* //
inline Mat4x4 perspectiveDX( float l, float r, float b, float t, float n, float f )
{
    return Mat4x4(2.0f*n / (r-l), 0.0f,           (l+r) / (l-r), 0.0f,
                  0.0f,           2.0f*n / (t-b), (b+t) / (b-t), 0.0f,
                  0.0f,           0.0f,           f / (f-n),     -n*f / (f-n),
                  0.0f,           0.0f,           1.0f,          0.0f);
}

// ************************************************************************* //
inline Mat4x4 perspectiveDX( float _fovY, float _aspectRatio, float _n, float _f )
{
    // cot(x) == tan(π/2 - x)
    float h = tan(π * 0.5f -_fovY / 2.0f);
    float w = h / _aspectRatio;
    return Mat4x4(w,    0.0f, 0.0f,         0.0f,
                  0.0f, h,    0.0f,         0.0f,
                  0.0f, 0.0f, _f / (_f-_n), -_n*_f/(_f-_n),
                  0.0f, 0.0f, 1.0f,         0.0f);
}

// ************************************************************************* //
inline Mat4x4 orthographicDX( float _l, float _r, float _b, float _t, float _n, float _f )
{
    return Mat4x4(2.0f / (_r-_l), 0.0f, 0.0f, -(_r+_l) / (_r-_l),
                  0.0f, 2.0f / (_t-_b), 0.0f, -(_t+_b) / (_t-_b),
                  0.0f, 0.0f, 1.0f / (_f-_n), -(_f+_n) / (_f-_n),
                  0.0f, 0.0f, 0.0f,           1.0f);
}