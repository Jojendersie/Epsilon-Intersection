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
    for(int i = 0; i < N * M; ++i)
        m_data[i] = _s;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1>
Matrix<T, M, N>::Matrix(const Matrix<T1,M,N>& _mat1)
{
    for(int i = 0; i < N * M; ++i)
        m_data[i] = static_cast<T>(_mat1[i]);
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
    static_assert(M * N == 12, "Constructor should not exist!");
    m_data[0] = _s0;    m_data[1] = _s1;    m_data[2] = _s2;    m_data[3] = _s3;
    m_data[4] = _s4;    m_data[5] = _s5;    m_data[6] = _s6;    m_data[7] = _s7;
    m_data[8] = _s8;    m_data[9] = _s9;    m_data[10] = _s10;  m_data[11] = _s11;
    m_data[12] = _s12; 	m_data[13] = _s13;  m_data[14] = _s14;  m_data[15] = _s15;
}


// ************************************************************************* //
//                               OPERATORS                                   //
// ************************************************************************* //

// ************************************************************************* //
template<typename T, uint M, uint N>
T& Matrix<T, M, N>::operator() (uint _row, uint _col)
{
    assertlvl2(_row < M && _col < N, "Index out of bounds!");
    return m_data[_row * N + _col];
}

template<typename T, uint M, uint N>
T Matrix<T, M, N>::operator() (uint _row, uint _col) const
{
    assertlvl2(_row < M && _col < N, "Index out of bounds!");
    return m_data[_row * N + _col];
}

// ************************************************************************* //
template<typename T, uint M, uint N>
T& Matrix<T, M, N>::operator[] (uint _index)
{
    assertlvl2(_index < N * M, "Index out of bounds!");
    return m_data[_index];
}

template<typename T, uint M, uint N>
T Matrix<T, M, N>::operator[] (uint _index) const
{
    assertlvl2(_index < N * M, "Index out of bounds!");
    return m_data[_index];
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1>
Matrix<RESULT_TYPE(+), M, N> Matrix<T, M, N>::operator+ (const Matrix<T1,M,N>& _mat1)
{
    Matrix<RESULT_TYPE(+), M, N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = (*this)[i] + _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1>
Matrix<RESULT_TYPE(-), M, N> Matrix<T, M, N>::operator- (const Matrix<T1,M,N>& _mat1)
{
    Matrix<RESULT_TYPE(-), M, N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = (*this)[i] - _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
Matrix<T, M, N> Matrix<T, M, N>::operator- ()
{
    Matrix<T, M, N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = -(*this)[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, uint O>
typename std::conditional<M * O == 1, RESULT_TYPE(*), Matrix<RESULT_TYPE(*), M, O>>::type
Matrix<T, M, N>::operator* (const Matrix<T1,N,O>& _mat1)
{
    typename std::conditional<M * O == 1, RESULT_TYPE(*), Matrix<RESULT_TYPE(*), M, O>>::type result;
    for(int m = 0; m < M; ++m)
    {
        for(int o = 0; o < O; ++o)
        {
            RESULT_TYPE(*) acc = (*this)(m,0) * _mat1(0,o);
            for(int n = 1; n < N; ++n)
                acc += (*this)(m,n) * _mat1(n,o);
            *(reinterpret_cast<RESULT_TYPE(*)*>(&result) + m * O + o) = acc;
        }
    }
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<RESULT_TYPE(*), M, 1> Matrix<T, M, N>::operator* (const Matrix<T1,M,1>& _mat1)
{
    Matrix<RESULT_TYPE(*), M, 1> result;
    for(int i = 0; i < M; ++i)
        result[i] = (*this)[i] * _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<RESULT_TYPE(*), 1, N> Matrix<T, M, N>::operator* (const Matrix<T1,1,N>& _mat1)
{
    Matrix<RESULT_TYPE(*), 1, N> result;
    for(int i = 0; i < N; ++i)
        result[i] = (*this)[i] * _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<RESULT_TYPE(/), M, 1> Matrix<T, M, N>::operator/ (const Matrix<T1,M,1>& _mat1)
{
    Matrix<RESULT_TYPE(*), M, 1> result;
    for(int i = 0; i < M; ++i)
        result[i] = (*this)[i] / _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<RESULT_TYPE(/), 1, N> Matrix<T, M, N>::operator/ (const Matrix<T1,1,N>& _mat1)
{
    Matrix<RESULT_TYPE(*), 1, N> result;
    for(int i = 0; i < N; ++i)
        result[i] = (*this)[i] / _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline bool Matrix<T, M, N>::operator== (const Matrix<T,M,N>& _mat1)
{
    for(int i = 0; i < N * M; ++i)
        if((*this)[i] != _mat1[i]) return false;
    return true;
}

// ********************************************************************* //
template<typename T, uint M, uint N, typename T1>
inline Matrix<RESULT_TYPE(+), M, N> operator+ (const Matrix<T,M,N>& _mat, T1 _s)
{
    Matrix<RESULT_TYPE(+), M, N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = _mat[i] + _s;
    return result;
}

template<typename T1, typename T, uint M, uint N>
inline Matrix<RESULT_TYPE(+), M, N> operator+ (T1 _s, const Matrix<T,M,N>& _mat)
{
    Matrix<RESULT_TYPE(+), M, N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = _s + _mat[i];
    return result;
}

// ********************************************************************* //
template<typename T, uint M, uint N, typename T1>
inline Matrix<RESULT_TYPE(-), M, N> operator- (const Matrix<T, M, N>& _mat, T1 _s)
{
    Matrix<RESULT_TYPE(-), M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _mat[i] - _s;
    return result;
}

template<typename T1, typename T, uint M, uint N>
inline Matrix<RESULT_TYPE(-), M, N> operator- (T1 _s, const Matrix<T, M, N>& _mat)
{
    Matrix<RESULT_TYPE(-), M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _s - _mat[i];
    return result;
}

// ********************************************************************* //
template<typename T, uint M, uint N, typename T1>
inline Matrix<RESULT_TYPE(/), M, N> operator/ (const Matrix<T, M, N>& _mat, T1 _s)
{
    Matrix<RESULT_TYPE(/), M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _mat[i] / _s;
    return result;
}

template<typename T1, typename T, uint M, uint N>
inline Matrix<RESULT_TYPE(/), M, N> operator/ (T1 _s, const Matrix<T, M, N>& _mat)
{
    Matrix<RESULT_TYPE(/), M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _s / _mat[i];
    return result;
}

// ********************************************************************* //
template<typename T, uint M, uint N, typename T1>
inline Matrix<RESULT_TYPE(*), M, N> operator* (const Matrix<T,M,N>& _mat, T1 _s)
{
    Matrix<RESULT_TYPE(*), M, N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = _mat[i] * _s;
    return result;
}

template<typename T1, typename T, uint M, uint N>
inline Matrix<RESULT_TYPE(*), M, N> operator* (T1 _s, const Matrix<T,M,N>& _mat)
{
    Matrix<RESULT_TYPE(*), M, N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = _s * _mat[i];
    return result;
}

// ************************************************************************* //
//                               FUNCTIONS                                   //
// ************************************************************************* //

// ************************************************************************* //
template<typename T, uint M, uint N>
inline bool approx(const Matrix<T,M,N>& _mat0,
                   const Matrix<T,M,N>& _mat1,
                   float _epsilon)
{
    for(int i = 0; i < N * M; ++i)
        if(abs(_mat1[i] - _mat0[i]) > _epsilon) return false;
    return true;
}

// ********************************************************************* //
template<typename T, unsigned M, unsigned N, typename T1>
inline RESULT_TYPE(*) dot(const Matrix<T,M,N>& _mat0,
                          const Matrix<T1,M,N>& _mat1)
{
    RESULT_TYPE(*) sum = _mat0[0] * _mat1[0];
    for(int i = 1; i < N * M; ++i)
        sum += _mat0[i] * _mat1[i];
    return sum;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline T lensq(const Matrix<T,M,N>& _mat0)
{
    return dot(_mat0, _mat0);
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline decltype(sqrt(std::declval<T>())) len(const Matrix<T,M,N>& _mat0)
{
    return sqrt(dot(_mat0, _mat0));
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<T,M,N> max(const Matrix<T,M,N>& _mat0,
                         const Matrix<T,M,N>& _mat1)
{
    Matrix<T,M,N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = max(_mat0[i], _mat1[i]);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<T,M,N> min(const Matrix<T,M,N>& _mat0,
                         const Matrix<T,M,N>& _mat1)
{
    Matrix<T,M,N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = min(_mat0[i], _mat1[i]);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline T max(const Matrix<T,M,N>& _mat0)
{
    T result = _mat0[0];
    for(int i = 1; i < N * M; ++i)
        result = max(_mat0[i], result);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline T min(const Matrix<T,M,N>& _mat0)
{
    T result = _mat0[0];
    for(int i = 1; i < N * M; ++i)
        result = min(_mat0[i], result);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline T sum(const Matrix<T,M,N>& _mat0)
{
    T result = _mat0[0];
    for(int i = 1; i < N * M; ++i)
        result += _mat0[i];
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
    for(int i = 0; i < N * M; ++i)
        result[i] = abs(_mat0[i]);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<T,M,N> sign(const Matrix<T,M,N>& _mat0)
{
    Matrix<T,M,N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = sign(_mat0[i]);
    return result;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
inline Matrix<T,M,N> sgn(const Matrix<T,M,N>& _mat0)
{
    Matrix<T,M,N> result;
    for(int i = 0; i < N * M; ++i)
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