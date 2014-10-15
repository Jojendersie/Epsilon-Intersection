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
    static_assert(M * N == 16, "Constructor should not exist!");
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
Matrix<RESULT_TYPE(+), M, N> Matrix<T, M, N>::operator+ (const Matrix<T1,M,N>& _mat1) const
{
    Matrix<RESULT_TYPE(+), M, N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = (*this)[i] + _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1>
Matrix<RESULT_TYPE(-), M, N> Matrix<T, M, N>::operator- (const Matrix<T1,M,N>& _mat1) const
{
    Matrix<RESULT_TYPE(-), M, N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = (*this)[i] - _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
Matrix<T, M, N> Matrix<T, M, N>::operator- () const
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
Matrix<T, M, N>::operator* (const Matrix<T1,N,O>& _mat1) const
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
Matrix<RESULT_TYPE(*), M, 1> Matrix<T, M, N>::operator* (const Matrix<T1,M,1>& _mat1) const
{
    Matrix<RESULT_TYPE(*), M, 1> result;
    for(int i = 0; i < M; ++i)
        result[i] = (*this)[i] * _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<RESULT_TYPE(*), 1, N> Matrix<T, M, N>::operator* (const Matrix<T1,1,N>& _mat1) const
{
    Matrix<RESULT_TYPE(*), 1, N> result;
    for(int i = 0; i < N; ++i)
        result[i] = (*this)[i] * _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<RESULT_TYPE(/), M, 1> Matrix<T, M, N>::operator/ (const Matrix<T1,M,1>& _mat1) const
{
    Matrix<RESULT_TYPE(*), M, 1> result;
    for(int i = 0; i < M; ++i)
        result[i] = (*this)[i] / _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<typename T1, class>
Matrix<RESULT_TYPE(/), 1, N> Matrix<T, M, N>::operator/ (const Matrix<T1,1,N>& _mat1) const
{
    Matrix<RESULT_TYPE(*), 1, N> result;
    for(int i = 0; i < N; ++i)
        result[i] = (*this)[i] / _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline Matrix<bool,M,N> Matrix<T, M, N>::operator== (const Matrix<T,M,N>& _mat1) const
{
    Matrix<bool,M,N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = (*this)[i] == _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline Matrix<bool,M,N> Matrix<T, M, N>::operator<= (const Matrix<T,M,N>& _mat1) const
{
    Matrix<bool,M,N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = (*this)[i] <= _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline Matrix<bool,M,N> Matrix<T, M, N>::operator< (const Matrix<T,M,N>& _mat1) const
{
    Matrix<bool,M,N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = (*this)[i] < _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline Matrix<bool,M,N> Matrix<T, M, N>::operator> (const Matrix<T,M,N>& _mat1) const
{
    Matrix<bool,M,N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = (*this)[i] > _mat1[i];
    return result;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
inline Matrix<bool,M,N> Matrix<T, M, N>::operator>= (const Matrix<T,M,N>& _mat1) const
{
    Matrix<bool,M,N> result;
    for(int i = 0; i < N * M; ++i)
        result[i] = (*this)[i] >= _mat1[i];
    return result;
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
inline Matrix<bool, M, N> operator< (const Matrix<T, M, N>& _mat, T1 _s)
{
    Matrix<bool, M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _mat[i] < _s;
    return result;
}

template<typename T1, typename T, uint M, uint N>
inline Matrix<bool, M, N> operator< (T1 _s, const Matrix<T, M, N>& _mat)
{
    Matrix<bool, M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _s < _mat[i];
    return result;
}

// ********************************************************************* //
template<typename T, uint M, uint N, typename T1>
inline Matrix<bool, M, N> operator<= (const Matrix<T, M, N>& _mat, T1 _s)
{
    Matrix<bool, M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _mat[i] <= _s;
    return result;
}

template<typename T1, typename T, uint M, uint N>
inline Matrix<bool, M, N> operator<= (T1 _s, const Matrix<T, M, N>& _mat)
{
    Matrix<bool, M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _s <= _mat[i];
    return result;
}

// ********************************************************************* //
template<typename T, uint M, uint N, typename T1>
inline Matrix<bool, M, N> operator>= (const Matrix<T, M, N>& _mat, T1 _s)
{
    Matrix<bool, M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _mat[i] >= _s;
    return result;
}

template<typename T1, typename T, uint M, uint N>
inline Matrix<bool, M, N> operator>= (T1 _s, const Matrix<T, M, N>& _mat)
{
    Matrix<bool, M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _s >= _mat[i];
    return result;
}

// ********************************************************************* //
template<typename T, uint M, uint N, typename T1>
inline Matrix<bool, M, N> operator> (const Matrix<T, M, N>& _mat, T1 _s)
{
    Matrix<bool, M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _mat[i] > _s;
    return result;
}

template<typename T1, typename T, uint M, uint N>
inline Matrix<bool, M, N> operator> (T1 _s, const Matrix<T, M, N>& _mat)
{
    Matrix<bool, M, N> result;
    for (int i = 0; i < N * M; ++i)
        result[i] = _s > _mat[i];
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
inline decltype(std::declval<T>() + std::declval<T>()) sum(const Matrix<T,M,N>& _mat0)
{
    decltype(std::declval<T>() + std::declval<T>()) result = _mat0[0];
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

// ************************************************************************* //
template<unsigned M, unsigned N>
bool any(const Matrix<bool,M,N>& _mat0)
{
    for(int i = 0; i < N * M; ++i)
        if(_mat0[i]) return true;
    return false;
}

// ************************************************************************* //
template<unsigned M, unsigned N>
bool none(const Matrix<bool,M,N>& _mat0)
{
    for(int i = 0; i < N * M; ++i)
        if(_mat0[i]) return false;
    return true;
}

// ************************************************************************* //
template<unsigned M, unsigned N>
bool all(const Matrix<bool,M,N>& _mat0)
{
    for(int i = 0; i < N * M; ++i)
        if(!_mat0[i]) return false;
    return true;
}

// ************************************************************************* //
template<typename T, unsigned M, unsigned N>
Matrix<T,N,M> transpose(const Matrix<T,M,N>& _mat0)
{
    Matrix<T,N,M> result;
    // This counter avoids one index computation y*N+x in the inner loop
    int i = 0;
    for(int x = 0; x < N; ++x)
        for(int y = 0; y < M; ++y)
            result[i++] = _mat0(y,x);
    return result;
}


// ************************************************************************* //
//                              TRANSFORMATIONS                              //
// ************************************************************************* //

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N,N> identity()
{
    Matrix<T,N,N> result(T(0));
    for(int n = 0; n < N; ++n)
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

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N+1,N+1> translation( const Matrix<T, N, 1>& _vector )
{
    Matrix<T,N+1,N+1> result = identity<T,N+1>();
    for(int i = 0; i < N; ++i)
        result[i * (N+1) + N] = _vector[i];
}

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N,N> scaling( const Matrix<T, N, 1>& _scale )
{
    Matrix<T,N,N> result(T(0));
    for(int n = 0; n < N; ++n)
        result[n*n] = _scale[n];
    return result;
}

// ************************************************************************* //
template<typename T, unsigned N>
Matrix<T,N,N> scaling( T _scale )
{
    Matrix<T,N,N> result(T(0));
    for(int n = 0; n < N; ++n)
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