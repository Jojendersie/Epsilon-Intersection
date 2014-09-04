// ************************************************************************* //
//								 CONSTRUCTORS   						     //
// ************************************************************************* //

// ************************************************************************* //
template<typename T, uint M, uint N>
Matrix<T, M, N>::Matrix()
{
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _v0, T _v1)
{
	static_assert(M * N == 2, "Constructor should not exist!");
	m_data[0] = _v0;
	m_data[1] = _v1;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _v0, T _v1, T _v2)
{
	static_assert(M * N == 3, "Constructor should not exist!");
	m_data[0] = _v0;
	m_data[1] = _v1;
	m_data[2] = _v2;
}

// ************************************************************************* //
template<typename T, uint M, uint N>
template<class>
Matrix<T, M, N>::Matrix(T _v0, T _v1, T _v2, T _v3)
{
	static_assert(M * N == 4, "Constructor should not exist!");
	m_data[0] = _v0;
	m_data[1] = _v1;
	m_data[2] = _v2;
	m_data[3] = _v3;
}


// ************************************************************************* //
//								 OPERATORS      						     //
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
//								 FUNCTIONS								     //
// ************************************************************************* //

// ************************************************************************* //
template<typename T, uint M, uint N>
inline bool approx(const Matrix<T,M,N>& _mat0, const Matrix<T,M,N>& _mat1, float _epsilon)
{
	for(int i = 0; i < N * M; ++i)
		if(abs(_mat1[i] - _mat0[i]) > _epsilon) return false;
	return true;
}
