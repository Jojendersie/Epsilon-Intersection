namespace details {
    // ********************************************************************* //
    //                        IMPLEMETATION UTILITIES                        //
    // ********************************************************************* //

    // ********************************************************************* //
    /// \brief Lift vector or squared matrix to homogeneous space.
    /// \details Adds a row and a column with zeros to a matrix and sets the
    ///    new diagonal element to 1.
    ///    A     =>    A 0
    ///                0 1
    ///
    ///   Appends 1 to vectors: v    =>   (v 1)
    template<typename T, unsigned N>
    ei::Matrix<T,N+1,N+1> incrementDims( const ei::Matrix<T,N,N>& _mat0 )
    {
        ei::Matrix<T,N+1,N+1> result;
        // Indices for _mat0 and result
        unsigned i = 0, j = 0;
        for(unsigned y = 0; y < N; ++y)
        {
            // Copy NxN part
            for(unsigned x = 0; x < N; ++x)
                result[j++] = _mat0[i++];
            // New element at the end of the row is 0
            result[j++] = T(0);
        }
        // Fill new row
        for(unsigned x = 0; x < N; ++x)
            result[j + x] = T(0);
        result[j + N] = T(1);
        return result;
    }

    template<typename T, unsigned N>
    ei::Matrix<T,N+1,1> incrementDims( const ei::Matrix<T,N,1>& _v0 )
    {
        ei::Matrix<T,N+1,1> result;
        for(unsigned i = 0; i < N; ++i)
            result[i] = _v0[i];
        result[N] = T(1);
        return result;
    }

    template<typename T, unsigned N>
    ei::Matrix<T,1,N+1> incrementDims( const ei::Matrix<T,1,N>& _v0 )
    {
        ei::Matrix<T,1,N+1> result;
        for(unsigned i = 0; i < N; ++i)
            result[i] = _v0[i];
        result[N] = T(1);
        return result;
    }
}