/// \brief Header to extend vectors for std:: functions like to_string and hash
#pragma once

#include "vector.hpp"
#include <ostream>
#include <iomanip>

namespace std {

    template <class T> struct hash;
    template <class T> struct equal_to;

    /// \brief Custom hash function for vectors.
    template <typename T, uint M, uint N>
    struct hash<ei::Matrix<T, M, N>>
    {
        using argument_type = ei::Matrix<T, M, N>;
        using result_type = std::size_t;

        std::size_t operator()(const ei::Matrix<T, M, N>& _key) const
        {
            hash<T> hashfunc;
            size_t h = 0;
            for( int i = 0; i < M*N; ++i )
            {
                // A sum of hashes is universal:
                // http://en.wikipedia.org/wiki/Universal_hashing#Hashing_vectors
                h += hashfunc(_key[i]);
            }
            return h;
        }
    };

    /// \brief Stronger hash function for 2D vectors.
    /// \details This function generates different hashes for permutations too.
    template <typename T>
    struct hash<ei::Vec<T, 2>>
    {
        using argument_type = ei::Vec<T, 2>;
        using result_type = std::size_t;

        std::size_t operator()(const ei::Vec<T, 2>& _key) const
        {
            hash<T> hashfunc;
            size_t h = 0xbd73a0fb;
            h += hashfunc(_key[0]) * 0xf445f0a9;
            h += hashfunc(_key[1]) * 0x5c23b2e1;
            return h;
        }
    };

    /// \brief Stronger hash function for 3D vectors.
    /// \details This function generates different hashes for permutations too.
    template <typename T>
    struct hash<ei::Vec<T, 3>>
    {
        using argument_type = ei::Vec<T, 3>;
        using result_type = std::size_t;

        std::size_t operator()(const ei::Vec<T, 3>& _key) const
        {
            hash<T> hashfunc;
            size_t h = 0xbd73a0fb;
            h += hashfunc(_key[0]) * 0xf445f0a9;
            h += hashfunc(_key[1]) * 0x5c23b2e1;
            h += hashfunc(_key[2]) * 0x7d25f695;
            return h;
        }
    };

    /// \brief Stronger hash function for 4D vectors.
    /// \details This function generates different hashes for permutations too.
    template <typename T>
    struct hash<ei::Vec<T, 4>>
    {
        using argument_type = ei::Vec<T, 4>;
        using result_type = std::size_t;

        std::size_t operator()(const ei::Vec<T, 4>& _key) const
        {
            hash<T> hashfunc;
            size_t h = 0xbd73a0fb;
            h += hashfunc(_key[0]) * 0xf445f0a9;
            h += hashfunc(_key[1]) * 0x5c23b2e1;
            h += hashfunc(_key[2]) * 0x7d25f695;
            h += hashfunc(_key[2]) * 0x11d6a9f3;
            return h;
        }
    };




    /// \brief Pretty printer for vectors and matrices.
    /// \details The output syntax differs for matrices and vectors by using different
    ///     brackets.
    ///     Column vector: (X, Y, Z) or row vector: (X, Y, Z)'
    ///
    ///     |m00 m01|
    ///     |m10 m11|
    ///     The printer does not add line breaks at the end. Only the matrix version
    ///     adds line breaks between rows.
    template <typename T, uint N> // Row vector
    std::ostream& operator << (std::ostream& _os, const ei::Matrix<T, 1, N>& _mat)  // TESTED
    {
        _os << '(';
        for(uint i = 0; i < N-1; ++i)
            _os << _mat[i] << ", ";
        _os << _mat[N-1] << ")'";

        return _os;
    }

    template <typename T, uint M> // Column vector
    std::ostream& operator << (std::ostream& _os, const ei::Matrix<T, M, 1>& _mat)  // TESTED
    {
        _os << '(';
        for(uint i = 0; i < M-1; ++i)
            _os << _mat[i] << ", ";
        _os << _mat[M-1] << ")";

        return _os;
    }

    template <typename T, uint M, uint N> // Matrix
    std::ostream& operator << (std::ostream& _os, const ei::Matrix<T, M, N>& _mat)  // TESTED
    {
        // Set to scientific to align columns
        //auto flags = _os.flags();
        //  _os << std::scientific;
        for(uint j = 0; j < M; ++j)
        {
            _os << '|';
            for(uint i = 0; i < N-1; ++i)
                _os << std::setw(11) << _mat[i + N * j] << ' ';
            _os << std::setw(11) << _mat[N-1 + N * j] << '|';
            if(j < M-1)
                _os << '\n';
        }
        // Revert stream format
        //_os.setf(flags);

        return _os;
    }
} // namespace std