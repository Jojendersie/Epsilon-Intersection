/// \brief Header to extend vectors for std:: functions like to_string and hash
#pragma once

#include "vector.hpp"

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


    /// \brief Container need to compare items for equality too
    template <typename T, uint M, uint N>
    struct equal_to<ei::Matrix<T, M, N>>
    {
        bool operator()(const ei::Matrix<T, M, N>& _lhs, const ei::Matrix<T, M, N>& _rhs) const
        {
            return all(_lhs == _rhs);
        }
    };
} // namespace std