#pragma once

#include "elementarytypes.hpp"
#include "details/primetest.hpp"

namespace ei {

    /// Test if an (positive) integer is prime.
    template<typename T>
    constexpr EIAPI bool isPrime(T _n)
    {
        // No negative primes and no 0 and 1.
        if(_n < 2) return false;
        if(_n == 2 || _n == 3) return true;
        return details::isPrime_MRTest(_n);
    }
    template<>
    constexpr EIAPI bool isPrime(uint64 _n)
    {
        // No negative primes and no 0 and 1.
        if(_n < 2) return false;
        if(_n == 2 || _n == 3) return true;
        // BAD: fall back to slower test for larger numbers, because MRtest is
        // only implemented for 32bit numbers (otherwise big-number libraries
        // would be necessary.
        return _n <= 0xffffffff ?
            details::isPrime_MRTest(static_cast<uint32>(_n)) :
            details::isPrime_DivisorBasedTest(_n);
    }
    template<>
    constexpr EIAPI bool isPrime(int64 _n) { return isPrime(static_cast<uint64>(_n)); }


    /// Get a prime number which is greater or equal than the given input.
    /// \returns 2 for all _numbers <= 2, the _number itself if it is prime
    ///    and the smallest prime greater than _number if it is not prime.
    template<typename T>
    constexpr EIAPI T nextPrimeGreaterOrEqual(T _number)
    {
        if(_number <= 2) return 2;
        if(_number == 3) return _number;
        // Find the first number 6k > number, but if there is a 6k(+1)==number
        // test 6k+1 separately.
        T r = _number % 6;
        if(r <= 1)
        {
            _number += 1-r;
            if(details::isPrime_DivisorBasedTest(_number)) return _number;
            _number += 5;
        } else _number += 6-r;
        while(true)
        {
            if(details::isPrime_DivisorBasedTest(_number - 1)) return _number - 1;
            if(details::isPrime_DivisorBasedTest(_number + 1)) return _number + 1;
            _number += 6;
        }
    }

    /// Get a prime number which is greater than the given input.
    /// \returns 2 for all _numbers < 2 and the smallest prime greater than
    ///    _number otherwise.
    template<typename T>
    constexpr EIAPI T nextPrime(T _number)
    {
        return nextPrimeGreaterOrEqual(_number + 1);
    }

} // namespace ei
