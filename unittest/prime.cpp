#include "ei/prime.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace ei;
using namespace std;

bool test_primes()
{
    bool result = true;

    TEST(!isPrime(1), "1 is not a prime!");
    TEST(!isPrime(-2), "-2 is not a prime!");
    TEST(isPrime(2), "2 is a prime!");
    TEST(isPrime(3), "3 is a prime!");
    TEST(!isPrime(4), "4 is not a prime!");
    TEST(isPrime(5), "5 is a prime!");
    TEST(isPrime(1009), "1009 is a prime!");
    TEST(!isPrime((1 << 30) - 1), "2^30-1 is not a prime!");

    TEST(nextPrime(2) == 3, "nextPrime(2) is 3!");
    TEST(nextPrime(3) == 5, "nextPrime(3) is 5!");
    TEST(nextPrime(1000) == 1009, "nextPrime(1000) is 1009!");
    TEST(nextPrime((1 << 30) - 1) == 1073741827, "nextPrime(1 << 30 - 1) is 1073741827!");

    TEST(nextPrimeGreaterOrEqual(3) == 3, "nextPrimeGreaterOrEqual(3) is 3!");
    TEST(nextPrimeGreaterOrEqual(4) == 5, "nextPrimeGreaterOrEqual(4) is 5!");
    TEST(nextPrimeGreaterOrEqual(9) == 11, "nextPrimeGreaterOrEqual(9) is 11!");
    TEST(nextPrimeGreaterOrEqual(1073741827) == 1073741827, "nextPrimeGreaterOrEqual(1073741827) is 1073741827!");

    // Test currently disabled because isPrime_MRTest uses divisor test for small numbers
    // as performance fallback.
    /*for(uint i = 4; i < 1000; ++i)
    {
        if(details::isPrime_DivisorBasedTest(i) != details::isPrime_MRTest(i))
        {
            TEST(false, "Miller-Rabin test provides a wrong answer for " << i);
            break;
        }
    }*/

    // Benchmark for small primes
    uint64 start = ticks();
    for(uint i = 0; i < 100000; ++i)
        eatMyDummy(isPrime(i));
    uint64 end = ticks();
    std::cerr << "[timing] isPrime test for small numbers: " << deltaTicksToMicroSeconds(end-start) << u8" µs.\n";

    // Benchmark for larger primes
    start = ticks();
    for(uint i = 1000000; i < 10000000; ++i)
        eatMyDummy(isPrime(i));
    end = ticks();
    std::cerr << "[timing] isPrime test for largerer numbers: " << deltaTicksToMilliSeconds(end-start) << " ms.\n";

    return result;
}