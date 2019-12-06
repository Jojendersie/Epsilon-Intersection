namespace ei { namespace details {

    // Test numbers >= 4 correctly testing all possible dividers up to the
    // square root. This test avoids some ifs for the trivial small primes.
    //
    // Based on: http://stackoverflow.com/questions/30052316/find-next-prime-number-algorithm
    // Main contribution: all primes are in the form of 6k+-1 -> only test those numbers.
    // https://www.quora.com/Is-every-prime-number-other-than-2-and-3-of-the-form-6k%C2%B11
    template<typename T>
    constexpr EIAPI bool isPrime_DivisorBasedTest(T _n)
    {
        eiAssertWeak(_n >= 4, "Wrong input for the reduced prime test");
        if(_n % 2 == 0 || _n % 3 == 0) return false;
        T divisor = 6;
        while(divisor * divisor - 2 * divisor + 1 <= _n) // (divisor-1)^2 <= n
        {
            if(_n % (divisor - 1) == 0)
                return false;

            if(_n % (divisor + 1) == 0)
                return false;

            divisor += 6;
        }

        return true;
    }

    // Test if _n is a strong pseudo-prime respective to base _a.
    // The test fails for _a^2 or (_n-1)^2 >= 2^64. I.e. it is only correct
    // for 32 bit integers.
    // Source: https://de.wikipedia.org/wiki/Miller-Rabin-Test
    // _n must be odd and 1 < _a < _n-1
    constexpr EIAPI bool millerRabinTest(const uint32 _n, const uint32 _a)
    {
        const uint32 n1 = _n - 1;
        uint32 d = n1 >> 1;
        int j = 1;
        while((d & 1) == 0) d >>= 1, ++j;
        /*if((d & 0xffff) == 0) d >>= 16, j += 16;
        if((d & 0xff) == 0) d >>= 8, j += 8;
        if((d & 0xf) == 0) d >>= 4, j += 4;
        if((d & 0x3) == 0) d >>= 2, j += 2;
        if((d & 0x1) == 0) d >>= 1, ++j;*/
        uint64 t = _a, p = _a; // needs to have 64 bit to avoid overflow of squares.
        while(d >>= 1) { // square and multiply: a^d mod n
            p = (p * p) % _n;
            if(d & 1) t = (t * p) % _n;
        }
        if(t == 1 || t == n1) return true; // _n is probably prime
        for(int k = 1; k < j; ++k) {
            t = (t * t) % _n;
            if(t == n1) return true; // _n is probably prime
            if(t <= 1) break;
        }
        return false; // _n is surely not prime
    }

    // Deterministic Miller-Rabin test for 32 bit integers
    constexpr EIAPI bool isPrime_MRTest(uint32 _n)
    {
        if(_n % 2 == 0 || _n % 3 == 0) return false;
        // _n must be larger than 61 for correctness. The current value is chosen
        // for performance reasons.
        if(_n < 4000000) return isPrime_DivisorBasedTest(_n);
        if(!millerRabinTest(_n, 2)) return false;
        if(!millerRabinTest(_n, 7)) return false;
        return millerRabinTest(_n, 61);
    }

}} // namespace details
