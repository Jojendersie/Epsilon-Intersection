#pragma once

// You must provide an eiconfig.hpp at a location where your build tool finds it.
// To create that file, copy the defaultconfig.hpp to your chosen location and
// configure the Epsilon-Intersection library by switching the defines on and off.
#include <eiconfig.hpp>

// Just show the namespace once to make sure it exists when aliasing.
namespace ei {}
#ifdef EI_USE_UNICODE_NAMES
    namespace ε = ei;
#endif
// An easter egg for german speaking people ;-).
namespace egg = ei;


// Define the assertions depending on the option level.
#if defined(DEBUG) || defined(_DEBUG)
#   ifdef EI_DEBUG_BREAK
#   undef EI_DEBUG_BREAK
#   endif

#   if defined(__CUDA_ARCH__)
#       include <cstdio>
#       define EI_DEBUG_BREAK(msg) printf("Assertion violated: %s", msg)
#   elif defined(_WIN32)
#       define EI_DEBUG_BREAK(msg) __debugbreak()
#   else
#       define EI_DEBUG_BREAK(msg) ::kill(0, SIGTRAP)
#   endif

    /// \brief Default assert macro for all our needs. Use this instead of <cassert>
    /// \details This form of assertion breaks at the very point of code
    ///    where you placed it. Further you are able to step over the assertion
    ///    in any condition.
    ///
    ///    The do while loop is the only way to guarantee robustness of the macro.
    ///    It forces you to set a semicolon and can be used within if, for, ...
    ///    without curled brackets:
    ///        if(...)
    ///            eiAssert(..., "Even if compiled out the if statement does not
    ///                             infer with the next line.");
    ///        somestatement;
#   define eiAssert(condition, errorMessage)   \
        do {                                   \
            if((condition) == false) {         \
                EI_DEBUG_BREAK(errorMessage);  \
            }                                  \
        } while((void)0,0)

    /// \brief Optional assert macro for highly frequented code.
    /// \details \see eiAssert()
#   if EI_ASSERTION_LEVEL > 1
#       define eiAssertWeak(condition, errorMessage) \
            do {                                   \
                if((condition) == false) {         \
                    EI_DEBUG_BREAK(errorMessage);  \
                }                                  \
            } while((void)0,0)
#   endif

#else

    /// \brief Empty replacement if switched off.
    /// \details \see assertlvl1()
#   define eiAssert(condition, errorMessage) do { } while((void)0,0)
#   if EI_ASSERTION_LEVEL > 1
#       define eiAssertWeak(condition, errorMessage) do { } while((void)0,0)
#   endif
#endif


// Dependend
namespace ei {
#if defined(DOUBLE_PRECISION)
	/// \brief Default threshold value for numerical instable tests
	constexpr double EPSILON = 1e-14;
#else
	/// \brief Default threshold value for numerical instable tests
	constexpr float EPSILON = 1e-6f;
#endif
}


// Make (partial) CUDA compatible
#ifdef __CUDACC__
#define EIAPI	__host__ __device__ inline
#else
#define EIAPI inline
#endif
