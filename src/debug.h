#ifndef DEBUG_H_
#define DEBUG_H_

#ifdef DEBUG
#if defined(__APPLE__)
    #if defined(__x86_64__) || defined(__i386__)
        #include <xmmintrin.h>
        void EnableFloatingPointExceptions() {
            _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &
                                   ~(_MM_MASK_INVALID |
                                     _MM_MASK_DIV_ZERO |
                                     _MM_MASK_OVERFLOW));
        }
    #elif defined(__aarch64__) || defined(__arm64__)
        #include <fenv.h>
        #pragma STDC FENV_ACCESS ON
        void EnableFloatingPointExceptions() {
            fenv_t env;
            fegetenv(&env);
            env.__fpcr = env.__fpcr & ~((1 << 8) | (1 << 9));
            fesetenv(&env);
        }
    #else
        void EnableFloatingPointExceptions() {}
    #endif
#elif defined(__linux__)
    #include <fenv.h>
    void EnableFloatingPointExceptions() {
        feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    }
#else
    void EnableFloatingPointExceptions() {}
#endif
#endif // DEBUG

#endif // DEBUG_H_