#pragma once

#if defined(__SSE2__)
    #include <xmmintrin.h>
#endif
#if defined(__AVX__)
    #include <immintrin.h>
#endif

namespace amal
{
    namespace internal
    {
#if defined(__SSE2__)
        inline __v4sf sqrt(__v4sf x) { return _mm_sqrt_ps(x); }

        inline __v4sf inverse_sqrt(__v4sf x) { return _mm_rsqrt_ps(x); }
#endif
#if defined(__AVX__)
        inline __v4df sqrt(__v4df x) { return _mm256_sqrt_pd(x); }

        inline __v4df inverse_sqrt(__v4df x)
        {
            __v4df sqrt_x = _mm256_sqrt_pd(x);
            __v4df one = _mm256_set1_pd(1.0);
            return _mm256_div_pd(one, sqrt_x);
        }
#endif
    } // namespace internal
} // namespace amal