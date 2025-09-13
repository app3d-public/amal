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
        inline __m128_u sqrt(__m128_u x) { return _mm_sqrt_ps(x); }

        inline __m128_u inverse_sqrt(__m128_u x) { return _mm_rsqrt_ps(x); }
#endif
#if defined(__AVX__)
        inline __m256d_u sqrt(__m256d_u x) { return _mm256_sqrt_pd(x); }

        inline __m256d_u inverse_sqrt(__m256d_u x)
        {
            __m256d_u sqrt_x = _mm256_sqrt_pd(x);
            __m256d_u one = _mm256_set1_pd(1.0);
            return _mm256_div_pd(one, sqrt_x);
        }
#endif
    } // namespace internal
} // namespace amal