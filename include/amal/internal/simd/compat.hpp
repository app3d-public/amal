#pragma once

#include <xmmintrin.h>
#if defined(__SSE4_1__) || defined(__AVX2__)
    #include <smmintrin.h>
#endif
#if defined(__AVX__)
    #include <immintrin.h>
#endif

namespace amal
{
    namespace internal
    {
        static inline __v4si mm_mullo_epi32_compat(__v4si a, __v4si b)
        {
#if defined(__SSE4_1__) || defined(__AVX2__)
            return _mm_mullo_epi32(a, b);
#else
            __m128i lo = _mm_mul_epu32(a, b);
            __m128i hi = _mm_mul_epu32(_mm_srli_si128(a, 4), _mm_srli_si128(b, 4));
            hi = _mm_shuffle_epi32(hi, _MM_SHUFFLE(0, 0, 2, 0));
            return _mm_unpacklo_epi64(lo, hi);
#endif
        }
    } // namespace internal
} // namespace amal