#pragma once

#include <xmmintrin.h>
#include "../type_info.hpp"

#if defined(__SSE4_1__) || defined(__AVX2__)
    #include <smmintrin.h>
#endif
#if defined(__AVX__)
    #include <immintrin.h>
#endif

#if defined(AMAL_FMA_ENABLE)
    #define AMAL_FMA_ADD(a, b, c) _mm_fmadd_ps((a), (b), (c))
    #define AMAL_FMA_SUB(a, b, c) _mm_fmsub_ps((a), (b), (c))
    #if defined(__AVX__)
        #define AMAL_FMA_ADD_PD(a, b, c) _mm256_fmadd_pd((a), (b), (c))
        #define AMAL_FMA_SUB_PD(a, b, c) _mm256_fmsub_pd((a), (b), (c))
    #endif
#else
    #define AMAL_FMA_ADD(a, b, c) _mm_add_ps(_mm_mul_ps((a), (b)), (c))
    #define AMAL_FMA_SUB(a, b, c) _mm_sub_ps(_mm_mul_ps((a), (b)), (c))
    #if defined(__AVX__)
        #define AMAL_FMA_ADD_PD(a, b, c) _mm256_add_pd(_mm256_mul_pd((a), (b)), (c))
        #define AMAL_FMA_SUB_PD(a, b, c) _mm256_sub_pd(_mm256_mul_pd((a), (b)), (c))
    #endif
#endif

namespace amal::internal
{
    static inline __v4si_u mm_mullo_epi32_compat(__v4si_u a, __v4si_u b)
    {
#if defined(__SSE4_1__) || defined(__AVX2__)
        return AMAL_V4SI(_mm_mullo_epi32(AMAL_M128I(a), AMAL_M128I(b)));
#else
        __m128i a128 = AMAL_M128I(a);
        __m128i b128 = AMAL_M128I(b);
        __m128i lo = _mm_mul_epu32(a128, b128);
        __m128i hi = _mm_mul_epu32(_mm_srli_si128(a128, 4), _mm_srli_si128(b128, 4));
        hi = _mm_shuffle_epi32(hi, _MM_SHUFFLE(0, 0, 2, 0));
        return AMAL_V4SI(_mm_unpacklo_epi64(lo, hi));
#endif
    }
} // namespace amal::internal
