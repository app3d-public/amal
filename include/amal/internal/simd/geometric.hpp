#pragma once

#include "compat.hpp"

namespace amal
{
    namespace internal
    {
#if defined(__SSE2__)
        inline __m128_u cross_yzx(__m128_u x, __m128_u y)
        {
            const __m128_u x_yzx = _mm_shuffle_ps(x, x, _MM_SHUFFLE(3, 0, 2, 1));
            const __m128_u y_yzx = _mm_shuffle_ps(y, y, _MM_SHUFFLE(3, 0, 2, 1));
            const __m128_u mul1 = _mm_mul_ps(x, y_yzx);
            const __m128_u mul2 = _mm_mul_ps(x_yzx, y);
            return _mm_sub_ps(mul1, mul2);
        }

        inline __m128_u cross(__m128_u x, __m128_u y)
        {
            __m128_u c = cross_yzx(x, y);
            return _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 0, 2, 1));
        }

        inline __v4si_u cross_yzx(__v4si_u a, __v4si_u b)
        {
            return _mm_sub_epi32(mm_mullo_epi32_compat(_mm_shuffle_epi32(a, _MM_SHUFFLE(3, 0, 2, 1)),
                                                       _mm_shuffle_epi32(b, _MM_SHUFFLE(3, 1, 0, 2))),
                                 mm_mullo_epi32_compat(_mm_shuffle_epi32(a, _MM_SHUFFLE(3, 1, 0, 2)),
                                                       _mm_shuffle_epi32(b, _MM_SHUFFLE(3, 0, 2, 1))));
        }

        inline __v4si_u cross(__v4si_u a, __v4si_u b)
        {
            return _mm_shuffle_epi32(cross_yzx(a, b), _MM_SHUFFLE(3, 0, 2, 1));
        }

        inline __v4si_u dot(__v4si_u a, __v4si_u b)
        {
            __v4si_u mul = mm_mullo_epi32_compat(a, b);
    #if defined(__SSE4_1__) || defined(__AVX2__)
            __v4si_u temp = _mm_add_epi32(mul, _mm_shuffle_epi32(mul, _MM_SHUFFLE(2, 2, 0, 1)));
            int sum = _mm_cvtsi128_si32(temp) + _mm_extract_epi32(temp, 1);
    #else
            __v4si_u t1 = _mm_add_epi32(mul, _mm_srli_si128(mul, 8));
            __v4si_u t2 = _mm_add_epi32(t1, _mm_srli_si128(t1, 4));
            int sum = _mm_cvtsi128_si32(t2);
    #endif
            return _mm_set1_epi32(sum);
        }

        inline __m128_u dot(__m128_u a, __m128_u b)
        {
            // here I Remove _mm_dp_ps due to slow performance
    #if defined(AMAL_ARCH_SSE3)
            __m128 mul = _mm_mul_ps(a, b);     // [p0 p1 p2 p3]
            __m128 h1 = _mm_hadd_ps(mul, mul); // [p0+p1, p2+p3, p0+p1, p2+p3]
            __m128 h2 = _mm_hadd_ps(h1, h1);   // [sum,   sum,   sum,   sum]
            return h2;
    #else // SSE2
            __m128 p = _mm_mul_ps(a, b);                       // [p0 p1 p2 p3]
            __m128 t = _mm_add_ps(p, _mm_movehl_ps(p, p));     // [p0+p2, p1+p3, *, *]
            __m128 s = _mm_add_ss(t, _mm_shuffle_ps(t, t, 1)); // [sum, 0, 0, 0]
            return _mm_shuffle_ps(s, s, 0);                    // [sum sum sum sum]
    #endif
        }
#endif
#if defined(__AVX__)
        inline __m256d_u dot(__m256d_u a, __m256d_u b)
        {
            __m256d_u m = _mm256_mul_pd(a, b);
            __m256d_u h1 = _mm256_hadd_pd(m, m);
            return _mm256_hadd_pd(h1, h1);
        }

        inline __m256d_u cross_yzx(__m256d_u a, __m256d_u b)
        {
            const __m256d_u a_yzx = _mm256_permute4x64_pd(a, _MM_SHUFFLE(3, 0, 2, 1));
            const __m256d_u b_yzx = _mm256_permute4x64_pd(b, _MM_SHUFFLE(3, 0, 2, 1));
            return _mm256_sub_pd(_mm256_mul_pd(a, b_yzx), _mm256_mul_pd(a_yzx, b));
        }

        inline __m256d_u cross(__m256d_u a, __m256d_u b)
        {
            return _mm256_permute4x64_pd(cross_yzx(a, b), _MM_SHUFFLE(3, 0, 2, 1));
        }
#endif
    } // namespace internal
} // namespace amal