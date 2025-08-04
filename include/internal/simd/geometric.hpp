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
        inline __v4sf cross_yzx(__v4sf x, __v4sf y)
        {
            const __v4sf x_yzx = _mm_shuffle_ps(x, x, _MM_SHUFFLE(3, 0, 2, 1));
            const __v4sf y_yzx = _mm_shuffle_ps(y, y, _MM_SHUFFLE(3, 0, 2, 1));
            const __v4sf mul1 = _mm_mul_ps(x, y_yzx);
            const __v4sf mul2 = _mm_mul_ps(x_yzx, y);
            return _mm_sub_ps(mul1, mul2);
        }

        inline __v4sf cross(__v4sf x, __v4sf y)
        {
            __v4sf c = cross_yzx(x, y);
            return _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 0, 2, 1));
        }

        inline __v4si cross_yzx(__v4si a, __v4si b)
        {
            return _mm_sub_epi32(_mm_mullo_epi32(_mm_shuffle_epi32(a, _MM_SHUFFLE(3, 0, 2, 1)),
                                                 _mm_shuffle_epi32(b, _MM_SHUFFLE(3, 1, 0, 2))),
                                 _mm_mullo_epi32(_mm_shuffle_epi32(a, _MM_SHUFFLE(3, 1, 0, 2)),
                                                 _mm_shuffle_epi32(b, _MM_SHUFFLE(3, 0, 2, 1))));
        }

        inline __v4si cross(__v4si a, __v4si b) { return _mm_shuffle_epi32(cross_yzx(a, b), _MM_SHUFFLE(3, 0, 2, 1)); }

        inline __v4si dot(__v4si a, __v4si b)
        {
            __v4si mul = _mm_mullo_epi32(a, b);
            __v4si temp = _mm_add_epi32(mul, _mm_shuffle_epi32(mul, _MM_SHUFFLE(2, 2, 0, 1)));
            int sum = _mm_cvtsi128_si32(temp) + _mm_extract_epi32(temp, 1);
            return _mm_set1_epi32(sum);
        }

    #if defined(__AVX__)
        inline __v4sf dot(__v4sf a, __v4sf b) { return _mm_dp_ps(a, b, 0xFF); }
    #elif defined(__SSE3__)
        inline __v4sf dot(__v4sf a, __v4sf b)
        {
            __v4sf mul = _mm_mul_ps(a, b);
            __v4sf hadd = _mm_hadd_ps(mul, mul);
            return _mm_hadd_ps(hadd, hadd);
        }
    #else
        inline __v4sf dot(__v4sf a, __v4sf b)
        {
            __v4sf mul = _mm_mul_ps(a, b);
            __v4sf temp = _mm_add_ps(mul, _mm_movehl_ps(mul, mul));
            __v4sf shuffled = _mm_shuffle_ps(temp, temp, 1);
            return _mm_add_ss(temp, shuffled);
        }
    #endif
#endif
#if defined(__AVX__)
        inline __v4df dot(__v4df a, __v4df b)
        {
            __v4df m = _mm256_mul_pd(a, b);
            __v4df h1 = _mm256_hadd_pd(m, m);
            return _mm256_hadd_pd(h1, h1);
        }

        inline __v4df cross_yzx(__v4df a, __v4df b)
        {
            const __v4df a_yzx = _mm256_permute4x64_pd(a, _MM_SHUFFLE(3, 0, 2, 1));
            const __v4df b_yzx = _mm256_permute4x64_pd(b, _MM_SHUFFLE(3, 0, 2, 1));
            return _mm256_sub_pd(_mm256_mul_pd(a, b_yzx), _mm256_mul_pd(a_yzx, b));
        }

        inline __v4df cross(__v4df a, __v4df b)
        {
            return _mm256_permute4x64_pd(cross_yzx(a, b), _MM_SHUFFLE(3, 0, 2, 1));
        }
#endif
    } // namespace internal
} // namespace amal