#pragma once

#include <cstdint>
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
        inline __v4sf abs(__v4sf x)
        {
            const __v4si abs_mask = {0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF};
            return reinterpret_cast<__v4sf>(reinterpret_cast<__v4si>(x) & abs_mask);
        }

        inline __v4sf mix(__v4sf x, __v4sf y, __v4si a)
        {
            __v4si mask = -a;
            __v4si xi = reinterpret_cast<__v4si &>(x);
            __v4si yi = reinterpret_cast<__v4si &>(y);
            __v4si result = (xi & ~mask) | (yi & mask);
            return reinterpret_cast<__v4sf &>(result);
        }

        inline __v4sf min(__v4sf a, __v4sf b) { return _mm_min_ps(a, b); }

        inline __v4sf max(__v4sf a, __v4sf b) { return _mm_max_ps(a, b); }

        inline __v4sf clamp(__v4sf x, __v4sf min_val, __v4sf max_val)
        {
            return _mm_min_ps(_mm_max_ps(x, min_val), max_val);
        }

        inline __v4si mix(__v4si x, __v4si y, __v4si a)
        {
            __v4si mask = -a;
            return (x & ~mask) | (y & mask);
        }

        inline float extract_scalar(__v4sf v) { return _mm_cvtss_f32(v); }

        inline int extract_scalar(__v4si v) { return _mm_cvtsi128_si32(v); }

    #if defined(__SSE4_1__)
        inline __v4sf round(__v4sf x) { return _mm_round_ps(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC); }

        inline __v4sf floor(__v4sf x) { return _mm_floor_ps(x); }

        inline __v4si abs(__v4si x)
        {
            __m128i xi = reinterpret_cast<__m128i &>(x);
            __m128i r = _mm_abs_epi32(xi);
            return reinterpret_cast<__v4si &>(r);
        }

        inline __v4sf trunc(__v4sf x) { return _mm_round_ps(x, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC); }

        inline __v4sf ceil(__v4sf x) { return _mm_ceil_ps(x); }

        inline __v4si min(__v4si a, __v4si b) { return _mm_min_epi32(a, b); }

        inline __v4si max(__v4si a, __v4si b) { return _mm_max_epi32(a, b); }
    #else
        inline __v4sf round(__v4sf x)
        {
            const __v4si sign_mask = {INT32_MIN, INT32_MIN, INT32_MIN, INT32_MIN};
            const __v4sf magic = _mm_set1_ps(8388608.0f);

            __v4sf sgn = reinterpret_cast<__v4sf>(reinterpret_cast<__v4si>(x) & sign_mask);
            __v4sf or0 = _mm_or_ps(sgn, magic);
            __v4sf add = _mm_add_ps(x, or0);
            return _mm_sub_ps(add, or0);
        }

        inline __v4sf floor(__v4sf x)
        {
            __v4sf r = round(x);
            __v4sf mask = _mm_cmplt_ps(x, r);
            __v4sf one = _mm_set1_ps(1.0f);
            __v4sf adj = _mm_and_ps(mask, one);
            return _mm_sub_ps(r, adj);
        }

        inline __v4si abs(__v4si x)
        {
            __v4si sign = x >> 31;
            return (x ^ sign) - sign;
        }

        inline __v4sf ceil(__v4sf x)
        {
            __v4sf r = round(x);
            __v4sf mask = _mm_cmplt_ps(r, x);
            __v4sf one = _mm_set1_ps(1.0f);
            __v4sf adj = _mm_and_ps(mask, one);
            return _mm_add_ps(r, adj);
        }

        inline __v4sf trunc(__v4sf x)
        {
            const __v4sf zero   = _mm_setzero_ps();
            __v4sf is_neg = _mm_cmplt_ps(x, zero);
            __v4sf floored = floor(x);
            __v4sf ceiled = ceil(x);
            return _mm_or_ps(_mm_and_ps(is_neg, ceiled), _mm_andnot_ps(is_neg, floored));
        }

        inline __v4si min(__v4si a, __v4si b)
        {
            __v4si mask = _mm_cmplt_epi32(a, b);
            return (a & mask) | (b & ~mask);
        }

        inline __v4si max(__v4si a, __v4si b)
        {
            __v4si mask = _mm_cmpgt_epi32(a, b);
            return (a & mask) | (b & ~mask);
        }
    #endif
#endif
#if defined(__AVX__)
        inline __v4df abs(__v4df x)
        {
            const __v4di abs_mask = {0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF};
            return reinterpret_cast<__v4df>(reinterpret_cast<__v4di>(x) & abs_mask);
        }

        inline __v4df mix(__v4df x, __v4df y, __v4si a)
        {
            __v4di mask64 = {static_cast<long long>(a[0]), static_cast<long long>(a[1]), static_cast<long long>(a[2]),
                             static_cast<long long>(a[3])};
            __v4di xi = reinterpret_cast<__v4di &>(x);
            __v4di yi = reinterpret_cast<__v4di &>(y);
            __v4di ri = (yi & mask64) | (xi & ~mask64);
            return reinterpret_cast<__v4df &>(ri);
        }

        inline double extract_scalar(__v4df v)
        {
            __v2df low = _mm256_castpd256_pd128(v);
            return _mm_cvtsd_f64(low);
        }

        inline __v4df round(__v4df x) { return _mm256_round_pd(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC); }

        inline __v4df floor(__v4df x) { return _mm256_floor_pd(x); }

        inline __v4df trunc(__v4df x) { return _mm256_round_pd(x, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC); }

        inline __v4df ceil(__v4df x) { return _mm256_ceil_pd(x); }

        inline __v4df min(__v4df a, __v4df b) { return _mm256_min_pd(a, b); }

        inline __v4df max(__v4df a, __v4df b) { return _mm256_max_pd(a, b); }
#endif
    } // namespace internal
} // namespace amal