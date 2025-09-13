#pragma once

#include "../type_info.hpp"
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
        inline __m128_u abs(__m128_u x)
        {
            const __v4si_u abs_mask = {0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF};
            return reinterpret_cast<__m128_u>(reinterpret_cast<__v4si_u>(x) & abs_mask);
        }

        inline __m128_u fma(__m128_u const &x, __m128_u const &y, __m128_u const &z)
        {
            return _mm_add_ps(_mm_mul_ps(x, y), z);
        }

        template <int c>
        inline __m128_u splat(__m128_u x)
        {
            const int s = _MM_SHUFFLE(c, c, c, c);
    #if defined(__AVX__)
            return _mm_permute_ps(x, s);
    #else
            return _mm_shuffle_ps(x, x, s);
    #endif
        }

        inline __m128_u mix(__m128_u x, __m128_u y, __v4si_u a)
        {
            __v4si_u mask = -a;
            __v4si_u xi = reinterpret_cast<__v4si_u &>(x);
            __v4si_u yi = reinterpret_cast<__v4si_u &>(y);
            __v4si_u result = (xi & ~mask) | (yi & mask);
            return reinterpret_cast<__m128_u &>(result);
        }

        inline __m128_u min(__m128_u a, __m128_u b) { return _mm_min_ps(a, b); }

        inline __m128_u max(__m128_u a, __m128_u b) { return _mm_max_ps(a, b); }

        inline __m128_u clamp(__m128_u x, __m128_u min_val, __m128_u max_val)
        {
            return _mm_min_ps(_mm_max_ps(x, min_val), max_val);
        }

        inline __v4si_u mix(__v4si_u x, __v4si_u y, __v4si_u a)
        {
            __v4si_u mask = -a;
            return (x & ~mask) | (y & mask);
        }

        inline float extract_scalar(__m128_u v) { return _mm_cvtss_f32(v); }

        inline int extract_scalar(__v4si_u v) { return _mm_cvtsi128_si32(v); }

    #if defined(__SSE4_1__)
        inline __m128_u round(__m128_u x) { return _mm_round_ps(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC); }

        inline __m128_u floor(__m128_u x) { return _mm_floor_ps(x); }

        inline __v4si_u abs(__v4si_u x)
        {
            __m128i xi = reinterpret_cast<__m128i &>(x);
            __m128i r = _mm_abs_epi32(xi);
            return reinterpret_cast<__v4si_u &>(r);
        }

        inline __m128_u trunc(__m128_u x) { return _mm_round_ps(x, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC); }

        inline __m128_u ceil(__m128_u x) { return _mm_ceil_ps(x); }

        inline __v4si_u min(__v4si_u a, __v4si_u b) { return _mm_min_epi32(a, b); }

        inline __v4si_u max(__v4si_u a, __v4si_u b) { return _mm_max_epi32(a, b); }
    #else
        inline __m128_u round(__m128_u x)
        {
            const __v4si_u sign_mask = {INT32_MIN, INT32_MIN, INT32_MIN, INT32_MIN};
            const __m128_u magic = _mm_set1_ps(8388608.0f);

            __m128_u sgn = reinterpret_cast<__m128_u>(reinterpret_cast<__v4si_u>(x) & sign_mask);
            __m128_u or0 = _mm_or_ps(sgn, magic);
            __m128_u add = _mm_add_ps(x, or0);
            return _mm_sub_ps(add, or0);
        }

        inline __m128_u floor(__m128_u x)
        {
            __m128_u r = round(x);
            __m128_u mask = _mm_cmplt_ps(x, r);
            __m128_u one = _mm_set1_ps(1.0f);
            __m128_u adj = _mm_and_ps(mask, one);
            return _mm_sub_ps(r, adj);
        }

        inline __v4si_u abs(__v4si_u x)
        {
            __v4si_u sign = x >> 31;
            return (x ^ sign) - sign;
        }

        inline __m128_u ceil(__m128_u x)
        {
            __m128_u r = round(x);
            __m128_u mask = _mm_cmplt_ps(r, x);
            __m128_u one = _mm_set1_ps(1.0f);
            __m128_u adj = _mm_and_ps(mask, one);
            return _mm_add_ps(r, adj);
        }

        inline __m128_u trunc(__m128_u x)
        {
            const __m128_u zero = _mm_setzero_ps();
            __m128_u is_neg = _mm_cmplt_ps(x, zero);
            __m128_u floored = floor(x);
            __m128_u ceiled = ceil(x);
            return _mm_or_ps(_mm_and_ps(is_neg, ceiled), _mm_andnot_ps(is_neg, floored));
        }

        inline __v4si_u min(__v4si_u a, __v4si_u b)
        {
            __v4si_u mask = _mm_cmplt_epi32(a, b);
            return (a & mask) | (b & ~mask);
        }

        inline __v4si_u max(__v4si_u a, __v4si_u b)
        {
            __v4si_u mask = _mm_cmpgt_epi32(a, b);
            return (a & mask) | (b & ~mask);
        }
    #endif
#endif
#if defined(__AVX__)
        inline __m256d_u abs(__m256d_u x)
        {
            const __v4di abs_mask = {0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF};
            return reinterpret_cast<__m256d_u>(reinterpret_cast<__v4di>(x) & abs_mask);
        }

        inline __m256d_u fma(__m256d_u const &x, __m256d_u const &y, __m256d_u const &z)
        {
            return _mm256_fmadd_pd(x, y, z);
        }

        template <int c>
        inline __m256d_u splat(__m256d_u x)
        {
            const int s = _MM_SHUFFLE(c, c, c, c);
            return _mm256_permute_ps(x, s);
        }

        inline __m256d_u mix(__m256d_u x, __m256d_u y, __v4si_u a)
        {
            __v4di mask64 = {static_cast<long long>(a[0]), static_cast<long long>(a[1]), static_cast<long long>(a[2]),
                             static_cast<long long>(a[3])};
            __v4di xi = reinterpret_cast<__v4di &>(x);
            __v4di yi = reinterpret_cast<__v4di &>(y);
            __v4di ri = (yi & mask64) | (xi & ~mask64);
            return reinterpret_cast<__m256d_u &>(ri);
        }

        inline double extract_scalar(__m256d_u v)
        {
            __v2df low = _mm256_castpd256_pd128(v);
            return _mm_cvtsd_f64(low);
        }

        inline __m256d_u round(__m256d_u x)
        {
            return _mm256_round_pd(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
        }

        inline __m256d_u floor(__m256d_u x) { return _mm256_floor_pd(x); }

        inline __m256d_u trunc(__m256d_u x) { return _mm256_round_pd(x, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC); }

        inline __m256d_u ceil(__m256d_u x) { return _mm256_ceil_pd(x); }

        inline __m256d_u min(__m256d_u a, __m256d_u b) { return _mm256_min_pd(a, b); }

        inline __m256d_u max(__m256d_u a, __m256d_u b) { return _mm256_max_pd(a, b); }
#endif
    } // namespace internal
} // namespace amal