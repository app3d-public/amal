#pragma once

#if defined(__SSE2__)
    #include <xmmintrin.h>
#endif
#if defined(__AVX__)
    #include <immintrin.h>
#endif
#include "compat.hpp"

namespace amal
{
    namespace internal
    {
#if defined(__SSE2__)
        inline __m128_u log2_poly_1pf(__m128_u f)
        {
            const __m128 c1 = _mm_set1_ps(1.4426899f);
            const __m128 c2 = _mm_set1_ps(-0.7211659f);
            const __m128 c3 = _mm_set1_ps(0.4808983f);
            const __m128 c4 = _mm_set1_ps(-0.36067376f);
            const __m128 c5 = _mm_set1_ps(0.28853902f);

            __m128 p = AMAL_FMA_ADD(c5, f, c4);
            p = AMAL_FMA_ADD(p, f, c3);
            p = AMAL_FMA_ADD(p, f, c2);
            p = AMAL_FMA_ADD(p, f, c1);
            return _mm_mul_ps(f, p);
        }

        inline __m128_u exp2_poly_fracf(__m128_u f)
        {
            const __m128 k1 = _mm_set1_ps(0.69314718f);
            const __m128 k2 = _mm_set1_ps(0.24022651f);
            const __m128 k3 = _mm_set1_ps(0.05550411f);
            const __m128 k4 = _mm_set1_ps(0.00961813f);

            __m128 p = AMAL_FMA_ADD(k4, f, k3);
            p = AMAL_FMA_ADD(p, f, k2);
            p = AMAL_FMA_ADD(p, f, k1);
            return AMAL_FMA_ADD(f, p, _mm_set1_ps(1.0f));
        }

        inline __m128_u sqrt(__m128_u x) { return _mm_sqrt_ps(x); }

        inline __m128_u inverse_sqrt(__m128_u x) { return _mm_rsqrt_ps(x); }

        inline __m128_u log2(__m128_u x)
        {
            // Keep strictly positive domain for log2.
            const __m128 min_norm = _mm_set1_ps(1.17549435e-38f);
            x = _mm_max_ps(x, min_norm);

            __m128i ix = _mm_castps_si128(x);
            __m128i exp_i = _mm_sub_epi32(_mm_and_si128(_mm_srli_epi32(ix, 23), _mm_set1_epi32(255)),
                                          _mm_set1_epi32(127));
            __m128 exp_f = _mm_cvtepi32_ps(exp_i);

            __m128i mant_i = _mm_or_si128(_mm_and_si128(ix, _mm_set1_epi32(0x007FFFFF)), _mm_set1_epi32(0x3F800000));
            __m128 m = _mm_castsi128_ps(mant_i);
            __m128 f = _mm_sub_ps(m, _mm_set1_ps(1.0f));
            return _mm_add_ps(exp_f, log2_poly_1pf(f));
        }

        inline __m128_u exp2(__m128_u x)
        {
            // Clamp to float exponent range.
            x = _mm_min_ps(_mm_max_ps(x, _mm_set1_ps(-126.0f)), _mm_set1_ps(127.0f));

            __m128 n_f;
#if defined(__SSE4_1__)
            // Fast path when SSE4.1 is available.
            n_f = _mm_floor_ps(x);
#else
            // SSE2 fallback: floor via truncation + correction for negative fractional values.
            n_f = _mm_cvtepi32_ps(_mm_cvttps_epi32(x));
            __m128 mask = _mm_cmplt_ps(x, n_f);
            n_f = _mm_sub_ps(n_f, _mm_and_ps(mask, _mm_set1_ps(1.0f)));
#endif
            __m128 f = _mm_sub_ps(x, n_f);

            __m128i n_i = _mm_cvttps_epi32(n_f);
            __m128i exp_bits = _mm_slli_epi32(_mm_add_epi32(n_i, _mm_set1_epi32(127)), 23);
            __m128 two_pow_n = _mm_castsi128_ps(exp_bits);

            return _mm_mul_ps(two_pow_n, exp2_poly_fracf(f));
        }

        inline __m128_u pow(__m128_u base, __m128_u exponent) { return exp2(_mm_mul_ps(log2(base), exponent)); }

        inline __m128_u log(__m128_u x)
        {
            static constexpr float inv_log2_e = 0.69314718055994530942f;
            return _mm_mul_ps(log2(x), _mm_set1_ps(inv_log2_e));
        }

        inline __m128_u log10(__m128_u x)
        {
            static constexpr float inv_log2_10 = 0.30102999566398119521f;
            return _mm_mul_ps(log2(x), _mm_set1_ps(inv_log2_10));
        }

        inline __m128_u exp(__m128_u x)
        {
            static constexpr float log2_e = 1.44269504088896340736f;
            return exp2(_mm_mul_ps(x, _mm_set1_ps(log2_e)));
        }
#endif
#if defined(__AVX__)
        inline __m128d fmadd_pd128(__m128d a, __m128d b, __m128d c)
        {
#if defined(AMAL_FMA_ENABLE)
            return _mm_fmadd_pd(a, b, c);
#else
            return _mm_add_pd(_mm_mul_pd(a, b), c);
#endif
        }

        inline __m128d log2_pd_sse2(__m128d x)
        {
            const __m128d min_norm = _mm_set1_pd(2.2250738585072014e-308);
            x = _mm_max_pd(x, min_norm);

            __m128i ix = _mm_castpd_si128(x);
            __m128i exp_i64 = _mm_sub_epi64(_mm_and_si128(_mm_srli_epi64(ix, 52), _mm_set1_epi64x(0x7FF)),
                                            _mm_set1_epi64x(1023));
            __m128i exp_i32_packed = _mm_shuffle_epi32(exp_i64, _MM_SHUFFLE(0, 0, 2, 0));
            __m128d exp_f = _mm_cvtepi32_pd(exp_i32_packed);

            __m128i mant_i =
                _mm_or_si128(_mm_and_si128(ix, _mm_set1_epi64x(0x000FFFFFFFFFFFFF)), _mm_set1_epi64x(0x3FF0000000000000));
            __m128d m = _mm_castsi128_pd(mant_i);
            __m128d f = _mm_sub_pd(m, _mm_set1_pd(1.0));

            const __m128d c1 = _mm_set1_pd(1.4426950408889634074);
            const __m128d c2 = _mm_set1_pd(-0.7213475204444817037);
            const __m128d c3 = _mm_set1_pd(0.4808983469629878025);
            const __m128d c4 = _mm_set1_pd(-0.3606737602222408518);
            const __m128d c5 = _mm_set1_pd(0.2885390081777926815);
            const __m128d c6 = _mm_set1_pd(-0.2404491734814939012);

            __m128d p = fmadd_pd128(c6, f, c5);
            p = fmadd_pd128(p, f, c4);
            p = fmadd_pd128(p, f, c3);
            p = fmadd_pd128(p, f, c2);
            p = fmadd_pd128(p, f, c1);
            __m128d log2_m = _mm_mul_pd(f, p);
            return _mm_add_pd(exp_f, log2_m);
        }

        inline __m128d exp2_pd(__m128d x)
        {
            x = _mm_min_pd(_mm_max_pd(x, _mm_set1_pd(-1022.0)), _mm_set1_pd(1023.0));

            __m128d n = _mm_floor_pd(x);
            __m128d f = _mm_sub_pd(x, n);

            const __m128d ln2 = _mm_set1_pd(0.6931471805599453094);
            const __m128d k1 = ln2;
            const __m128d k2 = _mm_set1_pd(0.2402265069591007123);
            const __m128d k3 = _mm_set1_pd(0.05550410866482157995);
            const __m128d k4 = _mm_set1_pd(0.009618129107628477161);
            const __m128d k5 = _mm_set1_pd(0.001333355814642844342);
            const __m128d k6 = _mm_set1_pd(0.0001540353039338161124);

            __m128d p = fmadd_pd128(k6, f, k5);
            p = fmadd_pd128(p, f, k4);
            p = fmadd_pd128(p, f, k3);
            p = fmadd_pd128(p, f, k2);
            p = fmadd_pd128(p, f, k1);
            __m128d two_pow_f = fmadd_pd128(f, p, _mm_set1_pd(1.0));

            __m128i n_i32 = _mm_cvttpd_epi32(n);
            __m128i sign = _mm_srai_epi32(n_i32, 31);
            __m128i n_i64 = _mm_unpacklo_epi32(n_i32, sign);
            __m128i exp_i64 = _mm_add_epi64(n_i64, _mm_set_epi32(0, 1023, 0, 1023));
            __m128i exp_bits = _mm_slli_epi64(exp_i64, 52);
            __m128d two_pow_n = _mm_castsi128_pd(exp_bits);
            return _mm_mul_pd(two_pow_n, two_pow_f);
        }

        inline __m256d_u log2(__m256d_u x)
        {
#if defined(__AVX2__)
            const __m256d min_norm = _mm256_set1_pd(2.2250738585072014e-308);
            x = _mm256_max_pd(x, min_norm);

            __m256i ix = _mm256_castpd_si256(x);
            __m256i exp_i64 =
                _mm256_sub_epi64(_mm256_and_si256(_mm256_srli_epi64(ix, 52), _mm256_set1_epi64x(0x7FF)),
                                 _mm256_set1_epi64x(1023));

            __m128i lo_i64 = _mm256_castsi256_si128(exp_i64);
            __m128i hi_i64 = _mm256_extractf128_si256(exp_i64, 1);
            __m128d lo_e = _mm_cvtepi32_pd(_mm_shuffle_epi32(lo_i64, _MM_SHUFFLE(0, 0, 2, 0)));
            __m128d hi_e = _mm_cvtepi32_pd(_mm_shuffle_epi32(hi_i64, _MM_SHUFFLE(0, 0, 2, 0)));
            __m256d exp_f = _mm256_castpd128_pd256(lo_e);
            exp_f = _mm256_insertf128_pd(exp_f, hi_e, 1);

            __m256i mant_i = _mm256_or_si256(_mm256_and_si256(ix, _mm256_set1_epi64x(0x000FFFFFFFFFFFFF)),
                                             _mm256_set1_epi64x(0x3FF0000000000000));
            __m256d m = _mm256_castsi256_pd(mant_i);
            __m256d f = _mm256_sub_pd(m, _mm256_set1_pd(1.0));

            const __m256d c1 = _mm256_set1_pd(1.4426950408889634074);
            const __m256d c2 = _mm256_set1_pd(-0.7213475204444817037);
            const __m256d c3 = _mm256_set1_pd(0.4808983469629878025);
            const __m256d c4 = _mm256_set1_pd(-0.3606737602222408518);
            const __m256d c5 = _mm256_set1_pd(0.2885390081777926815);
            const __m256d c6 = _mm256_set1_pd(-0.2404491734814939012);

            __m256d p = AMAL_FMA_ADD_PD(c6, f, c5);
            p = AMAL_FMA_ADD_PD(p, f, c4);
            p = AMAL_FMA_ADD_PD(p, f, c3);
            p = AMAL_FMA_ADD_PD(p, f, c2);
            p = AMAL_FMA_ADD_PD(p, f, c1);
            __m256d log2_m = _mm256_mul_pd(f, p);
            return _mm256_add_pd(exp_f, log2_m);
#else
            __m128d lo = _mm256_castpd256_pd128(x);
            __m128d hi = _mm256_extractf128_pd(x, 1);
            __m256d_u out = _mm256_castpd128_pd256(log2_pd_sse2(lo));
            return _mm256_insertf128_pd(out, log2_pd_sse2(hi), 1);
#endif
        }

        inline __m256d_u exp2(__m256d_u x)
        {
#if defined(__AVX2__)
            x = _mm256_min_pd(_mm256_max_pd(x, _mm256_set1_pd(-1022.0)), _mm256_set1_pd(1023.0));

            __m256d n = _mm256_cvtepi32_pd(_mm256_cvttpd_epi32(x));
            __m256d mask = _mm256_cmp_pd(x, n, _CMP_LT_OQ);
            n = _mm256_sub_pd(n, _mm256_and_pd(mask, _mm256_set1_pd(1.0)));
            __m256d f = _mm256_sub_pd(x, n);

            const __m256d k1 = _mm256_set1_pd(0.6931471805599453094);
            const __m256d k2 = _mm256_set1_pd(0.2402265069591007123);
            const __m256d k3 = _mm256_set1_pd(0.05550410866482157995);
            const __m256d k4 = _mm256_set1_pd(0.009618129107628477161);
            const __m256d k5 = _mm256_set1_pd(0.001333355814642844342);
            const __m256d k6 = _mm256_set1_pd(0.0001540353039338161124);

            __m256d p = AMAL_FMA_ADD_PD(k6, f, k5);
            p = AMAL_FMA_ADD_PD(p, f, k4);
            p = AMAL_FMA_ADD_PD(p, f, k3);
            p = AMAL_FMA_ADD_PD(p, f, k2);
            p = AMAL_FMA_ADD_PD(p, f, k1);
            __m256d two_pow_f = AMAL_FMA_ADD_PD(f, p, _mm256_set1_pd(1.0));

            __m128d n_lo = _mm256_castpd256_pd128(n);
            __m128d n_hi = _mm256_extractf128_pd(n, 1);

            __m128i n_lo_i32 = _mm_cvttpd_epi32(n_lo);
            __m128i n_hi_i32 = _mm_cvttpd_epi32(n_hi);
            __m128i lo_sign = _mm_srai_epi32(n_lo_i32, 31);
            __m128i hi_sign = _mm_srai_epi32(n_hi_i32, 31);
            __m128i lo_i64 = _mm_unpacklo_epi32(n_lo_i32, lo_sign);
            __m128i hi_i64 = _mm_unpacklo_epi32(n_hi_i32, hi_sign);

            __m256i n_i64 = _mm256_castsi128_si256(lo_i64);
            n_i64 = _mm256_insertf128_si256(n_i64, hi_i64, 1);
            __m256i exp_i64 = _mm256_add_epi64(n_i64, _mm256_set1_epi64x(1023));
            __m256i exp_bits = _mm256_slli_epi64(exp_i64, 52);
            __m256d two_pow_n = _mm256_castsi256_pd(exp_bits);
            return _mm256_mul_pd(two_pow_n, two_pow_f);
#else
            // AVX without AVX2: evaluate two SSE2 vectors.
            __m128d lo = _mm256_castpd256_pd128(x);
            __m128d hi = _mm256_extractf128_pd(x, 1);
            __m256d_u out = _mm256_castpd128_pd256(exp2_pd(lo));
            return _mm256_insertf128_pd(out, exp2_pd(hi), 1);
#endif
        }

        inline __m256d_u pow(__m256d_u base, __m256d_u exponent) { return exp2(_mm256_mul_pd(log2(base), exponent)); }

        inline __m256d_u sqrt(__m256d_u x) { return _mm256_sqrt_pd(x); }

        inline __m256d_u inverse_sqrt(__m256d_u x)
        {
            __m256d_u sqrt_x = _mm256_sqrt_pd(x);
            __m256d_u one = _mm256_set1_pd(1.0);
            return _mm256_div_pd(one, sqrt_x);
        }

        inline __m256d_u log(__m256d_u x)
        {
            static constexpr double inv_log2_e = 0.69314718055994530942;
            return _mm256_mul_pd(log2(x), _mm256_set1_pd(inv_log2_e));
        }

        inline __m256d_u log10(__m256d_u x)
        {
            static constexpr double inv_log2_10 = 0.30102999566398119521;
            return _mm256_mul_pd(log2(x), _mm256_set1_pd(inv_log2_10));
        }

        inline __m256d_u exp(__m256d_u x)
        {
            static constexpr double log2_e = 1.44269504088896340736;
            return exp2(_mm256_mul_pd(x, _mm256_set1_pd(log2_e)));
        }
#endif
    } // namespace internal
} // namespace amal
