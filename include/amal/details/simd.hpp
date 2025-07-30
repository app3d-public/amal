#pragma once

#if defined(__SSE2__)
    #include <emmintrin.h>
    #include <xmmintrin.h>

#endif
#if defined(__AVX__)
    #include <immintrin.h>
#endif

namespace amal
{
    namespace details
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

        inline __v4si less_than_mask(__v4sf x, __v4sf y)
        {
            __m128 cmp = _mm_cmplt_ps(x, y);
            return reinterpret_cast<__v4si &>(cmp);
        }

        inline __v4si less_than_equal_mask(__v4sf x, __v4sf y)
        {
            __m128 cmp = _mm_cmple_ps(x, y);
            return reinterpret_cast<__v4si &>(cmp);
        }

        inline __v4si greater_than_mask(__v4sf x, __v4sf y)
        {
            __m128 cmp = _mm_cmpgt_ps(x, y);
            return reinterpret_cast<__v4si &>(cmp);
        }

        inline __v4si greater_than_equal_mask(__v4sf x, __v4sf y)
        {
            __m128 cmp = _mm_cmpge_ps(x, y);
            return reinterpret_cast<__v4si &>(cmp);
        }

        inline __v4si equal_mask(__v4sf x, __v4sf y)
        {
            __m128 cmp = _mm_cmpeq_ps(x, y);
            return reinterpret_cast<__v4si &>(cmp);
        }

        inline __v4si not_equal_mask(__v4sf x, __v4sf y)
        {
            __m128 cmp = _mm_cmpneq_ps(x, y);
            return reinterpret_cast<__v4si &>(cmp);
        }

        inline __v4sf min(__v4sf a, __v4sf b) { return _mm_min_ps(a, b); }

        inline __v4sf max(__v4sf a, __v4sf b) { return _mm_max_ps(a, b); }

        inline __v4sf clamp(__v4sf x, __v4sf min_val, __v4sf max_val)
        {
            return _mm_min_ps(_mm_max_ps(x, min_val), max_val);
        }

        inline __v4sf sqrt(__v4sf x) { return _mm_sqrt_ps(x); }

        inline __v4sf inverse_sqrt(__v4sf x) { return _mm_rsqrt_ps(x); }

        inline __v4sf cross(__v4sf x, __v4sf y)
        {
            const __v4sf x_yzx = _mm_shuffle_ps(x, x, _MM_SHUFFLE(3, 0, 2, 1));
            const __v4sf y_yzx = _mm_shuffle_ps(y, y, _MM_SHUFFLE(3, 0, 2, 1));
            const __v4sf mul1 = _mm_mul_ps(x, y_yzx);
            const __v4sf mul2 = _mm_mul_ps(x_yzx, y);
            const __v4sf result = _mm_sub_ps(mul1, mul2);
            return _mm_shuffle_ps(result, result, _MM_SHUFFLE(3, 0, 2, 1));
        }

    #if defined(__SSE4_1__)
        inline __v4sf round(__v4sf x) { return _mm_round_ps(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC); }

        inline __v4sf floor(__v4sf x) { return _mm_floor_ps(x); }

        inline __v4si abs(__v4si x)
        {
            __m128i xi = reinterpret_cast<__m128i &>(x);
            __m128i r = _mm_abs_epi32(xi);
            return reinterpret_cast<__v4si &>(r);
        }

        inline __v4si less_than_mask(__v4si x, __v4si y)
        {
            __m128i xi = reinterpret_cast<__m128i &>(x);
            __m128i yi = reinterpret_cast<__m128i &>(y);
            __m128i cmp = _mm_cmplt_epi32(xi, yi);
            return reinterpret_cast<__v4si &>(cmp);
        }

        inline __v4si less_than_equal_mask(__v4si x, __v4si y)
        {
            __m128i xi = reinterpret_cast<__m128i &>(x);
            __m128i yi = reinterpret_cast<__m128i &>(y);
            __m128i lt = _mm_cmplt_epi32(xi, yi);
            __m128i eq = _mm_cmpeq_epi32(xi, yi);
            __m128i le = _mm_or_si128(lt, eq);
            return reinterpret_cast<__v4si &>(le);
        }

        inline __v4si greater_than_equal_mask(__v4si x, __v4si y)
        {
            __m128i cmp = _mm_cmplt_epi32(reinterpret_cast<__m128i &>(x), reinterpret_cast<__m128i &>(y));
            __m128i inv_cmp = _mm_andnot_si128(cmp, _mm_set1_epi32(-1));
            return reinterpret_cast<__v4si &>(inv_cmp);
        }

        inline __v4si greater_than_mask(__v4si x, __v4si y)
        {
            __m128i xi = reinterpret_cast<__m128i &>(x);
            __m128i yi = reinterpret_cast<__m128i &>(y);
            __m128i cmp = _mm_cmpgt_epi32(xi, yi);
            return reinterpret_cast<__v4si &>(cmp);
        }

        inline __v4si equal_mask(__v4si x, __v4si y)
        {
            __m128i cmp = _mm_cmpeq_epi32(reinterpret_cast<__m128i &>(x), reinterpret_cast<__m128i &>(y));
            return reinterpret_cast<__v4si &>(cmp);
        }

        inline __v4sf trunc(__v4sf x) { return _mm_round_ps(x, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC); }

        inline __v4sf ceil(__v4sf x) { return _mm_ceil_ps(x); }

        inline __v4si min(__v4si a, __v4si b) { return _mm_min_epi32(a, b); }

        inline __v4si max(__v4si a, __v4si b) { return _mm_max_epi32(a, b); }
    #elif defined(__SSE3__)
        inline __v4sf dot(__v4sf a, __v4sf b)
        {
            __v4sf mul = _mm_mul_ps(a, b);
            __v4sf hadd = _mm_hadd_ps(mul, mul);
            return _mm_hadd_ps(hadd, hadd);
        }
    #else
        inline __v4sf round(__v4sf x)
        {
            const __v4si sign_mask = {0x80000000, 0x80000000, 0x80000000, 0x80000000};
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

        inline __v4si less_than_mask(__v4si x, __v4si y) { return ((x - y) >> 31); }

        inline __v4si less_than_equal_mask(__v4si x, __v4si y)
        {
            __v4si lt = ((x - y) >> 31);
            __v4si eq = ~(x ^ y);
            return lt | eq;
        }

        inline __v4si greater_than_equal_mask(__v4si x, __v4si y) { return ~((x - y) >> 31); }
        inline __v4si greater_than_mask(__v4si x, __v4si y) { return ((y - x) >> 31); }

        inline __v4si equal_mask(__v4si x, __v4si y)
        {
            __v4si diff = x ^ y;
            return ~(diff | -diff) >> 31;
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
            __v4sf is_neg = _mm_cmplt_ps(x, 0.0f);
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

        inline __v4sf dot(__v4sf a, __v4sf b)
        {
            __v4sf mul = _mm_mul_ps(a, b);
            __v4sf temp = _mm_add_ps(mul, _mm_movehl_ps(mul, mul));
            __v4sf shuffled = _mm_shuffle_ps(temp, temp, 1);
            return _mm_add_ss(temp, shuffled);
        }
    #endif

        inline __v4si mix(__v4si x, __v4si y, __v4si a)
        {
            __v4si mask = -a;
            return (x & ~mask) | (y & mask);
        }

        inline __v4si not_equal_mask(__v4si x, __v4si y) { return (x ^ y); }
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

        inline __v4df sqrt(__v4df x) { return _mm256_sqrt_pd(x); }

        inline __v4df inverse_sqrt(__v4df x)
        {
            __v4df sqrt_x = _mm256_sqrt_pd(x);
            __v4df one = _mm256_set1_pd(1.0);
            return _mm256_div_pd(one, sqrt_x);
        }

        inline __v4sf dot(__v4sf a, __v4sf b) { return _mm_dp_ps(a, b, 0xFF); }

        inline __v4df dot(__v4df a, __v4df b)
        {
            __v4df mul = _mm256_mul_pd(a, b);
            __v2df lo = _mm256_castpd256_pd128(mul);
            __v2df hi = _mm256_extractf128_pd(mul, 1);
            __v2df sum2 = _mm_add_pd(lo, hi);
            __v2df swapped = _mm_shuffle_pd(sum2, sum2, 0x1);
            __v2df dot128 = _mm_add_sd(sum2, swapped);
            return _mm256_broadcast_sd(reinterpret_cast<const double *>(&dot128));
        }

        inline double extract_scalar(__v4df v)
        {
            __v2df low = _mm256_castpd256_pd128(v);
            return _mm_cvtsd_f64(low);
        }

        inline __v4df cross(__v4df a, __v4df b)
        {
            const __m256d a_yzx = _mm256_permute4x64_pd(a, _MM_SHUFFLE(3, 0, 2, 1));
            const __m256d b_yzx = _mm256_permute4x64_pd(b, _MM_SHUFFLE(3, 0, 2, 1));
            const __m256d c = _mm256_sub_pd(_mm256_mul_pd(a, b_yzx), _mm256_mul_pd(a_yzx, b));
            return _mm256_permute4x64_pd(c, _MM_SHUFFLE(3, 0, 2, 1));
        }

    #if defined(__AVX2__)
        inline __v4si less_than_mask(__v4df x, __v4df y)
        {
            __m256d cmp = _mm256_cmp_pd(x, y, _CMP_LT_OQ);
            __m256i cmp_i = _mm256_castpd_si256(cmp);
            __m128i lower = _mm256_castsi256_si128(cmp_i);
            return reinterpret_cast<__v4si &>(lower);
        }

        inline __v4si less_than_equal_mask(__v4df x, __v4df y)
        {
            __m256d cmp = _mm256_cmp_pd(x, y, _CMP_LE_OQ);
            __m256i cmp_i = _mm256_castpd_si256(cmp);
            __m128i lower = _mm256_castsi256_si128(cmp_i);
            return reinterpret_cast<__v4si &>(lower);
        }

        inline __v4si greater_than_mask(__v4df x, __v4df y)
        {
            __m256d cmp = _mm256_cmp_pd(x, y, _CMP_GT_OQ);
            __m256i cmp_i = _mm256_castpd_si256(cmp);
            __m128i lower = _mm256_castsi256_si128(cmp_i);
            return reinterpret_cast<__v4si &>(lower);
        }

        inline __v4si greater_than_equal_mask(__v4df x, __v4df y)
        {
            __m256d cmp = _mm256_cmp_pd(x, y, _CMP_GE_OQ);
            __m256i cmp_i = _mm256_castpd_si256(cmp);
            __m128i lower = _mm256_castsi256_si128(cmp_i);
            return reinterpret_cast<__v4si &>(lower);
        }

        inline __v4si equal_mask(__v4df x, __v4df y)
        {
            __m256d cmp = _mm256_cmp_pd(x, y, _CMP_EQ_OQ);
            __m256i cmp_i = _mm256_castpd_si256(cmp);
            __m128i lower = _mm256_castsi256_si128(cmp_i);
            return reinterpret_cast<__v4si &>(lower);
        }

        inline __v4si not_equal_mask(__v4df x, __v4df y)
        {
            __m256d cmp = _mm256_cmp_pd(x, y, _CMP_NEQ_OQ);
            __m256i cmp_i = _mm256_castpd_si256(cmp);
            __m128i lower = _mm256_castsi256_si128(cmp_i);
            return reinterpret_cast<__v4si &>(lower);
        }
    #else
        inline __v4si less_than_mask(__v4df x, __v4df y)
        {
            __m256d cmp_pd = _mm256_cmp_pd(x, y, _CMP_LT_OQ);
            int mask = _mm256_movemask_pd(cmp_pd);
            return __v4si{(mask & 0x1) ? -1 : 0, (mask & 0x2) ? -1 : 0, (mask & 0x4) ? -1 : 0, (mask & 0x8) ? -1 : 0};
        }

        inline __v4si less_than_equal_mask(__v4df x, __v4df y)
        {
            __m256d cmp = _mm256_cmp_pd(x, y, _CMP_LE_OQ);
            int mask = _mm256_movemask_pd(cmp);
            return __v4si{(mask & 0x1) ? -1 : 0, (mask & 0x2) ? -1 : 0, (mask & 0x4) ? -1 : 0, (mask & 0x8) ? -1 : 0};
        }

        inline __v4si greater_than_mask(__v4df x, __v4df y)
        {
            __m256d cmp = _mm256_cmp_pd(x, y, _CMP_GT_OQ);
            int mask = _mm256_movemask_pd(cmp);
            return __v4si{(mask & 0x1) ? -1 : 0, (mask & 0x2) ? -1 : 0, (mask & 0x4) ? -1 : 0, (mask & 0x8) ? -1 : 0};
        }

        inline __v4si greater_than_equal_mask(__v4df x, __v4df y)
        {
            __m256d cmp = _mm256_cmp_pd(x, y, _CMP_GE_OQ);
            int mask = _mm256_movemask_pd(cmp);
            return __v4si{(mask & 0x1) ? -1 : 0, (mask & 0x2) ? -1 : 0, (mask & 0x4) ? -1 : 0, (mask & 0x8) ? -1 : 0};
        }

        inline __v4si equal_mask(__v4df x, __v4df y)
        {
            __m256d cmp = _mm256_cmp_pd(x, y, _CMP_EQ_OQ);
            int mask = _mm256_movemask_pd(cmp);
            return __v4si{(mask & 0x1) ? -1 : 0, (mask & 0x2) ? -1 : 0, (mask & 0x4) ? -1 : 0, (mask & 0x8) ? -1 : 0};
        }

        inline __v4si not_equal_mask(__v4df x, __v4df y)
        {
            __m256d cmp = _mm256_cmp_pd(x, y, _CMP_NEQ_OQ);
            int mask = _mm256_movemask_pd(cmp);
            return __v4si{(mask & 0x1) ? -1 : 0, (mask & 0x2) ? -1 : 0, (mask & 0x4) ? -1 : 0, (mask & 0x8) ? -1 : 0};
        }
    #endif

        inline __v4df round(__v4df x) { return _mm256_round_pd(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC); }

        inline __v4df floor(__v4df x) { return _mm256_floor_pd(x); }

        inline __v4df trunc(__v4df x) { return _mm256_round_pd(x, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC); }

        inline __v4df ceil(__v4df x) { return _mm256_ceil_pd(x); }

        inline __v4df min(__v4df a, __v4df b) { return _mm256_min_pd(a, b); }

        inline __v4df max(__v4df a, __v4df b) { return _mm256_max_pd(a, b); }

        inline float extract_scalar(__v4sf v) { return _mm_cvtss_f32(v); }
#endif
    } // namespace details
} // namespace amal