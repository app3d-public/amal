#pragma once

#if defined(__SSE2__)
    #include <emmintrin.h>
#endif
#if defined(__AVX__)
    #include <immintrin.h>
#endif

namespace amal
{
    namespace internal
    {
#if defined(__SSE2__)
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

        inline __v4si not_equal_mask(__v4si x, __v4si y) { return (x ^ y); }

    #if defined(__SSE4_1__)
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
    #else
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

    #endif
#endif
#if defined(__AVX__)
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
#endif
    } // namespace internal
} // namespace amal