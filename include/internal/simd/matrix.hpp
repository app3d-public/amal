#pragma once

#include "../type_info.hpp"
#include "geometric.hpp"

namespace amal
{
    namespace internal
    {
#ifdef __SSE2__
        template <int C1, int R1, int C2, int R2, typename T>
        inline void multiply_matrix(T const (&m1)[C1], T const (&m2)[C2], T (&out)[C2]);

    #include <matrix_multiply_v4sf.hpp>

        template <length_t C, length_t R>
        inline __v4sf multiply_matrix(__v4sf const (&m)[C], __v4sf const &v)
        {
            __v4sf result;
    #if defined(AMAL_FMA_ENABLE)
            result = _mm_mul_ps(m[0], _mm_shuffle_ps(v, v, _MM_SHUFFLE(0, 0, 0, 0)));
            result = _mm_fmadd_ps(m[1], _mm_shuffle_ps(v, v, _MM_SHUFFLE(1, 1, 1, 1)), result);
            if constexpr (C > 2) result = _mm_fmadd_ps(m[2], _mm_shuffle_ps(v, v, _MM_SHUFFLE(2, 2, 2, 2)), result);
            if constexpr (C > 3) result = _mm_fmadd_ps(m[3], _mm_shuffle_ps(v, v, _MM_SHUFFLE(3, 3, 3, 3)), result);
    #else
            result = _mm_setzero_ps();
            result = _mm_add_ps(result, _mm_mul_ps(m[0], _mm_shuffle_ps(v, v, _MM_SHUFFLE(0, 0, 0, 0))));
            result = _mm_add_ps(result, _mm_mul_ps(m[1], _mm_shuffle_ps(v, v, _MM_SHUFFLE(1, 1, 1, 1))));
            if constexpr (C > 2)
                result = _mm_add_ps(result, _mm_mul_ps(m[2], _mm_shuffle_ps(v, v, _MM_SHUFFLE(2, 2, 2, 2))));
            if constexpr (C > 3)
                result = _mm_add_ps(result, _mm_mul_ps(m[3], _mm_shuffle_ps(v, v, _MM_SHUFFLE(3, 3, 3, 3))));
    #endif
            return result;
        }

        template <length_t N, length_t M>
        void transpose(__v4sf const (&in)[N], __v4sf (&out)[M]);

        template <>
        inline void transpose(__v4sf const (&in)[2], __v4sf (&out)[2])
        {
            __v4sf t0 = _mm_unpacklo_ps(in[0], in[1]);
            out[0] = _mm_movelh_ps(t0, t0);
            out[1] = _mm_movehl_ps(t0, t0);
        }

        template <>
        inline void transpose(__v4sf const (&in)[3], __v4sf (&out)[3])
        {
            __v4sf t0 = _mm_unpacklo_ps(in[0], in[1]);
            __v4sf t1 = _mm_unpackhi_ps(in[0], in[1]);
            out[0] = _mm_movelh_ps(t0, t1);
            out[1] = _mm_movehl_ps(t1, t0);
            out[2] = _mm_shuffle_ps(in[0], in[1], _MM_SHUFFLE(2, 2, 2, 2));
            out[2] = _mm_insert_ps(out[2], in[2], 0b00100000);
        }

        template <>
        inline void transpose(__v4sf const (&in)[2], __v4sf (&out)[3])
        {
            __v4sf t0 = _mm_unpacklo_ps(in[0], in[1]);
            __v4sf t1 = _mm_unpackhi_ps(in[0], in[1]);
            out[0] = _mm_movelh_ps(t0, t0);
            out[1] = _mm_movehl_ps(t0, t0);
            out[2] = t1;
        }

        template <>
        inline void transpose(__v4sf const (&in)[3], __v4sf (&out)[2])
        {
            __v4sf t0 = _mm_unpacklo_ps(in[0], in[1]);
            __v4sf t1 = _mm_unpacklo_ps(in[2], _mm_setzero_ps());
            out[0] = _mm_movelh_ps(t0, t1);
            out[1] = _mm_movehl_ps(t1, t0);
        }

        template <>
        inline void transpose(__v4sf const (&in)[4], __v4sf (&out)[4])
        {
            __v4sf t0 = _mm_unpacklo_ps(in[0], in[1]);
            __v4sf t1 = _mm_unpackhi_ps(in[0], in[1]);
            __v4sf t2 = _mm_unpacklo_ps(in[2], in[3]);
            __v4sf t3 = _mm_unpackhi_ps(in[2], in[3]);
            out[0] = _mm_movelh_ps(t0, t2);
            out[1] = _mm_movehl_ps(t2, t0);
            out[2] = _mm_movelh_ps(t1, t3);
            out[3] = _mm_movehl_ps(t3, t1);
        }

        template <>
        inline void transpose(__v4sf const (&in)[2], __v4sf (&out)[4])
        {
            __v4sf t0 = _mm_unpacklo_ps(in[0], in[1]);
            __v4sf t1 = _mm_unpackhi_ps(in[0], in[1]);
            out[0] = _mm_movelh_ps(t0, t0);
            out[1] = _mm_movehl_ps(t0, t0);
            out[2] = _mm_movelh_ps(t1, t1);
            out[3] = _mm_movehl_ps(t1, t1);
        }

        template <>
        inline void transpose(__v4sf const (&in)[3], __v4sf (&out)[4])
        {
            __v4sf t0 = _mm_unpacklo_ps(in[0], in[1]);
            __v4sf t1 = _mm_unpackhi_ps(in[0], in[1]);
            __v4sf t2 = _mm_unpacklo_ps(in[2], _mm_setzero_ps());
            __v4sf t3 = _mm_unpackhi_ps(in[2], _mm_setzero_ps());

            out[0] = _mm_movelh_ps(t0, t2);
            out[1] = _mm_movehl_ps(t2, t0);
            out[2] = _mm_movelh_ps(t1, t3);
            out[3] = _mm_movehl_ps(t3, t1);
        }

        template <>
        inline void transpose(__v4sf const (&in)[4], __v4sf (&out)[2])
        {
            __v4sf t0 = _mm_unpacklo_ps(in[0], in[1]);
            __v4sf t1 = _mm_unpacklo_ps(in[2], in[3]);
            out[0] = _mm_movelh_ps(t0, t1);
            out[1] = _mm_movehl_ps(t1, t0);
        }

        template <>
        inline void transpose(__v4sf const (&in)[4], __v4sf (&out)[3])
        {
            __v4sf t0 = _mm_unpacklo_ps(in[0], in[1]);
            __v4sf t1 = _mm_unpackhi_ps(in[0], in[1]);
            __v4sf t2 = _mm_unpacklo_ps(in[2], in[3]);
            __v4sf t3 = _mm_unpackhi_ps(in[2], in[3]);
            out[0] = _mm_movelh_ps(t0, t2);
            out[1] = _mm_movehl_ps(t2, t0);
            out[2] = _mm_movelh_ps(t1, t3);
        }

        inline __v4sf determinant(__v4sf const (&m)[2])
        {
            __v4sf a = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 0, 0));
            __v4sf d = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(1, 1, 1, 1));
            __v4sf bc = _mm_mul_ps(_mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(1, 1, 1, 1)),
                                   _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 0)));
    #if defined(AMAL_FMA_ENABLE)
            return _mm_fmsub_ps(a, d, bc);
    #else
            __v4sf ad = _mm_mul_ps(a, d);
            return _mm_sub_ps(ad, bc);
    #endif
        }

        inline __v4sf determinant(__v4sf const (&m)[3])
        {
            // m[0] = [a, b, c, _]
            // m[1] = [d, e, f, _]
            // m[2] = [g, h, i, _]

            // ei - fh
            __v4sf ei_fh = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 2, 2)),  // [f, f, f, f]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 1, 1))), // [h, h, h, h]
                                      _mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 1, 1)),  // [e, e, e, e]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 2, 2)))  // [i, i, i, i]
            ); // [ei - fh, ei - fh, ei - fh, ei - fh]

            __v4sf a = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 0, 0)); // [a, a, a, a]
            __v4sf term0 = _mm_mul_ps(a, ei_fh);

            // di - fg
            __v4sf di_fg = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 0)),  // [d, d, d, d]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 2, 2))), // [i, i, i, i]
                                      _mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 2, 2)),  // [f, f, f, f]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 0, 0)))  // [g, g, g, g]
            );
            __v4sf b = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 1, 1));
            __v4sf term1 = _mm_mul_ps(b, di_fg);

            // dh - eg
            __v4sf dh_eg = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 0)),  // [d, d, d, d]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 1, 1))), // [h, h, h, h]
                                      _mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 1, 1)),  // [e, e, e, e]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 0, 0)))  // [g, g, g, g]
            );
            __v4sf c = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 2, 2));
            __v4sf term2 = _mm_mul_ps(c, dh_eg);

            return _mm_sub_ps(_mm_add_ps(term0, term2), term1);
        }

        inline __v4sf determinant(__v4sf const (&m)[4])
        {
            __v4sf s0 = _mm_mul_ps(_mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 1, 1, 2)),
                                                         _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(3, 2, 3, 3))),
                                              _mm_mul_ps(_mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(3, 2, 3, 3)),
                                                         _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(0, 1, 1, 2)))),
                                   _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 1)));

            __v4sf s1 =
                _mm_mul_ps(_mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 1, 2)),
                                                 _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(1, 2, 0, 0))),
                                      _mm_movehl_ps(_mm_setzero_ps(),
                                                    _mm_mul_ps(_mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 1, 2)),
                                                               _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(1, 2, 0, 0))))),
                           _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(2, 3, 3, 3)));

            __v4sf det = _mm_add_ps(s0, s1);
            det = _mm_mul_ps(det, _mm_setr_ps(1.f, -1.f, 1.f, -1.f));
            return _mm_mul_ps(det, m[0]);
        }

        inline void inverse_matrix(__v4sf const (&m)[2], __v4sf (&out)[2])
        {
            __v4sf rdet = _mm_rcp_ps(determinant(m));
            out[0] = _mm_mul_ps(_mm_setr_ps(m[1][1], -m[0][1], 0.f, 0.f), rdet);
            out[1] = _mm_mul_ps(_mm_setr_ps(-m[1][0], m[0][0], 0.f, 0.f), rdet);
        }

        inline void inverse_matrix(__v4sf const (&m)[3], __v4sf (&out)[3])
        {
            const __v4sf a = m[0];
            const __v4sf b = m[1];
            const __v4sf c = m[2];

            const __v4sf i0 = cross_yzx(b, c);
            const __v4sf i1 = cross_yzx(c, a);
            const __v4sf i2 = cross_yzx(a, b);

            out[0] = _mm_setr_ps(i0[0], i1[0], i2[0], 0.f);
            out[1] = _mm_setr_ps(i0[1], i1[1], i2[1], 0.f);
            out[2] = _mm_setr_ps(i0[2], i1[2], i2[2], 0.f);

            const __v4sf rdet = _mm_rcp_ps(dot(a, cross_yzx(b, c)));

            out[0] = _mm_mul_ps(out[0], rdet);
            out[1] = _mm_mul_ps(out[1], rdet);
            out[2] = _mm_mul_ps(out[2], rdet);
        }

        inline void inverse_matrix(__v4sf const (&m)[4], __v4sf (&out)[4])
        {
            const __v4sf A = m[0]; // [a0 a1 a2 a3]
            const __v4sf B = m[1]; // [b0 b1 b2 b3]
            const __v4sf C = m[2]; // [c0 c1 c2 c3]
            const __v4sf D = m[3]; // [d0 d1 d2 d3]

            const __v4sf AB_lo = _mm_movelh_ps(A, B); // [a0 a1 b0 b1]
            const __v4sf AB_hi = _mm_movehl_ps(B, A); // [a2 a3 b2 b3]
            const __v4sf CD_lo = _mm_movelh_ps(C, D); // [c0 c1 d0 d1]
            const __v4sf CD_hi = _mm_movehl_ps(D, C); // [c2 c3 d2 d3]

            const __v4sf AB_lo_swapped = _mm_shuffle_ps(AB_lo, AB_lo, _MM_SHUFFLE(2, 3, 0, 1)); // [b0 b1 a0 a1]
            const __v4sf CD_hi_swapped = _mm_shuffle_ps(CD_hi, CD_hi, _MM_SHUFFLE(2, 3, 0, 1)); // [d2 d3 c2 c3]

            const __v4sf ab_cd = _mm_mul_ps(AB_lo, CD_hi);
            const __v4sf ab_cd_swapped = _mm_mul_ps(AB_lo_swapped, CD_hi_swapped);
            const __v4sf cof0 = _mm_sub_ps(ab_cd, ab_cd_swapped);

            const __v4sf AB_hi_swapped = _mm_shuffle_ps(AB_hi, AB_hi, _MM_SHUFFLE(2, 3, 0, 1));
            const __v4sf CD_lo_swapped = _mm_shuffle_ps(CD_lo, CD_lo, _MM_SHUFFLE(2, 3, 0, 1));

            const __v4sf ab2_cd2 = _mm_mul_ps(AB_hi, CD_lo);
            const __v4sf ab2_cd2_swapped = _mm_mul_ps(AB_hi_swapped, CD_lo_swapped);
            const __v4sf cof1 = _mm_sub_ps(ab2_cd2, ab2_cd2_swapped);

            out[0] = _mm_shuffle_ps(cof0, cof1, _MM_SHUFFLE(2, 0, 2, 0)); // [+, -, +, -]
            out[1] = _mm_shuffle_ps(cof0, cof1, _MM_SHUFFLE(3, 1, 3, 1)); // [-, +, -, +]

            const __v4sf sign0 = _mm_setr_ps(+1.f, -1.f, +1.f, -1.f);
            const __v4sf sign1 = _mm_setr_ps(-1.f, +1.f, -1.f, +1.f);

            out[0] = _mm_mul_ps(out[0], sign0);
            out[1] = _mm_mul_ps(out[1], sign1);

            const __v4sf CD_lo_swapped2 = _mm_shuffle_ps(CD_lo, CD_lo, _MM_SHUFFLE(2, 3, 0, 1));
            const __v4sf AB_hi_swapped2 = _mm_shuffle_ps(AB_hi, AB_hi, _MM_SHUFFLE(2, 3, 0, 1));

            const __v4sf cd_ab = _mm_mul_ps(CD_lo, AB_hi);
            const __v4sf cd_ab_swapped = _mm_mul_ps(CD_lo_swapped2, AB_hi_swapped2);
            const __v4sf cof2 = _mm_sub_ps(cd_ab, cd_ab_swapped);

            const __v4sf CD_hi_swapped2 = _mm_shuffle_ps(CD_hi, CD_hi, _MM_SHUFFLE(2, 3, 0, 1));
            const __v4sf AB_lo_swapped2 = _mm_shuffle_ps(AB_lo, AB_lo, _MM_SHUFFLE(2, 3, 0, 1));

            const __v4sf cd_ab2 = _mm_mul_ps(CD_hi, AB_lo);
            const __v4sf cd_ab2_swapped = _mm_mul_ps(CD_hi_swapped2, AB_lo_swapped2);
            const __v4sf cof3 = _mm_sub_ps(cd_ab2, cd_ab2_swapped);

            out[2] = _mm_shuffle_ps(cof2, cof3, _MM_SHUFFLE(2, 0, 2, 0)); // [+, -, +, -]
            out[3] = _mm_shuffle_ps(cof2, cof3, _MM_SHUFFLE(3, 1, 3, 1)); // [-, +, -, +]

            const __v4sf sign2 = sign0;
            const __v4sf sign3 = sign1;
            out[2] = _mm_mul_ps(out[2], sign2);
            out[3] = _mm_mul_ps(out[3], sign3);

            __v4sf detvec = dot(m[0], out[0]);
            __v4sf rdet = _mm_rcp_ps(detvec);

            out[0] = _mm_mul_ps(out[0], rdet);
            out[1] = _mm_mul_ps(out[1], rdet);
            out[2] = _mm_mul_ps(out[2], rdet);
            out[3] = _mm_mul_ps(out[3], rdet);
        }

        template <length_t N, length_t M>
        void transpose(__v4si const (&in)[N], __v4si (&out)[M]);

        inline void transpose(__v4si const (&in)[2], __v4si (&out)[2])
        {
            __m128i t0 = _mm_unpacklo_epi32(in[0], in[1]);
            out[0] = _mm_unpacklo_epi64(t0, t0); // a0 b0
            out[1] = _mm_unpackhi_epi64(t0, t0); // a1 b1
        }

        inline void transpose(__v4si const (&in)[3], __v4si (&out)[2])
        {
            __m128i t0 = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i t1 = _mm_unpacklo_epi32(in[2], _mm_setzero_si128());
            out[0] = _mm_unpacklo_epi64(t0, t1);
            out[1] = _mm_unpackhi_epi64(t0, t1);
        }

        inline void transpose(__v4si const (&in)[2], __v4si (&out)[3])
        {
            __m128i t0 = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i t1 = _mm_unpackhi_epi32(in[0], in[1]);
            out[0] = _mm_unpacklo_epi64(t0, t0);
            out[1] = _mm_unpackhi_epi64(t0, t0);
            out[2] = _mm_unpacklo_epi64(t1, t1);
        }

        inline void transpose(__v4si const (&in)[3], __v4si (&out)[3])
        {
            __m128i t0 = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i t1 = _mm_unpackhi_epi32(in[0], in[1]);
            out[0] = _mm_unpacklo_epi64(t0, t0);
            out[1] = _mm_unpackhi_epi64(t0, t0);
            out[2] = _mm_unpacklo_epi64(t1, in[2]);
        }

        inline void transpose(__v4si const (&in)[4], __v4si (&out)[2])
        {
            __m128i ab_lo = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i cd_lo = _mm_unpacklo_epi32(in[2], in[3]);

            out[0] = _mm_unpacklo_epi64(ab_lo, cd_lo);
            out[1] = _mm_unpackhi_epi64(ab_lo, cd_lo);
        }

        inline void transpose(__v4si const (&in)[2], __v4si (&out)[4])
        {
            __m128i lo = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i hi = _mm_unpackhi_epi32(in[0], in[1]);

            out[0] = _mm_unpacklo_epi64(lo, lo);
            out[1] = _mm_unpackhi_epi64(lo, lo);
            out[2] = _mm_unpacklo_epi64(hi, hi);
            out[3] = _mm_unpackhi_epi64(hi, hi);
        }

        inline void transpose(__v4si const (&in)[3], __v4si (&out)[4])
        {
            __m128i ab_lo = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i ab_hi = _mm_unpackhi_epi32(in[0], in[1]);
            __m128i c_lo = _mm_unpacklo_epi32(in[2], _mm_setzero_si128());
            __m128i c_hi = _mm_unpackhi_epi32(in[2], _mm_setzero_si128());

            out[0] = _mm_unpacklo_epi64(ab_lo, c_lo);
            out[1] = _mm_unpackhi_epi64(ab_lo, c_lo);
            out[2] = _mm_unpacklo_epi64(ab_hi, c_hi);
            out[3] = _mm_unpackhi_epi64(ab_hi, c_hi);
        }

        inline void transpose(__v4si const (&in)[4], __v4si (&out)[3])
        {
            __m128i ab_lo = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i ab_hi = _mm_unpackhi_epi32(in[0], in[1]);
            __m128i cd_lo = _mm_unpacklo_epi32(in[2], in[3]);
            __m128i cd_hi = _mm_unpackhi_epi32(in[2], in[3]);

            out[0] = _mm_unpacklo_epi64(ab_lo, cd_lo);
            out[1] = _mm_unpackhi_epi64(ab_lo, cd_lo);
            out[2] = _mm_unpacklo_epi64(ab_hi, cd_hi);
        }

        inline void transpose(__v4si const (&in)[4], __v4si (&out)[4])
        {
            __m128i ab_lo = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i ab_hi = _mm_unpackhi_epi32(in[0], in[1]);
            __m128i cd_lo = _mm_unpacklo_epi32(in[2], in[3]);
            __m128i cd_hi = _mm_unpackhi_epi32(in[2], in[3]);

            out[0] = _mm_unpacklo_epi64(ab_lo, cd_lo);
            out[1] = _mm_unpackhi_epi64(ab_lo, cd_lo);
            out[2] = _mm_unpacklo_epi64(ab_hi, cd_hi);
            out[3] = _mm_unpackhi_epi64(ab_hi, cd_hi);
        }

        template <>
        inline void transpose(__v4si const (&in)[4], __v4si (&out)[4])
        {
            __m128i a = in[0];
            __m128i b = in[1];
            __m128i c = in[2];
            __m128i d = in[3];

            // unpack 32-bit integers
            __m128i ab_lo = _mm_unpacklo_epi32(a, b); // a0 b0 a1 b1
            __m128i ab_hi = _mm_unpackhi_epi32(a, b); // a2 b2 a3 b3
            __m128i cd_lo = _mm_unpacklo_epi32(c, d); // c0 d0 c1 d1
            __m128i cd_hi = _mm_unpackhi_epi32(c, d); // c2 d2 c3 d3

            out[0] = _mm_unpacklo_epi64(ab_lo, cd_lo); // a0 b0 c0 d0
            out[1] = _mm_unpackhi_epi64(ab_lo, cd_lo); // a1 b1 c1 d1
            out[2] = _mm_unpacklo_epi64(ab_hi, cd_hi); // a2 b2 c2 d2
            out[3] = _mm_unpackhi_epi64(ab_hi, cd_hi); // a3 b3 c3 d3
        }

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

    #include <matrix_multiply_v4si.hpp>

        inline __v4si determinant(__v4si const (&m)[2])
        {
            __v4si ab = mm_mullo_epi32_compat(m[0], _mm_shuffle_epi32(m[1], _MM_SHUFFLE(0, 0, 1, 1)));
            __v4si det = _mm_sub_epi32(ab, _mm_shuffle_epi32(ab, _MM_SHUFFLE(1, 1, 1, 1)));
            return det; // .x = det
        }

        inline __v4si determinant(__v4si const (&m)[3])
        {
            __v4si a = m[0];
            __v4si b = m[1];
            __v4si c = m[2];

            __v4si bc = mm_mullo_epi32_compat(_mm_shuffle_epi32(b, _MM_SHUFFLE(1, 2, 0, 1)),
                                              _mm_shuffle_epi32(c, _MM_SHUFFLE(2, 0, 1, 2)));
            __v4si cb = mm_mullo_epi32_compat(_mm_shuffle_epi32(b, _MM_SHUFFLE(2, 0, 1, 2)),
                                              _mm_shuffle_epi32(c, _MM_SHUFFLE(1, 2, 0, 1)));

            __v4si cross = _mm_sub_epi32(bc, cb);
            __v4si det = mm_mullo_epi32_compat(a, cross);

            return _mm_add_epi32(_mm_shuffle_epi32(det, _MM_SHUFFLE(0, 0, 0, 0)),
                                 _mm_add_epi32(_mm_shuffle_epi32(det, _MM_SHUFFLE(1, 1, 1, 1)),
                                               _mm_shuffle_epi32(det, _MM_SHUFFLE(2, 2, 2, 2))));
        }

        inline __v4si determinant(__v4si const (&m)[4])
        {
            // s0 = (m2[0,1,1,2] * m3[3,2,3,3] - m2[3,2,3,3] * m3[0,1,1,2]) * m1[0,0,0,1]
            __v4si s0 = mm_mullo_epi32_compat(
                _mm_sub_epi32(mm_mullo_epi32_compat(_mm_shuffle_epi32(m[2], _MM_SHUFFLE(0, 1, 1, 2)),
                                                    _mm_shuffle_epi32(m[3], _MM_SHUFFLE(3, 2, 3, 3))),
                              mm_mullo_epi32_compat(_mm_shuffle_epi32(m[2], _MM_SHUFFLE(3, 2, 3, 3)),
                                                    _mm_shuffle_epi32(m[3], _MM_SHUFFLE(0, 1, 1, 2)))),
                _mm_shuffle_epi32(m[1], _MM_SHUFFLE(0, 0, 0, 1)));

            // prod = m2[0,0,1,2] * m3[1,2,0,0]
            __v4si prod = mm_mullo_epi32_compat(_mm_shuffle_epi32(m[2], _MM_SHUFFLE(0, 0, 1, 2)),
                                                _mm_shuffle_epi32(m[3], _MM_SHUFFLE(1, 2, 0, 0)));
            __v4si upper = _mm_castps_si128(_mm_movehl_ps(_mm_setzero_ps(), _mm_castsi128_ps(prod)));

            __v4si s1 =
                mm_mullo_epi32_compat(_mm_sub_epi32(prod, upper), _mm_shuffle_epi32(m[1], _MM_SHUFFLE(2, 3, 3, 3)));

            __v4si det = _mm_add_epi32(s0, s1);
            __v4si signs = _mm_setr_epi32(1, -1, 1, -1);
            det = mm_mullo_epi32_compat(det, signs);
            return mm_mullo_epi32_compat(det, m[0]);
        }

        template <length_t C, length_t R>
        inline __v4si multiply_matrix(__v4si const (&m)[C], __v4si const &v)
        {
            __v4si result = _mm_setzero_si128();
            result = _mm_add_epi32(result, mm_mullo_epi32_compat(m[0], _mm_shuffle_epi32(v, _MM_SHUFFLE(0, 0, 0, 0))));
            result = _mm_add_epi32(result, mm_mullo_epi32_compat(m[1], _mm_shuffle_epi32(v, _MM_SHUFFLE(1, 1, 1, 1))));
            if constexpr (C > 2)
                result =
                    _mm_add_epi32(result, mm_mullo_epi32_compat(m[2], _mm_shuffle_epi32(v, _MM_SHUFFLE(2, 2, 2, 2))));
            if constexpr (C > 3)
                result =
                    _mm_add_epi32(result, mm_mullo_epi32_compat(m[3], _mm_shuffle_epi32(v, _MM_SHUFFLE(3, 3, 3, 3))));
            return result;
        }

    #ifdef __AVX__
        #include <matrix_multiply_v4df.hpp>

        template <size_t C, size_t R>
        inline __v4df multiply_matrix(__v4df const (&m)[C], __v4df const &v)
        {
            __m256d result;
        #if defined(AMAL_FMA_ENABLE)
            result = _mm256_mul_pd(reinterpret_cast<const __m256d &>(m[0]),
                                   _mm256_permute4x64_pd(reinterpret_cast<const __m256d &>(v), 0x00));
            result = _mm256_fmadd_pd(reinterpret_cast<const __m256d &>(m[1]),
                                     _mm256_permute4x64_pd(reinterpret_cast<const __m256d &>(v), 0x55), result);
            if constexpr (C > 2)
                result = _mm256_fmadd_pd(reinterpret_cast<const __m256d &>(m[2]),
                                         _mm256_permute4x64_pd(reinterpret_cast<const __m256d &>(v), 0xAA), result);
            if constexpr (C > 3)
                result = _mm256_fmadd_pd(reinterpret_cast<const __m256d &>(m[3]),
                                         _mm256_permute4x64_pd(reinterpret_cast<const __m256d &>(v), 0xFF), result);
        #else
            result = _mm256_setzero_pd();
            result =
                _mm256_add_pd(result, _mm256_mul_pd(reinterpret_cast<const __m256d &>(m[0]),
                                                    _mm256_permute4x64_pd(reinterpret_cast<const __m256d &>(v), 0x00)));
            result =
                _mm256_add_pd(result, _mm256_mul_pd(reinterpret_cast<const __m256d &>(m[1]),
                                                    _mm256_permute4x64_pd(reinterpret_cast<const __m256d &>(v), 0x55)));
            if constexpr (C > 2)
                result = _mm256_add_pd(
                    result, _mm256_mul_pd(reinterpret_cast<const __m256d &>(m[2]),
                                          _mm256_permute4x64_pd(reinterpret_cast<const __m256d &>(v), 0xAA)));
            if constexpr (C > 3)
                result = _mm256_add_pd(
                    result, _mm256_mul_pd(reinterpret_cast<const __m256d &>(m[3]),
                                          _mm256_permute4x64_pd(reinterpret_cast<const __m256d &>(v), 0xFF)));
        #endif

            return reinterpret_cast<__v4df &>(result);
        }

        template <length_t N, length_t M>
        void transpose(__v4df const (&in)[N], __v4df (&out)[M]);

        template <>
        inline void transpose(__v4df const (&in)[2], __v4df (&out)[2])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            out[0] = _mm256_permute2f128_pd(t0, t0, 0x20);
            out[1] = _mm256_permute2f128_pd(t0, t0, 0x31);
        }

        template <>
        inline void transpose(__v4df const (&in)[3], __v4df (&out)[3])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            __m256d t1 = _mm256_unpackhi_pd(in[0], in[1]);
            out[0] = _mm256_permute2f128_pd(t0, t1, 0x20);
            out[1] = _mm256_permute2f128_pd(t0, t1, 0x31);
            out[2] =
                _mm256_blend_pd(_mm256_castpd128_pd256(_mm256_castpd256_pd128(in[2])),
                                _mm256_insertf128_pd(_mm256_setzero_pd(), _mm256_extractf128_pd(in[2], 1), 1), 0b1100);
        }

        template <>
        inline void transpose(__v4df const (&in)[2], __v4df (&out)[3])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            __m256d t1 = _mm256_unpackhi_pd(in[0], in[1]);
            out[0] = _mm256_permute2f128_pd(t0, t0, 0x20);
            out[1] = _mm256_permute2f128_pd(t0, t0, 0x31);
            out[2] = t1;
        }

        template <>
        inline void transpose(__v4df const (&in)[3], __v4df (&out)[2])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            __m256d t1 = _mm256_unpacklo_pd(in[2], _mm256_setzero_pd());
            out[0] = _mm256_permute2f128_pd(t0, t1, 0x20);
            out[1] = _mm256_permute2f128_pd(t0, t1, 0x31);
        }

        template <>
        inline void transpose(__v4df const (&in)[4], __v4df (&out)[4])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            __m256d t1 = _mm256_unpackhi_pd(in[0], in[1]);
            __m256d t2 = _mm256_unpacklo_pd(in[2], in[3]);
            __m256d t3 = _mm256_unpackhi_pd(in[2], in[3]);
            out[0] = _mm256_permute2f128_pd(t0, t2, 0x20);
            out[1] = _mm256_permute2f128_pd(t0, t2, 0x31);
            out[2] = _mm256_permute2f128_pd(t1, t3, 0x20);
            out[3] = _mm256_permute2f128_pd(t1, t3, 0x31);
        }

        template <>
        inline void transpose(__v4df const (&in)[4], __v4df (&out)[3])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            __m256d t1 = _mm256_unpackhi_pd(in[0], in[1]);
            __m256d t2 = _mm256_unpacklo_pd(in[2], in[3]);
            __m256d t3 = _mm256_unpackhi_pd(in[2], in[3]);

            out[0] = _mm256_permute2f128_pd(t0, t2, 0x20);
            out[1] = _mm256_permute2f128_pd(t0, t2, 0x31);
            out[2] = _mm256_permute2f128_pd(t1, t3, 0x20);
        }

        template <>
        inline void transpose(__v4df const (&in)[3], __v4df (&out)[4])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            __m256d t1 = _mm256_unpackhi_pd(in[0], in[1]);
            __m256d t2 = _mm256_unpacklo_pd(in[2], _mm256_setzero_pd());
            __m256d t3 = _mm256_unpackhi_pd(in[2], _mm256_setzero_pd());

            out[0] = _mm256_permute2f128_pd(t0, t2, 0x20);
            out[1] = _mm256_permute2f128_pd(t0, t2, 0x31);
            out[2] = _mm256_permute2f128_pd(t1, t3, 0x20);
            out[3] = _mm256_permute2f128_pd(t1, t3, 0x31);
        }

        template <>
        inline void transpose(__v4df const (&in)[2], __v4df (&out)[4])
        {
            __m128d a = _mm256_castpd256_pd128(in[0]);
            __m128d b = _mm256_castpd256_pd128(in[1]);
            __m128d c = _mm256_extractf128_pd(in[0], 1);
            __m128d d = _mm256_extractf128_pd(in[1], 1);
            out[0] = _mm256_insertf128_pd(_mm256_setzero_pd(), _mm_unpacklo_pd(a, b), 0);
            out[1] = _mm256_insertf128_pd(_mm256_setzero_pd(), _mm_unpackhi_pd(a, b), 0);
            out[2] = _mm256_insertf128_pd(_mm256_setzero_pd(), _mm_unpacklo_pd(c, d), 0);
            out[3] = _mm256_insertf128_pd(_mm256_setzero_pd(), _mm_unpackhi_pd(c, d), 0);
        }

        template <>
        inline void transpose(__v4df const (&in)[4], __v4df (&out)[2])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            __m256d t1 = _mm256_unpackhi_pd(in[0], in[1]);
            __m256d t2 = _mm256_unpacklo_pd(in[2], in[3]);
            __m256d t3 = _mm256_unpackhi_pd(in[2], in[3]);

            out[0] = _mm256_permute2f128_pd(t0, t2, 0x20);
            out[1] = _mm256_permute2f128_pd(t1, t3, 0x20);
        }

        inline __v4df determinant(__v4df const (&m)[2])
        {
            __m256d mul = _mm256_mul_pd(m[0], _mm256_permute4x64_pd(m[1], 0b01000000)); // m[1]: [a, a, b, b]
            return _mm256_sub_pd(mul, _mm256_permute4x64_pd(mul, 0b01010101));          // mul - [b*c, b*c, b*c, b*c]
        }

        inline __v4df determinant(__v4df const (&m)[3])
        {
            __m256d a = m[0];
            __m256d b = m[1];
            __m256d c = m[2];

            __m256d bc = _mm256_mul_pd(_mm256_permute4x64_pd(b, 0b11010010),  // [b1, b2, b0, b1]
                                       _mm256_permute4x64_pd(c, 0b10010011)); // [c2, c0, c1, c2]

            __m256d cb = _mm256_mul_pd(_mm256_permute4x64_pd(b, 0b10010011),  // [b2, b0, b1, b2]
                                       _mm256_permute4x64_pd(c, 0b11010010)); // [c1, c2, c0, c1]

            __m256d cross = _mm256_sub_pd(bc, cb);
            __m256d mul = _mm256_mul_pd(a, cross);

            __m256d tmp1 = _mm256_permute4x64_pd(mul, 0b00000000);
            __m256d tmp2 = _mm256_permute4x64_pd(mul, 0b01010101);
            __m256d tmp3 = _mm256_permute4x64_pd(mul, 0b10101010);

            return _mm256_add_pd(tmp1, _mm256_add_pd(tmp2, tmp3));
        }

        inline __v4df determinant(__v4df const (&m)[4])
        {
            // s0 = (m2[0,1,1,2] * m3[3,2,3,3] - m2[3,2,3,3] * m3[0,1,1,2]) * m1[0,0,0,1]
            __v4df m2a = _mm256_permute4x64_pd(m[2], _MM_SHUFFLE(0, 1, 1, 2));
            __v4df m3a = _mm256_permute4x64_pd(m[3], _MM_SHUFFLE(3, 2, 3, 3));
            __v4df m2b = _mm256_permute4x64_pd(m[2], _MM_SHUFFLE(3, 2, 3, 3));
            __v4df m3b = _mm256_permute4x64_pd(m[3], _MM_SHUFFLE(0, 1, 1, 2));

            __v4df s0 = _mm256_mul_pd(_mm256_sub_pd(_mm256_mul_pd(m2a, m3a), _mm256_mul_pd(m2b, m3b)),
                                      _mm256_permute4x64_pd(m[1], _MM_SHUFFLE(0, 0, 0, 1)));

            // s1 = ((m2[0,0,1,2] * m3[1,2,0,0]) - hi(m2*m3)) * m1[2,3,3,3]
            __v4df t1 = _mm256_mul_pd(_mm256_permute4x64_pd(m[2], _MM_SHUFFLE(0, 0, 1, 2)),
                                      _mm256_permute4x64_pd(m[3], _MM_SHUFFLE(1, 2, 0, 0)));

            __m128d hi128 = _mm256_extractf128_pd(t1, 1);
            __v4df upper = _mm256_castpd128_pd256(hi128);                  // [hi2, hi3, ?, ?]
            upper = _mm256_permute4x64_pd(upper, _MM_SHUFFLE(1, 1, 1, 1)); // [hi3, hi3, hi3, hi3]

            __v4df s1 = _mm256_mul_pd(_mm256_sub_pd(t1, upper), _mm256_permute4x64_pd(m[1], _MM_SHUFFLE(2, 3, 3, 3)));

            __v4df det = _mm256_add_pd(s0, s1);
            det = _mm256_mul_pd(det, _mm256_setr_pd(1.0, -1.0, 1.0, -1.0));
            return _mm256_mul_pd(det, m[0]);
        }

        inline void inverse_matrix(__v4df const (&m)[2], __v4df (&out)[2])
        {
            __v4df a = m[0], b = m[1];
            __v4df det = _mm256_sub_pd(_mm256_mul_pd(a, _mm256_permute4x64_pd(b, _MM_SHUFFLE(0, 0, 3, 3))),
                                       _mm256_mul_pd(_mm256_permute4x64_pd(a, _MM_SHUFFLE(0, 0, 3, 3)), b));
            __v4df invd = _mm256_div_pd(_mm256_set1_pd(1.0), det);
            __v4df adj0 =
                _mm256_setr_pd(+((double *)&b)[3], -((double *)&b)[2], -((double *)&a)[3], +((double *)&a)[2]);
            __v4df adj1 =
                _mm256_setr_pd(-((double *)&b)[1], +((double *)&b)[0], +((double *)&a)[1], -((double *)&a)[0]);
            out[0] = _mm256_mul_pd(adj0, invd);
            out[1] = _mm256_mul_pd(adj1, invd);
        }

        inline void inverse_matrix(__v4df const (&m)[3], __v4df (&out)[3])
        {
            __v4df a = m[0], b = m[1], c = m[2];
            __v4df i0 = cross_yzx(b, c);
            __v4df i1 = cross_yzx(c, a);
            __v4df i2 = cross_yzx(a, b);
            __v4df det = dot(a, i0);
            __v4df invd = _mm256_div_pd(_mm256_set1_pd(1.0), det);
            out[0] = _mm256_mul_pd(i0, invd);
            out[1] = _mm256_mul_pd(i1, invd);
            out[2] = _mm256_mul_pd(i2, invd);
        }

        inline void inverse_matrix(__v4df const (&m)[4], __v4df (&out)[4])
        {
            __v4df A = m[0], B = m[1], C = m[2], D = m[3];
        #if defined(__AVX2__)
            __v4df AB_lo = _mm256_unpacklo_pd(A, B);
            __v4df AB_hi = _mm256_unpackhi_pd(A, B);
            __v4df CD_lo = _mm256_unpacklo_pd(C, D);
            __v4df CD_hi = _mm256_unpackhi_pd(C, D);

            auto perm = [](const __v4df &v) { return _mm256_permute4x64_pd(v, _MM_SHUFFLE(2, 3, 0, 1)); };

            __v4df cof0 = _mm256_sub_pd(_mm256_mul_pd(AB_lo, CD_hi), _mm256_mul_pd(perm(AB_lo), perm(CD_hi)));
            __v4df cof1 = _mm256_sub_pd(_mm256_mul_pd(AB_hi, CD_lo), _mm256_mul_pd(perm(AB_hi), perm(CD_lo)));
            __v4df cof2 = _mm256_sub_pd(_mm256_mul_pd(CD_lo, AB_hi), _mm256_mul_pd(perm(CD_lo), perm(AB_hi)));
            __v4df cof3 = _mm256_sub_pd(_mm256_mul_pd(CD_hi, AB_lo), _mm256_mul_pd(perm(CD_hi), perm(AB_lo)));

            __v4df t0 = _mm256_unpacklo_pd(cof0, cof1);
            __v4df t1 = _mm256_unpacklo_pd(cof2, cof3);
            out[0] = _mm256_permute2f128_pd(t0, t1, 0x20);

            __v4df u0 = _mm256_unpackhi_pd(cof0, cof1);
            __v4df u1 = _mm256_unpackhi_pd(cof2, cof3);
            out[1] = _mm256_permute2f128_pd(u0, u1, 0x20);
            out[2] = _mm256_permute2f128_pd(t0, t1, 0x31);
            out[3] = _mm256_permute2f128_pd(u0, u1, 0x31);
        #else
            __v4df t0 = _mm256_shuffle_pd(A, B, 0x0);
            __v4df t1 = _mm256_shuffle_pd(C, D, 0x0);
            __v4df t2 = _mm256_shuffle_pd(A, B, 0xF);
            __v4df t3 = _mm256_shuffle_pd(C, D, 0xF);

            __v4df cof0 = _mm256_sub_pd(_mm256_mul_pd(t0, t3), _mm256_mul_pd(t2, t1));
            __v4df cof1 = _mm256_sub_pd(_mm256_mul_pd(t2, t0), _mm256_mul_pd(t1, t3));
            __v4df cof2 = _mm256_sub_pd(_mm256_mul_pd(t3, t1), _mm256_mul_pd(t0, t2));
            __v4df cof3 = _mm256_sub_pd(_mm256_mul_pd(t1, t3), _mm256_mul_pd(t2, t0));

            out[0] = _mm256_setr_pd(cof0[0], cof1[0], cof2[0], cof3[0]);
            out[1] = _mm256_setr_pd(cof0[1], cof1[1], cof2[1], cof3[1]);
            out[2] = _mm256_setr_pd(cof0[2], cof1[2], cof2[2], cof3[2]);
            out[3] = _mm256_setr_pd(cof0[3], cof1[3], cof2[3], cof3[3]);
        #endif
            const __v4df S0 = _mm256_setr_pd(+1, -1, +1, -1);
            const __v4df S1 = _mm256_setr_pd(-1, +1, -1, +1);
            out[0] = _mm256_mul_pd(out[0], S0);
            out[1] = _mm256_mul_pd(out[1], S1);
            out[2] = _mm256_mul_pd(out[2], S0);
            out[3] = _mm256_mul_pd(out[3], S1);

            __v4df det = dot(A, out[0]);
            __v4df invd = _mm256_div_pd(_mm256_set1_pd(1.0), det);
            out[0] = _mm256_mul_pd(out[0], invd);
            out[1] = _mm256_mul_pd(out[1], invd);
            out[2] = _mm256_mul_pd(out[2], invd);
            out[3] = _mm256_mul_pd(out[3], invd);
        }

    #endif
#endif
    } // namespace internal
} // namespace amal