#pragma once

#include <cmath>
#include "../type_info.hpp"
#include "geometric.hpp"

namespace amal
{
    namespace internal
    {
#ifdef __SSE2__
        template <int C1, int R1, int C2, int R2, typename T>
        inline void multiply_matrix(T const (&m1)[C1], T const (&m2)[C2], T (&out)[C2]);

        template <length_t C>
        inline __m128_u multiply_matrix(__m128_u const (&m)[C], __m128_u const &v)
        {
            __m128_u result;
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
        void transpose(__m128_u const (&in)[N], __m128_u (&out)[M]);

        template <>
        inline void transpose(__m128_u const (&in)[2], __m128_u (&out)[2])
        {
            __m128_u t0 = _mm_unpacklo_ps(in[0], in[1]);
            out[0] = _mm_movelh_ps(t0, t0);
            out[1] = _mm_movehl_ps(t0, t0);
        }

        template <>
        inline void transpose(__m128_u const (&in)[3], __m128_u (&out)[3])
        {
            __m128_u t0 = _mm_unpacklo_ps(in[0], in[1]);
            __m128_u t1 = _mm_unpackhi_ps(in[0], in[1]);
            out[0] = _mm_movelh_ps(t0, t1);
            out[1] = _mm_movehl_ps(t1, t0);
    #if defined(__SSE4_1__)
            out[2] = _mm_shuffle_ps(in[0], in[1], _MM_SHUFFLE(2, 2, 2, 2));
            out[2] = _mm_insert_ps(out[2], in[2], 0b00100000);
    #else
            __m128_u z = _mm_xor_ps(in[2], in[2]);
            __m128_u t3 = _mm_unpackhi_ps(in[2], z);
            out[2] = _mm_movelh_ps(t1, t3);
    #endif
        }

        template <>
        inline void transpose(__m128_u const (&in)[2], __m128_u (&out)[3])
        {
            __m128_u t0 = _mm_unpacklo_ps(in[0], in[1]);
            __m128_u t1 = _mm_unpackhi_ps(in[0], in[1]);
            out[0] = _mm_movelh_ps(t0, t0);
            out[1] = _mm_movehl_ps(t0, t0);
            out[2] = t1;
        }

        template <>
        inline void transpose(__m128_u const (&in)[3], __m128_u (&out)[2])
        {
            __m128_u t0 = _mm_unpacklo_ps(in[0], in[1]);
            __m128_u t1 = _mm_unpacklo_ps(in[2], _mm_setzero_ps());
            out[0] = _mm_movelh_ps(t0, t1);
            out[1] = _mm_movehl_ps(t1, t0);
        }

        template <>
        inline void transpose(__m128_u const (&in)[4], __m128_u (&out)[4])
        {
            __m128_u t0 = _mm_unpacklo_ps(in[0], in[1]);
            __m128_u t1 = _mm_unpackhi_ps(in[0], in[1]);
            __m128_u t2 = _mm_unpacklo_ps(in[2], in[3]);
            __m128_u t3 = _mm_unpackhi_ps(in[2], in[3]);
            out[0] = _mm_movelh_ps(t0, t2);
            out[1] = _mm_movehl_ps(t2, t0);
            out[2] = _mm_movelh_ps(t1, t3);
            out[3] = _mm_movehl_ps(t3, t1);
        }

        template <>
        inline void transpose(__m128_u const (&in)[2], __m128_u (&out)[4])
        {
            __m128_u t0 = _mm_unpacklo_ps(in[0], in[1]);
            __m128_u t1 = _mm_unpackhi_ps(in[0], in[1]);
            out[0] = _mm_movelh_ps(t0, t0);
            out[1] = _mm_movehl_ps(t0, t0);
            out[2] = _mm_movelh_ps(t1, t1);
            out[3] = _mm_movehl_ps(t1, t1);
        }

        template <>
        inline void transpose(__m128_u const (&in)[3], __m128_u (&out)[4])
        {
            __m128_u t0 = _mm_unpacklo_ps(in[0], in[1]);
            __m128_u t1 = _mm_unpackhi_ps(in[0], in[1]);
            __m128_u t2 = _mm_unpacklo_ps(in[2], _mm_setzero_ps());
            __m128_u t3 = _mm_unpackhi_ps(in[2], _mm_setzero_ps());

            out[0] = _mm_movelh_ps(t0, t2);
            out[1] = _mm_movehl_ps(t2, t0);
            out[2] = _mm_movelh_ps(t1, t3);
            out[3] = _mm_movehl_ps(t3, t1);
        }

        template <>
        inline void transpose(__m128_u const (&in)[4], __m128_u (&out)[2])
        {
            __m128_u t0 = _mm_unpacklo_ps(in[0], in[1]);
            __m128_u t1 = _mm_unpacklo_ps(in[2], in[3]);
            out[0] = _mm_movelh_ps(t0, t1);
            out[1] = _mm_movehl_ps(t1, t0);
        }

        template <>
        inline void transpose(__m128_u const (&in)[4], __m128_u (&out)[3])
        {
            __m128_u t0 = _mm_unpacklo_ps(in[0], in[1]);
            __m128_u t1 = _mm_unpackhi_ps(in[0], in[1]);
            __m128_u t2 = _mm_unpacklo_ps(in[2], in[3]);
            __m128_u t3 = _mm_unpackhi_ps(in[2], in[3]);
            out[0] = _mm_movelh_ps(t0, t2);
            out[1] = _mm_movehl_ps(t2, t0);
            out[2] = _mm_movelh_ps(t1, t3);
        }

        inline __m128_u determinant(__m128_u const (&m)[2])
        {
            __m128_u a = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 0, 0));
            __m128_u d = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(1, 1, 1, 1));
            __m128_u bc = _mm_mul_ps(_mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(1, 1, 1, 1)),
                                   _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 0)));
    #if defined(AMAL_FMA_ENABLE)
            return _mm_fmsub_ps(a, d, bc);
    #else
            __m128_u ad = _mm_mul_ps(a, d);
            return _mm_sub_ps(ad, bc);
    #endif
        }

        inline __m128_u determinant(__m128_u const (&m)[3])
        {
            // m[0] = [a, b, c, _]
            // m[1] = [d, e, f, _]
            // m[2] = [g, h, i, _]

            // ei - fh
            __m128_u ei_fh = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 2, 2)),  // [f, f, f, f]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 1, 1))), // [h, h, h, h]
                                      _mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 1, 1)),  // [e, e, e, e]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 2, 2)))  // [i, i, i, i]
            ); // [ei - fh, ei - fh, ei - fh, ei - fh]

            __m128_u a = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 0, 0)); // [a, a, a, a]
            __m128_u term0 = _mm_mul_ps(a, ei_fh);

            // di - fg
            __m128_u di_fg = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 0)),  // [d, d, d, d]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 2, 2))), // [i, i, i, i]
                                      _mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 2, 2)),  // [f, f, f, f]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 0, 0)))  // [g, g, g, g]
            );
            __m128_u b = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 1, 1));
            __m128_u term1 = _mm_mul_ps(b, di_fg);

            // dh - eg
            __m128_u dh_eg = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 0)),  // [d, d, d, d]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 1, 1))), // [h, h, h, h]
                                      _mm_mul_ps(_mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 1, 1)),  // [e, e, e, e]
                                                 _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 0, 0)))  // [g, g, g, g]
            );
            __m128_u c = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 2, 2));
            __m128_u term2 = _mm_mul_ps(c, dh_eg);

            return _mm_sub_ps(_mm_add_ps(term0, term2), term1);
        }

        inline __m128_u determinant(__m128_u const (&m)[4])
        {
            __m128_u s0 = _mm_mul_ps(_mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 1, 1, 2)),
                                                         _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(3, 2, 3, 3))),
                                              _mm_mul_ps(_mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(3, 2, 3, 3)),
                                                         _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(0, 1, 1, 2)))),
                                   _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 1)));

            __m128_u s1 =
                _mm_mul_ps(_mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 1, 2)),
                                                 _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(1, 2, 0, 0))),
                                      _mm_movehl_ps(_mm_setzero_ps(),
                                                    _mm_mul_ps(_mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 1, 2)),
                                                               _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(1, 2, 0, 0))))),
                           _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(2, 3, 3, 3)));

            __m128_u det = _mm_add_ps(s0, s1);
            det = _mm_mul_ps(det, _mm_setr_ps(1.f, -1.f, 1.f, -1.f));
            return _mm_mul_ps(det, m[0]);
        }

        inline void inverse_matrix(__m128_u const (&m)[2], __m128_u (&out)[2])
        {
            __m128_u rdet = _mm_rcp_ps(determinant(m));
            out[0] = _mm_mul_ps(_mm_setr_ps(m[1][1], -m[0][1], 0.f, 0.f), rdet);
            out[1] = _mm_mul_ps(_mm_setr_ps(-m[1][0], m[0][0], 0.f, 0.f), rdet);
        }

        inline void inverse_matrix(__m128_u const (&m)[3], __m128_u (&out)[3])
        {
            const __m128_u a = m[0];
            const __m128_u b = m[1];
            const __m128_u c = m[2];

            const __m128_u i0 = cross_yzx(b, c);
            const __m128_u i1 = cross_yzx(c, a);
            const __m128_u i2 = cross_yzx(a, b);

            out[0] = _mm_setr_ps(i0[0], i1[0], i2[0], 0.f);
            out[1] = _mm_setr_ps(i0[1], i1[1], i2[1], 0.f);
            out[2] = _mm_setr_ps(i0[2], i1[2], i2[2], 0.f);

            const __m128_u rdet = _mm_rcp_ps(dot(a, cross_yzx(b, c)));

            out[0] = _mm_mul_ps(out[0], rdet);
            out[1] = _mm_mul_ps(out[1], rdet);
            out[2] = _mm_mul_ps(out[2], rdet);
        }

        inline void inverse_matrix(__m128_u const (&m)[4], __m128_u (&out)[4])
        {
            __m128 v3 = m[3], v2 = m[2], v1 = m[1];

            __m128 A3_3 = _mm_shuffle_ps(v3, v2, _MM_SHUFFLE(3, 3, 3, 3));
            __m128 A2_2 = _mm_shuffle_ps(v3, v2, _MM_SHUFFLE(2, 2, 2, 2));
            __m128 A1_1 = _mm_shuffle_ps(v3, v2, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 A0_0 = _mm_shuffle_ps(v3, v2, _MM_SHUFFLE(0, 0, 0, 0));

            __m128 B3_3 = _mm_shuffle_ps(v2, v1, _MM_SHUFFLE(3, 3, 3, 3));
            __m128 B2_2 = _mm_shuffle_ps(v2, v1, _MM_SHUFFLE(2, 2, 2, 2));
            __m128 B1_1 = _mm_shuffle_ps(v2, v1, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 B0_0 = _mm_shuffle_ps(v2, v1, _MM_SHUFFLE(0, 0, 0, 0));

            __m128 Fac0, Fac1, Fac2, Fac3, Fac4, Fac5;
            {
                __m128 a = _mm_shuffle_ps(A3_3, A3_3, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 b = _mm_shuffle_ps(A2_2, A2_2, _MM_SHUFFLE(2, 0, 0, 0));
    #ifdef AMAL_FMA_ENABLE
                Fac0 = _mm_fmsub_ps(B2_2, a, _mm_mul_ps(b, B3_3));
    #else
                Fac0 = _mm_sub_ps(_mm_mul_ps(B2_2, a), _mm_mul_ps(b, B3_3));
    #endif
            }
            {
                __m128 a = _mm_shuffle_ps(A3_3, A3_3, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 b = _mm_shuffle_ps(A1_1, A1_1, _MM_SHUFFLE(2, 0, 0, 0));
    #ifdef AMAL_FMA_ENABLE
                Fac1 = _mm_fmsub_ps(B1_1, a, _mm_mul_ps(b, B3_3));
    #else
                Fac1 = _mm_sub_ps(_mm_mul_ps(B1_1, a), _mm_mul_ps(b, B3_3));
    #endif
            }
            // Fac2 = B1_1*swp(A2_2) - B2_2*swp(A1_1)
            {
                __m128 a = _mm_shuffle_ps(A2_2, A2_2, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 b = _mm_shuffle_ps(A1_1, A1_1, _MM_SHUFFLE(2, 0, 0, 0));
    #ifdef AMAL_FMA_ENABLE
                Fac2 = _mm_fmsub_ps(B1_1, a, _mm_mul_ps(b, B2_2));
    #else
                Fac2 = _mm_sub_ps(_mm_mul_ps(B1_1, a), _mm_mul_ps(b, B2_2));
    #endif
            }
            // Fac3 = B0_0*swp(A3_3) - B3_3*swp(A0_0)
            {
                __m128 a = _mm_shuffle_ps(A3_3, A3_3, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 b = _mm_shuffle_ps(A0_0, A0_0, _MM_SHUFFLE(2, 0, 0, 0));
    #ifdef AMAL_FMA_ENABLE
                Fac3 = _mm_fmsub_ps(B0_0, a, _mm_mul_ps(b, B3_3));
    #else
                Fac3 = _mm_sub_ps(_mm_mul_ps(B0_0, a), _mm_mul_ps(b, B3_3));
    #endif
            }
            // Fac4 = B0_0*swp(A2_2) - B2_2*swp(A0_0)
            {
                __m128 a = _mm_shuffle_ps(A2_2, A2_2, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 b = _mm_shuffle_ps(A0_0, A0_0, _MM_SHUFFLE(2, 0, 0, 0));
    #ifdef AMAL_FMA_ENABLE
                Fac4 = _mm_fmsub_ps(B0_0, a, _mm_mul_ps(b, B2_2));
    #else
                Fac4 = _mm_sub_ps(_mm_mul_ps(B0_0, a), _mm_mul_ps(b, B2_2));
    #endif
            }
            // Fac5 = B0_0*swp(A1_1) - B1_1*swp(A0_0)
            {
                __m128 a = _mm_shuffle_ps(A1_1, A1_1, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 b = _mm_shuffle_ps(A0_0, A0_0, _MM_SHUFFLE(2, 0, 0, 0));
    #ifdef AMAL_FMA_ENABLE
                Fac5 = _mm_fmsub_ps(B0_0, a, _mm_mul_ps(b, B1_1));
    #else
                Fac5 = _mm_sub_ps(_mm_mul_ps(B0_0, a), _mm_mul_ps(b, B1_1));
    #endif
            }

            const __m128 sign_a = _mm_set_ps(1.f, -1.f, 1.f, -1.f);
            const __m128 sign_b = _mm_set_ps(-1.f, 1.f, -1.f, 1.f);

            __m128 T0 = _mm_shuffle_ps(m[1], m[0], _MM_SHUFFLE(0, 0, 0, 0));
            __m128 T1 = _mm_shuffle_ps(m[1], m[0], _MM_SHUFFLE(1, 1, 1, 1));
            __m128 T2 = _mm_shuffle_ps(m[1], m[0], _MM_SHUFFLE(2, 2, 2, 2));
            __m128 T3 = _mm_shuffle_ps(m[1], m[0], _MM_SHUFFLE(3, 3, 3, 3));

            __m128 V0 = _mm_shuffle_ps(T0, T0, _MM_SHUFFLE(2, 2, 2, 0));
            __m128 V1 = _mm_shuffle_ps(T1, T1, _MM_SHUFFLE(2, 2, 2, 0));
            __m128 V2 = _mm_shuffle_ps(T2, T2, _MM_SHUFFLE(2, 2, 2, 0));
            __m128 V3 = _mm_shuffle_ps(T3, T3, _MM_SHUFFLE(2, 2, 2, 0));

            __m128 inv0 = _mm_mul_ps(
                sign_b, _mm_add_ps(_mm_sub_ps(_mm_mul_ps(V1, Fac0), _mm_mul_ps(V2, Fac1)), _mm_mul_ps(V3, Fac2)));
            __m128 inv1 = _mm_mul_ps(
                sign_a, _mm_add_ps(_mm_sub_ps(_mm_mul_ps(V0, Fac0), _mm_mul_ps(V2, Fac3)), _mm_mul_ps(V3, Fac4)));
            __m128 inv2 = _mm_mul_ps(
                sign_b, _mm_add_ps(_mm_sub_ps(_mm_mul_ps(V0, Fac1), _mm_mul_ps(V1, Fac3)), _mm_mul_ps(V3, Fac5)));
            __m128 inv3 = _mm_mul_ps(
                sign_a, _mm_add_ps(_mm_sub_ps(_mm_mul_ps(V0, Fac2), _mm_mul_ps(V1, Fac4)), _mm_mul_ps(V2, Fac5)));

            __m128 row0 = _mm_shuffle_ps(inv0, inv1, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 row1 = _mm_shuffle_ps(inv2, inv3, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 row2 = _mm_shuffle_ps(row0, row1, _MM_SHUFFLE(2, 0, 2, 0));

            __m128 det = dot(m[0], row2);
            __m128 rdet = _mm_div_ps(_mm_set1_ps(1.f), det);

            out[0] = _mm_mul_ps(inv0, rdet);
            out[1] = _mm_mul_ps(inv1, rdet);
            out[2] = _mm_mul_ps(inv2, rdet);
            out[3] = _mm_mul_ps(inv3, rdet);
        }

    #if defined(AMAL_FMA_ENABLE)
        #define AMAL_FMA_ADD(a, b, c) _mm_fmadd_ps((a), (b), (c))
    #else
        #define AMAL_FMA_ADD(a, b, c) _mm_add_ps(_mm_mul_ps((a), (b)), (c))
    #endif

        inline void rotate(__m128_u const (&m)[4], float angle, __m128_u const &axis4, __m128_u (&out)[4])
        {
            float c = cosf(angle);
            float s = sinf(angle);

            __m128 one = _mm_set1_ps(1.0f);
            __m128 vc = _mm_set1_ps(c);
            __m128 vs = _mm_set1_ps(s);

            __m128 axis = axis4;
    #if defined(__SSE4_1__)
            __m128 len2 = _mm_dp_ps(axis, axis, 0x7F);
            axis = _mm_mul_ps(axis, _mm_rsqrt_ps(len2));
    #else
            __m128 t = _mm_mul_ps(axis, axis);
            __m128 shuf = _mm_shuffle_ps(t, t, _MM_SHUFFLE(2, 3, 0, 1));
            __m128 sum2 = _mm_add_ps(t, shuf);
            __m128 dot = _mm_add_ss(sum2, _mm_shuffle_ps(sum2, sum2, _MM_SHUFFLE(1, 1, 1, 1)));
            dot = _mm_shuffle_ps(dot, dot, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 invlen = _mm_rsqrt_ps(dot);
            const __m128 mask_xyz = _mm_castsi128_ps(_mm_set_epi32(0, -1, -1, -1));
            invlen = _mm_and_ps(invlen, mask_xyz);
            axis = _mm_mul_ps(axis, invlen);
    #endif

            __m128 ic = _mm_sub_ps(one, vc);
            __m128 zero = _mm_setzero_ps();

            __m128 ax = _mm_shuffle_ps(axis, axis, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 ay = _mm_shuffle_ps(axis, axis, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 az = _mm_shuffle_ps(axis, axis, _MM_SHUFFLE(2, 2, 2, 2));

            __m128 ax_ax = _mm_mul_ps(ax, ax);
            __m128 ax_ay = _mm_mul_ps(ax, ay);
            __m128 ax_az = _mm_mul_ps(ax, az);
            __m128 ay_ay = _mm_mul_ps(ay, ay);
            __m128 ay_az = _mm_mul_ps(ay, az);
            __m128 az_az = _mm_mul_ps(az, az);

            __m128 ic_ax_ax = _mm_mul_ps(ic, ax_ax);
            __m128 ic_ax_ay = _mm_mul_ps(ic, ax_ay);
            __m128 ic_ax_az = _mm_mul_ps(ic, ax_az);
            __m128 ic_ay_ay = _mm_mul_ps(ic, ay_ay);
            __m128 ic_ay_az = _mm_mul_ps(ic, ay_az);
            __m128 ic_az_az = _mm_mul_ps(ic, az_az);

            __m128 r00 = AMAL_FMA_ADD(ic_ax_ax, one, vc);
            __m128 r01 = AMAL_FMA_ADD(ic_ax_ay, one, AMAL_FMA_ADD(vs, az, zero));
            __m128 r02 = AMAL_FMA_ADD(ic_ax_az, one, AMAL_FMA_ADD(_mm_sub_ps(zero, vs), ay, zero));

            __m128 r10 = AMAL_FMA_ADD(ic_ax_ay, one, AMAL_FMA_ADD(_mm_sub_ps(zero, vs), az, zero));
            __m128 r11 = AMAL_FMA_ADD(ic_ay_ay, one, vc);
            __m128 r12 = AMAL_FMA_ADD(ic_ay_az, one, AMAL_FMA_ADD(vs, ax, zero));

            __m128 r20 = AMAL_FMA_ADD(ic_ax_az, one, AMAL_FMA_ADD(vs, ay, zero));
            __m128 r21 = AMAL_FMA_ADD(ic_ay_az, one, AMAL_FMA_ADD(_mm_sub_ps(zero, vs), ax, zero));
            __m128 r22 = AMAL_FMA_ADD(ic_az_az, one, vc);

            out[0] = AMAL_FMA_ADD(m[2], r02, AMAL_FMA_ADD(m[1], r01, _mm_mul_ps(m[0], r00)));
            out[1] = AMAL_FMA_ADD(m[2], r12, AMAL_FMA_ADD(m[1], r11, _mm_mul_ps(m[0], r10)));
            out[2] = AMAL_FMA_ADD(m[2], r22, AMAL_FMA_ADD(m[1], r21, _mm_mul_ps(m[0], r20)));

            out[3] = m[3];
        }

        inline void scale(__m128_u const (&in)[4], __m128_u const &v, __m128_u (&out)[4])
        {
            out[0] = _mm_mul_ps(in[0], _mm_shuffle_ps(v, v, _MM_SHUFFLE(0, 0, 0, 0)));
            out[1] = _mm_mul_ps(in[1], _mm_shuffle_ps(v, v, _MM_SHUFFLE(1, 1, 1, 1)));
            out[2] = _mm_mul_ps(in[2], _mm_shuffle_ps(v, v, _MM_SHUFFLE(2, 2, 2, 2)));
            out[3] = in[3];
        }

        inline void shear(__m128_u const (&in)[4], __m128_u const &point, __m128_u const &lxy_lxz, __m128_u const &lyx_lyz,
                          __m128_u const &lzx_lzy, __m128_u (&out)[4])
        {
            float lambda_xy = lxy_lxz[0], lambda_xz = lxy_lxz[1];
            float lambda_yx = lyx_lyz[0], lambda_yz = lyx_lyz[1];
            float lambda_zx = lzx_lzy[0], lambda_zy = lzx_lzy[1];

            float px = point[0], py = point[1], pz = point[2];

            __m128_u col0 = _mm_set_ps(0.0f, lambda_xz, lambda_xy, 1.0f);
            __m128_u col1 = _mm_set_ps(0.0f, lambda_yz, 1.0f, lambda_yx);
            __m128_u col2 = _mm_set_ps(0.0f, 1.0f, lambda_zy, lambda_zx);
            __m128_u col3 = _mm_set_ps(1.0f, -pz * (lambda_zx + lambda_zy), -py * (lambda_yx + lambda_yz),
                                     -px * (lambda_xy + lambda_xz));

            out[0] = _mm_add_ps(
                _mm_add_ps(_mm_mul_ps(in[0], _mm_set1_ps(col0[0])), _mm_mul_ps(in[1], _mm_set1_ps(col1[0]))),
                _mm_add_ps(_mm_mul_ps(in[2], _mm_set1_ps(col2[0])), _mm_mul_ps(in[3], _mm_set1_ps(col3[0]))));

            out[1] = _mm_add_ps(
                _mm_add_ps(_mm_mul_ps(in[0], _mm_set1_ps(col0[1])), _mm_mul_ps(in[1], _mm_set1_ps(col1[1]))),
                _mm_add_ps(_mm_mul_ps(in[2], _mm_set1_ps(col2[1])), _mm_mul_ps(in[3], _mm_set1_ps(col3[1]))));

            out[2] = _mm_add_ps(
                _mm_add_ps(_mm_mul_ps(in[0], _mm_set1_ps(col0[2])), _mm_mul_ps(in[1], _mm_set1_ps(col1[2]))),
                _mm_add_ps(_mm_mul_ps(in[2], _mm_set1_ps(col2[2])), _mm_mul_ps(in[3], _mm_set1_ps(col3[2]))));

            out[3] = _mm_add_ps(
                _mm_add_ps(_mm_mul_ps(in[0], _mm_set1_ps(col0[3])), _mm_mul_ps(in[1], _mm_set1_ps(col1[3]))),
                _mm_add_ps(_mm_mul_ps(in[2], _mm_set1_ps(col2[3])), _mm_mul_ps(in[3], _mm_set1_ps(col3[3]))));
        }

        template <length_t N, length_t M>
        void transpose(__v4si_u const (&in)[N], __v4si_u (&out)[M]);

        inline void transpose(__v4si_u const (&in)[2], __v4si_u (&out)[2])
        {
            __m128i t0 = _mm_unpacklo_epi32(in[0], in[1]);
            out[0] = _mm_unpacklo_epi64(t0, t0); // a0 b0
            out[1] = _mm_unpackhi_epi64(t0, t0); // a1 b1
        }

        inline void transpose(__v4si_u const (&in)[3], __v4si_u (&out)[2])
        {
            __m128i t0 = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i t1 = _mm_unpacklo_epi32(in[2], _mm_setzero_si128());
            out[0] = _mm_unpacklo_epi64(t0, t1);
            out[1] = _mm_unpackhi_epi64(t0, t1);
        }

        inline void transpose(__v4si_u const (&in)[2], __v4si_u (&out)[3])
        {
            __m128i t0 = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i t1 = _mm_unpackhi_epi32(in[0], in[1]);
            out[0] = _mm_unpacklo_epi64(t0, t0);
            out[1] = _mm_unpackhi_epi64(t0, t0);
            out[2] = _mm_unpacklo_epi64(t1, t1);
        }

        inline void transpose(__v4si_u const (&in)[3], __v4si_u (&out)[3])
        {
            __m128i t0 = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i t1 = _mm_unpackhi_epi32(in[0], in[1]);
            out[0] = _mm_unpacklo_epi64(t0, t0);
            out[1] = _mm_unpackhi_epi64(t0, t0);
            out[2] = _mm_unpacklo_epi64(t1, in[2]);
        }

        inline void transpose(__v4si_u const (&in)[4], __v4si_u (&out)[2])
        {
            __m128i ab_lo = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i cd_lo = _mm_unpacklo_epi32(in[2], in[3]);

            out[0] = _mm_unpacklo_epi64(ab_lo, cd_lo);
            out[1] = _mm_unpackhi_epi64(ab_lo, cd_lo);
        }

        inline void transpose(__v4si_u const (&in)[2], __v4si_u (&out)[4])
        {
            __m128i lo = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i hi = _mm_unpackhi_epi32(in[0], in[1]);

            out[0] = _mm_unpacklo_epi64(lo, lo);
            out[1] = _mm_unpackhi_epi64(lo, lo);
            out[2] = _mm_unpacklo_epi64(hi, hi);
            out[3] = _mm_unpackhi_epi64(hi, hi);
        }

        inline void transpose(__v4si_u const (&in)[3], __v4si_u (&out)[4])
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

        inline void transpose(__v4si_u const (&in)[4], __v4si_u (&out)[3])
        {
            __m128i ab_lo = _mm_unpacklo_epi32(in[0], in[1]);
            __m128i ab_hi = _mm_unpackhi_epi32(in[0], in[1]);
            __m128i cd_lo = _mm_unpacklo_epi32(in[2], in[3]);
            __m128i cd_hi = _mm_unpackhi_epi32(in[2], in[3]);

            out[0] = _mm_unpacklo_epi64(ab_lo, cd_lo);
            out[1] = _mm_unpackhi_epi64(ab_lo, cd_lo);
            out[2] = _mm_unpacklo_epi64(ab_hi, cd_hi);
        }

        inline void transpose(__v4si_u const (&in)[4], __v4si_u (&out)[4])
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
        inline void transpose(__v4si_u const (&in)[4], __v4si_u (&out)[4])
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

        inline __v4si_u determinant(__v4si_u const (&m)[2])
        {
            __v4si_u ab = mm_mullo_epi32_compat(m[0], _mm_shuffle_epi32(m[1], _MM_SHUFFLE(0, 0, 1, 1)));
            __v4si_u det = _mm_sub_epi32(ab, _mm_shuffle_epi32(ab, _MM_SHUFFLE(1, 1, 1, 1)));
            return det; // .x = det
        }

        inline __v4si_u determinant(__v4si_u const (&m)[3])
        {
            __v4si_u a = m[0];
            __v4si_u b = m[1];
            __v4si_u c = m[2];

            __v4si_u bc = mm_mullo_epi32_compat(_mm_shuffle_epi32(b, _MM_SHUFFLE(1, 2, 0, 1)),
                                              _mm_shuffle_epi32(c, _MM_SHUFFLE(2, 0, 1, 2)));
            __v4si_u cb = mm_mullo_epi32_compat(_mm_shuffle_epi32(b, _MM_SHUFFLE(2, 0, 1, 2)),
                                              _mm_shuffle_epi32(c, _MM_SHUFFLE(1, 2, 0, 1)));

            __v4si_u cross = _mm_sub_epi32(bc, cb);
            __v4si_u det = mm_mullo_epi32_compat(a, cross);

            return _mm_add_epi32(_mm_shuffle_epi32(det, _MM_SHUFFLE(0, 0, 0, 0)),
                                 _mm_add_epi32(_mm_shuffle_epi32(det, _MM_SHUFFLE(1, 1, 1, 1)),
                                               _mm_shuffle_epi32(det, _MM_SHUFFLE(2, 2, 2, 2))));
        }

        inline __v4si_u determinant(__v4si_u const (&m)[4])
        {
            // s0 = (m2[0,1,1,2] * m3[3,2,3,3] - m2[3,2,3,3] * m3[0,1,1,2]) * m1[0,0,0,1]
            __v4si_u s0 = mm_mullo_epi32_compat(
                _mm_sub_epi32(mm_mullo_epi32_compat(_mm_shuffle_epi32(m[2], _MM_SHUFFLE(0, 1, 1, 2)),
                                                    _mm_shuffle_epi32(m[3], _MM_SHUFFLE(3, 2, 3, 3))),
                              mm_mullo_epi32_compat(_mm_shuffle_epi32(m[2], _MM_SHUFFLE(3, 2, 3, 3)),
                                                    _mm_shuffle_epi32(m[3], _MM_SHUFFLE(0, 1, 1, 2)))),
                _mm_shuffle_epi32(m[1], _MM_SHUFFLE(0, 0, 0, 1)));

            // prod = m2[0,0,1,2] * m3[1,2,0,0]
            __v4si_u prod = mm_mullo_epi32_compat(_mm_shuffle_epi32(m[2], _MM_SHUFFLE(0, 0, 1, 2)),
                                                _mm_shuffle_epi32(m[3], _MM_SHUFFLE(1, 2, 0, 0)));
            __v4si_u upper = _mm_castps_si128(_mm_movehl_ps(_mm_setzero_ps(), _mm_castsi128_ps(prod)));

            __v4si_u s1 =
                mm_mullo_epi32_compat(_mm_sub_epi32(prod, upper), _mm_shuffle_epi32(m[1], _MM_SHUFFLE(2, 3, 3, 3)));

            __v4si_u det = _mm_add_epi32(s0, s1);
            __v4si_u signs = _mm_setr_epi32(1, -1, 1, -1);
            det = mm_mullo_epi32_compat(det, signs);
            return mm_mullo_epi32_compat(det, m[0]);
        }

        template <length_t C>
        inline __v4si_u multiply_matrix(__v4si_u const (&m)[C], __v4si_u const &v)
        {
            __v4si_u result = _mm_setzero_si128();
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
        template <size_t C>
        inline __m256d_u multiply_matrix(__m256d_u const (&m)[C], __m256d_u const &v)
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

            return reinterpret_cast<__m256d_u &>(result);
        }

        template <length_t N, length_t M>
        void transpose(__m256d_u const (&in)[N], __m256d_u (&out)[M]);

        template <>
        inline void transpose(__m256d_u const (&in)[2], __m256d_u (&out)[2])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            out[0] = _mm256_permute2f128_pd(t0, t0, 0x20);
            out[1] = _mm256_permute2f128_pd(t0, t0, 0x31);
        }

        template <>
        inline void transpose(__m256d_u const (&in)[3], __m256d_u (&out)[3])
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
        inline void transpose(__m256d_u const (&in)[2], __m256d_u (&out)[3])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            __m256d t1 = _mm256_unpackhi_pd(in[0], in[1]);
            out[0] = _mm256_permute2f128_pd(t0, t0, 0x20);
            out[1] = _mm256_permute2f128_pd(t0, t0, 0x31);
            out[2] = t1;
        }

        template <>
        inline void transpose(__m256d_u const (&in)[3], __m256d_u (&out)[2])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            __m256d t1 = _mm256_unpacklo_pd(in[2], _mm256_setzero_pd());
            out[0] = _mm256_permute2f128_pd(t0, t1, 0x20);
            out[1] = _mm256_permute2f128_pd(t0, t1, 0x31);
        }

        template <>
        inline void transpose(__m256d_u const (&in)[4], __m256d_u (&out)[4])
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
        inline void transpose(__m256d_u const (&in)[4], __m256d_u (&out)[3])
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
        inline void transpose(__m256d_u const (&in)[3], __m256d_u (&out)[4])
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
        inline void transpose(__m256d_u const (&in)[2], __m256d_u (&out)[4])
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
        inline void transpose(__m256d_u const (&in)[4], __m256d_u (&out)[2])
        {
            __m256d t0 = _mm256_unpacklo_pd(in[0], in[1]);
            __m256d t1 = _mm256_unpackhi_pd(in[0], in[1]);
            __m256d t2 = _mm256_unpacklo_pd(in[2], in[3]);
            __m256d t3 = _mm256_unpackhi_pd(in[2], in[3]);

            out[0] = _mm256_permute2f128_pd(t0, t2, 0x20);
            out[1] = _mm256_permute2f128_pd(t1, t3, 0x20);
        }

        inline __m256d_u determinant(__m256d_u const (&m)[2])
        {
            __m256d mul = _mm256_mul_pd(m[0], _mm256_permute4x64_pd(m[1], 0b01000000)); // m[1]: [a, a, b, b]
            return _mm256_sub_pd(mul, _mm256_permute4x64_pd(mul, 0b01010101));          // mul - [b*c, b*c, b*c, b*c]
        }

        inline __m256d_u determinant(__m256d_u const (&m)[3])
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

        inline __m256d_u determinant(__m256d_u const (&m)[4])
        {
            // s0 = (m2[0,1,1,2] * m3[3,2,3,3] - m2[3,2,3,3] * m3[0,1,1,2]) * m1[0,0,0,1]
            __m256d_u m2a = _mm256_permute4x64_pd(m[2], _MM_SHUFFLE(0, 1, 1, 2));
            __m256d_u m3a = _mm256_permute4x64_pd(m[3], _MM_SHUFFLE(3, 2, 3, 3));
            __m256d_u m2b = _mm256_permute4x64_pd(m[2], _MM_SHUFFLE(3, 2, 3, 3));
            __m256d_u m3b = _mm256_permute4x64_pd(m[3], _MM_SHUFFLE(0, 1, 1, 2));

            __m256d_u s0 = _mm256_mul_pd(_mm256_sub_pd(_mm256_mul_pd(m2a, m3a), _mm256_mul_pd(m2b, m3b)),
                                      _mm256_permute4x64_pd(m[1], _MM_SHUFFLE(0, 0, 0, 1)));

            // s1 = ((m2[0,0,1,2] * m3[1,2,0,0]) - hi(m2*m3)) * m1[2,3,3,3]
            __m256d_u t1 = _mm256_mul_pd(_mm256_permute4x64_pd(m[2], _MM_SHUFFLE(0, 0, 1, 2)),
                                      _mm256_permute4x64_pd(m[3], _MM_SHUFFLE(1, 2, 0, 0)));

            __m128d hi128 = _mm256_extractf128_pd(t1, 1);
            __m256d_u upper = _mm256_castpd128_pd256(hi128);                  // [hi2, hi3, ?, ?]
            upper = _mm256_permute4x64_pd(upper, _MM_SHUFFLE(1, 1, 1, 1)); // [hi3, hi3, hi3, hi3]

            __m256d_u s1 = _mm256_mul_pd(_mm256_sub_pd(t1, upper), _mm256_permute4x64_pd(m[1], _MM_SHUFFLE(2, 3, 3, 3)));

            __m256d_u det = _mm256_add_pd(s0, s1);
            det = _mm256_mul_pd(det, _mm256_setr_pd(1.0, -1.0, 1.0, -1.0));
            return _mm256_mul_pd(det, m[0]);
        }

        inline void inverse_matrix(__m256d_u const (&m)[2], __m256d_u (&out)[2])
        {
            // m[0] = {a00, a10, x, x}
            // m[1] = {a01, a11, x, x}
            __m256d a = m[0];
            __m256d b = m[1];

            // ---------------- det = a00*a11 − a10*a01 ----------------
            __m128d a_low = _mm256_castpd256_pd128(a); // [a00 a10]
            __m128d b_low = _mm256_castpd256_pd128(b); // [a01 a11]

            double a00 = _mm_cvtsd_f64(a_low);                         // lane-0
            double a10 = _mm_cvtsd_f64(_mm_unpackhi_pd(a_low, a_low)); // lane-1
            double a01 = _mm_cvtsd_f64(b_low);
            double a11 = _mm_cvtsd_f64(_mm_unpackhi_pd(b_low, b_low));

            double det = a00 * a11 - a10 * a01;
            __m256d invd = _mm256_set1_pd(1.0 / det);

            // ---------------- adjugate ----------------
            __m256d adj0 = _mm256_setr_pd(+a11, -a01, 0.0, 0.0); // [ d, −b, 0,0 ]
            __m256d adj1 = _mm256_setr_pd(-a10, +a00, 0.0, 0.0); // [−c,  a, 0,0 ]

            out[0] = _mm256_mul_pd(adj0, invd);
            out[1] = _mm256_mul_pd(adj1, invd);
        }

        inline void inverse_matrix(__m256d_u const (&m)[3], __m256d_u (&out)[3])
        {
            __m256d_u a = m[0], b = m[1], c = m[2];
            __m256d_u i0 = cross_yzx(b, c);
            __m256d_u i1 = cross_yzx(c, a);
            __m256d_u i2 = cross_yzx(a, b);
            __m256d_u det = dot(a, i0);
            __m256d_u invd = _mm256_div_pd(_mm256_set1_pd(1.0), det);
            out[0] = _mm256_mul_pd(i0, invd);
            out[1] = _mm256_mul_pd(i1, invd);
            out[2] = _mm256_mul_pd(i2, invd);
        }

        inline void inverse_matrix(__m256d_u const (&m)[4], __m256d_u (&out)[4])
        {
            __m256d_u A = m[0], B = m[1], C = m[2], D = m[3];
        #if defined(__AVX2__)
            __m256d_u AB_lo = _mm256_unpacklo_pd(A, B);
            __m256d_u AB_hi = _mm256_unpackhi_pd(A, B);
            __m256d_u CD_lo = _mm256_unpacklo_pd(C, D);
            __m256d_u CD_hi = _mm256_unpackhi_pd(C, D);

            auto perm = [](const __m256d_u &v) { return _mm256_permute4x64_pd(v, _MM_SHUFFLE(2, 3, 0, 1)); };

            __m256d_u cof0 = _mm256_sub_pd(_mm256_mul_pd(AB_lo, CD_hi), _mm256_mul_pd(perm(AB_lo), perm(CD_hi)));
            __m256d_u cof1 = _mm256_sub_pd(_mm256_mul_pd(AB_hi, CD_lo), _mm256_mul_pd(perm(AB_hi), perm(CD_lo)));
            __m256d_u cof2 = _mm256_sub_pd(_mm256_mul_pd(CD_lo, AB_hi), _mm256_mul_pd(perm(CD_lo), perm(AB_hi)));
            __m256d_u cof3 = _mm256_sub_pd(_mm256_mul_pd(CD_hi, AB_lo), _mm256_mul_pd(perm(CD_hi), perm(AB_lo)));

            __m256d_u t0 = _mm256_unpacklo_pd(cof0, cof1);
            __m256d_u t1 = _mm256_unpacklo_pd(cof2, cof3);
            out[0] = _mm256_permute2f128_pd(t0, t1, 0x20);

            __m256d_u u0 = _mm256_unpackhi_pd(cof0, cof1);
            __m256d_u u1 = _mm256_unpackhi_pd(cof2, cof3);
            out[1] = _mm256_permute2f128_pd(u0, u1, 0x20);
            out[2] = _mm256_permute2f128_pd(t0, t1, 0x31);
            out[3] = _mm256_permute2f128_pd(u0, u1, 0x31);
        #else
            __m256d_u t0 = _mm256_shuffle_pd(A, B, 0x0);
            __m256d_u t1 = _mm256_shuffle_pd(C, D, 0x0);
            __m256d_u t2 = _mm256_shuffle_pd(A, B, 0xF);
            __m256d_u t3 = _mm256_shuffle_pd(C, D, 0xF);

            __m256d_u cof0 = _mm256_sub_pd(_mm256_mul_pd(t0, t3), _mm256_mul_pd(t2, t1));
            __m256d_u cof1 = _mm256_sub_pd(_mm256_mul_pd(t2, t0), _mm256_mul_pd(t1, t3));
            __m256d_u cof2 = _mm256_sub_pd(_mm256_mul_pd(t3, t1), _mm256_mul_pd(t0, t2));
            __m256d_u cof3 = _mm256_sub_pd(_mm256_mul_pd(t1, t3), _mm256_mul_pd(t2, t0));

            out[0] = _mm256_setr_pd(cof0[0], cof1[0], cof2[0], cof3[0]);
            out[1] = _mm256_setr_pd(cof0[1], cof1[1], cof2[1], cof3[1]);
            out[2] = _mm256_setr_pd(cof0[2], cof1[2], cof2[2], cof3[2]);
            out[3] = _mm256_setr_pd(cof0[3], cof1[3], cof2[3], cof3[3]);
        #endif
            const __m256d_u S0 = _mm256_setr_pd(+1, -1, +1, -1);
            const __m256d_u S1 = _mm256_setr_pd(-1, +1, -1, +1);
            out[0] = _mm256_mul_pd(out[0], S0);
            out[1] = _mm256_mul_pd(out[1], S1);
            out[2] = _mm256_mul_pd(out[2], S0);
            out[3] = _mm256_mul_pd(out[3], S1);

            __m256d_u det = dot(A, out[0]);
            __m256d_u invd = _mm256_div_pd(_mm256_set1_pd(1.0), det);
            out[0] = _mm256_mul_pd(out[0], invd);
            out[1] = _mm256_mul_pd(out[1], invd);
            out[2] = _mm256_mul_pd(out[2], invd);
            out[3] = _mm256_mul_pd(out[3], invd);
        }

        inline void rotate(__m256d_u const (&m)[4], double angle, __m256d_u const &axis4, __m256d_u (&out)[4])
        {
            __m256d one = _mm256_set1_pd(1.0);
            __m256d c = _mm256_set1_pd(cos(angle));
            __m256d s = _mm256_set1_pd(sin(angle));

            // normalize axis
            __m256d dot = _mm256_add_pd(_mm256_mul_pd(axis4, axis4),
                                        _mm256_permute4x64_pd(_mm256_mul_pd(axis4, axis4), _MM_SHUFFLE(1, 2, 0, 3)));
            __m256d len = _mm256_sqrt_pd(dot);
            __m256d axis = _mm256_div_pd(axis4, len);

            __m256d ic = _mm256_sub_pd(one, c);

            // splat
            __m256d ax = _mm256_permute4x64_pd(axis, _MM_SHUFFLE(0, 0, 0, 0));
            __m256d ay = _mm256_permute4x64_pd(axis, _MM_SHUFFLE(1, 1, 1, 1));
            __m256d az = _mm256_permute4x64_pd(axis, _MM_SHUFFLE(2, 2, 2, 2));

            // build rotation coefficients
            __m256d r00 = _mm256_add_pd(c, _mm256_mul_pd(ic, _mm256_mul_pd(ax, ax)));
            __m256d r01 = _mm256_add_pd(_mm256_mul_pd(ic, _mm256_mul_pd(ax, ay)), _mm256_mul_pd(s, az));
            __m256d r02 = _mm256_sub_pd(_mm256_mul_pd(ic, _mm256_mul_pd(ax, az)), _mm256_mul_pd(s, ay));

            __m256d r10 = _mm256_sub_pd(_mm256_mul_pd(ic, _mm256_mul_pd(ay, ax)), _mm256_mul_pd(s, az));
            __m256d r11 = _mm256_add_pd(c, _mm256_mul_pd(ic, _mm256_mul_pd(ay, ay)));
            __m256d r12 = _mm256_add_pd(_mm256_mul_pd(ic, _mm256_mul_pd(ay, az)), _mm256_mul_pd(s, ax));

            __m256d r20 = _mm256_add_pd(_mm256_mul_pd(ic, _mm256_mul_pd(az, ax)), _mm256_mul_pd(s, ay));
            __m256d r21 = _mm256_sub_pd(_mm256_mul_pd(ic, _mm256_mul_pd(az, ay)), _mm256_mul_pd(s, ax));
            __m256d r22 = _mm256_add_pd(c, _mm256_mul_pd(ic, _mm256_mul_pd(az, az)));

            out[0] = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(m[0], r00), _mm256_mul_pd(m[1], r01)),
                                   _mm256_mul_pd(m[2], r02));

            out[1] = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(m[0], r10), _mm256_mul_pd(m[1], r11)),
                                   _mm256_mul_pd(m[2], r12));

            out[2] = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(m[0], r20), _mm256_mul_pd(m[1], r21)),
                                   _mm256_mul_pd(m[2], r22));

            out[3] = m[3]; // preserve translation
        }

        inline void scale(__m256d_u const (&in)[4], __m256d_u const &v, __m256d_u (&out)[4])
        {
            out[0] = _mm256_mul_pd(in[0], _mm256_permute4x64_pd(v, _MM_SHUFFLE(0, 0, 0, 0)));
            out[1] = _mm256_mul_pd(in[1], _mm256_permute4x64_pd(v, _MM_SHUFFLE(1, 1, 1, 1)));
            out[2] = _mm256_mul_pd(in[2], _mm256_permute4x64_pd(v, _MM_SHUFFLE(2, 2, 2, 2)));
            out[3] = in[3];
        }

        inline void shear(__m256d_u const (&in)[4], __m256d_u const &point, __m256d_u const &lxy_lxz, __m256d_u const &lyx_lyz,
                          __m256d_u const &lzx_lzy, __m256d_u (&out)[4])
        {
            double lambda_xy = lxy_lxz[0], lambda_xz = lxy_lxz[1];
            double lambda_yx = lyx_lyz[0], lambda_yz = lyx_lyz[1];
            double lambda_zx = lzx_lzy[0], lambda_zy = lzx_lzy[1];

            double px = point[0], py = point[1], pz = point[2];

            __m256d_u col0 = _mm256_set_pd(0.0, lambda_xz, lambda_xy, 1.0);
            __m256d_u col1 = _mm256_set_pd(0.0, lambda_yz, 1.0, lambda_yx);
            __m256d_u col2 = _mm256_set_pd(0.0, 1.0, lambda_zy, lambda_zx);
            __m256d_u col3 = _mm256_set_pd(1.0, -pz * (lambda_zx + lambda_zy), -py * (lambda_yx + lambda_yz),
                                        -px * (lambda_xy + lambda_xz));

            out[0] = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(in[0], _mm256_set1_pd(col0[0])),
                                                 _mm256_mul_pd(in[1], _mm256_set1_pd(col1[0]))),
                                   _mm256_add_pd(_mm256_mul_pd(in[2], _mm256_set1_pd(col2[0])),
                                                 _mm256_mul_pd(in[3], _mm256_set1_pd(col3[0]))));

            out[1] = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(in[0], _mm256_set1_pd(col0[1])),
                                                 _mm256_mul_pd(in[1], _mm256_set1_pd(col1[1]))),
                                   _mm256_add_pd(_mm256_mul_pd(in[2], _mm256_set1_pd(col2[1])),
                                                 _mm256_mul_pd(in[3], _mm256_set1_pd(col3[1]))));

            out[2] = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(in[0], _mm256_set1_pd(col0[2])),
                                                 _mm256_mul_pd(in[1], _mm256_set1_pd(col1[2]))),
                                   _mm256_add_pd(_mm256_mul_pd(in[2], _mm256_set1_pd(col2[2])),
                                                 _mm256_mul_pd(in[3], _mm256_set1_pd(col3[2]))));

            out[3] = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(in[0], _mm256_set1_pd(col0[3])),
                                                 _mm256_mul_pd(in[1], _mm256_set1_pd(col1[3]))),
                                   _mm256_add_pd(_mm256_mul_pd(in[2], _mm256_set1_pd(col2[3])),
                                                 _mm256_mul_pd(in[3], _mm256_set1_pd(col3[3]))));
        }
    #endif
#endif
    } // namespace internal
} // namespace amal