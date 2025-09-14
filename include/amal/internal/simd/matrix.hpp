#pragma once

#include <cmath>
#include "../type_info.hpp"
#include "geometric.hpp"

namespace amal
{
    namespace internal
    {
#ifdef __SSE2__
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
            return AMAL_FMA_SUB(a, d, bc);
        }

        inline __m128_u determinant(__m128_u const (&m)[3])
        {
            const __m128 a = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 0, 0)); // a
            const __m128 b = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 1, 1)); // b
            const __m128 c = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 2, 2)); // c

            const __m128 d = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 0)); // d
            const __m128 e = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 1, 1)); // e
            const __m128 f = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 2, 2)); // f

            const __m128 g = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 0, 0)); // g
            const __m128 h = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 1, 1)); // h
            const __m128 i = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 2, 2)); // i

            const __m128 ei_fh = AMAL_FMA_SUB(e, i, _mm_mul_ps(f, h));
            const __m128 di_fg = AMAL_FMA_SUB(d, i, _mm_mul_ps(f, g));
            const __m128 dh_eg = AMAL_FMA_SUB(d, h, _mm_mul_ps(e, g));

            const __m128 term0 = _mm_mul_ps(a, ei_fh);
            const __m128 term1 = _mm_mul_ps(b, di_fg);
            const __m128 term2 = _mm_mul_ps(c, dh_eg);

            return _mm_sub_ps(_mm_add_ps(term0, term2), term1);
        }

        inline __m128_u determinant(__m128_u const (&m)[4])
        {
            const __m128 m2_a = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 1, 1, 2));
            const __m128 m3_a = _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(3, 2, 3, 3));
            const __m128 m2_b = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(3, 2, 3, 3));
            const __m128 m3_b = _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(0, 1, 1, 2));
            const __m128 w0 = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 1));

            const __m128 s0_inner = AMAL_FMA_SUB(m2_a, m3_a, _mm_mul_ps(m2_b, m3_b));
            const __m128 s0 = _mm_mul_ps(s0_inner, w0);

            const __m128 m2_c = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 1, 2));
            const __m128 m3_c = _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(1, 2, 0, 0));
            const __m128 prod = _mm_mul_ps(m2_c, m3_c);
            const __m128 high = _mm_movehl_ps(_mm_setzero_ps(), prod);
            const __m128 w1 = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(2, 3, 3, 3));

            const __m128 s1_inner = AMAL_FMA_SUB(m2_c, m3_c, high);
            const __m128 s1 = _mm_mul_ps(s1_inner, w1);

            __m128 det = _mm_add_ps(s0, s1);
            const __m128 sign_mask = _mm_set_ps(-0.0f, +0.0f, -0.0f, +0.0f);
            det = _mm_xor_ps(det, sign_mask);
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
            // column-major: a=m[0], b=m[1], c=m[2]
            __m128 a = (__m128)m[0];
            __m128 b = (__m128)m[1];
            __m128 c = (__m128)m[2];

            __m128 i0 = (__m128)cross_yzx((__m128_u)b, (__m128_u)c); // cross(b,c)
            __m128 i1 = (__m128)cross_yzx((__m128_u)c, (__m128_u)a); // cross(c,a)
            __m128 i2 = (__m128)cross_yzx((__m128_u)a, (__m128_u)b); // cross(a,b)

            __m128 det = dot((__m128_u)a, (__m128_u)i0);
            det = _mm_shuffle_ps(det, det, 0);

            __m128 rcp = _mm_rcp_ps(det);

            i0 = _mm_mul_ps(i0, rcp);
            i1 = _mm_mul_ps(i1, rcp);
            i2 = _mm_mul_ps(i2, rcp);

            i0 = _mm_move_ss(i0, _mm_setzero_ps());
            i1 = _mm_move_ss(i1, _mm_setzero_ps());
            i2 = _mm_move_ss(i2, _mm_setzero_ps());

            out[0] = (__m128_u)i0;
            out[1] = (__m128_u)i1;
            out[2] = (__m128_u)i2;
        }

        inline void inverse_matrix(__m128_u const (&m)[4], __m128_u (&out)[4]) noexcept
        {
            __m128 c0 = (__m128)m[0];
            __m128 c1 = (__m128)m[1];
            __m128 c2 = (__m128)m[2];
            __m128 c3 = (__m128)m[3];

            __m128 fac0;
            {
                __m128 swp0a = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(3, 3, 3, 3));
                __m128 swp0b = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(2, 2, 2, 2));
                __m128 swp00 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(2, 2, 2, 2));
                __m128 swp01 = _mm_shuffle_ps(swp0a, swp0a, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp02 = _mm_shuffle_ps(swp0b, swp0b, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp03 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(3, 3, 3, 3));
                fac0 = AMAL_FMA_SUB(swp00, swp01, _mm_mul_ps(swp02, swp03));
            }
            __m128 fac1;
            {
                __m128 swp0a = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(3, 3, 3, 3));
                __m128 swp0b = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(1, 1, 1, 1));
                __m128 swp00 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(1, 1, 1, 1));
                __m128 swp01 = _mm_shuffle_ps(swp0a, swp0a, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp02 = _mm_shuffle_ps(swp0b, swp0b, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp03 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(3, 3, 3, 3));
                fac1 = AMAL_FMA_SUB(swp00, swp01, _mm_mul_ps(swp02, swp03));
            }
            __m128 fac2;
            {
                __m128 swp0a = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(2, 2, 2, 2));
                __m128 swp0b = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(1, 1, 1, 1));
                __m128 swp00 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(1, 1, 1, 1));
                __m128 swp01 = _mm_shuffle_ps(swp0a, swp0a, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp02 = _mm_shuffle_ps(swp0b, swp0b, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp03 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(2, 2, 2, 2));
                fac2 = AMAL_FMA_SUB(swp00, swp01, _mm_mul_ps(swp02, swp03));
            }
            __m128 fac3;
            {
                __m128 swp0a = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(3, 3, 3, 3));
                __m128 swp0b = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(0, 0, 0, 0));
                __m128 swp00 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(0, 0, 0, 0));
                __m128 swp01 = _mm_shuffle_ps(swp0a, swp0a, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp02 = _mm_shuffle_ps(swp0b, swp0b, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp03 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(3, 3, 3, 3));
                fac3 = AMAL_FMA_SUB(swp00, swp01, _mm_mul_ps(swp02, swp03));
            }
            __m128 fac4;
            {
                __m128 swp0a = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(2, 2, 2, 2));
                __m128 swp0b = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(0, 0, 0, 0));
                __m128 swp00 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(0, 0, 0, 0));
                __m128 swp01 = _mm_shuffle_ps(swp0a, swp0a, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp02 = _mm_shuffle_ps(swp0b, swp0b, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp03 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(2, 2, 2, 2));
                fac4 = AMAL_FMA_SUB(swp00, swp01, _mm_mul_ps(swp02, swp03));
            }
            __m128 fac5;
            {
                __m128 swp0a = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(1, 1, 1, 1));
                __m128 swp0b = _mm_shuffle_ps(c3, c2, _MM_SHUFFLE(0, 0, 0, 0));
                __m128 swp00 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(0, 0, 0, 0));
                __m128 swp01 = _mm_shuffle_ps(swp0a, swp0a, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp02 = _mm_shuffle_ps(swp0b, swp0b, _MM_SHUFFLE(2, 0, 0, 0));
                __m128 swp03 = _mm_shuffle_ps(c2, c1, _MM_SHUFFLE(1, 1, 1, 1));
                fac5 = AMAL_FMA_SUB(swp00, swp01, _mm_mul_ps(swp02, swp03));
            }

            __m128 tmp0 = _mm_shuffle_ps(c1, c0, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 v0 = _mm_shuffle_ps(tmp0, tmp0, _MM_SHUFFLE(2, 2, 2, 0));

            __m128 tmp1 = _mm_shuffle_ps(c1, c0, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 v1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(2, 2, 2, 0));

            __m128 tmp2 = _mm_shuffle_ps(c1, c0, _MM_SHUFFLE(2, 2, 2, 2));
            __m128 v2 = _mm_shuffle_ps(tmp2, tmp2, _MM_SHUFFLE(2, 2, 2, 0));

            __m128 tmp3 = _mm_shuffle_ps(c1, c0, _MM_SHUFFLE(3, 3, 3, 3));
            __m128 v3 = _mm_shuffle_ps(tmp3, tmp3, _MM_SHUFFLE(2, 2, 2, 0));

            static const __m128 mask_b = _mm_set_ps(+0.0f, -0.0f, +0.0f, -0.0f); // [-,+,-,+] по [x,y,z,w]
            static const __m128 mask_a = _mm_set_ps(-0.0f, +0.0f, -0.0f, +0.0f); // [+,-,+,-]

            __m128 ta = AMAL_FMA_SUB(v1, fac0, _mm_mul_ps(v2, fac1));
            __m128 col0 = AMAL_FMA_ADD(v3, fac2, ta);
            __m128 inv0 = _mm_xor_ps(col0, mask_b);

            __m128 tb = AMAL_FMA_SUB(v0, fac0, _mm_mul_ps(v2, fac3));
            __m128 col1 = AMAL_FMA_ADD(v3, fac4, tb);
            __m128 inv1 = _mm_xor_ps(col1, mask_a);

            __m128 tc = AMAL_FMA_SUB(v0, fac1, _mm_mul_ps(v1, fac3));
            __m128 col2 = AMAL_FMA_ADD(v3, fac5, tc);
            __m128 inv2 = _mm_xor_ps(col2, mask_b);

            __m128 tD = AMAL_FMA_SUB(v0, fac2, _mm_mul_ps(v1, fac4));
            __m128 col3 = AMAL_FMA_ADD(v2, fac5, tD);
            __m128 inv3 = _mm_xor_ps(col3, mask_a);

            __m128 r0 = _mm_shuffle_ps(inv0, inv1, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 r1 = _mm_shuffle_ps(inv2, inv3, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 r2 = _mm_shuffle_ps(r0, r1, _MM_SHUFFLE(2, 0, 2, 0));

            __m128 rcp = _mm_rcp_ps(dot(c0, r2));

            out[0] = (__m128_u)_mm_mul_ps(inv0, rcp);
            out[1] = (__m128_u)_mm_mul_ps(inv1, rcp);
            out[2] = (__m128_u)_mm_mul_ps(inv2, rcp);
            out[3] = (__m128_u)_mm_mul_ps(inv3, rcp);
        }

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

        inline void shear(__m128_u const (&in)[4], __m128_u const &point, __m128_u const &lxy_lxz,
                          __m128_u const &lyx_lyz, __m128_u const &lzx_lzy, __m128_u (&out)[4])
        {
            const float l_xy = lxy_lxz[0], l_xz = lxy_lxz[1];
            const float l_yx = lyx_lyz[0], l_yz = lyx_lyz[1];
            const float l_zx = lzx_lzy[0], l_zy = lzx_lzy[1];

            const float px = point[0], py = point[1], pz = point[2];

            // четвертый столбец (трансляция) заранее:
            const float t_x = -px * (l_xy + l_xz);
            const float t_y = -py * (l_yx + l_yz);
            const float t_z = -pz * (l_zx + l_zy);

            // broadcasts, переиспользуются
            const __m128 s_yx = _mm_set1_ps(l_yx);
            const __m128 s_zy = _mm_set1_ps(l_zy);
            const __m128 s_zx = _mm_set1_ps(l_zx);
            const __m128 s_xy = _mm_set1_ps(l_xy);
            const __m128 s_xz = _mm_set1_ps(l_xz);
            const __m128 s_yz = _mm_set1_ps(l_yz);

            const __m128 s_tx = _mm_set1_ps(t_x);
            const __m128 s_ty = _mm_set1_ps(t_y);
            const __m128 s_tz = _mm_set1_ps(t_z);

            // читаем входные столбцы один раз
            const __m128 c0 = (__m128)in[0];
            const __m128 c1 = (__m128)in[1];
            const __m128 c2 = (__m128)in[2];
            const __m128 c3 = (__m128)in[3];

            // out[0] = in0*1 + in1*l_yx + in2*l_zx + in3*0
            __m128 acc0 = c0; // *1
            acc0 = AMAL_FMA_ADD(c1, s_yx, acc0);
            acc0 = AMAL_FMA_ADD(c2, s_zx, acc0);
            out[0] = (__m128_u)acc0;

            // out[1] = in0*l_xy + in1*1 + in2*l_zy + in3*0
            __m128 acc1 = AMAL_FMA_ADD(c0, s_xy, c1); // c1 *1 + c0*l_xy
            acc1 = AMAL_FMA_ADD(c2, s_zy, acc1);
            out[1] = (__m128_u)acc1;

            // out[2] = in0*l_xz + in1*l_yz + in2*1 + in3*0
            __m128 acc2 = AMAL_FMA_ADD(c0, s_xz, c2); // c2 *1 + c0*l_xz
            acc2 = AMAL_FMA_ADD(c1, s_yz, acc2);
            out[2] = (__m128_u)acc2;

            // out[3] = in0*t_x + in1*t_y + in2*t_z + in3*1
            __m128 acc3 = AMAL_FMA_ADD(c0, s_tx, c3); // c3 *1 + c0*t_x
            acc3 = AMAL_FMA_ADD(c1, s_ty, acc3);
            acc3 = AMAL_FMA_ADD(c2, s_tz, acc3);
            out[3] = (__m128_u)acc3;
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

    #ifdef __AVX__
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
            __m256d_u upper = _mm256_castpd128_pd256(hi128);               // [hi2, hi3, ?, ?]
            upper = _mm256_permute4x64_pd(upper, _MM_SHUFFLE(1, 1, 1, 1)); // [hi3, hi3, hi3, hi3]

            __m256d_u s1 =
                _mm256_mul_pd(_mm256_sub_pd(t1, upper), _mm256_permute4x64_pd(m[1], _MM_SHUFFLE(2, 3, 3, 3)));

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

        inline void inverse_matrix(const __m256d_u (&m)[3], __m256d_u (&out)[3])
        {
            __m256d a = m[0], b = m[1], c = m[2];

            __m256d i0 = (__m256d)cross_yzx((__m256d_u)b, (__m256d_u)c); // col0 adj(a)
            __m256d i1 = (__m256d)cross_yzx((__m256d_u)c, (__m256d_u)a); // col1
            __m256d i2 = (__m256d)cross_yzx((__m256d_u)a, (__m256d_u)b); // col2

            __m256d det = (__m256d)dot((__m256d_u)a, (__m256d_u)i0);
            det = _mm256_permute4x64_pd(det, _MM_SHUFFLE(0, 0, 0, 0));

            __m256d invd = _mm256_div_pd(_mm256_set1_pd(1.0), det);

            i0 = _mm256_mul_pd(i0, invd);
            i1 = _mm256_mul_pd(i1, invd);
            i2 = _mm256_mul_pd(i2, invd);

            const __m256d mask_w0 = _mm256_castsi256_pd(_mm256_set_epi64x(0, -1LL, -1LL, -1LL));
            i0 = _mm256_and_pd(i0, mask_w0);
            i1 = _mm256_and_pd(i1, mask_w0);
            i2 = _mm256_and_pd(i2, mask_w0);

            out[0] = i0;
            out[1] = i1;
            out[2] = i2;
        }

        static inline __m256d mat2_mul(__m256d x, __m256d y)
        {
            __m256d xc0 = _mm256_permute2f128_pd(x, x, 0x00);
            __m256d xc1 = _mm256_permute2f128_pd(x, x, 0x11);
            __m256d yrow0 = _mm256_shuffle_pd(y, y, 0x0);
            __m256d yrow1 = _mm256_shuffle_pd(y, y, 0xF);
            return AMAL_FMA_ADD_PD(xc1, yrow1, _mm256_mul_pd(xc0, yrow0));
        }

        static inline __m256d mat2_inv(__m256d V)
        {
            __m128d lo = _mm256_castpd256_pd128(V);       // [a00 a10]
            __m128d hi = _mm256_extractf128_pd(V, 1);     // [a01 a11]
            __m128d hi_sw = _mm_shuffle_pd(hi, hi, 0b01); // [a11 a01]
            __m128d mul = _mm_mul_pd(lo, hi_sw);          // [a00*a11, a10*a01]
            __m128d det128 = _mm_hsub_pd(mul, mul);       // [det, det]
            __m128d invd128 = _mm_div_sd(_mm_set_sd(1.0), det128);
            __m256d invd = _mm256_broadcast_sd((const double *)&invd128);
            __m256d rev = _mm256_permute4x64_pd(V, 0b00011011);
            static const __m256i sgn = _mm256_set_epi64x(0x0, 0x8000000000000000ull, 0x8000000000000000ull, 0x0);
            __m256d adj = _mm256_xor_pd(rev, _mm256_castsi256_pd(sgn));
            return _mm256_mul_pd(adj, invd);
        }

        static inline __m256d block2x2_top(__m256d c0, __m256d c1)
        {
            __m128d a_lo = _mm256_castpd256_pd128(c0); // [r00 r10]
            __m128d b_lo = _mm256_castpd256_pd128(c1); // [r01 r11]
            return _mm256_set_m128d(b_lo, a_lo);       // [r00,r10,r01,r11]
        }
        static inline __m256d block2x2_bottom(__m256d c0, __m256d c1)
        {
            __m128d a_hi = _mm256_extractf128_pd(c0, 1); // [r20 r30]
            __m128d b_hi = _mm256_extractf128_pd(c1, 1); // [r21 r31]
            return _mm256_set_m128d(b_hi, a_hi);         // [r20,r30,r21,r31]
        }

        static inline __m256d make_col_from_blocks_0(__m256d TL, __m256d BL)
        {
            __m128d tl_lo = _mm256_castpd256_pd128(TL);
            __m128d bl_lo = _mm256_castpd256_pd128(BL);
            return _mm256_set_m128d(bl_lo, tl_lo);
        }
        static inline __m256d make_col_from_blocks_1(__m256d TL, __m256d BL)
        {
            __m128d tl_hi = _mm256_extractf128_pd(TL, 1);
            __m128d bl_hi = _mm256_extractf128_pd(BL, 1);
            return _mm256_set_m128d(bl_hi, tl_hi);
        }
        static inline __m256d make_col_from_blocks_2(__m256d TR, __m256d BR)
        {
            __m128d tr_lo = _mm256_castpd256_pd128(TR);
            __m128d br_lo = _mm256_castpd256_pd128(BR);
            return _mm256_set_m128d(br_lo, tr_lo);
        }
        static inline __m256d make_col_from_blocks_3(__m256d TR, __m256d BR)
        {
            __m128d tr_hi = _mm256_extractf128_pd(TR, 1);
            __m128d br_hi = _mm256_extractf128_pd(BR, 1);
            return _mm256_set_m128d(br_hi, tr_hi);
        }

        inline void inverse_matrix(const __m256d_u (&m)[4], __m256d_u (&out)[4]) noexcept
        {
            __m256d c0 = (__m256d)m[0];
            __m256d c1 = (__m256d)m[1];
            __m256d c2 = (__m256d)m[2];
            __m256d c3 = (__m256d)m[3];

            __m256d a = block2x2_top(c0, c1);
            __m256d b = block2x2_top(c2, c3);
            __m256d c = block2x2_bottom(c0, c1);
            __m256d d = block2x2_bottom(c2, c3);

            // det(a)
            __m128d a_lo = _mm256_castpd256_pd128(a);           // [a00 a10]
            __m128d a_hi = _mm256_extractf128_pd(a, 1);         // [a01 a11]
            __m128d a_hi_sw = _mm_shuffle_pd(a_hi, a_hi, 0b01); // [a11 a01]
            __m128d a_mul = _mm_mul_pd(a_lo, a_hi_sw);          // [a00*a11, a10*a01]
            __m128d det_a128 = _mm_hsub_pd(a_mul, a_mul);       // [deta, deta]
            double det_a;
            _mm_store_sd(&det_a, det_a128);

            if (fabs(det_a) > 1e-18)
            {
                __m256d inv_a = mat2_inv(a);
                __m256d a_inv_b = mat2_mul(inv_a, b);
                __m256d c_inv_a = mat2_mul(c, inv_a);
                __m256d s = _mm256_sub_pd(d, mat2_mul(c_inv_a, b)); // S = D - c*inv_a*B
                __m256d inv_s = mat2_inv(s);

                __m256d tl =
                    _mm256_add_pd(inv_a, mat2_mul(a_inv_b, mat2_mul(inv_s, c_inv_a))); // inv_a + inv_a*b*inv_s*c*inv_a
                __m256d tr = _mm256_sub_pd(_mm256_setzero_pd(), mat2_mul(a_inv_b, inv_s)); // -inv_a*b*inv_s
                __m256d bl = _mm256_sub_pd(_mm256_setzero_pd(), mat2_mul(inv_s, c_inv_a)); // -inv_s*c*inv_a
                __m256d br = inv_s;

                out[0] = (__m256d_u)make_col_from_blocks_0(tl, bl);
                out[1] = (__m256d_u)make_col_from_blocks_1(tl, bl);
                out[2] = (__m256d_u)make_col_from_blocks_2(tr, br);
                out[3] = (__m256d_u)make_col_from_blocks_3(tr, br);
                return;
            }

            __m256d inv_d = mat2_inv(d);
            __m256d b_inv_d = mat2_mul(b, inv_d);
            __m256d inv_dc = mat2_mul(inv_d, c);
            __m256d t = _mm256_sub_pd(a, mat2_mul(b_inv_d, c)); // T = A - B*inv_d*c
            __m256d inv_t = mat2_inv(t);

            __m256d tl = inv_t;
            __m256d tr = _mm256_sub_pd(_mm256_setzero_pd(), mat2_mul(inv_t, b_inv_d));
            __m256d bl = _mm256_sub_pd(_mm256_setzero_pd(), mat2_mul(inv_dc, inv_t));
            __m256d br = _mm256_add_pd(inv_d, mat2_mul(inv_dc, mat2_mul(inv_t, b_inv_d)));

            out[0] = (__m256d_u)make_col_from_blocks_0(tl, bl);
            out[1] = (__m256d_u)make_col_from_blocks_1(tl, bl);
            out[2] = (__m256d_u)make_col_from_blocks_2(tr, br);
            out[3] = (__m256d_u)make_col_from_blocks_3(tr, br);
        }

        inline void rotate(__m256d_u const (&m)[4], double angle, __m256d_u const &axis4, __m256d_u (&out)[4])
        {
            const __m256d one = _mm256_set1_pd(1.0);
            const __m256d c = _mm256_set1_pd(std::cos(angle));
            const __m256d s = _mm256_set1_pd(std::sin(angle));

            const __m256d sq = _mm256_mul_pd(axis4, axis4);
            const __m256d dot = _mm256_add_pd(sq, _mm256_permute4x64_pd(sq, _MM_SHUFFLE(1, 2, 0, 3)));
            const __m256d len = _mm256_sqrt_pd(dot);
            const __m256d axis = _mm256_div_pd(axis4, len);

            const __m256d ic = _mm256_sub_pd(one, c);

            // splat
            const __m256d ax = _mm256_permute4x64_pd(axis, _MM_SHUFFLE(0, 0, 0, 0));
            const __m256d ay = _mm256_permute4x64_pd(axis, _MM_SHUFFLE(1, 1, 1, 1));
            const __m256d az = _mm256_permute4x64_pd(axis, _MM_SHUFFLE(2, 2, 2, 2));

            const __m256d xx = _mm256_mul_pd(ax, ax);
            const __m256d yy = _mm256_mul_pd(ay, ay);
            const __m256d zz = _mm256_mul_pd(az, az);

            const __m256d xy = _mm256_mul_pd(ax, ay);
            const __m256d xz = _mm256_mul_pd(ax, az);
            const __m256d yz = _mm256_mul_pd(ay, az);

            const __m256d saz = _mm256_mul_pd(s, az);
            const __m256d say = _mm256_mul_pd(s, ay);
            const __m256d sax = _mm256_mul_pd(s, ax);

            const __m256d r00 = AMAL_FMA_ADD_PD(ic, xx, c); // c + ic*xx
            const __m256d r11 = AMAL_FMA_ADD_PD(ic, yy, c); // c + ic*yy
            const __m256d r22 = AMAL_FMA_ADD_PD(ic, zz, c); // c + ic*zz

            const __m256d r01 = AMAL_FMA_ADD_PD(ic, xy, saz); // ic*xy + s*az
            const __m256d r02 = AMAL_FMA_SUB_PD(ic, xz, say); // ic*xz - s*ay

            const __m256d r10 = AMAL_FMA_SUB_PD(ic, xy, saz); // ic*xy - s*az
            const __m256d r12 = AMAL_FMA_ADD_PD(ic, yz, sax); // ic*yz + s*ax

            const __m256d r20 = AMAL_FMA_ADD_PD(ic, xz, say); // ic*xz + s*ay
            const __m256d r21 = AMAL_FMA_SUB_PD(ic, yz, sax); // ic*yz - s*ax

            const __m256d t0 = AMAL_FMA_ADD_PD(m[0], r00, _mm256_mul_pd(m[1], r01));
            out[0] = AMAL_FMA_ADD_PD(m[2], r02, t0);

            const __m256d t1 = AMAL_FMA_ADD_PD(m[0], r10, _mm256_mul_pd(m[1], r11));
            out[1] = AMAL_FMA_ADD_PD(m[2], r12, t1);

            const __m256d t2 = AMAL_FMA_ADD_PD(m[0], r20, _mm256_mul_pd(m[1], r21));
            out[2] = AMAL_FMA_ADD_PD(m[2], r22, t2);

            out[3] = m[3]; // preserve translation
        }

        inline void scale(__m256d_u const (&in)[4], __m256d_u const &v, __m256d_u (&out)[4])
        {
            out[0] = _mm256_mul_pd(in[0], _mm256_permute4x64_pd(v, _MM_SHUFFLE(0, 0, 0, 0)));
            out[1] = _mm256_mul_pd(in[1], _mm256_permute4x64_pd(v, _MM_SHUFFLE(1, 1, 1, 1)));
            out[2] = _mm256_mul_pd(in[2], _mm256_permute4x64_pd(v, _MM_SHUFFLE(2, 2, 2, 2)));
            out[3] = in[3];
        }

        inline void shear(const __m256d_u (&in)[4], const __m256d_u &point, const __m256d_u &lxy_lxz,
                          const __m256d_u &lyx_lyz, const __m256d_u &lzx_lzy, __m256d_u (&out)[4])
        {
            const double lambda_xy = lxy_lxz[0], lambda_xz = lxy_lxz[1];
            const double lambda_yx = lyx_lyz[0], lambda_yz = lyx_lyz[1];
            const double lambda_zx = lzx_lzy[0], lambda_zy = lzx_lzy[1];

            const double px = point[0], py = point[1], pz = point[2];

            const __m256d b00 = _mm256_set1_pd(1.0);
            const __m256d b01 = _mm256_set1_pd(lambda_xy);
            const __m256d b02 = _mm256_set1_pd(lambda_xz);
            const __m256d b03 = _mm256_set1_pd(0.0);

            __m256d acc0 = AMAL_FMA_ADD_PD(in[1], b01, _mm256_mul_pd(in[0], b00)); // in0*b00 + in1*b01
            acc0 = AMAL_FMA_ADD_PD(in[2], b02, acc0);                              // + in2*b02
            out[0] = AMAL_FMA_ADD_PD(in[3], b03, acc0);                            // + in3*b03

            const __m256d b10 = _mm256_set1_pd(lambda_yx);
            const __m256d b11 = _mm256_set1_pd(1.0);
            const __m256d b12 = _mm256_set1_pd(lambda_yz);
            const __m256d b13 = _mm256_set1_pd(0.0);

            __m256d acc1 = AMAL_FMA_ADD_PD(in[1], b11, _mm256_mul_pd(in[0], b10));
            acc1 = AMAL_FMA_ADD_PD(in[2], b12, acc1);
            out[1] = AMAL_FMA_ADD_PD(in[3], b13, acc1);

            const __m256d b20 = _mm256_set1_pd(lambda_zx);
            const __m256d b21 = _mm256_set1_pd(lambda_zy);
            const __m256d b22 = _mm256_set1_pd(1.0);
            const __m256d b23 = _mm256_set1_pd(0.0);

            __m256d acc2 = AMAL_FMA_ADD_PD(in[1], b21, _mm256_mul_pd(in[0], b20));
            acc2 = AMAL_FMA_ADD_PD(in[2], b22, acc2);
            out[2] = AMAL_FMA_ADD_PD(in[3], b23, acc2);

            const __m256d b30 = _mm256_set1_pd(-px * (lambda_xy + lambda_xz));
            const __m256d b31 = _mm256_set1_pd(-py * (lambda_yx + lambda_yz));
            const __m256d b32 = _mm256_set1_pd(-pz * (lambda_zx + lambda_zy));
            const __m256d b33 = _mm256_set1_pd(1.0);

            __m256d acc3 = AMAL_FMA_ADD_PD(in[1], b31, _mm256_mul_pd(in[0], b30));
            acc3 = AMAL_FMA_ADD_PD(in[2], b32, acc3);
            out[3] = AMAL_FMA_ADD_PD(in[3], b33, acc3);
        }

    #endif
#endif
    } // namespace internal
} // namespace amal