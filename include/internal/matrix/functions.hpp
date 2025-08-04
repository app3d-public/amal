#pragma once

#include <cmath>
#include <vec4.hpp>
#include "../fwd/matrix.hpp"
#include "../simd/common.hpp"
#include "../simd/matrix.hpp"

namespace amal
{
    template <length_t C, length_t R, typename T, bool aligned>
    inline AMAL_NMAT_VAL_SIMD(R, C) transpose(AMAL_NMAT(C, R) const &m)
    {
        typename AMAL_NMAT(R, C)::simd_type s[C];
        using simd_t = typename AMAL_NMAT(R, C)::simd_type::value_type;
        simd_t out[C];
        if constexpr (R < 4) __builtin_memset(out, 0, sizeof(out));
        internal::inverse_matrix(*reinterpret_cast<simd_t const(*)[C]>(m.data), out);
        internal::transpose(s, m.data.s);
        return AMAL_NMAT(R, C)(s);
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_NMAT_VAL_NOSIMD(2, 2) transpose(AMAL_NMAT(2, 2) const &m)
    {
        return AMAL_NMAT(2, 2)(m[0][0], m[1][0], m[0][1], m[1][1]);
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_NMAT_VAL_NOSIMD(3, 2) transpose(AMAL_NMAT(2, 3) const &m)
    {
        return AMAL_NMAT(3, 2)(m[0][0], m[1][0], m[0][1], m[1][1], m[0][2], m[1][2]);
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_NMAT_VAL_NOSIMD(4, 2) transpose(AMAL_NMAT(2, 4) const &m)
    {
        return AMAL_NMAT(4, 2)(m[0][0], m[1][0], m[0][1], m[1][1], m[0][2], m[1][2], m[0][3], m[1][3]);
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_NMAT_VAL_NOSIMD(2, 3) transpose(AMAL_NMAT(3, 2) const &m)
    {
        return AMAL_NMAT(2, 3)(m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1]);
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_NMAT_VAL_NOSIMD(3, 3) transpose(AMAL_NMAT(3, 3) const &m)
    {
        return AMAL_NMAT(3, 3)(m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2]);
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_NMAT_VAL_NOSIMD(4, 3) transpose(AMAL_NMAT(3, 4) const &m)
    {
        return AMAL_NMAT(4, 3)(m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2], m[0][3],
                               m[1][3], m[2][3]);
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_NMAT_VAL_NOSIMD(2, 4) transpose(AMAL_NMAT(4, 2) const &m)
    {
        return AMAL_NMAT(2, 4)(m[0][0], m[1][0], m[2][0], m[3][0], m[0][1], m[1][1], m[2][1], m[3][1]);
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_NMAT_VAL_NOSIMD(3, 4) transpose(AMAL_NMAT(4, 3) const &m)
    {
        return AMAL_NMAT(3, 4)(m[0][0], m[1][0], m[2][0], m[3][0], m[0][1], m[1][1], m[2][1], m[3][1], m[0][2], m[1][2],
                               m[2][2], m[3][2]);
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_NMAT_VAL_NOSIMD(4, 4) transpose(AMAL_NMAT(4, 4) const &m)
    {
        return AMAL_NMAT(4, 4)(m[0][0], m[1][0], m[2][0], m[3][0], m[0][1], m[1][1], m[2][1], m[3][1], m[0][2], m[1][2],
                               m[2][2], m[3][2], m[0][3], m[1][3], m[2][3], m[3][3]);
    }

    template <length_t R, length_t C, typename T, bool aligned>
    inline AMAL_TYPE_SIMD(AMAL_NMAT(R, C), T) determinant(AMAL_NMAT(R, C) const &m)
    {
        static_assert(R == C, "Matrix must be square");
        return internal::extract_scalar(internal::determinant(m.data));
    }

    template <length_t R, length_t C, typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(R, C), T) determinant(AMAL_NMAT(R, C) const &m)
    {
        static_assert(R == C, "Matrix must be square");
    }

#if defined(AMAL_FMA_ENABLE)
    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(2, 2), T) determinant(AMAL_NMAT(2, 2) const &m)
    {
        return std::fma(m[0][0], m[1][1], -m[0][1] * m[1][0]);
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(3, 3), T) determinant(AMAL_NMAT(3, 3) const &m)
    {
        return std::fma(m[0][0], std::fma(m[1][1], m[2][2], -m[2][1] * m[1][2]),
                        std::fma(-m[1][0], std::fma(m[0][1], m[2][2], -m[2][1] * m[0][2]),
                                 m[2][0] * (m[0][1] * m[1][2] - m[1][1] * m[0][2])));
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(4, 4), T) determinant(AMAL_NMAT(4, 4) const &m)
    {
        T s0 = std::fma(m[2][2], m[3][3], -m[3][2] * m[2][3]);
        T s1 = std::fma(m[2][1], m[3][3], -m[3][1] * m[2][3]);
        T s2 = std::fma(m[2][1], m[3][2], -m[3][1] * m[2][2]);
        T s3 = std::fma(m[2][0], m[3][3], -m[3][0] * m[2][3]);
        T s4 = std::fma(m[2][0], m[3][2], -m[3][0] * m[2][2]);
        T s5 = std::fma(m[2][0], m[3][1], -m[3][0] * m[2][1]);

        T c0 = std::fma(m[1][1], s0, std::fma(-m[1][2], s1, m[1][3] * s2));
        T c1 = std::fma(-m[1][0], s0, std::fma(m[1][2], s3, -m[1][3] * s4));
        T c2 = std::fma(m[1][0], s1, std::fma(-m[1][1], s3, m[1][3] * s5));
        T c3 = std::fma(-m[1][0], s2, std::fma(m[1][1], s4, -m[1][2] * s5));

        return std::fma(m[0][0], c0, std::fma(m[0][1], c1, std::fma(m[0][2], c2, m[0][3] * c3)));
    }

#else
    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(2, 2), T) determinant(AMAL_NMAT(2, 2) const &m)
    {
        return (m[0][0] * m[1][1]) - (m[0][1] * m[1][0]);
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(3, 3), T) determinant(AMAL_NMAT(3, 3) const &m)
    {
        return (m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])) -
               (m[1][0] * (m[0][1] * m[2][2] - m[2][1] * m[0][2])) +
               (m[2][0] * (m[0][1] * m[1][2] - m[1][1] * m[0][2]));
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(4, 4), T) determinant(AMAL_NMAT(4, 4) const &m)
    {
        T s0 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
        T s1 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
        T s2 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
        T s3 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
        T s4 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
        T s5 = m[2][0] * m[3][1] - m[3][0] * m[2][1];

        AMAL_NVEC(4)
        coff(+(m[1][1] * s0 - m[1][2] * s1 + m[1][3] * s2), -(m[1][0] * s0 - m[1][2] * s3 + m[1][3] * s4),
             +(m[1][0] * s1 - m[1][1] * s3 + m[1][3] * s5), -(m[1][0] * s2 - m[1][1] * s4 + m[1][2] * s5));
        return m[0][0] * coff[0] + m[0][1] * coff[1] + m[0][2] * coff[2] + m[0][3] * coff[3];
    }
#endif

    template <length_t R, length_t C, typename T, bool aligned>
    inline AMAL_NMAT_VAL_SIMD(C, R) inverse_matrix(AMAL_NMAT(R, C) const &m)
    {
        static_assert(is_floating_point_v<T>, "inverse_matrix only supports floating point types");
        static_assert(R == C, "Matrix must be square");
        using simd_t = typename AMAL_NMAT(R, C)::simd_type::value_type;
        simd_t out[C];
        if constexpr (R < 4) __builtin_memset(out, 0, sizeof(out));
        internal::inverse_matrix(*reinterpret_cast<simd_t const(*)[C]>(m.data), out);
        return AMAL_NMAT(R, C)(out);
    }

    template <typename T, bool aligned>
    AMAL_NMAT_VAL_NOSIMD(2, 2)
    inverse_matrix(AMAL_NMAT(2, 2) const &m)
    {
        static_assert(is_floating_point_v<T>, "inverse_matrix only supports floating point types");
        T rdet = static_cast<T>(1) / (+m[0][0] * m[1][1] - m[1][0] * m[0][1]);
        return AMAL_NMAT(2, 2)(+m[1][1] * rdet, -m[0][1] * rdet, -m[1][0] * rdet, +m[0][0] * rdet);
    }

    template <typename T, bool aligned>
    AMAL_NMAT_VAL_NOSIMD(3, 3)
    inverse_matrix(AMAL_NMAT(3, 3) const &m)
    {
        static_assert(is_floating_point_v<T>, "inverse_matrix only supports floating point types");

        T rdet = static_cast<T>(1) / (+m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
                                      m[1][0] * (m[0][1] * m[2][2] - m[2][1] * m[0][2]) +
                                      m[2][0] * (m[0][1] * m[1][2] - m[1][1] * m[0][2]));

        return AMAL_NMAT(3, 3){
                   {+(m[1][1] * m[2][2] - m[2][1] * m[1][2]),  //
                    -(m[0][1] * m[2][2] - m[2][1] * m[0][2]),  //
                    +(m[0][1] * m[1][2] - m[1][1] * m[0][2])}, //

                   {-(m[1][0] * m[2][2] - m[2][0] * m[1][2]),  //
                    +(m[0][0] * m[2][2] - m[2][0] * m[0][2]),  //
                    -(m[0][0] * m[1][2] - m[1][0] * m[0][2])}, //

                   {+(m[1][0] * m[2][1] - m[2][0] * m[1][1]), //
                    -(m[0][0] * m[2][1] - m[2][0] * m[0][1]), //
                    +(m[0][0] * m[1][1] - m[1][0] * m[0][1])} //
               } *
               rdet;
    }

    template <typename T, bool aligned>
    AMAL_NMAT_VAL_NOSIMD(4, 4)
    inverse_matrix(AMAL_NMAT(4, 4) const &m)
    {
        static_assert(is_floating_point_v<T>, "inverse_matrix only supports floating point types");

        const T                                          //
            c00 = m[2][2] * m[3][3] - m[3][2] * m[2][3], //
            c02 = m[1][2] * m[3][3] - m[3][2] * m[1][3], //
            c03 = m[1][2] * m[2][3] - m[2][2] * m[1][3], //

            c04 = m[2][1] * m[3][3] - m[3][1] * m[2][3], //
            c06 = m[1][1] * m[3][3] - m[3][1] * m[1][3], //
            c07 = m[1][1] * m[2][3] - m[2][1] * m[1][3], //

            c08 = m[2][1] * m[3][2] - m[3][1] * m[2][2], //
            c10 = m[1][1] * m[3][2] - m[3][1] * m[1][2], //
            c11 = m[1][1] * m[2][2] - m[2][1] * m[1][2], //

            c12 = m[2][0] * m[3][3] - m[3][0] * m[2][3], //
            c14 = m[1][0] * m[3][3] - m[3][0] * m[1][3], //
            c15 = m[1][0] * m[2][3] - m[2][0] * m[1][3], //

            c16 = m[2][0] * m[3][2] - m[3][0] * m[2][2], //
            c18 = m[1][0] * m[3][2] - m[3][0] * m[1][2], //
            c19 = m[1][0] * m[2][2] - m[2][0] * m[1][2], //

            c20 = m[2][0] * m[3][1] - m[3][0] * m[2][1], //
            c22 = m[1][0] * m[3][1] - m[3][0] * m[1][1], //
            c23 = m[1][0] * m[2][1] - m[2][0] * m[1][1]; //

        const AMAL_NVEC(4) f0(c00, c00, c02, c03);
        const AMAL_NVEC(4) f1(c04, c04, c06, c07);
        const AMAL_NVEC(4) f2(c08, c08, c10, c11);
        const AMAL_NVEC(4) f3(c12, c12, c14, c15);
        const AMAL_NVEC(4) f4(c16, c16, c18, c19);
        const AMAL_NVEC(4) f5(c20, c20, c22, c23);

        const AMAL_NVEC(4) v0(m[1][0], m[0][0], m[0][0], m[0][0]);
        const AMAL_NVEC(4) v1(m[1][1], m[0][1], m[0][1], m[0][1]);
        const AMAL_NVEC(4) v2(m[1][2], m[0][2], m[0][2], m[0][2]);
        const AMAL_NVEC(4) v3(m[1][3], m[0][3], m[0][3], m[0][3]);

        const AMAL_NVEC(4) inv0 = v1 * f0 - v2 * f1 + v3 * f2;
        const AMAL_NVEC(4) inv1 = v0 * f0 - v2 * f3 + v3 * f4;
        const AMAL_NVEC(4) inv2 = v0 * f1 - v1 * f3 + v3 * f5;
        const AMAL_NVEC(4) inv3 = v0 * f2 - v1 * f4 + v2 * f5;

        const AMAL_NVEC(4) sign_a(+1, -1, +1, -1);
        const AMAL_NVEC(4) sign_b(-1, +1, -1, +1);

        const AMAL_NMAT(4, 4) inv{inv0 * sign_a, inv1 * sign_b, inv2 * sign_a, inv3 * sign_b};

        T dot = m[0][0] * inv[0][0] + m[0][1] * inv[1][0] + m[0][2] * inv[2][0] + m[0][3] * inv[3][0];
        return inv * (static_cast<T>(1) / dot);
    }

} // namespace amal