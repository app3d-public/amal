#pragma once

#include "internal/fwd/matrix.hpp"
#include "internal/fwd/vector.hpp"
#include "internal/simd/compare.hpp"

namespace amal
{
    template <length_t N, typename T, bool aligned>
    inline AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC) less_than(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = internal::less_than_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC) less_than(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] < y[i];
        return r;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC) less_than_equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = internal::less_than_equal_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC)
        less_than_equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] <= y[i];
        return r;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC) greater_than(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = internal::greater_than_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC) greater_than(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] > y[i];
        return r;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC) greater_thanEqual(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = internal::greater_than_equal_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC)
        greater_than_equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] >= y[i];
        return r;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC) equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = internal::equal_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC) equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] == y[i];
        return r;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC) not_equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = internal::equal_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC) not_equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] == y[i];
        return r;
    }

    template <length_t N, bool aligned>
    inline constexpr bool any(AMAL_BVEC const &v)
    {
        bool r = false;
        for (length_t i = 0; i < N; ++i) r = r || v[i];
        return r;
    }

    template <length_t N, bool aligned>
    inline constexpr bool all(AMAL_BVEC const &v)
    {
        bool r = false;
        for (length_t i = 0; i < N; ++i) r = r && v[i];
        return r;
    }

    template <length_t N, bool aligned>
    inline constexpr bool invert(AMAL_BVEC const &v)
    {
        bool r = false;
        for (length_t i = 0; i < N; ++i) !v[i];
        return r;
    }

    template <length_t C, length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_VEC(C, bool, aligned) equal(AMAL_MAT_SELF const &a, AMAL_MAT_SELF const &b)
    {
        AMAL_VEC(C, bool, aligned) r(true);
        for (length_t i = 0; i < C; ++i) r[i] = all(equal(a[i], b[i]));
        return r;
    }

    template <length_t C, length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_VEC(C, bool, aligned) equal(AMAL_MAT_SELF const &a, AMAL_MAT_SELF const &b, T epsilon)
    {
        return equal(a, b, AMAL_NVEC(C)(epsilon));
    }

    template <length_t C, length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_VEC(C, bool, aligned)
        equal(AMAL_MAT_SELF const &a, AMAL_MAT_SELF const &b, AMAL_NVEC(C) const &epsilon)
    {
        AMAL_VEC(C, bool, aligned) r(true);
        for (length_t i = 0; i < C; ++i) r[i] = all(equal(a[i], b[i], epsilon[i]));
        return r;
    }

    template <length_t C, length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_VEC(C, bool, aligned) not_equal(AMAL_MAT_SELF const &a, AMAL_MAT_SELF const &b)
    {
        AMAL_VEC(C, bool, aligned) r(true);
        for (length_t i = 0; i < C; ++i) r[i] = any(not_equal(a[i], b[i]));
        return r;
    }

    template <length_t C, length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_VEC(C, bool, aligned)
        not_equal(AMAL_MAT_SELF const &a, AMAL_MAT_SELF const &b, T epsilon)
    {
        return not_equal(a, b, AMAL_NVEC(C)(epsilon));
    }

    template <length_t C, length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_VEC(C, bool, aligned)
        not_equal(AMAL_MAT_SELF const &a, AMAL_MAT_SELF const &b, AMAL_NVEC(C) const &epsilon)
    {
        AMAL_VEC(C, bool, aligned) r(true);
        for (length_t i = 0; i < C; ++i) r[i] = any(not_equal(a[i], b[i], epsilon[i]));
        return r;
    }

    template <length_t C, length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_VEC(C, bool, aligned) equal(AMAL_MAT_SELF const &a, AMAL_MAT_SELF const &b, int max_ulps)
    {
        return equal(a, b, AMAL_VEC(C, int, aligned)(max_ulps));
    }

    template <length_t C, length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_VEC(C, bool, aligned)
        equal(AMAL_MAT_SELF const &a, AMAL_MAT_SELF const &b, AMAL_VEC(C, int, aligned) const &max_ulps)
    {
        AMAL_VEC(C, bool, aligned) r(true);
        for (length_t i = 0; i < C; ++i) r[i] = all(equal(a[i], b[i], max_ulps[i]));
        return r;
    }

    template <length_t C, length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_VEC(C, bool, aligned)
        not_equal(AMAL_MAT_SELF const &a, AMAL_MAT_SELF const &b, int max_ulps)
    {
        return not_equal(a, b, AMAL_VEC(C, int, aligned)(max_ulps));
    }

    template <length_t C, length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_VEC(C, bool, aligned)
        not_equal(AMAL_MAT_SELF const &a, AMAL_MAT_SELF const &b, AMAL_VEC(C, int, aligned) const &max_ulps)
    {
        AMAL_VEC(C, bool, aligned) r(true);
        for (length_t i = 0; i < C; ++i) r[i] = any(not_equal(a[i], b[i], max_ulps[i]));
        return r;
    }
} // namespace amal