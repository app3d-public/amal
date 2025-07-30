#pragma once

#include "details/fwd/vector.hpp"
#include "details/simd.hpp"

namespace amal
{
    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC) less_than(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = details::less_than_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC) less_than(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] < y[i];
        return r;
    }

    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC)
        less_than_equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = details::less_than_equal_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC)
        less_than_equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] <= y[i];
        return r;
    }

    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC) greater_than(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = details::greater_than_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC) greater_than(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] > y[i];
        return r;
    }

    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC)
        greater_thanEqual(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = details::greater_than_equal_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC)
        greater_than_equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] >= y[i];
        return r;
    }

    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC) equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = details::equal_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC) equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] == y[i];
        return r;
    }

    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_SIMD(AMAL_BVEC, AMAL_BVEC) not_equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        auto mask = details::equal_mask(x.s, y.s);
        if constexpr (N == 4)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0, mask[3] != 0);
        else if constexpr (N == 3)
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0, mask[2] != 0);
        else
            return AMAL_BVEC(mask[0] != 0, mask[1] != 0);
    }

    template <length_t N, typename T, Pack P>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_BVEC, AMAL_BVEC) not_equal(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        AMAL_BVEC r(true);
        for (length_t i = 0; i < N; ++i) r[i] = x[i] == y[i];
        return r;
    }

    template <length_t N, Pack P>
    inline constexpr bool any(AMAL_BVEC const &v)
    {
        bool r = false;
        for (length_t i = 0; i < N; ++i) r = r || v[i];
        return r;
    }

    template <length_t N, Pack P>
    inline constexpr bool all(AMAL_BVEC const &v)
    {
        bool r = false;
        for (length_t i = 0; i < N; ++i) r = r && v[i];
        return r;
    }

    template <length_t N, Pack P>
    inline constexpr bool invert(AMAL_BVEC const &v)
    {
        bool r = false;
        for (length_t i = 0; i < N; ++i) !v[i];
        return r;
    }
} // namespace amal