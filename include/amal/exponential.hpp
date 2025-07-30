#pragma once

#include "details/simd.hpp"
#include "vector.hpp"

namespace amal
{
    using ::pow;
    template <length_t N, typename T, Pack P>
    inline AMAL_VEC_SELF pow(AMAL_VEC_SELF const &base, AMAL_VEC_SELF const &exponent)
    {
        return details::create_by_call(base, exponent, std::pow);
    }

    using ::exp;
    template <length_t N, typename T, Pack P>
    inline AMAL_VEC_SELF exp(AMAL_VEC_SELF const &v)
    {
        return details::create_by_call(v, std::exp);
    }

    using ::exp2;
    template <length_t N, typename T, Pack P>
    inline AMAL_VEC_SELF exp2(AMAL_VEC_SELF const &v)
    {
        return details::create_by_call(v, std::exp2);
    }

    using ::log;
    template <length_t N, typename T, Pack P>
    inline AMAL_VEC_SELF log(AMAL_VEC_SELF const &v)
    {
        return details::create_by_call(v, std::log);
    }

    using ::log2;
    template <length_t N, typename T, Pack P>
    inline AMAL_VEC_SELF log2(AMAL_VEC_SELF const &v)
    {
        return details::create_by_call(v, std::log2);
    }

    using ::log10;
    template <length_t N, typename T, Pack P>
    inline AMAL_VEC_SELF log10(AMAL_VEC_SELF const &v)
    {
        return details::create_by_call(v, std::log10);
    }

    using ::sqrt;

    template <length_t N, typename T, Pack P>
    inline AMAL_VEC_VAL_SIMD sqrt(AMAL_VEC_SELF const &v)
    {
        return AMAL_VEC_SELF(details::sqrt(v.s));
    }

    template <length_t N, typename T, Pack P>
    inline AMAL_VEC_VAL_NOSIMD sqrt(AMAL_VEC_SELF const &v)
    {
        static_assert(is_floating_point_v<T>, "sqrt only supports floating point types");
        return details::create_by_call(v, std::sqrt);
    }

    template <typename T>
    inline T inverse_sqrt(T x)
    {
        return T(1) / sqrt(x);
    }

    template <length_t N, typename T, Pack P>
    inline AMAL_VEC_VAL_SIMD inverse_sqrt(AMAL_VEC_SELF const &v)
    {
        return AMAL_VEC_SELF(details::inverse_sqrt(v.s));
    }

    template <length_t N, typename T, Pack P>
    inline AMAL_VEC_VAL_NOSIMD inverse_sqrt(AMAL_VEC_SELF const &v)
    {
        return details::create_by_call(v, inverse_sqrt<T>);
    }
} // namespace amal