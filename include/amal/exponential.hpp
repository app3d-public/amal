#pragma once

#include "vector.hpp"
#ifdef AMAL_SIMD_ENABLE
    #include "internal/simd/exponential.hpp"
#endif

namespace amal
{
#if !defined(__cpp_lib_constexpr_cmath) || (__cpp_lib_constexpr_cmath < 202306L)
    namespace detail
    {
        template <typename T>
        inline constexpr T cx_log(T x)
        {
            if (x <= static_cast<T>(0)) return -static_cast<T>(1e30);

            constexpr T ln2 = static_cast<T>(0.6931471805599453094);
            int k = 0;
            while (x > static_cast<T>(2))
            {
                x *= static_cast<T>(0.5);
                ++k;
            }
            while (x < static_cast<T>(1))
            {
                x *= static_cast<T>(2);
                --k;
            }

            const T y = (x - static_cast<T>(1)) / (x + static_cast<T>(1));
            const T y2 = y * y;
            T term = y;
            T sum = term;
            for (int n = 3; n <= 21; n += 2)
            {
                term *= y2;
                sum += term / static_cast<T>(n);
            }
            return static_cast<T>(2) * sum + static_cast<T>(k) * ln2;
        }

        template <typename T>
        inline constexpr T cx_exp(T x)
        {
            constexpr T ln2 = static_cast<T>(0.6931471805599453094);
            int k = 0;
            while (x > ln2)
            {
                x -= ln2;
                ++k;
            }
            while (x < -ln2)
            {
                x += ln2;
                --k;
            }

            T term = static_cast<T>(1);
            T sum = static_cast<T>(1);
            for (int n = 1; n <= 16; ++n)
            {
                term *= x / static_cast<T>(n);
                sum += term;
            }

            while (k > 0)
            {
                sum *= static_cast<T>(2);
                --k;
            }
            while (k < 0)
            {
                sum *= static_cast<T>(0.5);
                ++k;
            }
            return sum;
        }

        template <typename T>
        inline constexpr T cx_pow(T base, T exponent)
        {
            return cx_exp(exponent * cx_log(base));
        }
    } // namespace detail
#endif

    using std::pow;
#ifdef AMAL_SIMD_ENABLE
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD pow(AMAL_VEC_SELF const &base, AMAL_VEC_SELF const &exponent)
    {
        if constexpr (is_floating_point_v<T>)
            return AMAL_VEC_SELF(internal::pow(base.s, exponent.s));
        else
            return internal::create_by_call(base, exponent, static_cast<T (*)(T, T)>(std::pow));
    }
#endif

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD pow(AMAL_VEC_SELF const &base, AMAL_VEC_SELF const &exponent)
    {
        return internal::create_by_call(base, exponent, static_cast<T (*)(T, T)>(std::pow));
    }

    using std::exp;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD exp(AMAL_VEC_SELF const &v)
    {
        if constexpr (is_floating_point_v<T>)
            return AMAL_VEC_SELF(internal::exp(v.s));
        else
            return internal::create_by_call(v, static_cast<T (*)(T)>(std::exp));
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD exp(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::exp);
    }

    using std::exp2;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD exp2(AMAL_VEC_SELF const &v)
    {
        if constexpr (std::is_floating_point_v<T>)
            return AMAL_VEC_SELF(internal::exp2(v.s));
        else
            return internal::create_by_call(v, static_cast<T (*)(T)>(std::exp2));
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD exp2(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::exp2);
    }

    using std::log;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD log(AMAL_VEC_SELF const &v)
    {
        if constexpr (is_floating_point_v<T>)
            return AMAL_VEC_SELF(internal::log(v.s));
        else
            return internal::create_by_call(v, static_cast<T (*)(T)>(std::log));
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD log(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::log);
    }

    using std::log2;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD log2(AMAL_VEC_SELF const &v)
    {
        if constexpr (std::is_floating_point_v<T>)
            return AMAL_VEC_SELF(internal::log2(v.s));
        else
            return internal::create_by_call(v, static_cast<T (*)(T)>(std::log2));
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD log2(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::log2);
    }

    using std::log10;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD log10(AMAL_VEC_SELF const &v)
    {
        if constexpr (std::is_floating_point_v<T>)
            return AMAL_VEC_SELF(internal::log10(v.s));
        else
            return internal::create_by_call(v, static_cast<T (*)(T)>(std::log10));
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD log10(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::log10);
    }

    using std::sqrt;
#ifdef AMAL_SIMD_ENABLE
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD sqrt(AMAL_VEC_SELF const &v)
    {
        return AMAL_VEC_SELF(internal::sqrt(v.s));
    }
#endif

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD sqrt(AMAL_VEC_SELF const &v)
    {
        static_assert(is_floating_point_v<T>, "sqrt only supports floating point types");
        return internal::create_by_call(v, std::sqrt);
    }

    template <typename T>
    inline T inverse_sqrt(T x)
    {
        return T(1) / sqrt(x);
    }

#ifdef AMAL_SIMD_ENABLE
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD inverse_sqrt(AMAL_VEC_SELF const &v)
    {
        return AMAL_VEC_SELF(internal::inverse_sqrt(v.s));
    }
#endif

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD inverse_sqrt(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, inverse_sqrt<T>);
    }
} // namespace amal
