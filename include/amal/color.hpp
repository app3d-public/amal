#pragma once

#include "common.hpp"
#include "exponential.hpp"

namespace amal
{
    template <typename T>
    inline constexpr std::enable_if_t<is_floating_point_v<T>, T> srgb_to_linear(T c)
    {
        if (c <= static_cast<T>(0.04045)) return c / static_cast<T>(12.92);
        const T base = (c + static_cast<T>(0.055)) / static_cast<T>(1.055);
        if consteval
        {
#if defined(__cpp_lib_constexpr_cmath) && (__cpp_lib_constexpr_cmath >= 202306L)
            return static_cast<T>(std::pow(static_cast<double>(base), 2.4));
#else
            return detail::cx_pow(base, static_cast<T>(2.4));
#endif
        }
        else
        {
            return static_cast<T>(std::pow(static_cast<double>(base), 2.4));
        }
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_NVEC_VAL_SIMD(4) srgb_to_linear(AMAL_NVEC(4) const &v)
    {
        static_assert(is_floating_point_v<T>, "srgb_to_linear only supports floating point vectors");
        if consteval { return internal::create_by_call(v, srgb_to_linear<T>); }
        else
        {
            AMAL_NVEC(4) const lo = v / static_cast<T>(12.92);
            AMAL_NVEC(4) const hi = pow((v + static_cast<T>(0.055)) / static_cast<T>(1.055), AMAL_NVEC(4)(2.4));
            __v4si_u mask = internal::less_than_equal_mask(v.s, AMAL_NVEC(4)(static_cast<T>(0.04045)).s);
            AMAL_NVEC(4) out(internal::mix(hi.s, lo.s, mask));
            out[3] = v[3];
            return out;
        }
    }
    template <typename T, bool aligned>
    inline constexpr AMAL_NVEC_VAL_NOSIMD(4) srgb_to_linear(AMAL_NVEC(4) const &v)
    {
        static_assert(is_floating_point_v<T>, "srgb_to_linear only supports floating point vectors");
        if consteval { return internal::create_by_call(v, srgb_to_linear<T>); }
        else
        {
            AMAL_NVEC(4) out{};
            for (length_t i = 0; i < 3; ++i) out[i] = srgb_to_linear<T>(v[i]);
            out[3] = v[3];
            return out;
        }
    }

    template <typename T>
    inline constexpr std::enable_if_t<is_floating_point_v<T>, T> srgb8_to_linear(T c)
    {
        return srgb_to_linear(c / static_cast<T>(255));
    }

    template <typename T>
    inline constexpr std::enable_if_t<std::is_integral_v<T>, float> srgb8_to_linear(T c)
    {
        return srgb_to_linear(static_cast<float>(c) / 255.0f);
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF srgb8_to_linear(AMAL_VEC_SELF const &v)
    {
        if consteval { return internal::create_by_call(v, srgb8_to_linear<T>); }
        else
        {
            return srgb_to_linear(v / static_cast<T>(255));
        }
    }
} // namespace amal
