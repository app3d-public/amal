#pragma once

#include "common.hpp"
#include "exponential.hpp"

namespace amal
{
    inline constexpr vec4 rgba8_to_vec4(uint8_t r, uint8_t g, uint8_t b, uint8_t a = 255u)
    {
        return {static_cast<float>(r) / 255.0f, static_cast<float>(g) / 255.0f, static_cast<float>(b) / 255.0f,
                static_cast<float>(a) / 255.0f};
    }

    template <typename T>
    inline constexpr std::enable_if_t<is_floating_point_v<T>, T> hue_to_rgb(T p, T q, T t)
    {
        if (t < static_cast<T>(0)) t += static_cast<T>(1);
        if (t > static_cast<T>(1)) t -= static_cast<T>(1);
        if (t < static_cast<T>(1) / static_cast<T>(6)) return p + (q - p) * static_cast<T>(6) * t;
        if (t < static_cast<T>(1) / static_cast<T>(2)) return q;
        if (t < static_cast<T>(2) / static_cast<T>(3))
            return p + (q - p) * (static_cast<T>(2) / static_cast<T>(3) - t) * static_cast<T>(6);
        return p;
    }

    template <typename T>
    inline constexpr std::enable_if_t<is_floating_point_v<T>, vec<3, T, true>> hsl_to_rgb(T hue_deg, T s, T l)
    {
        const T h = hue_deg / static_cast<T>(360);
        if (s <= static_cast<T>(0)) return {l, l, l};

        const T q = l < static_cast<T>(0.5) ? l * (static_cast<T>(1) + s) : l + s - l * s;
        const T p = static_cast<T>(2) * l - q;
        return {hue_to_rgb(p, q, h + static_cast<T>(1) / static_cast<T>(3)), hue_to_rgb(p, q, h),
                hue_to_rgb(p, q, h - static_cast<T>(1) / static_cast<T>(3))};
    }

    template <typename T>
    inline constexpr std::enable_if_t<is_floating_point_v<T>, vec<4, T, true>> hsl_to_rgba(T hue_deg, T s, T l,
                                                                                           T alpha = static_cast<T>(1))
    {
        const vec<3, T, true> rgb = hsl_to_rgb(hue_deg, s, l);
        return {rgb.x, rgb.y, rgb.z, alpha};
    }

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
