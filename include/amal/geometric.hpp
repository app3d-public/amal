#pragma once

#include "exponential.hpp"
#include "internal/simd/common.hpp"
#include "internal/simd/geometric.hpp"
#include "vector.hpp"

namespace amal
{
    template <length_t N, typename T, bool aligned>
    inline AMAL_TYPE_SIMD(AMAL_NVEC(N), T) dot(AMAL_NVEC(N) const &v1, AMAL_NVEC(N) const &v2)
    {
        return internal::extract_scalar(internal::dot(v1.s, v2.s));
    }

#if defined(AMAL_FMA_ENABLE)
    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NVEC(2), T) dot(AMAL_NVEC(2) const &v1, AMAL_NVEC(2) const &v2)
    {
        return std::fma(v1.y, v2.y, v1.x * v2.x);
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NVEC(3), T) dot(AMAL_NVEC(3) const &v1, AMAL_NVEC(3) const &v2)
    {
        return std::fma(v1.z, v2.z, std::fma(v1.y, v2.y, v1.x * v2.x));
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NVEC(4), T) dot(AMAL_NVEC(4) const &v1, AMAL_NVEC(4) const &v2)
    {
        return std::fma(v1.w, v2.w, std::fma(v1.z, v2.z, std::fma(v1.y, v2.y, v1.x * v2.x)));
    }
#else
    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NVEC(2), T) dot(AMAL_NVEC(2) const &v1, AMAL_NVEC(2) const &v2)
    {
        AMAL_NVEC(2) tmp(v1 * v2);
        return tmp.x + tmp.y;
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NVEC(3), T) dot(AMAL_NVEC(3) const &v1, AMAL_NVEC(3) const &v2)
    {
        AMAL_NVEC(3) tmp(v1 * v2);
        return tmp.x + tmp.y + tmp.z;
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NVEC(4), T) dot(AMAL_NVEC(4) const &v1, AMAL_NVEC(4) const &v2)
    {
        AMAL_NVEC(4) tmp(v1 * v2);
        return tmp.x + tmp.y + tmp.z + tmp.w;
    }
#endif

    template <typename T>
    inline T length(T x)
    {
        return abs(x);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD length(AMAL_VEC_SELF const &v)
    {
        return internal::sqrt(internal::dot(v.s, v.s));
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD length(AMAL_VEC_SELF const &v)
    {
        return sqrt(dot(v, v));
    }

    template <typename T>
    inline T distance(T x, T y)
    {
        return length(x - y);
    }

    template <length_t N, typename T, bool aligned>
    inline T distance(AMAL_VEC_SELF const &v1, AMAL_VEC_SELF const &v2)
    {
        return length(v1 - v2);
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_SIMD(AMAL_NVEC(3), AMAL_NVEC(3)) cross(AMAL_NVEC(3) const &v1, AMAL_NVEC(3) const &v2)
    {
        return AMAL_NVEC(3)(internal::cross(v1.s, v2.s));
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NVEC(3), AMAL_NVEC(3)) cross(AMAL_NVEC(3) const &v1, AMAL_NVEC(3) const &v2)
    {
        return AMAL_NVEC(3)(v1.y * v2.z - v2.y * v1.z, v1.z * v2.x - v2.z * v1.x, v1.x * v2.y - v2.x * v1.y);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD normalize(AMAL_VEC_SELF const &v)
    {
        static_assert(is_floating_point_v<T>, "Normalize only valid for floating-point types");
        return AMAL_VEC_SELF(v.s * internal::inverse_sqrt(internal::dot(v.s, v.s)));
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD normalize(AMAL_VEC_SELF const &v)
    {
        static_assert(is_floating_point_v<T>, "Normalize only valid for floating-point types");
        return v * inverse_sqrt(dot(v, v));
    }

    template <typename T>
    inline T face_forward(T const &n, T const &i, T const &ref)
    {
        return dot(ref, i) < static_cast<T>(0) ? n : -n;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF face_forward(AMAL_VEC_SELF const &n, AMAL_VEC_SELF const &i, AMAL_VEC_SELF const &ref)
    {
        return face_forward<AMAL_VEC_SELF>(n, i, ref);
    }

    template <typename T>
    inline T reflect(T const &i, T const &n)
    {
        return i - n * dot(n, i) * static_cast<T>(2);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF reflect(AMAL_VEC_SELF const &i, AMAL_VEC_SELF const &n)
    {
        return reflect<AMAL_VEC_SELF>(i, n);
    }

    template <typename T>
    inline T refract(T const &i, T const &n, T const &eta)
    {
        T const dot_product(dot(n, i));
        T const k(static_cast<T>(1) - eta * eta * (static_cast<T>(1) - dot_product * dot_product));
        return (eta * i - (eta * dot_product + sqrt(k)) * n) * static_cast<T>(k >= static_cast<T>(0));
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF refract(AMAL_VEC_SELF const &i, AMAL_VEC_SELF const &n, T const &eta)
    {
        return refract<AMAL_VEC_SELF>(i, n, eta);
    }

    enum class axis
    {
        x,
        y,
        z
    };

    template <typename T, bool aligned>
    inline AMAL_NVEC(3) normal_by_axis(axis axis)
    {
        switch (axis)
        {
            case axis::x:
                return AMAL_NVEC(3)(1, 0, 0);
            case axis::y:
                return AMAL_NVEC(3)(0, 1, 0);
            case axis::z:
                return AMAL_NVEC(3)(0, 0, 1);
            default:
                return AMAL_NVEC(3)(0);
        }
    }
} // namespace amal