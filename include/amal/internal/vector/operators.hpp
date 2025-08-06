#pragma once

#include <cmath>
#include "../fwd/vector.hpp"
#include "../type_info.hpp"

namespace amal
{
    template <length_t N, typename T, bool aligned, typename U>
    inline AMAL_VEC_REF_SIMD operator+=(AMAL_VEC_SELF &lhs, U scalar)
    {
        lhs.s += static_cast<T>(scalar);
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline constexpr AMAL_VEC_REF_NOSIMD operator+=(AMAL_VEC_SELF &lhs, U scalar)
    {
        lhs.x += static_cast<T>(scalar);
        lhs.y += static_cast<T>(scalar);
        if constexpr (N > 2) lhs.z += static_cast<T>(scalar);
        if constexpr (N > 3) lhs.w += static_cast<T>(scalar);
        return lhs;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_REF_SIMD operator+=(AMAL_VEC_SELF &lhs, AMAL_VEC_SELF const &rhs)
    {
        lhs.s += rhs.s;
        return lhs;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_REF_NOSIMD operator+=(AMAL_VEC_SELF &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        lhs.x += static_cast<T>(rhs.x);
        lhs.y += static_cast<T>(rhs.y);
        if constexpr (NB > 2)
        {
            if constexpr (N > 2) lhs.z += static_cast<T>(rhs.z);
            if constexpr (N > 3) lhs.w += static_cast<T>(rhs.w);
        }
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline AMAL_VEC_REF_SIMD operator-=(AMAL_VEC_SELF &lhs, U scalar)
    {
        lhs.s -= static_cast<T>(scalar);
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline constexpr AMAL_VEC_REF_NOSIMD operator-=(AMAL_VEC_SELF &lhs, U scalar)
    {
        lhs.x -= static_cast<T>(scalar);
        lhs.y -= static_cast<T>(scalar);
        if constexpr (N > 2) lhs.z -= static_cast<T>(scalar);
        if constexpr (N > 3) lhs.w -= static_cast<T>(scalar);
        return lhs;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_REF_SIMD operator-=(AMAL_VEC_SELF &lhs, AMAL_VEC_SELF const &rhs)
    {
        lhs.s -= rhs.s;
        return lhs;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_REF_NOSIMD operator-=(AMAL_VEC_SELF &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        lhs.x -= static_cast<T>(rhs.x);
        lhs.y -= static_cast<T>(rhs.y);
        if constexpr (NB > 2)
        {
            if constexpr (N > 2) lhs.z -= static_cast<T>(rhs.z);
            if constexpr (N > 3) lhs.w -= static_cast<T>(rhs.w);
        }
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline AMAL_VEC_REF_SIMD operator/=(AMAL_VEC_SELF &lhs, U scalar)
    {
        lhs.s /= static_cast<T>(scalar);
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline constexpr AMAL_VEC_REF_NOSIMD operator/=(AMAL_VEC_SELF &lhs, U scalar)
    {
        lhs.x /= static_cast<T>(scalar);
        lhs.y /= static_cast<T>(scalar);
        if constexpr (N > 2) lhs.z /= static_cast<T>(scalar);
        if constexpr (N > 3) lhs.w /= static_cast<T>(scalar);
        return lhs;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_REF_SIMD operator/=(AMAL_VEC_SELF &lhs, AMAL_VEC_SELF const &rhs)
    {
        lhs.s /= rhs.s;
        return lhs;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_REF_NOSIMD operator/=(AMAL_VEC_SELF &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        lhs.x /= static_cast<T>(rhs.x);
        lhs.y /= static_cast<T>(rhs.y);
        if constexpr (NB > 2)
        {
            if constexpr (N > 2) lhs.z /= static_cast<T>(rhs.z);
            if constexpr (N > 3) lhs.w /= static_cast<T>(rhs.w);
        }
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline AMAL_VEC_REF_SIMD operator*=(AMAL_VEC_SELF &lhs, U scalar)
    {
        lhs.s *= static_cast<T>(scalar);
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline constexpr AMAL_VEC_REF_NOSIMD operator*=(AMAL_VEC_SELF &lhs, U scalar)
    {
        lhs.x *= static_cast<T>(scalar);
        lhs.y *= static_cast<T>(scalar);
        if constexpr (N > 2) lhs.z *= static_cast<T>(scalar);
        if constexpr (N > 3) lhs.w *= static_cast<T>(scalar);
        return lhs;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_REF_SIMD operator*=(AMAL_VEC_SELF &lhs, AMAL_VEC_SELF const &rhs)
    {
        lhs.s *= rhs.s;
        return lhs;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_REF_NOSIMD operator*=(AMAL_VEC_SELF &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        lhs.x *= static_cast<T>(rhs.x);
        lhs.y *= static_cast<T>(rhs.y);
        if constexpr (NB > 2)
        {
            if constexpr (N > 2) lhs.z *= static_cast<T>(rhs.z);
            if constexpr (N > 3) lhs.w *= static_cast<T>(rhs.w);
        }
        return lhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF &operator++(AMAL_VEC_SELF &lhs)
    {
        return lhs += 1;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF &operator--(AMAL_VEC_SELF &lhs)
    {
        return lhs -= 1;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator++(AMAL_VEC_SELF &lhs, int)
    {
        vec tmp{lhs};
        ++lhs;
        return tmp;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator--(AMAL_VEC_SELF &lhs, int)
    {
        vec tmp{lhs};
        ++lhs;
        return tmp;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline constexpr AMAL_VEC_SELF &operator%=(AMAL_VEC_SELF &lhs, U scalar)
    {
        if constexpr (std::is_integral_v<T>)
        {
            lhs.x %= static_cast<T>(scalar);
            lhs.y %= static_cast<T>(scalar);
            if constexpr (N > 2) lhs.z %= static_cast<T>(scalar);
            if constexpr (N > 3) lhs.w %= static_cast<T>(scalar);
        }
        else
        {
            lhs.x = std::fmod(lhs.x, static_cast<T>(scalar));
            lhs.y = std::fmod(lhs.y, static_cast<T>(scalar));
            if constexpr (N > 2) lhs.z = std::fmod(lhs.z, static_cast<T>(scalar));
            if constexpr (N > 3) lhs.w = std::fmod(lhs.w, static_cast<T>(scalar));
        }
        return lhs;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF &operator%=(AMAL_VEC_SELF &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        if constexpr (std::is_integral_v<T> && std::is_integral_v<TB>)
        {
            lhs.x %= static_cast<T>(rhs.x);
            lhs.y %= static_cast<T>(rhs.y);
            if constexpr (NB > 2)
            {
                if constexpr (N > 2) lhs.z %= static_cast<T>(rhs.z);
                if constexpr (N > 3) lhs.w %= static_cast<T>(rhs.w);
            }
        }
        else
        {
            lhs.x = std::fmod(lhs.x, static_cast<T>(rhs.x));
            lhs.y = std::fmod(lhs.y, static_cast<T>(rhs.y));
            if constexpr (NB > 2)
            {
                if constexpr (N > 2) lhs.z = std::fmod(lhs.z, static_cast<T>(rhs.z));
                if constexpr (N > 3) lhs.w = std::fmod(lhs.w, static_cast<T>(rhs.w));
            }
        }
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline constexpr AMAL_VEC_SELF &operator&=(AMAL_VEC_SELF &lhs, U scalar)
    {
        if constexpr (std::is_integral_v<T>)
        {
            lhs.x &= static_cast<T>(scalar);
            lhs.y &= static_cast<T>(scalar);
            if constexpr (N > 2) lhs.z &= static_cast<T>(scalar);
            if constexpr (N > 3) lhs.w &= static_cast<T>(scalar);
        }
        else
            static_assert(std::is_integral_v<T>, "Bitwise AND not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_REF_SIMD operator&=(AMAL_VEC_SELF &lhs, AMAL_VEC_SELF const &rhs)
    {
        if constexpr (std::is_integral_v<T>)
            lhs.s &= rhs.s;
        else
            static_assert(std::is_integral_v<T>, "Bitwise AND not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF &operator&=(AMAL_VEC_SELF &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        if constexpr (std::is_integral_v<T> && std::is_integral_v<TB>)
        {
            lhs.x &= static_cast<T>(rhs.x);
            lhs.y &= static_cast<T>(rhs.y);
            if constexpr (NB > 2)
            {
                if constexpr (N > 2) lhs.z &= static_cast<T>(rhs.z);
                if constexpr (N > 3) lhs.w &= static_cast<T>(rhs.w);
            }
        }
        else
            static_assert(std::is_integral_v<T> && std::is_integral_v<TB>,
                          "Bitwise AND not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline constexpr AMAL_VEC_SELF &operator|=(AMAL_VEC_SELF &lhs, U scalar)
    {
        if constexpr (std::is_integral_v<T>)
        {
            lhs.x |= static_cast<T>(scalar);
            lhs.y |= static_cast<T>(scalar);
            if constexpr (N > 2) lhs.z |= static_cast<T>(scalar);
            if constexpr (N > 3) lhs.w |= static_cast<T>(scalar);
        }
        else
            static_assert(std::is_integral_v<T>, "Bitwise OR not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_REF_SIMD operator|=(AMAL_VEC_SELF &lhs, AMAL_VEC_SELF const &rhs)
    {
        if constexpr (std::is_integral_v<T> && std::is_integral_v<T>)
            lhs.s |= rhs.s;
        else
            static_assert(std::is_integral_v<T>, "Bitwise OR not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF &operator|=(AMAL_VEC_SELF &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        if constexpr (std::is_integral_v<T> && std::is_integral_v<T>)
        {
            lhs.x |= static_cast<T>(rhs.x);
            lhs.y |= static_cast<T>(rhs.y);
            if constexpr (NB > 2)
            {
                if constexpr (N > 2) lhs.z |= static_cast<T>(rhs.z);
                if constexpr (N > 3) lhs.w |= static_cast<T>(rhs.w);
            }
        }
        else
            static_assert(std::is_integral_v<T> && std::is_integral_v<TB>,
                          "Bitwise OR not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline constexpr AMAL_VEC_SELF &operator^=(AMAL_VEC_SELF &lhs, U scalar)
    {
        if constexpr (std::is_integral_v<T>)
        {
            lhs.x ^= static_cast<T>(scalar);
            lhs.y ^= static_cast<T>(scalar);
            if constexpr (N > 2) lhs.z ^= static_cast<T>(scalar);
            if constexpr (N > 3) lhs.w ^= static_cast<T>(scalar);
        }
        else
            static_assert(std::is_integral_v<T>, "Bitwise XOR not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_REF_SIMD operator^=(AMAL_VEC_SELF &lhs, AMAL_VEC_SELF const &rhs)
    {
        if constexpr (std::is_integral_v<T> && std::is_integral_v<T>)
            lhs.s ^= rhs.s;
        else
            static_assert(std::is_integral_v<T>, "Bitwise XOR not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF &operator^=(AMAL_VEC_SELF &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        if constexpr (std::is_integral_v<T> && std::is_integral_v<T>)
        {
            lhs.x |= static_cast<T>(rhs.x);
            lhs.y |= static_cast<T>(rhs.y);
            if constexpr (NB > 2)
            {
                if constexpr (N > 2) lhs.z ^= static_cast<T>(rhs.z);
                if constexpr (N > 3) lhs.w ^= static_cast<T>(rhs.w);
            }
        }
        else
            static_assert(std::is_integral_v<T> && std::is_integral_v<TB>,
                          "Bitwise XOR not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline constexpr AMAL_VEC_SELF &operator<<=(AMAL_VEC_SELF &lhs, U scalar)
    {
        if constexpr (std::is_integral_v<T>)
        {
            lhs.x <<= static_cast<T>(scalar);
            lhs.y <<= static_cast<T>(scalar);
            if constexpr (N > 2) lhs.z <<= static_cast<T>(scalar);
            if constexpr (N > 3) lhs.w <<= static_cast<T>(scalar);
        }
        else
            static_assert(std::is_integral_v<T>, "Bitwise SHIFT LEFT not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_REF_SIMD operator<<=(AMAL_VEC_SELF &lhs, AMAL_VEC_SELF const &rhs)
    {
        if constexpr (std::is_integral_v<T> && std::is_integral_v<T>)
            lhs.s <<= rhs.s;
        else
            static_assert(std::is_integral_v<T>, "Bitwise SHIFT LEFT not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF &operator<<=(AMAL_VEC_SELF &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        if constexpr (std::is_integral_v<T> && std::is_integral_v<T>)
        {
            lhs.x <<= static_cast<T>(rhs.x);
            lhs.y <<= static_cast<T>(rhs.y);
            if constexpr (NB > 2)
            {
                if constexpr (N > 2) lhs.z <<= static_cast<T>(rhs.z);
                if constexpr (N > 3) lhs.w <<= static_cast<T>(rhs.w);
            }
        }
        else
            static_assert(std::is_integral_v<T> && std::is_integral_v<TB>,
                          "Bitwise OR not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline constexpr AMAL_VEC_SELF &operator>>=(AMAL_VEC_SELF &lhs, U scalar)
    {
        if constexpr (std::is_integral_v<T>)
        {
            lhs.x >>= static_cast<T>(scalar);
            lhs.y >>= static_cast<T>(scalar);
            if constexpr (N > 2) lhs.z >>= static_cast<T>(scalar);
            if constexpr (N > 3) lhs.w >>= static_cast<T>(scalar);
        }
        else
            static_assert(std::is_integral_v<T>, "Bitwise SHIFT RIGHT not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_REF_SIMD operator>>=(AMAL_VEC_SELF &lhs, AMAL_VEC_SELF const &rhs)
    {
        if constexpr (std::is_integral_v<T> && std::is_integral_v<T>)
            lhs.s >>= rhs.s;
        else
            static_assert(std::is_integral_v<T>, "Bitwise SHIFT RIGHT not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF &operator>>=(AMAL_VEC_SELF &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        if constexpr (std::is_integral_v<T> && std::is_integral_v<T>)
        {
            lhs.x >>= static_cast<T>(rhs.x);
            lhs.y >>= static_cast<T>(rhs.y);
            if constexpr (NB > 2)
            {
                if constexpr (N > 2) lhs.z >>= static_cast<T>(rhs.z);
                if constexpr (N > 3) lhs.w >>= static_cast<T>(rhs.w);
            }
        }
        else
            static_assert(std::is_integral_v<T> && std::is_integral_v<TB>,
                          "Bitwise SHIFT RIGHT not valid for floating-point types");
        return lhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator+(AMAL_VEC_SELF const &rhs)
    {
        return rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator-(AMAL_VEC_SELF const &rhs)
    {
        AMAL_VEC_SELF tmp(0);
        return tmp -= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator+(AMAL_VEC_SELF const &v, T scalar)
    {
        AMAL_VEC_SELF tmp(v);
        return tmp += scalar;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF operator+(AMAL_VEC_SELF const &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        AMAL_VEC_SELF tmp(lhs);
        return tmp += rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator+(T scalar, AMAL_VEC_SELF const &rhs)
    {
        AMAL_VEC_SELF tmp(scalar);
        return scalar += rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator-(AMAL_VEC_SELF const &v, T scalar)
    {
        AMAL_VEC_SELF tmp(v);
        return tmp -= scalar;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF operator-(AMAL_VEC_SELF const &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        AMAL_VEC_SELF tmp(lhs);
        return tmp -= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator-(T scalar, AMAL_VEC_SELF const &rhs)
    {
        AMAL_VEC_SELF tmp(scalar);
        return tmp -= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator/(AMAL_VEC_SELF const &v, T scalar)
    {
        AMAL_VEC_SELF tmp(v);
        return v /= scalar;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF operator/(AMAL_VEC_SELF const &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        AMAL_VEC_SELF tmp(lhs);
        return tmp /= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator/(T scalar, AMAL_VEC_SELF const &rhs)
    {
        AMAL_VEC_SELF tmp(scalar);
        return tmp /= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator*(AMAL_VEC_SELF const &v, T scalar)
    {
        AMAL_VEC_SELF tmp(v);
        return tmp *= scalar;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF operator*(AMAL_VEC_SELF const &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        AMAL_VEC_SELF tmp(lhs);
        return tmp *= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator*(T scalar, AMAL_VEC_SELF const &rhs)
    {
        AMAL_VEC_SELF tmp(scalar);
        return tmp *= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator%(AMAL_VEC_SELF const &v, T scalar)
    {
        AMAL_VEC_SELF tmp(v);
        return tmp %= scalar;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF operator%(AMAL_VEC_SELF const &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        AMAL_VEC_SELF tmp(lhs);
        return tmp %= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator%(T scalar, AMAL_VEC_SELF const &rhs)
    {
        AMAL_VEC_SELF tmp(scalar);
        return tmp %= scalar;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator&(AMAL_VEC_SELF const &v, T scalar)
    {
        AMAL_VEC_SELF tmp(v);
        return tmp &= scalar;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF operator&(AMAL_VEC_SELF const &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        AMAL_VEC_SELF tmp(lhs);
        return lhs &= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator&(T scalar, AMAL_VEC_SELF const &rhs)
    {
        AMAL_VEC_SELF tmp(scalar);
        return tmp &= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator|(AMAL_VEC_SELF const &v, T scalar)
    {
        AMAL_VEC_SELF tmp(v);
        return tmp |= scalar;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF operator|(AMAL_VEC_SELF const &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        AMAL_VEC_SELF tmp(lhs);
        return tmp |= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator|(T scalar, AMAL_VEC_SELF const &rhs)
    {
        AMAL_VEC_SELF tmp(scalar);
        return tmp |= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator^(AMAL_VEC_SELF const &v, T scalar)
    {
        AMAL_VEC_SELF tmp(v);
        return tmp ^= scalar;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF operator^(AMAL_VEC_SELF const &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        AMAL_VEC_SELF tmp(lhs);
        return tmp ^= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator^(T scalar, AMAL_VEC_SELF const &rhs)
    {
        AMAL_VEC_SELF tmp(scalar);
        return tmp ^= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator<<(AMAL_VEC_SELF const &v, T scalar)
    {
        AMAL_VEC_SELF tmp(v);
        return tmp %= scalar;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF operator<<(AMAL_VEC_SELF const &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        AMAL_VEC_SELF tmp(lhs);
        return tmp %= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator<<(T scalar, AMAL_VEC_SELF const &rhs)
    {
        AMAL_VEC_SELF tmp(scalar);
        return tmp %= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator>>(AMAL_VEC_SELF const &v, T scalar)
    {
        AMAL_VEC_SELF tmp(v);
        return tmp %= scalar;
    }

    template <length_t N, typename T, bool aligned, length_t NB, typename TB, bool AB>
    inline constexpr AMAL_VEC_SELF operator>>(AMAL_VEC_SELF const &lhs, AMAL_VEC(NB, TB, AB) const &rhs)
    {
        AMAL_VEC_SELF tmp(lhs);
        return tmp %= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF operator>>(T scalar, AMAL_VEC_SELF const &rhs)
    {
        AMAL_VEC_SELF tmp(scalar);
        return tmp %= rhs;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD operator~(AMAL_VEC_SELF const &v)
    {
        return AMAL_VEC_SELF(~v.s);
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_VAL_NOSIMD operator~(AMAL_VEC_SELF const &v)
    {
        AMAL_VEC_SELF r;
        r.x = ~v.x;
        r.y = ~v.y;
        if constexpr (N > 2) r.z = ~v.z;
        if constexpr (N > 3) r.w = ~v.w;
        return r;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr bool operator==(AMAL_VEC_SELF const &lhs, AMAL_VEC_SELF const &rhs)
    {
        bool r = lhs.x == rhs.x && lhs.y == rhs.y;
        if constexpr (N > 2) r = r && lhs.z == rhs.z;
        if constexpr (N > 3) r = r && lhs.w == rhs.w;
        return r;
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr bool operator!=(AMAL_VEC_SELF const &lhs, AMAL_VEC_SELF const &rhs)
    {
        return !(lhs == rhs);
    }

    template <length_t N, bool aligned>
    inline constexpr AMAL_BVEC_VAL_SIMD operator||(AMAL_BVEC const &lhs, AMAL_BVEC const &rhs)
    {
        return AMAL_BVEC(lhs.s | rhs.s);
    }

    template <length_t N, bool aligned>
    inline constexpr AMAL_BVEC_VAL_NOSIMD operator||(AMAL_BVEC const &lhs, AMAL_BVEC const &rhs)
    {
        AMAL_BVEC r;
        r.x = lhs.x || rhs.x;
        r.y = lhs.y || rhs.y;
        if constexpr (N > 2) r.z = lhs.z || rhs.z;
        if constexpr (N > 3) r.w = lhs.w || rhs.w;
        return r;
    }
} // namespace amal