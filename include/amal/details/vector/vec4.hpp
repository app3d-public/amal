#pragma once

#include <cassert>
#include "../fwd/vector.hpp"

namespace amal
{
    template <typename T, enum Pack P>
    struct vec<4, T, P>
    {
        using value_type = T;
        static inline constexpr const length_t vec_dimension = 4;
        using simd_type = details::simd_type<4, T, P>;

        union
        {
            struct
            {
                T x, y, z, w;
            };
            struct
            {
                T r, g, b, a;
            };
            typename simd_type::value_type s;
            T data[simd_type::max_size];
        };

#include "definition.inl"

        constexpr explicit vec(T scalar) noexcept : data{scalar, scalar, scalar, scalar} {}

        constexpr vec(T a, T b, T c, T d) : data{a, b, c, d} {}

        template <typename X, typename Y, typename Z, typename W>
        inline constexpr vec(X x, Y y, Z z, W w)
            : data{static_cast<T>(x), static_cast<T>(y), static_cast<T>(z), static_cast<T>(w)}
        {
        }

        constexpr vec(const AMAL_VEC(2, T, P) & v, T z = 0, T w = 0) : data{v.x, v.y, z, w} {}

        template <typename U, Pack Q>
        constexpr vec(const AMAL_VEC(3, U, Q) & v, T w)
            : data{static_cast<T>(v.x), static_cast<T>(v.y), static_cast<T>(v.z), w}
        {
        }

        template <typename U, Pack Q>
        constexpr vec(T x, const AMAL_VEC(3, U, Q) & v)
            : data{x, static_cast<T>(v.x), static_cast<T>(v.y), static_cast<T>(v.z)}
        {
        }

        template <typename U, Pack Q>
        constexpr vec(const AMAL_VEC(2, U, Q) & a, const AMAL_VEC(2, U, Q) & b)
            : data{static_cast<T>(a.x), static_cast<T>(a.y), static_cast<T>(b.x), static_cast<T>(b.y)}
        {
        }

        template <typename U, Pack Q>
        constexpr vec(T x, const AMAL_VEC(2, U, Q) & v, T w = 0) : data{x, static_cast<T>(v.x), static_cast<T>(v.y), w}
        {
        }

        template <typename U, Pack Q>
        constexpr vec(const AMAL_VEC(2, U, Q) & v, T z = 0, T w = 0)
            : data{static_cast<T>(v.x), static_cast<T>(v.y), z, w}
        {
        }

        template <typename U, Pack Q>
        constexpr vec(T x, T y, const AMAL_VEC(2, U, Q) & v) : data{x, y, static_cast<T>(v.x), static_cast<T>(v.y)}
        {
        }

        template <typename U, Pack Q, typename = std::enable_if_t<!std::is_same_v<U, T>>>
        constexpr explicit vec(const AMAL_VEC(4, U, Q) & v)
            : data{static_cast<T>(v.x), static_cast<T>(v.y), static_cast<T>(v.z), static_cast<T>(v.w)}
        {
        }
    };

    namespace details
    {
        template <typename T, enum Pack P>
        inline constexpr AMAL_NVEC(4) create_by_call(AMAL_NVEC(4) const &v, T (*call)(T))
        {
            return AMAL_NVEC(4)(call(v.x), call(v.y), call(v.z), call(v.w));
        }

        template <typename T, enum Pack P>
        inline constexpr AMAL_NVEC(4)
            create_by_call(AMAL_NVEC(4) const &v1, AMAL_NVEC(4) const &v2, T (*call)(T, T)) noexcept
        {
            return AMAL_NVEC(4)(call(v1.x, v2.x), call(v1.y, v2.y), call(v1.z, v2.z), call(v1.w, v2.w));
        }
    } // namespace details
} // namespace amal