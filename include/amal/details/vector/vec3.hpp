#pragma once

#include <cassert>
#include "../fwd/vector.hpp"

namespace amal
{
    template <typename T, Pack P>
    struct vec<3, T, P>
    {
        using value_type = T;
        static inline constexpr const length_t vec_dimension = 3;
        using simd_type = details::simd_type<3, T, P>;

        union
        {
            struct
            {
                T x, y, z;
            };
            struct
            {
                T r, g, b;
            };
            struct
            {
                T u, v, w;
            };
            typename simd_type::value_type s;
            T data[simd_type::max_size];
        };

#include "definition.inl"

        template <AMAL_CONSTRUCT_SIMD>
        constexpr explicit vec(T scalar) noexcept : data{scalar, scalar, scalar, 0}
        {
        }

        template <AMAL_CONSTRUCT_NOSIMD>
        constexpr explicit vec(T scalar) noexcept : data{scalar, scalar, scalar}
        {
        }

        template <AMAL_CONSTRUCT_SIMD>
        constexpr vec(T a, T b, T c) : data{a, b, c, 0}
        {
        }

        template <AMAL_CONSTRUCT_NOSIMD>
        constexpr vec(T a, T b, T c) : data{a, b, c}
        {
        }

        template <typename X, typename Y, typename Z, AMAL_CONSTRUCT_SIMD>
        inline constexpr vec(X x, Y y, Z z) : data{static_cast<T>(x), static_cast<T>(y), static_cast<T>(z), 0}
        {
        }

        template <typename X, typename Y, typename Z, AMAL_CONSTRUCT_NOSIMD>
        inline constexpr vec(X x, Y y, Z z) : data{static_cast<T>(x), static_cast<T>(y), static_cast<T>(z)}
        {
        }

        template <AMAL_CONSTRUCT_SIMD>
        constexpr vec(const AMAL_VEC(2, T, P) & v, T z) : data{v.x, v.y, z, 0}
        {
        }

        template <AMAL_CONSTRUCT_NOSIMD>
        constexpr vec(const AMAL_VEC(2, T, P) & v, T z) : data{v.x, v.y, z}
        {
        }

        template <AMAL_CONSTRUCT_NOSIMD>
        constexpr explicit vec(const AMAL_VEC(4, T, P) & v) : data{v.x, v.y, v.z, 0}
        {
        }

        template <typename U, enum Pack Q, AMAL_CONSTRUCT_SIMD>
        constexpr explicit vec(const AMAL_VEC(3, U, Q) & v)
            : data{static_cast<T>(v.x), static_cast<T>(v.y), static_cast<T>(v.z), 0}
        {
        }

        template <typename U, enum Pack Q, AMAL_CONSTRUCT_NOSIMD>
        constexpr explicit vec(const AMAL_VEC(3, U, Q) & v)
            : data{static_cast<T>(v.x), static_cast<T>(v.y), static_cast<T>(v.z)}
        {
        }
    };

    namespace details
    {
        template <typename T, enum Pack P>
        inline constexpr AMAL_NVEC(3) create_by_call(AMAL_NVEC(3) const &v, T (*call)(T)) noexcept
        {
            return AMAL_NVEC(3)(call(v.x), call(v.y), call(v.z));
        }

        template <typename T, enum Pack P, typename Func>
        inline constexpr AMAL_NVEC(3)
            create_by_call(AMAL_NVEC(3) const &v1, AMAL_NVEC(3) const &v2, Func&& call) noexcept
        {
            return AMAL_NVEC(3)(call(v1.x, v2.x), call(v1.y, v2.y), call(v1.z, v2.z));
        }
    } // namespace details
} // namespace amal