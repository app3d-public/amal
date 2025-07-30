#pragma once

#include <cassert>
#include "../fwd/vector.hpp"

namespace amal
{
    template <typename T, Pack P>
    struct vec<2, T, P>
    {
        using value_type = T;
        static inline constexpr const length_t vec_dimension = 2;
        using simd_type = details::simd_type<2, T, P>;

        union
        {
            struct
            {
                T x, y;
            };
            struct
            {
                T u, v;
            };
            typename simd_type::value_type s;
            T data[simd_type::max_size];
        };

#include "definition.inl"

        template <AMAL_CONSTRUCT_SIMD>
        constexpr explicit vec(T scalar) noexcept : data{scalar, scalar, 0, 0}
        {
        }

        template <AMAL_CONSTRUCT_NOSIMD>
        constexpr explicit vec(T scalar) noexcept : data{scalar, scalar}
        {
        }

        template <AMAL_CONSTRUCT_SIMD>
        constexpr vec(T a, T b) : data{a, b, 0, 0}
        {
        }

        template <AMAL_CONSTRUCT_NOSIMD>
        constexpr vec(T a, T b) : data{a, b}
        {
        }

        template <typename X, typename Y, AMAL_CONSTRUCT_SIMD>
        inline constexpr vec(X x, Y y) : data{static_cast<T>(x), static_cast<T>(y), 0, 0}
        {
        }

        template <typename X, typename Y, typename Z, AMAL_CONSTRUCT_NOSIMD>
        inline constexpr vec(X x, Y y) : data{static_cast<T>(x), static_cast<T>(y)}
        {
        }

        template <AMAL_CONSTRUCT_NOSIMD>
        constexpr explicit vec(const AMAL_VEC(3, T, P) & v) : data{v.x, v.y, 0, 0}
        {
        }

        template <AMAL_CONSTRUCT_NOSIMD>
        constexpr explicit vec(const AMAL_VEC(4, T, P) & v) : data{v.x, v.y, 0, 0}
        {
        }

        template <typename U, enum Pack Q, AMAL_CONSTRUCT_SIMD>
        constexpr explicit vec(const AMAL_VEC(2, U, Q) & v) : data{static_cast<T>(v.x), static_cast<T>(v.y), 0, 0}
        {
        }

        template <typename U, enum Pack Q, AMAL_CONSTRUCT_NOSIMD>
        constexpr explicit vec(const AMAL_VEC(2, U, Q) & v) : data{static_cast<T>(v.x), static_cast<T>(v.y)}
        {
        }
    };

    namespace details
    {
        template <typename T, enum Pack P>
        inline constexpr AMAL_NVEC(2) create_by_call(AMAL_NVEC(2) const &v, T (*call)(T))
        {
            return AMAL_NVEC(2)(call(v.x), call(v.y));
        }

        template <typename T, enum Pack P>
        inline constexpr AMAL_NVEC(2)
            create_by_call(AMAL_NVEC(2) const &v1, AMAL_NVEC(2) const &v2, T (*call)(T, T)) noexcept
        {
            return AMAL_NVEC(2)(call(v1.x, v2.x), call(v1.y, v2.y));
        }
    } // namespace details
} // namespace amal