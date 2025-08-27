#pragma once

#include <amal/internal/vec2.hpp>
#include <amal/internal/vec3.hpp>
#include <amal/internal/vec4.hpp>
#include <amal/internal/vector/operators.hpp>
#include <functional>
#include <utility>

#ifndef AMAL_DETAIL_SPLAT_VEC_DEFINED
    #define AMAL_DETAIL_SPLAT_VEC_DEFINED

namespace amal::detail
{
    template <length_t N, typename T, bool A, std::size_t... I>
    constexpr amal::vec<N, T, A> splat_vec_impl(const T &v, std::index_sequence<I...>)
    {
        (void)sizeof...(I);
        return amal::vec<N, T, A>{(static_cast<void>(I), v)...};
    }

    template <length_t N, typename T, bool A>
    constexpr amal::vec<N, T, A> splat_vec(const T &v)
    {
        return splat_vec_impl<N, T, A>(v, std::make_index_sequence<N>{});
    }
} // namespace amal::detail
#endif // AMAL_DETAIL_SPLAT_VEC_DEFINED

namespace std
{
    template <amal::length_t N, typename T, bool aligned>
    class numeric_limits<amal::vec<N, T, aligned>>
    {
        using base = std::numeric_limits<T>;

    public:
        static constexpr bool is_specialized = base::is_specialized;
        static constexpr bool is_signed = base::is_signed;
        static constexpr bool is_integer = base::is_integer;
        static constexpr bool is_exact = base::is_exact;
        static constexpr bool is_modulo = base::is_modulo;
        static constexpr bool is_bounded = base::is_bounded;
        static constexpr bool is_iec559 = base::is_iec559;
        static constexpr bool has_infinity = base::has_infinity;
        static constexpr bool has_quiet_NaN = base::has_quiet_NaN;
        static constexpr bool has_signaling_NaN = base::has_signaling_NaN;
        static constexpr bool has_denorm_loss = base::has_denorm_loss;
        static constexpr bool traps = base::traps;
        static constexpr bool tinyness_before = base::tinyness_before;

        static constexpr std::float_round_style round_style = base::round_style;

        static constexpr int digits = base::digits;
        static constexpr int digits10 = base::digits10;
        static constexpr int max_digits10 = base::max_digits10;
        static constexpr int radix = base::radix;
        static constexpr int min_exponent = base::min_exponent;
        static constexpr int min_exponent10 = base::min_exponent10;
        static constexpr int max_exponent = base::max_exponent;
        static constexpr int max_exponent10 = base::max_exponent10;

        static constexpr amal::vec<N, T, aligned> min() noexcept
        {
            return amal::detail::splat_vec<N, T, aligned>(base::min());
        }
        static constexpr amal::vec<N, T, aligned> lowest() noexcept
        {
            return amal::detail::splat_vec<N, T, aligned>(base::lowest());
        }
        static constexpr amal::vec<N, T, aligned> max() noexcept
        {
            return amal::detail::splat_vec<N, T, aligned>(base::max());
        }
        static constexpr amal::vec<N, T, aligned> epsilon() noexcept
        {
            return amal::detail::splat_vec<N, T, aligned>(base::epsilon());
        }
        static constexpr amal::vec<N, T, aligned> round_error() noexcept
        {
            return amal::detail::splat_vec<N, T, aligned>(base::round_error());
        }
        static constexpr amal::vec<N, T, aligned> infinity() noexcept
        {
            if constexpr (base::has_infinity) { return amal::detail::splat_vec<N, T, aligned>(base::infinity()); }
            else { return amal::detail::splat_vec<N, T, aligned>(T{}); }
        }
        static constexpr amal::vec<N, T, aligned> quiet_NaN() noexcept
        {
            if constexpr (base::has_quiet_NaN) { return amal::detail::splat_vec<N, T, aligned>(base::quiet_NaN()); }
            else { return amal::detail::splat_vec<N, T, aligned>(T{}); }
        }
        static constexpr amal::vec<N, T, aligned> signaling_NaN() noexcept
        {
            if constexpr (base::has_signaling_NaN)
                return amal::detail::splat_vec<N, T, aligned>(base::signaling_NaN());
            else
                return amal::detail::splat_vec<N, T, aligned>(T{});
        }
        static constexpr amal::vec<N, T, aligned> denorm_min() noexcept
        {
            return amal::detail::splat_vec<N, T, aligned>(base::denorm_min());
        }
    };

    template <typename T, bool aligned>
    struct hash<amal::vec<2, T, aligned>>
    {
        [[nodiscard]] std::size_t operator()(const amal::vec<2, T, aligned> &k) const noexcept
        {
            if constexpr (std::is_integral_v<T>)
                return ((std::hash<T>()(k.x) ^ (std::hash<T>()(k.y) << 1)) >> 1);
            else
            {
                std::hash<T> hasher{};
                size_t h1 = hasher(k.x);
                size_t h2 = hasher(k.y);
                return h1 ^ (h2 << 1);
            }
        }
    };

    template <typename T, bool aligned>
    struct hash<amal::vec<3, T, aligned>>
    {
        [[nodiscard]] std::size_t operator()(const amal::vec<3, T, aligned> &k) const noexcept
        {
            if constexpr (std::is_integral_v<T>)
                return ((std::hash<int>()(k.x) ^ (std::hash<int>()(k.y) << 1)) >> 1) ^ (std::hash<int>()(k.z) << 1);
            else
            {
                std::hash<T> hasher{};
                size_t h1 = hasher(k.x);
                size_t h2 = hasher(k.y);
                size_t h3 = hasher(k.z);
                return h1 ^ (h2 << 1) ^ (h3 << 2);
            }
        }
    };

    template <typename T, bool aligned>
    struct hash<amal::vec<4, T, aligned>>
    {
        [[nodiscard]] std::size_t operator()(const amal::vec<4, T, aligned> &k) const noexcept
        {
            if constexpr (std::is_integral_v<T>)
                return ((std::hash<int>()(k.x) ^ (std::hash<int>()(k.y) << 1)) >> 1) ^ (std::hash<int>()(k.z) << 1) ^
                       (std::hash<int>()(k.w) << 2);
            else
            {
                std::hash<T> hasher{};
                size_t h1 = hasher(k.x);
                size_t h2 = hasher(k.y);
                size_t h3 = hasher(k.z);
                size_t h4 = hasher(k.w);
                return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
            }
        }
    };
} // namespace std