#pragma once

#include <limits>
#include <utility>

#include <amal/internal/mat2x2.hpp>
#include <amal/internal/mat2x3.hpp>
#include <amal/internal/mat2x4.hpp>
#include <amal/internal/mat3x2.hpp>
#include <amal/internal/mat3x3.hpp>
#include <amal/internal/mat3x4.hpp>
#include <amal/internal/mat4x2.hpp>
#include <amal/internal/mat4x3.hpp>
#include <amal/internal/mat4x4.hpp>
#include <amal/internal/matrix/operators.hpp>

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
#endif

#ifndef AMAL_DETAIL_SPLAT_MAT_DEFINED
    #define AMAL_DETAIL_SPLAT_MAT_DEFINED
namespace amal::detail
{
    template <length_t C, length_t R, typename T, bool A, std::size_t... J>
    constexpr amal::mat<C, R, T, A> splat_mat_impl(const T &v, std::index_sequence<J...>)
    {
        (void)sizeof...(J);
        return amal::mat<C, R, T, A>{(static_cast<void>(J), splat_vec<R, T, A>(v))...};
    }

    template <length_t C, length_t R, typename T, bool A>
    constexpr amal::mat<C, R, T, A> splat_mat(const T &v)
    {
        return splat_mat_impl<C, R, T, A>(v, std::make_index_sequence<C>{});
    }
} // namespace amal::detail
#endif // AMAL_DETAIL_SPLAT_MAT_DEFINED

#ifndef AMAL_STD_NUMERIC_LIMITS_MAT_DEFINED
    #define AMAL_STD_NUMERIC_LIMITS_MAT_DEFINED
namespace std
{

    template <amal::length_t C, amal::length_t R, typename T, bool aligned>
    class numeric_limits<amal::mat<C, R, T, aligned>>
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

        static constexpr amal::mat<C, R, T, aligned> min() noexcept
        {
            return amal::detail::splat_mat<C, R, T, aligned>(base::min());
        }
        static constexpr amal::mat<C, R, T, aligned> lowest() noexcept
        {
            return amal::detail::splat_mat<C, R, T, aligned>(base::lowest());
        }
        static constexpr amal::mat<C, R, T, aligned> max() noexcept
        {
            return amal::detail::splat_mat<C, R, T, aligned>(base::max());
        }
        static constexpr amal::mat<C, R, T, aligned> epsilon() noexcept
        {
            return amal::detail::splat_mat<C, R, T, aligned>(base::epsilon());
        }
        static constexpr amal::mat<C, R, T, aligned> round_error() noexcept
        {
            return amal::detail::splat_mat<C, R, T, aligned>(base::round_error());
        }
        static constexpr amal::mat<C, R, T, aligned> infinity() noexcept
        {
            if constexpr (base::has_infinity) { return amal::detail::splat_mat<C, R, T, aligned>(base::infinity()); }
            else { return amal::detail::splat_mat<C, R, T, aligned>(T{}); }
        }
        static constexpr amal::mat<C, R, T, aligned> quiet_NaN() noexcept
        {
            if constexpr (base::has_quiet_NaN) { return amal::detail::splat_mat<C, R, T, aligned>(base::quiet_NaN()); }
            else { return amal::detail::splat_mat<C, R, T, aligned>(T{}); }
        }
        static constexpr amal::mat<C, R, T, aligned> signaling_NaN() noexcept
        {
            if constexpr (base::has_signaling_NaN)
            {
                return amal::detail::splat_mat<C, R, T, aligned>(base::signaling_NaN());
            }
            else { return amal::detail::splat_mat<C, R, T, aligned>(T{}); }
        }
        static constexpr amal::mat<C, R, T, aligned> denorm_min() noexcept
        {
            return amal::detail::splat_mat<C, R, T, aligned>(base::denorm_min());
        }
    };

} // namespace std
#endif // AMAL_STD_NUMERIC_LIMITS_MAT_DEFINED
