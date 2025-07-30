#pragma once

#include "details/half.hpp"

namespace amal
{
    struct half
    {
        details::uint16 data;

        constexpr half() : data(0) {};

        constexpr explicit half(details::uint16 rhs) : data(rhs) {}

        template <
            class T,
            std::enable_if_t<std::is_arithmetic_v<T> && !std::is_same_v<std::remove_cv_t<T>, details::uint16>, int> = 0>
        constexpr half(T v) noexcept : data(static_cast<details::uint16>(details::float_to_half(static_cast<float>(v))))
        {
        }

        constexpr operator float() const { return details::half_to_float(data); }

        constexpr half &operator=(float rhs)
        {
            data = static_cast<details::uint16>(details::float_to_half(rhs));
            return *this;
        }

        constexpr half &operator+=(half rhs) { return *this = *this + rhs; }
        constexpr half &operator-=(half rhs) { return *this = *this - rhs; }
        constexpr half &operator*=(half rhs) { return *this = *this * rhs; }
        constexpr half &operator/=(half rhs) { return *this = *this / rhs; }
        constexpr half &operator+=(float rhs) { return *this = *this + rhs; }
        constexpr half &operator-=(float rhs) { return *this = *this - rhs; }
        constexpr half &operator*=(float rhs) { return *this = *this * rhs; }
        constexpr half &operator/=(float rhs) { return *this = *this / rhs; }
        half &operator++() { return *this = *this + half(static_cast<details::uint16>(0x3C00)); }
        half &operator--() { return *this = *this + half(static_cast<details::uint16>(0xBC00)); }

        half operator++(int)
        {
            half out(*this);
            ++*this;
            return out;
        }

        half operator--(int)
        {
            half out(*this);
            --*this;
            return out;
        }
    };

    inline constexpr bool operator==(half x, half y)
    {
        return !details::compsignal(x.data, y.data) && (x.data == y.data || !((x.data | y.data) & 0x7FFF));
    }

    inline constexpr bool operator!=(half x, half y)
    {
        return details::compsignal(x.data, y.data) || (x.data != y.data && ((x.data | y.data) & 0x7FFF));
    }

    inline constexpr bool operator<(half x, half y)
    {
        return !details::compsignal(x.data, y.data) &&
               ((x.data ^ (0x8000 | (0x8000 - (x.data >> 15)))) + (x.data >> 15)) <
                   ((y.data ^ (0x8000 | (0x8000 - (y.data >> 15)))) + (y.data >> 15));
    }

    inline constexpr bool operator>(half x, half y)
    {
        return !details::compsignal(x.data, y.data) &&
               ((x.data ^ (0x8000 | (0x8000 - (x.data >> 15)))) + (x.data >> 15)) >
                   ((y.data ^ (0x8000 | (0x8000 - (y.data >> 15)))) + (y.data >> 15));
    }

    inline constexpr bool operator<=(half x, half y)
    {
        return !details::compsignal(x.data, y.data) &&
               ((x.data ^ (0x8000 | (0x8000 - (x.data >> 15)))) + (x.data >> 15)) <=
                   ((y.data ^ (0x8000 | (0x8000 - (y.data >> 15)))) + (y.data >> 15));
    }

    inline constexpr bool operator>=(half x, half y)
    {
        return !details::compsignal(x.data, y.data) &&
               ((x.data ^ (0x8000 | (0x8000 - (x.data >> 15)))) + (x.data >> 15)) >=
                   ((y.data ^ (0x8000 | (0x8000 - (y.data >> 15)))) + (y.data >> 15));
    }

    inline constexpr half operator+(half arg) { return arg; }
    inline constexpr half operator-(half arg) { return half(static_cast<details::uint16>(arg.data ^ 0x8000)); }
    inline constexpr half operator+(half x, half y) { return half(static_cast<float>(x) + static_cast<float>(y)); }
    inline constexpr half operator-(half x, half y) { return half(static_cast<float>(x) - static_cast<float>(y)); }
    inline constexpr half operator*(half x, half y) { return half(static_cast<float>(x) * static_cast<float>(y)); }
    inline constexpr half operator/(half x, half y) { return half(static_cast<float>(x) / static_cast<float>(y)); }

    inline constexpr half fabs(half arg) { return half(static_cast<details::uint16>(arg.data & 0x7FFF)); }
    inline constexpr half abs(half arg) { return fabs(arg); }

    inline half fmod(half x, half y)
    {
        unsigned int absx = x.data & 0x7FFF, absy = y.data & 0x7FFF, sign = x.data & 0x8000;
        if (absx >= 0x7C00 || absy >= 0x7C00)
            return half(static_cast<details::uint16>((absx > 0x7C00 || absy > 0x7C00) ? details::signal(x.data, y.data)
                                                     : (absx == 0x7C00)               ? AMAL_HALF_INVALID
                                                                                      : x.data));
        if (!absy) return half(static_cast<details::uint16>(AMAL_HALF_INVALID));
        if (!absx) return x;
        if (absx == absy) return half(static_cast<details::uint16>(sign));
        return half(static_cast<details::uint16>(sign | details::mod<false, false>(absx, absy)));
    }

    inline constexpr half exp(half arg) { return half(std::exp(static_cast<float>(arg))); }
    inline constexpr half exp2(half arg) { return half(std::exp2(static_cast<float>(arg))); }
    inline constexpr half expm1(half arg) { return half(std::expm1(static_cast<float>(arg))); }
    inline constexpr half log(half arg) { return half(std::log(static_cast<float>(arg))); }
    inline constexpr half log2(half arg) { return half(std::log2(static_cast<float>(arg))); }
    inline constexpr half log10(half arg) { return half(std::log10(static_cast<float>(arg))); }
    inline constexpr half log1p(half arg) { return half(std::log1p(static_cast<float>(arg))); }
    inline constexpr half sqrt(half arg) { return half(std::sqrt(static_cast<float>(arg))); }
    inline constexpr half cbrt(half arg) { return half(std::cbrt(static_cast<float>(arg))); }
    inline constexpr half pow(half x, half y) { return half(std::pow(static_cast<float>(x), static_cast<float>(y))); }
    inline constexpr half sin(half arg) { return half(std::sin(static_cast<float>(arg))); }
    inline constexpr half cos(half arg) { return half(std::cos(static_cast<float>(arg))); }
    inline constexpr half tan(half arg) { return half(std::tan(static_cast<float>(arg))); }
    inline constexpr half asin(half arg) { return half(std::asin(static_cast<float>(arg))); }
    inline constexpr half acos(half arg) { return half(std::acos(static_cast<float>(arg))); }
    inline constexpr half atan(half arg) { return half(std::atan(static_cast<float>(arg))); }
    inline constexpr half atan2(half y, half x)
    {
        return half(std::atan2(static_cast<float>(y), static_cast<float>(x)));
    }
    inline constexpr half sinh(half arg) { return half(std::sinh(static_cast<float>(arg))); }
    inline constexpr half cosh(half arg) { return half(std::cosh(static_cast<float>(arg))); }
    inline constexpr half tanh(half arg) { return half(std::tanh(static_cast<float>(arg))); }
    inline constexpr half asinh(half arg) { return half(std::asinh(static_cast<float>(arg))); }
    inline constexpr half acosh(half arg) { return half(std::acosh(static_cast<float>(arg))); }
    inline constexpr half atanh(half arg) { return half(std::atanh(static_cast<float>(arg))); }
    inline constexpr half floor(half arg) { return half(std::floor(static_cast<float>(arg))); }
    inline constexpr half ceil(half arg) { return half(std::ceil(static_cast<float>(arg))); }
    inline constexpr half trunc(half arg) { return half(std::trunc(static_cast<float>(arg))); }
    inline constexpr half round(half arg) { return half(std::round(static_cast<float>(arg))); }
    inline constexpr half lgamma(half arg) { return half(std::lgamma(static_cast<float>(arg))); }
    inline constexpr half tgamma(half arg) { return half(std::tgamma(static_cast<float>(arg))); }

    inline constexpr bool isfinite(half arg) { return (arg.data & 0x7C00) != 0x7C00; }
    inline constexpr bool isinf(half arg) { return (arg.data & 0x7FFF) == 0x7C00; }
    inline constexpr bool isnan(half arg) { return (arg.data & 0x7FFF) > 0x7C00; }
    inline constexpr bool isnormal(half arg) { return ((arg.data & 0x7C00) != 0) & ((arg.data & 0x7C00) != 0x7C00); }
    inline constexpr bool signbit(half arg) { return (arg.data & 0x8000) != 0; }
    inline constexpr bool isgreater(half x, half y)
    {
        return ((x.data ^ (0x8000 | (0x8000 - (x.data >> 15)))) + (x.data >> 15)) >
                   ((y.data ^ (0x8000 | (0x8000 - (y.data >> 15)))) + (y.data >> 15)) &&
               !isnan(x) && !isnan(y);
    }
    inline constexpr bool isgreaterequal(half x, half y)
    {
        return ((x.data ^ (0x8000 | (0x8000 - (x.data >> 15)))) + (x.data >> 15)) >=
                   ((y.data ^ (0x8000 | (0x8000 - (y.data >> 15)))) + (y.data >> 15)) &&
               !isnan(x) && !isnan(y);
    }
    inline constexpr bool isless(half x, half y)
    {
        return ((x.data ^ (0x8000 | (0x8000 - (x.data >> 15)))) + (x.data >> 15)) <
                   ((y.data ^ (0x8000 | (0x8000 - (y.data >> 15)))) + (y.data >> 15)) &&
               !isnan(x) && !isnan(y);
    }
    inline constexpr bool islessequal(half x, half y)
    {
        return ((x.data ^ (0x8000 | (0x8000 - (x.data >> 15)))) + (x.data >> 15)) <=
                   ((y.data ^ (0x8000 | (0x8000 - (y.data >> 15)))) + (y.data >> 15)) &&
               !isnan(x) && !isnan(y);
    }
    inline constexpr bool islessgreater(half x, half y)
    {
        return x.data != y.data && ((x.data | y.data) & 0x7FFF) && !isnan(x) && !isnan(y);
    }
    inline constexpr bool isunordered(half x, half y) { return isnan(x) || isnan(y); }

    inline constexpr half fmax(half x, half y)
    {
        return (!isnan(y) && (isnan(x) || (x.data ^ (0x8000 | (0x8000 - (x.data >> 15)))) <
                                              (y.data ^ (0x8000 | (0x8000 - (y.data >> 15))))))
                   ? y
                   : x;
    }

    inline constexpr half fmin(half x, half y)
    {
        return (!isnan(y) && (isnan(x) || (x.data ^ (0x8000 | (0x8000 - (x.data >> 15)))) >
                                              (y.data ^ (0x8000 | (0x8000 - (y.data >> 15))))))
                   ? y
                   : x;
    }
} // namespace amal

namespace std
{
    template <>
    class numeric_limits<amal::half>
    {
    public:
        static constexpr bool is_specialized = true;
        static constexpr bool is_signed = true;
        static constexpr bool is_integer = false;
        static constexpr bool is_exact = false;
        static constexpr bool is_modulo = false;
        static constexpr bool is_bounded = true;
        static constexpr bool is_iec559 = true;
        static constexpr bool has_infinity = true;
        static constexpr bool has_quiet_NaN = true;
        static constexpr bool has_signaling_NaN = true;
        static constexpr bool has_denorm_loss = false;
        static constexpr bool traps = false;
        static constexpr bool tinyness_before = false;
        static constexpr float_round_style round_style = std::float_round_style::round_to_nearest;
        static constexpr int digits = 11;
        static constexpr int digits10 = 3;
        static constexpr int max_digits10 = 5;
        static constexpr int radix = 2;
        static constexpr int min_exponent = -13;
        static constexpr int min_exponent10 = -4;
        static constexpr int max_exponent = 16;
        static constexpr int max_exponent10 = 4;
        static constexpr amal::half min() noexcept { return amal::half(static_cast<amal::details::uint16>(0x0400)); }
        static constexpr amal::half lowest() noexcept { return amal::half(static_cast<amal::details::uint16>(0xFBFF)); }
        static constexpr amal::half max() noexcept { return amal::half(static_cast<amal::details::uint16>(0x7BFF)); }
        static constexpr amal::half epsilon() noexcept
        {
            return amal::half(static_cast<amal::details::uint16>(0x1400));
        }
        static constexpr amal::half round_error() noexcept
        {
            return amal::half(static_cast<amal::details::uint16>(0x3800));
        }
        static constexpr amal::half infinity() noexcept
        {
            return amal::half(static_cast<amal::details::uint16>(0x7C00));
        }
        static constexpr amal::half quiet_NaN() noexcept
        {
            return amal::half(static_cast<amal::details::uint16>(0x7FFF));
        }
        static constexpr amal::half signaling_NaN() noexcept
        {
            return amal::half(static_cast<amal::details::uint16>(0x7DFF));
        }
        static constexpr amal::half denorm_min() noexcept
        {
            return amal::half(static_cast<amal::details::uint16>(0x0001));
        }
    };
} // namespace std

typedef amal::half f16;