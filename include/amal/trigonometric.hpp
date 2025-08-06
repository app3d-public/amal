#pragma once

#include "vector.hpp"

namespace amal
{
    template <typename T>
    inline constexpr T pi()
    {
        static_assert(std::numeric_limits<T>::is_iec559, "pi only accepts floating-point inputs");
        return static_cast<T>(3.14159265358979323846264338327950288);
    }

    template <typename T>
    inline constexpr T half_pi()
    {
        return T(1.57079632679489661923132169163975144);
    }

    template <typename T>
    inline constexpr T quarter_pi()
    {
        return T(0.785398163397448309615660845819875721);
    }

    template <typename T>
    inline constexpr T e()
    {
        return T(2.71828182845904523536);
    }

    template <typename T>
    inline constexpr T radians(T degrees)
    {
        static_assert(std::numeric_limits<T>::is_iec559, "'radians' only accept floating-point input");
        return degrees * static_cast<T>(0.01745329251994329576923690768489);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF radians(AMAL_VEC_SELF const &degrees)
    {
        return internal::create_by_call(degrees, radians);
    }

    template <typename T>
    inline constexpr T degrees(T radians)
    {
        static_assert(std::numeric_limits<T>::is_iec559, "'degrees' only accept floating-point input");
        return degrees * static_cast<T>(57.295779513082320876798154814105);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF degrees(AMAL_VEC_SELF const &radians)
    {
        return internal::create_by_call(radians, degrees);
    }

    using std::sin;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF sin(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::sin);
    }

    using std::cos;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF cos(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::cos);
    }

    using std::tan;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF tan(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::tan);
    }

    using std::asin;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF asin(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::asin);
    }

    using std::acos;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF acos(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::acos);
    }

    using std::atan;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF atan(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::atan);
    }

    using std::atan2;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF atan2(AMAL_VEC_SELF const &y, AMAL_VEC_SELF const &x)
    {
        return internal::create_by_call(y, x, std::atan2<T>);
    }

    using std::sinh;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF sinh(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::sinh);
    }

    using std::cosh;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF cosh(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::cosh);
    }

    using std::tanh;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF tanh(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::tanh);
    }

    using std::asinh;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF asinh(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::asinh);
    }

    using std::acosh;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF acosh(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::acosh);
    }

    using std::atanh;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF atanh(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::atanh);
    }
} // namespace amal