#pragma once

#include <algorithm>
#ifdef AMAL_SIMD_ENABLE
    #include "internal/simd/common.hpp"
#endif
#include "compare.hpp"
#include "vector.hpp"

namespace amal
{
    using std::abs;
#ifdef AMAL_SIMD_ENABLE
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD abs(AMAL_VEC_SELF const &v)
    {
        return AMAL_VEC_SELF(internal::abs(v.s));
    }
#endif

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_VAL_NOSIMD abs(AMAL_NVEC(N) const &v)
    {
        if constexpr (is_floating_point_v<T>)
            return internal::create_by_call(v, std::fabs);
        else
            return internal::create_by_call(v, std::abs);
    }

#ifdef AMAL_SIMD_ENABLE
    template <length_t N, typename T, typename U, bool aligned>
    inline AMAL_VEC_VAL_SIMD mix(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y, AMAL_VEC(N, U, aligned) const &a)
    {
        return AMAL_VEC_SELF(internal::mix(x.s, y.s, a.s));
    }
#endif

    template <length_t N, typename T, typename U, bool aligned>
    inline constexpr AMAL_VEC_VAL_NOSIMD mix(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y,
                                             AMAL_VEC(N, U, aligned) const &a)
    {
        if constexpr (std::is_same_v<U, bool>)
        {
            AMAL_VEC_SELF r;
            for (length_t i = 0; i < x.length(); ++i) r[i] = a[i] ? y[i] : x[i];
            return r;
        }
        return AMAL_VEC_SELF(x * (static_cast<U>(1) - a) + y * a);
    }

    template <length_t N, typename T, typename U, bool aligned>
    inline constexpr std::enable_if_t<is_arithmetic_v<U>, AMAL_VEC_SELF> mix(AMAL_VEC_SELF const &x,
                                                                             AMAL_VEC_SELF const &y, U a)
    {
        if constexpr (std::is_same_v<U, bool>) return a ? y : x;
        return mix(x, y, AMAL_VEC_SELF(a));
    }

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF sign(AMAL_VEC_SELF const &x)
    {
        if constexpr (is_floating_point_v<T>)
            return AMAL_VEC_SELF(less_than(AMAL_VEC_SELF(0), x)) - AMAL_VEC_SELF(less_than(x, AMAL_VEC_SELF(0)));
        else
        {
            T const shift(static_cast<T>(sizeof(T) * 8 - 1));
            AMAL_VEC_SELF const y(vec<N, typename std::make_unsigned<T>::type, aligned>(-x) >>
                                  typename std::make_unsigned<T>::type(shift));
            return (x >> shift) | y;
        }
    }

    template <typename T>
    inline constexpr std::enable_if_t<is_floating_point_v<T>, T> sign(T value)
    {
        return (value > 0) ? T(1) : ((value < 0) ? T(-1) : T(0));
    }

    using std::floor;
#ifdef AMAL_SIMD_ENABLE
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD floor(AMAL_VEC_SELF const &v)
    {
        if constexpr (is_floating_point_v<T>)
            return AMAL_VEC_SELF(internal::floor(v.s));
        else
            return AMAL_VEC_SELF(v.s);
    }
#endif

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD floor(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::floor);
    }

    using std::round;
#ifdef AMAL_SIMD_ENABLE
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD round(AMAL_VEC_SELF const &v)
    {
        if constexpr (is_floating_point_v<T>)
            return AMAL_VEC_SELF(internal::round(v.s));
        else
            return AMAL_VEC_SELF(v.s);
    }
#endif

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD round(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::round);
    }
    // Rounds a floating-point number to the nearest power of 10.
    template <typename T>
    constexpr T round10(T value)
    {
        if (value == 0) return 0;
        float value_abs = abs(value);
        int exponent = static_cast<int>(round(log10(value_abs)));
        float base = pow(10, exponent);
        return base * sign(value);
    }

    using std::ceil;
#ifdef AMAL_SIMD_ENABLE
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD ceil(AMAL_VEC_SELF const &v)
    {
        if constexpr (is_floating_point_v<T>)
            return AMAL_VEC_SELF(internal::ceil(v.s));
        else
            return AMAL_VEC_SELF(v.s);
    }
#endif

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD ceil(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::ceil);
    }

    using std::trunc;
#ifdef AMAL_SIMD_ENABLE
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD trunc(AMAL_VEC_SELF const &v)
    {
        if constexpr (is_floating_point_v<T>)
            return AMAL_VEC_SELF(internal::trunc(v.s));
        else
            return AMAL_VEC_SELF(v.s);
    }
#endif

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_NOSIMD trunc(AMAL_VEC_SELF const &v)
    {
        return internal::create_by_call(v, std::trunc);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF fract(AMAL_VEC_SELF const &v)
    {
        return v - floor(v);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF mod(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
#if defined(__FMA__)
        return fma(-y, floor(x / y), x);
#else
        return x - y * floor(x / y);
#endif
    }

    template <length_t N, typename T, typename U, bool aligned>
    inline std::enable_if_t<is_arithmetic_v<U>, AMAL_VEC_SELF> mod(AMAL_VEC_SELF const &x, U scalar)
    {
        return mod(x, AMAL_VEC_SELF(scalar));
    }

    using std::modf;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF modf(AMAL_VEC_SELF const &x, AMAL_VEC_SELF &i)
    {
        return internal::create_by_call(x, i, std::modf);
    }

    using std::fma;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF fma(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y, AMAL_VEC_SELF const &z)
    {
        return x * y + z;
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF splat_x(AMAL_VEC_SELF const &v)
    {
        return AMAL_VEC_SELF(v.x);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF splat_y(AMAL_VEC_SELF const &v)
    {
        return AMAL_VEC_SELF(v.y);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF splat_z(AMAL_VEC_SELF const &v)
    {
        return AMAL_VEC_SELF(v.z);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF splat_w(AMAL_VEC_SELF const &v)
    {
        return AMAL_VEC_SELF(v.w);
    }

    using std::min;
#ifdef AMAL_SIMD_ENABLE
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD min(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        return AMAL_VEC_SELF(internal::min(x.s, y.s));
    }
#endif

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_VAL_NOSIMD min(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        using PFN_min = const T &(*)(const T &, const T &);
        return internal::create_by_call(x, y, (PFN_min)std::min<T>);
    }

    using std::max;
#ifdef AMAL_SIMD_ENABLE
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_VAL_SIMD max(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        return AMAL_VEC_SELF(internal::max(x.s, y.s));
    }
#endif

    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_VAL_NOSIMD max(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        using PFN_max = const T &(*)(const T &, const T &);
        return internal::create_by_call(x, y, (PFN_max)std::max<T>);
    }

    using std::minmax;
    template <length_t N, typename T, bool aligned>
    inline std::pair<AMAL_VEC_SELF, AMAL_VEC_SELF> minmax(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &y)
    {
        return std::make_pair(min(x, y), max(x, y));
    }

    using std::clamp;
    template <length_t N, typename T, bool aligned>
    inline constexpr AMAL_VEC_SELF clamp(AMAL_VEC_SELF const &x, AMAL_VEC_SELF const &min_val,
                                         AMAL_VEC_SELF const &max_val)
    {
        return min(max(x, min_val), max_val);
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline std::enable_if_t<is_arithmetic_v<U>, AMAL_VEC_SELF> clamp(AMAL_VEC_SELF const &x, U min_val, U max_val)
    {
        return clamp(x, AMAL_VEC_SELF(min_val), AMAL_VEC_SELF(max_val));
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF step(AMAL_VEC_SELF const &edge, AMAL_VEC_SELF const &x)
    {
        return mix(AMAL_VEC_SELF(1), AMAL_VEC_SELF(0), less_than(x, edge));
    }

    template <length_t N, typename T, bool aligned, typename U>
    inline std::enable_if_t<is_arithmetic_v<U>, AMAL_VEC_SELF> step(U edge, AMAL_VEC_SELF const &x)
    {
        return step(AMAL_VEC_SELF(edge), x);
    }

    template <typename T>
    inline T smoothstep(T edge0, T edge1, T x)
    {
        T const tmp(std::clamp((x - edge0) / (edge1 - edge0), static_cast<T>(0), static_cast<T>(1)));
        return tmp * tmp * (static_cast<T>(3) - static_cast<T>(2) * tmp);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF smoothstep(AMAL_VEC_SELF const &edge0, AMAL_VEC_SELF const &edge1, AMAL_VEC_SELF const &x)
    {
        AMAL_VEC_SELF const tmp(clamp((x - edge0) / (edge1 - edge0), static_cast<T>(0), static_cast<T>(1)));
        return tmp * tmp * (static_cast<T>(3) - static_cast<T>(2) * tmp);
    }

    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF smoothstep(T edge0, T edge1, AMAL_VEC_SELF const &x)
    {
        return smoothstep(AMAL_VEC_SELF(edge0), AMAL_VEC_SELF(edge1), x);
    }

    using std::isnan;
    template <length_t N, typename T, bool aligned>
    inline AMAL_BVEC isnan(AMAL_VEC_SELF const &v)
    {
        AMAL_BVEC r(0);
        for (length_t l = 0; l < v.length(); ++l) r[l] = std::isnan(v[l]);
        return r;
    }

    using std::isinf;
    template <length_t N, typename T, bool aligned>
    inline AMAL_BVEC isinf(AMAL_VEC_SELF const &v)
    {
        AMAL_BVEC r(0);
        for (length_t l = 0; l < v.length(); ++l) r[l] = std::isinf(v[l]);
        return r;
    }

    using std::frexp;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF frexp(AMAL_VEC_SELF const &v, AMAL_IVEC &exp)
    {
        AMAL_VEC_SELF r(0);
        for (length_t l = 0; l < v.length(); ++l) r[l] = std::frexp(v[l], &exp[l]);
        return r;
    }

    using std::ldexp;
    template <length_t N, typename T, bool aligned>
    inline AMAL_VEC_SELF ldexp(AMAL_VEC_SELF const &v, AMAL_IVEC const &exp)
    {
        AMAL_VEC_SELF r(0);
        for (length_t l = 0; l < v.length(); ++l) r[l] = std::ldexp(v[l], exp[l]);
        return r;
    }
} // namespace amal