#pragma once

#include <cstring>
#include <type_traits>
#if defined(__SSE2__)
    #include <emmintrin.h>
    #include <xmmintrin.h>

#endif
#if defined(__AVX__)
    #include <immintrin.h>
#endif

#if defined(AMAL_FORCE_HALF_PRECISION)
    #define AMAL_FLOAT_TYPE half
#else
    #define AMAL_FLOAT_TYPE float
#endif

#if defined(__FMA__) && !defined(AMAL_FMA_DISABLE)
    #define AMAL_FMA_ENABLE
#endif

#if defined(__SSE2__) && !defined(AMAL_SIMD_DISABLE)
    #define AMAL_SIMD_ENABLE
#endif

#if defined(AMAL_FMA_ENABLE) or defined(AMAL_SIMD_ENABLE)
    #define AMAL_CONSTEXPR
#else
    #define AMAL_CONSTEXPR constexpr
#endif

#define AMAL_CEIL(a, b)    (((a % b) == 0) ? a / b : (a / b) + 1)
#define AMAL_ALIGN_SIZE(a) (int(AMAL_CEIL(a, 4)) * 4)

typedef int __v4si_u __attribute__((__vector_size__(16), __aligned__(1)));

namespace amal
{
    using length_t = int;

    typedef struct half half;

    namespace internal
    {
        struct simd_disabled
        {
        };

        template <length_t L, typename T, bool aligned>
        struct simd_type
        {
            using value_type = simd_disabled;
            static inline constexpr const length_t max_size = L;
        };

#define AMAL_SIMD_TYPE_BUILD(L, T, V)                                   \
    template <>                                                         \
    struct simd_type<L, T, true>                                        \
    {                                                                   \
        using value_type = V;                                           \
        static inline constexpr length_t max_size = AMAL_ALIGN_SIZE(L); \
    };

#define AMAL_SIMD_TYPE_BUILD_GROUP(T, V) \
    AMAL_SIMD_TYPE_BUILD(2, T, V);       \
    AMAL_SIMD_TYPE_BUILD(3, T, V);       \
    AMAL_SIMD_TYPE_BUILD(4, T, V);

#if defined(__SSE2__)
        AMAL_SIMD_TYPE_BUILD_GROUP(float, __m128_u);
        AMAL_SIMD_TYPE_BUILD_GROUP(int, __v4si_u);
#endif

#if defined(__AVX__)
        AMAL_SIMD_TYPE_BUILD_GROUP(double, __m256d_u);
#endif

        template <typename V>
        inline constexpr bool is_simd_enabled_v = !std::is_same_v<V, simd_disabled>;

        template <class T>
        struct is_floating_point : public std::false_type
        {
        };

        template <>
        struct is_floating_point<float> : public std::true_type
        {
        };
        template <>
        struct is_floating_point<double> : public std::true_type
        {
        };
        template <>
        struct is_floating_point<long double> : public std::true_type
        {
        };
        template <>
        struct is_floating_point<half> : public std::true_type
        {
        };
    } // namespace internal

    template <class T>
    struct is_floating_point : public internal::is_floating_point<std::__remove_cv_t<T>>
    {
    };

    template <class T>
    inline constexpr bool is_floating_point_v = is_floating_point<T>::value;

    template <class T>
    struct is_arithmetic
        : public std::integral_constant<bool, std::is_integral<T>::value || is_floating_point<T>::value>
    {
    };

    template <class T>
    inline constexpr bool is_arithmetic_v = is_arithmetic<T>::value;

#define AMAL_TYPE_SIMD(C, ...) \
    std::enable_if_t<internal::is_simd_enabled_v<typename C::simd_type::value_type>, __VA_ARGS__>
#define AMAL_TYPE_NOSIMD(C, ...) \
    std::enable_if_t<!internal::is_simd_enabled_v<typename C::simd_type::value_type>, __VA_ARGS__>
#define AMAL_CONSTRUCT_SIMD \
    typename V = typename simd_type::value_type, std::enable_if_t<internal::is_simd_enabled_v<V>, int> = 0
#define AMAL_CONSTRUCT_NOSIMD \
    typename V = typename simd_type::value_type, std::enable_if_t<!internal::is_simd_enabled_v<V>, int> = 0

    typedef unsigned int uint;
} // namespace amal

#ifndef AMAL_NO_GLOBAL_ALIASES
using uint = amal::uint;
#endif