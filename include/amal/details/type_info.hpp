#pragma once

#include <type_traits>
#if defined(__SSE2__)
    #include <emmintrin.h>
    #include <xmmintrin.h>

#endif
#if defined(__AVX__)
    #include <immintrin.h>
#endif

#define AMAL_PRECISION_MEDIUM 0
#define AMAL_PRECISION_HIGH   1

#ifndef AMAL_PRECISION
    #define AMAL_PRECISION AMAL_PRECISION_HIGH
#endif

#define AMAL_CEIL(a, b)    (((a % b) == 0) ? a / b : (a / b) + 1)
#define AMAL_ALIGN_SIZE(a) (int(AMAL_CEIL(a, 4)) * 4)

namespace amal
{
    using length_t = int;

    enum class Pack
    {
        aligned,
        packed
    };

    typedef struct half half;

    namespace details
    {
        struct simd_disabled
        {
        };

        template <length_t L, typename T, Pack P>
        struct simd_type
        {
            using value_type = simd_disabled;
            static inline constexpr const length_t max_size = L;
        };

#ifndef _MSC_VER
    #define AMAL_SIMD_TYPE_BUILD(L, T, P, V)                                                                 \
        template <>                                                                                          \
        struct simd_type<L, T, P>                                                                            \
        {                                                                                                    \
            using value_type = V;                                                                            \
            static inline constexpr length_t max_size = (P != amal::Pack::aligned ? L : AMAL_ALIGN_SIZE(L)); \
        };

    #define AMAL_SIMD_TYPE_BUILD_GROUP(T, V)          \
        AMAL_SIMD_TYPE_BUILD(2, T, Pack::aligned, V); \
        AMAL_SIMD_TYPE_BUILD(3, T, Pack::aligned, V); \
        AMAL_SIMD_TYPE_BUILD(4, T, Pack::aligned, V);

    #if defined(__SSE2__)
        AMAL_SIMD_TYPE_BUILD_GROUP(float, __v4sf);
        AMAL_SIMD_TYPE_BUILD_GROUP(int, __v4si);
    #endif

    #if defined(__AVX__)
        AMAL_SIMD_TYPE_BUILD_GROUP(double, __v4df);
    #endif
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
    } // namespace details

    template <class T>
    struct is_floating_point : public details::is_floating_point<std::__remove_cv_t<T>>
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
    std::enable_if_t<details::is_simd_enabled_v<typename C::simd_type::value_type>, __VA_ARGS__>
#define AMAL_TYPE_NOSIMD(C, ...) \
    std::enable_if_t<!details::is_simd_enabled_v<typename C::simd_type::value_type>, __VA_ARGS__>
#define AMAL_CONSTRUCT_SIMD \
    typename V = typename simd_type::value_type, std::enable_if_t<details::is_simd_enabled_v<V>, int> = 0
#define AMAL_CONSTRUCT_NOSIMD \
    typename V = typename simd_type::value_type, std::enable_if_t<!details::is_simd_enabled_v<V>, int> = 0
} // namespace amal

typedef unsigned int uint;