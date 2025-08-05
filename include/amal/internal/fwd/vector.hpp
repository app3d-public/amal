#pragma once

#include "../type_info.hpp"

#if defined(AMAL_FORCE_HALF_PRECISION)
    #include "half.hpp"
#endif

namespace amal
{
    template <length_t L, typename T, bool aligned>
    struct vec;
    using vec2_aligned = vec<2, AMAL_FLOAT_TYPE, true>;
    using vec2_packed = vec<2, AMAL_FLOAT_TYPE, false>;
    using dvec2_aligned = vec<2, double, true>;
    using dvec2_packed = vec<2, double, false>;
    using ivec2_aligned = vec<2, int, true>;
    using ivec2_packed = vec<2, int, false>;
    using vec3_aligned = vec<3, AMAL_FLOAT_TYPE, true>;
    using vec3_packed = vec<3, AMAL_FLOAT_TYPE, false>;
    using dvec3_aligned = vec<3, double, true>;
    using dvec3_packed = vec<3, double, false>;
    using ivec3_aligned = vec<3, int, true>;
    using ivec3_packed = vec<3, int, false>;

#ifdef AMAL_FORCE_ALIGNED_TYPES
    using vec2 = vec2_aligned;
    using dvec2 = dvec2_aligned;
    using ivec2 = ivec2_aligned;

    using vec3 = vec3_aligned;
    using dvec3 = dvec3_aligned;
    using ivec3 = ivec3_aligned;
#else
    using vec2 = vec2_packed;
    using dvec2 = dvec2_packed;
    using ivec2 = ivec2_packed;

    using vec3 = vec3_packed;
    using dvec3 = dvec3_packed;
    using ivec3 = ivec3_packed;
#endif

    using bvec2 = vec<2, bool, false>;
    using bvec3 = vec<3, bool, false>;
    using bvec4 = vec<4, bool, true>;

    using uvec2 = vec<2, uint, false>;
    using uvec3 = vec<3, uint, false>;
    using uvec4 = vec<4, uint, true>;

    using vec4 = vec<4, AMAL_FLOAT_TYPE, true>;
    using dvec4 = vec<4, double, true>;
    using ivec4 = vec<4, int, true>;

    template <typename T>
    struct is_vector : std::false_type
    {
    };

    template <length_t N, typename T, bool aligned>
    struct is_vector<vec<N, T, aligned>> : std::true_type
    {
    };

    template <class T>
    inline constexpr bool is_vector_v = is_vector<T>::value;

#define AMAL_VEC(N, T, aligned) vec<N, T, aligned>
#define AMAL_NVEC(N)            AMAL_VEC(N, T, aligned)
#define AMAL_VEC_SELF           AMAL_VEC(N, T, aligned)
#define AMAL_VEC_REF_SIMD       AMAL_TYPE_SIMD(AMAL_VEC_SELF, AMAL_VEC_SELF &)
#define AMAL_VEC_REF_NOSIMD     AMAL_TYPE_NOSIMD(AMAL_VEC_SELF, AMAL_VEC_SELF &)
#define AMAL_VEC_VAL_SIMD       AMAL_TYPE_SIMD(AMAL_VEC_SELF, AMAL_VEC_SELF)
#define AMAL_VEC_VAL_NOSIMD     AMAL_TYPE_NOSIMD(AMAL_VEC_SELF, AMAL_VEC_SELF)
#define AMAL_BVEC               AMAL_VEC(N, bool, aligned)
#define AMAL_BVEC_VAL_SIMD      AMAL_TYPE_SIMD(AMAL_BVEC, bool)
#define AMAL_BVEC_VAL_NOSIMD    AMAL_TYPE_NOSIMD(AMAL_BVEC, bool)
#define AMAL_IVEC               AMAL_VEC(N, int, aligned)
} // namespace amal