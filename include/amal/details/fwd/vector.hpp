#pragma once

#include "../type_info.hpp"

#if AMAL_PRECISION == AMAL_PRECISION_MEDIUM
    #include "half.hpp"
#endif

namespace amal
{
    template <length_t L, typename T, Pack P>
    struct vec;

#if AMAL_PRECISION == AMAL_PRECISION_MEDIUM
    #define AMAL_VEC_FLOAT_TYPE half
#else
    #define AMAL_VEC_FLOAT_TYPE float
#endif
    using vec2_aligned = vec<2, AMAL_VEC_FLOAT_TYPE, Pack::aligned>;
    using vec2_packed = vec<2, AMAL_VEC_FLOAT_TYPE, Pack::packed>;
    using dvec2_aligned = vec<2, double, Pack::aligned>;
    using dvec2_packed = vec<2, double, Pack::packed>;
    using ivec2_aligned = vec<2, int, Pack::aligned>;
    using ivec2_packed = vec<2, int, Pack::packed>;
    using vec3_aligned = vec<3, AMAL_VEC_FLOAT_TYPE, Pack::aligned>;
    using vec3_packed = vec<3, AMAL_VEC_FLOAT_TYPE, Pack::packed>;
    using dvec3_aligned = vec<3, double, Pack::aligned>;
    using dvec3_packed = vec<3, double, Pack::packed>;
    using ivec3_aligned = vec<3, int, Pack::aligned>;
    using ivec3_packed = vec<3, int, Pack::packed>;

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

    using bvec2 = vec<2, bool, Pack::packed>;
    using bvec3 = vec<3, bool, Pack::packed>;
    using bvec4 = vec<4, bool, Pack::aligned>;

    using uvec2 = vec<2, uint, Pack::packed>;
    using uvec3 = vec<3, uint, Pack::packed>;
    using uvec4 = vec<4, uint, Pack::aligned>;

    using vec4 = vec<4, AMAL_VEC_FLOAT_TYPE, Pack::aligned>;
    using dvec4 = vec<4, double, Pack::aligned>;
    using ivec4 = vec<4, int, Pack::aligned>;

    template <typename T>
    struct is_vector : std::false_type
    {
    };

    template <length_t N, typename T, Pack P>
    struct is_vector<vec<N, T, P>> : std::true_type
    {
    };

    template <class T>
    inline constexpr bool is_vector_v = is_vector<T>::value;

#define AMAL_VEC(N, T, P)    vec<N, T, P>
#define AMAL_NVEC(N)         AMAL_VEC(N, T, P)
#define AMAL_VEC_SELF        AMAL_VEC(N, T, P)
#define AMAL_VEC_REF_SIMD    AMAL_TYPE_SIMD(AMAL_VEC_SELF, AMAL_VEC_SELF &)
#define AMAL_VEC_REF_NOSIMD  AMAL_TYPE_NOSIMD(AMAL_VEC_SELF, AMAL_VEC_SELF &)
#define AMAL_VEC_VAL_SIMD    AMAL_TYPE_SIMD(AMAL_VEC_SELF, AMAL_VEC_SELF)
#define AMAL_VEC_VAL_NOSIMD  AMAL_TYPE_NOSIMD(AMAL_VEC_SELF, AMAL_VEC_SELF)
#define AMAL_BVEC            AMAL_VEC(N, bool, P)
#define AMAL_BVEC_VAL_SIMD   AMAL_TYPE_SIMD(AMAL_BVEC, bool)
#define AMAL_BVEC_VAL_NOSIMD AMAL_TYPE_NOSIMD(AMAL_BVEC, bool)
#define AMAL_IVEC            AMAL_VEC(N, int, P)
} // namespace amal