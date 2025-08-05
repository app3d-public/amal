#pragma once

#include "../type_info.hpp"

#if defined(AMAL_FORCE_HALF_PRECISION)
    #include "half.hpp"
#endif

namespace amal
{
    template <length_t C, length_t R, typename T, bool aligned = false>
    struct mat;

    using mat2x2_aligned = mat<2, 2, AMAL_FLOAT_TYPE, true>;
    using mat2x2_packed = mat<2, 2, AMAL_FLOAT_TYPE, false>;
    using dmat2x2_aligned = mat<2, 2, double, true>;
    using dmat2x2_packed = mat<2, 2, double, false>;
    using imat2x2_aligned = mat<2, 2, int, true>;
    using imat2x2_packed = mat<2, 2, int, false>;

    using mat2x3_aligned = mat<2, 3, AMAL_FLOAT_TYPE, true>;
    using mat2x3_packed = mat<2, 3, AMAL_FLOAT_TYPE, false>;
    using dmat2x3_aligned = mat<2, 3, double, true>;
    using dmat2x3_packed = mat<2, 3, double, false>;
    using imat2x3_aligned = mat<2, 3, int, true>;
    using imat2x3_packed = mat<2, 3, int, false>;

    using mat2x4 = mat<2, 4, AMAL_FLOAT_TYPE, true>;
    using dmat2x4 = mat<2, 4, double, true>;
    using imat2x4 = mat<2, 4, int, true>;

    using mat3x2_aligned = mat<3, 2, AMAL_FLOAT_TYPE, true>;
    using mat3x2_packed = mat<3, 2, AMAL_FLOAT_TYPE, false>;
    using dmat3x2_aligned = mat<3, 2, double, true>;
    using dmat3x2_packed = mat<3, 2, double, false>;
    using imat3x2_aligned = mat<3, 2, int, true>;
    using imat3x2_packed = mat<3, 2, int, false>;

    using mat3x3_aligned = mat<3, 3, AMAL_FLOAT_TYPE, true>;
    using mat3x3_packed = mat<3, 3, AMAL_FLOAT_TYPE, false>;
    using dmat3x3_aligned = mat<3, 3, double, true>;
    using dmat3x3_packed = mat<3, 3, double, false>;
    using imat3x3_aligned = mat<3, 3, int, true>;
    using imat3x3_packed = mat<3, 3, int, false>;

    using mat3x4 = mat<3, 4, AMAL_FLOAT_TYPE, true>;
    using dmat3x4 = mat<3, 4, double, true>;
    using imat3x4 = mat<3, 4, int, true>;

    using mat4x2_aligned = mat<4, 2, AMAL_FLOAT_TYPE, true>;
    using mat4x2_packed = mat<4, 2, AMAL_FLOAT_TYPE, false>;
    using dmat4x2_aligned = mat<4, 2, double, true>;
    using dmat4x2_packed = mat<4, 2, double, false>;
    using imat4x2_aligned = mat<4, 2, int, true>;
    using imat4x2_packed = mat<4, 2, int, false>;

    using mat4x3_aligned = mat<4, 3, AMAL_FLOAT_TYPE, true>;
    using mat4x3_packed = mat<4, 3, AMAL_FLOAT_TYPE, false>;
    using dmat4x3_aligned = mat<4, 3, double, true>;
    using dmat4x3_packed = mat<4, 3, double, false>;
    using imat4x3_aligned = mat<4, 3, int, true>;
    using imat4x3_packed = mat<4, 3, int, false>;

    using mat4x4 = mat<4, 4, AMAL_FLOAT_TYPE, true>;
    using dmat4x4 = mat<4, 4, double, true>;
    using imat4x4 = mat<4, 4, int, true>;

#ifdef AMAL_FORCE_ALIGNED_TYPES
    using mat2x2 = mat2x2_aligned;
    using mat2x3 = mat2x3_aligned;
    using mat3x2 = mat3x2_aligned;
    using mat3x3 = mat3x3_aligned;
    using mat3x4 = mat3x4_aligned;
    using mat4x2 = mat4x2_aligned;
    using mat4x3 = mat4x3_aligned;
#else
    using mat2x2 = mat2x2_packed;
    using mat2x3 = mat2x3_packed;
    using mat3x2 = mat3x2_packed;
    using mat3x3 = mat3x3_packed;
    using mat4x2 = mat4x2_packed;
    using mat4x3 = mat4x3_packed;
#endif

    using mat2 = mat2x2;
    using mat3 = mat3x3;
    using mat4 = mat4x4;

    template <typename T>
    struct is_matrix : std::false_type
    {
    };

    template <length_t C, length_t R, typename T, bool aligned>
    struct is_matrix<mat<C, R, T, aligned>> : std::true_type
    {
    };

    template <class T>
    inline constexpr bool is_matrix_v = is_matrix<T>::value;

#define AMAL_MAT(C, R, T, aligned) mat<C, R, T, aligned>
#define AMAL_MAT_SELF              AMAL_MAT(C, R, T, aligned)
#define AMAL_NMAT(C, R)            AMAL_MAT(C, R, T, aligned)
#define AMAL_NMAT_VAL_SIMD(C, R)   AMAL_TYPE_SIMD(AMAL_NMAT(C, R), AMAL_NMAT(C, R))
#define AMAL_NMAT_VAL_NOSIMD(C, R) AMAL_TYPE_NOSIMD(AMAL_NMAT(C, R), AMAL_NMAT(C, R))

    namespace internal
    {
        template <length_t C1, length_t R1, length_t C2, length_t R2, typename T, bool aligned>
        inline constexpr bool is_matrix_multiply_simdable =
            is_simd_enabled_v<typename mat<C1, R1, T, aligned>::simd_type::value_type> &&
            is_simd_enabled_v<typename mat<C2, R2, T, aligned>::simd_type::value_type> &&
            is_simd_enabled_v<typename mat<C2, R1, T, aligned>::simd_type::value_type>;
    }

#define AMAL_MAT_MUL_SIMD \
    std::enable_if_t<internal::is_matrix_multiply_simdable<C1, R1, C2, R2, T, aligned>, AMAL_NMAT(C2, R1)>
#define AMAL_NMAT_MUL_NOSIMD(C1, R1, C2, R2) \
    std::enable_if_t<!internal::is_matrix_multiply_simdable<C1, R1, C2, R2, T, aligned>, AMAL_NMAT(C2, R1)>
} // namespace amal