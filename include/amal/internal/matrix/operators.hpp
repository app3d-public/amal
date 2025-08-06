#pragma once

#include "../fwd/matrix.hpp"
#include "../type_info.hpp"
#include "functions.hpp"

#define AMAL_MATRIX_TEMP_DECL length_t C, length_t R, typename T, bool aligned

namespace amal
{
    template <AMAL_MATRIX_TEMP_DECL, typename U>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF &operator+=(AMAL_MAT_SELF &self, U scalar)
    {
        self.data[0] += scalar;
        self.data[1] += scalar;
        if constexpr (R > 2) self.data[2] += scalar;
        if constexpr (R > 3) self.data[3] += scalar;
        return self;
    }

    template <AMAL_MATRIX_TEMP_DECL, typename U>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF &operator+=(AMAL_MAT_SELF &self, AMAL_MAT(C, R, U, aligned) const &m)
    {
        self.data[0] += m.data[0];
        self.data[1] += m.data[1];
        if constexpr (R > 2) self.data[2] += m.data[2];
        if constexpr (R > 3) self.data[3] += m.data[3];
        return self;
    }

    template <AMAL_MATRIX_TEMP_DECL, typename U>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF &operator-=(AMAL_MAT_SELF &self, U scalar)
    {
        self.data[0] -= scalar;
        self.data[1] -= scalar;
        if constexpr (R > 2) self.data[2] -= scalar;
        if constexpr (R > 3) self.data[3] -= scalar;
        return self;
    }

    template <AMAL_MATRIX_TEMP_DECL, typename U>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF &operator-=(AMAL_MAT_SELF &self, AMAL_MAT(C, R, U, aligned) const &m)
    {
        self.data[0] -= m.data[0];
        self.data[1] -= m.data[1];
        if constexpr (R > 2) self.data[2] -= m.data[2];
        if constexpr (R > 3) self.data[3] -= m.data[3];
        return self;
    }

    template <AMAL_MATRIX_TEMP_DECL, typename U>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF &operator*=(AMAL_MAT_SELF &self, U scalar)
    {
        self.data[0] *= scalar;
        self.data[1] *= scalar;
        if constexpr (R > 2) self.data[2] *= scalar;
        if constexpr (R > 3) self.data[3] *= scalar;
        return self;
    }

    template <AMAL_MATRIX_TEMP_DECL, typename U>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF &operator*=(AMAL_MAT_SELF &self, AMAL_MAT(C, R, U, aligned) const &m)
    {
        return self = self * m;
    }

    template <AMAL_MATRIX_TEMP_DECL, typename U>
    AMAL_CONSTEXPR AMAL_MAT_SELF &operator/=(AMAL_MAT_SELF &self, U scalar)
    {
        self.data[0] /= scalar;
        self.data[1] /= scalar;
        if constexpr (R > 2) self.data[2] /= scalar;
        if constexpr (R > 3) self.data[3] /= scalar;
        return self;
    }

    template <AMAL_MATRIX_TEMP_DECL, typename U>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF &operator/=(AMAL_MAT_SELF &self, AMAL_MAT(C, R, U, aligned) const &m)
    {
        return self = self * inverse_matrix(m);
    }

    template <AMAL_MATRIX_TEMP_DECL>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF &operator++(AMAL_MAT_SELF &self)
    {
        ++self.data[0];
        ++self.data[1];
        if constexpr (R > 2) ++self.data[2];
        if constexpr (R > 3) ++self.data[3];
        return self;
    }

    template <AMAL_MATRIX_TEMP_DECL>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF &operator--(AMAL_MAT_SELF &self)
    {
        --self.data[0];
        --self.data[1];
        if constexpr (R > 2) --self.data[2];
        if constexpr (R > 3) --self.data[3];
        return self;
    }

    template <AMAL_MATRIX_TEMP_DECL>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF operator++(AMAL_MAT_SELF const &self, int)
    {
        AMAL_MAT_SELF ret(self);
        ++self;
        return ret;
    }

    template <AMAL_MATRIX_TEMP_DECL>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF operator--(AMAL_MAT_SELF const &self, int)
    {
        AMAL_MAT_SELF ret(self);
        --self;
        return ret;
    }

    template <AMAL_MATRIX_TEMP_DECL>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF operator+(AMAL_MAT_SELF const &self)
    {
        return self;
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(2, R) operator-(AMAL_NMAT(2, R) const &self)
    {
        return AMAL_NMAT(2, R)(-self.data[0], -self.data[1]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(3, R) operator-(AMAL_NMAT(3, R) const &self)
    {
        return AMAL_NMAT(3, R)(-self.data[0], -self.data[1], -self.data[2]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(4, R) operator-(AMAL_NMAT(4, R) const &self)
    {
        return AMAL_NMAT(4, R)(-self.data[0], -self.data[1], -self.data[2], -self.data[3]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(2, R) operator+(AMAL_NMAT(2, R) const &self, T scalar)
    {
        return AMAL_NMAT(2, R)(self.data[0] + scalar, self.data[1] + scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(3, R) operator+(AMAL_NMAT(3, R) const &self, T scalar)
    {
        return AMAL_NMAT(3, R)(self.data[0] + scalar, self.data[1] + scalar, self.data[2] + scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(4, R) operator+(AMAL_NMAT(4, R) const &self, T scalar)
    {
        return AMAL_NMAT(4, R)(self.data[0] + scalar, self.data[1] + scalar, self.data[2] + scalar,
                               self.data[3] + scalar);
    }

    template <AMAL_MATRIX_TEMP_DECL>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF operator+(T scalar, AMAL_MAT_SELF const &self)
    {
        return self + scalar;
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(2, R) operator-(AMAL_NMAT(2, R) const &self, T scalar)
    {
        return AMAL_NMAT(2, R)(self.data[0] - scalar, self.data[1] - scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(3, R) operator-(AMAL_NMAT(3, R) const &self, T scalar)
    {
        return AMAL_NMAT(3, R)(self.data[0] - scalar, self.data[1] - scalar, self.data[2] - scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(4, R) operator-(AMAL_NMAT(4, R) const &self, T scalar)
    {
        return AMAL_NMAT(4, R)(self.data[0] - scalar, self.data[1] - scalar, self.data[2] - scalar,
                               self.data[3] - scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(2, R) operator-(T scalar, AMAL_NMAT(2, R) const &self)
    {
        return AMAL_NMAT(2, R)(scalar - self.data[0], scalar - self.data[1]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(3, R) operator-(T scalar, AMAL_NMAT(3, R) const &self)
    {
        return AMAL_NMAT(3, R)(scalar - self.data[0], scalar - self.data[1], scalar - self.data[2]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(4, R) operator-(T scalar, AMAL_NMAT(4, R) const &self)
    {
        return AMAL_NMAT(4, R)(scalar - self.data[0], scalar - self.data[1], scalar - self.data[2],
                               scalar - self.data[3]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(2, R) operator+(AMAL_NMAT(2, R) const &m1, AMAL_NMAT(2, R) const &m2)
    {
        return AMAL_NMAT(2, R)(m1.data[0] + m2.data[0], m1.data[1] + m2.data[1]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(3, R) operator+(AMAL_NMAT(3, R) const &m1, AMAL_NMAT(3, R) const &m2)
    {
        return AMAL_NMAT(3, R)(m1.data[0] + m2.data[0], m1.data[1] + m2.data[1], m1.data[2] + m2.data[2]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(4, R) operator+(AMAL_NMAT(4, R) const &m1, AMAL_NMAT(4, R) const &m2)
    {
        return AMAL_NMAT(4, R)(m1.data[0] + m2.data[0], m1.data[1] + m2.data[1], m1.data[2] + m2.data[2],
                               m1.data[3] + m2.data[3]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(2, R) operator-(AMAL_NMAT(2, R) const &m1, AMAL_NMAT(2, R) const &m2)
    {
        return AMAL_NMAT(2, R)(m1.data[0] - m2.data[0], m1.data[1] - m2.data[1]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(3, R) operator-(AMAL_NMAT(3, R) const &m1, AMAL_NMAT(3, R) const &m2)
    {
        return AMAL_NMAT(3, R)(m1.data[0] - m2.data[0], m1.data[1] - m2.data[1], m1.data[2] - m2.data[2]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(4, R) operator-(AMAL_NMAT(4, R) const &m1, AMAL_NMAT(4, R) const &m2)
    {
        return AMAL_NMAT(4, R)(m1.data[0] - m2.data[0], m1.data[1] - m2.data[1], m1.data[2] - m2.data[2],
                               m1.data[3] - m2.data[3]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(2, R) operator*(AMAL_NMAT(2, R) const &self, T scalar)
    {
        return AMAL_NMAT(2, R)(self.data[0] * scalar, self.data[1] * scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(3, R) operator*(AMAL_NMAT(3, R) const &self, T scalar)
    {
        return AMAL_NMAT(3, R)(self.data[0] * scalar, self.data[1] * scalar, self.data[2] * scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(4, R) operator*(AMAL_NMAT(4, R) const &self, T scalar)
    {
        return AMAL_NMAT(4, R)(self.data[0] * scalar, self.data[1] * scalar, self.data[2] * scalar,
                               self.data[3] * scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(2, R) operator*(T scalar, AMAL_NMAT(2, R) const &self)
    {
        return AMAL_NMAT(2, R)(self.data[0] * scalar, self.data[1] * scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(3, R) operator*(T scalar, AMAL_NMAT(3, R) const &self)
    {
        return AMAL_NMAT(3, R)(self.data[0] * scalar, self.data[1] * scalar, self.data[2] * scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(4, R) operator*(T scalar, AMAL_NMAT(4, R) const &self)
    {
        return AMAL_NMAT(4, R)(self.data[0] * scalar, self.data[1] * scalar, self.data[2] * scalar,
                               self.data[3] * scalar);
    }

    template <AMAL_MATRIX_TEMP_DECL>
    inline AMAL_TYPE_SIMD(AMAL_MAT_SELF, typename AMAL_MAT_SELF::column_type)
    operator*(AMAL_MAT_SELF const &self, typename AMAL_MAT_SELF::row_type const &v)
    {
        using simd_type = typename AMAL_MAT_SELF::simd_type::value_type;
        simd_type s = internal::multiply_matrix(*reinterpret_cast<simd_type const(*)[C]>(self.data), v.s);
        return typename AMAL_MAT_SELF::column_type(s);
    }

#if defined(AMAL_FMA_ENABLE)
    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(2, 2), AMAL_NVEC(2)) operator*(AMAL_NMAT(2, 2) const &self, AMAL_NVEC(2) const &v)
    {
        return AMAL_NVEC(2)(std::fma(self[1][0], v.y, self[0][0] * v.x),  //
                            std::fma(self[1][1], v.y, self[0][1] * v.x)); //
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(2, 3), AMAL_NVEC(3)) operator*(AMAL_NMAT(2, 3) const &self, AMAL_NVEC(2) const &v)
    {
        return AMAL_NVEC(3)(std::fma(self[1][0], v.y, self[0][0] * v.x),  //
                            std::fma(self[1][1], v.y, self[0][1] * v.x),  //
                            std::fma(self[1][2], v.y, self[0][2] * v.x)); //
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(2, 4), AMAL_NVEC(4)) operator*(AMAL_NMAT(2, 4) const &self, AMAL_NVEC(2) const &v)
    {
        return AMAL_NVEC(4)(std::fma(self[1][0], v.y, self[0][0] * v.x),  //
                            std::fma(self[1][1], v.y, self[0][1] * v.x),  //
                            std::fma(self[1][2], v.y, self[0][2] * v.x),  //
                            std::fma(self[1][3], v.y, self[0][3] * v.x)); //
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(3, 2), AMAL_NVEC(2)) operator*(AMAL_NMAT(3, 2) const &self, AMAL_NVEC(3) const &v)
    {
        return AMAL_NVEC(2)(std::fma(self[2][0], v.z, std::fma(self[1][0], v.y, self[0][0] * v.x)),  //
                            std::fma(self[2][1], v.z, std::fma(self[1][1], v.y, self[0][1] * v.x))); //
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(3, 3), AMAL_NVEC(3)) operator*(AMAL_NMAT(3, 3) const &self, AMAL_NVEC(3) const &v)
    {
        return AMAL_NVEC(3)(std::fma(self[2][0], v.z, std::fma(self[1][0], v.y, self[0][0] * v.x)),  //
                            std::fma(self[2][1], v.z, std::fma(self[1][1], v.y, self[0][1] * v.x)),  //
                            std::fma(self[2][2], v.z, std::fma(self[1][2], v.y, self[0][2] * v.x))); //
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(3, 4), AMAL_NVEC(4)) operator*(AMAL_NMAT(3, 4) const &self, AMAL_NVEC(3) const &v)
    {
        return AMAL_NVEC(4)(std::fma(self[2][0], v.z, std::fma(self[1][0], v.y, self[0][0] * v.x)),  //
                            std::fma(self[2][1], v.z, std::fma(self[1][1], v.y, self[0][1] * v.x)),  //
                            std::fma(self[2][2], v.z, std::fma(self[1][2], v.y, self[0][2] * v.x)),  //
                            std::fma(self[2][3], v.z, std::fma(self[1][3], v.y, self[0][3] * v.x))); //
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(4, 2), AMAL_NVEC(2)) operator*(AMAL_NMAT(4, 2) const &self, AMAL_NVEC(4) const &v)
    {
        return AMAL_NVEC(2)(
            std::fma(self[3][0], v.w, std::fma(self[2][0], v.z, std::fma(self[1][0], v.y, self[0][0] * v.x))),  //
            std::fma(self[3][1], v.w, std::fma(self[2][1], v.z, std::fma(self[1][1], v.y, self[0][1] * v.x)))); //
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(4, 3), AMAL_NVEC(3)) operator*(AMAL_NMAT(4, 3) const &self, AMAL_NVEC(4) const &v)
    {
        return AMAL_NVEC(3)(
            std::fma(self[3][0], v.w, std::fma(self[2][0], v.z, std::fma(self[1][0], v.y, self[0][0] * v.x))),  //
            std::fma(self[3][1], v.w, std::fma(self[2][1], v.z, std::fma(self[1][1], v.y, self[0][1] * v.x))),  //
            std::fma(self[3][2], v.w, std::fma(self[2][2], v.z, std::fma(self[1][2], v.y, self[0][2] * v.x)))); //
    }

    template <typename T, bool aligned>
    inline AMAL_TYPE_NOSIMD(AMAL_NMAT(4, 4), AMAL_NVEC(4)) operator*(AMAL_NMAT(4, 4) const &self, AMAL_NVEC(4) const &v)
    {
        return AMAL_NVEC(4)(
            std::fma(self[3][0], v.w, std::fma(self[2][0], v.z, std::fma(self[1][0], v.y, self[0][0] * v.x))),
            std::fma(self[3][1], v.w, std::fma(self[2][1], v.z, std::fma(self[1][1], v.y, self[0][1] * v.x))),
            std::fma(self[3][2], v.w, std::fma(self[2][2], v.z, std::fma(self[1][2], v.y, self[0][2] * v.x))),
            std::fma(self[3][3], v.w, std::fma(self[2][3], v.z, std::fma(self[1][3], v.y, self[0][3] * v.x))));
    }

#else
    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(2, 2), AMAL_NVEC(2)) operator*(AMAL_NMAT(2, 2) const &self,
                                                                               AMAL_NVEC(2) const &v)
    {
        return AMAL_NVEC(2)(self[0][0] * v.x + self[1][0] * v.y, self[0][1] * v.x + self[1][1] * v.y);
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(2, 3), AMAL_NVEC(3)) operator*(AMAL_NMAT(2, 3) const &self,
                                                                               AMAL_NVEC(2) const &v)
    {
        return AMAL_NVEC(3)(self[0][0] * v.x + self[1][0] * v.y,  //
                            self[0][1] * v.x + self[1][1] * v.y,  //
                            self[0][2] * v.x + self[1][2] * v.y); //
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(2, 4), AMAL_NVEC(4)) operator*(AMAL_NMAT(2, 4) const &self,
                                                                               AMAL_NVEC(2) const &v)
    {
        return AMAL_NVEC(4)(self[0][0] * v.x + self[1][0] * v.y,  //
                            self[0][1] * v.x + self[1][1] * v.y,  //
                            self[0][2] * v.x + self[1][2] * v.y,  //
                            self[0][3] * v.x + self[1][3] * v.y); //
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(3, 2), AMAL_NVEC(2)) operator*(AMAL_NMAT(3, 2) const &self,
                                                                               AMAL_NVEC(3) const &v)
    {
        return AMAL_NVEC(2)(self[0][0] * v.x + self[1][0] * v.y + self[2][0] * v.z,  //
                            self[0][1] * v.x + self[1][1] * v.y + self[2][1] * v.z); //
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(3, 3), AMAL_NVEC(3)) operator*(AMAL_NMAT(3, 3) const &self,
                                                                               AMAL_NVEC(3) const &v)
    {
        return AMAL_NVEC(3)(self[0][0] * v.x + self[1][0] * v.y + self[2][0] * v.z,  //
                            self[0][1] * v.x + self[1][1] * v.y + self[2][1] * v.z,  //
                            self[0][2] * v.x + self[1][2] * v.y + self[2][2] * v.z); //
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(3, 4), AMAL_NVEC(4)) operator*(AMAL_NMAT(3, 4) const &self,
                                                                               AMAL_NVEC(3) const &v)
    {
        return AMAL_NVEC(4)(self[0][0] * v.x + self[1][0] * v.y + self[2][0] * v.z,  //
                            self[0][1] * v.x + self[1][1] * v.y + self[2][1] * v.z,  //
                            self[0][2] * v.x + self[1][2] * v.y + self[2][2] * v.z,  //
                            self[0][3] * v.x + self[1][3] * v.y + self[2][3] * v.z); //
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(4, 2), AMAL_NVEC(2)) operator*(AMAL_NMAT(4, 2) const &self,
                                                                               AMAL_NVEC(4) const &v)
    {
        return AMAL_NVEC(2)(self[0][0] * v.x + self[1][0] * v.y + self[2][0] * v.z + self[3][0] * v.w,  //
                            self[0][1] * v.x + self[1][1] * v.y + self[2][1] * v.z + self[3][1] * v.w); //
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(4, 3), AMAL_NVEC(3)) operator*(AMAL_NMAT(4, 3) const &self,
                                                                               AMAL_NVEC(4) const &v)
    {
        return AMAL_NVEC(3)(self[0][0] * v.x + self[1][0] * v.y + self[2][0] * v.z + self[3][0] * v.w,  //
                            self[0][1] * v.x + self[1][1] * v.y + self[2][1] * v.z + self[3][1] * v.w,  //
                            self[0][2] * v.x + self[1][2] * v.y + self[2][2] * v.z + self[3][2] * v.w); //
    }

    template <typename T, bool aligned>
    inline constexpr AMAL_TYPE_NOSIMD(AMAL_NMAT(4, 4), AMAL_NVEC(4)) operator*(AMAL_NMAT(4, 4) const &self,
                                                                               AMAL_NVEC(4) const &v)
    {
        return AMAL_NVEC(4)(self[0][0] * v.x + self[1][0] * v.y + self[2][0] * v.z + self[3][0] * v.w,  //
                            self[0][1] * v.x + self[1][1] * v.y + self[2][1] * v.z + self[3][1] * v.w,  //
                            self[0][2] * v.x + self[1][2] * v.y + self[2][2] * v.z + self[3][2] * v.w,  //
                            self[0][3] * v.x + self[1][3] * v.y + self[2][3] * v.z + self[3][3] * v.w); //
    }
#endif

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR typename AMAL_NMAT(2, R)::row_type operator*(typename AMAL_NMAT(2, R)::column_type const &v,
                                                                       AMAL_NMAT(2, R) const &m)
    {
        return typename AMAL_NMAT(2, R)::row_type(dot(m.data[0], v), dot(m.data[1], v));
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR typename AMAL_NMAT(3, R)::row_type operator*(typename AMAL_NMAT(3, R)::column_type const &v,
                                                                       AMAL_NMAT(3, R) const &m)
    {
        return typename AMAL_NMAT(3, R)::row_type(dot(m.data[0], v), dot(m.data[1], v), dot(m.data[2], v));
    }

    template <length_t R, typename T, bool aligned>

    inline AMAL_CONSTEXPR typename AMAL_NMAT(4, R)::row_type operator*(typename AMAL_NMAT(4, R)::column_type const &v,
                                                                       AMAL_NMAT(4, R) const &m)
    {
        return typename AMAL_NMAT(4, R)::row_type(dot(m.data[0], v), dot(m.data[1], v), dot(m.data[2], v),
                                                  dot(m.data[3], v));
    }

    template <length_t C1, length_t R1, length_t C2, length_t R2, typename T, bool aligned>
    inline AMAL_MAT_MUL_SIMD operator*(AMAL_NMAT(C1, R1) const &m1, AMAL_NMAT(C2, R2) const &m2)
    {
        static_assert(R1 == C2, "Invalid matrix dimensions for multiplication");
        using simd_t = typename AMAL_NMAT(C2, R1)::simd_type::value_type;
        simd_t out[C2];
        if constexpr (R1 < 4) __builtin_memset(out, 0, sizeof(out));
        internal::multiply_matrix<C1, R1, C2, R2>(*reinterpret_cast<simd_t const(*)[C1]>(m1.data),
                                                  *reinterpret_cast<simd_t const(*)[C2]>(m2.data), out);
        return AMAL_NMAT(C2, R1)(out);
    }

#if defined(AMAL_FMA_ENABLE)
    #include <amal/internal/matrix_multiply_fma.hpp>
#else
    #include <amal/internal/matrix_multiply_linear.hpp>
#endif

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(2, R) operator/(AMAL_NMAT(2, R) const &self, T scalar)
    {
        return AMAL_NMAT(2, R)(self.data[0] / scalar, self.data[1] / scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(3, R) operator/(AMAL_NMAT(3, R) const &self, T scalar)
    {
        return AMAL_NMAT(3, R)(self.data[0] / scalar, self.data[1] / scalar, self.data[2] / scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(4, R) operator/(AMAL_NMAT(4, R) const &self, T scalar)
    {
        return AMAL_NMAT(4, R)(self.data[0] / scalar, self.data[1] / scalar, self.data[2] / scalar,
                               self.data[3] / scalar);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(2, R) operator/(T scalar, AMAL_NMAT(2, R) const &self)
    {
        return AMAL_NMAT(2, R)(scalar / self.data[0], scalar / self.data[1]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(3, R) operator/(T scalar, AMAL_NMAT(3, R) const &self)
    {
        return AMAL_NMAT(3, R)(scalar / self.data[0], scalar / self.data[1], scalar / self.data[2]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NMAT(4, R) operator/(T scalar, AMAL_NMAT(4, R) const &self)
    {
        return AMAL_NMAT(4, R)(scalar / self.data[0], scalar / self.data[1], scalar / self.data[2],
                               scalar / self.data[3]);
    }

    template <AMAL_MATRIX_TEMP_DECL>
    inline AMAL_CONSTEXPR typename  AMAL_MAT_SELF::column_type operator/(AMAL_MAT_SELF const &self,
                                                               typename AMAL_MAT_SELF::row_type const &v)
    {
        return inverse_matrix(self) * v;
    }

    template <AMAL_MATRIX_TEMP_DECL>
    inline AMAL_CONSTEXPR typename AMAL_MAT_SELF::row_type operator/(typename AMAL_MAT_SELF::column_type const &v,
                                                            AMAL_MAT_SELF const &self)
    {
        return v * inverse_matrix(self);
    }

    template <AMAL_MATRIX_TEMP_DECL>
    inline AMAL_CONSTEXPR AMAL_MAT_SELF operator/(AMAL_MAT_SELF const &m1, AMAL_MAT_SELF const &m2)
    {
        AMAL_MAT_SELF tmp(m1);
        return tmp /= m2;
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR bool operator==(AMAL_NMAT(2, R) const &m1, AMAL_NMAT(2, R) const &m2)
    {
        return (m1.data[0] == m2.data[0]) && (m1.data[1] == m2.data[1]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR bool operator==(AMAL_NMAT(3, R) const &m1, AMAL_NMAT(3, R) const &m2)
    {
        return (m1.data[0] == m2.data[0]) && (m1.data[1] == m2.data[1]) && (m1.data[2] == m2.data[2]);
    }

    template <length_t R, typename T, bool aligned>
    inline AMAL_CONSTEXPR bool operator==(AMAL_NMAT(4, R) const &m1, AMAL_NMAT(4, R) const &m2)
    {
        return (m1.data[0] == m2.data[0]) && (m1.data[1] == m2.data[1]) && (m1.data[2] == m2.data[2]) &&
               (m1.data[3] == m2.data[3]);
    }
} // namespace amal