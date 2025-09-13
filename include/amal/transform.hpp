#pragma once

#include "geometric.hpp"
#include "matrix.hpp"

namespace amal
{
    namespace internal
    {
        template <length_t N, typename T, bool aligned>
        inline constexpr bool is_mat4_nvec_simdable =
            is_simd_enabled_v<typename mat<4, 4, T, true>::simd_type::value_type> &&
            is_simd_enabled_v<typename vec<N, T, aligned>::simd_type::value_type>;

#define AMAL_MAT4_SIMD(N)   std::enable_if_t<internal::is_mat4_nvec_simdable<N, T, aligned>, mat4>
#define AMAL_MAT4_NOSIMD(N) std::enable_if_t<!internal::is_mat4_nvec_simdable<N, T, aligned>, mat4>
    } // namespace internal

    template <typename T, bool aligned>
    inline AMAL_NVEC(2) screen_to_ndc(const AMAL_NVEC(2) & screen, int x, int y)
    {
        T x_ndc = (screen.x / x) * static_cast<T>(2) - static_cast<T>(1);
        T y_ndc = 1.0f - (screen.y / y) * static_cast<T>(2);
        return AMAL_NVEC(2)(x_ndc, y_ndc);
    }

    template <typename T, bool aligned>
    inline AMAL_CONSTEXPR AMAL_NVEC(3)
        project(AMAL_NVEC(3) const &obj, mat4 const &model, mat4 const &proj, AMAL_VEC(4, T, true) const &viewport)
    {
        AMAL_VEC(4, T, true) tmp(obj, static_cast<T>(1));
        tmp = model * tmp;
        tmp = proj * tmp;

        tmp /= tmp.w;
#ifdef AMAL_CLIP_SPACE_NO
        tmp = tmp * static_cast<T>(0.5) + static_cast<T>(0.5);
#else
        tmp.x = tmp.x * static_cast<T>(0.5) + static_cast<T>(0.5);
        tmp.y = tmp.y * static_cast<T>(0.5) + static_cast<T>(0.5);
#endif
#ifdef AMAL_FMA_ENABLE
        tmp[0] = fma(tmp[0], T(viewport[2]), T(viewport[0]));
        tmp[1] = fma(tmp[1], T(viewport[3]), T(viewport[1]));
#else
        tmp[0] = tmp[0] * T(viewport[2]) + T(viewport[0]);
        tmp[1] = tmp[1] * T(viewport[3]) + T(viewport[1]);
#endif

        return AMAL_NVEC(3)(tmp);
    }

    template <typename T, bool aligned>
    inline AMAL_NVEC(3)
        unproject(AMAL_NVEC(3) const &win, mat4 const &model, mat4 const &proj, AMAL_VEC(4, T, true) const &viewport)
    {
        mat4 inverse = inverse_matrix(proj * model);
        AMAL_VEC(4, T, true) tmp(win, T(1));
        tmp.x = (tmp.x - T(viewport[0])) / T(viewport[2]);
        tmp.y = (tmp.y - T(viewport[1])) / T(viewport[3]);
#ifdef AMAL_CLIP_SPACE_NO
        tmp = tmp * static_cast<T>(2) - static_cast<T>(1);
#else
        tmp.x = tmp.x * static_cast<T>(2) - static_cast<T>(1);
        tmp.y = tmp.y * static_cast<T>(2) - static_cast<T>(1);
#endif

        AMAL_VEC(4, T, true) obj = inverse * tmp;
        obj /= obj.w;
        return AMAL_NVEC(3)(obj);
    }

    template <typename T, bool aligned>
    inline constexpr mat4 translate(mat4 const &m, AMAL_VEC(3, T, aligned) const &v)
    {
        mat4 r(m);
        r[3] = m[0] * v[0] + m[1] * v[1] + m[2] * v[2] + m[3];
        return r;
    }

    template<typename T, bool aligned>
    inline constexpr mat4 translate(AMAL_VEC(3, T, aligned) const &v)
    {
        return translate(mat4(static_cast<T>(1)), v);
    }

    template <typename T, bool aligned>
    inline AMAL_MAT4_SIMD(3) rotate(mat4 const &m, T angle, AMAL_NVEC(3) const &axis)
    {
        static_assert(is_floating_point_v<T>, "rotate only supports floating point types");
        using simd_t = typename mat4::simd_type::value_type;
        simd_t out[4];
        __builtin_memset(out, 0, sizeof(out));
        internal::rotate(*reinterpret_cast<simd_t const(*)[4]>(m.data), angle, axis.s, out);
        return mat4(out);
    }

    template <typename T, bool aligned>
    inline AMAL_MAT4_SIMD(4) rotate(mat4 const &m, T angle, AMAL_NVEC(4) const &axis)
    {
        static_assert(is_floating_point_v<T>, "rotate only supports floating point types");
        using simd_t = typename mat4::simd_type::value_type;
        simd_t out[4];
        internal::rotate(*reinterpret_cast<simd_t const(*)[4]>(m.data), angle, axis.s, out);
        return mat4(out);
    }

#if defined(AMAL_FMA_ENABLE)
    template <typename T, bool aligned>
    inline AMAL_MAT4_NOSIMD(3) rotate(mat4 const &m, T angle, AMAL_NVEC(3) const &v)
    {
        static_assert(is_floating_point_v<T>, "rotate only supports floating point types");
        T const c = cos(angle);
        T const s = sin(angle);
        AMAL_NVEC(3) axis = normalize(v);
        AMAL_NVEC(3) temp = (T(1) - c) * axis;

        mat4 rotate;

        rotate[0][0] = fma(temp[0], axis[0], c);
        rotate[0][1] = fma(temp[0], axis[1], s * axis[2]);
        rotate[0][2] = fma(temp[0], axis[2], -s * axis[1]);

        rotate[1][0] = fma(temp[1], axis[0], -s * axis[2]);
        rotate[1][1] = fma(temp[1], axis[1], c);
        rotate[1][2] = fma(temp[1], axis[2], s * axis[0]);

        rotate[2][0] = fma(temp[2], axis[0], s * axis[1]);
        rotate[2][1] = fma(temp[2], axis[1], -s * axis[0]);
        rotate[2][2] = fma(temp[2], axis[2], c);

        mat4 r;

        r[0] = m[0] * rotate[0][0] + m[1] * rotate[0][1] + m[2] * rotate[0][2];
        r[1] = m[0] * rotate[1][0] + m[1] * rotate[1][1] + m[2] * rotate[1][2];
        r[2] = m[0] * rotate[2][0] + m[1] * rotate[2][1] + m[2] * rotate[2][2];
        r[3] = m[3];

        return r;
    }
#else
    template <typename T, bool aligned>
    inline AMAL_MAT4_NOSIMD(4) rotate(mat4 const &m, T angle, AMAL_NVEC(3) const &v)
    {
        static_assert(is_floating_point_v<T>, "rotate only supports floating point types");
        T const a = angle;
        T const c = cos(a);
        T const s = sin(a);

        AMAL_NVEC(3) axis(normalize(v));
        AMAL_NVEC(3) temp((T(1) - c) * axis);

        mat4 rotate;
        rotate[0][0] = c + temp[0] * axis[0];
        rotate[0][1] = temp[0] * axis[1] + s * axis[2];
        rotate[0][2] = temp[0] * axis[2] - s * axis[1];

        rotate[1][0] = temp[1] * axis[0] - s * axis[2];
        rotate[1][1] = c + temp[1] * axis[1];
        rotate[1][2] = temp[1] * axis[2] + s * axis[0];

        rotate[2][0] = temp[2] * axis[0] + s * axis[1];
        rotate[2][1] = temp[2] * axis[1] - s * axis[0];
        rotate[2][2] = c + temp[2] * axis[2];

        mat4 r;
        r[0] = m[0] * rotate[0][0] + m[1] * rotate[0][1] + m[2] * rotate[0][2];
        r[1] = m[0] * rotate[1][0] + m[1] * rotate[1][1] + m[2] * rotate[1][2];
        r[2] = m[0] * rotate[2][0] + m[1] * rotate[2][1] + m[2] * rotate[2][2];
        r[3] = m[3];
        return r;
    }
#endif

    template <typename T, bool aligned>
    inline AMAL_MAT4_SIMD(3) scale(mat4 const &m, AMAL_NVEC(3) const &axis)
    {
        static_assert(is_floating_point_v<T>, "scale only supports floating point types");
        using simd_t = typename mat4::simd_type::value_type;
        simd_t out[4];
        __builtin_memset(out, 0, sizeof(out));
        internal::scale(*reinterpret_cast<simd_t const(*)[4]>(m.data), axis.s, out);
        return mat4(out);
    }

    template <typename T, bool aligned>
    inline AMAL_MAT4_SIMD(4) scale(mat4 const &m, AMAL_NVEC(4) const &axis)
    {
        static_assert(is_floating_point_v<T>, "scale only supports floating point types");
        using simd_t = typename mat4::simd_type::value_type;
        simd_t out[4];
        internal::scale(*reinterpret_cast<simd_t const(*)[4]>(m.data), axis.s, out);
        return mat4(out);
    }

    template <typename T, bool aligned>
    inline AMAL_MAT4_NOSIMD(3) scale(mat4 const &m, AMAL_NVEC(3) const &axis)
    {
        static_assert(is_floating_point_v<T>, "scale only supports floating point types");
        mat4 out;
        out[0] = m[0] * axis[0];
        out[1] = m[1] * axis[1];
        out[2] = m[2] * axis[2];
        out[3] = m[3];
        return out;
    }

    template <typename T, bool aligned>
    inline AMAL_MAT4_SIMD(2) shear(mat4 const &m, AMAL_VEC(3, T, true) const &p, AMAL_VEC(2, T, true) const &l_x,
                                   AMAL_VEC(2, T, true) const &l_y, AMAL_VEC(2, T, true) const &l_z)
    {
        static_assert(is_floating_point_v<T>, "shear only supports floating point types");
        using simd_t = typename mat4::simd_type::value_type;
        simd_t out[4];
        __builtin_memset(out, 0, sizeof(out));
        internal::shear(*reinterpret_cast<simd_t const(*)[4]>(m.data), p.s, l_x.s, l_y.s, l_z.s, out);
        return mat4(out);
    }

    template <typename T, bool aligned>
    inline AMAL_MAT4_SIMD(2) shear(mat4 const &m, AMAL_NVEC(4) const &p, AMAL_VEC(2, T, true) const &l_x,
                                   AMAL_VEC(2, T, true) const &l_y, AMAL_VEC(2, T, true) const &l_z)
    {
        static_assert(is_floating_point_v<T>, "shear only supports floating point types");
        using simd_t = typename mat4::simd_type::value_type;
        simd_t out[4];
        internal::shear(*reinterpret_cast<simd_t const(*)[4]>(m.data), p.s, l_x.s, l_y.s, l_z.s, out);
        return mat4(out);
    }

    template <typename T, bool aligned>
    inline AMAL_MAT4_NOSIMD(2) shear(mat4 const &m, AMAL_NVEC(3) const &p, AMAL_VEC(2, T, true) const &l_x,
                                     AMAL_VEC(2, T, true) const &l_y, AMAL_VEC(2, T, true) const &l_z)
    {
        static_assert(is_floating_point_v<T>, "shear only supports floating point types");
        T const lambda_xy = l_x[0];
        T const lambda_xz = l_x[1];
        T const lambda_yx = l_y[0];
        T const lambda_yz = l_y[1];
        T const lambda_zx = l_z[0];
        T const lambda_zy = l_z[1];

        AMAL_NVEC(3) point_lambda((lambda_xy + lambda_xz), (lambda_yx + lambda_yz), (lambda_zx + lambda_zy));

        mat4 shear(1, lambda_yx, lambda_zx, 0, lambda_xy, 1, lambda_zy, 0, lambda_xz, lambda_yz, 1, 0,
                   -point_lambda[0] * p[0], -point_lambda[1] * p[1], -point_lambda[2] * p[2], 1);

        mat4 r(1);
        r[0] = m[0] * shear[0][0] + m[1] * shear[0][1] + m[2] * shear[0][2] + m[3] * shear[0][3];
        r[1] = m[0] * shear[1][0] + m[1] * shear[1][1] + m[2] * shear[1][2] + m[3] * shear[1][3];
        r[2] = m[0] * shear[2][0] + m[1] * shear[2][1] + m[2] * shear[2][2] + m[3] * shear[2][3];
        r[3] = m[0] * shear[3][0] + m[1] * shear[3][1] + m[2] * shear[3][2] + m[3] * shear[3][3];
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_CONSTEXPR mat4 look_at_rh(AMAL_NVEC(3) const &eye, AMAL_NVEC(3) const &center, AMAL_NVEC(3) const &up)
    {
        AMAL_NVEC(3) const f(normalize(center - eye));
        AMAL_NVEC(3) const s(normalize(cross(f, up)));
        AMAL_NVEC(3) const u(cross(s, f));

        mat4 r(1);
        r[0][0] = s.x;
        r[1][0] = s.y;
        r[2][0] = s.z;
        r[0][1] = u.x;
        r[1][1] = u.y;
        r[2][1] = u.z;
        r[0][2] = -f.x;
        r[1][2] = -f.y;
        r[2][2] = -f.z;
        r[3][0] = -dot(s, eye);
        r[3][1] = -dot(u, eye);
        r[3][2] = dot(f, eye);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_CONSTEXPR mat4 look_at_lh(AMAL_NVEC(3) const &eye, AMAL_NVEC(3) const &center, AMAL_NVEC(3) const &up)
    {
        AMAL_NVEC(3) const f(normalize(center - eye));
        AMAL_NVEC(3) const s(normalize(cross(up, f)));
        AMAL_NVEC(3) const u(cross(f, s));

        mat4 r;
        r[0][0] = s.x;
        r[1][0] = s.y;
        r[2][0] = s.z;
        r[0][1] = u.x;
        r[1][1] = u.y;
        r[2][1] = u.z;
        r[0][2] = f.x;
        r[1][2] = f.y;
        r[2][2] = f.z;
        r[3][0] = -dot(s, eye);
        r[3][1] = -dot(u, eye);
        r[3][2] = -dot(f, eye);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_CONSTEXPR mat4 look_at(AMAL_NVEC(3) const &eye, AMAL_NVEC(3) const &center, AMAL_NVEC(3) const &up)
    {
#ifdef AMAL_CLIP_SPACE_LH
        return look_at_lh<T, aligned>(eye, center, up);
#else
        return look_at_rh<T, aligned>(eye, center, up);
#endif
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 orthographic(T left, T right, T bottom, T top)
    {
        mat4 r(static_cast<T>(1));
        r[0][0] = static_cast<T>(2) / (right - left);
        r[1][1] = static_cast<T>(2) / (top - bottom);
        r[2][2] = -static_cast<T>(1);
        r[3][0] = -(right + left) / (right - left);
        r[3][1] = -(top + bottom) / (top - bottom);
        return r;
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 orthographic_lh_zo(T left, T right, T bottom, T top, T z_near, T z_far)
    {
        mat4 r(1);
        r[0][0] = static_cast<T>(2) / (right - left);
        r[1][1] = static_cast<T>(2) / (top - bottom);
        r[2][2] = static_cast<T>(1) / (z_far - z_near);
        r[3][0] = -(right + left) / (right - left);
        r[3][1] = -(top + bottom) / (top - bottom);
        r[3][2] = -z_near / (z_far - z_near);
        return r;
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 orthographic_lh_no(T left, T right, T bottom, T top, T z_near, T z_far)
    {
        mat4 r(1);
        r[0][0] = static_cast<T>(2) / (right - left);
        r[1][1] = static_cast<T>(2) / (top - bottom);
        r[2][2] = static_cast<T>(2) / (z_far - z_near);
        r[3][0] = -(right + left) / (right - left);
        r[3][1] = -(top + bottom) / (top - bottom);
        r[3][2] = -(z_far + z_near) / (z_far - z_near);
        return r;
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 orthographic_rh_zo(T left, T right, T bottom, T top, T z_near, T z_far)
    {
        mat4 r(1);
        r[0][0] = static_cast<T>(2) / (right - left);
        r[1][1] = static_cast<T>(2) / (top - bottom);
        r[2][2] = -static_cast<T>(1) / (z_far - z_near);
        r[3][0] = -(right + left) / (right - left);
        r[3][1] = -(top + bottom) / (top - bottom);
        r[3][2] = -z_near / (z_far - z_near);
        return r;
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 orthographic_rh_no(T left, T right, T bottom, T top, T z_near, T z_far)
    {
        mat4 r(1);
        r[0][0] = static_cast<T>(2) / (right - left);
        r[1][1] = static_cast<T>(2) / (top - bottom);
        r[2][2] = -static_cast<T>(2) / (z_far - z_near);
        r[3][0] = -(right + left) / (right - left);
        r[3][1] = -(top + bottom) / (top - bottom);
        r[3][2] = -(z_far + z_near) / (z_far - z_near);
        return r;
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 orthographic(T left, T right, T bottom, T top, T z_near, T z_far)
    {
#if defined(AMAL_CLIP_SPACE_LH)
    #if defined(AMAL_CLIP_SPACE_NO)
        return orthographic_lh_no(left, right, bottom, top, z_near, z_far);
    #else
        return orthographic_lh_zo(left, right, bottom, top, z_near, z_far);
    #endif
#else
    #if defined(AMAL_CLIP_SPACE_NO)
        return orthographic_rh_no(left, right, bottom, top, z_near, z_far);
    #else
        return orthographic_rh_zo(left, right, bottom, top, z_near, z_far);
    #endif
#endif
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 perspective_rh_zo(T fovy, T aspect, T z_near, T z_far)
    {
        assert(abs(aspect - std::numeric_limits<T>::epsilon()) > static_cast<T>(0));
        T const tan_half_fovy = tan(fovy / static_cast<T>(2));

        mat4 r(static_cast<T>(0));
        r[0][0] = static_cast<T>(1) / (aspect * tan_half_fovy);
        r[1][1] = static_cast<T>(1) / (tan_half_fovy);
        r[2][2] = z_far / (z_near - z_far);
        r[2][3] = -static_cast<T>(1);
        r[3][2] = -(z_far * z_near) / (z_far - z_near);
        return r;
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 perspective_rh_no(T fovy, T aspect, T z_near, T z_far)
    {
        assert(abs(aspect - std::numeric_limits<T>::epsilon()) > static_cast<T>(0));
        T const tan_half_fovy = tan(fovy / static_cast<T>(2));

        mat4 r(static_cast<T>(0));
        r[0][0] = static_cast<T>(1) / (aspect * tan_half_fovy);
        r[1][1] = static_cast<T>(1) / (tan_half_fovy);
        r[2][2] = -(z_far + z_near) / (z_far - z_near);
        r[2][3] = -static_cast<T>(1);
        r[3][2] = -(static_cast<T>(2) * z_far * z_near) / (z_far - z_near);
        return r;
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 perspective_lh_zo(T fovy, T aspect, T z_near, T z_far)
    {
        assert(abs(aspect - std::numeric_limits<T>::epsilon()) > static_cast<T>(0));
        T const tan_half_fovy = tan(fovy / static_cast<T>(2));

        mat4 r(static_cast<T>(0));
        r[0][0] = static_cast<T>(1) / (aspect * tan_half_fovy);
        r[1][1] = static_cast<T>(1) / (tan_half_fovy);
        r[2][2] = z_far / (z_far - z_near);
        r[2][3] = static_cast<T>(1);
        r[3][2] = -(z_far * z_near) / (z_far - z_near);
        return r;
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 perspective_lh_no(T fovy, T aspect, T z_near, T z_far)
    {
        assert(abs(aspect - std::numeric_limits<T>::epsilon()) > static_cast<T>(0));
        T const tan_half_fovy = tan(fovy / static_cast<T>(2));

        mat4 r(static_cast<T>(0));
        r[0][0] = static_cast<T>(1) / (aspect * tan_half_fovy);
        r[1][1] = static_cast<T>(1) / (tan_half_fovy);
        r[2][2] = (z_far + z_near) / (z_far - z_near);
        r[2][3] = static_cast<T>(1);
        r[3][2] = -(static_cast<T>(2) * z_far * z_near) / (z_far - z_near);
        return r;
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 perspective(T fovy, T aspect, T z_near, T z_far)
    {
#if defined(AMAL_CLIP_SPACE_LH)
    #if defined(AMAL_CLIP_SPACE_NO)
        return perspective_lh_no<T>(fovy, aspect, z_near, z_far);
    #else
        return perspective_lh_zo<T>(fovy, aspect, z_near, z_far);
    #endif
#else
    #if defined(AMAL_CLIP_SPACE_NO)
        return perspective_rh_no<T>(fovy, aspect, z_near, z_far);
    #else
        return perspective_rh_zo<T>(fovy, aspect, z_near, z_far);
    #endif
#endif
    }

    template <typename T>
    inline AMAL_CONSTEXPR mat4 perspective_infinite_rh_no(T fovy, T aspect, T z_near, T ep)
    {
        const T range = tan(fovy * T(0.5)) * z_near;
        const T left = -range * aspect;
        const T right = range * aspect;
        const T bottom = -range;
        const T top = range;

        mat4 r(T(0));
        r[0][0] = (T(2) * z_near) / (right - left);
        r[1][1] = (T(2) * z_near) / (top - bottom);
        r[2][2] = ep - T(1);
        r[2][3] = T(-1);
        r[3][2] = (ep - T(2)) * z_near;
        return r;
    }

    template <typename T>
    inline constexpr mat4 perspective_infinite_rh_zo(T fovy, T aspect, T z_near, T ep)
    {
        const T range = tan(fovy * T(0.5)) * z_near;
        const T left = -range * aspect;
        const T right = range * aspect;
        const T bottom = -range;
        const T top = range;

        mat4 r(T(0));
        r[0][0] = (T(2) * z_near) / (right - left);
        r[1][1] = (T(2) * z_near) / (top - bottom);
        r[2][2] = T(-1) * (ep - T(1));
        r[2][3] = T(-1);
        r[3][2] = T(-1) * (ep - T(2)) * z_near;
        return r;
    }

    template <typename T>
    inline constexpr mat4 perspective_infinite_lh_no(T fovy, T aspect, T z_near, T ep)
    {
        const T range = tan(fovy * T(0.5)) * z_near;
        const T left = -range * aspect;
        const T right = range * aspect;
        const T bottom = -range;
        const T top = range;

        mat4 r(T(0));
        r[0][0] = (T(2) * z_near) / (right - left);
        r[1][1] = (T(2) * z_near) / (top - bottom);
        r[2][2] = ep - T(1);
        r[2][3] = T(1);
        r[3][2] = T(-1) * (ep - T(2)) * z_near;
        return r;
    }

    template <typename T>
    inline constexpr mat4 perspective_infinite_lh_zo(T fovy, T aspect, T z_near, T ep)
    {
        const T range = tan(fovy * T(0.5)) * z_near;
        const T left = -range * aspect;
        const T right = range * aspect;
        const T bottom = -range;
        const T top = range;

        mat4 r(T(0));
        r[0][0] = (T(2) * z_near) / (right - left);
        r[1][1] = (T(2) * z_near) / (top - bottom);
        r[2][2] = (ep - T(1));
        r[2][3] = T(1);
        r[3][2] = T(-1) * (ep - T(2)) * z_near;
        return r;
    }

    // @ref: Infinite projection matrix: http://www.terathon.com/gdc07_lengyel.pdf
    template <typename T>
    inline constexpr mat4 perspective_infinite(T fovy, T aspect, T z_near, T ep = std::numeric_limits<T>::epsilon())
    {
#if defined(AMAL_CLIP_SPACE_LH)
    #if defined(AMAL_CLIP_SPACE_NO)
        return perspective_infinite_lh_no<T>(fovy, aspect, z_near, ep);
    #else
        return perspective_infinite_lh_zo<T>(fovy, aspect, z_near, ep);
    #endif
#else
    #if defined(AMAL_CLIP_SPACE_NO)
        return perspective_infinite_rh_no<T>(fovy, aspect, z_near, ep);
    #else
        return perspective_infinite_rh_zo<T>(fovy, aspect, z_near, ep);
    #endif
#endif
    }

} // namespace amal