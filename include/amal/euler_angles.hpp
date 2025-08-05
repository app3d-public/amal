#pragma once

#include "matrix.hpp"

namespace amal
{
    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_x(T const &angle_x)
    {
        T cos_x = cos(angle_x);
        T sin_x = sin(angle_x);
        return AMAL_NMAT(4, 4)(1, 0, 0, 0, 0, cos_x, sin_x, 0, 0, -sin_x, cos_x, 0, 0, 0, 0, 1);
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_y(T const &angle_y)
    {
        T cos_y = cos(angle_y);
        T sin_y = sin(angle_y);
        return AMAL_NMAT(4, 4)(cos_y, 0, -sin_y, 0, 0, 1, 0, 0, sin_y, 0, cos_y, 0, 0, 0, 0, 1);
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_z(T const &angle_z)
    {
        T cos_z = cos(angle_z);
        T sin_z = sin(angle_z);
        return AMAL_NMAT(4, 4)(cos_z, sin_z, 0, 0, -sin_z, cos_z, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_xy(T const &angle_x, T const &angle_y)
    {
        T cos_x = cos(angle_x);
        T sin_x = sin(angle_x);
        T cos_y = cos(angle_y);
        T sin_y = sin(angle_y);
        return AMAL_NMAT(4, 4)(cos_y, -sin_x * -sin_y, cos_x * -sin_y, 0, 0, cos_x, sin_x, 0, sin_y, -sin_x * cos_y,
                               cos_x * cos_y, 0, 0, 0, 0, 1);
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_yx(T const &angle_y, T const &angle_x)
    {
        T cos_x = cos(angle_x);
        T sin_x = sin(angle_x);
        T cos_y = cos(angle_y);
        T sin_y = sin(angle_y);
        return AMAL_NMAT(4, 4)(cos_y, 0, -sin_y, 0, sin_y * sin_x, cos_x, cos_y * sin_x, 0, sin_y * cos_x, -sin_x,
                               cos_y * cos_x, 0, 0, 0, 0, 1);
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_xz(T const &angle_x, T const &angle_z)
    {
        return euler_angle_x(angle_x) * euler_angle_z(angle_z);
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_zx(T const &angle_z, T const &angle_x)
    {
        return euler_angle_z(angle_z) * euler_angle_x(angle_x);
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_yz(T const &angle_y, T const &angle_z)
    {
        return euler_angle_y(angle_y) * euler_angle_z(angle_z);
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_zy(T const &angle_z, T const &angle_y)
    {
        return euler_angle_z(angle_z) * euler_angle_y(angle_y);
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_xyz(T const &t1, T const &t2, T const &t3)
    {
        T c1 = cos(-t1);
        T c2 = cos(-t2);
        T c3 = cos(-t3);
        T s1 = sin(-t1);
        T s2 = sin(-t2);
        T s3 = sin(-t3);

        AMAL_NMAT(4, 4) r;
        r[0][0] = c2 * c3;
        r[0][1] = -c1 * s3 + s1 * s2 * c3;
        r[0][2] = s1 * s3 + c1 * s2 * c3;
        r[0][3] = static_cast<T>(0);
        r[1][0] = c2 * s3;
        r[1][1] = c1 * c3 + s1 * s2 * s3;
        r[1][2] = -s1 * c3 + c1 * s2 * s3;
        r[1][3] = static_cast<T>(0);
        r[2][0] = -s2;
        r[2][1] = s1 * c2;
        r[2][2] = c1 * c2;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_yxz(T const &yaw, T const &pitch, T const &roll)
    {
        T tmp_ch = cos(yaw);
        T tmp_sh = sin(yaw);
        T tmp_cp = cos(pitch);
        T tmp_sp = sin(pitch);
        T tmp_cb = cos(roll);
        T tmp_sb = sin(roll);

        AMAL_NMAT(4, 4) r;
        r[0][0] = tmp_ch * tmp_cb + tmp_sh * tmp_sp * tmp_sb;
        r[0][1] = tmp_sb * tmp_cp;
        r[0][2] = -tmp_sh * tmp_cb + tmp_ch * tmp_sp * tmp_sb;
        r[0][3] = static_cast<T>(0);
        r[1][0] = -tmp_ch * tmp_sb + tmp_sh * tmp_sp * tmp_cb;
        r[1][1] = tmp_cb * tmp_cp;
        r[1][2] = tmp_sb * tmp_sh + tmp_ch * tmp_sp * tmp_cb;
        r[1][3] = static_cast<T>(0);
        r[2][0] = tmp_sh * tmp_cp;
        r[2][1] = -tmp_sp;
        r[2][2] = tmp_ch * tmp_cp;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_xzx(T const &t1, T const &t2, T const &t3)
    {
        T c1 = cos(t1);
        T s1 = sin(t1);
        T c2 = cos(t2);
        T s2 = sin(t2);
        T c3 = cos(t3);
        T s3 = sin(t3);

        AMAL_NMAT(4, 4) r;
        r[0][0] = c2;
        r[0][1] = c1 * s2;
        r[0][2] = s1 * s2;
        r[0][3] = static_cast<T>(0);
        r[1][0] = -c3 * s2;
        r[1][1] = c1 * c2 * c3 - s1 * s3;
        r[1][2] = c1 * s3 + c2 * c3 * s1;
        r[1][3] = static_cast<T>(0);
        r[2][0] = s2 * s3;
        r[2][1] = -c3 * s1 - c1 * c2 * s3;
        r[2][2] = c1 * c3 - c2 * s1 * s3;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_xyx(T const &t1, T const &t2, T const &t3)
    {
        T c1 = cos(t1);
        T s1 = sin(t1);
        T c2 = cos(t2);
        T s2 = sin(t2);
        T c3 = cos(t3);
        T s3 = sin(t3);

        AMAL_NMAT(4, 4) r;
        r[0][0] = c2;
        r[0][1] = s1 * s2;
        r[0][2] = -c1 * s2;
        r[0][3] = static_cast<T>(0);
        r[1][0] = s2 * s3;
        r[1][1] = c1 * c3 - c2 * s1 * s3;
        r[1][2] = c3 * s1 + c1 * c2 * s3;
        r[1][3] = static_cast<T>(0);
        r[2][0] = c3 * s2;
        r[2][1] = -c1 * s3 - c2 * c3 * s1;
        r[2][2] = c1 * c2 * c3 - s1 * s3;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_yxy(T const &t1, T const &t2, T const &t3)
    {
        T c1 = cos(t1);
        T s1 = sin(t1);
        T c2 = cos(t2);
        T s2 = sin(t2);
        T c3 = cos(t3);
        T s3 = sin(t3);

        AMAL_NMAT(4, 4) r;
        r[0][0] = c1 * c3 - c2 * s1 * s3;
        r[0][1] = s2 * s3;
        r[0][2] = -c3 * s1 - c1 * c2 * s3;
        r[0][3] = static_cast<T>(0);
        r[1][0] = s1 * s2;
        r[1][1] = c2;
        r[1][2] = c1 * s2;
        r[1][3] = static_cast<T>(0);
        r[2][0] = c1 * s3 + c2 * c3 * s1;
        r[2][1] = -c3 * s2;
        r[2][2] = c1 * c2 * c3 - s1 * s3;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_yzy(T const &t1, T const &t2, T const &t3)
    {
        T c1 = cos(t1);
        T s1 = sin(t1);
        T c2 = cos(t2);
        T s2 = sin(t2);
        T c3 = cos(t3);
        T s3 = sin(t3);

        AMAL_NMAT(4, 4) r;
        r[0][0] = c1 * c2 * c3 - s1 * s3;
        r[0][1] = c3 * s2;
        r[0][2] = -c1 * s3 - c2 * c3 * s1;
        r[0][3] = static_cast<T>(0);
        r[1][0] = -c1 * s2;
        r[1][1] = c2;
        r[1][2] = s1 * s2;
        r[1][3] = static_cast<T>(0);
        r[2][0] = c3 * s1 + c1 * c2 * s3;
        r[2][1] = s2 * s3;
        r[2][2] = c1 * c3 - c2 * s1 * s3;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_zyz(T const &t1, T const &t2, T const &t3)
    {
        T c1 = cos(t1);
        T s1 = sin(t1);
        T c2 = cos(t2);
        T s2 = sin(t2);
        T c3 = cos(t3);
        T s3 = sin(t3);

        AMAL_NMAT(4, 4) r;
        r[0][0] = c1 * c2 * c3 - s1 * s3;
        r[0][1] = c1 * s3 + c2 * c3 * s1;
        r[0][2] = -c3 * s2;
        r[0][3] = static_cast<T>(0);
        r[1][0] = -c3 * s1 - c1 * c2 * s3;
        r[1][1] = c1 * c3 - c2 * s1 * s3;
        r[1][2] = s2 * s3;
        r[1][3] = static_cast<T>(0);
        r[2][0] = c1 * s2;
        r[2][1] = s1 * s2;
        r[2][2] = c2;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_zxz(T const &t1, T const &t2, T const &t3)
    {
        T c1 = cos(t1);
        T s1 = sin(t1);
        T c2 = cos(t2);
        T s2 = sin(t2);
        T c3 = cos(t3);
        T s3 = sin(t3);

        AMAL_NMAT(4, 4) r;
        r[0][0] = c1 * c3 - c2 * s1 * s3;
        r[0][1] = c3 * s1 + c1 * c2 * s3;
        r[0][2] = s2 * s3;
        r[0][3] = static_cast<T>(0);
        r[1][0] = -c1 * s3 - c2 * c3 * s1;
        r[1][1] = c1 * c2 * c3 - s1 * s3;
        r[1][2] = c3 * s2;
        r[1][3] = static_cast<T>(0);
        r[2][0] = s1 * s2;
        r[2][1] = -c1 * s2;
        r[2][2] = c2;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_xzy(T const &t1, T const &t2, T const &t3)
    {
        T c1 = cos(t1);
        T s1 = sin(t1);
        T c2 = cos(t2);
        T s2 = sin(t2);
        T c3 = cos(t3);
        T s3 = sin(t3);

        AMAL_NMAT(4, 4) r;
        r[0][0] = c2 * c3;
        r[0][1] = s1 * s3 + c1 * c3 * s2;
        r[0][2] = c3 * s1 * s2 - c1 * s3;
        r[0][3] = static_cast<T>(0);
        r[1][0] = -s2;
        r[1][1] = c1 * c2;
        r[1][2] = c2 * s1;
        r[1][3] = static_cast<T>(0);
        r[2][0] = c2 * s3;
        r[2][1] = c1 * s2 * s3 - c3 * s1;
        r[2][2] = c1 * c3 + s1 * s2 * s3;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_yzx(T const &t1, T const &t2, T const &t3)
    {
        T c1 = cos(t1);
        T s1 = sin(t1);
        T c2 = cos(t2);
        T s2 = sin(t2);
        T c3 = cos(t3);
        T s3 = sin(t3);

        AMAL_NMAT(4, 4) r;
        r[0][0] = c1 * c2;
        r[0][1] = s2;
        r[0][2] = -c2 * s1;
        r[0][3] = static_cast<T>(0);
        r[1][0] = s1 * s3 - c1 * c3 * s2;
        r[1][1] = c2 * c3;
        r[1][2] = c1 * s3 + c3 * s1 * s2;
        r[1][3] = static_cast<T>(0);
        r[2][0] = c3 * s1 + c1 * s2 * s3;
        r[2][1] = -c2 * s3;
        r[2][2] = c1 * c3 - s1 * s2 * s3;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_ZYX(T const &t1, T const &t2, T const &t3)
    {
        T c1 = cos(t1);
        T s1 = sin(t1);
        T c2 = cos(t2);
        T s2 = sin(t2);
        T c3 = cos(t3);
        T s3 = sin(t3);

        AMAL_NMAT(4, 4) r;
        r[0][0] = c1 * c2;
        r[0][1] = c2 * s1;
        r[0][2] = -s2;
        r[0][3] = static_cast<T>(0);
        r[1][0] = c1 * s2 * s3 - c3 * s1;
        r[1][1] = c1 * c3 + s1 * s2 * s3;
        r[1][2] = c2 * s3;
        r[1][3] = static_cast<T>(0);
        r[2][0] = s1 * s3 + c1 * c3 * s2;
        r[2][1] = c3 * s1 * s2 - c1 * s3;
        r[2][2] = c2 * c3;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline AMAL_NMAT(4, 4) euler_angle_zxy(T const &t1, T const &t2, T const &t3)
    {
        T c1 = cos(t1);
        T s1 = sin(t1);
        T c2 = cos(t2);
        T s2 = sin(t2);
        T c3 = cos(t3);
        T s3 = sin(t3);

        AMAL_NMAT(4, 4) r;
        r[0][0] = c1 * c3 - s1 * s2 * s3;
        r[0][1] = c3 * s1 + c1 * s2 * s3;
        r[0][2] = -c2 * s3;
        r[0][3] = static_cast<T>(0);
        r[1][0] = -c2 * s1;
        r[1][1] = c1 * c2;
        r[1][2] = s2;
        r[1][3] = static_cast<T>(0);
        r[2][0] = c1 * s3 + c3 * s1 * s2;
        r[2][1] = s1 * s3 - c1 * c3 * s2;
        r[2][2] = c2 * c3;
        r[2][3] = static_cast<T>(0);
        r[3][0] = static_cast<T>(0);
        r[3][1] = static_cast<T>(0);
        r[3][2] = static_cast<T>(0);
        r[3][3] = static_cast<T>(1);
        return r;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_xyz(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(M[2][1], M[2][2]);
        T C2 = sqrt(M[0][0] * M[0][0] + M[1][0] * M[1][0]);
        T T2 = atan2(-M[2][0], C2);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(S1 * M[0][2] - C1 * M[0][1], C1 * M[1][1] - S1 * M[1][2]);
        t1 = -T1;
        t2 = -T2;
        t3 = -T3;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_yxz(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(M[2][0], M[2][2]);
        T C2 = sqrt(M[0][1] * M[0][1] + M[1][1] * M[1][1]);
        T T2 = atan2(-M[2][1], C2);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(S1 * M[1][2] - C1 * M[1][0], C1 * M[0][0] - S1 * M[0][2]);
        t1 = T1;
        t2 = T2;
        t3 = T3;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_xzx(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(M[0][2], M[0][1]);
        T S2 = sqrt(M[1][0] * M[1][0] + M[2][0] * M[2][0]);
        T T2 = atan2(S2, M[0][0]);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(C1 * M[1][2] - S1 * M[1][1], C1 * M[2][2] - S1 * M[2][1]);
        t1 = T1;
        t2 = T2;
        t3 = T3;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_xyx(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(M[0][1], -M[0][2]);
        T S2 = sqrt(M[1][0] * M[1][0] + M[2][0] * M[2][0]);
        T T2 = atan2(S2, M[0][0]);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(-C1 * M[2][1] - S1 * M[2][2], C1 * M[1][1] + S1 * M[1][2]);
        t1 = T1;
        t2 = T2;
        t3 = T3;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_yxy(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(M[1][0], M[1][2]);
        T S2 = sqrt(M[0][1] * M[0][1] + M[2][1] * M[2][1]);
        T T2 = atan2(S2, M[1][1]);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(C1 * M[2][0] - S1 * M[2][2], C1 * M[0][0] - S1 * M[0][2]);
        t1 = T1;
        t2 = T2;
        t3 = T3;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_yzy(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(M[1][2], -M[1][0]);
        T S2 = sqrt(M[0][1] * M[0][1] + M[2][1] * M[2][1]);
        T T2 = atan2(S2, M[1][1]);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(-S1 * M[0][0] - C1 * M[0][2], S1 * M[2][0] + C1 * M[2][2]);
        t1 = T1;
        t2 = T2;
        t3 = T3;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_zyz(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(M[2][1], M[2][0]);
        T S2 = sqrt(M[0][2] * M[0][2] + M[1][2] * M[1][2]);
        T T2 = atan2(S2, M[2][2]);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(C1 * M[0][1] - S1 * M[0][0], C1 * M[1][1] - S1 * M[1][0]);
        t1 = T1;
        t2 = T2;
        t3 = T3;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_zxz(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(M[2][0], -M[2][1]);
        T S2 = sqrt(M[0][2] * M[0][2] + M[1][2] * M[1][2]);
        T T2 = atan2(S2, M[2][2]);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(-C1 * M[1][0] - S1 * M[1][1], C1 * M[0][0] + S1 * M[0][1]);
        t1 = T1;
        t2 = T2;
        t3 = T3;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_xzy(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(M[1][2], M[1][1]);
        T C2 = sqrt(M[0][0] * M[0][0] + M[2][0] * M[2][0]);
        T T2 = atan2(-M[1][0], C2);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(S1 * M[0][1] - C1 * M[0][2], C1 * M[2][2] - S1 * M[2][1]);
        t1 = T1;
        t2 = T2;
        t3 = T3;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_yzx(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(-M[0][2], M[0][0]);
        T C2 = sqrt(M[1][1] * M[1][1] + M[2][1] * M[2][1]);
        T T2 = atan2(M[0][1], C2);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(S1 * M[1][0] + C1 * M[1][2], S1 * M[2][0] + C1 * M[2][2]);
        t1 = T1;
        t2 = T2;
        t3 = T3;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_zyx(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(M[0][1], M[0][0]);
        T C2 = sqrt(M[1][2] * M[1][2] + M[2][2] * M[2][2]);
        T T2 = atan2(-M[0][2], C2);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(S1 * M[2][0] - C1 * M[2][1], C1 * M[1][1] - S1 * M[1][0]);
        t1 = T1;
        t2 = T2;
        t3 = T3;
    }

    template <typename T, bool aligned>
    inline void extract_euler_angle_zxy(AMAL_NMAT(4, 4) const &M, T &t1, T &t2, T &t3)
    {
        T T1 = atan2(-M[1][0], M[1][1]);
        T C2 = sqrt(M[0][2] * M[0][2] + M[2][2] * M[2][2]);
        T T2 = atan2(M[1][2], C2);
        T S1 = sin(T1);
        T C1 = cos(T1);
        T T3 = atan2(C1 * M[2][0] + S1 * M[2][1], C1 * M[0][0] + S1 * M[0][1]);
        t1 = T1;
        t2 = T2;
        t3 = T3;
    }
} // namespace amal