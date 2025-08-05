#pragma once

#include "matrix.hpp"
#include "vector.hpp"

namespace amal
{
    template <typename T, bool aligned>
    inline AMAL_NVEC(2) screen_to_ndc(const AMAL_NVEC(2) & screen, int x, int y)
    {
        T x_ndc = (screen.x / x) * static_cast<T>(2) - static_cast<T>(1);
        T y_ndc = 1.0f - (screen.y / y) * static_cast<T>(2);
        return AMAL_NVEC(2)(x_ndc, y_ndc);
    }

    template <typename T, typename U, bool aligned>
    inline AMAL_NVEC(3) project(AMAL_NVEC(3) const &obj, AMAL_NMAT(4, 4) const &model, AMAL_NMAT(4, 4) const &proj,
                                AMAL_NVEC(4) const &viewport)
    {
        AMAL_NVEC(4) tmp(obj, static_cast<T>(1));
        tmp = model * tmp;
        tmp = proj * tmp;

        tmp /= tmp.w;
        tmp = tmp * static_cast<T>(0.5) + static_cast<T>(0.5);
#ifdef AMAL_FMA_ENABLE
        tmp[0] = fma(tmp[0], T(viewport[2]), T(viewport[0]));
        tmp[1] = fma(tmp[1], T(viewport[3]), T(viewport[1]));
#else
        tmp[0] = tmp[0] * T(viewport[2]) + T(viewport[0]);
        tmp[1] = tmp[1] * T(viewport[3]) + T(viewport[1]);
#endif

        return AMAL_NVEC(3)(tmp);
    }

    template <typename T, typename U, bool aligned>
    inline AMAL_NVEC(3) unproject(AMAL_NVEC(3) const &win, AMAL_NMAT(4, 4) const &model, AMAL_NMAT(4, 4) const &proj,
                                  AMAL_NVEC(4) const &viewport)
    {
        AMAL_NMAT(4, 4) inverse = inverse(proj * model);

        AMAL_NVEC(4) tmp(win, T(1));
        tmp.x = (tmp.x - T(viewport[0])) / T(viewport[2]);
        tmp.y = (tmp.y - T(viewport[1])) / T(viewport[3]);
        tmp = tmp * static_cast<T>(2) - static_cast<T>(1);

        AMAL_NVEC(4) obj = inverse * tmp;
        obj /= obj.w;

        return AMAL_NVEC(3)(obj);
    }
} // namespace amal