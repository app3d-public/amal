#pragma once

#include "vector.hpp"

namespace amal
{
    template <typename T, enum Pack P>
    inline AMAL_NVEC(2) screen_to_ndc(const AMAL_NVEC(2) & screen, int x, int y)
    {
        T x_ndc = (screen.x / x) * static_cast<T>(2) - static_cast<T>(1);
        T y_ndc = 1.0f - (screen.y / y) * static_cast<T>(2);
        return AMAL_NVEC(2)(x_ndc, y_ndc);
    }
} // namespace amal