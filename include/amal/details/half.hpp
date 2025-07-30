#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>

#if defined(__F16C__)
    #include <immintrin.h>
#endif

#define AMAL_HALF_OVERFLOW(x)  x | 0x7C00
#define AMAL_HALF_UNDERFLOW(x) x
#define AMAL_HALF_INVALID      0x7FFF

namespace amal
{
    namespace details
    {
        typedef std::uint_least16_t uint16;
        typedef std::uint_fast32_t uint32;
        typedef std::int_fast32_t int32;

        inline constexpr unsigned int round(unsigned int value, int g, int s) { return value + (g & (s | value)); }

        template <typename T>
        unsigned int float_to_half(T value)
        {
            unsigned int hbits = static_cast<unsigned>(builtin_signbit(value)) << 15;
            if (value == T()) return hbits;
            if (builtin_isnan(value)) return hbits | AMAL_HALF_INVALID;
            if (builtin_isinf(value)) return hbits | 0x7C00;
            int exp;
            std::frexp(value, &exp);
            if (exp > 15) return AMAL_HALF_OVERFLOW(hbits);
            if (exp < -13)
                value = std::ldexp(value, 25);
            else
            {
                value = std::ldexp(value, 12 - exp);
                hbits |= ((exp + 13) << 10);
            }
            T ival, frac = std::modf(value, &ival);
            int m = std::abs(static_cast<int>(ival));
            return round(hbits + (m >> 1), m & 1, frac != T());
        }

        template <>
        inline unsigned int float_to_half<float>(float value)
        {
#if defined(__F16C__)
            return _mm_cvtsi128_si32(_mm_cvtps_ph(_mm_set_ss(value), _MM_FROUND_TO_NEAREST_INT));
#else
            static_assert(sizeof(float) == sizeof(std::uint32_t), "float must be 32-bit IEEE 754");
            uint32_t fbits;
            memcpy(&fbits, &value, sizeof(float));

            unsigned int sign = (fbits >> 16) & 0x8000;
            fbits &= 0x7FFFFFFF;
            if (fbits >= 0x7F800000)
                return sign | 0x7C00 | ((fbits > 0x7F800000) ? (0x200 | ((fbits >> 13) & 0x3FF)) : 0);
            if (fbits >= 0x47800000) return AMAL_HALF_OVERFLOW(sign);
            if (fbits >= 0x38800000)
                return round(sign | (((fbits >> 23) - 112) << 10) | ((fbits >> 13) & 0x3FF), (fbits >> 12) & 1,
                             (fbits & 0xFFF) != 0);
            if (fbits >= 0x33000000)
            {
                int i = 125 - (fbits >> 23);
                fbits = (fbits & 0x7FFFFF) | 0x800000;
                return round(sign | (fbits >> (i + 1)), (fbits >> i) & 1,
                             (fbits & ((static_cast<uint32>(1) << i) - 1)) != 0);
            }
            if (fbits != 0) return AMAL_HALF_UNDERFLOW(sign);
            return sign;
#endif
        }
        inline float half_to_float(unsigned int value)
        {
#if defined(__F16C__)
            return _mm_cvtss_f32(_mm_cvtph_ps(_mm_cvtsi32_si128(value)));
#else
            const uint32 sign = uint32(value & 0x8000) << 16;
            uint32 exp_h = (value & 0x7C00) >> 10;
            uint32 mant_h = value & 0x03FF;
            uint32 fbits;

            if (exp_h == 0) // zero / denorm
            {
                if (mant_h == 0)
                    fbits = sign; // ±0
                else
                {
                    while ((mant_h & 0x0400) == 0)
                    {
                        mant_h <<= 1;
                        --exp_h;
                    }
                    mant_h &= 0x03FF;
                    ++exp_h;
                    fbits = sign | ((exp_h + 112) << 23) | (mant_h << 13);
                }
            }
            else if (exp_h == 0x1F) // Inf / NaN
                fbits = sign | 0x7F800000 | (mant_h << 13);
            else
                fbits = sign | ((exp_h + 112) << 23) | (mant_h << 13);

            float out;
            memcpy(&out, &fbits, sizeof(out));
            return out;
#endif
        }

        inline constexpr unsigned int signal(unsigned int x, unsigned int y)
        {
            return ((x & 0x7FFF) > 0x7C00) ? (x | 0x200) : (y | 0x200);
        }

        inline constexpr bool compsignal(uint16 x, uint16 y) { return (x & 0x7FFF) > 0x7C00 || (y & 0x7FFF) > 0x7C00; }

        template <bool Q, bool R>
        unsigned int mod(unsigned int x, unsigned int y, int *quo = NULL)
        {
            unsigned int q = 0;
            if (x > y)
            {
                int absx = x, absy = y, expx = 0, expy = 0;
                for (; absx < 0x400; absx <<= 1, --expx);
                for (; absy < 0x400; absy <<= 1, --expy);
                expx += absx >> 10;
                expy += absy >> 10;
                int mx = (absx & 0x3FF) | 0x400, my = (absy & 0x3FF) | 0x400;
                for (int d = expx - expy; d; --d)
                {
                    if (!Q && mx == my) return 0;
                    if (mx >= my)
                    {
                        mx -= my;
                        q += Q;
                    }
                    mx <<= 1;
                    q <<= static_cast<int>(Q);
                }
                if (!Q && mx == my) return 0;
                if (mx >= my)
                {
                    mx -= my;
                    ++q;
                }
                if (Q)
                {
                    q &= (1 << (std::numeric_limits<int>::digits - 1)) - 1;
                    if (!mx) return *quo = q, 0;
                }
                for (; mx < 0x400; mx <<= 1, --expy);
                x = (expy > 0) ? ((expy << 10) | (mx & 0x3FF)) : (mx >> (1 - expy));
            }
            if (R)
            {
                unsigned int a, b;
                if (y < 0x800)
                {
                    a = (x < 0x400) ? (x << 1) : (x + 0x400);
                    b = y;
                }
                else
                {
                    a = x;
                    b = y - 0x400;
                }
                if (a > b || (a == b && (q & 1)))
                {
                    int exp = (y >> 10) + (y <= 0x3FF), d = exp - (x >> 10) - (x <= 0x3FF);
                    int m =
                        (((y & 0x3FF) | ((y > 0x3FF) << 10)) << 1) - (((x & 0x3FF) | ((x > 0x3FF) << 10)) << (1 - d));
                    for (; m < 0x800 && exp > 1; m <<= 1, --exp);
                    x = 0x8000 + ((exp - 1) << 10) + (m >> 1);
                    q += Q;
                }
            }
            if (Q) *quo = q;
            return x;
        }

        inline uint32 sign_mask(uint32 arg)
        {
            static const int N = std::numeric_limits<uint32>::digits - 1;
            return static_cast<int32>(arg) >> N;
        }

        inline uint32 arithmetic_shift(uint32 arg, int i) { return static_cast<int32>(arg) >> i; }
    } // namespace details
} // namespace amal