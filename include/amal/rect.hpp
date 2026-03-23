#pragma once

#include "vector.hpp"

namespace amal
{
    namespace internal
    {
        template <typename T>
        struct rect_base
        {
            using value_type = T::value_type;

            T offset;
            T size;

            rect_base(T offset, T size) : offset(offset), size(size) {}
            rect_base(value_type x, value_type y, value_type w, value_type h) : offset(x, y), size(w, h) {}
            rect_base(value_type offset, value_type size) : offset(offset), size(size) {}
            rect_base() : offset(static_cast<value_type>(0)), size(static_cast<value_type>(0)) {}
        };
    } // namespace internal

    using rect = internal::rect_base<vec2>;
    using irect = internal::rect_base<ivec2>;

    template <typename R>
    inline R::value_type get_rect_left(const R &r)
    {
        return r.offset.x;
    }

    template <typename R>
    inline R::value_type get_rect_top(const R &r)
    {
        return r.offset.y;
    }

    template <typename R>
    inline R::value_type get_rect_right(const R &r)
    {
        return get_rect_left(r) + r.size.x;
    }

    template <typename R>
    inline R::value_type get_rect_bottom(const R &r)
    {
        return get_rect_top(r) + r.size.y;
    }

    template <typename R>
    inline R::value_type get_rect_area(const R &r)
    {
        return r.size.x * r.size.y;
    }

    template <typename R>
    inline bool is_rects_overlap(const R &a, const R &b)
    {
        return get_rect_left(a) < get_rect_right(b) && get_rect_right(a) > get_rect_left(b) &&
               get_rect_top(a) < get_rect_bottom(b) && get_rect_bottom(a) > get_rect_top(b);
    }

    template <typename R>
    inline bool is_rect_empty(const R &rect)
    {
        return rect.size.x <= 0 || rect.size.y <= 0;
    }

    template <typename R>
    inline bool is_rect_contains(const R &outer, const R &inner)
    {
        return get_rect_left(inner) >= get_rect_left(outer) && get_rect_top(inner) >= get_rect_top(outer) &&
               get_rect_right(inner) <= get_rect_right(outer) && get_rect_bottom(inner) <= get_rect_bottom(outer);
    }

} // namespace amal