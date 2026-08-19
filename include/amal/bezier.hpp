#pragma once

#include <acul/vector.hpp>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include "common.hpp"
#include "geometric.hpp"

namespace amal
{
    // Incoming handles control the segment ending at a knot; outgoing handles
    // control the segment starting at it.
    enum class bezier_handle : std::uint8_t
    {
        in,
        out
    };

    enum class bezier_handle_mode : std::uint8_t
    {
        // Editing either handle leaves the other one untouched.
        independent,
        // Handles stay collinear while retaining independent lengths.
        aligned,
        // Handles stay collinear with equal lengths.
        mirrored
    };

    enum class bezier_interpolation : std::uint8_t
    {
        linear,
        smooth
    };

    template <typename Point>
    struct cubic_bezier
    {
        Point p0{};
        Point p1{};
        Point p2{};
        Point p3{};
    };

    template <typename Point>
    struct cubic_bezier_split
    {
        cubic_bezier<Point> left{};
        cubic_bezier<Point> right{};
    };

    template <typename Point, typename Scalar>
    struct bezier_sample
    {
        Point position{};
        Scalar t = Scalar(0);
    };

    template <typename Point, typename Scalar>
    struct bezier_arc_length_table
    {
        struct entry
        {
            Point position{};
            std::size_t segment_index = 0u;
            Scalar t = Scalar(0);
            Scalar distance = Scalar(0);
        };

        using point_type = Point;
        using scalar_type = Scalar;
        using entry_type = entry;

        acul::vector<entry> entries;
        Scalar total_length = Scalar(0);

        bool empty() const { return entries.empty(); }
        void clear()
        {
            entries.clear();
            total_length = Scalar(0);
        }
    };

    using bezier_arc_length_table2 = bezier_arc_length_table<vec2, AMAL_FLOAT_TYPE>;
    using bezier_arc_length_table3 = bezier_arc_length_table<vec3, AMAL_FLOAT_TYPE>;

    template <typename Point>
    struct bezier_frame
    {
        Point position{};
        Point tangent{};
        Point normal{};
        Point binormal{};
    };

    template <typename Point>
    struct bezier_point
    {
        Point position{};
        Point handle_in{};
        Point handle_out{};
        bezier_handle_mode handle_mode = bezier_handle_mode::mirrored;

        bezier_point() = default;
        explicit bezier_point(const Point &value) : position(value), handle_in(value), handle_out(value) {}
        bezier_point(const Point &value, const Point &in, const Point &out,
                     bezier_handle_mode mode = bezier_handle_mode::mirrored)
            : position(value), handle_in(in), handle_out(out), handle_mode(mode)
        {
        }
    };

    template <typename Point>
    struct bezier_curve
    {
        using point_type = Point;
        using knot_type = bezier_point<Point>;

        acul::vector<knot_type> points;
        bool closed = false;

        bool empty() const { return points.empty(); }
        std::size_t size() const { return points.size(); }
        std::size_t segment_count() const
        {
            if (points.size() < 2u) return 0u;
            return closed ? points.size() : points.size() - 1u;
        }
    };

    using bezier_curve2 = bezier_curve<vec2>;
    using bezier_point2 = bezier_point<vec2>;

    using bezier_curve3 = bezier_curve<vec3>;
    using bezier_point3 = bezier_point<vec3>;

    namespace detail
    {
        template <typename Point, bool = is_vector_v<Point>>
        struct point_scalar
        {
            using type = Point;
        };

        template <typename Point>
        struct point_scalar<Point, true>
        {
            using type = typename Point::value_type;
        };

        template <typename Point>
        using point_scalar_t = typename point_scalar<Point>::type;

        template <typename Point, typename Scalar>
        struct bezier_flatten_node
        {
            cubic_bezier<Point> curve{};
            Scalar t0 = Scalar(0);
            Scalar t1 = Scalar(1);
            std::uint32_t depth = 0u;
        };

        template <typename Point>
        struct bezier_coefficients
        {
            Point a{};
            Point b{};
            Point c{};
            Point d{};
        };

        template <typename Point>
        inline bezier_coefficients<Point> make_bezier_coefficients(const cubic_bezier<Point> &curve)
        {
            using scalar = point_scalar_t<Point>;
            return {curve.p3 - curve.p2 * scalar(3) + curve.p1 * scalar(3) - curve.p0,
                    curve.p2 * scalar(3) - curve.p1 * scalar(6) + curve.p0 * scalar(3),
                    (curve.p1 - curve.p0) * scalar(3), curve.p0};
        }

        template <typename Point, typename Scalar>
        inline Point bezier_evaluate_unchecked(const bezier_coefficients<Point> &coefficients, Scalar t)
        {
            using scalar = point_scalar_t<Point>;
            const scalar value = static_cast<scalar>(t);
            return ((coefficients.a * value + coefficients.b) * value + coefficients.c) * value + coefficients.d;
        }

        template <typename Point, typename Scalar>
        inline Point bezier_derivative_unchecked(const bezier_coefficients<Point> &coefficients, Scalar t)
        {
            using scalar = point_scalar_t<Point>;
            const scalar value = static_cast<scalar>(t);
            return (coefficients.a * (scalar(3) * value) + coefficients.b * scalar(2)) * value + coefficients.c;
        }

        template <typename Point, typename Scalar>
        inline Point bezier_second_derivative_unchecked(const bezier_coefficients<Point> &coefficients, Scalar t)
        {
            using scalar = point_scalar_t<Point>;
            return coefficients.a * (scalar(6) * static_cast<scalar>(t)) + coefficients.b * scalar(2);
        }

        template <typename Point>
        inline point_scalar_t<Point> point_dot(const Point &a, const Point &b)
        {
            if constexpr (is_vector_v<Point>) return dot(a, b);
            else return a * b;
        }

        template <typename Point>
        inline point_scalar_t<Point> point_distance_squared(const Point &a, const Point &b)
        {
            const Point delta = a - b;
            return point_dot(delta, delta);
        }

        template <typename Point>
        inline point_scalar_t<Point> point_segment_distance_squared(const Point &point, const Point &a, const Point &b)
        {
            using scalar = point_scalar_t<Point>;
            const Point ab = b - a;
            const scalar denominator = point_dot(ab, ab);
            if (denominator <= std::numeric_limits<scalar>::epsilon()) return point_distance_squared(point, a);
            const scalar t = amal::clamp(point_dot(point - a, ab) / denominator, scalar(0), scalar(1));
            return point_distance_squared(point, a + ab * t);
        }

        template <typename Point>
        inline Point safe_normalize(const Point &value, const Point &fallback = Point{})
        {
            using scalar = point_scalar_t<Point>;
            const scalar value_length = length(value);
            if (value_length <= std::numeric_limits<scalar>::epsilon()) return fallback;
            return value / value_length;
        }

        template <typename Point>
        inline Point opposite_handle_position(const bezier_point<Point> &point, bezier_handle edited,
                                              const Point &value, bezier_handle_mode mode)
        {
            using scalar = point_scalar_t<Point>;
            if (mode == bezier_handle_mode::mirrored) return point.position * scalar(2) - value;

            const Point &opposite = edited == bezier_handle::in ? point.handle_out : point.handle_in;
            const scalar opposite_length = distance(point.position, opposite);
            const Point direction = point.position - value;
            const scalar edited_length = length(direction);
            if (edited_length <= std::numeric_limits<scalar>::epsilon()) return opposite;
            return point.position + direction * (opposite_length / edited_length);
        }

        template <typename Curve>
        inline void append_point(Curve &curve, const typename Curve::knot_type &point)
        {
            curve.points.push_back(point);
        }

        struct identity_projection
        {
            template <typename Point>
            const Point &operator()(const Point &point) const
            {
                return point;
            }
        };
    } // namespace detail

    template <typename Point>
    using bezier_scalar_t = detail::point_scalar_t<Point>;

    template <typename Point>
    using bezier_flatten_workspace = acul::vector<detail::bezier_flatten_node<Point, bezier_scalar_t<Point>>>;

    using bezier_flatten_workspace2 = bezier_flatten_workspace<vec2>;
    using bezier_flatten_workspace3 = bezier_flatten_workspace<vec3>;

    template <typename Point>
    struct bezier_curve_location
    {
        using scalar_type = bezier_scalar_t<Point>;

        std::size_t segment_index = 0u;
        scalar_type t = scalar_type(0);
        scalar_type distance_squared = std::numeric_limits<scalar_type>::max();
        bool valid = false;
    };

    template <typename Point, typename Scalar>
    inline Point bezier_evaluate_unchecked(const cubic_bezier<Point> &curve, Scalar t)
    {
        return detail::bezier_evaluate_unchecked(detail::make_bezier_coefficients(curve), t);
    }

    template <typename Point, typename Scalar>
    inline Point bezier_derivative_unchecked(const cubic_bezier<Point> &curve, Scalar t)
    {
        return detail::bezier_derivative_unchecked(detail::make_bezier_coefficients(curve), t);
    }

    template <typename Point, typename Scalar>
    inline Point bezier_second_derivative_unchecked(const cubic_bezier<Point> &curve, Scalar t)
    {
        return detail::bezier_second_derivative_unchecked(detail::make_bezier_coefficients(curve), t);
    }

    template <typename Point, typename Scalar>
    inline Point bezier_evaluate(const cubic_bezier<Point> &curve, Scalar t)
    {
        using scalar = bezier_scalar_t<Point>;
        const scalar value = clamp(static_cast<scalar>(t), scalar(0), scalar(1));
        return detail::bezier_evaluate_unchecked(detail::make_bezier_coefficients(curve), value);
    }

    template <typename Point, typename Scalar>
    inline Point bezier_derivative(const cubic_bezier<Point> &curve, Scalar t)
    {
        using scalar = bezier_scalar_t<Point>;
        const scalar value = clamp(static_cast<scalar>(t), scalar(0), scalar(1));
        return detail::bezier_derivative_unchecked(detail::make_bezier_coefficients(curve), value);
    }

    template <typename Point, typename Scalar>
    inline Point bezier_second_derivative(const cubic_bezier<Point> &curve, Scalar t)
    {
        using scalar = bezier_scalar_t<Point>;
        const scalar value = clamp(static_cast<scalar>(t), scalar(0), scalar(1));
        return detail::bezier_second_derivative_unchecked(detail::make_bezier_coefficients(curve), value);
    }

    template <typename Point, typename Scalar>
    inline cubic_bezier_split<Point> bezier_split(const cubic_bezier<Point> &curve, Scalar t)
    {
        using scalar = bezier_scalar_t<Point>;
        const scalar value = clamp(static_cast<scalar>(t), scalar(0), scalar(1));
        const Point p01 = mix(curve.p0, curve.p1, value);
        const Point p12 = mix(curve.p1, curve.p2, value);
        const Point p23 = mix(curve.p2, curve.p3, value);
        const Point p012 = mix(p01, p12, value);
        const Point p123 = mix(p12, p23, value);
        const Point point = mix(p012, p123, value);
        return {{curve.p0, p01, p012, point}, {point, p123, p23, curve.p3}};
    }

    namespace detail
    {
        template <typename Point, typename Workspace, typename Projection, typename Tolerance, typename Visitor>
        inline void flatten_projected_visit(const cubic_bezier<Point> &curve, Workspace &workspace,
                                            Projection projection, Tolerance tolerance, std::uint32_t max_depth,
                                            bool include_start, Visitor visitor)
        {
            using scalar = point_scalar_t<Point>;
            using projected_point = std::decay_t<decltype(projection(curve.p0))>;
            using projected_scalar = point_scalar_t<projected_point>;
            const projected_scalar safe_tolerance =
                amal::max(static_cast<projected_scalar>(tolerance), std::numeric_limits<projected_scalar>::epsilon());
            const projected_scalar tolerance_squared = safe_tolerance * safe_tolerance;

            workspace.clear();
            workspace.push_back(bezier_flatten_node<Point, scalar>{curve, scalar(0), scalar(1), 0u});
            if (include_start) visitor(curve.p0, scalar(0));

            while (!workspace.empty())
            {
                const auto node = workspace.back();
                workspace.pop_back();
                const projected_point p0 = projection(node.curve.p0);
                const projected_point p1 = projection(node.curve.p1);
                const projected_point p2 = projection(node.curve.p2);
                const projected_point p3 = projection(node.curve.p3);
                const projected_scalar d1 = point_segment_distance_squared(p1, p0, p3);
                const projected_scalar d2 = point_segment_distance_squared(p2, p0, p3);
                if (node.depth >= max_depth || amal::max(d1, d2) <= tolerance_squared)
                {
                    visitor(node.curve.p3, node.t1);
                    continue;
                }

                const auto halves = bezier_split(node.curve, scalar(0.5));
                const scalar middle_t = (node.t0 + node.t1) * scalar(0.5);
                const std::uint32_t next_depth = node.depth + 1u;
                workspace.push_back(bezier_flatten_node<Point, scalar>{halves.right, middle_t, node.t1, next_depth});
                workspace.push_back(bezier_flatten_node<Point, scalar>{halves.left, node.t0, middle_t, next_depth});
            }
        }
    } // namespace detail

    template <typename Point, typename OutputIt, typename Workspace, typename Projection, typename Tolerance>
    inline OutputIt bezier_flatten_projected_with_workspace(const cubic_bezier<Point> &curve, OutputIt output,
                                                            Workspace &workspace, Projection projection,
                                                            Tolerance tolerance, std::uint32_t max_depth = 12u,
                                                            bool include_start = true)
    {
        detail::flatten_projected_visit(curve, workspace, projection, tolerance, max_depth, include_start,
                                        [&](const Point &point, auto) { *output++ = point; });
        return output;
    }

    template <typename Point, typename OutputIt, typename Workspace, typename Tolerance>
    inline OutputIt bezier_flatten_with_workspace(const cubic_bezier<Point> &curve, OutputIt output,
                                                  Workspace &workspace, Tolerance tolerance,
                                                  std::uint32_t max_depth = 12u, bool include_start = true)
    {
        return bezier_flatten_projected_with_workspace(curve, output, workspace, detail::identity_projection{},
                                                       tolerance, max_depth, include_start);
    }

    template <typename Point, typename OutputIt, typename Workspace, typename Projection, typename Tolerance>
    inline OutputIt bezier_flatten_samples_projected_with_workspace(const cubic_bezier<Point> &curve, OutputIt output,
                                                                    Workspace &workspace, Projection projection,
                                                                    Tolerance tolerance, std::uint32_t max_depth = 12u,
                                                                    bool include_start = true)
    {
        using scalar = bezier_scalar_t<Point>;
        detail::flatten_projected_visit(
            curve, workspace, projection, tolerance, max_depth, include_start,
            [&](const Point &point, scalar t) { *output++ = bezier_sample<Point, scalar>{point, t}; });
        return output;
    }

    template <typename Point>
    inline bezier_scalar_t<Point> bezier_closest_parameter(const cubic_bezier<Point> &curve, const Point &point,
                                                           std::uint32_t coarse_steps = 32u,
                                                           std::uint32_t refine_steps = 8u)
    {
        using scalar = bezier_scalar_t<Point>;
        coarse_steps = coarse_steps < 2u ? 2u : coarse_steps;
        const auto coefficients = detail::make_bezier_coefficients(curve);
        scalar best_t = scalar(0);
        scalar best_distance = detail::point_distance_squared(curve.p0, point);
        for (std::uint32_t i = 1u; i <= coarse_steps; ++i)
        {
            const scalar t = static_cast<scalar>(i) / static_cast<scalar>(coarse_steps);
            const scalar candidate =
                detail::point_distance_squared(detail::bezier_evaluate_unchecked(coefficients, t), point);
            if (candidate < best_distance)
            {
                best_distance = candidate;
                best_t = t;
            }
        }

        for (std::uint32_t i = 0u; i < refine_steps; ++i)
        {
            const Point delta = detail::bezier_evaluate_unchecked(coefficients, best_t) - point;
            const Point first = detail::bezier_derivative_unchecked(coefficients, best_t);
            const scalar denominator =
                detail::point_dot(first, first) +
                detail::point_dot(delta, detail::bezier_second_derivative_unchecked(coefficients, best_t));
            if (abs(denominator) <= std::numeric_limits<scalar>::epsilon()) break;
            const scalar next = clamp(best_t - detail::point_dot(delta, first) / denominator, scalar(0), scalar(1));
            if (abs(next - best_t) <= std::numeric_limits<scalar>::epsilon()) break;
            best_t = next;
        }
        return best_t;
    }

    template <typename Point>
    inline Point bezier_closest_point(const cubic_bezier<Point> &curve, const Point &point,
                                      std::uint32_t coarse_steps = 32u, std::uint32_t refine_steps = 8u)
    {
        return bezier_evaluate(curve, bezier_closest_parameter(curve, point, coarse_steps, refine_steps));
    }

    template <typename Point, typename OutputIt>
    inline OutputIt bezier_flatten(const cubic_bezier<Point> &curve, OutputIt output,
                                   bezier_scalar_t<Point> tolerance = bezier_scalar_t<Point>(0.5),
                                   std::uint32_t max_depth = 12u, bool include_start = true)
    {
        using scalar = bezier_scalar_t<Point>;
        const scalar safe_tolerance = max(tolerance, std::numeric_limits<scalar>::epsilon());
        const scalar tolerance_squared = safe_tolerance * safe_tolerance;
        if (include_start) *output++ = curve.p0;

        auto recurse = [&](auto &&self, const cubic_bezier<Point> &part, std::uint32_t depth) -> void {
            const scalar d1 = detail::point_segment_distance_squared(part.p1, part.p0, part.p3);
            const scalar d2 = detail::point_segment_distance_squared(part.p2, part.p0, part.p3);
            if (depth >= max_depth || max(d1, d2) <= tolerance_squared)
            {
                *output++ = part.p3;
                return;
            }
            const auto halves = bezier_split(part, scalar(0.5));
            self(self, halves.left, depth + 1u);
            self(self, halves.right, depth + 1u);
        };
        recurse(recurse, curve, 0u);
        return output;
    }

    template <typename Point>
    inline bezier_scalar_t<Point> bezier_approximate_length(const cubic_bezier<Point> &curve, std::uint32_t steps = 32u)
    {
        using scalar = bezier_scalar_t<Point>;
        steps = steps == 0u ? 1u : steps;
        const auto coefficients = detail::make_bezier_coefficients(curve);
        scalar result = scalar(0);
        Point previous = curve.p0;
        for (std::uint32_t i = 1u; i <= steps; ++i)
        {
            const Point current =
                detail::bezier_evaluate_unchecked(coefficients, static_cast<scalar>(i) / static_cast<scalar>(steps));
            result += distance(previous, current);
            previous = current;
        }
        return result;
    }

    template <typename Curve>
    inline typename Curve::knot_type &curve_point(Curve &curve, std::size_t index)
    {
        return curve.points[index];
    }

    template <typename Curve>
    inline const typename Curve::knot_type &curve_point(const Curve &curve, std::size_t index)
    {
        return curve.points[index];
    }

    template <typename Curve>
    inline cubic_bezier<typename Curve::point_type> curve_segment(const Curve &curve, std::size_t index)
    {
        const std::size_t next = (index + 1u) % curve.points.size();
        const auto &a = curve.points[index];
        const auto &b = curve.points[next];
        return {a.position, a.handle_out, b.handle_in, b.position};
    }

    template <typename Curve, typename Scalar>
    inline typename Curve::point_type curve_evaluate(const Curve &curve, std::size_t segment_index, Scalar t)
    {
        return bezier_evaluate(curve_segment(curve, segment_index), t);
    }

    template <typename Curve>
    inline bezier_curve_location<typename Curve::point_type>
    curve_closest_location(const Curve &curve, const typename Curve::point_type &point,
                           std::uint32_t coarse_steps = 32u, std::uint32_t refine_steps = 8u)
    {
        using point_type = typename Curve::point_type;
        bezier_curve_location<point_type> result{};
        for (std::size_t i = 0u; i < curve.segment_count(); ++i)
        {
            const auto segment = curve_segment(curve, i);
            const auto t = bezier_closest_parameter(segment, point, coarse_steps, refine_steps);
            const auto distance_squared = detail::point_distance_squared(bezier_evaluate(segment, t), point);
            if (!result.valid || distance_squared < result.distance_squared)
            {
                result.segment_index = i;
                result.t = t;
                result.distance_squared = distance_squared;
                result.valid = true;
            }
        }
        return result;
    }

    template <typename Curve, typename Table, typename Workspace, typename Projection, typename Tolerance>
    inline bool curve_build_arc_length_table_projected(const Curve &curve, Table &table, Workspace &workspace,
                                                       Projection projection, Tolerance tolerance,
                                                       std::uint32_t max_depth = 12u)
    {
        using point_type = typename Curve::point_type;
        using scalar = typename Table::scalar_type;
        table.clear();
        if (curve.segment_count() == 0u) return false;

        point_type previous{};
        bool has_previous = false;
        scalar accumulated = scalar(0);
        for (std::size_t segment_index = 0u; segment_index < curve.segment_count(); ++segment_index)
        {
            detail::flatten_projected_visit(curve_segment(curve, segment_index), workspace, projection, tolerance,
                                            max_depth, true, [&](const point_type &position, auto t) {
                                                if (has_previous)
                                                    accumulated += static_cast<scalar>(distance(previous, position));
                                                table.entries.push_back(typename Table::entry_type{
                                                    position, segment_index, static_cast<scalar>(t), accumulated});
                                                previous = position;
                                                has_previous = true;
                                            });
        }
        table.total_length = accumulated;
        return true;
    }

    template <typename Curve, typename Table, typename Workspace, typename Tolerance>
    inline bool curve_build_arc_length_table(const Curve &curve, Table &table, Workspace &workspace,
                                             Tolerance tolerance, std::uint32_t max_depth = 12u)
    {
        return curve_build_arc_length_table_projected(curve, table, workspace, detail::identity_projection{}, tolerance,
                                                      max_depth);
    }

    template <typename Table>
    inline bezier_curve_location<typename Table::point_type>
    curve_resolve_arc_length(const Table &table, typename Table::scalar_type distance)
    {
        using point_type = typename Table::point_type;
        using scalar = typename Table::scalar_type;
        bezier_curve_location<point_type> result{};
        if (table.entries.empty()) return result;

        const scalar target = clamp(distance, scalar(0), table.total_length);
        std::size_t low = 0u;
        std::size_t high = table.entries.size();
        while (low < high)
        {
            const std::size_t middle = low + (high - low) / 2u;
            if (table.entries[middle].distance <= target) low = middle + 1u;
            else high = middle;
        }

        const std::size_t left_index = low == 0u ? 0u : low - 1u;
        const auto &left = table.entries[left_index];
        if (low >= table.entries.size())
        {
            result.segment_index = left.segment_index;
            result.t = static_cast<typename decltype(result)::scalar_type>(left.t);
            result.valid = true;
            return result;
        }

        std::size_t right_index = low;
        while (right_index + 1u < table.entries.size() && table.entries[right_index].distance <= left.distance)
            ++right_index;
        const auto &right = table.entries[right_index];
        if (right.segment_index != left.segment_index || right.distance <= left.distance)
        {
            result.segment_index = right.segment_index;
            result.t = static_cast<typename decltype(result)::scalar_type>(right.t);
            result.valid = true;
            return result;
        }

        const scalar factor = (target - left.distance) / (right.distance - left.distance);
        result.segment_index = left.segment_index;
        result.t = static_cast<typename decltype(result)::scalar_type>(mix(left.t, right.t, factor));
        result.valid = true;
        return result;
    }

    template <typename Curve, typename Table>
    inline typename Curve::point_type curve_evaluate_at_distance(const Curve &curve, const Table &table,
                                                                 typename Table::scalar_type distance)
    {
        const auto location = curve_resolve_arc_length(table, distance);
        if (!location.valid) return typename Curve::point_type{};
        return curve_evaluate(curve, location.segment_index, location.t);
    }

    template <typename Curve, typename Scalar>
    inline typename Curve::point_type
    curve_tangent(const Curve &curve, std::size_t segment_index, Scalar t,
                  const typename Curve::point_type &fallback = typename Curve::point_type{})
    {
        return detail::safe_normalize(bezier_derivative(curve_segment(curve, segment_index), t), fallback);
    }

    template <typename Curve, typename Scalar>
    inline typename Curve::point_type
    curve_normal(const Curve &curve, std::size_t segment_index, Scalar t,
                 const typename Curve::point_type &fallback = typename Curve::point_type{})
    {
        const auto segment = curve_segment(curve, segment_index);
        const auto tangent = detail::safe_normalize(bezier_derivative(segment, t));
        const auto acceleration = bezier_second_derivative(segment, t);
        return detail::safe_normalize(acceleration - tangent * detail::point_dot(acceleration, tangent), fallback);
    }

    template <typename Point>
    inline bezier_frame<Point> make_bezier_frame(const Point &position, const Point &direction,
                                                 const Point &reference_up)
    {
        using scalar = bezier_scalar_t<Point>;
        Point tangent = detail::safe_normalize(direction, Point(scalar(0), scalar(0), scalar(1)));
        Point binormal = cross(tangent, reference_up);
        if (length(binormal) <= std::numeric_limits<scalar>::epsilon())
        {
            const Point fallback_axis = abs(tangent.x) < scalar(0.9) ? Point(scalar(1), scalar(0), scalar(0))
                                                                     : Point(scalar(0), scalar(1), scalar(0));
            binormal = cross(tangent, fallback_axis);
        }
        binormal = detail::safe_normalize(binormal);
        const Point normal = detail::safe_normalize(cross(binormal, tangent));
        return {position, tangent, normal, binormal};
    }

    template <typename Curve, typename Scalar>
    inline bezier_frame<typename Curve::point_type> curve_frame(const Curve &curve, std::size_t segment_index, Scalar t,
                                                                const typename Curve::point_type &reference_up)
    {
        const auto segment = curve_segment(curve, segment_index);
        return make_bezier_frame(bezier_evaluate(segment, t), bezier_derivative(segment, t), reference_up);
    }

    template <typename Point>
    inline bezier_frame<Point> transport_bezier_frame(const bezier_frame<Point> &previous, const Point &position,
                                                      const Point &direction)
    {
        using scalar = bezier_scalar_t<Point>;
        const Point tangent = detail::safe_normalize(direction, previous.tangent);
        Point normal = previous.normal - tangent * detail::point_dot(previous.normal, tangent);
        if (length(normal) <= std::numeric_limits<scalar>::epsilon())
            return make_bezier_frame(position, tangent, previous.binormal);
        normal = detail::safe_normalize(normal);
        const Point binormal = detail::safe_normalize(cross(tangent, normal), previous.binormal);
        normal = detail::safe_normalize(cross(binormal, tangent), normal);
        return {position, tangent, normal, binormal};
    }

    template <typename Curve, typename Table>
    inline bezier_frame<typename Curve::point_type>
    curve_frame_at_distance(const Curve &curve, const Table &table, typename Table::scalar_type distance,
                            const typename Curve::point_type &reference_up)
    {
        const auto location = curve_resolve_arc_length(table, distance);
        if (!location.valid) return {};
        return curve_frame(curve, location.segment_index, location.t, reference_up);
    }

    template <typename Curve, typename Table, typename OutputIt>
    inline OutputIt curve_sample_uniform(const Curve &curve, const Table &table, OutputIt output,
                                         std::size_t sample_count, bool include_end = true)
    {
        using scalar = typename Table::scalar_type;
        if (sample_count == 0u || table.entries.empty()) return output;
        if (sample_count == 1u)
        {
            *output++ = curve_evaluate_at_distance(curve, table, scalar(0));
            return output;
        }

        const std::size_t intervals = include_end ? sample_count - 1u : sample_count;
        for (std::size_t i = 0u; i < sample_count; ++i)
        {
            const scalar factor = static_cast<scalar>(i) / static_cast<scalar>(intervals);
            *output++ = curve_evaluate_at_distance(curve, table, table.total_length * factor);
        }
        return output;
    }

    template <typename Curve, typename Table, typename OutputIt>
    inline OutputIt curve_sample_by_spacing(const Curve &curve, const Table &table, OutputIt output,
                                            typename Table::scalar_type spacing, bool include_end = true)
    {
        using scalar = typename Table::scalar_type;
        if (table.entries.empty() || spacing <= std::numeric_limits<scalar>::epsilon()) return output;
        scalar current = scalar(0);
        while (current < table.total_length)
        {
            *output++ = curve_evaluate_at_distance(curve, table, current);
            current += spacing;
        }
        if (include_end) *output++ = curve_evaluate_at_distance(curve, table, table.total_length);
        return output;
    }

    template <typename Point>
    inline void bezier_point_set_position(bezier_point<Point> &point, const Point &position)
    {
        const Point delta = position - point.position;
        point.position = position;
        point.handle_in += delta;
        point.handle_out += delta;
    }

    template <typename Point>
    inline void bezier_point_set_handle(bezier_point<Point> &point, bezier_handle handle, const Point &position,
                                        bezier_handle_mode mode)
    {
        if (handle == bezier_handle::in) point.handle_in = position;
        else point.handle_out = position;

        point.handle_mode = mode;
        if (mode == bezier_handle_mode::independent) return;

        const Point opposite = detail::opposite_handle_position(point, handle, position, mode);
        if (handle == bezier_handle::in) point.handle_out = opposite;
        else point.handle_in = opposite;
    }

    template <typename Point>
    inline void bezier_point_set_handle(bezier_point<Point> &point, bezier_handle handle, const Point &position)
    {
        bezier_point_set_handle(point, handle, position, point.handle_mode);
    }

    template <typename Point>
    inline void bezier_point_set_handle_independent(bezier_point<Point> &point, bezier_handle handle,
                                                    const Point &position)
    {
        bezier_point_set_handle(point, handle, position, bezier_handle_mode::independent);
    }

    template <typename Point>
    inline void bezier_point_break_handles(bezier_point<Point> &point)
    {
        point.handle_mode = bezier_handle_mode::independent;
    }

    template <typename Point>
    inline void bezier_point_set_handle_mode(bezier_point<Point> &point, bezier_handle_mode mode,
                                             bezier_handle driving_handle = bezier_handle::out)
    {
        const Point position = driving_handle == bezier_handle::in ? point.handle_in : point.handle_out;
        bezier_point_set_handle(point, driving_handle, position, mode);
    }

    template <typename Curve>
    inline bool curve_set_point_position(Curve &curve, std::size_t point_index,
                                         const typename Curve::point_type &position)
    {
        if (point_index >= curve.points.size()) return false;
        bezier_point_set_position(curve.points[point_index], position);
        return true;
    }

    template <typename Curve>
    inline bool curve_set_handle(Curve &curve, std::size_t point_index, bezier_handle handle,
                                 const typename Curve::point_type &position, bezier_handle_mode mode)
    {
        if (point_index >= curve.points.size()) return false;
        bezier_point_set_handle(curve.points[point_index], handle, position, mode);
        return true;
    }

    template <typename Curve>
    inline bool curve_set_handle(Curve &curve, std::size_t point_index, bezier_handle handle,
                                 const typename Curve::point_type &position)
    {
        if (point_index >= curve.points.size()) return false;
        bezier_point_set_handle(curve.points[point_index], handle, position);
        return true;
    }

    template <typename Curve>
    inline bool curve_set_handle_independent(Curve &curve, std::size_t point_index, bezier_handle handle,
                                             const typename Curve::point_type &position)
    {
        return curve_set_handle(curve, point_index, handle, position, bezier_handle_mode::independent);
    }

    template <typename Curve>
    inline bool
    curve_is_segment_linear(const Curve &curve, std::size_t segment_index,
                            bezier_scalar_t<typename Curve::point_type> tolerance =
                                std::numeric_limits<bezier_scalar_t<typename Curve::point_type>>::epsilon() *
                                bezier_scalar_t<typename Curve::point_type>(64))
    {
        using scalar = bezier_scalar_t<typename Curve::point_type>;
        if (segment_index >= curve.segment_count()) return false;
        const auto segment = curve_segment(curve, segment_index);
        const scalar scale = max(distance(segment.p0, segment.p3), scalar(1));
        const scalar epsilon = max(tolerance, scalar(0)) * scale;
        return distance(segment.p1, mix(segment.p0, segment.p3, scalar(1) / scalar(3))) <= epsilon &&
               distance(segment.p2, mix(segment.p0, segment.p3, scalar(2) / scalar(3))) <= epsilon;
    }

    template <typename Curve>
    inline bool curve_make_segment_linear(Curve &curve, std::size_t segment_index)
    {
        using scalar = bezier_scalar_t<typename Curve::point_type>;
        if (segment_index >= curve.segment_count()) return false;
        const std::size_t next = (segment_index + 1u) % curve.points.size();
        auto &a = curve.points[segment_index];
        auto &b = curve.points[next];
        a.handle_out = mix(a.position, b.position, scalar(1) / scalar(3));
        b.handle_in = mix(a.position, b.position, scalar(2) / scalar(3));
        a.handle_mode = bezier_handle_mode::independent;
        b.handle_mode = bezier_handle_mode::independent;
        if (!curve.closed && segment_index == 0u) a.handle_in = a.position;
        if (!curve.closed && next + 1u == curve.points.size()) b.handle_out = b.position;
        return true;
    }

    template <typename Curve>
    inline bool curve_make_range_linear(Curve &curve, std::size_t first_segment, std::size_t segment_count)
    {
        if (first_segment > curve.segment_count() || segment_count > curve.segment_count() - first_segment)
            return false;
        for (std::size_t i = 0u; i < segment_count; ++i) curve_make_segment_linear(curve, first_segment + i);
        return true;
    }

    template <typename Curve>
    inline bool curve_interpolate_segment(
        Curve &curve, std::size_t segment_index, bezier_interpolation interpolation,
        bezier_scalar_t<typename Curve::point_type> tension = bezier_scalar_t<typename Curve::point_type>(0))
    {
        using scalar = bezier_scalar_t<typename Curve::point_type>;
        if (segment_index >= curve.segment_count()) return false;
        if (interpolation == bezier_interpolation::linear) return curve_make_segment_linear(curve, segment_index);

        const std::size_t count = curve.points.size();
        const std::size_t i0 = segment_index;
        const std::size_t i1 = (segment_index + 1u) % count;
        const std::size_t prev = i0 == 0u ? (curve.closed ? count - 1u : i0) : i0 - 1u;
        const std::size_t next = i1 + 1u < count ? i1 + 1u : (curve.closed ? 0u : i1);
        const scalar scale = (scalar(1) - clamp(tension, scalar(0), scalar(1))) / scalar(6);
        curve.points[i0].handle_out =
            curve.points[i0].position + (curve.points[i1].position - curve.points[prev].position) * scale;
        curve.points[i1].handle_in =
            curve.points[i1].position - (curve.points[next].position - curve.points[i0].position) * scale;
        curve.points[i0].handle_mode = bezier_handle_mode::independent;
        curve.points[i1].handle_mode = bezier_handle_mode::independent;
        return true;
    }

    template <typename Curve>
    inline bool curve_interpolate_range(
        Curve &curve, std::size_t first_segment, std::size_t segment_count, bezier_interpolation interpolation,
        bezier_scalar_t<typename Curve::point_type> tension = bezier_scalar_t<typename Curve::point_type>(0))
    {
        if (first_segment > curve.segment_count() || segment_count > curve.segment_count() - first_segment)
            return false;
        for (std::size_t i = 0u; i < segment_count; ++i)
            curve_interpolate_segment(curve, first_segment + i, interpolation, tension);
        if (interpolation == bezier_interpolation::smooth)
        {
            if (curve.closed && segment_count == curve.segment_count())
            {
                for (auto &point : curve.points) point.handle_mode = bezier_handle_mode::mirrored;
            }
            else
            {
                const std::size_t first_inner = first_segment + 1u;
                const std::size_t last_inner = first_segment + segment_count;
                for (std::size_t point_index = first_inner; point_index < last_inner; ++point_index)
                    curve.points[point_index].handle_mode = bezier_handle_mode::mirrored;
            }
        }
        return true;
    }

    template <typename Curve, typename Scalar>
    inline bool curve_insert_point_preserving_shape(Curve &curve, std::size_t segment_index, Scalar t)
    {
        using scalar = bezier_scalar_t<typename Curve::point_type>;
        if (segment_index >= curve.segment_count()) return false;
        const scalar value = static_cast<scalar>(t);
        if (!(value > scalar(0) && value < scalar(1))) return false;

        const std::size_t next = (segment_index + 1u) % curve.points.size();
        const auto halves = bezier_split(curve_segment(curve, segment_index), value);
        curve.points[segment_index].handle_out = halves.left.p1;
        curve.points[segment_index].handle_mode = bezier_handle_mode::independent;
        curve.points[next].handle_in = halves.right.p2;
        curve.points[next].handle_mode = bezier_handle_mode::independent;
        typename Curve::knot_type inserted{halves.left.p3, halves.left.p2, halves.right.p1,
                                           bezier_handle_mode::aligned};

        if (curve.closed && next == 0u) curve.points.push_back(inserted);
        else curve.points.insert(curve.points.begin() + static_cast<std::ptrdiff_t>(next), inserted);
        return true;
    }

    template <typename Curve>
    inline void curve_append_point(
        Curve &curve, const typename Curve::point_type &position,
        bezier_interpolation interpolation = bezier_interpolation::linear,
        bezier_scalar_t<typename Curve::point_type> tension = bezier_scalar_t<typename Curve::point_type>(0))
    {
        curve.points.push_back(typename Curve::knot_type(position));
        if (curve.segment_count() == 0u) return;
        if (interpolation == bezier_interpolation::linear)
        {
            curve_make_segment_linear(curve, curve.segment_count() - 1u);
            return;
        }

        const std::size_t first = curve.segment_count() > 1u ? curve.segment_count() - 2u : 0u;
        curve_interpolate_range(curve, first, curve.segment_count() - first, interpolation, tension);
        if (first > 0u) curve.points[first].handle_mode = bezier_handle_mode::mirrored;
    }

    template <typename Curve, typename InputIt>
    inline void curve_append_points(
        Curve &curve, InputIt first, InputIt last, bezier_interpolation interpolation = bezier_interpolation::linear,
        bezier_scalar_t<typename Curve::point_type> tension = bezier_scalar_t<typename Curve::point_type>(0))
    {
        const std::size_t first_new_point = curve.points.size();
        for (; first != last; ++first) curve.points.push_back(typename Curve::knot_type(*first));
        if (curve.segment_count() == 0u || first_new_point == curve.points.size()) return;
        const std::size_t first_segment = first_new_point == 0u ? 0u : first_new_point - 1u;
        curve_interpolate_range(curve, first_segment, curve.segment_count() - first_segment, interpolation, tension);
        if (interpolation == bezier_interpolation::smooth && first_segment > 0u)
        {
            curve_interpolate_segment(curve, first_segment - 1u, interpolation, tension);
            curve.points[first_segment].handle_mode = bezier_handle_mode::mirrored;
        }
    }

    template <typename Curve, typename Scalar>
    inline bool curve_insert_point_linear(Curve &curve, std::size_t segment_index, Scalar t)
    {
        using scalar = bezier_scalar_t<typename Curve::point_type>;
        if (segment_index >= curve.segment_count()) return false;
        const scalar value = static_cast<scalar>(t);
        if (!(value > scalar(0) && value < scalar(1))) return false;
        const auto source = curve_segment(curve, segment_index);
        const typename Curve::point_type position = mix(source.p0, source.p3, value);
        const std::size_t next = (segment_index + 1u) % curve.points.size();
        typename Curve::knot_type inserted(position);
        inserted.handle_mode = bezier_handle_mode::aligned;
        if (curve.closed && next == 0u) curve.points.push_back(inserted);
        else curve.points.insert(curve.points.begin() + static_cast<std::ptrdiff_t>(next), inserted);
        curve_make_segment_linear(curve, segment_index);
        curve_make_segment_linear(curve, segment_index + 1u);
        return true;
    }

    template <typename Curve, typename Scalar>
    inline bool curve_insert_point_auto(Curve &curve, std::size_t segment_index, Scalar t)
    {
        if (segment_index >= curve.segment_count()) return false;
        if (curve_is_segment_linear(curve, segment_index)) return curve_insert_point_linear(curve, segment_index, t);
        return curve_insert_point_preserving_shape(curve, segment_index, t);
    }

    template <typename Curve>
    inline bool curve_insert_point_approximated(
        Curve &curve, std::size_t segment_index, const typename Curve::point_type &position,
        bezier_scalar_t<typename Curve::point_type> handle_scale = bezier_scalar_t<typename Curve::point_type>(1) /
                                                                   bezier_scalar_t<typename Curve::point_type>(3))
    {
        using scalar = bezier_scalar_t<typename Curve::point_type>;
        using point_type = typename Curve::point_type;
        if (segment_index >= curve.segment_count()) return false;
        const std::size_t next = (segment_index + 1u) % curve.points.size();
        const point_type &a = curve.points[segment_index].position;
        const point_type &b = curve.points[next].position;
        point_type direction = b - a;
        const scalar direction_length = length(direction);
        if (direction_length > std::numeric_limits<scalar>::epsilon()) direction /= direction_length;
        else direction = point_type{};
        const scalar scale = max(handle_scale, scalar(0));
        typename Curve::knot_type inserted{position, position - direction * (distance(a, position) * scale),
                                           position + direction * (distance(position, b) * scale),
                                           bezier_handle_mode::aligned};
        if (curve.closed && next == 0u) curve.points.push_back(inserted);
        else curve.points.insert(curve.points.begin() + static_cast<std::ptrdiff_t>(next), inserted);
        return true;
    }

    template <typename Curve>
    struct bezier_curve_split
    {
        Curve left{};
        Curve right{};
        bool valid = false;
    };

    template <typename Curve, typename Scalar>
    inline bezier_curve_split<Curve> curve_split(const Curve &source, std::size_t segment_index, Scalar t)
    {
        bezier_curve_split<Curve> result{};
        if (source.closed || segment_index >= source.segment_count()) return result;
        Curve divided = source;
        if (!curve_insert_point_preserving_shape(divided, segment_index, t)) return result;
        const std::size_t split_index = segment_index + 1u;
        for (std::size_t i = 0u; i <= split_index; ++i) detail::append_point(result.left, divided.points[i]);
        for (std::size_t i = split_index; i < divided.points.size(); ++i)
            detail::append_point(result.right, divided.points[i]);
        result.valid = true;
        return result;
    }

    template <typename Curve>
    inline Curve curve_join(
        const Curve &left, const Curve &right,
        bezier_scalar_t<typename Curve::point_type> weld_tolerance = bezier_scalar_t<typename Curve::point_type>(0))
    {
        Curve result = left;
        result.closed = false;
        if (right.points.empty()) return result;
        if (result.points.empty())
        {
            result = right;
            result.closed = false;
            return result;
        }

        std::size_t first = 0u;
        if (distance(result.points.back().position, right.points.front().position) <=
            max(weld_tolerance, bezier_scalar_t<typename Curve::point_type>(0)))
        {
            result.points.back().handle_out = right.points.front().handle_out;
            result.points.back().handle_mode = bezier_handle_mode::independent;
            first = 1u;
        }
        for (std::size_t i = first; i < right.points.size(); ++i) detail::append_point(result, right.points[i]);
        return result;
    }

    template <typename Curve>
    inline void curve_reverse(Curve &curve)
    {
        std::reverse(curve.points.begin(), curve.points.end());
        for (auto &point : curve.points) std::swap(point.handle_in, point.handle_out);
    }

    template <typename Curve, typename Transform>
    inline void curve_transform(Curve &curve, Transform transform)
    {
        for (auto &point : curve.points)
        {
            point.position = transform(point.position);
            point.handle_in = transform(point.handle_in);
            point.handle_out = transform(point.handle_out);
        }
    }

    template <typename Curve>
    inline void curve_translate(Curve &curve, const typename Curve::point_type &delta)
    {
        curve_transform(curve, [&](const typename Curve::point_type &point) { return point + delta; });
    }

    template <typename Curve>
    inline bool curve_remove_point_approximated(
        Curve &curve, std::size_t point_index,
        bezier_scalar_t<typename Curve::point_type> tension = bezier_scalar_t<typename Curve::point_type>(0))
    {
        if (point_index >= curve.points.size()) return false;
        if (!curve.closed && (point_index == 0u || point_index + 1u == curve.points.size())) return false;
        if (curve.points.size() <= (curve.closed ? 3u : 2u)) return false;

        const std::size_t previous = point_index == 0u ? curve.points.size() - 1u : point_index - 1u;
        curve.points.erase(curve.points.begin() + static_cast<std::ptrdiff_t>(point_index));
        const std::size_t segment = previous >= curve.points.size() ? curve.points.size() - 1u : previous;
        return curve_interpolate_segment(curve, segment, bezier_interpolation::smooth, tension);
    }

    template <typename Curve, typename OutputIt>
    inline bool curve_sample(const Curve &curve, OutputIt output, std::uint32_t samples_per_segment = 32u)
    {
        using scalar = bezier_scalar_t<typename Curve::point_type>;
        if (curve.segment_count() == 0u || samples_per_segment == 0u) return false;
        for (std::size_t segment_index = 0u; segment_index < curve.segment_count(); ++segment_index)
        {
            const auto segment = curve_segment(curve, segment_index);
            const auto coefficients = detail::make_bezier_coefficients(segment);
            const std::uint32_t first = segment_index == 0u ? 0u : 1u;
            for (std::uint32_t i = first; i <= samples_per_segment; ++i)
                *output++ = detail::bezier_evaluate_unchecked(
                    coefficients, static_cast<scalar>(i) / static_cast<scalar>(samples_per_segment));
        }
        return true;
    }

    template <typename Curve, typename OutputIt>
    inline bool curve_flatten(
        const Curve &curve, OutputIt output,
        bezier_scalar_t<typename Curve::point_type> tolerance = bezier_scalar_t<typename Curve::point_type>(0.5),
        std::uint32_t max_depth = 12u)
    {
        if (curve.segment_count() == 0u) return false;
        for (std::size_t i = 0u; i < curve.segment_count(); ++i)
            output = bezier_flatten(curve_segment(curve, i), output, tolerance, max_depth, i == 0u);
        return true;
    }

    template <typename Curve, typename OutputIt, typename Workspace, typename Projection, typename Tolerance>
    inline OutputIt curve_flatten_projected_with_workspace(const Curve &curve, OutputIt output, Workspace &workspace,
                                                           Projection projection, Tolerance tolerance,
                                                           std::uint32_t max_depth = 12u)
    {
        for (std::size_t i = 0u; i < curve.segment_count(); ++i)
            output = bezier_flatten_projected_with_workspace(curve_segment(curve, i), output, workspace, projection,
                                                             tolerance, max_depth, i == 0u);
        return output;
    }

    template <typename Curve, typename OutputIt, typename Workspace, typename Tolerance>
    inline OutputIt curve_flatten_with_workspace(const Curve &curve, OutputIt output, Workspace &workspace,
                                                 Tolerance tolerance, std::uint32_t max_depth = 12u)
    {
        return curve_flatten_projected_with_workspace(curve, output, workspace, detail::identity_projection{},
                                                      tolerance, max_depth);
    }

    template <typename Curve>
    inline bezier_scalar_t<typename Curve::point_type> curve_approximate_length(const Curve &curve,
                                                                                std::uint32_t steps = 32u)
    {
        using scalar = bezier_scalar_t<typename Curve::point_type>;
        scalar result = scalar(0);
        for (std::size_t i = 0u; i < curve.segment_count(); ++i)
            result += bezier_approximate_length(curve_segment(curve, i), steps);
        return result;
    }
} // namespace amal
