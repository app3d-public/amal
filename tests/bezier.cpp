#include <amal/bezier.hpp>

#ifdef near
    #undef near
#endif

using namespace amal;

namespace
{
    bool near(float a, float b, float epsilon = 1e-4f) { return abs(a - b) <= epsilon; }
    bool near(const vec2 &a, const vec2 &b, float epsilon = 1e-4f)
    {
        return near(a.x, b.x, epsilon) && near(a.y, b.y, epsilon);
    }
    bool near(const vec3 &a, const vec3 &b, float epsilon = 1e-4f)
    {
        return near(a.x, b.x, epsilon) && near(a.y, b.y, epsilon) && near(a.z, b.z, epsilon);
    }

    bezier_curve2 make_test_curve()
    {
        bezier_curve2 curve;
        curve.points = {
            {{0.0f, 0.0f}, {0.0f, 0.0f}, {1.0f, 2.0f}, bezier_handle_mode::independent},
            {{3.0f, 2.0f}, {2.0f, 3.0f}, {4.0f, 1.0f}, bezier_handle_mode::independent},
            {{6.0f, 0.0f}, {5.0f, -1.0f}, {6.0f, 0.0f}, bezier_handle_mode::independent},
        };
        return curve;
    }
} // namespace

void test_bezier()
{
    {
        const cubic_bezier<vec2> curve{{0.0f, 0.0f}, {1.0f, 2.0f}, {2.0f, 2.0f}, {3.0f, 0.0f}};
        assert(near(bezier_evaluate(curve, 0.0f), curve.p0));
        assert(near(bezier_evaluate(curve, 1.0f), curve.p3));
        assert(near(bezier_evaluate(curve, 0.5f), {1.5f, 1.5f}));
        assert(near(bezier_derivative(curve, 0.0f), (curve.p1 - curve.p0) * 3.0f));

        for (std::uint32_t i = 0u; i <= 20u; ++i)
        {
            const float t = static_cast<float>(i) / 20.0f;
            assert(near(bezier_evaluate_unchecked(curve, t), bezier_evaluate(curve, t)));
            assert(near(bezier_derivative_unchecked(curve, t), bezier_derivative(curve, t)));
            assert(near(bezier_second_derivative_unchecked(curve, t), bezier_second_derivative(curve, t)));
        }

        const auto halves = bezier_split(curve, 0.35f);
        assert(near(halves.left.p3, halves.right.p0));
        for (std::uint32_t i = 0u; i <= 20u; ++i)
        {
            const float t = static_cast<float>(i) / 20.0f;
            if (t <= 0.35f) assert(near(bezier_evaluate(curve, t), bezier_evaluate(halves.left, t / 0.35f), 2e-4f));
            else assert(near(bezier_evaluate(curve, t), bezier_evaluate(halves.right, (t - 0.35f) / 0.65f), 2e-4f));
        }
    }

    {
        bezier_point<vec2> point{{2.0f, 2.0f}, {1.0f, 2.0f}, {3.0f, 2.0f}, bezier_handle_mode::mirrored};
        const vec2 preserved_in = point.handle_in;
        bezier_point_set_handle_independent(point, bezier_handle::out, {4.0f, 3.0f});
        assert(near(point.handle_in, preserved_in));
        assert(near(point.handle_out, {4.0f, 3.0f}));
        assert(point.handle_mode == bezier_handle_mode::independent);

        bezier_point_set_handle(point, bezier_handle::out, {3.0f, 2.0f}, bezier_handle_mode::mirrored);
        assert(near(point.handle_in, {1.0f, 2.0f}));
        bezier_point_set_handle(point, bezier_handle::out, {2.0f, 4.0f}, bezier_handle_mode::aligned);
        assert(near(point.handle_in, {2.0f, 1.0f}));
    }

    {
        auto curve = make_test_curve();
        const auto original_first = curve_segment(curve, 0u);
        const auto original_second = curve_segment(curve, 1u);
        const vec2 preserved_in = curve.points[1].handle_in;
        assert(curve_set_handle_independent(curve, 1u, bezier_handle::out, {5.0f, 2.0f}));
        assert(near(curve.points[1].handle_in, preserved_in));
        assert(!curve_set_handle_independent(curve, 99u, bezier_handle::out, {0.0f, 0.0f}));
        curve.points[1].handle_out = original_second.p1;
        assert(curve_insert_point_preserving_shape(curve, 0u, 0.4f));
        assert(curve.points.size() == 4u);
        for (std::uint32_t i = 0u; i <= 20u; ++i)
        {
            const float t = static_cast<float>(i) / 20.0f;
            const vec2 expected = bezier_evaluate(original_first, t);
            const vec2 actual = t <= 0.4f ? bezier_evaluate(curve_segment(curve, 0u), t / 0.4f)
                                          : bezier_evaluate(curve_segment(curve, 1u), (t - 0.4f) / 0.6f);
            assert(near(expected, actual, 3e-4f));
        }
        assert(near(curve_segment(curve, 2u).p0, original_second.p0));
        assert(near(curve_segment(curve, 2u).p3, original_second.p3));
    }

    {
        auto curve = make_test_curve();
        const auto divided = curve_split(curve, 0u, 0.25f);
        assert(divided.valid);
        assert(divided.left.points.size() == 2u);
        assert(divided.right.points.size() == 3u);
        assert(near(divided.left.points.back().position, divided.right.points.front().position));
        const auto joined = curve_join(divided.left, divided.right, 1e-5f);
        assert(joined.points.size() == 4u);
        assert(joined.segment_count() == 3u);
    }

    {
        bezier_curve2 linear;
        linear.points = {bezier_point<vec2>({0.0f, 0.0f}), bezier_point<vec2>({9.0f, 0.0f})};
        assert(curve_make_segment_linear(linear, 0u));
        assert(near(linear.points[0].handle_out, {3.0f, 0.0f}));
        assert(near(linear.points[1].handle_in, {6.0f, 0.0f}));
        assert(curve_is_segment_linear(linear, 0u));
        assert(curve_insert_point_auto(linear, 0u, 0.5f));
        assert(linear.points.size() == 3u);
        assert(near(linear.points[1].position, {4.5f, 0.0f}));
        assert(curve_is_segment_linear(linear, 0u));
        assert(curve_is_segment_linear(linear, 1u));
    }

    {
        auto interpolated = make_test_curve();
        const auto original = curve_segment(interpolated, 0u);
        assert(!curve_is_segment_linear(interpolated, 0u));
        assert(curve_insert_point_auto(interpolated, 0u, 0.4f));
        for (std::uint32_t i = 0u; i <= 20u; ++i)
        {
            const float t = static_cast<float>(i) / 20.0f;
            const vec2 expected = bezier_evaluate(original, t);
            const vec2 actual = t <= 0.4f ? bezier_evaluate(curve_segment(interpolated, 0u), t / 0.4f)
                                          : bezier_evaluate(curve_segment(interpolated, 1u), (t - 0.4f) / 0.6f);
            assert(near(expected, actual, 3e-4f));
        }
    }

    {
        const cubic_bezier<vec2> line{{0.0f, 0.0f}, {1.0f, 0.0f}, {2.0f, 0.0f}, {3.0f, 0.0f}};
        const float t = bezier_closest_parameter(line, vec2{1.2f, 2.0f});
        assert(near(t, 0.4f, 2e-3f));
        assert(near(bezier_closest_point(line, vec2{1.2f, 2.0f}), {1.2f, 0.0f}, 2e-3f));
        assert(near(bezier_approximate_length(line), 3.0f));

        acul::vector<vec2> flattened;
        bezier_flatten(line, std::back_inserter(flattened), 0.01f);
        assert(flattened.size() == 2u);
        assert(near(flattened.front(), line.p0));
        assert(near(flattened.back(), line.p3));
    }

    {
        auto curve = make_test_curve();
        const vec2 original = curve.points[0].position;
        assert(curve_remove_point_approximated(curve, 1u));
        assert(curve.points.size() == 2u);
        assert(near(curve.points[0].position, original));
        curve_reverse(curve);
        assert(near(curve.points.back().position, original));

        acul::vector<vec2> samples;
        assert(curve_sample(curve, std::back_inserter(samples), 8u));
        assert(samples.size() == 9u);
    }

    {
        bezier_curve2 curve;
        const vec2 points[] = {{0.0f, 0.0f}, {2.0f, 1.0f}, {4.0f, 0.0f}};
        curve_append_points(curve, std::begin(points), std::end(points), bezier_interpolation::smooth);
        assert(curve.points.size() == 3u);
        assert(curve.segment_count() == 2u);

        const auto location = curve_closest_location(curve, vec2{2.0f, 1.2f});
        assert(location.valid);
        assert(location.segment_index < curve.segment_count());

        acul::vector<vec2> flattened;
        assert(curve_flatten(curve, std::back_inserter(flattened), 0.05f));
        assert(flattened.size() >= 3u);
        const float old_length = curve_approximate_length(curve);
        curve_translate(curve, vec2{10.0f, -3.0f});
        assert(near(curve.points.front().position, {10.0f, -3.0f}));
        assert(near(curve_approximate_length(curve), old_length));
    }

    {
        bezier_curve2 curve;
        curve_append_point(curve, vec2{0.0f, 0.0f});
        curve_append_point(curve, vec2{1.0f, 1.0f});
        assert(curve.segment_count() == 1u);
        assert(curve_insert_point_preserving_shape(curve, 0u, 0.5f));
        assert(curve.points.size() == 3u);
    }

    {
        auto curve = make_test_curve();
        bezier_flatten_workspace2 workspace;
        acul::vector<vec2> world_points;
        curve_flatten_with_workspace(curve, std::back_inserter(world_points), workspace, 0.05f);

        acul::vector<vec2> screen_points;
        const auto projection = [](const vec2 &point) { return vec2{point.x, point.y * 100.0f}; };
        curve_flatten_projected_with_workspace(curve, std::back_inserter(screen_points), workspace, projection, 0.05f);
        assert(screen_points.size() > world_points.size());
        assert(near(screen_points.front(), curve.points.front().position));
        assert(near(screen_points.back(), curve.points.back().position));
    }

    {
        bezier_curve3 curve;
        curve_append_point(curve, vec3{0.0f, 0.0f, 0.0f});
        curve_append_point(curve, vec3{3.0f, 0.0f, 0.0f});
        curve_append_point(curve, vec3{3.0f, 4.0f, 0.0f});

        bezier_flatten_workspace3 workspace;
        bezier_arc_length_table3 table;
        assert(curve_build_arc_length_table(curve, table, workspace, 1e-4f));
        assert(near(table.total_length, 7.0f, 1e-3f));

        const auto first_location = curve_resolve_arc_length(table, 1.5f);
        assert(first_location.segment_index == 0u);
        assert(near(first_location.t, 0.5f, 1e-3f));
        const auto second_location = curve_resolve_arc_length(table, 5.0f);
        assert(second_location.segment_index == 1u);
        assert(near(second_location.t, 0.5f, 1e-3f));
        assert(near(curve_evaluate_at_distance(curve, table, 5.0f), {3.0f, 2.0f, 0.0f}, 1e-3f));

        acul::vector<vec3> uniform_samples;
        curve_sample_uniform(curve, table, std::back_inserter(uniform_samples), 1000u);
        assert(uniform_samples.size() == 1000u);
        assert(near(uniform_samples.front(), curve.points.front().position));
        assert(near(uniform_samples.back(), curve.points.back().position));

        acul::vector<vec3> spaced_samples;
        curve_sample_by_spacing(curve, table, std::back_inserter(spaced_samples), 2.0f);
        assert(spaced_samples.size() == 5u);
        assert(near(spaced_samples.back(), curve.points.back().position));

        const auto frame = curve_frame_at_distance(curve, table, 1.0f, vec3{0.0f, 0.0f, 1.0f});
        assert(near(length(frame.tangent), 1.0f));
        assert(near(length(frame.normal), 1.0f));
        assert(near(length(frame.binormal), 1.0f));
        assert(near(dot(frame.tangent, frame.normal), 0.0f));
        assert(near(dot(frame.tangent, frame.binormal), 0.0f));
        assert(near(dot(frame.normal, frame.binormal), 0.0f));

        const auto transported = transport_bezier_frame(frame, vec3{3.0f, 1.0f, 0.0f}, vec3{0.0f, 1.0f, 0.0f});
        assert(near(transported.tangent, {0.0f, 1.0f, 0.0f}));
        assert(near(dot(transported.tangent, transported.normal), 0.0f));
    }

    {
        bezier_curve3 curve;
        for (std::uint32_t i = 0u; i < 1000u; ++i) curve_append_point(curve, vec3{static_cast<float>(i), 0.0f, 0.0f});
        assert(curve.points.size() == 1000u);
        assert(curve.segment_count() == 999u);

        bezier_flatten_workspace3 workspace;
        bezier_arc_length_table3 table;
        assert(curve_build_arc_length_table(curve, table, workspace, 1e-4f));
        assert(near(table.total_length, 999.0f, 1e-3f));
        assert(near(curve_evaluate_at_distance(curve, table, 999.0f), {999.0f, 0.0f, 0.0f}));
    }

    printf("test_bezier passed!\n");
}
