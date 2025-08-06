#include <amal/transform.hpp>
#include <amal/trigonometric.hpp>
#include <cstdio>

using namespace amal;

void test_transform()
{
    {
        vec2 screen(400.0f, 300.0f);
        vec2 ndc = screen_to_ndc<float, false>(screen, 800, 600);
        assert(abs(ndc.x) < 1e-5f);
        assert(abs(ndc.y) < 1e-5f);
        printf("screen_to_ndc passed\n");
    }

    {
        // project + unproject
        vec3 obj(1.0f, 2.0f, 3.0f);
        mat4 model = translate(mat4(1.0f), vec3(0.0f, 0.0f, 0.0f));
        mat4 view = look_at(vec3(0.0f, 0.0f, 10.0f), vec3(0.0f), vec3(0.0f, 1.0f, 0.0f));
        mat4 proj = perspective(radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
        vec4 viewport(0.0f, 0.0f, 800.0f, 600.0f);

        vec3 win = project(obj, view * model, proj, viewport);
        vec3 obj2 = unproject(win, view * model, proj, viewport);
        assert(length(obj - obj2) < 1e-3f);
        printf("project/unproject passed\n");
    }

    {
        // translate
        mat4 m(1.0f);
        m = translate(m, vec3(2.0f, 3.0f, 4.0f));
        vec4 p = m * vec4(1.0f);
        assert(abs(p.x - 3.0f) < 1e-5f);
        assert(abs(p.y - 4.0f) < 1e-5f);
        assert(abs(p.z - 5.0f) < 1e-5f);
        printf("translate passed\n");
    }

    {
        // rotate (around Z)
        float angle = radians(90.0f);
        mat4 m = rotate(mat4(1.0f), angle, vec3(0.0f, 0.0f, 1.0f));
        vec4 p = m * vec4(1.0f, 0.0f, 0.0f, 1.0f);
        assert(abs(p.x) < 1e-5f);
        assert(abs(p.y - 1.0f) < 1e-5f);
        printf("rotate passed\n");
    }

    {
        // scale
        mat4 m = scale(mat4(1.0f), vec3(2.0f, 3.0f, 4.0f));
        vec4 p = m * vec4(1.0f);
        assert(abs(p.x - 2.0f) < 1e-5f);
        assert(abs(p.y - 3.0f) < 1e-5f);
        assert(abs(p.z - 4.0f) < 1e-5f);
        printf("scale passed\n");
    }

    {
        // look_at (basic direction)
        vec3 eye(0.0f, 0.0f, 3.0f);
        vec3 center(0.0f, 0.0f, 0.0f);
        vec3 up(0.0f, 1.0f, 0.0f);
        mat4 view = look_at(eye, center, up);
        vec4 forward = view * vec4(0.0f, 0.0f, -1.0f, 0.0f);
        assert(abs(forward.z + 1.0f) < 1e-5f);
        printf("look_at passed\n");
    }

    {
        // orthographic
        mat4 proj = orthographic(-1.0f, 1.0f, -1.0f, 1.0f);
        vec4 p = proj * vec4(1.0f, 1.0f, 0.0f, 1.0f);
        assert(abs(p.x - 1.0f) < 1e-5f);
        assert(abs(p.y - 1.0f) < 1e-5f);
        printf("orthographic passed\n");
    }

    {
        mat4 proj = perspective(radians(90.0f), 1.0f, 0.1f, 10.0f);
        vec4 p = proj * vec4(0.0f, 0.0f, -1.0f, 1.0f);
        printf("%f %f %f %f\n", p.x, p.y, p.z, p.w);
#if defined(AMAL_CLIP_SPACE_NO)
        assert(p.w > 0.0f && p.w < 1.0f);
        assert(p.z >= -1.0f && p.z <= 1.0f);
#else
        assert(p.w == 1.0f);
        assert(p.z >= 0.0f && p.z <= 1.0f);
#endif
        printf("perspective passed\n");
    }

    printf("All transform tests passed!\n");
}