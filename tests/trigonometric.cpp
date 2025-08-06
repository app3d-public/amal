#include <amal/trigonometric.hpp>
#include <cstdio>

using namespace amal;

void test_trigonometric()
{
    {
        vec3 v(0.0f, half_pi<float>(), pi<float>());
        vec3 s = sin(v);
        assert(abs(s.x - 0.0f) < 1e-5f);
        assert(abs(s.y - 1.0f) < 1e-5f);
        assert(abs(s.z - 0.0f) < 1e-5f);
        printf("sin passed\n");
    }

    {
        vec3 v(0.0f, half_pi<float>(), pi<float>());
        vec3 c = cos(v);
        assert(abs(c.x - 1.0f) < 1e-5f);
        assert(abs(c.y - 0.0f) < 1e-5f);
        assert(abs(c.z + 1.0f) < 1e-5f);
        printf("cos passed\n");
    }

    {
        vec3 v(0.0f, 1.0f, -1.0f);
        vec3 t = tan(v);
        assert(abs(t.x - 0.0f) < 1e-5f);
        assert(t.y > 1.5f && t.y < 2.0f);
        assert(t.z < -1.5f && t.z > -2.0f);
        printf("tan passed\n");
    }

    {
        vec3 v(0.0f, 0.5f, 1.0f);
        vec3 r = asin(v);
        assert(abs(r.x - 0.0f) < 1e-5f);
        assert(abs(r.y - std::asin(0.5f)) < 1e-5f);
        assert(abs(r.z - std::asin(1.0f)) < 1e-5f);
        printf("asin passed\n");
    }

    {
        vec3 v(1.0f, 0.5f, 0.0f);
        vec3 r = acos(v);
        assert(abs(r.x - 0.0f) < 1e-5f);
        assert(abs(r.y - std::acos(0.5f)) < 1e-5f);
        assert(abs(r.z - half_pi<float>()) < 1e-5f);
        printf("acos passed\n");
    }

    {
        vec3 v(0.0f, 1.0f, -1.0f);
        vec3 r = atan(v);
        assert(abs(r.x - 0.0f) < 1e-5f);
        assert(abs(r.y - std::atan(1.0f)) < 1e-5f);
        assert(abs(r.z - std::atan(-1.0f)) < 1e-5f);
        printf("atan(vec) passed\n");
    }

    {
        vec3 y(1.0f, 0.0f, -1.0f);
        vec3 x(1.0f, 1.0f, -1.0f);
        vec3 r = atan2(y, x);
        assert(abs(r.x - quarter_pi<float>()) < 1e-5f);
        assert(abs(r.y - 0.0f) < 1e-5f);
        assert(abs(r.z + 3 * quarter_pi<float>()) < 1e-5f);
        printf("atan2(vec, vec) passed\n");
    }

    printf("test_trigonometric passed!\n");
}