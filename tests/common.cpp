#include <amal/common.hpp>
#include <amal/vector.hpp>
#include <cassert>
#include <cmath>
#include <cstdio>

using namespace amal;

void test_common()
{
    {
        vec3 v(-1.0f, 0.0f, 2.0f);
        vec3 r = abs(v);
        assert(r.x == 1.0f && r.y == 0.0f && r.z == 2.0f);
        printf("abs passed\n");
    }

    {
        vec3 v(-5.0f, 0.0f, 7.0f);
        vec3 r = sign(v);
        assert(r.x == -1.0f && r.y == 0.0f && r.z == 1.0f);
        printf("sign passed\n");
    }

    {
        vec3 v(1.2f, -3.7f, 4.9f);
        vec3 r = floor(v);
        assert(r.x == 1.0f && r.y == -4.0f && r.z == 4.0f);
        printf("floor passed\n");
    }

    {
        vec3 v(1.2f, -3.7f, 4.9f);
        vec3 r = ceil(v);
        assert(r.x == 2.0f && r.y == -3.0f && r.z == 5.0f);
        printf("ceil passed\n");
    }

    {
        vec3 v(1.75f, -1.25f, 0.5f);
        vec3 r = fract(v);
        assert(abs(r.x - 0.75f) < 1e-5f);
        assert(abs(r.y - 0.75f) < 1e-5f); // fract(-1.25) = 0.75
        assert(abs(r.z - 0.5f) < 1e-5f);
        printf("fract passed\n");
    }

    {
        vec3 a(1.0f, 2.0f, 3.0f);
        vec3 b(2.0f, 1.0f, 5.0f);
        vec3 r = min(a, b);
        assert(r.x == 1.0f && r.y == 1.0f && r.z == 3.0f);
        printf("min passed\n");
    }

    {
        vec3 a(1.0f, 2.0f, 3.0f);
        vec3 b(2.0f, 1.0f, 5.0f);
        vec3 r = max(a, b);
        assert(r.x == 2.0f && r.y == 2.0f && r.z == 5.0f);
        printf("max passed\n");
    }

    {
        vec3 v(2.0f, -1.0f, 5.0f);
        vec3 r = clamp(v, 0.0f, 3.0f);
        assert(r.x == 2.0f && r.y == 0.0f && r.z == 3.0f);
        printf("clamp passed\n");
    }

    {
        vec3 a(0.0f, 2.0f, 4.0f);
        vec3 b(1.0f, 4.0f, 2.0f);
        vec3 r = mix(a, b, 0.5f); // expected: (0.5, 3.0, 3.0)
        assert(abs(r.x - 0.5f) < 1e-5f);
        assert(abs(r.y - 3.0f) < 1e-5f);
        assert(abs(r.z - 3.0f) < 1e-5f);
        printf("mix passed\n");
    }

    {
        vec3 v(0.2f, 0.5f, 0.8f);
        vec3 r = step(0.5f, v);
        assert(r.x == 0.0f && r.y == 1.0f && r.z == 1.0f);
        printf("step passed\n");
    }

    {
        vec3 v(0.2f, 0.5f, 0.8f);
        vec3 r = smoothstep(0.0f, 1.0f, v); // smoothstep on [0..1]
        assert(r.x > 0.0f && r.x < 1.0f);
        assert(abs(r.y - 0.5f) < 0.01f); // mid point
        assert(r.z > 0.5f && r.z < 1.0f);
        printf("smoothstep passed\n");
    }

    {
        vec3 a(5.5f, -3.5f, 8.0f);
        vec3 r = mod(a, 3.0f);
        assert(abs(r.x - 2.5f) < 1e-5f);
        assert(abs(r.y - 2.5f) < 1e-5f);
        assert(abs(r.z - 2.0f) < 1e-5f);
        printf("mod passed\n");
    }

    printf("All tests passed!\n");
}
