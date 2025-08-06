#include <amal/common.hpp>
#include <amal/compare.hpp>
#include <amal/matrix.hpp>
#include <amal/vector.hpp>
#include <cassert>
#include <cstdio>

using namespace amal;

void test_compare()
{
    {
        vec3 a(1.0f, 2.0f, 3.0f);
        vec3 b(1.0f, 2.0f, 3.0f);
        bvec3 r = equal(a, b);
        assert(r.x && r.y && r.z);
        printf("equal(vec3) passed\n");
    }

    {
        vec4 a(1.0f, 2.0f, 3.0f, 4.0f);
        vec4 b(1.0f, 0.0f, 3.0f, 5.0f);
        bvec4 r = equal(a, b);
        assert(r.x && !r.y && r.z && !r.w);
        printf("equal(vec4 partial) passed\n");
    }

    {
        vec2 a(1.0f, 2.0f);
        vec2 b(1.0f, 3.0f);
        bvec2 r = not_equal(a, b);
        assert(!r.x && r.y);
        printf("not_equal(vec2) passed\n");
    }

    {
        mat3 a(1.0f);
        mat3 b(1.0f);
        bvec3 r = equal(a, b);
        assert(r.x && r.y && r.z);
        printf("matrix_equal(mat3 identity) passed\n");
    }

    {
        mat2 a(vec2(1.0f, 2.0f), vec2(3.0f, 4.0f));
        mat2 b(vec2(1.0f, 2.1f), vec2(3.0f, 4.0f));
        bvec2 r = equal(a, b);
        assert(!r.x && r.y);
        printf("matrix_equal(mat2 partial) passed\n");
    }

    {
        mat4 a(1.0f);
        mat4 b(1.0f);
        bvec4 r = not_equal(a, b);
        assert(!r.x && !r.y && !r.z && !r.w);
        printf("matrix_not_equal(mat4 identical) passed\n");
    }

    {
        mat3 a(vec3(1.0f, 0.0f, 0.0f),
               vec3(0.0f, 1.0f, 0.0f),
               vec3(0.0f, 0.0f, 1.0f));
        mat3 b = a;
        b[1].y = 0.0f;
        bvec3 r = not_equal(a, b);
        assert(!r.x && r.y && !r.z);
        printf("matrix_not_equal(mat3 partial) passed\n");
    }

    printf("test_compare passed!\n");
}
