#include <amal/exponential.hpp>
#include <cstdio>

using namespace amal;

void test_exponential()
{
    {
        vec3 v(1.0f, 4.0f, 9.0f);
        vec3 r = sqrt(v);
        assert(abs(r.x - 1.0f) < 1e-5f);
        assert(abs(r.y - 2.0f) < 1e-5f);
        assert(abs(r.z - 3.0f) < 1e-5f);
        printf("sqrt passed\n");
    }

    {
        vec3 v(1.0f, 4.0f, 9.0f);
        vec3 r = inverse_sqrt(v);
        assert(abs(r.x - 1.0f) < 1e-5f);
        assert(abs(r.y - 0.5f) < 1e-5f);
        assert(abs(r.z - 1.0f / 3.0f) < 1e-5f);
        printf("inversesqrt passed\n");
    }

    {
        vec3 v(1.0f, 2.0f, 3.0f);
        vec3 r = exp(v);
        assert(abs(r.x - exp(1.0f)) < 1e-4f);
        assert(abs(r.y - exp(2.0f)) < 1e-4f);
        assert(abs(r.z - exp(3.0f)) < 1e-4f);
        printf("exp passed\n");
    }

    {
        vec3 v(1.0f, 2.0f, 4.0f);
        vec3 r = log(v);
        assert(abs(r.x - log(1.0f)) < 1e-5f);
        assert(abs(r.y - log(2.0f)) < 1e-5f);
        assert(abs(r.z - log(4.0f)) < 1e-5f);
        printf("log passed\n");
    }

    {
        vec3 v(1.0f, 2.0f, 3.0f);
        vec3 r = exp2(v);
        assert(abs(r.x - 2.0f) < 1e-4f);
        assert(abs(r.y - 4.0f) < 1e-4f);
        assert(abs(r.z - 8.0f) < 1e-4f);
        printf("exp2 passed\n");
    }

    {
        vec3 v(1.0f, 2.0f, 4.0f);
        vec3 r = log2(v);
        assert(abs(r.x - 0.0f) < 1e-5f);
        assert(abs(r.y - 1.0f) < 1e-5f);
        assert(abs(r.z - 2.0f) < 1e-5f);
        printf("log2 passed\n");
    }

    {
        vec3 base(2.0f, 3.0f, 10.0f);
        vec3 exponent(3.0f, 2.0f, 1.0f);
        vec3 r = pow(base, exponent);
        assert(abs(r.x - 8.0f) < 1e-4f);
        assert(abs(r.y - 9.0f) < 1e-4f);
        assert(abs(r.z - 10.0f) < 1e-4f);
        printf("pow passed\n");
    }

    printf("test_exponential passed!\n");
}
