#include <amal/geometric.hpp>
#include <cstdio>

using namespace amal;

void test_geometric()
{
    {
        vec3 v(3.0f, 4.0f, 0.0f);
        float l = length(v);
        assert(abs(l - 5.0f) < 1e-5f);
        printf("length passed\n");
    }

    {
        vec3 a(1.0f, 2.0f, 3.0f);
        vec3 b(4.0f, 6.0f, 3.0f);
        float d = distance(a, b);
        assert(abs(d - 5.0f) < 1e-5f);
        printf("distance passed\n");
    }

    {
        vec3 a(1.0f, 2.0f, 3.0f);
        vec3 b(4.0f, -5.0f, 6.0f);
        float d = dot(a, b);
        assert(abs(d - 12.0f) < 1e-5f);
        printf("dot passed\n");
    }

    {
        vec3 a(1.0f, 0.0f, 0.0f);
        vec3 b(0.0f, 1.0f, 0.0f);
        vec3 r = cross(a, b);
        assert(abs(r.x - 0.0f) < 1e-5f);
        assert(abs(r.y - 0.0f) < 1e-5f);
        assert(abs(r.z - 1.0f) < 1e-5f);
        printf("cross passed\n");
    }

    {
        vec3 v(0.0f, 3.0f, 4.0f);
        vec3 r = normalize(v);
        assert(abs(r.x - 0.0f) < 1e-5f);
        assert(abs(r.y - 0.6f) < 1e-5f);
        assert(abs(r.z - 0.8f) < 1e-5f);
        printf("normalize passed\n");
    }

    {
        vec3 n(0.0f, 0.0f, 1.0f);
        vec3 i(0.0f, 0.0f, -1.0f);
        vec3 nref(0.0f, 0.0f, 1.0f);
        vec3 r = face_forward(n, i, nref);
        assert(abs(r.z - 1.0f) < 1e-5f);
        printf("face_forward passed\n");
    }

    {
        vec3 n(0.0f, 1.0f, 0.0f);
        vec3 i(1.0f, -1.0f, 0.0f);
        vec3 r = reflect(i, n);
        assert(abs(r.x - 1.0f) < 1e-5f);
        assert(abs(r.y - 1.0f) < 1e-5f);
        assert(abs(r.z - 0.0f) < 1e-5f);
        printf("reflect passed\n");
    }

    {
        vec3 n(0.0f, 1.0f, 0.0f);
        vec3 i(0.0f, -1.0f, 0.0f);
        float eta = 0.5f;
        vec3 r = refract(i, n, eta);
        assert(abs(r.x) < 1e-5f);
        assert(abs(r.y + 1.0f) < 1e-5f);
        assert(abs(r.z) < 1e-5f);
        printf("refract (normal incidence) passed\n");
    }

    printf("test_geometric passed!\n");
}
