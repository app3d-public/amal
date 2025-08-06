#include <amal/matrix.hpp>
#include <cstdio>

using namespace amal;

void test_matrix()
{
    {
        mat2 m(vec2(1.0f, 2.0f), vec2(3.0f, 4.0f));
        mat2 t = transpose(m);
        assert(t[0].x == 1.0f && t[0].y == 3.0f);
        assert(t[1].x == 2.0f && t[1].y == 4.0f);
        printf("transpose passed\n");
    }

    {
        mat2 m(vec2(4.0f, 7.0f), vec2(2.0f, 6.0f));
        float d = determinant(m);
        assert(abs(d - 10.0f) < 1e-5f);
        printf("determinant passed\n");
    }

    {
        mat2 m(vec2(4.0f, 7.0f), vec2(2.0f, 6.0f));
        mat2 inv = inverse_matrix(m);
        mat2 identity = m * inv;

        assert(abs(identity[0].x - 1.0f) < 1e-3f);
        assert(abs(identity[0].y - 0.0f) < 1e-3f);
        assert(abs(identity[1].x - 0.0f) < 1e-3f);
        assert(abs(identity[1].y - 1.0f) < 1e-3f);
        printf("inverse passed\n");
    }

    printf("test_matrix passed!\n");
}
