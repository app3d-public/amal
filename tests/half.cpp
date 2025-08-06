#include <amal/half.hpp>
#include <cassert>
#include <cstdio>

using namespace amal;

void test_half()
{
    half h1 = 1.5f;
    half h2 = 2.0f;
    half h3 = static_cast<half>(3.0);

    float f1 = static_cast<float>(h1);
    assert(f1 > 1.4f && f1 < 1.6f);

    half sum = h1 + h2;
    half mul = h1 * h3;

    assert(static_cast<float>(sum) > 3.4f && static_cast<float>(sum) < 3.6f);
    assert(static_cast<float>(mul) > 4.4f && static_cast<float>(mul) < 4.6f);

    assert(h2 > h1);
    assert(h1 != h3);
    assert(h1 < h3);

    ++h1;
    assert(h1 > half(2.4f) && h1 < half(2.6f));

    h1 += half(1.0f);
    assert(h1 > half(3.4f) && h1 < half(3.6f));

    half div = h3 / h2;
    assert(static_cast<float>(div) > 1.4f && static_cast<float>(div) < 1.6f);

    assert(h2 == 2.0f);
    assert(h1 > 3.0f);

    printf("All half tests passed.\n");
}