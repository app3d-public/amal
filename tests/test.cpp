#include <amal/matrix.hpp>
#include <cstdio>
#include <typeinfo>

int main()
{
    amal::mat2x2 m1{1, 2, 3, 4};
    printf("simd type: %s\n", typeid(typename amal::mat2x2_aligned::simd_type::value_type).name());
    amal::mat2x2 m2{5, 6, 7, 8};
    auto m3 = m1 * m2;
    printf("%f %f %f %f\n", m3[0][0], m3[0][1], m3[1][0], m3[1][1]);
    return 0;
}
