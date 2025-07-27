#include <cstdio>
#include <glm/glm.hpp>
#include <typeinfo>
#include "../include/amal/half.hpp"
#include "../include/amal/vector.hpp"


// #if GLM_HAS_ALIGNOF && (GLM_LANG & GLM_LANG_CXXMS_FLAG) && (GLM_ARCH & GLM_ARCH_SIMD_BIT)
// #	define GLM_CONFIG_SIMD GLM_ENABLE
// #else
// #	define GLM_CONFIG_SIMD GLM_DISABLE
// #endif

int main()
{

    amal::vec4 d{false};
    printf("simd type: %s\n", typeid(d.s).name());
    printf("__v4sf: %s\n", typeid(__v4sf).name());
    // // bool a[3];
    // printf("sizeof(d) = %zu\n", sizeof(amal::bvec2));
    // f16 d{15.0};
    // printf("sizeof(d) = %zu, %f\n", sizeof(d), static_cast<bool>(d));
    // amal::vec3 v1{4.0f};
    // printf("type: %s\n", typeid(v1.x).name());
    // amal::vec3 v2{2.0f, 2.0f, 2.0f};
    // amal::vec3 v3 = v1 / v2;
    // printf("%f %f %f\n", v3.data[0], v3.data[1], v3.data[2]);
    // amal::vec3 v4{1.0f, 2.0f, 3.0f};
    // v4 += v3;
    // printf("%f %f %f\n", v4.data[0], v4.data[1], v4.data[2]);
    // printf("sizeof(amal::vec3) = %zu\n", sizeof(amal::vec3));
    return 0;
}