#include <benchmark/benchmark.h>
#include <glm/glm.hpp>
#include "amal/geometric.hpp"
#include "amal/vector.hpp"

int main()
{
    amal::vec4 v{1.0f, 2.0f, 3.0f, 4.0f};
    amal::vec4 r = amal::normalize(v);
    printf("r.x = %f, r.y = %f, r.z = %f, r.w = %f\n", r.x, r.y, r.z, r.w);

    return 0;
}