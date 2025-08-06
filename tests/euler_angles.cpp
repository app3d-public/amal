#include <amal/euler_angles.hpp>
#include <cstdio>

using namespace amal;

void test_euler_angles()
{
    {
        float ax = 0.37f;
        float ay = -1.15f;
        float az = 2.47f;

        mat4 m1 = euler_angle_xyz<float, false>(ax, ay, az);

        float rx, ry, rz;
        extract_euler_angle_xyz<float, false>(m1, rx, ry, rz);

        mat4 m2 = euler_angle_xyz<float, false>(rx, ry, rz);

        for (int c = 0; c < 4; ++c)
        {
            for (int r = 0; r < 4; ++r)
            {
                float a = m1[c][r];
                float b = m2[c][r];
                assert(abs(a - b) < 1e-4f);
            }
        }

        printf("euler_angle_xyz/extract_xyz matrix comparison passed\n");
    }

    {
        float ax = -0.9f;
        float ay = 0.4f;
        float az = 1.2f;

        mat4 m1 = euler_angle_xyz<float, false>(ax, ay, az);
        float rx, ry, rz;
        extract_euler_angle_xyz<float, false>(m1, rx, ry, rz);
        mat4 m2 = euler_angle_xyz<float, false>(rx, ry, rz);

        for (int c = 0; c < 4; ++c)
        {
            for (int r = 0; r < 4; ++r)
            {
                float a = m1[c][r];
                float b = m2[c][r];
                assert(abs(a - b) < 1e-4f);
            }
        }

        printf("euler_angle_xyz/extract_xyz (2nd case) passed\n");
    }

    printf("test_euler_angles passed!\n");
}
