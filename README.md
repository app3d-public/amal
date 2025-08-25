# App3D Math Library

**Amal** is a C++ header-only math library, developed as part of the App3D project. It is designed specifically for high-performance real-time applications, such as 3D rendering and interactive tools.

> [!NOTE]
> Amal is not a general-purpose or scientific math library.
> All algorithms are optimized for speed over numerical precision.  
> Do not use Amal for scientific computing, high-precision modeling, or anywhere strict numerical guarantees are required.

Structures are aligned in accordance with std140.\
The amal is not a fully constexpr-ready library by default settings in accordance with using SIMD and FMA instructions, which can't be calculated at compile time. By default settings, you can access constexpr base constructors. Functions and operators are not supported. To enable the constexpr functions, you need to disable FMA and SIMD using flags described below.

## Limitations
- SIMD limitations:
    - ARM/NEON is not supported.
    - The minimal instruction sets is SSE2.
    - SIMD on `bool`, `half` and `uint` types is not supported.
- Half-precision arithmetic is not IEEE 754-compliant

## Dependencies:
- [acbt](https://git.homedatasrv.ru/app3d/acbt)
- [jinja2 Python Package](https://pypi.org/project/Jinja2/)

## Benchmark

<details>
<summary>amal/glm comparison</summary>

OS: Microsoft Windows 10.0.26100\
CPU: Intel(R) Core(TM) Ultra 9 185H\
Types size: packed (default)\
glm flags:
- `GLM_ENABLE_CXX_20`
- `GLM_FORCE_INTRINSICS`

amal flags: none (default)

| Benchmark | Bandwidth (MiB/s) | avg (M) | min (M) | max (M) |
|-----------|-------|---------|---------|---------|
| BM_glm_vec3_add              | 771.616     | 67.425  | 45.327  | 75.241  |
| BM_amal_vec3_add             | 1283.06     | 112.116 | 79.481  | 123.873 |
| BM_glm_vec3_mul_scalar       | 766.925     | 67.015  | 47.356  | 73.586  |
| BM_amal_vec3_mul_scalar      | 1387.01     | 121.199 | 70.189  | 140.083 |
| BM_glm_vec3_dot              | 636.984     | 55.661  | 36.614  | 63.068  |
| BM_amal_vec3_dot             | 2802.34     | 244.873 | 118.581 | 271.861 |
| BM_glm_vec3_normalize        | 271.45      | 23.720  | 12.224  | 26.682  |
| BM_amal_vec3_normalize       | 914.421     | 79.903  | 45.708  | 88.700  |
| BM_glm_vec3_cross            | 1258.46     | 109.966 | 58.241  | 124.444 |
| BM_amal_vec3_cross           | 1459.45     | 127.529 | 85.467  | 141.000 |
| BM_glm_vec4_add              | 939.148     | 61.548  | 21.338  | 68.758  |
| BM_amal_vec4_add             | 1560.72     | 102.283 | 62.117  | 116.070 |
| BM_glm_vec4_mul_scalar       | 953.040     | 62.458  | 37.240  | 69.028  |
| BM_amal_vec4_mul_scalar      | 1704.75     | 111.722 | 58.317  | 125.275 |
| BM_glm_vec4_dot              | 827.122     | 54.206  | 31.559  | 59.183  |
| BM_amal_vec4_dot             | 3106.97     | 203.618 | 120.995 | 230.727 |
| BM_glm_vec4_normalize        | 320.867     | 21.028  | 10.499  | 24.086  |
| BM_amal_vec4_normalize       | 2167.96     | 142.079 | 93.828  | 162.673 |
| BM_glm_mat3_mat_add          | 528.763     | 15.401  | 9.708   | 17.123  |
| BM_amal_mat3_mat_add         | 1161.57     | 33.833  | 18.944  | 37.670  |
| BM_glm_mat3_mat_mul_scalar   | 568.616     | 16.562  | 9.963   | 18.692  |
| BM_amal_mat3_mat_mul_scalar  | 1189.93     | 34.659  | 23.161  | 38.266  |
| BM_glm_mat3_mat_mul_vec      | 252.550     | 7.356   | 4.976   | 8.229   |
| BM_amal_mat3_mat_mul_vec     | 980.659     | 28.564  | 17.414  | 33.030  |
| BM_glm_mat3_mat_mul_mat      | 162.555     | 4.735   | 2.374   | 5.249   |
| BM_amal_mat3_mat_mul_mat     | 200.660     | 5.845   | 3.385   | 6.639   |
| BM_glm_mat3_mat_transpose    | 301.828     | 8.791   | 6.313   | 9.782   |
| BM_amal_mat3_mat_transpose   | 1048.97     | 30.553  | 18.849  | 34.288  |
| BM_glm_mat3_inverse          | 92.7561     | 2.702   | 1.854   | 3.028   |
| BM_amal_mat3_inverse         | 285.721     | 8.322   | 4.566   | 9.247   |
| BM_glm_mat4_mat_add          | 714.624     | 11.708  | 6.334   | 13.605  |
| BM_amal_mat4_mat_add         | 1631.48     | 26.730  | 16.922  | 29.744  |
| BM_glm_mat4_mat_mul_scalar   | 724.140     | 11.864  | 6.939   | 13.355  |
| BM_amal_mat4_mat_mul_scalar  | 1630.42     | 26.713  | 17.678  | 29.987  |
| BM_glm_mat4_mat_mul_vec      | 494.231     | 8.097   | 4.449   | 9.265   |
| BM_amal_mat4_mat_mul_vec     | 5518.89     | 90.421  | 56.932  | 102.866 |
| BM_glm_mat4_mat_mul_mat      | 148.886     | 2.439   | 1.315   | 2.725   |
| BM_amal_mat4_mat_mul_mat     | 2080.79     | 34.092  | 26.447  | 37.365  |
| BM_glm_mat4_mat_transpose    | 294.063     | 4.818   | 3.098   | 5.651   |
| BM_amal_mat4_mat_transpose   | 3338.26     | 54.694  | 37.848  | 61.600  |
| BM_glm_mat4_translate        | 506.278     | 8.295   | 4.414   | 9.351   |
| BM_amal_mat4_translate       | 1029.54     | 16.868  | 9.373   | 18.756  |
| BM_glm_mat4_inverse          | 61.7332     | 1.011   | 0.613   | 1.141   |
| BM_amal_mat4_inverse         | 1202.08     | 19.695  | 14.687  | 21.641  |

</details>

## Building
### Supported compilers:
- GNU GCC
- Clang

### Requirements:
- Cmake (3.17 or newer)
- Python3

Amal is provided as CMake interface library.
### Cmake options:
#### Defines:
- `AMAL_FORCE_HALF_PRECISION`: Force half-precision arithmetic instead of float32
- `AMAL_FORCE_ALIGNED_TYPES`: Type aliases are forced to be aligned
- `AMAL_FMA_DISABLE`: Disable FMA instructions
- `AMAL_SIMD_DISABLE`: Disable SIMD instructions
- `AMAL_CLIP_SPACE_LH`: Use left-handed coordinate system. Default is right-handed
- `AMAL_CLIP_SPACE_NO`: Use Negative One to One clip space. Default is Zero to One
- `AMAL_NO_GLOBAL_ALIASES`: Disable global type aliases

#### Options:
- `BUILD_TESTS`: Enable testing
- `ENABLE_COVERAGE`: Enable code coverage

## License
This project is licensed under the [MIT License](LICENSE).

## Contacts
For any questions or feedback, you can reach out via [email](mailto:wusikijeronii@gmail.com) or open a new issue.