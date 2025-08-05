# App3D Math Library

**Amal** is a C++ header-only math library, developed as part of the App3D project. It is designed specifically for high-performance real-time applications, such as 3D rendering and interactive tools.

**_NOTE:_** Amal is not a general-purpose or scientific math library.
All algorithms are optimized for speed over numerical precision.  
Do not use Amal for scientific computing, high-precision modeling, or anywhere strict numerical guarantees are required.

Structures are aligned in accordance with std140.\
The amal is not a fully constexpr-ready library by default settings in accordance with using SIMD and FMA instructions, which can't be calculated at compile time. By default settings, you can access constexpr base constructors. Functions and operators are not supported. To enable the constexpr functions, you need to disable FMA and SIMD using flags described below.

## Limitations
- SIMD limitations:
    - ARM/NEON is not supported.
    - The minimal instruction sets is SSE2.
    - SIMD on `bool`, `half` and `uint` types is not supported.
- Half-precision arithmetic is not IEEE 754-compliant

## Dependencies:
- [acbt](https://github.com/app3d/acbt) (Testing Only)
- [jinja2 Python Package](https://pypi.org/project/Jinja2/)

## Benchmark

<details>
<summary>amal/glm comparison</summary>

OS: Microsoft Windows 11 10.0.26100\
CPU: Intel(R) Core(TM) Ultra 9 185H\
Types size: packed (default)\
glm flags:
- `GLM_ENABLE_CXX_20`
- `GLM_FORCE_INTRINSICS`

amal flags: none (default)

| Benchmark                               |   Time (ns) |   Bandwidth (GiB/s) |   cps_avg (M) |   cps_max (M) |   cps_min (M) |
|:----------------------------------------|------------:|--------------------:|--------------:|--------------:|--------------:|
| BM_glm_vec3_add             |       14804 |               0.755 |        70.808 |        80     |         1.134 |
| BM_amal_vec3_add            |        8974 |               1.245 |       117.791 |       133.333 |         1.655 |
| BM_glm_vec3_mul_scalar      |       15257 |               0.733 |        68.209 |        78.125 |         0.657 |
| BM_amal_vec3_mul_scalar     |        8143 |               1.373 |       128.309 |       144.928 |         2.706 |
| BM_glm_vec3_dot             |       17533 |               0.637 |        59.315 |        67.568 |         1.351 |
| BM_amal_vec3_dot            |        4087 |               2.735 |       252.177 |       303.030 |         3.394 |
| BM_glm_vec3_normalize       |       40701 |               0.275 |        25.355 |        27.855 |         0.875 |
| BM_amal_vec3_normalize      |       12529 |               0.892 |        83.066 |        97.087 |         2.355 |
| BM_glm_vec3_cross           |        9081 |               1.231 |       114.932 |       129.870 |         2.203 |
| BM_amal_vec3_cross          |        7874 |               1.419 |       132.273 |       149.254 |         0.687 |
| BM_glm_vec4_add             |       15606 |               0.955 |        65.826 |        75.188 |         4.728 |
| BM_amal_vec4_add            |        9321 |               1.599 |       110.555 |       120.482 |         2.355 |
| BM_glm_vec4_mul_scalar      |       15785 |               0.944 |        65.976 |        75.188 |         1.094 |
| BM_amal_vec4_mul_scalar     |        8684 |               1.716 |       119.388 |       131.579 |         1.732 |
| BM_glm_vec4_dot             |       18201 |               0.838 |        56.684 |        63.291 |         2.690 |
| BM_amal_vec4_dot            |        5138 |               2.900 |       202.256 |       256.410 |         2.470 |
| BM_glm_vec4_normalize       |       45298 |               0.337 |        22.825 |        25.445 |         1.071 |
| BM_amal_vec4_normalize      |        7193 |               2.072 |       146.167 |       178.571 |         1.241 |
| BM_glm_mat3_mat_add         |       62716 |               0.535 |        16.374 |        18.519 |         2.473 |
| BM_amal_mat3_mat_add        |       28724 |               1.167 |        35.945 |        39.682 |         2.311 |
| BM_glm_mat3_mat_mul_scalar  |       55663 |               0.603 |        18.346 |        20.534 |         2.345 |
| BM_amal_mat3_mat_mul_scalar |       27275 |               1.229 |        37.783 |        41.494 |         1.199 |
| BM_glm_mat3_mat_mul_vec     |      134977 |               0.248 |         7.638 |         8.795 |         1.253 |
| BM_amal_mat3_mat_mul_vec    |       30802 |               1.088 |        33.581 |        42.918 |         1.677 |
| BM_glm_mat3_mat_mul_mat     |      209121 |               0.164 |         4.937 |         5.587 |         1.033 |
| BM_amal_mat3_mat_mul_mat    |      168716 |               0.203 |         6.130 |         7.994 |         1.007 |
| BM_glm_mat3_mat_transpose   |      116172 |               0.289 |         8.764 |        10.428 |         1.074 |
| BM_amal_mat3_mat_transpose  |       30984 |               1.082 |        33.350 |        42.194 |         1.850 |
| BM_glm_mat4_mat_add         |       79906 |               0.746 |        12.915 |        14.749 |         1.551 |
| BM_amal_mat4_mat_add        |       37269 |               1.599 |        27.800 |        30.212 |         2.444 |
| BM_glm_mat4_mat_mul_scalar  |       78847 |               0.756 |        13.057 |        14.706 |         0.822 |
| BM_amal_mat4_mat_mul_scalar |       35909 |               1.660 |        28.741 |        30.864 |         0.784 |
| BM_glm_mat4_mat_mul_vec     |      116827 |               0.522 |         8.840 |        10.081 |         1.771 |
| BM_amal_mat4_mat_mul_vec    |       10268 |               5.805 |        99.521 |       112.360 |         1.339 |
| BM_glm_mat4_mat_mul_mat     |      424352 |               0.144 |         2.435 |         2.755 |         0.822 |
| BM_amal_mat4_mat_mul_mat    |       28936 |               2.060 |        35.493 |        40.161 |         2.141 |
| BM_glm_mat4_mat_transpose   |      203053 |               0.295 |         5.077 |         6.177 |         1.430 |
| BM_amal_mat4_mat_transpose  |       17040 |               3.498 |        60.892 |        69.930 |         1.549 |
| BM_glm_mat4_translate       |      114098 |               0.522 |         8.975 |        10.132 |         1.217 |
| BM_amal_mat4_translate      |       67315 |               0.887 |        15.341 |        16.892 |         2.488 |

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
#### Flags:
- `AMAL_FORCE_HALF_PRECISION`: Force half-precision arithmetic instead of float32
- `AMAL_FORCE_ALIGNED_TYPES`: Type aliases are forced to be aligned
- `AMAL_FMA_DISABLE`: Disable FMA instructions
- `AMAL_SIMD_DISABLE`: Disable SIMD instructions
- `AMAL_RIGHT_HANDED`: Use right-handed coordinate system. Default is left-handed

#### Options:
- `BUILD_TESTS`: Enable testing
- `ENABLE_COVERAGE`: Enable code coverage