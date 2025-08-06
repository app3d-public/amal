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

OS: Microsoft Windows 10.0.26100\
CPU: Intel(R) Core(TM) Ultra 9 185H\
Types size: packed (default)\
glm flags:
- `GLM_ENABLE_CXX_20`
- `GLM_FORCE_INTRINSICS`

amal flags: none (default)

| Benchmark                               |   Time (ns) |   Bandwidth (GiB/s) |   cps_avg (M) |   cps_max (M) |   cps_min (M) |
|:----------------------------------------|------------:|--------------------:|--------------:|--------------:|--------------:|
| BM_glm_vec3_add            |      14379 |             0.7772 |       70.7756 |       80.0000 |        0.9924 |
| BM_amal_vec3_add           |       8754 |             1.2767 |      119.3210 |      133.3330 |        0.5678 |
| BM_glm_vec3_mul_scalar     |      15199 |             0.7353 |       68.1126 |       78.1250 |        1.6176 |
| BM_amal_vec3_mul_scalar    |       8133 |             1.3741 |      129.2630 |      144.9280 |        1.9286 |
| BM_glm_vec3_dot            |      17368 |             0.6435 |       60.0295 |       67.5676 |        0.8031 |
| BM_amal_vec3_dot           |       4236 |             2.6381 |      248.4100 |      303.0300 |        5.6818 |
| BM_glm_vec3_normalize      |      41741 |             0.2677 |       24.9370 |       27.8552 |        0.9571 |
| BM_amal_vec3_normalize     |      12069 |             0.9260 |       85.5178 |       98.0392 |        2.1758 |
| BM_glm_vec3_cross          |       8959 |             1.2474 |      115.9240 |      129.8700 |        2.6911 |
| BM_amal_vec3_cross         |       7728 |             1.4461 |      132.4150 |      149.2540 |        1.1120 |
| BM_glm_vec4_add            |      16081 |             0.9266 |       65.0168 |       75.1880 |        1.4861 |
| BM_amal_vec4_add           |       9549 |             1.5605 |      109.6090 |      120.4820 |        1.6633 |
| BM_glm_vec4_mul_scalar     |      15999 |             0.9314 |       65.2765 |       74.6269 |        2.0338 |
| BM_amal_vec4_mul_scalar    |       8524 |             1.7482 |      120.2810 |      131.5790 |        3.6603 |
| BM_glm_vec4_dot            |      18075 |             0.8244 |       56.7989 |       63.2911 |        1.7259 |
| BM_amal_vec4_dot           |       5144 |             2.8968 |      202.7530 |      256.4100 |        1.6319 |
| BM_glm_vec4_normalize      |      45676 |             0.3262 |       22.7955 |       25.5102 |        3.0609 |
| BM_amal_vec4_normalize     |       6861 |             2.1718 |      151.1290 |      178.5710 |        2.9904 |
| BM_glm_mat3_mat_add        |      62754 |             0.5343 |       16.5438 |       18.5529 |        1.9759 |
| BM_amal_mat3_mat_add       |      29376 |             1.1413 |       35.6115 |       39.6825 |        1.9693 |
| BM_glm_mat3_mat_mul_scalar |      56283 |             0.5957 |       18.4367 |       20.5761 |        1.3275 |
| BM_amal_mat3_mat_mul_scalar|      27358 |             1.2255 |       37.6666 |       41.4938 |        2.8209 |
| BM_glm_mat3_mat_mul_vec    |     134919 |             0.2485 |        7.7070 |        8.7951 |        0.9537 |
| BM_amal_mat3_mat_mul_vec   |      35424 |             0.9465 |       29.2337 |       36.3636 |        1.2819 |
| BM_glm_mat3_mat_mul_mat    |     205935 |             0.1628 |        5.0027 |        5.5866 |        1.8077 |
| BM_amal_mat3_mat_mul_mat   |     187013 |             0.1793 |        5.4869 |        7.0175 |        1.4443 |
| BM_glm_mat3_mat_transpose  |     121522 |             0.2759 |        8.4939 |       10.1317 |        1.4843 |
| BM_amal_mat3_mat_transpose |      35072 |             0.9560 |       29.7086 |       36.9004 |        2.3691 |
| BM_glm_mat3_inverse        |     388949 |             0.0862 |        2.6307 |        3.1250 |        0.7598 |
| BM_amal_mat3_inverse       |     128956 |             0.2600 |        7.8872 |        9.9404 |        2.1268 |
| BM_glm_mat4_mat_add        |      78899 |             0.7555 |       12.9791 |       14.5349 |        1.0396 |
| BM_amal_mat4_mat_add       |      38167 |             1.5617 |       27.1162 |       29.5858 |        3.6140 |
| BM_glm_mat4_mat_mul_scalar |      83260 |             0.7159 |       12.5723 |       14.7493 |        1.9944 |
| BM_amal_mat4_mat_mul_scalar|      36613 |             1.6280 |       28.2258 |       30.4878 |        3.1017 |
| BM_glm_mat4_mat_mul_vec    |     117386 |             0.5078 |        8.8260 |        9.9404 |        0.7092 |
| BM_amal_mat4_mat_mul_vec   |      10452 |             5.7027 |       99.4458 |      112.3600 |        3.4674 |
| BM_glm_mat4_mat_mul_mat    |     401056 |             0.1486 |        2.5858 |        2.8977 |        0.8968 |
| BM_amal_mat4_mat_mul_mat   |      29076 |             2.0499 |       35.8042 |       40.3226 |        1.2721 |
| BM_glm_mat4_mat_transpose  |     205584 |             0.2899 |        5.0495 |        5.9067 |        0.9642 |
| BM_amal_mat4_mat_transpose |      17102 |             3.4853 |       60.6989 |       68.0272 |        2.3026 |
| BM_glm_mat4_translate      |     114782 |             0.5193 |        9.0026 |       10.1112 |        2.2502 |
| BM_amal_mat4_translate     |      58879 |             1.0123 |       17.7301 |       19.6850 |        0.6772 |
| BM_glm_mat4_inverse        |     933741 |             0.0638 |        1.0825 |        1.2145 |        0.4827 |
| BM_amal_mat4_inverse       |      49699 |             1.1993 |       20.6486 |       22.6757 |        1.9205 |

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
- `AMAL_CLIP_SPACE_LH`: Use left-handed coordinate system. Default is right-handed
- `AMAL_CLIP_SPACE_NO`: Use Negative One to One clip space. Default is Zero to One
- `AMAL_NO_GLOBAL_ALIASES`: Disable global type aliases

#### Options:
- `BUILD_TESTS`: Enable testing
- `ENABLE_COVERAGE`: Enable code coverage