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
    - SIMD on `bool`, `half` and `uint` types is not supported.
- Half-precision arithmetic is not IEEE 754-compliant

## Dependencies:
- [acbt](https://github.com/app3d-public/acbt)
- [acul](https://github.com/app3d-public/acul) - Optional. Used for acul integration only
- [jinja2 Python Package](https://pypi.org/project/Jinja2/)

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
- `AMAL_INTEGRATE_ACUL`: Specialize acul functions
- `BUILD_TESTS`: Enable testing
- `ENABLE_COVERAGE`: Enable code coverage

## License
This project is licensed under the [MIT License](LICENSE).

## Contacts
For any questions or feedback, you can reach out via [email](mailto:wusikijeronii@gmail.com) or open a new issue.