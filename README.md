# App3D Math Library

**Amal** is a header-only math library, developed as part of the App3D project. It is designed specifically for high-performance real-time applications, such as 3D rendering and interactive tools.

**_NOTE:_** Amal is not a general-purpose or scientific math library.
All algorithms are optimized for speed over numerical precision.  
Do not use Amal for scientific computing, high-precision modeling, or anywhere strict numerical guarantees are required.

## Limitations
- SIMD limitations:
    - ARM/NEON is not supported.
    - The minimal instruction sets is SSE2.
    - Not supported with MSVC compiler.
    - SIMD on `bool`, `half` and `uint` types is not supported.
- Half-precision arithmetic is not IEEE 754-compliant

## Building
Amal is integrated via CMake.\
Preprocessor flags:
- `AMAL_PRECISION_MEDIUM`: The default floating type becomes `amal::half`
- `AMAL_PRECISION_HIGH` (default): The default floating type becomes `float`