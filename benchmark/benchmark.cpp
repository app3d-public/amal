// benchmark_math_ops.cpp

#include <algorithm>
#include <amal/geometric.hpp>
#include <amal/matrix.hpp>
#include <amal/vector.hpp>
#include <benchmark/benchmark.h>
#include <cfloat>
#include <chrono>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <random>
#include <vector>

constexpr int N = 1000;
static std::vector<glm::vec3> glm3_a(N), glm3_b(N), glm3_out(N);
static std::vector<amal::vec3> amal3_a(N), amal3_b(N), amal3_out(N);
static std::vector<glm::vec4> glm4_a(N), glm4_b(N), glm4_out(N);
static std::vector<amal::vec4> amal4_a(N), amal4_b(N), amal4_out(N);
static std::vector<glm::mat3> glm_mat3_a(N), glm_mat3_b(N), glm_mat3_out(N);
static std::vector<amal::mat3> amal_mat3_a(N), amal_mat3_b(N), amal_mat3_out(N);
static std::vector<glm::vec3> glm_mat3_v(N), glm_mat3_vout(N);
static std::vector<amal::vec3> amal_mat3_v(N), amal_mat3_vout(N);

static std::vector<glm::mat4> glm_mat4_a(N), glm_mat4_b(N), glm_mat4_out(N);
static std::vector<amal::mat4> amal_mat4_a(N), amal_mat4_b(N), amal_mat4_out(N);
static std::vector<glm::vec4> glm_mat4_v(N), glm_mat4_vout(N);
static std::vector<amal::vec4> amal_mat4_v(N), amal_mat4_vout(N);

static void InitData()
{
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(0.f, 1.f);
    for (int i = 0; i < N; ++i)
    {
        glm3_a[i] = glm::vec3(dist(rng), dist(rng), dist(rng));
        glm3_b[i] = glm::vec3(dist(rng), dist(rng), dist(rng));
        glm3_out[i] = glm::vec3(0);
        amal3_a[i] = amal::vec3(dist(rng), dist(rng), dist(rng));
        amal3_b[i] = amal::vec3(dist(rng), dist(rng), dist(rng));
        amal3_out[i] = amal::vec3(0);

        glm4_a[i] = glm::vec4(dist(rng), dist(rng), dist(rng), dist(rng));
        glm4_b[i] = glm::vec4(dist(rng), dist(rng), dist(rng), dist(rng));
        glm4_out[i] = glm::vec4(0);
        amal4_a[i] = amal::vec4(dist(rng), dist(rng), dist(rng), dist(rng));
        amal4_b[i] = amal::vec4(dist(rng), dist(rng), dist(rng), dist(rng));
        amal4_out[i] = amal::vec4(0);

        // mat3
        glm_mat3_a[i] = glm::mat3(1.0f);
        glm_mat3_b[i] = glm::mat3(2.0f);
        glm_mat3_out[i] = glm::mat3(0.0f);
        glm_mat3_v[i] = glm::vec3(1.0f);
        glm_mat3_vout[i] = glm::vec3(0.0f);

        amal_mat3_a[i] = amal::mat3(1.0f);
        amal_mat3_b[i] = amal::mat3(2.0f);
        amal_mat3_out[i] = amal::mat3(0.0f);
        amal_mat3_v[i] = amal::vec3(1.0f);
        amal_mat3_vout[i] = amal::vec3(0.0f);

        // mat4
        glm_mat4_a[i] = glm::mat4(1.0f);
        glm_mat4_b[i] = glm::mat4(2.0f);
        glm_mat4_out[i] = glm::mat4(0.0f);
        glm_mat4_v[i] = glm::vec4(1.0f);
        glm_mat4_vout[i] = glm::vec4(0.0f);

        amal_mat4_a[i] = amal::mat4(1.0f);
        amal_mat4_b[i] = amal::mat4(2.0f);
        amal_mat4_out[i] = amal::mat4(0.0f);
        amal_mat4_v[i] = amal::vec4(1.0f);
        amal_mat4_vout[i] = amal::vec4(0.0f);
    }
}

// Macro to report stats
#define REPORT_STATS(bytes_per_iter)                                                                             \
    do {                                                                                                         \
        double avg = sum_cps / state.iterations();                                                               \
        state.counters["cps_min"] = min_cps;                                                                     \
        state.counters["cps_avg"] = avg;                                                                         \
        state.counters["cps_max"] = max_cps;                                                                     \
        state.counters["MiB/s"] = benchmark::Counter(static_cast<double>((bytes_per_iter) * state.iterations()), \
                                                     benchmark::Counter::kIsRate, benchmark::Counter::kIs1024);  \
    } while (0)

#define RUN_BENCHMARK(C, FUNC)                                      \
    InitData();                                                     \
    double min_cps = DBL_MAX, max_cps = 0.0, sum_cps = 0.0;         \
    const int64_t bytes_per_iter = N * sizeof(C);                   \
    for (auto _ : state)                                            \
    {                                                               \
        auto t0 = std::chrono::high_resolution_clock::now();        \
        for (int i = 0; i < N; ++i) benchmark::DoNotOptimize(FUNC); \
        auto t1 = std::chrono::high_resolution_clock::now();        \
        double dt = std::chrono::duration<double>(t1 - t0).count(); \
        double cps = N / dt;                                        \
        sum_cps += cps;                                             \
        min_cps = std::min(min_cps, cps);                           \
        max_cps = std::max(max_cps, cps);                           \
        state.SetIterationTime(dt);                                 \
    }                                                               \
    REPORT_STATS(bytes_per_iter);

// vec3 tests

// --- vec3 + vec3 ---
#define FUNC_VEC_ADD(C) C##_out[i] = C##_a[i] + C##_b[i]
static void BM_glm_vec3_add(benchmark::State &state) { RUN_BENCHMARK(glm::vec3, FUNC_VEC_ADD(glm3)); }
static void BM_amal_vec3_add(benchmark::State &state) { RUN_BENCHMARK(amal::vec3, FUNC_VEC_ADD(amal3)); }

// --- vec3 * scalar ---
#define FUNC_VEC_MUL_SCALAR(C) C##_out[i] = C##_a[i] * 0.5f
static void BM_glm_vec3_mul_scalar(benchmark::State &state) { RUN_BENCHMARK(glm::vec3, FUNC_VEC_MUL_SCALAR(glm3)); }
static void BM_amal_vec3_mul_scalar(benchmark::State &state) { RUN_BENCHMARK(amal::vec3, FUNC_VEC_MUL_SCALAR(amal3)); }

// --- vec3 dot vec3 ---
#define FUNC_VEC_DOT(N, C) N::dot(C##_a[i], C##_b[i])
static void BM_glm_vec3_dot(benchmark::State &state) { RUN_BENCHMARK(glm::vec3, FUNC_VEC_DOT(glm, glm3)); }
static void BM_amal_vec3_dot(benchmark::State &state) { RUN_BENCHMARK(amal::vec3, FUNC_VEC_DOT(amal, amal3)); }

// --- vec3 normalize ---
#define FUNC_VEC_NORMALIZE(N, C) C##_out[i] = N::normalize(C##_a[i])
static void BM_glm_vec3_normalize(benchmark::State &state) { RUN_BENCHMARK(glm::vec3, FUNC_VEC_NORMALIZE(glm, glm3)); }
static void BM_amal_vec3_normalize(benchmark::State &state)
{
    RUN_BENCHMARK(amal::vec3, FUNC_VEC_NORMALIZE(amal, amal3));
}

// --- vec3 cross vec3 ---
#define FUNC_VEC_CROSS(N, C) C##_out[i] = N::cross(C##_a[i], C##_b[i])
static void BM_glm_vec3_cross(benchmark::State &state) { RUN_BENCHMARK(glm::vec3, FUNC_VEC_CROSS(glm, glm3)); }
static void BM_amal_vec3_cross(benchmark::State &state) { RUN_BENCHMARK(amal::vec3, FUNC_VEC_CROSS(amal, amal3)); }

// vec4 tests

// --- vec4 + vec4 ---
static void BM_glm_vec4_add(benchmark::State &state) { RUN_BENCHMARK(glm::vec4, FUNC_VEC_ADD(glm4)); }
static void BM_amal_vec4_add(benchmark::State &state) { RUN_BENCHMARK(amal::vec4, FUNC_VEC_ADD(amal4)); }

// --- vec4 * scalar ---
static void BM_glm_vec4_mul_scalar(benchmark::State &state) { RUN_BENCHMARK(glm::vec4, FUNC_VEC_MUL_SCALAR(glm4)); }
static void BM_amal_vec4_mul_scalar(benchmark::State &state) { RUN_BENCHMARK(amal::vec4, FUNC_VEC_MUL_SCALAR(amal4)); }

// --- vec4 dot vec4 ---
static void BM_glm_vec4_dot(benchmark::State &state) { RUN_BENCHMARK(glm::vec4, FUNC_VEC_DOT(glm, glm4)); }
static void BM_amal_vec4_dot(benchmark::State &state) { RUN_BENCHMARK(amal::vec4, FUNC_VEC_DOT(amal, amal4)); }

// --- vec4 normalize ---
static void BM_glm_vec4_normalize(benchmark::State &state) { RUN_BENCHMARK(glm::vec4, FUNC_VEC_NORMALIZE(glm, glm4)); }
static void BM_amal_vec4_normalize(benchmark::State &state)
{
    RUN_BENCHMARK(amal::vec4, FUNC_VEC_NORMALIZE(amal, amal4));
}

// mat3 tests

// --- mat3 + mat3 ---
static void BM_glm_mat3_mat_add(benchmark::State &state)
{
    RUN_BENCHMARK(glm::mat3, glm_mat3_out[i] = glm_mat3_a[i] + glm_mat3_b[i]);
}
static void BM_amal_mat3_mat_add(benchmark::State &state)
{
    RUN_BENCHMARK(amal::mat3, amal_mat3_out[i] = amal_mat3_a[i] + amal_mat3_b[i]);
}

// --- mat3 * scalar ---
static void BM_glm_mat3_mat_mul_scalar(benchmark::State &state)
{
    RUN_BENCHMARK(glm::mat3, glm_mat3_out[i] = glm_mat3_a[i] * 0.5f);
}
static void BM_amal_mat3_mat_mul_scalar(benchmark::State &state)
{
    RUN_BENCHMARK(amal::mat3, amal_mat3_out[i] = amal_mat3_a[i] * 0.5f);
}

// --- mat3 * vec3 ---
static void BM_glm_mat3_mat_mul_vec(benchmark::State &state)
{
    RUN_BENCHMARK(glm::mat3, glm_mat3_vout[i] = glm_mat3_a[i] * glm_mat3_v[i]);
}
static void BM_amal_mat3_mat_mul_vec(benchmark::State &state)
{
    RUN_BENCHMARK(amal::mat3, amal_mat3_vout[i] = amal_mat3_a[i] * amal_mat3_v[i]);
}

// --- mat3 * mat3 ---
static void BM_glm_mat3_mat_mul_mat(benchmark::State &state)
{
    RUN_BENCHMARK(glm::mat3, glm_mat3_out[i] = glm_mat3_a[i] * glm_mat3_b[i]);
}
static void BM_amal_mat3_mat_mul_mat(benchmark::State &state)
{
    RUN_BENCHMARK(amal::mat3, amal_mat3_out[i] = amal_mat3_a[i] * amal_mat3_b[i]);
}

// --- mat3 transpose ---
static void BM_glm_mat3_mat_transpose(benchmark::State &state)
{
    RUN_BENCHMARK(glm::mat3, glm_mat3_out[i] = glm::transpose(glm_mat3_a[i]));
}
static void BM_amal_mat3_mat_transpose(benchmark::State &state)
{
    RUN_BENCHMARK(amal::mat3, amal_mat3_out[i] = amal::transpose(amal_mat3_a[i]));
}

// mat4 tests

// --- mat4 + mat4 ---
static void BM_glm_mat4_mat_add(benchmark::State &state)
{
    RUN_BENCHMARK(glm::mat4, glm_mat4_out[i] = glm_mat4_a[i] + glm_mat4_b[i]);
}
static void BM_amal_mat4_mat_add(benchmark::State &state)
{
    RUN_BENCHMARK(amal::mat4, amal_mat4_out[i] = amal_mat4_a[i] + amal_mat4_b[i]);
}

// --- mat4 * scalar ---
static void BM_glm_mat4_mat_mul_scalar(benchmark::State &state)
{
    RUN_BENCHMARK(glm::mat4, glm_mat4_out[i] = glm_mat4_a[i] * 0.5f);
}
static void BM_amal_mat4_mat_mul_scalar(benchmark::State &state)
{
    RUN_BENCHMARK(amal::mat4, amal_mat4_out[i] = amal_mat4_a[i] * 0.5f);
}

// --- mat4 * vec4 ---
static void BM_glm_mat4_mat_mul_vec(benchmark::State &state)
{
    RUN_BENCHMARK(glm::mat4, glm_mat4_vout[i] = glm_mat4_a[i] * glm_mat4_v[i]);
}
static void BM_amal_mat4_mat_mul_vec(benchmark::State &state)
{
    RUN_BENCHMARK(amal::mat4, amal_mat4_vout[i] = amal_mat4_a[i] * amal_mat4_v[i]);
}

// --- mat4 * mat4 ---
static void BM_glm_mat4_mat_mul_mat(benchmark::State &state)
{
    RUN_BENCHMARK(glm::mat4, glm_mat4_out[i] = glm_mat4_a[i] * glm_mat4_b[i]);
}
static void BM_amal_mat4_mat_mul_mat(benchmark::State &state)
{
    RUN_BENCHMARK(amal::mat4, amal_mat4_out[i] = amal_mat4_a[i] * amal_mat4_b[i]);
}

// --- mat4 transpose ---
static void BM_glm_mat4_mat_transpose(benchmark::State &state)
{
    RUN_BENCHMARK(glm::mat4, glm_mat4_out[i] = glm::transpose(glm_mat4_a[i]));
}
static void BM_amal_mat4_mat_transpose(benchmark::State &state)
{
    RUN_BENCHMARK(amal::mat4, amal_mat4_out[i] = amal::transpose(amal_mat4_a[i]));
}

// --- mat4 translate ---
static void BM_glm_mat4_translate(benchmark::State &state)
{
    RUN_BENCHMARK(glm::mat4, glm_mat4_out[i] = glm::translate(glm_mat4_a[i], glm3_a[i]));
}

static void BM_amal_mat4_translate(benchmark::State &state)
{
    RUN_BENCHMARK(amal::mat4, amal_mat4_out[i] = amal::translate(amal_mat4_a[i], amal3_a[i]));
}


// Register benchmarks
BENCHMARK(BM_glm_vec3_add)->UseManualTime();
BENCHMARK(BM_amal_vec3_add)->UseManualTime();
BENCHMARK(BM_glm_vec3_mul_scalar)->UseManualTime();
BENCHMARK(BM_amal_vec3_mul_scalar)->UseManualTime();
BENCHMARK(BM_glm_vec3_dot)->UseManualTime();
BENCHMARK(BM_amal_vec3_dot)->UseManualTime();
BENCHMARK(BM_glm_vec3_normalize)->UseManualTime();
BENCHMARK(BM_amal_vec3_normalize)->UseManualTime();
BENCHMARK(BM_glm_vec3_cross)->UseManualTime();
BENCHMARK(BM_amal_vec3_cross)->UseManualTime();

BENCHMARK(BM_glm_vec4_add)->UseManualTime();
BENCHMARK(BM_amal_vec4_add)->UseManualTime();
BENCHMARK(BM_glm_vec4_mul_scalar)->UseManualTime();
BENCHMARK(BM_amal_vec4_mul_scalar)->UseManualTime();
BENCHMARK(BM_glm_vec4_dot)->UseManualTime();
BENCHMARK(BM_amal_vec4_dot)->UseManualTime();
BENCHMARK(BM_glm_vec4_normalize)->UseManualTime();
BENCHMARK(BM_amal_vec4_normalize)->UseManualTime();

BENCHMARK(BM_glm_mat3_mat_add)->UseManualTime();
BENCHMARK(BM_amal_mat3_mat_add)->UseManualTime();
BENCHMARK(BM_glm_mat3_mat_mul_scalar)->UseManualTime();
BENCHMARK(BM_amal_mat3_mat_mul_scalar)->UseManualTime();
BENCHMARK(BM_glm_mat3_mat_mul_vec)->UseManualTime();
BENCHMARK(BM_amal_mat3_mat_mul_vec)->UseManualTime();
BENCHMARK(BM_glm_mat3_mat_mul_mat)->UseManualTime();
BENCHMARK(BM_amal_mat3_mat_mul_mat)->UseManualTime();
// BENCHMARK(BM_glm_mat3_mat_transpose)->UseManualTime();
// BENCHMARK(BM_amal_mat3_mat_transpose)->UseManualTime();

// BENCHMARK(BM_glm_mat4_mat_add)->UseManualTime();
// BENCHMARK(BM_amal_mat4_mat_add)->UseManualTime();
// BENCHMARK(BM_glm_mat4_mat_mul_scalar)->UseManualTime();
// BENCHMARK(BM_amal_mat4_mat_mul_scalar)->UseManualTime();
// BENCHMARK(BM_glm_mat4_mat_mul_vec)->UseManualTime();
// BENCHMARK(BM_amal_mat4_mat_mul_vec)->UseManualTime();
// BENCHMARK(BM_glm_mat4_mat_mul_mat)->UseManualTime();
// BENCHMARK(BM_amal_mat4_mat_mul_mat)->UseManualTime();
// BENCHMARK(BM_glm_mat4_mat_transpose)->UseManualTime();
// BENCHMARK(BM_amal_mat4_mat_transpose)->UseManualTime();
// BENCHMARK(BM_glm_mat4_translate)->UseManualTime();
// BENCHMARK(BM_amal_mat4_translate)->UseManualTime();

BENCHMARK_MAIN();
