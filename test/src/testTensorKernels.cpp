/**
 * @file testTensorKernels.cpp
 * @brief Equality and lane-coverage tests for the FEM tensor-product kernels.
 *
 * The SIMD paths in FEM/src/tensor.cpp are dispatched only for M in {5,7,9};
 * every other M runs the scalar fallback. That gives a free oracle check: at a
 * non-dispatched M the library IS the scalar path, so agreement there pins the
 * reference below to the real semantics before it is used to judge M=5,7,9.
 *
 * Inputs are multiples of 1/64 with bounded magnitude, so every product and sum
 * is exact in double and the FMA-based SIMD path must agree with scalar mul+add
 * bit for bit. Outputs are poisoned with NaN so a lane that is never written is
 * distinguishable from a lane that is written wrongly.
 */

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <mpi.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "tensor.h"

namespace {

// M values FEM/src/tensor.cpp dispatches to a SIMD kernel.
const std::vector<int> kSimdM     = {5, 7, 9};
// M values with no SIMD specialization: the library runs its scalar fallback.
const std::vector<int> kScalarM   = {3, 4, 6, 8, 10};

struct Rng {
    std::uint64_t s;
    explicit Rng(std::uint64_t seed) : s(seed) {}
    // multiples of 1/64 in [-4, 4]: products and 11-term sums stay exact
    double next() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return (double)((int)((s >> 33) % 513) - 256) / 64.0;
    }
};

// The five references mirror the accumulation order of the scalar fallbacks in
// FEM/src/tensor.cpp, so bitwise comparison is meaningful.

void refAIIX(int M, const double* A, const double* X, double* Y) {
    const int MM = M * M;
    for (int i = 0; i < M; ++i) {
        double d = A[i];
        for (int j = 0; j < MM; ++j) Y[i * MM + j] = d * X[j];
        for (int k = 1; k < M; ++k) {
            d = A[i + k * M];
            for (int j = 0; j < MM; ++j) Y[i * MM + j] += d * X[MM * k + j];
        }
    }
}

void refIIAX(int M, const double* A, const double* X, double* Y) {
    for (int i = 0; i < M * M; ++i)
        for (int j = 0; j < M; ++j) {
            double e = 0;
            for (int k = 0; k < M; ++k) e += X[i * M + k] * A[k * M + j];
            Y[i * M + j] = e;
        }
}

void refIAIX(int M, const double* A, const double* X, double* Y) {
    const int MM = M * M;
    for (int ib = 0; ib < M; ++ib)
        for (int i = 0; i < M; ++i) {
            double d = A[i];
            for (int j = 0; j < M; ++j)
                Y[ib * MM + i * M + j] = d * X[ib * MM + j];
            for (int k = 1; k < M; ++k) {
                d = A[i + k * M];
                for (int j = 0; j < M; ++j)
                    Y[ib * MM + i * M + j] += d * X[ib * MM + M * k + j];
            }
        }
}

void refIAX2D(int M, const double* A, const double* X, double* Y) {
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) {
            double e = 0;
            for (int k = 0; k < M; ++k) e += X[i * M + k] * A[k * M + j];
            Y[i * M + j] = e;
        }
}

void refAIX2D(int M, const double* A, const double* X, double* Y) {
    for (int i = 0; i < M; ++i) {
        double d = A[i];
        for (int j = 0; j < M; ++j) Y[i * M + j] = d * X[j];
        for (int k = 1; k < M; ++k) {
            d = A[i + k * M];
            for (int j = 0; j < M; ++j) Y[i * M + j] += d * X[M * k + j];
        }
    }
}

typedef void (*Kernel)(const int, const double* __restrict__,
                       const double* __restrict__, double* __restrict__);
typedef void (*Ref)(int, const double*, const double*, double*);

struct KernelUnderTest {
    const char* name;
    Kernel lib;
    Ref ref;
    bool cubic;  // true: X and Y are M^3; false: M^2 (2D face kernels)
};

const std::vector<KernelUnderTest>& kernels() {
    static const std::vector<KernelUnderTest> k = {
        {"AIIX", DENDRO_TENSOR_AIIX_APPLY_ELEM, refAIIX, true},
        {"IIAX", DENDRO_TENSOR_IIAX_APPLY_ELEM, refIIAX, true},
        {"IAIX", DENDRO_TENSOR_IAIX_APPLY_ELEM, refIAIX, true},
        {"IAX_2D", DENDRO_TENSOR_IAX_APPLY_ELEM_2D, refIAX2D, false},
        {"AIX_2D", DENDRO_TENSOR_AIX_APPLY_ELEM_2D, refAIX2D, false},
    };
    return k;
}

struct Result {
    int unwritten = 0;  // lanes the kernel never touched (still NaN)
    int wrong     = 0;  // lanes written with a value the reference disagrees on
};

Result compare(const KernelUnderTest& k, int M, std::uint64_t seed) {
    const std::size_t MM  = (std::size_t)M * M;
    const std::size_t n   = k.cubic ? MM * (std::size_t)M : MM;
    Rng rng(seed);
    std::vector<double> A(MM), X(n);
    for (double& v : A) v = rng.next();
    for (double& v : X) v = rng.next();

    const double nan = std::nan("");
    std::vector<double> got(n, nan), want(n, nan);
    k.lib(M, A.data(), X.data(), got.data());
    k.ref(M, A.data(), X.data(), want.data());

    Result r;
    for (std::size_t i = 0; i < n; ++i) {
        if (std::isnan(got[i]) && !std::isnan(want[i]))
            r.unwritten++;
        else if (got[i] != want[i])
            r.wrong++;
    }
    return r;
}

}  // namespace

TEST_CASE("scalar reference matches the library where SIMD does not dispatch") {
    for (int M : kScalarM)
        for (const KernelUnderTest& k : kernels())
            for (std::uint64_t seed : {1u, 2u, 3u}) {
                const Result r            = compare(k, M, seed);
                const std::string kernel  = k.name;
                CAPTURE(kernel);
                CAPTURE(M);
                CAPTURE(seed);
                CHECK(r.unwritten == 0);
                CHECK(r.wrong == 0);
            }
}

TEST_CASE("tensor kernels write every output lane at M = 5, 7, 9") {
    for (int M : kSimdM)
        for (const KernelUnderTest& k : kernels())
            for (std::uint64_t seed : {1u, 2u, 3u}) {
                const Result r            = compare(k, M, seed);
                const std::string kernel  = k.name;
                CAPTURE(kernel);
                CAPTURE(M);
                CAPTURE(seed);
                CHECK(r.unwritten == 0);
            }
}

TEST_CASE("tensor kernels equal the scalar reference at M = 5, 7, 9") {
    for (int M : kSimdM)
        for (const KernelUnderTest& k : kernels())
            for (std::uint64_t seed : {1u, 2u, 3u}) {
                const Result r            = compare(k, M, seed);
                const std::string kernel  = k.name;
                CAPTURE(kernel);
                CAPTURE(M);
                CAPTURE(seed);
                CHECK(r.wrong == 0);
            }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    doctest::Context ctx;
    ctx.applyCommandLine(argc, argv);
    const int res = ctx.run();
    MPI_Finalize();
    return res;
}
