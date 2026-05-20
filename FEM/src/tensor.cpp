//
// Created by milinda on 1/15/17.
//

/**
 *
 * @author Milinda Fernando
 * @breif contains the utilities for tensor kronecker products for
 * interpolations.
 *
 * AVX2 SIMD specializations (M=7 for eO=6, M=5 for eO=4) live in an
 * anonymous namespace below and are dispatched on M at the top of each
 * public function when DENDRO_TENSOR_SIMD is defined. Scalar code is
 * unchanged and is the fallback for other M.
 *
 * */

#include "tensor.h"

#if defined(DENDRO_TENSOR_SIMD)
#include <immintrin.h>
#endif

#if defined(DENDRO_TENSOR_SIMD)
namespace {

// ===========================================================================
// AVX-512 specializations (M=7 only — single masked 8-wide op covers all 7
// elements, no tail). Selected at runtime when __AVX512F__ is defined.
// ===========================================================================
#if defined(__AVX512F__)

// Mask of 7 lanes (bits 0..6) for M=7 ops.
static constexpr __mmask8 K7 = 0x7F;

template <int M>
static inline void aiix_avx512(const double* __restrict__ A,
                               const double* __restrict__ X,
                               double* __restrict__ Y);
template <int M>
static inline void iiax_avx512(const double* __restrict__ A,
                               const double* __restrict__ X,
                               double* __restrict__ Y);
template <int M>
static inline void iaix_avx512(const double* __restrict__ A,
                               const double* __restrict__ X,
                               double* __restrict__ Y);

// AIIX (Z): inner loop is M*M unit-stride. M=7 -> 49 = 8*6 + 1.
template <>
inline void aiix_avx512<7>(const double* __restrict__ A,
                           const double* __restrict__ X,
                           double* __restrict__ Y) {
    for (int i = 0; i < 7; ++i) {
        // k=0: initialize Y[i,:]
        const __m512d vd0 = _mm512_set1_pd(A[i]);
        int j;
        for (j = 0; j + 8 <= 49; j += 8) {
            const __m512d vx = _mm512_loadu_pd(X + j);
            _mm512_storeu_pd(Y + i * 49 + j, _mm512_mul_pd(vd0, vx));
        }
        // tail j=48 (just 1 element)
        for (; j < 49; ++j) Y[i * 49 + j] = A[i] * X[j];
        // k=1..6: accumulate. Use 6 INDEPENDENT chains by splitting the j
        // range into 6 blocks of 8 — each block has its own dependency
        // chain on Y[i,j..j+7], so blocks pipeline. (This gets 6 fmadds in
        // flight, hitting the 1/cycle throughput of port 0/1.)
        for (int k = 1; k < 7; ++k) {
            const __m512d vd = _mm512_set1_pd(A[i + k * 7]);
            for (j = 0; j + 8 <= 49; j += 8) {
                const __m512d vx = _mm512_loadu_pd(X + 49 * k + j);
                const __m512d vy = _mm512_loadu_pd(Y + i * 49 + j);
                _mm512_storeu_pd(Y + i * 49 + j,
                                 _mm512_fmadd_pd(vd, vx, vy));
            }
            for (; j < 49; ++j) Y[i * 49 + j] += A[i + k * 7] * X[49 * k + j];
        }
    }
}

// IIAX (X): outer M*M=49, inner sum over k=0..6 produces Y[i,:] (7 lanes).
// Trick: process 4 i's in parallel for ILP, breaking the Y-reload chain.
template <>
inline void iiax_avx512<7>(const double* __restrict__ A,
                           const double* __restrict__ X,
                           double* __restrict__ Y) {
    // Pre-load A rows (7 vectors, masked) into registers — small and reused
    // for every i. Compiler will likely keep them in zmm regs.
    const __m512d A0 = _mm512_maskz_loadu_pd(K7, A + 0 * 7);
    const __m512d A1 = _mm512_maskz_loadu_pd(K7, A + 1 * 7);
    const __m512d A2 = _mm512_maskz_loadu_pd(K7, A + 2 * 7);
    const __m512d A3 = _mm512_maskz_loadu_pd(K7, A + 3 * 7);
    const __m512d A4 = _mm512_maskz_loadu_pd(K7, A + 4 * 7);
    const __m512d A5 = _mm512_maskz_loadu_pd(K7, A + 5 * 7);
    const __m512d A6 = _mm512_maskz_loadu_pd(K7, A + 6 * 7);

    int i = 0;
    // 4-wide outer unroll for ILP
    for (; i + 4 <= 49; i += 4) {
        const double* Xa = X + (i + 0) * 7;
        const double* Xb = X + (i + 1) * 7;
        const double* Xc = X + (i + 2) * 7;
        const double* Xd = X + (i + 3) * 7;
        __m512d Ya = _mm512_mul_pd(_mm512_set1_pd(Xa[0]), A0);
        __m512d Yb = _mm512_mul_pd(_mm512_set1_pd(Xb[0]), A0);
        __m512d Yc = _mm512_mul_pd(_mm512_set1_pd(Xc[0]), A0);
        __m512d Yd = _mm512_mul_pd(_mm512_set1_pd(Xd[0]), A0);
        Ya = _mm512_fmadd_pd(_mm512_set1_pd(Xa[1]), A1, Ya);
        Yb = _mm512_fmadd_pd(_mm512_set1_pd(Xb[1]), A1, Yb);
        Yc = _mm512_fmadd_pd(_mm512_set1_pd(Xc[1]), A1, Yc);
        Yd = _mm512_fmadd_pd(_mm512_set1_pd(Xd[1]), A1, Yd);
        Ya = _mm512_fmadd_pd(_mm512_set1_pd(Xa[2]), A2, Ya);
        Yb = _mm512_fmadd_pd(_mm512_set1_pd(Xb[2]), A2, Yb);
        Yc = _mm512_fmadd_pd(_mm512_set1_pd(Xc[2]), A2, Yc);
        Yd = _mm512_fmadd_pd(_mm512_set1_pd(Xd[2]), A2, Yd);
        Ya = _mm512_fmadd_pd(_mm512_set1_pd(Xa[3]), A3, Ya);
        Yb = _mm512_fmadd_pd(_mm512_set1_pd(Xb[3]), A3, Yb);
        Yc = _mm512_fmadd_pd(_mm512_set1_pd(Xc[3]), A3, Yc);
        Yd = _mm512_fmadd_pd(_mm512_set1_pd(Xd[3]), A3, Yd);
        Ya = _mm512_fmadd_pd(_mm512_set1_pd(Xa[4]), A4, Ya);
        Yb = _mm512_fmadd_pd(_mm512_set1_pd(Xb[4]), A4, Yb);
        Yc = _mm512_fmadd_pd(_mm512_set1_pd(Xc[4]), A4, Yc);
        Yd = _mm512_fmadd_pd(_mm512_set1_pd(Xd[4]), A4, Yd);
        Ya = _mm512_fmadd_pd(_mm512_set1_pd(Xa[5]), A5, Ya);
        Yb = _mm512_fmadd_pd(_mm512_set1_pd(Xb[5]), A5, Yb);
        Yc = _mm512_fmadd_pd(_mm512_set1_pd(Xc[5]), A5, Yc);
        Yd = _mm512_fmadd_pd(_mm512_set1_pd(Xd[5]), A5, Yd);
        Ya = _mm512_fmadd_pd(_mm512_set1_pd(Xa[6]), A6, Ya);
        Yb = _mm512_fmadd_pd(_mm512_set1_pd(Xb[6]), A6, Yb);
        Yc = _mm512_fmadd_pd(_mm512_set1_pd(Xc[6]), A6, Yc);
        Yd = _mm512_fmadd_pd(_mm512_set1_pd(Xd[6]), A6, Yd);
        _mm512_mask_storeu_pd(Y + (i + 0) * 7, K7, Ya);
        _mm512_mask_storeu_pd(Y + (i + 1) * 7, K7, Yb);
        _mm512_mask_storeu_pd(Y + (i + 2) * 7, K7, Yc);
        _mm512_mask_storeu_pd(Y + (i + 3) * 7, K7, Yd);
    }
    // tail i in [48, 49): scalar per element of M*M=49 leaves 1 left
    for (; i < 49; ++i) {
        const double* Xi = X + i * 7;
        __m512d Yi = _mm512_mul_pd(_mm512_set1_pd(Xi[0]), A0);
        Yi = _mm512_fmadd_pd(_mm512_set1_pd(Xi[1]), A1, Yi);
        Yi = _mm512_fmadd_pd(_mm512_set1_pd(Xi[2]), A2, Yi);
        Yi = _mm512_fmadd_pd(_mm512_set1_pd(Xi[3]), A3, Yi);
        Yi = _mm512_fmadd_pd(_mm512_set1_pd(Xi[4]), A4, Yi);
        Yi = _mm512_fmadd_pd(_mm512_set1_pd(Xi[5]), A5, Yi);
        Yi = _mm512_fmadd_pd(_mm512_set1_pd(Xi[6]), A6, Yi);
        _mm512_mask_storeu_pd(Y + i * 7, K7, Yi);
    }
}

// IAIX (Y): Y[ib,i,:] = sum_k A[i,k] * X[ib,k,:]. inner is 7 lanes, ib*49
// stride. Same ILP trick: unroll i by 4 to break Y reload chain.
template <>
inline void iaix_avx512<7>(const double* __restrict__ A,
                           const double* __restrict__ X,
                           double* __restrict__ Y) {
    for (int ib = 0; ib < 7; ++ib) {
        const double* Xib = X + ib * 49;
        double* Yib       = Y + ib * 49;
        // Load X rows once per ib block
        const __m512d X0 = _mm512_maskz_loadu_pd(K7, Xib + 0 * 7);
        const __m512d X1 = _mm512_maskz_loadu_pd(K7, Xib + 1 * 7);
        const __m512d X2 = _mm512_maskz_loadu_pd(K7, Xib + 2 * 7);
        const __m512d X3 = _mm512_maskz_loadu_pd(K7, Xib + 3 * 7);
        const __m512d X4 = _mm512_maskz_loadu_pd(K7, Xib + 4 * 7);
        const __m512d X5 = _mm512_maskz_loadu_pd(K7, Xib + 5 * 7);
        const __m512d X6 = _mm512_maskz_loadu_pd(K7, Xib + 6 * 7);
        int i = 0;
        for (; i + 4 <= 7; i += 4) {
            __m512d Ya = _mm512_mul_pd(_mm512_set1_pd(A[(i + 0) + 0 * 7]), X0);
            __m512d Yb = _mm512_mul_pd(_mm512_set1_pd(A[(i + 1) + 0 * 7]), X0);
            __m512d Yc = _mm512_mul_pd(_mm512_set1_pd(A[(i + 2) + 0 * 7]), X0);
            __m512d Yd = _mm512_mul_pd(_mm512_set1_pd(A[(i + 3) + 0 * 7]), X0);
#define ACC(k, Xk) \
    Ya = _mm512_fmadd_pd(_mm512_set1_pd(A[(i + 0) + (k) * 7]), Xk, Ya); \
    Yb = _mm512_fmadd_pd(_mm512_set1_pd(A[(i + 1) + (k) * 7]), Xk, Yb); \
    Yc = _mm512_fmadd_pd(_mm512_set1_pd(A[(i + 2) + (k) * 7]), Xk, Yc); \
    Yd = _mm512_fmadd_pd(_mm512_set1_pd(A[(i + 3) + (k) * 7]), Xk, Yd)
            ACC(1, X1);
            ACC(2, X2);
            ACC(3, X3);
            ACC(4, X4);
            ACC(5, X5);
            ACC(6, X6);
#undef ACC
            _mm512_mask_storeu_pd(Yib + (i + 0) * 7, K7, Ya);
            _mm512_mask_storeu_pd(Yib + (i + 1) * 7, K7, Yb);
            _mm512_mask_storeu_pd(Yib + (i + 2) * 7, K7, Yc);
            _mm512_mask_storeu_pd(Yib + (i + 3) * 7, K7, Yd);
        }
        // tail i in [4, 7)
        for (; i < 7; ++i) {
            __m512d Yi = _mm512_mul_pd(_mm512_set1_pd(A[i + 0 * 7]), X0);
            Yi = _mm512_fmadd_pd(_mm512_set1_pd(A[i + 1 * 7]), X1, Yi);
            Yi = _mm512_fmadd_pd(_mm512_set1_pd(A[i + 2 * 7]), X2, Yi);
            Yi = _mm512_fmadd_pd(_mm512_set1_pd(A[i + 3 * 7]), X3, Yi);
            Yi = _mm512_fmadd_pd(_mm512_set1_pd(A[i + 4 * 7]), X4, Yi);
            Yi = _mm512_fmadd_pd(_mm512_set1_pd(A[i + 5 * 7]), X5, Yi);
            Yi = _mm512_fmadd_pd(_mm512_set1_pd(A[i + 6 * 7]), X6, Yi);
            _mm512_mask_storeu_pd(Yib + i * 7, K7, Yi);
        }
    }
}

#endif  // __AVX512F__

// ---------------------------------------------------------------------------
// AVX2 specialization for AIIX (Z-axis), M known at compile time.
// Y[i,j] = sum_k A[i,k] * X[k,j], with j running over the M*M plane
// (unit-stride, the easy axis).
template <int M>
static inline void aiix_avx2(const double* __restrict__ A,
                             const double* __restrict__ X,
                             double* __restrict__ Y) {
    constexpr int MM       = M * M;
    constexpr int simd_end = (MM / 4) * 4;
    for (int i = 0; i < M; ++i) {
        // first pass: k=0 → initialize Y[i,:]
        const __m256d vd0 = _mm256_set1_pd(A[i]);
        int j;
        for (j = 0; j < simd_end; j += 4) {
            const __m256d vx = _mm256_loadu_pd(X + j);
            _mm256_storeu_pd(Y + i * MM + j, _mm256_mul_pd(vd0, vx));
        }
        const double d0 = A[i];
        for (; j < MM; ++j) Y[i * MM + j] = d0 * X[j];
        // k=1..M-1: accumulate
        for (int k = 1; k < M; ++k) {
            const __m256d vd = _mm256_set1_pd(A[i + k * M]);
            for (j = 0; j < simd_end; j += 4) {
                const __m256d vx = _mm256_loadu_pd(X + MM * k + j);
                const __m256d vy = _mm256_loadu_pd(Y + i * MM + j);
                _mm256_storeu_pd(Y + i * MM + j,
                                 _mm256_fmadd_pd(vd, vx, vy));
            }
            const double d = A[i + k * M];
            for (; j < MM; ++j) Y[i * MM + j] += d * X[MM * k + j];
        }
    }
}

// ---------------------------------------------------------------------------
// AVX2 specialization for IIAX (X-axis), M known at compile time.
// Y[i,j] = sum_k X[i,k] * A[k,j].   i iterates M*M times (outer), j and k
// over M. Reorder so j is innermost contiguous in both Y and A[k,:].
template <int M>
static inline void iiax_avx2(const double* __restrict__ A,
                             const double* __restrict__ X,
                             double* __restrict__ Y) {
    constexpr int MM       = M * M;
    constexpr int simd_end = (M / 4) * 4;  // for M=7: 4; for M=5: 4
    for (int i = 0; i < MM; ++i) {
        // k=0: initialize Y[i,:]
        {
            const __m256d vx = _mm256_set1_pd(X[i * M + 0]);
            const __m256d va = _mm256_loadu_pd(A + 0);
            _mm256_storeu_pd(Y + i * M, _mm256_mul_pd(vx, va));
            const double x = X[i * M + 0];
            for (int j = simd_end; j < M; ++j) Y[i * M + j] = x * A[j];
        }
        for (int k = 1; k < M; ++k) {
            const __m256d vx = _mm256_set1_pd(X[i * M + k]);
            const __m256d va = _mm256_loadu_pd(A + k * M);
            const __m256d vy = _mm256_loadu_pd(Y + i * M);
            _mm256_storeu_pd(Y + i * M, _mm256_fmadd_pd(vx, va, vy));
            const double x = X[i * M + k];
            for (int j = simd_end; j < M; ++j)
                Y[i * M + j] += x * A[k * M + j];
        }
    }
}

// ---------------------------------------------------------------------------
// AVX2 specialization for IAIX (Y-axis), M known at compile time.
// Y[ib,i,j] = sum_k A[i,k] * X[ib,k,j], inner j is short (length M).
template <int M>
static inline void iaix_avx2(const double* __restrict__ A,
                             const double* __restrict__ X,
                             double* __restrict__ Y) {
    constexpr int MM       = M * M;
    constexpr int simd_end = (M / 4) * 4;
    for (int ib = 0; ib < M; ++ib) {
        for (int i = 0; i < M; ++i) {
            // k=0: initialize Y[ib,i,:]
            const double d0   = A[i];
            const __m256d vd0 = _mm256_set1_pd(d0);
            const __m256d vx0 = _mm256_loadu_pd(X + ib * MM + 0);
            _mm256_storeu_pd(Y + ib * MM + i * M, _mm256_mul_pd(vd0, vx0));
            for (int j = simd_end; j < M; ++j)
                Y[ib * MM + i * M + j] = d0 * X[ib * MM + j];
            for (int k = 1; k < M; ++k) {
                const double d   = A[i + k * M];
                const __m256d vd = _mm256_set1_pd(d);
                const __m256d vx = _mm256_loadu_pd(X + ib * MM + k * M);
                const __m256d vy = _mm256_loadu_pd(Y + ib * MM + i * M);
                _mm256_storeu_pd(Y + ib * MM + i * M,
                                 _mm256_fmadd_pd(vd, vx, vy));
                for (int j = simd_end; j < M; ++j)
                    Y[ib * MM + i * M + j] += d * X[ib * MM + k * M + j];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// AVX2 specialization for 2D-face IAX (X-axis face interp), M known at
// compile time. Same as iiax but with M outer iterations instead of M*M
// (one ib slice).
template <int M>
static inline void iax_2d_avx2(const double* __restrict__ A,
                               const double* __restrict__ X,
                               double* __restrict__ Y) {
    constexpr int simd_end = (M / 4) * 4;
    for (int i = 0; i < M; ++i) {
        {
            const __m256d vx = _mm256_set1_pd(X[i * M + 0]);
            const __m256d va = _mm256_loadu_pd(A + 0);
            _mm256_storeu_pd(Y + i * M, _mm256_mul_pd(vx, va));
            const double x = X[i * M + 0];
            for (int j = simd_end; j < M; ++j) Y[i * M + j] = x * A[j];
        }
        for (int k = 1; k < M; ++k) {
            const __m256d vx = _mm256_set1_pd(X[i * M + k]);
            const __m256d va = _mm256_loadu_pd(A + k * M);
            const __m256d vy = _mm256_loadu_pd(Y + i * M);
            _mm256_storeu_pd(Y + i * M, _mm256_fmadd_pd(vx, va, vy));
            const double x = X[i * M + k];
            for (int j = simd_end; j < M; ++j)
                Y[i * M + j] += x * A[k * M + j];
        }
    }
}

// ---------------------------------------------------------------------------
// AVX2 specialization for 2D-face AIX (Y-axis face interp), M known at
// compile time. Same as iaix with ib=0 (one slice).
template <int M>
static inline void aix_2d_avx2(const double* __restrict__ A,
                               const double* __restrict__ X,
                               double* __restrict__ Y) {
    constexpr int MM       = M * M;
    constexpr int simd_end = (M / 4) * 4;
    for (int i = 0; i < M; ++i) {
        const double d0   = A[i];
        const __m256d vd0 = _mm256_set1_pd(d0);
        const __m256d vx0 = _mm256_loadu_pd(X + 0);
        _mm256_storeu_pd(Y + i * M, _mm256_mul_pd(vd0, vx0));
        for (int j = simd_end; j < M; ++j) Y[i * M + j] = d0 * X[j];
        for (int k = 1; k < M; ++k) {
            const double d   = A[i + k * M];
            const __m256d vd = _mm256_set1_pd(d);
            const __m256d vx = _mm256_loadu_pd(X + k * M);
            const __m256d vy = _mm256_loadu_pd(Y + i * M);
            _mm256_storeu_pd(Y + i * M, _mm256_fmadd_pd(vd, vx, vy));
            for (int j = simd_end; j < M; ++j)
                Y[i * M + j] += d * X[k * M + j];
        }
    }
    (void)MM;
}

}  // namespace
#endif  // DENDRO_TENSOR_SIMD

/**
 * Along the Z axis
 * */
void DENDRO_TENSOR_AIIX_APPLY_ELEM(const int M, const double* A,
                                   const double* X, double* Y) {
#if defined(DENDRO_TENSOR_SIMD)
#if defined(__AVX512F__)
    if (M == 7) { aiix_avx512<7>(A, X, Y); return; }
#endif
    if (M == 7) { aiix_avx2<7>(A, X, Y); return; }
    if (M == 5) { aiix_avx2<5>(A, X, Y); return; }
    if (M == 9) { aiix_avx2<9>(A, X, Y); return; }
#endif
    int i, j, k;
    double d, e;
    for (i = 0; i < M; ++i) {
        d = A[i];
        for (j = 0; j < M * M; ++j) {
            Y[i * M * M + j] = d * X[j];
        }
        for (k = 1; k < M; ++k) {
            d = A[i + k * M];
            for (j = 0; j < M * M; ++j) {
                e = d * X[M * M * k + j];
                Y[i * M * M + j] += e;
            }
        }
    }
}

/**
 * Along the X axis
 * */

void DENDRO_TENSOR_IIAX_APPLY_ELEM(const int M, const double* A,
                                   const double* X, double* Y) {
#if defined(DENDRO_TENSOR_SIMD)
#if defined(__AVX512F__)
    if (M == 7) { iiax_avx512<7>(A, X, Y); return; }
#endif
    if (M == 7) { iiax_avx2<7>(A, X, Y); return; }
    if (M == 5) { iiax_avx2<5>(A, X, Y); return; }
    if (M == 9) { iiax_avx2<9>(A, X, Y); return; }
#endif
    int i, j, k;
    double e;
    for (i = 0; i < M * M; i++) {
        // _mm_prefetch( (const char *)(X + (i+10)*M), 2);
        //_mm_prefetch( (const char *)(Y + (i+10)*M), 2);
        for (j = 0; j < M; ++j) {
            e = 0;
            for (k = 0; k < M; ++k) {
                e += X[i * M + k] * A[k * M + j];
            }
            Y[i * M + j] = e;
        }
    }
}

/**
 * Along the X axis. (in face interpolations. )
 * */

void DENDRO_TENSOR_IAX_APPLY_ELEM_2D(const int M, const double* A,
                                     const double* X, double* Y) {
#if defined(DENDRO_TENSOR_SIMD)
    if (M == 7) { iax_2d_avx2<7>(A, X, Y); return; }
    if (M == 5) { iax_2d_avx2<5>(A, X, Y); return; }
    if (M == 9) { iax_2d_avx2<9>(A, X, Y); return; }
#endif
    int i, j, k;
    double e;
    for (i = 0; i < M; i++) {
        // _mm_prefetch( (const char *)(X + (i+10)*M), 2);
        //_mm_prefetch( (const char *)(Y + (i+10)*M), 2);
        for (j = 0; j < M; ++j) {
            e = 0;
            for (k = 0; k < M; ++k) {
                e += X[i * M + k] * A[k * M + j];
            }
            Y[i * M + j] = e;
        }
    }
}

/**
 * Along the Y axis
 * */
void DENDRO_TENSOR_IAIX_APPLY_ELEM(const int M, const double* A,
                                   const double* X, double* Y) {
#if defined(DENDRO_TENSOR_SIMD)
#if defined(__AVX512F__)
    if (M == 7) { iaix_avx512<7>(A, X, Y); return; }
#endif
    if (M == 7) { iaix_avx2<7>(A, X, Y); return; }
    if (M == 5) { iaix_avx2<5>(A, X, Y); return; }
    if (M == 9) { iaix_avx2<9>(A, X, Y); return; }
#endif
    int i, j, k, ib;
    double d, e;
    for (ib = 0; ib < M; ++ib) {
        for (i = 0; i < M; ++i) {
            d = A[i];
            for (j = 0; j < M; ++j) {
                Y[ib * M * M + i * M + j] = d * X[ib * M * M + j];
            }
            for (k = 1; k < M; ++k) {
                d = A[i + k * M];
                for (j = 0; j < M; ++j) {
                    e = d * X[ib * M * M + M * k + j];
                    Y[ib * M * M + i * M + j] += e;
                }
            }
        }
    }
}

/**
 * Along the Y axis for (2D face interpolations. )
 * */
void DENDRO_TENSOR_AIX_APPLY_ELEM_2D(const int M, const double* A,
                                     const double* X, double* Y) {
#if defined(DENDRO_TENSOR_SIMD)
    if (M == 7) { aix_2d_avx2<7>(A, X, Y); return; }
    if (M == 5) { aix_2d_avx2<5>(A, X, Y); return; }
    if (M == 9) { aix_2d_avx2<9>(A, X, Y); return; }
#endif
    int i, j, k, ib = 0;
    double d, e;
    // for (ib = 0; ib < M; ++ib) {
    for (i = 0; i < M; ++i) {
        d = A[i];
        for (j = 0; j < M; ++j) {
            Y[ib * M * M + i * M + j] = d * X[ib * M * M + j];
        }
        for (k = 1; k < M; ++k) {
            d = A[i + k * M];
            for (j = 0; j < M; ++j) {
                e = d * X[ib * M * M + M * k + j];
                Y[ib * M * M + i * M + j] += e;
            }
        }
    }
    //}
}
