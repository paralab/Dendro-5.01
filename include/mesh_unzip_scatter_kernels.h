// Same-level scatter kernels for Mesh::unzip_scatter.
// When an element and its target block share the same refinement level,
// scattering DG element values into the block is a deterministic integer
// reindex plus contiguous row copies — no floating-point coordinate math
// or rounding is required.

#ifndef DENDRO_MESH_UNZIP_SCATTER_KERNELS_H
#define DENDRO_MESH_UNZIP_SCATTER_KERNELS_H

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <type_traits>

#if defined(DENDRO_TENSOR_SIMD)
#include <immintrin.h>
#endif

// Opt-in build banner (-DDENDRO_REPORT_BUILD_CONFIG). Fires wherever this
// header compiles, including the downstream solver build — where these
// header-template flags actually take effect, not just at libdendro5 build.
#if defined(DENDRO_REPORT_BUILD_CONFIG)
#if defined(DENDRO_UNZIP_SCATTER_FAST)
#pragma message("[dendro] unzip integer-index fast path : ENABLED")
#else
#pragma message( \
    "[dendro] unzip integer-index fast path : disabled (build with -DDENDRO_UNZIP_SCATTER_FAST=ON)")
#endif
#if defined(DENDRO_TENSOR_SIMD)
#pragma message("[dendro] wavelet SIMD tensor kernels   : ENABLED")
#else
#pragma message( \
    "[dendro] wavelet SIMD tensor kernels   : disabled (build with -DDENDRO_TENSOR_SIMD=ON)")
#endif
#if defined(DENDRO_UNZIP_OMP)
#pragma message("[dendro] OpenMP block-parallel unzip   : ENABLED")
#else
#pragma message( \
    "[dendro] OpenMP block-parallel unzip   : disabled (build with -DDENDRO_UNZIP_OMP=ON)")
#endif
#endif  // DENDRO_REPORT_BUILD_CONFIG

namespace dendro {
namespace unzip {

// Compile-time-sized row copy used when an element is fully inside the
// block (the common case). gcc can sometimes do this from std::memcpy
// with a fixed size, but observed asm shows it emitting 7 scalar vmovsd
// per row for eO=6. Explicit intrinsics force the right pattern.
template <int N, typename T>
inline void scatter_copy_row(T* __restrict__ dst, const T* __restrict__ src) {
    std::memcpy(dst, src, (std::size_t)N * sizeof(T));
}

#if defined(DENDRO_TENSOR_SIMD) && defined(__AVX512F__)
// double, N=7 (eOrder=6): single masked 8-wide AVX-512 load+store.
template <>
inline void scatter_copy_row<7, double>(double* __restrict__ dst,
                                        const double* __restrict__ src) {
    const __mmask8 k = 0x7F;
    _mm512_mask_storeu_pd(dst, k, _mm512_maskz_loadu_pd(k, src));
}
// double, N=5 (eOrder=4): one 4-wide AVX2 load+store + 1 scalar tail.
template <>
inline void scatter_copy_row<5, double>(double* __restrict__ dst,
                                        const double* __restrict__ src) {
    _mm256_storeu_pd(dst, _mm256_loadu_pd(src));
    dst[4] = src[4];
}
// double, N=9 (eOrder=8): one 8-wide AVX-512 load+store + 1 scalar tail.
template <>
inline void scatter_copy_row<9, double>(double* __restrict__ dst,
                                        const double* __restrict__ src) {
    _mm512_storeu_pd(dst, _mm512_loadu_pd(src));
    dst[8] = src[8];
}
// double, N=3 (eOrder=2): 2 scalar + 1 scalar.
template <>
inline void scatter_copy_row<3, double>(double* __restrict__ dst,
                                        const double* __restrict__ src) {
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
}
#endif

// Inputs are precomputed once per (element, block) pair in unzip_scatter.
//   dgWVec       : pointer to DG element values, layout dof * (eOrder+1)^3
//   uzWVec       : pointer to unzipped output buffer (full dof * unSz vector
//   base) dof          : number of variables (typically 1 in BSSN per-call)
//   unSz, dgSz   : per-variable strides in uzWVec and dgWVec
//   offset       : block base offset within uzWVec
//   lx, ly, lz   : block allocation extents (with padding)
//   i0, j0, k0   : signed start indices of this element inside the block array
//                  (may be negative when the element overlaps a ghost padding
//                  region)
//   The clip-to-block bounds [i_lo,i_hi) etc. are computed inside the kernel.

template <typename T, unsigned int EORDER>
__attribute__((always_inline)) inline void scatter_same_level_specialized(
    const T* __restrict__ dgWVec, T* __restrict__ uzWVec, unsigned int dof,
    std::size_t unSz, std::size_t dgSz, std::size_t offset, unsigned int lx,
    unsigned int ly, unsigned int lz, int i0, int j0, int k0) {
    static_assert(std::is_trivially_copyable<T>::value,
                  "scatter kernel requires trivially copyable T");
    constexpr int eOp1   = (int)EORDER + 1;
    constexpr int eOp1Sq = eOp1 * eOp1;

    const int i_lo       = std::max(0, -i0);
    const int i_hi       = std::min(eOp1, (int)lx - i0);
    const int j_lo       = std::max(0, -j0);
    const int j_hi       = std::min(eOp1, (int)ly - j0);
    const int k_lo       = std::max(0, -k0);
    const int k_hi       = std::min(eOp1, (int)lz - k0);
    if (i_hi <= i_lo || j_hi <= j_lo || k_hi <= k_lo) return;

    const std::size_t lxy        = (std::size_t)lx * (std::size_t)ly;
    const bool full_row          = (i_lo == 0 && i_hi == eOp1);
    const std::size_t row_bytes_ = (std::size_t)(i_hi - i_lo) * sizeof(T);

    for (int k = k_lo; k < k_hi; ++k) {
        const std::size_t dst_zoff = (std::size_t)(k0 + k) * lxy;
        const std::size_t src_zoff = (std::size_t)k * (std::size_t)eOp1Sq;
        for (int j = j_lo; j < j_hi; ++j) {
            const std::size_t dst_off =
                offset + dst_zoff + (std::size_t)(j0 + j) * (std::size_t)lx +
                (std::size_t)(i0 + i_lo);
            const std::size_t src_off = src_zoff +
                                        (std::size_t)j * (std::size_t)eOp1 +
                                        (std::size_t)i_lo;
            if (full_row) {
                for (unsigned int v = 0; v < dof; ++v) {
                    scatter_copy_row<eOp1, T>(
                        uzWVec + (std::size_t)v * unSz + dst_off,
                        dgWVec + (std::size_t)v * dgSz + src_off);
                }
            } else {
                for (unsigned int v = 0; v < dof; ++v) {
                    std::memcpy(uzWVec + (std::size_t)v * unSz + dst_off,
                                dgWVec + (std::size_t)v * dgSz + src_off,
                                row_bytes_);
                }
            }
        }
    }
}

template <typename T>
__attribute__((always_inline)) inline void scatter_same_level_generic(
    const T* __restrict__ dgWVec, T* __restrict__ uzWVec, unsigned int eOrder,
    unsigned int dof, std::size_t unSz, std::size_t dgSz, std::size_t offset,
    unsigned int lx, unsigned int ly, unsigned int lz, int i0, int j0, int k0) {
    static_assert(std::is_trivially_copyable<T>::value,
                  "scatter kernel requires trivially copyable T");
    const int eOp1   = (int)eOrder + 1;
    const int eOp1Sq = eOp1 * eOp1;

    const int i_lo   = std::max(0, -i0);
    const int i_hi   = std::min(eOp1, (int)lx - i0);
    const int j_lo   = std::max(0, -j0);
    const int j_hi   = std::min(eOp1, (int)ly - j0);
    const int k_lo   = std::max(0, -k0);
    const int k_hi   = std::min(eOp1, (int)lz - k0);
    if (i_hi <= i_lo || j_hi <= j_lo || k_hi <= k_lo) return;

    const std::size_t row_bytes = (std::size_t)(i_hi - i_lo) * sizeof(T);
    const std::size_t lxy       = (std::size_t)lx * (std::size_t)ly;

    for (int k = k_lo; k < k_hi; ++k) {
        const std::size_t dst_zoff = (std::size_t)(k0 + k) * lxy;
        const std::size_t src_zoff = (std::size_t)k * (std::size_t)eOp1Sq;
        for (int j = j_lo; j < j_hi; ++j) {
            const std::size_t dst_off =
                offset + dst_zoff + (std::size_t)(j0 + j) * (std::size_t)lx +
                (std::size_t)(i0 + i_lo);
            const std::size_t src_off = src_zoff +
                                        (std::size_t)j * (std::size_t)eOp1 +
                                        (std::size_t)i_lo;
            for (unsigned int v = 0; v < dof; ++v) {
                std::memcpy(uzWVec + (std::size_t)v * unSz + dst_off,
                            dgWVec + (std::size_t)v * dgSz + src_off,
                            row_bytes);
            }
        }
    }
}

template <typename T>
__attribute__((always_inline)) inline void scatter_same_level_dispatch(
    const T* __restrict__ dgWVec, T* __restrict__ uzWVec, unsigned int eOrder,
    unsigned int dof, std::size_t unSz, std::size_t dgSz, std::size_t offset,
    unsigned int lx, unsigned int ly, unsigned int lz, int i0, int j0, int k0) {
    switch (eOrder) {
        case 2:
            scatter_same_level_specialized<T, 2>(dgWVec, uzWVec, dof, unSz,
                                                 dgSz, offset, lx, ly, lz, i0,
                                                 j0, k0);
            return;
        case 4:
            scatter_same_level_specialized<T, 4>(dgWVec, uzWVec, dof, unSz,
                                                 dgSz, offset, lx, ly, lz, i0,
                                                 j0, k0);
            return;
        case 6:
            scatter_same_level_specialized<T, 6>(dgWVec, uzWVec, dof, unSz,
                                                 dgSz, offset, lx, ly, lz, i0,
                                                 j0, k0);
            return;
        case 8:
            scatter_same_level_specialized<T, 8>(dgWVec, uzWVec, dof, unSz,
                                                 dgSz, offset, lx, ly, lz, i0,
                                                 j0, k0);
            return;
        default:
            scatter_same_level_generic<T>(dgWVec, uzWVec, eOrder, dof, unSz,
                                          dgSz, offset, lx, ly, lz, i0, j0, k0);
            return;
    }
}

// ---------------------------------------------------------------------------
// Fine -> Coarse (element at level bLev+1, block at bLev).
//
// Only every other element CG node aligns with a block CG node. With even
// eOrder (the BSSN case eO=6) every (i,j,k) with even index aligns; with odd
// eOrder, every (i,j,k) with odd index aligns. The mapping is
//   iix = ei*(eOrder/2) + PW + i/2   (for even eOrder)
// where ei = (eleX - blkX) / (S_blk/2) is the element's index within the
// block's fine-coordinate grid. PRECONDITION: eOrder must be even — the
// caller must check and fall back to the FP path otherwise (only relevant
// for odd-order solvers; BSSN eO=6 is even).
//
// half_eOrder == eOrder/2. The (i_h, j_h, k_h) loop visits half_eOrder+1
// points per dim — for eO=6, that's 4^3 = 64 writes per element-block pair.
// Reads from dgWVec are stride-2, so unlike the same-level case we cannot
// use memcpy — the writes are point-by-point. The win comes from removing
// the std::round / fabs / tolerance / FP coord math.

template <typename T, unsigned int EORDER>
__attribute__((always_inline)) inline void scatter_fine_to_coarse_specialized(
    const T* __restrict__ dgWVec, T* __restrict__ uzWVec, unsigned int dof,
    std::size_t unSz, std::size_t dgSz, std::size_t offset, unsigned int lx,
    unsigned int ly, unsigned int lz, int i0, int j0, int k0) {
    static_assert(EORDER % 2 == 0,
                  "fine_to_coarse_specialized requires even EORDER");
    static_assert(std::is_trivially_copyable<T>::value,
                  "scatter kernel requires trivially copyable T");
    constexpr int eOp1    = (int)EORDER + 1;
    constexpr int eOp1Sq  = eOp1 * eOp1;
    constexpr int half_eO = (int)EORDER / 2;
    constexpr int n_h     = half_eO + 1;  // points per dim
    const int i_hi        = std::min(n_h, (int)lx - i0);
    const int j_hi        = std::min(n_h, (int)ly - j0);
    const int k_hi        = std::min(n_h, (int)lz - k0);
    const int i_lo        = std::max(0, -i0);
    const int j_lo        = std::max(0, -j0);
    const int k_lo        = std::max(0, -k0);
    if (i_hi <= i_lo || j_hi <= j_lo || k_hi <= k_lo) return;
    const std::size_t lxy = (std::size_t)lx * (std::size_t)ly;
    for (int k_h = k_lo; k_h < k_hi; ++k_h) {
        const int k                = 2 * k_h;
        const std::size_t dst_zoff = (std::size_t)(k0 + k_h) * lxy;
        const std::size_t src_zoff = (std::size_t)k * (std::size_t)eOp1Sq;
        for (int j_h = j_lo; j_h < j_hi; ++j_h) {
            const int j = 2 * j_h;
            const std::size_t dst_yoff =
                (std::size_t)(j0 + j_h) * (std::size_t)lx;
            const std::size_t src_yoff = (std::size_t)j * (std::size_t)eOp1;
            for (unsigned int v = 0; v < dof; ++v) {
                T* outRow = uzWVec + (std::size_t)v * unSz + offset + dst_zoff +
                            dst_yoff;
                const T* inRow =
                    dgWVec + (std::size_t)v * dgSz + src_zoff + src_yoff;
                for (int i_h = i_lo; i_h < i_hi; ++i_h) {
                    outRow[(std::size_t)(i0 + i_h)] =
                        inRow[(std::size_t)(2 * i_h)];
                }
            }
        }
    }
}

template <typename T>
__attribute__((always_inline)) inline void scatter_fine_to_coarse_generic_even(
    const T* __restrict__ dgWVec, T* __restrict__ uzWVec, unsigned int eOrder,
    unsigned int dof, std::size_t unSz, std::size_t dgSz, std::size_t offset,
    unsigned int lx, unsigned int ly, unsigned int lz, int i0, int j0, int k0) {
    // generic runtime version, requires eOrder even
    const int eOp1    = (int)eOrder + 1;
    const int eOp1Sq  = eOp1 * eOp1;
    const int half_eO = (int)eOrder / 2;
    const int n_h     = half_eO + 1;
    const int i_hi    = std::min(n_h, (int)lx - i0);
    const int j_hi    = std::min(n_h, (int)ly - j0);
    const int k_hi    = std::min(n_h, (int)lz - k0);
    const int i_lo    = std::max(0, -i0);
    const int j_lo    = std::max(0, -j0);
    const int k_lo    = std::max(0, -k0);
    if (i_hi <= i_lo || j_hi <= j_lo || k_hi <= k_lo) return;
    const std::size_t lxy = (std::size_t)lx * (std::size_t)ly;
    for (int k_h = k_lo; k_h < k_hi; ++k_h) {
        const int k                = 2 * k_h;
        const std::size_t dst_zoff = (std::size_t)(k0 + k_h) * lxy;
        const std::size_t src_zoff = (std::size_t)k * (std::size_t)eOp1Sq;
        for (int j_h = j_lo; j_h < j_hi; ++j_h) {
            const int j = 2 * j_h;
            const std::size_t dst_yoff =
                (std::size_t)(j0 + j_h) * (std::size_t)lx;
            const std::size_t src_yoff = (std::size_t)j * (std::size_t)eOp1;
            for (unsigned int v = 0; v < dof; ++v) {
                T* outRow = uzWVec + (std::size_t)v * unSz + offset + dst_zoff +
                            dst_yoff;
                const T* inRow =
                    dgWVec + (std::size_t)v * dgSz + src_zoff + src_yoff;
                for (int i_h = i_lo; i_h < i_hi; ++i_h) {
                    outRow[(std::size_t)(i0 + i_h)] =
                        inRow[(std::size_t)(2 * i_h)];
                }
            }
        }
    }
}

// Returns true if the integer fast path handled the case; false means
// caller must use the FP fallback (odd eOrder).
template <typename T>
__attribute__((always_inline)) inline bool scatter_fine_to_coarse_dispatch(
    const T* __restrict__ dgWVec, T* __restrict__ uzWVec, unsigned int eOrder,
    unsigned int dof, std::size_t unSz, std::size_t dgSz, std::size_t offset,
    unsigned int lx, unsigned int ly, unsigned int lz, int i0, int j0, int k0) {
    if (eOrder % 2 != 0) return false;  // odd eOrder: half-integer offsets
    switch (eOrder) {
        case 2:
            scatter_fine_to_coarse_specialized<T, 2>(dgWVec, uzWVec, dof, unSz,
                                                     dgSz, offset, lx, ly, lz,
                                                     i0, j0, k0);
            return true;
        case 4:
            scatter_fine_to_coarse_specialized<T, 4>(dgWVec, uzWVec, dof, unSz,
                                                     dgSz, offset, lx, ly, lz,
                                                     i0, j0, k0);
            return true;
        case 6:
            scatter_fine_to_coarse_specialized<T, 6>(dgWVec, uzWVec, dof, unSz,
                                                     dgSz, offset, lx, ly, lz,
                                                     i0, j0, k0);
            return true;
        case 8:
            scatter_fine_to_coarse_specialized<T, 8>(dgWVec, uzWVec, dof, unSz,
                                                     dgSz, offset, lx, ly, lz,
                                                     i0, j0, k0);
            return true;
        default:
            scatter_fine_to_coarse_generic_even<T>(dgWVec, uzWVec, eOrder, dof,
                                                   unSz, dgSz, offset, lx, ly,
                                                   lz, i0, j0, k0);
            return true;
    }
}

}  // namespace unzip
}  // namespace dendro

#endif  // DENDRO_MESH_UNZIP_SCATTER_KERNELS_H
