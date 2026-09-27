// Checks the integer-index unzip scatter kernels against the FP fallback over
// every element/block placement, comparing the written set as well as the
// values. No mesh, no MPI, no initial data.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "mesh_unzip_scatter_kernels.h"

namespace {

constexpr double D_COMPAR_TOL = 1e-10;  // as in mesh.tcc
constexpr double SENTINEL     = -1e30;  // marks "never written"

struct Geom {
    unsigned int maxDepth;  // m_uiMaxDepth
    unsigned int bLev;      // block regular grid level
    unsigned int nEle;      // block elements per dim (1 << (regLev - blkLev))
    unsigned int eOrder;
    unsigned int PW;   // ot::Block ctor: eOrder >> 1
    unsigned int lx;   // getAllocationSzX()
    uint64_t blkX;     // blkNode origin (aligned)
};

Geom makeGeom(unsigned int eOrder, unsigned int nEle) {
    Geom g;
    g.maxDepth = 13;
    g.bLev     = 11;
    g.nEle     = nEle;
    g.eOrder   = eOrder;
    g.PW       = eOrder >> 1u;
    g.lx       = eOrder * nEle + 1 + 2 * g.PW;
    g.blkX     = 4096;
    return g;
}

// FP fallback, same-level branch (element level == block level).
void fpSameLevel(const Geom& g, int ei, int ej, int ek, const double* dg,
                 double* uz, unsigned int dof, std::size_t unSz,
                 std::size_t dgSz) {
    const unsigned int eOrder = g.eOrder;
    const unsigned int lx = g.lx, ly = g.lx, lz = g.lx;
    const uint64_t szEle   = (uint64_t)1u << (g.maxDepth - g.bLev);
    const uint64_t blkSz   = szEle * g.nEle;
    const double hx        = szEle / (double)eOrder;
    const double xmin = g.blkX - g.PW * hx, xmax = g.blkX + blkSz + g.PW * hx;
    const double ymin = xmin, ymax = xmax, zmin = xmin, zmax = xmax;
    const double eX = (double)(g.blkX + (int64_t)ei * (int64_t)szEle);
    const double eY = (double)(g.blkX + (int64_t)ej * (int64_t)szEle);
    const double eZ = (double)(g.blkX + (int64_t)ek * (int64_t)szEle);
    const double hh = szEle / (double)eOrder, invhh = 1.0 / hh;
    const unsigned int eOp1 = eOrder + 1;

    for (unsigned int k = 0; k < eOp1; k++) {
        double zz = eZ + k * hh;
        if (std::fabs(zz - zmin) < D_COMPAR_TOL) zz = zmin;
        if (std::fabs(zz - zmax) < D_COMPAR_TOL) zz = zmax;
        if (zz < zmin || zz > zmax) continue;
        const int kkz = (int)std::round((zz - zmin) * invhh);
        for (unsigned int j = 0; j < eOp1; j++) {
            double yy = eY + j * hh;
            if (std::fabs(yy - ymin) < D_COMPAR_TOL) yy = ymin;
            if (std::fabs(yy - ymax) < D_COMPAR_TOL) yy = ymax;
            if (yy < ymin || yy > ymax) continue;
            const int jjy = (int)std::round((yy - ymin) * invhh);
            for (unsigned int i = 0; i < eOp1; i++) {
                double xx = eX + i * hh;
                if (std::fabs(xx - xmin) < D_COMPAR_TOL) xx = xmin;
                if (std::fabs(xx - xmax) < D_COMPAR_TOL) xx = xmax;
                if (xx < xmin || xx > xmax) continue;
                const int iix = (int)std::round((xx - xmin) * invhh);
                if (iix < 0 || iix >= (int)lx || jjy < 0 || jjy >= (int)ly ||
                    kkz < 0 || kkz >= (int)lz)
                    continue;
                for (unsigned int v = 0; v < dof; v++)
                    uz[v * unSz + (std::size_t)kkz * lx * ly +
                       (std::size_t)jjy * lx + iix] =
                        dg[v * dgSz + k * eOp1 * eOp1 + j * eOp1 + i];
            }
        }
    }
}

// FP fallback, fine->coarse branch (element level == bLev + 1).
void fpFineToCoarse(const Geom& g, int ei, int ej, int ek, const double* dg,
                    double* uz, unsigned int dof, std::size_t unSz,
                    std::size_t dgSz) {
    const unsigned int eOrder = g.eOrder;
    const unsigned int lx = g.lx, ly = g.lx, lz = g.lx;
    const uint64_t szBlkEle = (uint64_t)1u << (g.maxDepth - g.bLev);
    const uint64_t szEle    = (uint64_t)1u << (g.maxDepth - (g.bLev + 1));
    const uint64_t blkSz    = szBlkEle * g.nEle;
    const double hx         = szBlkEle / (double)eOrder;
    const double xmin = g.blkX - g.PW * hx, xmax = g.blkX + blkSz + g.PW * hx;
    const double ymin = xmin, ymax = xmax, zmin = xmin, zmax = xmax;
    const double eX = (double)(g.blkX + (int64_t)ei * (int64_t)szEle);
    const double eY = (double)(g.blkX + (int64_t)ej * (int64_t)szEle);
    const double eZ = (double)(g.blkX + (int64_t)ek * (int64_t)szEle);
    const double hh       = szEle / (double)eOrder;
    const double invhh    = 1.0 / (2 * hh);
    const unsigned int cb = (eOrder % 2 == 0) ? 0 : 1;

    for (unsigned int k = cb; k < eOrder + 1; k += 2) {
        double zz = eZ + k * hh;
        if (std::fabs(zz - zmin) < D_COMPAR_TOL) zz = zmin;
        if (std::fabs(zz - zmax) < D_COMPAR_TOL) zz = zmax;
        if (zz < zmin || zz > zmax) continue;
        const int kkz = (int)std::round((zz - zmin) * invhh);
        for (unsigned int j = cb; j < eOrder + 1; j += 2) {
            double yy = eY + j * hh;
            if (std::fabs(yy - ymin) < D_COMPAR_TOL) yy = ymin;
            if (std::fabs(yy - ymax) < D_COMPAR_TOL) yy = ymax;
            if (yy < ymin || yy > ymax) continue;
            const int jjy = (int)std::round((yy - ymin) * invhh);
            for (unsigned int i = cb; i < eOrder + 1; i += 2) {
                double xx = eX + i * hh;
                if (std::fabs(xx - xmin) < D_COMPAR_TOL) xx = xmin;
                if (std::fabs(xx - xmax) < D_COMPAR_TOL) xx = xmax;
                if (xx < xmin || xx > xmax) continue;
                const int iix = (int)std::round((xx - xmin) * invhh);
                if (iix < 0 || iix >= (int)lx || jjy < 0 || jjy >= (int)ly ||
                    kkz < 0 || kkz >= (int)lz)
                    continue;
                for (unsigned int v = 0; v < dof; v++)
                    uz[v * unSz + (std::size_t)kkz * lx * ly +
                       (std::size_t)jjy * lx + iix] =
                        dg[v * dgSz + k * (eOrder + 1) * (eOrder + 1) +
                           j * (eOrder + 1) + i];
            }
        }
    }
}

enum Branch { SAME_LEVEL = 1, FINE_TO_COARSE = 2 };

int runCase(unsigned int eOrder, unsigned int nEle, Branch branch,
            unsigned int dof) {
    const Geom g        = makeGeom(eOrder, nEle);
    const std::size_t nPe =
        (std::size_t)(eOrder + 1) * (eOrder + 1) * (eOrder + 1);
    const std::size_t unSz = (std::size_t)g.lx * g.lx * g.lx;

    std::vector<double> dg(nPe * dof);
    for (std::size_t i = 0; i < dg.size(); i++) dg[i] = (double)i + 1.0;

    // past the block on both sides: pad, fully-outside, and partial overlaps
    const int lo = (branch == SAME_LEVEL) ? -3 : -6;
    const int hi = (branch == SAME_LEVEL) ? (int)nEle + 3 : 2 * (int)nEle + 6;

    int missing = 0, extra = 0, wrong = 0, placements = 0;
    for (int ek = lo; ek <= hi; ek++) {
        for (int ej = lo; ej <= hi; ej++) {
            for (int ei = lo; ei <= hi; ei++) {
                std::vector<double> ref(unSz * dof, SENTINEL);
                std::vector<double> got(unSz * dof, SENTINEL);

                if (branch == SAME_LEVEL) {
                    fpSameLevel(g, ei, ej, ek, dg.data(), ref.data(), dof, unSz,
                                nPe);
                    const int i0 = ei * (int)eOrder + (int)g.PW;
                    const int j0 = ej * (int)eOrder + (int)g.PW;
                    const int k0 = ek * (int)eOrder + (int)g.PW;
                    dendro::unzip::scatter_same_level_dispatch<double>(
                        dg.data(), got.data(), eOrder, dof, unSz, nPe, 0, g.lx,
                        g.lx, g.lx, i0, j0, k0);
                } else {
                    fpFineToCoarse(g, ei, ej, ek, dg.data(), ref.data(), dof,
                                   unSz, nPe);
                    const int halfEO = (int)eOrder / 2;
                    const int i0     = ei * halfEO + (int)g.PW;
                    const int j0     = ej * halfEO + (int)g.PW;
                    const int k0     = ek * halfEO + (int)g.PW;
                    dendro::unzip::scatter_fine_to_coarse_dispatch<double>(
                        dg.data(), got.data(), eOrder, dof, unSz, nPe, 0, g.lx,
                        g.lx, g.lx, i0, j0, k0);
                }
                placements++;

                for (std::size_t t = 0; t < ref.size(); t++) {
                    const bool wroteRef = (ref[t] != SENTINEL);
                    const bool wroteGot = (got[t] != SENTINEL);
                    if (wroteRef && !wroteGot) {
                        if (missing < 3)
                            std::printf(
                                "    MISSING  ei=%d ej=%d ek=%d flat=%zu\n", ei,
                                ej, ek, t);
                        missing++;
                    } else if (!wroteRef && wroteGot) {
                        if (extra < 3)
                            std::printf(
                                "    EXTRA    ei=%d ej=%d ek=%d flat=%zu\n", ei,
                                ej, ek, t);
                        extra++;
                    } else if (wroteRef && got[t] != ref[t]) {
                        if (wrong < 3)
                            std::printf(
                                "    VALUE    ei=%d ej=%d ek=%d flat=%zu "
                                "ref=%g got=%g\n",
                                ei, ej, ek, t, ref[t], got[t]);
                        wrong++;
                    }
                }
            }
        }
    }

    const int bad = missing + extra + wrong;
    std::printf(
        "  eO=%u nEle=%u PW=%u lx=%2u dof=%u %-13s %5d placements  "
        "missing=%d extra=%d value=%d  %s\n",
        eOrder, nEle, g.PW, g.lx, dof,
        branch == SAME_LEVEL ? "same-level" : "fine->coarse", placements,
        missing, extra, wrong, bad ? "FAIL" : "ok");
    return bad;
}

}  // namespace

int main() {
    std::printf("unzip scatter kernels vs FP fallback\n");
#if defined(DENDRO_TENSOR_SIMD)
    std::printf("  DENDRO_TENSOR_SIMD : ON\n");
#else
    std::printf("  DENDRO_TENSOR_SIMD : off\n");
#endif
#if defined(__AVX512F__)
    std::printf("  __AVX512F__        : yes\n");
#else
    std::printf("  __AVX512F__        : no\n");
#endif

    int bad = 0;
    for (unsigned int eOrder : {2u, 4u, 6u, 8u})
        for (unsigned int nEle : {1u, 2u, 4u})
            for (Branch br : {SAME_LEVEL, FINE_TO_COARSE})
                for (unsigned int dof : {1u, 3u})
                    bad += runCase(eOrder, nEle, br, dof);

    if (bad) {
        std::printf("\nFAILED: %d mismatching samples\n", bad);
        return 1;
    }
    std::printf("\nPASSED: fast path is bit-identical to the FP fallback\n");
    return 0;
}
