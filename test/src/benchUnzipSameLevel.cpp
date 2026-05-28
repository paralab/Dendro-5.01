// Focused micro-benchmark for the same-level branch of Mesh::unzip_scatter.
//
// Mirrors the BSSN solver's per-variable calling pattern: BSSN_NUM_VARS_SIM
// sequential dof=1 unzip_scatter calls per "step". The bench builds a mesh
// representative of typical BSSN-GR usage (depth=8, wavelet_tol=1e-4,
// eleOrder=4) and times the call.
//
// For bit-exact verification across builds, after the timed loop the bench
// writes the last unzipped buffer to a binary file. Compile once with
// DENDRO_UNZIP_SCATTER_FAST=OFF, run, save the output, then compile with
// DENDRO_UNZIP_SCATTER_FAST=ON, run, and `cmp` the two output files.
//
// Usage:
//   mpirun -n <N> ./benchUnzipSameLevel <maxDepth> <wavelet_tol>
//   <partition_tol>
//                                       <eleOrder> [iters] [num_vars]
//                                       [out_file]

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "TreeNode.h"
#include "dendro.h"
#include "mesh.h"
#include "meshTestUtils.h"
#include "meshUtils.h"
#include "mpi.h"
#include "profiler.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if (argc < 5) {
        if (!rank) {
            std::fprintf(stderr,
                         "Usage: %s maxDepth wavelet_tol partition_tol "
                         "eleOrder [iters=32] [num_vars=24] [out_file=]\n",
                         argv[0]);
        }
        MPI_Abort(comm, 1);
    }

    m_uiMaxDepth         = std::atoi(argv[1]);
    double wavelet_tol   = std::atof(argv[2]);
    double partition_tol = std::atof(argv[3]);
    unsigned int eOrder  = (unsigned int)std::atoi(argv[4]);
    unsigned int iter    = (argc > 5) ? (unsigned int)std::atoi(argv[5]) : 32u;
    unsigned int n_vars  = (argc > 6) ? (unsigned int)std::atoi(argv[6]) : 24u;
    std::string out_file = (argc > 7) ? std::string(argv[7]) : std::string();

#if defined(DENDRO_UNZIP_SCATTER_FAST)
    const char* fast_tag = "FAST=ON";
#else
    const char* fast_tag = "FAST=OFF";
#endif

    if (!rank) {
        std::printf(
            "benchUnzipSameLevel [%s] maxDepth=%d wavelet_tol=%g "
            "partition_tol=%g eleOrder=%u iter=%u n_vars=%u\n",
            fast_tag, m_uiMaxDepth, wavelet_tol, partition_tol, eOrder, iter,
            n_vars);
    }

    _InitializeHcurve(m_uiDim);

    const double d_min = -10.0;
    const double d_max = 10.0;
    Point pt_min(d_min, d_min, d_min);
    Point pt_max(d_max, d_max, d_max);

    // Gaussian-sum (same shape as benchUtils.cpp) so refinement places
    // elements at multiple levels — exercises both same-level and refinement
    // branches of unzip_scatter.
    std::function<void(double, double, double, double*)> func =
        [](double x, double y, double z, double* var) {
            const double ca[] = {-2.0, 0.0, 0.0};
            const double cb[] = {2.0, 0.0, 0.0};
            const double rra  = (x - ca[0]) * (x - ca[0]) +
                                (y - ca[1]) * (y - ca[1]) +
                                (z - ca[2]) * (z - ca[2]);
            const double rrb  = (x - cb[0]) * (x - cb[0]) +
                                (y - cb[1]) * (y - cb[1]) +
                                (z - cb[2]) * (z - cb[2]);
            var[0]            = std::exp(-rra) + std::exp(-rrb);
        };

    std::function<double(double, double, double)> fr = [func, d_min, d_max](
                                                           double x, double y,
                                                           double z) {
        const double xx = (x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min;
        const double yy = (y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min;
        const double zz = (z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min;
        double v;
        func(xx, yy, zz, &v);
        return v;
    };

    std::vector<ot::TreeNode> tmpNodes;
    function2Octree(fr, tmpNodes, m_uiMaxDepth, wavelet_tol, eOrder, comm);

    ot::Mesh* mesh = ot::createMesh(
        tmpNodes.data(), tmpNodes.size(), eOrder, comm, 1, ot::SM_TYPE::FDM,
        DENDRO_DEFAULT_GRAIN_SZ, partition_tol, DENDRO_DEFAULT_SF_K);
    mesh->setDomainBounds(pt_min, pt_max);

    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin, lmax);
    if (!rank) std::printf("mesh lev_min=%u lev_max=%u\n", lmin, lmax);

    // One CG vector + one unzipped vector for measurement.
    double* u               = mesh->createCGVector<double>(func, 1);
    double* u_unzip         = mesh->createUnZippedVector<double>(1);

    const unsigned int unSz = mesh->getDegOfFreedomUnZip();

    mesh->readFromGhostBegin(u, 1);
    mesh->readFromGhostEnd(u, 1);

    // correctness sanity (analytic-field tolerance check)
    std::memset(u_unzip, 0, (size_t)unSz * sizeof(double));
    mesh->unzip_scatter(u, u_unzip, 1);
    const bool valid = ot::test::isUnzipValid(mesh, u_unzip, fr, 1e-3);
    if (!rank)
        std::printf("[%s] unzip_scatter correctness vs analytic: %s\n",
                    fast_tag, valid ? "PASS" : "FAIL");

    // Warm up.
    for (unsigned int w = 0; w < 4; ++w) {
        std::memset(u_unzip, 0, (size_t)unSz * sizeof(double));
        for (unsigned int v = 0; v < n_vars; ++v)
            mesh->unzip_scatter(u, u_unzip, 1);
    }

    profiler_t t_unzip;
    t_unzip.start();
    for (unsigned int i = 0; i < iter; ++i) {
        // Zero out the buffer at the top of each iteration so the comparison
        // file is deterministic — untouched padding cells stay zero.
        std::memset(u_unzip, 0, (size_t)unSz * sizeof(double));
        for (unsigned int v = 0; v < n_vars; ++v)
            mesh->unzip_scatter(u, u_unzip, 1);
    }
    t_unzip.stop();

    double t_local = t_unzip.seconds / (double)iter;
    double t_stat[3];
    par::computeOverallStats(&t_local, t_stat, comm,
                             "unzip_scatter per step (n_vars calls)");

    // Per-call (per-variable) average for direct comparison to other benches.
    double t_local_per_call = t_local / (double)n_vars;
    par::computeOverallStats(&t_local_per_call, t_stat, comm,
                             "unzip_scatter per call (dof=1)");

    // Write the last unzipped buffer to disk for bit-exact compare across
    // builds. Only rank 0 writes; the test should run with -n 1 for the
    // bit-exact check (other ranks have different partitions).
    if (!out_file.empty() && rank == 0) {
        std::ofstream ofs(out_file, std::ios::binary | std::ios::trunc);
        ofs.write(reinterpret_cast<const char*>(u_unzip),
                  (std::streamsize)((size_t)unSz * sizeof(double)));
        ofs.close();
        std::printf("[%s] wrote %zu bytes to %s\n", fast_tag,
                    (size_t)unSz * sizeof(double), out_file.c_str());
    }

    mesh->destroyVector(u);
    mesh->destroyVector(u_unzip);
    delete mesh;

    MPI_Finalize();
    return 0;
}
