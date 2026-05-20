// Focused test/bench for Mesh::isReMeshUnzip — the AMR refinement
// criterion that's called every K timesteps to decide where to refine
// or coarsen the octree. Computes wavelet coefficients per element
// across all evolved variables and flags elements that exceed
// tolerance.
//
// This test:
//   1. Builds a mesh on a Gaussian-sum field.
//   2. Creates an array of n_vars unzipped vectors (each = same field).
//   3. Calls isReMeshUnzip many times in a tight loop.
//   4. Reports return value (a global "mesh changed" bool) and timing.
//
// Compare across OMP=OFF, OMP=ON, and OMP thread counts. The return
// value must be identical (refinement decision is deterministic).

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
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
        if (!rank)
            std::fprintf(stderr,
                         "Usage: %s maxDepth wavelet_tol partition_tol "
                         "eleOrder [iters=16] [n_vars=24]\n",
                         argv[0]);
        MPI_Abort(comm, 1);
    }

    m_uiMaxDepth         = std::atoi(argv[1]);
    double wavelet_tol_v = std::atof(argv[2]);
    double partition_tol = std::atof(argv[3]);
    unsigned int eOrder  = (unsigned int)std::atoi(argv[4]);
    unsigned int iter    = (argc > 5) ? (unsigned int)std::atoi(argv[5]) : 16u;
    unsigned int n_vars  = (argc > 6) ? (unsigned int)std::atoi(argv[6]) : 24u;

#if defined(DENDRO_UNZIP_OMP)
    const char* omp_tag = "OMP=ON";
#else
    const char* omp_tag = "OMP=OFF";
#endif
    if (!rank)
        std::printf("testIsReMeshUnzip [%s] maxDepth=%d wavelet_tol=%g "
                    "eleOrder=%u iter=%u n_vars=%u\n",
                    omp_tag, m_uiMaxDepth, wavelet_tol_v, eOrder, iter,
                    n_vars);

    _InitializeHcurve(m_uiDim);
    const double d_min = -10.0;
    const double d_max = 10.0;
    Point pt_min(d_min, d_min, d_min);
    Point pt_max(d_max, d_max, d_max);

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
            var[0] = std::exp(-rra) + std::exp(-rrb);
        };
    std::function<double(double, double, double)> fr =
        [func, d_min, d_max](double x, double y, double z) {
            const double xx =
                (x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min;
            const double yy =
                (y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min;
            const double zz =
                (z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min;
            double v;
            func(xx, yy, zz, &v);
            return v;
        };

    std::vector<ot::TreeNode> tmpNodes;
    function2Octree(fr, tmpNodes, m_uiMaxDepth, wavelet_tol_v, eOrder, comm);
    ot::Mesh* mesh = ot::createMesh(
        tmpNodes.data(), tmpNodes.size(), eOrder, comm, 1, ot::SM_TYPE::FDM,
        DENDRO_DEFAULT_GRAIN_SZ, partition_tol, DENDRO_DEFAULT_SF_K);
    mesh->setDomainBounds(pt_min, pt_max);

    const size_t unSz = mesh->getDegOfFreedomUnZip();

    // Build n_vars unzipped vectors (each = same field for simplicity).
    std::vector<double*> all_unzip(n_vars);
    std::vector<double*> u_cg(n_vars);
    for (unsigned int v = 0; v < n_vars; v++) {
        u_cg[v]      = mesh->createCGVector<double>(func, 1);
        all_unzip[v] = mesh->createUnZippedVector<double>(1);
        mesh->readFromGhostBegin(u_cg[v], 1);
        mesh->readFromGhostEnd(u_cg[v], 1);
        mesh->unzip(u_cg[v], all_unzip[v], 1);
    }
    std::vector<const double*> unzippedVec(n_vars);
    for (unsigned int v = 0; v < n_vars; v++) unzippedVec[v] = all_unzip[v];

    std::vector<unsigned int> varIds(n_vars);
    for (unsigned int v = 0; v < n_vars; v++) varIds[v] = v;

    // wavelet tolerance callback (constant tolerance).
    std::function<double(double, double, double, double*)> wtol =
        [wavelet_tol_v](double, double, double, double*) {
            return wavelet_tol_v;
        };

    // Warm-up + correctness check (single call).
    bool ret0 = mesh->isReMeshUnzip(unzippedVec.data(), varIds.data(), n_vars,
                                    wtol);
    if (rank == 0)
        std::printf("[%s] isReMeshUnzip returned: %s\n", omp_tag,
                    ret0 ? "true" : "false");

    // Timed loop.
    profiler_t t_amr;
    t_amr.start();
    bool any_changed = false;
    for (unsigned int i = 0; i < iter; i++) {
        bool r = mesh->isReMeshUnzip(unzippedVec.data(), varIds.data(),
                                     n_vars, wtol);
        any_changed = any_changed || r;
    }
    t_amr.stop();

    double t_local = t_amr.seconds / (double)iter;
    double t_stat[3];
    par::computeOverallStats(&t_local, t_stat, comm,
                             "isReMeshUnzip per call");

    if (rank == 0)
        std::printf("[%s] any iteration changed mesh? %s\n", omp_tag,
                    any_changed ? "true" : "false");

    // Also output a checksum of the input data to verify the same input was
    // processed. (isReMeshUnzip itself returns only a bool, so we can't
    // directly hash its internal eleWMax. We can however confirm setup is
    // deterministic.)
    if (rank == 0) {
        double cs = 0.0;
        for (unsigned int v = 0; v < n_vars; v++)
            for (size_t i = 0; i < unSz; i++)
                cs += all_unzip[v][i] * (double)(i + 1);
        std::printf("[%s] unzipped-input checksum=%23.17e\n", omp_tag, cs);
    }

    for (unsigned int v = 0; v < n_vars; v++) {
        mesh->destroyVector(u_cg[v]);
        mesh->destroyVector(all_unzip[v]);
    }
    delete mesh;
    MPI_Finalize();
    return 0;
}
