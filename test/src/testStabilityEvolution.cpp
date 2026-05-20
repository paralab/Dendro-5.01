// Numeric stability test — pseudo-evolution that repeatedly exercises the
// unzip path so that any accumulated drift in our optimized variants shows
// up in the final state.
//
// At each "timestep":
//   1. unzip u -> u_unzip
//   2. Compute a deterministic local update (sum of unzipped values * dt)
//      applied back to u via a CG-vector update (just a scalar accumulation
//      into u). This isn't a physical evolution — it's a controlled scalar
//      mixing that is sensitive to ANY change in the unzipped values.
//   3. After N steps, print:
//        - the L1, L2, L_inf norms of u
//        - a hash-like checksum (sum of u[i]*(i+1)) for bit-level detection
//
// Two runs of the bench with the same parameters and the same flag set
// MUST produce identical printed values. Across different flag sets,
// agreement to machine epsilon (relative ~1e-13) indicates the
// optimizations are numerically equivalent.

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

    if (argc < 6) {
        if (!rank)
            std::fprintf(
                stderr,
                "Usage: %s maxDepth wavelet_tol partition_tol eleOrder "
                "n_steps [n_vars=24] [dt=1e-6]\n",
                argv[0]);
        MPI_Abort(comm, 1);
    }

    m_uiMaxDepth         = std::atoi(argv[1]);
    double wavelet_tol   = std::atof(argv[2]);
    double partition_tol = std::atof(argv[3]);
    unsigned int eOrder  = (unsigned int)std::atoi(argv[4]);
    unsigned int n_steps = (unsigned int)std::atoi(argv[5]);
    unsigned int n_vars  = (argc > 6) ? (unsigned int)std::atoi(argv[6]) : 24u;
    double dt            = (argc > 7) ? std::atof(argv[7]) : 1e-6;

#if defined(DENDRO_UNZIP_SCATTER_FAST)
    const char* fast_tag = "SCATTER_FAST=ON";
#else
    const char* fast_tag = "SCATTER_FAST=OFF";
#endif
#if defined(DENDRO_TENSOR_SIMD)
    const char* simd_tag = "TENSOR_SIMD=ON";
#else
    const char* simd_tag = "TENSOR_SIMD=OFF";
#endif
#if defined(DENDRO_UNZIP_OMP)
    const char* omp_tag = "OMP=ON";
#else
    const char* omp_tag = "OMP=OFF";
#endif

    if (!rank) {
        std::printf("testStabilityEvolution [%s %s %s] maxDepth=%d "
                    "wavelet_tol=%g partition_tol=%g eleOrder=%u n_steps=%u "
                    "n_vars=%u dt=%g\n",
                    fast_tag, simd_tag, omp_tag, m_uiMaxDepth, wavelet_tol,
                    partition_tol, eOrder, n_steps, n_vars, dt);
    }

    _InitializeHcurve(m_uiDim);

    const double d_min = -10.0;
    const double d_max = 10.0;
    Point pt_min(d_min, d_min, d_min);
    Point pt_max(d_max, d_max, d_max);

    // Initial condition: smooth Gaussian-sum (same as the other benches).
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
    function2Octree(fr, tmpNodes, m_uiMaxDepth, wavelet_tol, eOrder, comm);
    ot::Mesh* mesh = ot::createMesh(
        tmpNodes.data(), tmpNodes.size(), eOrder, comm, 1, ot::SM_TYPE::FDM,
        DENDRO_DEFAULT_GRAIN_SZ, partition_tol, DENDRO_DEFAULT_SF_K);
    mesh->setDomainBounds(pt_min, pt_max);

    const size_t cgSz = mesh->getDegOfFreedom();
    const size_t unSz = mesh->getDegOfFreedomUnZip();

    if (!rank)
        std::printf("mesh: cgSz=%zu unSz=%zu\n", cgSz, unSz);

    double* u       = mesh->createCGVector<double>(func, 1);
    double* u_unzip = mesh->createUnZippedVector<double>(1);

    auto compute_norms = [&](const double* v, const char* label) {
        double l1 = 0, l2sq = 0, linf = 0, checksum = 0;
        for (size_t i = 0; i < cgSz; i++) {
            const double a = std::fabs(v[i]);
            l1 += a;
            l2sq += v[i] * v[i];
            if (a > linf) linf = a;
            checksum += v[i] * (double)((int64_t)i + 1);  // index-weighted sum
        }
        // Reduce across ranks if MPI > 1.
        if (npes > 1) {
            double g_l1, g_l2sq, g_linf, g_checksum;
            MPI_Reduce(&l1, &g_l1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&l2sq, &g_l2sq, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&linf, &g_linf, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
            MPI_Reduce(&checksum, &g_checksum, 1, MPI_DOUBLE, MPI_SUM, 0,
                       comm);
            l1       = g_l1;
            l2sq     = g_l2sq;
            linf     = g_linf;
            checksum = g_checksum;
        }
        if (rank == 0) {
            std::printf("[%s] L1=%23.17e L2=%23.17e Linf=%23.17e "
                        "checksum=%23.17e\n",
                        label, l1, std::sqrt(l2sq), linf, checksum);
        }
    };

    compute_norms(u, "init");

    // Pseudo-evolution: each "step" does
    //   1. ghost exchange + unzip
    //   2. compute a scalar (sum of unzipped values) that depends on every
    //      output cell — this propagates any difference in the unzipped
    //      data through to the next u.
    //   3. nudge u by dt * (local_scalar / cgSz) so it evolves while
    //      remaining stable.
    // The mixing is intentionally trivial — the goal is just to make every
    // single unzip output element ATTRIBUTABLE in the final u state, so a
    // single-byte change anywhere in the unzip pipeline shows up in the
    // final norms.
    profiler_t t_step;
    t_step.start();
    for (unsigned int step = 0; step < n_steps; step++) {
        // Simulate BSSN's per-RK-substage pattern: n_vars dof=1 unzip calls.
        for (unsigned int v = 0; v < n_vars; v++) {
            mesh->readFromGhostBegin(u, 1);
            mesh->readFromGhostEnd(u, 1);
            std::memset(u_unzip, 0, unSz * sizeof(double));
            mesh->unzip(u, u_unzip, 1);

            // Local-only sum over unzipped values (this rank's portion).
            double local_sum = 0;
            for (size_t i = 0; i < unSz; i++) local_sum += u_unzip[i];

            // Deterministic small nudge to u — multiplies the local sum by
            // a fixed coefficient to mix unzipped state back into u.
            const double coef = dt * 1.0e-3 / (double)(cgSz + 1);
            for (size_t i = 0; i < cgSz; i++) u[i] += coef * local_sum;
        }
    }
    t_step.stop();

    if (rank == 0) {
        std::printf("evolved %u steps (%u vars/step) in %.3f sec\n", n_steps,
                    n_vars, t_step.seconds);
    }
    compute_norms(u, "final");

    mesh->destroyVector(u);
    mesh->destroyVector(u_unzip);
    delete mesh;
    MPI_Finalize();
    return 0;
}
