// Microbench for Mesh::unzip() (the GATHER variant) — what BSSN actually
// calls per evolved variable per RK substage.
//
// Mirrors test/src/benchUnzipSameLevel.cpp but exercises unzip() instead of
// unzip_scatter(). Use to (a) establish a baseline on the BSSN-relevant
// path and (b) profile with perf to find the next hotspots.

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
    const char* fast_tag = "SCATTER_FAST=ON";
#else
    const char* fast_tag = "SCATTER_FAST=OFF";
#endif

    if (!rank) {
        std::printf("benchUnzipGather [%s] maxDepth=%d wavelet_tol=%g "
                    "partition_tol=%g eleOrder=%u iter=%u n_vars=%u\n",
                    fast_tag, m_uiMaxDepth, wavelet_tol, partition_tol, eOrder,
                    iter, n_vars);
    }

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
    function2Octree(fr, tmpNodes, m_uiMaxDepth, wavelet_tol, eOrder, comm);

    ot::Mesh* mesh = ot::createMesh(
        tmpNodes.data(), tmpNodes.size(), eOrder, comm, 1, ot::SM_TYPE::FDM,
        DENDRO_DEFAULT_GRAIN_SZ, partition_tol, DENDRO_DEFAULT_SF_K);
    mesh->setDomainBounds(pt_min, pt_max);

    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin, lmax);
    if (!rank)
        std::printf("mesh lev_min=%u lev_max=%u\n", lmin, lmax);

    double* u       = mesh->createCGVector<double>(func, 1);
    double* u_unzip = mesh->createUnZippedVector<double>(1);

    const unsigned int unSz = mesh->getDegOfFreedomUnZip();

    mesh->readFromGhostBegin(u, 1);
    mesh->readFromGhostEnd(u, 1);

    // sanity (analytic-field tolerance)
    std::memset(u_unzip, 0, (size_t)unSz * sizeof(double));
    mesh->unzip(u, u_unzip, 1);
    const bool valid = ot::test::isUnzipValid(mesh, u_unzip, fr, 1e-3);
    if (!rank)
        std::printf("[%s] unzip correctness vs analytic: %s\n", fast_tag,
                    valid ? "PASS" : "FAIL");

    // warm-up
    for (unsigned int w = 0; w < 4; ++w) {
        std::memset(u_unzip, 0, (size_t)unSz * sizeof(double));
        for (unsigned int v = 0; v < n_vars; ++v)
            mesh->unzip(u, u_unzip, 1);
    }

    profiler_t t_unzip;
    t_unzip.start();
    for (unsigned int i = 0; i < iter; ++i) {
        std::memset(u_unzip, 0, (size_t)unSz * sizeof(double));
        for (unsigned int v = 0; v < n_vars; ++v)
            mesh->unzip(u, u_unzip, 1);
    }
    t_unzip.stop();

    double t_local = t_unzip.seconds / (double)iter;
    double t_stat[3];
    par::computeOverallStats(&t_local, t_stat, comm,
                             "unzip (gather) per step");

    double t_local_per_call = t_local / (double)n_vars;
    par::computeOverallStats(&t_local_per_call, t_stat, comm,
                             "unzip (gather) per call (dof=1)");

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
