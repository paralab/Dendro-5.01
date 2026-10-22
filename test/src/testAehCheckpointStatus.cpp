// AEH_BHaHAHA::create_checkpoint must report a checkpoint that was not fully
// written, on every active rank, so the caller can keep the previous one.
//
// Exit code 0 = pass.

#include <mpi.h>
#include <unistd.h>

#include <climits>
#include <cstdio>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "aeh_bhahaha.h"
#include "mesh.h"
#include "meshUtils.h"
#include "octUtils.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    m_uiMaxDepth = 6;
    _InitializeHcurve(m_uiDim);
    std::vector<ot::TreeNode> oct;
    createRegularOctree(oct, 3, m_uiDim, m_uiMaxDepth, comm);
    // small grain so up to 8 ranks stay active and rank 3 can be the writer
    ot::Mesh* mesh = ot::createMesh(oct.data(), oct.size(), 4, comm, 0,
                                    ot::SM_TYPE::FDM, 64);

    const Point grid_limits[2]   = {Point(0.0, 0.0, 0.0),
                                    Point(64.0, 64.0, 64.0)};
    const Point domain_limits[2] = {Point(-10.0, -10.0, -10.0),
                                    Point(10.0, 10.0, 10.0)};
    std::function<std::vector<double>(const std::vector<double>&)> transform =
        [](const std::vector<double>& in) { return in; };
    dendro_aeh::AEH_BHaHAHA aeh(
        3, true, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 3, {}, {},
        {}, {}, {}, {}, {}, {}, {}, 32, 64, "./",
        {{2.0, 0.0, 0.0, 0.5}, {-2.0, 0.0, 0.0, 0.5}},
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, transform, grid_limits,
        domain_limits, 10, 3, {8, 16, 32}, {16, 32, 64});

    int failures = 0;
    // every active rank must see the same status
    auto check   = [&](const char* path, bool want_fail, const char* what) {
        const int st = aeh.create_checkpoint(mesh, path);
        int lo       = mesh->isActive() ? st : INT_MAX;
        int hi       = mesh->isActive() ? st : INT_MIN;
        MPI_Allreduce(MPI_IN_PLACE, &lo, 1, MPI_INT, MPI_MIN, comm);
        MPI_Allreduce(MPI_IN_PLACE, &hi, 1, MPI_INT, MPI_MAX, comm);
        const bool ok = (lo == hi) && ((lo != 0) == want_fail);
        if (!ok) {
            if (!rank) printf("  FAIL: %s (status %d..%d)\n", what, lo, hi);
            failures++;
        }
    };

    const std::string good = "testAehCheckpointStatus_ok.json";
    check(good.c_str(), false, "a writable path must return 0");
    if (!rank) {
        std::ifstream in(good);
        if (!in) {
            printf("  FAIL: %s was not written\n", good.c_str());
            failures++;
        }
    }
    check("testAehCheckpointStatus_no_such_dir/chkpt.json", true,
          "an unopenable path must return nonzero");
    if (access("/dev/full", W_OK) == 0)
        check("/dev/full", true, "a failed write must return nonzero");

    MPI_Bcast(&failures, 1, MPI_INT, 0, comm);
    if (!rank) {
        std::remove(good.c_str());
        printf("testAehCheckpointStatus npes=%d: %s\n", npes,
               failures ? "FAIL" : "PASS");
    }

    delete mesh;
    MPI_Finalize();
    return failures ? 1 : 0;
}
