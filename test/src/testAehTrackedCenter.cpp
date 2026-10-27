// After a failed find, a BBH inspiral horizon must be searched for around its
// tracked puncture rather than its stale last-found center.
//
// Brill-Lindquist data (two unit-mass punctures, AH radius ~0.5, search
// radius 1). BH1 jumps by 2 between finds, out of its search sphere: that
// find fails, and the next one must recover it from the tracked location.
//
// Usage: testAehTrackedCenter [maxDepth=9]
// Exit code 0 = pass.

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "aeh_bhahaha.h"
#include "json.hpp"
#include "mesh.h"
#include "meshUtils.h"
#include "octUtils.h"

namespace {

constexpr double L      = 6.0;  // domain [-L, L]^3
constexpr double MASS   = 1.0;
constexpr double R_SRCH = 1.0;

using P3                = std::array<double, 3>;
// off the grid nodes, so no node sits on a puncture
const P3 BH1_OLD        = {2.1, 0.13, 0.07};
const P3 BH1_NEW        = {2.1, 2.13, 0.07};
const P3 BH2            = {-2.1, -0.11, 0.05};

double rad(const P3& p, double x, double y, double z) {
    return std::sqrt((x - p[0]) * (x - p[0]) + (y - p[1]) * (y - p[1]) +
                     (z - p[2]) * (z - p[2]));
}

struct State {
    bool failed;
    double t, x, y, z;
};

// horizon state as the checkpoint records it, on rank 0
std::vector<State> read_state(dendro_aeh::AEH_BHaHAHA& aeh,
                              const ot::Mesh* mesh, int rank) {
    const std::string fn = "testAehTrackedCenter_state.json";
    aeh.create_checkpoint(mesh, fn);
    std::vector<State> out;
    if (rank) return out;
    nlohmann::json j;
    std::ifstream(fn) >> j;
    std::remove(fn.c_str());
    const std::vector<bool> failed = j["FAILED_LAST_FIND"];
    const auto t = dendro_aeh::b91_decode<double>(j["T_M1"].get<std::string>());
    const auto x =
        dendro_aeh::b91_decode<double>(j["X_CENTER_M1"].get<std::string>());
    const auto y =
        dendro_aeh::b91_decode<double>(j["Y_CENTER_M1"].get<std::string>());
    const auto z =
        dendro_aeh::b91_decode<double>(j["Z_CENTER_M1"].get<std::string>());
    for (size_t h = 0; h < failed.size(); h++)
        out.push_back({failed[h], t[h], x[h], y[h], z[h]});
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    m_uiMaxDepth              = (argc > 1) ? std::atoi(argv[1]) : 9;
    const unsigned int eOrder = 6;
    _InitializeHcurve(m_uiDim);

    // graded toward every puncture position the test uses
    const double sc = 2.0 * L / (double)(1u << m_uiMaxDepth);
    std::vector<ot::TreeNode> oct, tmp;
    createRegularOctree(oct, 3, m_uiDim, m_uiMaxDepth, comm);
    for (unsigned int lev = 3; lev < m_uiMaxDepth; lev++) {
        tmp.clear();
        for (const ot::TreeNode& t : oct) {
            const double sz = (t.maxX() - t.minX()) * sc;
            const double cx = 0.5 * (t.minX() + t.maxX()) * sc - L;
            const double cy = 0.5 * (t.minY() + t.maxY()) * sc - L;
            const double cz = 0.5 * (t.minZ() + t.maxZ()) * sc - L;
            bool near       = false;
            for (const P3& p : {BH1_OLD, BH1_NEW, BH2})
                near = near || rad(p, cx, cy, cz) < 2.0 * sz;
            if (near && t.getLevel() == lev)
                t.addChildren(tmp);
            else
                tmp.push_back(t);
        }
        std::swap(oct, tmp);
    }
    ot::Mesh* mesh = ot::createMesh(oct.data(), oct.size(), eOrder, comm, 0);
    mesh->setDomainBounds(Point(-L, -L, -L), Point(L, L, L));
    DendroIntL nel = mesh->getNumLocalMeshElements(), nel_g = 0;
    MPI_Allreduce(&nel, &nel_g, 1, MPI_LONG_LONG, MPI_SUM, comm);

    // gxx..gzz, Kxx..Kzz with BH1 at bh1 (conformally flat, K = 0)
    P3 bh1          = BH1_OLD;
    const int NV    = 12;
    double* metric  = nullptr;
    auto set_metric = [&](const P3& p1) {
        bh1 = p1;
        if (metric) mesh->destroyVector(metric);
        std::function<void(double, double, double, double*)> f =
            [&bh1](double x, double y, double z, double* v) {
                const double psi = 1.0 + 0.5 * MASS / rad(bh1, x, y, z) +
                                   0.5 * MASS / rad(BH2, x, y, z);
                const double g   = psi * psi * psi * psi;
                for (int i = 0; i < NV; i++) v[i] = 0.0;
                v[0] = v[3] = v[5] = g;
            };
        metric = mesh->createCGVector<double>(f, NV);
        if (mesh->isActive()) {
            mesh->readFromGhostBegin(metric, NV);
            mesh->readFromGhostEnd(metric, NV);
        }
    };

    const Point grid_limits[2] = {
        Point(0.0, 0.0, 0.0),
        Point((double)(1u << m_uiMaxDepth), (double)(1u << m_uiMaxDepth),
              (double)(1u << m_uiMaxDepth))};
    const Point domain_limits[2] = {Point(-L, -L, -L), Point(L, L, L)};
    std::function<std::vector<double>(const std::vector<double>&)> transform =
        [](const std::vector<double>& in) { return in; };
    dendro_aeh::AEH_BHaHAHA aeh(
        3, true, {BH1_OLD[0], BH2[0], 0.0}, {BH1_OLD[1], BH2[1], 0.0},
        {BH1_OLD[2], BH2[2], 0.0}, 3, {}, {}, {10000, 10000, 10000}, {}, {}, {},
        {}, {R_SRCH, R_SRCH, 1.5}, {}, 32, 64, "./",
        {{BH1_OLD[0], BH1_OLD[1], BH1_OLD[2], MASS},
         {BH2[0], BH2[1], BH2[2], MASS}},
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, transform, grid_limits,
        domain_limits, /*file_output_freq*/ 0, 3, {8, 16, 32}, {16, 32, 64}, 0,
        /*verbosity*/ 1);

    int failures = 0;
    auto check   = [&](bool ok, const std::string& what) {
        if (!ok) {
            printf("  FAIL: %s\n", what.c_str());
            fflush(stdout);
            failures++;
        }
    };
    auto near = [](const State& s, const P3& p) {
        return std::fabs(s.x - p[0]) < 0.1 && std::fabs(s.y - p[1]) < 0.1 &&
               std::fabs(s.z - p[2]) < 0.1;
    };
    auto find = [&](unsigned int step, const P3& tracked1) {
        std::vector<const double*> v(NV);
        for (int i = 0; i < NV; i++)
            v[i] = metric + (size_t)i * mesh->getDegOfFreedom();
        aeh.find_horizons(mesh, v.data(), step, (double)step,
                          {Point(tracked1[0], tracked1[1], tracked1[2]),
                           Point(BH2[0], BH2[1], BH2[2])});
        return read_state(aeh, mesh, rank);
    };

    if (!rank)
        printf("testAehTrackedCenter npes=%d maxDepth=%u elements=%lld\n", npes,
               m_uiMaxDepth, (long long)nel_g);

    set_metric(BH1_OLD);
    auto s = find(0, BH1_OLD);
    if (!rank) {
        check(!s[0].failed && s[0].t == 0.0 && near(s[0], BH1_OLD),
              "BH1 must be found at its initial position");
        check(!s[1].failed && s[1].t == 0.0 && near(s[1], BH2),
              "BH2 must be found at t=0");
    }

    set_metric(BH1_NEW);
    s = find(1, BH1_NEW);
    if (!rank)
        check(s[0].failed && s[0].t == 0.0,
              "BH1 outside its search sphere must fail (test premise)");

    s = find(2, BH1_NEW);
    if (!rank) {
        check(!s[0].failed && s[0].t == 2.0 && near(s[0], BH1_NEW),
              "BH1 must be found at its tracked puncture after a failure");
        check(!s[1].failed && s[1].t == 2.0 && near(s[1], BH2),
              "BH2 must still be found");
        printf("testAehTrackedCenter: %s\n", failures ? "FAIL" : "PASS");
    }

    MPI_Bcast(&failures, 1, MPI_INT, 0, comm);
    mesh->destroyVector(metric);
    delete mesh;
    MPI_Finalize();
    return failures ? 1 : 0;
}
