// Gates for the AEH_BHaHAHA input guards, none of which need a mesh or a solve.
//
// checkpoint_root_rank(): the old `npesActive < 3 ? 0 : 3` named rank 3 when
// only ranks 0-2 exist, and MPI aborts on an invalid root.
//
// Constructor validation: each case below used to be an out-of-bounds read or
// write inside the solver instead of a diagnosable error. The two "must
// construct" cases are what stop the guards rejecting working setups.
//
// Exit code 0 = pass.

#include <mpi.h>
#include <unistd.h>

#include <cmath>
#include <csignal>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "aeh_bhahaha.h"

static int failures = 0;

static void check(bool ok, const std::string& what) {
    if (!ok) {
        // flushed: without the guards a later case segfaults and would
        // otherwise take the earlier output with it
        printf("  FAIL: %s\n", what.c_str());
        fflush(stdout);
        failures++;
    }
}

namespace {

// Valid by default; each case below breaks one field.
struct Cfg {
    unsigned int n_horizons          = 3;
    bool is_bbh                      = true;
    std::vector<double> xc           = {0.0, 0.0, 0.0};
    std::vector<double> yc           = {0.0, 0.0, 0.0};
    std::vector<double> zc           = {0.0, 0.0, 0.0};
    unsigned int n_res_multigrid     = 3;
    int n_res_after_find             = 3;
    std::vector<int> ntheta          = {8, 16, 32};
    std::vector<int> nphi            = {16, 32, 64};
    int ntheta_max                   = 32;
    int nphi_max                     = 64;
    std::vector<dendro_aeh::SimpleBlackHoleData> bhs = {
        {2.0, 0.0, 0.0, 0.5}, {-2.0, 0.0, 0.0, 0.5}};
};

bool construction_throws(const Cfg& c) {
    const Point grid_limits[2]   = {Point(0.0, 0.0, 0.0),
                                    Point(256.0, 256.0, 256.0)};
    const Point domain_limits[2] = {Point(-10.0, -10.0, -10.0),
                                    Point(10.0, 10.0, 10.0)};

    std::function<std::vector<double>(const std::vector<double>&)> transform =
        [](const std::vector<double>& in) { return in; };

    try {
        dendro_aeh::AEH_BHaHAHA aeh(
            c.n_horizons, c.is_bbh, c.xc, c.yc, c.zc, c.n_res_multigrid,
            /*m_scale*/ {}, /*cfl*/ {}, /*max_iter*/ {}, /*l2_tol*/ {},
            /*linf_tol*/ {}, /*eta_damp*/ {}, /*ko*/ {}, /*max_search_r*/ {},
            /*nr_interp_max*/ {}, c.ntheta_max, c.nphi_max, "./",
            c.bhs, /*indices_extract*/ {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
            transform, grid_limits, domain_limits, /*file_output_freq*/ 10,
            c.n_res_after_find, c.ntheta, c.nphi);
        return false;
    } catch (const std::runtime_error&) {
        return true;
    }
}

// BHaHAHA's radial grid set-up loops forever on some guesses
void on_alarm(int) {
    const char msg[] = "  FAIL: radial grid set-up hung\n";
    (void)!write(1, msg, sizeof(msg) - 1);
    _exit(1);
}

// what sanitize_radius_guess() hands to BHaHAHA must give a usable grid
bool radial_grid_ok(double r_min, double r_max) {
    constexpr int nr_max = 48;
    constexpr double R   = 1.0;
    dendro_aeh::AEH_BHaHAHA::sanitize_radius_guess(r_min, r_max, R);
    int nr    = -1;
    double r0 = NAN, dr = NAN;
    std::vector<double> radii(nr_max, NAN);
    alarm(5);
    bah_radial_grid_cell_centered_set_up(nr_max, R, r_min, r_max, &nr, &r0, &dr,
                                         radii.data());
    alarm(0);
    if (nr < 1 || nr > nr_max || !std::isfinite(r0) || !(dr > 0.0) ||
        !std::isfinite(dr))
        return false;
    for (int i = 0; i < nr; i++)
        if (!std::isfinite(radii[i]) || radii[i] < 0.0) return false;
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        // 3 was the only wrong value, but pin its neighbours too
        for (unsigned int npes = 1; npes <= 8; npes++) {
            const unsigned int root =
                dendro_aeh::AEH_BHaHAHA::checkpoint_root_rank(npes);
            check(root < npes, "checkpoint root " + std::to_string(root) +
                                   " is outside a communicator of size " +
                                   std::to_string(npes));
        }
        check(dendro_aeh::AEH_BHaHAHA::checkpoint_root_rank(3) == 0,
              "3 ranks must collect on rank 0");
        check(dendro_aeh::AEH_BHaHAHA::checkpoint_root_rank(4) == 3,
              "4+ ranks must still collect on rank 3 (unchanged behaviour)");

        {
            Cfg c;
            c.ntheta_max = 16;
            check(construction_throws(c),
                  "final ntheta != NTHETA_MAX must throw");
        }
        {
            Cfg c;
            c.nphi_max = 32;
            check(construction_throws(c), "final nphi != NPHI_MAX must throw");
        }
        {  // equal totals, transposed dims
            Cfg c;
            c.ntheta     = {8, 16, 64};
            c.nphi       = {16, 32, 32};
            c.ntheta_max = 64;
            c.nphi_max   = 32;
            check(!construction_throws(c),
                  "a consistent 64x32 ladder must construct");
            c.ntheta_max = 32;
            c.nphi_max   = 64;
            check(construction_throws(c),
                  "transposed theta/phi with equal totals must throw");
        }
        {
            Cfg c;
            c.n_horizons = 2;
            c.xc = c.yc = c.zc = {0.0, 0.0};
            check(construction_throws(c), "BBH with n_horizons < 3 must throw");
        }
        {
            Cfg c;
            c.n_horizons = 4;
            c.ntheta_max = 32;
            check(construction_throws(c),
                  "n_horizons beyond the center-list length must throw");
        }
        {
            Cfg c;
            c.bhs = {{2.0, 0.0, 0.0, 0.5}};
            check(construction_throws(c),
                  "BBH with fewer than 2 black holes must throw");
        }
        {
            Cfg c;
            c.n_res_multigrid  = 20;
            c.n_res_after_find = 20;
            c.ntheta.assign(20, 32);
            c.nphi.assign(20, 64);
            check(construction_throws(c),
                  "more than MAX_RESOLUTIONS levels must throw");
        }
        {
            Cfg c;
            c.n_res_multigrid  = 0;
            c.n_res_after_find = 0;
            c.ntheta.clear();
            c.nphi.clear();
            check(construction_throws(c), "zero multigrid levels must throw");
        }
        {
            Cfg c;
            c.n_res_multigrid = 1;
            check(construction_throws(c),
                  "n_resolutions_multigrid contradicting "
                  "num_resolutions_after_find must throw");
        }
        {
            Cfg c;
            c.ntheta = {8, 0, 32};
            check(construction_throws(c),
                  "a non-positive ntheta entry must throw");
        }
        {  // MAX_RESOLUTIONS itself must be accepted, one past it rejected
            Cfg c;
            c.n_res_multigrid = c.n_res_after_find = MAX_RESOLUTIONS;
            c.ntheta.assign(MAX_RESOLUTIONS, 32);
            c.nphi.assign(MAX_RESOLUTIONS, 64);
            check(!construction_throws(c),
                  "exactly MAX_RESOLUTIONS levels must construct");
            c.n_res_multigrid = c.n_res_after_find = MAX_RESOLUTIONS + 1;
            c.ntheta.assign(MAX_RESOLUTIONS + 1, 32);
            c.nphi.assign(MAX_RESOLUTIONS + 1, 64);
            check(construction_throws(c),
                  "MAX_RESOLUTIONS + 1 levels must throw");
        }
        {  // a single level is legal, so the bound must not over-reject
            Cfg c;
            c.n_res_multigrid = c.n_res_after_find = 1;
            c.ntheta                               = {32};
            c.nphi                                 = {64};
            check(!construction_throws(c),
                  "a single multigrid level must construct");
        }
        {  // non-BBH is still indexed to n_horizons
            Cfg c;
            c.is_bbh     = false;
            c.n_horizons = 3;
            c.xc = c.yc = c.zc = {0.0, 0.0};
            c.bhs              = {};
            check(construction_throws(c),
                  "non-BBH short center lists must throw");
        }
        {
            Cfg c;
            c.n_horizons = 0;
            check(construction_throws(c), "zero horizons must throw");
        }
        {
            Cfg c;
            check(!construction_throws(c),
                  "the default valid configuration must construct");
        }
        {  // the CCZ4/Z4c case: BBH rules must not reject it
            Cfg c;
            c.is_bbh     = false;
            c.n_horizons = 1;
            c.xc = c.yc = c.zc = {0.0};
            c.bhs              = {};
            check(!construction_throws(c),
                  "a single-horizon non-BBH configuration must construct");
        }

        {  // radius guesses: kept when usable, else the full search sphere
            const double nan = NAN, inf = INFINITY;
            const double kept[][2] = {
                {0.3, 0.8}, {-0.1, 0.5}, {0.5, 0.5}, {0.2, 5.0}, {0.0, 0.3}};
            for (const auto& g : kept) {
                double r0 = g[0], r1 = g[1];
                check(!dendro_aeh::AEH_BHaHAHA::sanitize_radius_guess(r0, r1,
                                                                      1.0) &&
                          r0 == g[0] && r1 == g[1],
                      "usable radius guess (" + std::to_string(g[0]) + ", " +
                          std::to_string(g[1]) + ") must be kept");
            }
            const double reset[][2] = {{nan, 1.0},  {0.2, nan},   {inf, 1.0},
                                       {0.2, inf},  {2.0, 0.0},   {1e17, 1.0},
                                       {0.5, -1.0}, {-1.0, -0.5}, {0.9, 0.1},
                                       {1.0, 2.0},  {0.0, 0.0}};
            for (const auto& g : reset) {
                double r0 = g[0], r1 = g[1];
                const std::string tag = "(" + std::to_string(g[0]) + ", " +
                                        std::to_string(g[1]) + ")";
                check(dendro_aeh::AEH_BHaHAHA::sanitize_radius_guess(r0, r1,
                                                                     1.0) &&
                          r0 == 0.0 && r1 == 1.0,
                      "radius guess " + tag + " must become the full sphere");
            }
            std::signal(SIGALRM, on_alarm);
            for (const auto& g : kept)
                check(radial_grid_ok(g[0], g[1]),
                      "radial grid from a kept guess");
            for (const auto& g : reset)
                check(radial_grid_ok(g[0], g[1]),
                      "radial grid from a sanitized guess (" +
                          std::to_string(g[0]) + ", " + std::to_string(g[1]) +
                          ")");
        }

        if (failures == 0) {
            printf("testAehBhahahaGuards: PASS\n");
        } else {
            printf("testAehBhahahaGuards: FAIL (%d checks failed)\n", failures);
        }
    }

    MPI_Bcast(&failures, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    return failures == 0 ? 0 : 1;
}
