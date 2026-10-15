// Embedded-pair error estimate and adaptive step control for ts::ETS.
//
// Tableau level (scalar ODEs, no mesh), per embedded pair:
//   - error row e = b - bhat sums to 0
//   - propagated row b converges at order 5, embedded row bhat at order 4
//   - one-step estimate h*sum(e_i k_i) scales like h^5 and matches the true
//     local error of bhat
// ETS level (real octree, real ETS::evolve, any rank count), on
//   du/dt = -K (u - a(x) cos t), u(0) = 0  (fast transient, then smooth):
//   - ESTIMATE leaves the solution bit-identical to OFF; norm scales like h^5
//   - ADAPTIVE meets the tolerance, rejects an oversized first step, grows dt
//     after the transient, counts accepted steps, advances t by dt taken
//   - the dt cap holds and is reached
//   - unsupported schemes (RK4, MSRK) and settings are refused
// Run: ./testETSErrorControl   (exit 0 = all pass, 1 = failure)

#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

#include "TreeNode.h"
#include "ctx.h"
#include "dendro.h"
#include "dvec.h"
#include "ets.h"
#include "ets_msrk.h"
#include "hcurvedata.h"
#include "mesh.h"
#include "meshUtils.h"
#include "mpi.h"
#include "octUtils.h"

using DVec  = ot::DVector<DendroScalar, unsigned int>;
using RhsFn = std::function<double(double, double)>;

struct Pair {
    ts::ETSType type;
    const char* name;
};

static const std::vector<Pair> PAIRS = {
    {ts::ETSType::RK45_CASH_KARP, "RK45_CASH_KARP"},
    {ts::ETSType::RKF45, "RKF45"},
    {ts::ETSType::RK45_DORMAND_PRINCE, "RK45_DORMAND_PRINCE"},
};

static int g_fail = 0;
static int g_rank = 0;

static void check(bool ok, const char* what, const char* name, double v) {
    if (!ok) g_fail++;
    if (!g_rank)
        std::printf("  %-20s %-44s %12.4e  %s\n", name, what, v,
                    ok ? "PASS" : "FAIL");
}

// ---------------------------------------------------------------------------
// ETS level

static const double KRATE = 50.0;

class ErrCtx : public ts::Ctx<ErrCtx, DendroScalar, unsigned int> {
    DVec m_evar;
    std::vector<double> m_amp;

   public:
    explicit ErrCtx(ot::Mesh* pMesh) {
        m_uiMesh = pMesh;
        m_evar.create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                             ot::DVEC_LOC::HOST, 1, true);
        m_uiTinfo = {0.0, 0.0, 0, 0.0, 0.0};
        m_amp.assign(m_uiMesh->getDegOfFreedom(), 0.0);
        const ot::TreeNode* pNodes = &(*(m_uiMesh->getAllElements().begin()));
        const unsigned int eo      = m_uiMesh->getElementOrder();
        const unsigned int* e2n_cg = &(*(m_uiMesh->getE2NMapping().begin()));
        const unsigned int* e2n_dg = &(*(m_uiMesh->getE2NMapping_DG().begin()));
        const unsigned int nPe     = m_uiMesh->getNumNodesPerElement();
        const double maxd          = (double)(1u << m_uiMaxDepth);
        for (unsigned int el = m_uiMesh->getElementLocalBegin();
             el < m_uiMesh->getElementLocalEnd(); el++)
            for (unsigned int q = 0; q < nPe; q++) {
                const unsigned int cg = e2n_cg[el * nPe + q];
                if (cg < m_uiMesh->getNodeLocalBegin() ||
                    cg >= m_uiMesh->getNodeLocalEnd())
                    continue;
                unsigned int own, ix, jy, kz;
                m_uiMesh->dg2eijk(e2n_dg[el * nPe + q], own, ix, jy, kz);
                const double len =
                    (double)(1u << (m_uiMaxDepth - pNodes[own].getLevel()));
                const double x = (pNodes[own].getX() + ix * len / eo) / maxd;
                const double y = (pNodes[own].getY() + jy * len / eo) / maxd;
                const double z = (pNodes[own].getZ() + kz * len / eo) / maxd;
                m_amp[cg]      = 1.0 + 0.5 * std::sin(2.0 * M_PI * x) *
                                           std::sin(2.0 * M_PI * y) *
                                           std::sin(2.0 * M_PI * z);
            }
    }
    ~ErrCtx() { m_evar.destroy_vector(); }

    int initialize() {
        DendroScalar* u = m_evar.get_vec_ptr();
        for (unsigned int n = 0; n < m_evar.get_size(); n++) u[n] = 0.0;
        return 0;
    }

    int rhs(DVec* in, DVec* out, unsigned int, DendroScalar time) {
        const DendroScalar* u = in->get_vec_ptr();
        DendroScalar* v       = out->get_vec_ptr();
        const double ct       = std::cos(time);
        for (unsigned int n = m_uiMesh->getNodeLocalBegin();
             n < m_uiMesh->getNodeLocalEnd(); n++)
            v[n] = -KRATE * (u[n] - m_amp[n] * ct);
        return 0;
    }

    double exact(double a, double t) const {
        const double k2 = KRATE * KRATE;
        return KRATE * a * (KRATE * std::cos(t) + std::sin(t)) / (k2 + 1.0) -
               k2 * a / (k2 + 1.0) * std::exp(-KRATE * t);
    }

    double max_error() {
        const DendroScalar* u = m_evar.get_vec_ptr();
        double e              = 0.0;
        for (unsigned int n = m_uiMesh->getNodeLocalBegin();
             n < m_uiMesh->getNodeLocalEnd(); n++)
            e = std::max(e,
                         std::fabs(u[n] - exact(m_amp[n], m_uiTinfo._m_uiT)));
        double eg = e;
        MPI_Allreduce(&e, &eg, 1, MPI_DOUBLE, MPI_MAX,
                      m_uiMesh->getMPIGlobalCommunicator());
        return eg;
    }

    std::vector<double> owned_state() {
        const DendroScalar* u = m_evar.get_vec_ptr();
        return std::vector<double>(u + m_uiMesh->getNodeLocalBegin(),
                                   u + m_uiMesh->getNodeLocalEnd());
    }

    DVec& get_evolution_vars() { return m_evar; }
    DVec& get_constraint_vars() { return m_evar; }
    DVec& get_primitive_vars() { return m_evar; }

    int pre_stage(DVec) { return 0; }
    int post_stage(DVec) { return 0; }
    int pre_timestep(DVec) { return 0; }
    int post_timestep(DVec) { return 0; }
    int rhs_blk(const DendroScalar*, DendroScalar*, unsigned int, unsigned int,
                DendroScalar) {
        return 0;
    }
    int pre_stage_blk(DendroScalar*, unsigned int, unsigned int, DendroScalar) {
        return 0;
    }
    int post_stage_blk(DendroScalar*, unsigned int, unsigned int,
                       DendroScalar) {
        return 0;
    }
    int pre_timestep_blk(DendroScalar*, unsigned int, unsigned int,
                         DendroScalar) {
        return 0;
    }
    int post_timestep_blk(DendroScalar*, unsigned int, unsigned int,
                          DendroScalar) {
        return 0;
    }
    bool is_remesh() { return false; }
    int write_vtu() { return 0; }
    int write_checkpt() { return 0; }
    int restore_checkpt() { return 0; }
    int finalize() { return 0; }
    int terminal_output() { return 0; }
    unsigned int get_async_batch_sz() { return 1; }
    unsigned int get_num_refine_vars() { return 1; }
    const unsigned int* get_refine_var_ids() {
        static unsigned int ids[1] = {0};
        return ids;
    }
    int grid_transfer(const ot::Mesh*) { return 0; }
    static unsigned int getBlkTimestepFac(unsigned int, unsigned int,
                                          unsigned int) {
        return 1;
    }
};

using ETSE = ts::ETS<DendroScalar, ErrCtx>;

// ---------------------------------------------------------------------------
// tableau level

struct Step {
    double y, yhat, est;
};

static Step rk_step(const RhsFn& f, unsigned int ns, const DendroScalar* c,
                    const DendroScalar* b, const DendroScalar* a,
                    const DendroScalar* e, double t, double y, double h) {
    std::vector<double> k(ns);
    for (unsigned int i = 0; i < ns; i++) {
        double acc = 0.0;
        for (unsigned int j = 0; j < i; j++) acc += a[i * ns + j] * k[j];
        k[i] = f(t + c[i] * h, y + h * acc);
    }
    double sb = 0.0, se = 0.0;
    for (unsigned int i = 0; i < ns; i++) {
        sb += b[i] * k[i];
        se += e[i] * k[i];
    }
    return {y + h * sb, y + h * (sb - se), h * se};
}

struct Ode {
    const char* name;
    RhsFn f;
    std::function<double(double)> exact;
};

static const std::vector<Ode> ODES = {
    {"y'=-y^2", [](double, double y) { return -y * y; },
     [](double t) { return 1.0 / (1.0 + t); }},
    {"y'=cos(t) y", [](double t, double y) { return std::cos(t) * y; },
     [](double t) { return std::exp(std::sin(t)); }},
};

// global-order problems: asymptotic over N = 20..80 for all three pairs
static const Ode ODE_B = {"y'=-2t y^2",
                          [](double t, double y) { return -2.0 * t * y * y; },
                          [](double t) { return 1.0 / (1.0 + t * t); }};
static const Ode ODE_BHAT = {
    "y'=y(1-y)", [](double, double y) { return y * (1.0 - y); },
    [](double t) { return 1.0 / (1.0 + std::exp(-t)); }};

// global error at T=1 for N steps, using b (use_hat=false) or bhat
static double global_err(const Ode& o, unsigned int ns, const DendroScalar* c,
                         const DendroScalar* b, const DendroScalar* a,
                         const DendroScalar* e, int N, bool use_hat) {
    double y = o.exact(0.0), t = 0.0;
    const double h = 1.0 / N;
    for (int s = 0; s < N; s++) {
        const Step st = rk_step(o.f, ns, c, b, a, e, t, y, h);
        y             = use_hat ? st.yhat : st.y;
        t += h;
    }
    return std::fabs(y - o.exact(1.0));
}

static double rate(double e0, double e1) { return std::log2(e0 / e1); }

// least-squares slope of log2(err) vs log2(N) over N = 20, 40, 80
static double global_order(const Ode& o, unsigned int ns, const DendroScalar* c,
                           const DendroScalar* b, const DendroScalar* a,
                           const DendroScalar* e, bool use_hat) {
    const double e0 = std::log2(global_err(o, ns, c, b, a, e, 20, use_hat));
    const double e2 = std::log2(global_err(o, ns, c, b, a, e, 80, use_hat));
    return (e0 - e2) / 2.0;
}

// exposes the tableau that ETS::set_ets_coefficients loads
class TableauETS : public ts::ETS<DendroScalar, ErrCtx> {
   public:
    using ts::ETS<DendroScalar, ErrCtx>::ETS;
    unsigned int ns() const { return m_uiNumStages; }
    const DendroScalar* c() const { return m_uiCi; }
    const DendroScalar* b() const { return m_uiBi; }
    const DendroScalar* a() const { return m_uiAij; }
    const DendroScalar* e() const { return m_uiEi; }
};

static void test_tableaux(ErrCtx& ctx) {
    if (!g_rank) std::printf("\n=== embedded pairs: tableau level ===\n");
    for (const Pair& p : PAIRS) {
        TableauETS tab(&ctx);
        if (tab.set_ets_coefficients(p.type) != 0 || tab.e() == nullptr) {
            check(false, "embedded error weights", p.name, 0.0);
            continue;
        }
        const unsigned int ns = tab.ns();
        const DendroScalar *c = tab.c(), *b = tab.b(), *a = tab.a(),
                           *e = tab.e();

        double se = 0.0, sb = 0.0;
        for (unsigned int i = 0; i < ns; i++) {
            se += e[i];
            sb += b[i];
        }
        check(std::fabs(se) < 1e-15, "sum(e) == 0", p.name, se);
        check(std::fabs(sb - 1.0) < 1e-15, "sum(b) == 1", p.name, sb - 1.0);

        char w[96];
        const double o5 = global_order(ODE_B, ns, c, b, a, e, false);
        std::snprintf(w, sizeof(w), "order(b) ~ 5    [%s]", ODE_B.name);
        check(o5 > 4.6 && o5 < 5.6, w, p.name, o5);
        const double o4 = global_order(ODE_BHAT, ns, c, b, a, e, true);
        std::snprintf(w, sizeof(w), "order(bhat) ~ 4 [%s]", ODE_BHAT.name);
        check(o4 > 3.6 && o4 < 4.4, w, p.name, o4);

        for (const Ode& o : ODES) {
            // one step from the exact state at t0
            const double t0 = 0.3, y0 = o.exact(t0);
            const Step s1  = rk_step(o.f, ns, c, b, a, e, t0, y0, 0.02);
            const Step s2  = rk_step(o.f, ns, c, b, a, e, t0, y0, 0.01);
            const double r = rate(std::fabs(s1.est), std::fabs(s2.est));
            std::snprintf(w, sizeof(w), "estimate ~ h^5  [%s]", o.name);
            check(r > 4.6 && r < 5.4, w, p.name, r);
            const double lte = o.exact(t0 + 0.01) - s2.yhat;
            const double rel = std::fabs(s2.est / lte - 1.0);
            std::snprintf(w, sizeof(w), "estimate/LTE(bhat) - 1 [%s]", o.name);
            check(rel < 0.05, w, p.name, rel);
        }
    }
}

// fixed-dt run of n steps from u=0; returns the owned state
static std::vector<double> run_fixed(ot::Mesh* mesh, ts::ETSType type,
                                     ts::ETSErrCtrl mode, double dt, int n,
                                     std::vector<double>* errs = nullptr) {
    ErrCtx ctx(mesh);
    ctx.set_ts_info({0.0, 1.0, 0, 0.0, dt});
    ETSE ets(&ctx);
    ets.set_ets_coefficients(type);
    ets.init();
    ets.set_error_control(mode, 1e-6, 1e-6);
    for (int i = 0; i < n; i++) {
        ets.evolve();
        if (errs) errs->push_back(ets.get_err_stats().err);
    }
    return ctx.owned_state();
}

struct AdaptRun {
    double err_final = 0.0, err_max_accepted = 0.0, t = 0.0, t_sum = 0.0;
    double dt_first = 0.0, dt_min = 1e300, dt_last = 0.0, dt_max_taken = 0.0;
    unsigned long long rejects = 0, forced = 0;
    long steps = 0, calls = 0;
    unsigned int first_rejects = 0;
};

static AdaptRun run_adaptive(ot::Mesh* mesh, ts::ETSType type, double tol,
                             double dt0, double dtmax, double T) {
    ErrCtx ctx(mesh);
    ctx.set_ts_info({0.0, T, 0, 0.0, dt0});
    ETSE ets(&ctx);
    ets.set_ets_coefficients(type);
    ets.init();
    ets.set_error_control(ts::ETSErrCtrl::ADAPTIVE, tol, tol, 0.9, 0.2, 5.0,
                          10);
    ets.set_dt_max(dtmax);
    AdaptRun r;
    while (ets.curr_time() < T) {
        ets.evolve();
        const ts::ETSErrStats& s = ets.get_err_stats();
        if (r.calls == 0) {
            r.dt_first      = s.dt_taken;
            r.first_rejects = s.rejects;
        }
        r.calls++;
        r.t_sum += s.dt_taken;
        r.err_max_accepted = std::max(r.err_max_accepted, s.err);
        r.dt_min           = std::min(r.dt_min, s.dt_taken);
        r.dt_max_taken     = std::max(r.dt_max_taken, s.dt_taken);
        r.dt_last          = s.dt_taken;
    }
    const ts::ETSErrStats& s = ets.get_err_stats();
    r.rejects                = s.total_rejects;
    r.forced                 = s.forced_accepts;
    r.steps                  = ets.curr_step();
    r.t                      = ets.curr_time();
    r.err_final              = ctx.max_error();
    return r;
}

static void test_ets(ot::Mesh* mesh) {
    if (!g_rank)
        std::printf(
            "\n=== ETS: du/dt = -%.0f (u - a(x) cos t) on the mesh ===\n",
            KRATE);

    for (const Pair& p : PAIRS) {
        // ESTIMATE must not touch the solution
        std::vector<double> eo;
        const std::vector<double> u_off =
            run_fixed(mesh, p.type, ts::ETSErrCtrl::OFF, 0.01, 5);
        const std::vector<double> u_est =
            run_fixed(mesh, p.type, ts::ETSErrCtrl::ESTIMATE, 0.01, 5, &eo);
        int same = (u_off == u_est) ? 1 : 0, same_g = 0;
        MPI_Allreduce(&same, &same_g, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        check(same_g == 1, "ESTIMATE bit-identical to OFF", p.name, 0.0);

        // first-step norm vs dt: inside the transient, local error ~ dt^5
        std::vector<double> e1, e2;
        run_fixed(mesh, p.type, ts::ETSErrCtrl::ESTIMATE, 0.004, 1, &e1);
        run_fixed(mesh, p.type, ts::ETSErrCtrl::ESTIMATE, 0.002, 1, &e2);
        const double r = rate(e1[0], e2[0]);
        check(r > 4.5 && r < 5.5, "ETS error norm ~ dt^5", p.name, r);

        // ADAPTIVE over the transient, three tolerances
        const double T = 2.0;
        double prev    = 0.0;
        long prev_n    = 1L << 40;
        for (double tol : {1e-4, 1e-6, 1e-8}) {
            const AdaptRun a = run_adaptive(mesh, p.type, tol, 0.1, 0.5, T);
            char w[96];
            if (!g_rank)
                std::printf(
                    "  %-20s tol=%.0e steps=%ld rejects=%llu dt[first,min,"
                    "last]=[%.2e,%.2e,%.2e] err=%.2e\n",
                    p.name, tol, a.steps, a.rejects, a.dt_first, a.dt_min,
                    a.dt_last, a.err_final);
            std::snprintf(w, sizeof(w), "tol=%.0e: global err / tol < 10", tol);
            check(a.err_final < 10.0 * tol, w, p.name, a.err_final / tol);
            std::snprintf(w, sizeof(w), "tol=%.0e: accepted norms <= 1", tol);
            check(a.err_max_accepted <= 1.0 && a.forced == 0, w, p.name,
                  a.err_max_accepted);
            std::snprintf(w, sizeof(w), "tol=%.0e: oversized dt0 rejected",
                          tol);
            check(a.first_rejects > 0 && a.dt_first < 0.1, w, p.name,
                  a.dt_first);
            std::snprintf(w, sizeof(w), "tol=%.0e: dt grows after transient",
                          tol);
            check(a.dt_max_taken > 5.0 * a.dt_min, w, p.name,
                  a.dt_max_taken / a.dt_min);
            std::snprintf(w, sizeof(w), "tol=%.0e: step = accepted steps", tol);
            check(a.steps == a.calls, w, p.name, (double)a.steps);
            std::snprintf(w, sizeof(w), "tol=%.0e: t == sum(dt taken)", tol);
            check(a.t == a.t_sum && a.t >= T, w, p.name, a.t - a.t_sum);
            if (prev > 0.0) {
                std::snprintf(w, sizeof(w),
                              "tol=%.0e: tighter tol, smaller err", tol);
                check(a.err_final < prev, w, p.name, a.err_final / prev);
                std::snprintf(w, sizeof(w), "tol=%.0e: tighter tol, more steps",
                              tol);
                check(a.steps > prev_n, w, p.name, (double)a.steps);
            }
            prev   = a.err_final;
            prev_n = a.steps;
        }

        // the cap holds and is reached once the transient has decayed
        const double cap = 0.01;
        const AdaptRun c = run_adaptive(mesh, p.type, 1e-6, 0.1, cap, 1.0);
        check(c.dt_max_taken <= cap && c.dt_last == cap,
              "dt cap holds and is reached", p.name, c.dt_max_taken);
    }

    // refusals
    {
        ErrCtx ctx(mesh);
        ctx.set_ts_info({0.0, 1.0, 0, 0.0, 0.01});
        ETSE rk4(&ctx);
        rk4.set_ets_coefficients(ts::ETSType::RK4);
        rk4.init();
        check(rk4.set_error_control(ts::ETSErrCtrl::ESTIMATE, 1e-6, 1e-6) != 0,
              "refused: RK4 has no embedded pair", "RK4", 0.0);
        check(rk4.set_error_control(ts::ETSErrCtrl::OFF, 0.0, 0.0) == 0,
              "OFF always accepted", "RK4", 0.0);

        ETSE ck(&ctx);
        ck.set_ets_coefficients(ts::ETSType::RKF45);
        ck.init();
        check(ck.set_error_control(ts::ETSErrCtrl::ADAPTIVE, 0.0, 1e-6) != 0,
              "refused: atol = 0", "RKF45", 0.0);
        check(ck.set_error_control(ts::ETSErrCtrl::ADAPTIVE, 1e-6, 1e-6, 1.5) !=
                  0,
              "refused: safety > 1", "RKF45", 0.0);
        check(ck.set_error_control(ts::ETSErrCtrl::ADAPTIVE, 1e-6, 1e-6) == 0,
              "accepted: RKF45 adaptive", "RKF45", 0.0);
        check(ck.set_ets_coefficients(ts::ETSType::RK4) != 0,
              "refused: switch to RK4 while on", "RKF45", 0.0);

        ts::ETS_MSRK<DendroScalar, ErrCtx> ms(&ctx, ts::ETSType::RK4_MSRK3);
        ms.set_ets_coefficients(ts::ETSType::RKF45);
        check(ms.set_error_control(ts::ETSErrCtrl::ESTIMATE, 1e-6, 1e-6) != 0,
              "refused: MSRK", "RK4_MSRK3", 0.0);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &g_rank);

    m_uiMaxDepth              = (argc > 1) ? atoi(argv[1]) : 4;
    const unsigned int eOrder = (argc > 2) ? atoi(argv[2]) : 4;
    _InitializeHcurve(m_uiDim);

    std::function<double(double, double, double)> fr = [](double x, double y,
                                                          double z) {
        const double L = (double)(1u << m_uiMaxDepth);
        return std::sin(2.0 * M_PI * x / L) * std::sin(2.0 * M_PI * y / L) *
               std::sin(2.0 * M_PI * z / L);
    };
    std::vector<ot::TreeNode> tmp;
    function2Octree(fr, tmp, m_uiMaxDepth, 1e-3, eOrder, comm);
    ot::Mesh* mesh = ot::createMesh(tmp.data(), tmp.size(), eOrder, comm, 1,
                                    ot::SM_TYPE::FDM, DENDRO_DEFAULT_GRAIN_SZ,
                                    0.3, DENDRO_DEFAULT_SF_K);
    mesh->setDomainBounds(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0));

    {
        ErrCtx ctx(mesh);
        test_tableaux(ctx);
    }
    test_ets(mesh);

    if (!g_rank)
        std::printf("\n%s (%d check%s failed)\n\n",
                    g_fail == 0 ? "ALL CHECKS PASS" : "FAILURES DETECTED",
                    g_fail, g_fail == 1 ? "" : "s");

    delete mesh;
    MPI_Finalize();
    return g_fail == 0 ? 0 : 1;
}
