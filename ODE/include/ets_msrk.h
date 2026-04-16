/**
 * @file ets_msrk.h
 * @brief Multistep Runge-Kutta (MSRK) explicit time stepper for Dendro.
 *
 * Implements the multistep RK4 methods from:
 *   L.T. Sanches, S.R. Brandt, J. Kalinani, L. Ji, E. Schnetter,
 *   "Accelerating Numerical Relativity Simulations with New Multistep
 *    Fourth-Order Runge-Kutta Methods", Class. Quantum Grav. (2026).
 *   arXiv:2603.05763
 *
 * Standard RK4 requires 4 RHS evaluations per time step.  MSRK methods
 * reuse RHS evaluations from previous time steps, so fewer new evaluations
 * are needed while maintaining 4th-order accuracy:
 *
 *   RK4-2(1), RK4-2(2): 2-step methods — 3 fresh RHS evals per step (25%
 *       savings).  Reuse 1 evaluation from the previous step.
 *   RK4-3:  3-step method — 2 fresh RHS evals per step (50% savings).
 *       Reuse 2 evaluations from the two preceding steps.
 *
 * Because each RHS evaluation involves ghost exchange, unzip, compute, and
 * zip, reducing the evaluation count also reduces the number of
 * synchronization points per step.
 *
 * Requires constant time step size (UTS or UTS_ADAP mode).  Falls back to
 * standard RK4 for the first few "bootstrap" steps and after remeshing,
 * following the same strategy described in the paper.
 *
 * @version 0.1
 * @date 2026-04-16
 */

#pragma once
#include "ets.h"

namespace ts {

/**
 * @brief Explicit time stepper using Multistep Runge-Kutta (MSRK) methods.
 *
 * Inherits from ETS and overrides evolve() and sync_with_mesh() to
 * support state-reuse across time steps.  The Ctx interface (rhs,
 * pre_stage, post_stage, etc.) is unchanged — applications only need to
 * replace ETS with ETS_MSRK at construction.
 *
 * @tparam T   Scalar type (typically DendroScalar / double).
 * @tparam Ctx Application context type (CRTP leaf of ts::Ctx).
 */
template <typename T, typename Ctx>
class ETS_MSRK : public ETS<T, Ctx> {
    /*---------------------------------------------------------------
     * Bring parent members into scope (same pattern as ExplicitNUTS).
     *-------------------------------------------------------------*/
    using ETS<T, Ctx>::m_uiAppCtx;
    using ETS<T, Ctx>::m_uiType;
    using ETS<T, Ctx>::m_uiNumStages;
    using ETS<T, Ctx>::m_uiEVar;
    using ETS<T, Ctx>::m_uiStVec;
    using ETS<T, Ctx>::m_uiEVecTmp;
    using ETS<T, Ctx>::m_uiTimeInfo;
    using ETS<T, Ctx>::m_uiAij;
    using ETS<T, Ctx>::m_uiBi;
    using ETS<T, Ctx>::m_uiCi;
    using ETS<T, Ctx>::m_uiIsInternalAlloc;

   protected:
    /**@brief: MSRK variant being used. */
    ETSType m_uiMSRKType;

    /**@brief: number of history slots (1 for 2-step methods, 2 for 3-step). */
    unsigned int m_uiNumHistorySlots;

    /**@brief: index of the first stage that requires a fresh RHS evaluation.
     * Stages [0, m_uiFirstFreshStage) are filled from history. */
    unsigned int m_uiFirstFreshStage;

    /**@brief: history stage vectors — stored RHS evaluations from previous
     * time steps, in zipped form (same layout as m_uiStVec entries). */
    std::vector<DVec> m_uiHistVec;

    /**@brief: bootstrap countdown.  When > 0, evolve() uses standard RK4 to
     * build up the required history.  Decremented each step. */
    unsigned int m_uiBootstrapRemaining;

    /**@brief: time step size from the previous step, used to detect dt
     * changes that would invalidate the history. */
    DendroScalar m_uiPrevDt;

    /**@brief: true when history vectors have been allocated. */
    bool m_uiIsHistoryAlloc = false;

    /*---------------------------------------------------------------
     * Owned coefficient storage.
     *
     * The parent class stores raw pointers (m_uiAij, m_uiBi, m_uiCi)
     * that we point at the arrays below.  During bootstrap we swap to
     * the standard RK4 arrays, then swap back to the MSRK arrays.
     *-------------------------------------------------------------*/

    /**@brief: MSRK Butcher-like tableau (4x4 row-major).
     * Rows for history stages are all zeros (unused).
     * See arXiv:2603.05763, Eq. 9–18, Table 1. */
    DendroScalar m_uiMSRK_Aij[16];

    /**@brief: MSRK final weights b_i.
     * y_{n+1} = y_n + h * sum(b_i * k_i). */
    DendroScalar m_uiMSRK_Bi[4];

    /**@brief: MSRK time coefficients c_i.
     * Stage i evaluates at t_n + c_i * h.  History stage c values are
     * set to 0 but never used (those stages skip the RHS call). */
    DendroScalar m_uiMSRK_Ci[4];

    /**@brief: standard RK4 stage weight matrix for bootstrap steps. */
    static constexpr DendroScalar RK4_STD_Aij[16] = {
        0.0, 0.0,       0.0, 0.0,
        0.5, 0.0,       0.0, 0.0,
        0.0, 0.5,       0.0, 0.0,
        0.0, 0.0,       1.0, 0.0};

    /**@brief: standard RK4 final weights for bootstrap steps. */
    static constexpr DendroScalar RK4_STD_Bi[4] = {
        1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};

    /**@brief: standard RK4 time coefficients for bootstrap steps. */
    static constexpr DendroScalar RK4_STD_Ci[4] = {0.0, 0.5, 0.5, 1.0};

    /*---------------------------------------------------------------
     * Internal helpers.
     *-------------------------------------------------------------*/

    /**@brief: allocate history DVec storage. */
    int allocate_history_vars();

    /**@brief: deallocate history DVec storage. */
    int deallocate_history_vars();

    /**@brief: invalidate history and reset the bootstrap counter. */
    void invalidate_history();

    /**@brief: perform a standard RK4 step and extract the f(t_n, y_n)
     * evaluation into the history buffer for future MSRK steps. */
    void evolve_bootstrap();

    /**@brief: perform a full MSRK step using stored history. */
    void evolve_msrk();

    /**@brief: after an MSRK step, rotate the history vectors so that the
     * most recently computed f(t_n, y_n) is saved for the next step. */
    void rotate_history();

    /**@brief: populate the MSRK coefficient arrays for the given variant.
     * Coefficients are from arXiv:2603.05763 Table 1. */
    void set_msrk_coefficients(ETSType type);

    /**@brief: point the parent coefficient pointers at the MSRK arrays. */
    void use_msrk_coefficients();

    /**@brief: point the parent coefficient pointers at the standard RK4
     * arrays (used during bootstrap). */
    void use_rk4_coefficients();

   public:
    /**
     * @brief Construct a new ETS_MSRK time stepper.
     *
     * @param appCtx   Application context (same as ETS).
     * @param msrkType Which MSRK variant to use.  Must be one of
     *                 RK4_MSRK2_1, RK4_MSRK2_2, or RK4_MSRK3.
     */
    ETS_MSRK(Ctx* appCtx, ETSType msrkType);

    /**@brief: destructor. */
    ~ETS_MSRK();

    /**@brief: initialize the time stepper (calls parent init, then sets
     * up history storage and bootstrap state). */
    void init();

    /**@brief: advance one time step.  Dispatches to evolve_bootstrap()
     * or evolve_msrk() depending on whether history is available. */
    void evolve();

    /**@brief: synchronize internal storage with a new mesh after
     * remeshing.  Invalidates history (triggers re-bootstrap). */
    int sync_with_mesh();
};

/*===================================================================
 * Static member definitions.
 *=================================================================*/

template <typename T, typename Ctx>
constexpr DendroScalar ETS_MSRK<T, Ctx>::RK4_STD_Aij[16];

template <typename T, typename Ctx>
constexpr DendroScalar ETS_MSRK<T, Ctx>::RK4_STD_Bi[4];

template <typename T, typename Ctx>
constexpr DendroScalar ETS_MSRK<T, Ctx>::RK4_STD_Ci[4];

/*===================================================================
 * Constructor / Destructor.
 *=================================================================*/

template <typename T, typename Ctx>
ETS_MSRK<T, Ctx>::ETS_MSRK(Ctx* appCtx, ETSType msrkType)
    : ETS<T, Ctx>(appCtx) {
    m_uiMSRKType = msrkType;

    // All MSRK variants use 4 logical stages.
    m_uiNumStages = 4;

    switch (msrkType) {
        case ETSType::RK4_MSRK2_1:
        case ETSType::RK4_MSRK2_2:
            m_uiNumHistorySlots = 1;
            m_uiFirstFreshStage = 1;
            break;
        case ETSType::RK4_MSRK3:
            m_uiNumHistorySlots = 2;
            m_uiFirstFreshStage = 2;
            break;
        default:
            dendro::logger::error(
                dendro::logger::Scope{"ETS_MSRK"},
                "Invalid MSRK type.  Use RK4_MSRK2_1, RK4_MSRK2_2, or "
                "RK4_MSRK3.");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    set_msrk_coefficients(msrkType);
    use_msrk_coefficients();

    m_uiBootstrapRemaining = m_uiNumHistorySlots;
    m_uiPrevDt             = -1.0;

    dendro::logger::info(
        dendro::logger::Scope{"ETS_MSRK"},
        "Multistep RK time stepper created (variant={}, history_slots={}, "
        "fresh_stages={}).  arXiv:2603.05763",
        static_cast<int>(msrkType), m_uiNumHistorySlots,
        m_uiNumStages - m_uiFirstFreshStage);
}

template <typename T, typename Ctx>
ETS_MSRK<T, Ctx>::~ETS_MSRK() {
    deallocate_history_vars();
}

/*===================================================================
 * Coefficient setup.
 *
 * All coefficient values are from arXiv:2603.05763, Table 1.
 * The 4x4 Aij tableau uses the same row-major layout as the parent
 * ETS class: Aij[stage * 4 + p] is the weight of stage p when
 * forming the intermediate state for stage `stage`.
 *
 * History stage rows are all zeros (those stages are pre-filled from
 * the history buffer, not computed via the Aij accumulation).
 *=================================================================*/

template <typename T, typename Ctx>
void ETS_MSRK<T, Ctx>::set_msrk_coefficients(ETSType type) {
    // Zero-initialize all arrays.
    std::fill(m_uiMSRK_Aij, m_uiMSRK_Aij + 16, 0.0);
    std::fill(m_uiMSRK_Bi, m_uiMSRK_Bi + 4, 0.0);
    std::fill(m_uiMSRK_Ci, m_uiMSRK_Ci + 4, 0.0);

    if (type == ETSType::RK4_MSRK2_1) {
        /*-----------------------------------------------------------
         * RK4-2(1): 2-step method, variant 1.
         * arXiv:2603.05763, Eq. 9–13, Table 1.
         *
         * Best imaginary-axis stability intercept among MSRK variants
         * (2.54 vs 2.83 for standard RK4).  Recommended for
         * wave-dominated systems such as BSSNOK.
         *
         * k0 = f(t_{n-1}, y_{n-1})                   [from history]
         * k1 = f(t_n, y_n)                            [fresh]
         * k2 = f(t_n + c2*h, y_n + h*(a20*k0+a21*k1)) [fresh]
         * k3 = f(t_n + c3*h, y_n + h*(a30*k0+a31*k1+a32*k2)) [fresh]
         * y_{n+1} = y_n + h*(b0*k0 + b1*k1 + b2*k2 + b3*k3)
         *---------------------------------------------------------*/

        // Final weights b_i.
        m_uiMSRK_Bi[0] = -643.0 / 1536.0;
        m_uiMSRK_Bi[1] = -4237.0 / 1092.0;
        m_uiMSRK_Bi[2] = 38125.0 / 10752.0;
        m_uiMSRK_Bi[3] = 4375.0 / 2496.0;

        // Time coefficients c_i.  c0 is unused (history stage).
        // c1 = 0 (evaluation at t_n).
        m_uiMSRK_Ci[0] = 0.0;
        m_uiMSRK_Ci[1] = 0.0;
        m_uiMSRK_Ci[2] = 7.0 / 25.0;
        m_uiMSRK_Ci[3] = -13.0 / 25.0;

        // Stage weight matrix a_ij.
        // Row 0: history stage (all zeros).
        // Row 1: f(t_n, y_n) — no previous-stage contributions.
        // Row 2:
        m_uiMSRK_Aij[2 * 4 + 0] = -49.0 / 1250.0;    // a20
        m_uiMSRK_Aij[2 * 4 + 1] = 399.0 / 1250.0;     // a21
        // Row 3:
        m_uiMSRK_Aij[3 * 4 + 0] = 7033.0 / 960000.0;  // a30
        m_uiMSRK_Aij[3 * 4 + 1] = -217633.0 / 210000.0;  // a31
        m_uiMSRK_Aij[3 * 4 + 2] = 5473.0 / 10752.0;   // a32

        dendro::logger::debug(dendro::logger::Scope{"ETS_MSRK"},
                              "Coefficients set for RK4-2(1)");

    } else if (type == ETSType::RK4_MSRK2_2) {
        /*-----------------------------------------------------------
         * RK4-2(2): 2-step method, variant 2.
         * arXiv:2603.05763, Eq. 9–13, Table 1.
         *
         * Same structure as RK4-2(1) but with alternative coefficients
         * that may perform better for certain PDE systems.
         * Imaginary-axis stability intercept: 2.46.
         *---------------------------------------------------------*/

        m_uiMSRK_Bi[0] = -191.0 / 882.0;
        m_uiMSRK_Bi[1] = 48241.0 / 59994.0;
        m_uiMSRK_Bi[2] = 193750.0 / 4351347.0;
        m_uiMSRK_Bi[3] = 100000.0 / 271791.0;

        m_uiMSRK_Ci[0] = 0.0;
        m_uiMSRK_Ci[1] = 0.0;
        m_uiMSRK_Ci[2] = -99.0 / 50.0;
        m_uiMSRK_Ci[3] = 101.0 / 100.0;

        m_uiMSRK_Aij[2 * 4 + 0] = 1309.0 / 15500.0;
        m_uiMSRK_Aij[2 * 4 + 1] = -31999.0 / 15500.0;

        m_uiMSRK_Aij[3 * 4 + 0] = -241289.0 / 5880000.0;
        m_uiMSRK_Aij[3 * 4 + 1] = 22846301.0 / 16170000.0;
        m_uiMSRK_Aij[3 * 4 + 2] = -936169.0 / 2587200.0;

        dendro::logger::debug(dendro::logger::Scope{"ETS_MSRK"},
                              "Coefficients set for RK4-2(2)");

    } else if (type == ETSType::RK4_MSRK3) {
        /*-----------------------------------------------------------
         * RK4-3: 3-step method.
         * arXiv:2603.05763, Eq. 14–18, Table 1.
         *
         * Reuses evaluations from the two preceding time steps,
         * requiring only 2 fresh RHS evaluations per step (50% savings).
         * Imaginary-axis stability intercept is lower (1.31), so a
         * smaller dt may be required compared to standard RK4.
         *
         * k0 = f(t_{n-2}, y_{n-2})                     [from history]
         * k1 = f(t_{n-1}, y_{n-1})                     [from history]
         * k2 = f(t_n, y_n)                              [fresh]
         * k3 = f(t_n + c3*h, y_n + h*(a30*k0+a31*k1+a32*k2)) [fresh]
         * y_{n+1} = y_n + h*(b0*k0 + b1*k1 + b2*k2 + b3*k3)
         *---------------------------------------------------------*/

        m_uiMSRK_Bi[0] = -85.0 / 1416.0;
        m_uiMSRK_Bi[1] = 131.0 / 408.0;
        m_uiMSRK_Bi[2] = -29.0 / 24.0;
        m_uiMSRK_Bi[3] = 15625.0 / 8024.0;

        // c0, c1 unused (history stages).  c2 = 0 (evaluation at t_n).
        m_uiMSRK_Ci[0] = 0.0;
        m_uiMSRK_Ci[1] = 0.0;
        m_uiMSRK_Ci[2] = 0.0;
        m_uiMSRK_Ci[3] = 9.0 / 25.0;

        // Only stage 3 has non-trivial Aij entries.
        m_uiMSRK_Aij[3 * 4 + 0] = 2511.0 / 62500.0;
        m_uiMSRK_Aij[3 * 4 + 1] = -2268.0 / 15625.0;
        m_uiMSRK_Aij[3 * 4 + 2] = 29061.0 / 62500.0;

        dendro::logger::debug(dendro::logger::Scope{"ETS_MSRK"},
                              "Coefficients set for RK4-3");
    }
}

template <typename T, typename Ctx>
void ETS_MSRK<T, Ctx>::use_msrk_coefficients() {
    m_uiAij = m_uiMSRK_Aij;
    m_uiBi  = m_uiMSRK_Bi;
    m_uiCi  = m_uiMSRK_Ci;
}

template <typename T, typename Ctx>
void ETS_MSRK<T, Ctx>::use_rk4_coefficients() {
    m_uiAij = const_cast<DendroScalar*>(RK4_STD_Aij);
    m_uiBi  = const_cast<DendroScalar*>(RK4_STD_Bi);
    m_uiCi  = const_cast<DendroScalar*>(RK4_STD_Ci);
}

/*===================================================================
 * History variable management.
 *=================================================================*/

template <typename T, typename Ctx>
int ETS_MSRK<T, Ctx>::allocate_history_vars() {
    if (m_uiIsHistoryAlloc) return 0;

    m_uiHistVec.resize(m_uiNumHistorySlots);
    for (unsigned int i = 0; i < m_uiNumHistorySlots; i++)
        m_uiHistVec[i].create_vector(m_uiAppCtx->get_mesh(),
                                     m_uiEVar.get_type(), m_uiEVar.get_loc(),
                                     m_uiEVar.get_dof(),
                                     m_uiEVar.is_ghost_allocated());

    m_uiIsHistoryAlloc = true;
    return 0;
}

template <typename T, typename Ctx>
int ETS_MSRK<T, Ctx>::deallocate_history_vars() {
    if (!m_uiIsHistoryAlloc) return 0;

    for (unsigned int i = 0; i < m_uiNumHistorySlots; i++)
        m_uiHistVec[i].destroy_vector();

    m_uiHistVec.clear();
    m_uiIsHistoryAlloc = false;
    return 0;
}

template <typename T, typename Ctx>
void ETS_MSRK<T, Ctx>::invalidate_history() {
    m_uiBootstrapRemaining = m_uiNumHistorySlots;
    m_uiPrevDt             = -1.0;

    dendro::logger::info(
        dendro::logger::Scope{"ETS_MSRK"},
        "History invalidated — will bootstrap with {} standard RK4 step(s)",
        m_uiBootstrapRemaining);
}

/*===================================================================
 * Initialization.
 *=================================================================*/

template <typename T, typename Ctx>
void ETS_MSRK<T, Ctx>::init() {
    dendro::logger::info(
        dendro::logger::Scope{"ETS_MSRK"},
        "Initializing ETS_MSRK (app context, internal vars, history)");

    m_uiAppCtx->initialize();
    m_uiTimeInfo = m_uiAppCtx->get_ts_info();
    this->allocate_internal_vars();
    allocate_history_vars();

    m_uiAppCtx->set_ets_synced(false);
    this->sync_with_mesh();
}

/*===================================================================
 * Mesh synchronization.
 *
 * After remeshing the old history vectors are on a stale mesh.
 * Rather than attempting an intergrid transfer (which could
 * introduce interpolation error into the time integrator), we
 * simply invalidate the history and re-bootstrap with standard RK4.
 * This is the same strategy recommended in arXiv:2603.05763 §4.
 *=================================================================*/

template <typename T, typename Ctx>
int ETS_MSRK<T, Ctx>::sync_with_mesh() {
    if (m_uiAppCtx->is_ets_synced()) return 0;

    dendro::logger::debug(
        dendro::logger::Scope{"ETS_MSRK"},
        "Syncing with mesh (reallocating internal + history variables)");

    m_uiEVar = m_uiAppCtx->get_evolution_vars();

    this->deallocate_internal_vars();
    deallocate_history_vars();

    this->allocate_internal_vars();
    allocate_history_vars();

    invalidate_history();
    m_uiAppCtx->set_ets_synced(true);

    dendro::logger::debug(dendro::logger::Scope{"ETS_MSRK"},
                          "Finished syncing ETS_MSRK with mesh");
    return 0;
}

/*===================================================================
 * Time stepping — evolve() dispatcher.
 *=================================================================*/

template <typename T, typename Ctx>
void ETS_MSRK<T, Ctx>::evolve() {
#ifdef __PROFILE_ETS__
    this->m_uiCtxpt[ETSPROFILE::EVOLVE].start();
#endif

    m_uiTimeInfo       = m_uiAppCtx->get_ts_info();
    const double dt    = m_uiTimeInfo._m_uiTh;
    const double dtEps = 1.0e-12 * std::abs(dt);

    // If the time step size changed, history is no longer valid.
    if (m_uiPrevDt > 0.0 && std::abs(dt - m_uiPrevDt) > dtEps) {
        dendro::logger::info(
            dendro::logger::Scope{"ETS_MSRK"},
            "Time step size changed ({} -> {}), invalidating history",
            m_uiPrevDt, dt);
        invalidate_history();
    }

    if (m_uiBootstrapRemaining > 0) {
        evolve_bootstrap();
    } else {
        evolve_msrk();
    }

    m_uiPrevDt = dt;

#ifdef __PROFILE_ETS__
    this->m_uiCtxpt[ETSPROFILE::EVOLVE].stop();
#endif
}

/*===================================================================
 * Bootstrap step — standard RK4 with history extraction.
 *
 * Runs a full 4-stage classical RK4 step using the standard Butcher
 * tableau.  After the step completes, extracts f(t_n, y_n) (which is
 * stage 0 of standard RK4, since c0 = 0) and stores it in the
 * history buffer for use by future MSRK steps.
 *=================================================================*/

template <typename T, typename Ctx>
void ETS_MSRK<T, Ctx>::evolve_bootstrap() {
    dendro::logger::info(
        dendro::logger::Scope{"ETS_MSRK"},
        "Bootstrap step ({} remaining) — using standard RK4",
        m_uiBootstrapRemaining);

    const ot::Mesh* pMesh  = m_uiAppCtx->get_mesh();
    const double current_t = m_uiTimeInfo._m_uiT;
    double current_t_adv   = current_t;
    const double dt        = m_uiTimeInfo._m_uiTh;

    m_uiAppCtx->pre_timestep(m_uiEVar);

    // Use standard RK4 coefficients for this step.
    use_rk4_coefficients();

    if (pMesh->isActive()) {
        for (unsigned int stage = 0; stage < m_uiNumStages; stage++) {
            dendro::logger::debug(
                dendro::logger::Scope{"ETS_MSRK"},
                "Bootstrap RK4 stage {}/{}", stage + 1, m_uiNumStages);

            m_uiEVecTmp[0].copy_data(m_uiEVar);

            for (unsigned int p = 0; p < stage; p++) {
                const DendroScalar aip = m_uiAij[stage * m_uiNumStages + p];
                if (aip != 0.0)
                    DVec::axpy(pMesh, aip * dt, m_uiStVec[p], m_uiEVecTmp[0]);
            }

            m_uiAppCtx->post_timestep(m_uiEVecTmp[0]);

            current_t_adv = current_t + m_uiCi[stage] * dt;
            m_uiAppCtx->pre_stage(m_uiStVec[stage]);
            m_uiAppCtx->rhs(&m_uiEVecTmp[0], &m_uiStVec[stage], 1,
                            current_t_adv);
            m_uiAppCtx->post_stage(m_uiStVec[stage]);
        }

        // Final update: y_{n+1} = y_n + sum(b_i * k_i * dt).
        for (unsigned int k = 0; k < m_uiNumStages; k++)
            DVec::axpy(pMesh, m_uiBi[k] * dt, m_uiStVec[k],
                       m_uiEVar);
    }

    m_uiAppCtx->post_timestep(m_uiEVar);

    /*---------------------------------------------------------------
     * Extract history from this bootstrap step.
     *
     * In standard RK4, stage 0 evaluates f(t_n, y_n) (since c0 = 0
     * and the intermediate state is just y_n).  This is exactly the
     * value needed by the MSRK methods as a "reused" evaluation.
     *
     * For the 2-step methods (1 history slot):
     *   After each bootstrap step, save stage 0 as history[0].
     *   Only 1 bootstrap step is needed.
     *
     * For the 3-step method (2 history slots):
     *   We need f(t_{n-2}, y_{n-2}) and f(t_{n-1}, y_{n-1}).
     *   On the first bootstrap step:  save stage 0 → history[0].
     *   On the second bootstrap step: shift history[0] → history[1]
     *     is wrong — we want oldest first.  Instead:
     *   We fill from the back: history is indexed so that
     *   history[0] = oldest, history[1] = most recent.
     *   Bootstrap step 1: save into history[0].
     *   Bootstrap step 2: shift history[0]→[0] stays, save new→[1].
     *   Wait — that's not right either.  Let's think carefully:
     *
     *   m_uiBootstrapRemaining starts at 2.
     *   Step with remaining=2: this is the first step ever.
     *     Save f(t_0, y_0) → history[0].
     *   Step with remaining=1: this is the second step.
     *     Shift: history[0] (=f(t_0,y_0)) is now the oldest.
     *     Save  f(t_1, y_1) → history[1].
     *   After bootstrap, at step 2 using MSRK:
     *     k0 = history[0] = f(t_0, y_0)  ← oldest
     *     k1 = history[1] = f(t_1, y_1)  ← previous step
     *     k2 = f(t_2, y_2)               ← fresh
     *   This is correct per Eq. 14–18.
     *-------------------------------------------------------------*/

    if (m_uiNumHistorySlots == 1) {
        // 2-step methods: just save the current f(t_n, y_n).
        m_uiHistVec[0].copy_data(m_uiStVec[0]);

    } else if (m_uiNumHistorySlots == 2) {
        if (m_uiBootstrapRemaining == 2) {
            // First bootstrap step: save into slot 0 (oldest).
            m_uiHistVec[0].copy_data(m_uiStVec[0]);
        } else {
            // Second bootstrap step: shift old [0] stays, save new to [1].
            // Actually, we need to shift: [0] already has f(t_0, y_0).
            // Now save f(t_1, y_1) into [1].
            m_uiHistVec[1].copy_data(m_uiStVec[0]);
        }
    }

    m_uiBootstrapRemaining--;

    // Restore MSRK coefficients for the next step.
    use_msrk_coefficients();

    m_uiAppCtx->increment_ts_info();
    m_uiTimeInfo = m_uiAppCtx->get_ts_info();
    pMesh->waitAll();

    dendro::logger::debug(dendro::logger::Scope{"ETS_MSRK"},
                          "Bootstrap step finished");
}

/*===================================================================
 * MSRK step — the core multistep Runge-Kutta algorithm.
 *
 * Stages [0, m_uiFirstFreshStage) are filled from the history buffer
 * (no RHS evaluation needed).  Stages [m_uiFirstFreshStage, 4) are
 * computed normally using the MSRK Butcher tableau.
 *
 * The Aij accumulation loop is identical to the parent ETS class,
 * except we skip RHS calls for history stages and avoid zero-weight
 * axpy operations.
 *=================================================================*/

template <typename T, typename Ctx>
void ETS_MSRK<T, Ctx>::evolve_msrk() {
    dendro::logger::debug(
        dendro::logger::Scope{"ETS_MSRK"},
        "MSRK evolve step (step={}, fresh_stages={})",
        m_uiTimeInfo._m_uiStep, m_uiNumStages - m_uiFirstFreshStage);

    const ot::Mesh* pMesh  = m_uiAppCtx->get_mesh();
    const double current_t = m_uiTimeInfo._m_uiT;
    double current_t_adv   = current_t;
    const double dt        = m_uiTimeInfo._m_uiTh;

    m_uiAppCtx->pre_timestep(m_uiEVar);

    if (pMesh->isActive()) {
        // Fill history stages from the stored evaluations.
        for (unsigned int s = 0; s < m_uiFirstFreshStage; s++) {
            m_uiStVec[s].copy_data(m_uiHistVec[s]);
            dendro::logger::debug(dendro::logger::Scope{"ETS_MSRK"},
                                  "Stage {}/{} filled from history", s + 1,
                                  m_uiNumStages);
        }

        // Compute fresh stages.
        for (unsigned int stage = m_uiFirstFreshStage;
             stage < m_uiNumStages; stage++) {
            dendro::logger::debug(dendro::logger::Scope{"ETS_MSRK"},
                                  "Computing fresh stage {}/{}",
                                  stage + 1, m_uiNumStages);

            m_uiEVecTmp[0].copy_data(m_uiEVar);

            // Accumulate contributions from all previous stages (including
            // history stages) via the Aij tableau.  Skip zero coefficients
            // to avoid unnecessary axpy work on large vectors.
            for (unsigned int p = 0; p < stage; p++) {
                const DendroScalar aip =
                    m_uiAij[stage * m_uiNumStages + p];
                if (aip != 0.0)
                    DVec::axpy(pMesh, aip * dt,
                               m_uiStVec[p], m_uiEVecTmp[0]);
            }

            m_uiAppCtx->post_timestep(m_uiEVecTmp[0]);

            current_t_adv = current_t + m_uiCi[stage] * dt;
            m_uiAppCtx->pre_stage(m_uiStVec[stage]);
            m_uiAppCtx->rhs(&m_uiEVecTmp[0], &m_uiStVec[stage], 1,
                            current_t_adv);
            m_uiAppCtx->post_stage(m_uiStVec[stage]);
        }

        // Final update: y_{n+1} = y_n + h * sum(b_i * k_i).
        dendro::logger::debug(dendro::logger::Scope{"ETS_MSRK"},
                              "Computing final update");
        for (unsigned int k = 0; k < m_uiNumStages; k++)
            DVec::axpy(pMesh, m_uiBi[k] * dt,
                       m_uiStVec[k], m_uiEVar);
    }

    m_uiAppCtx->post_timestep(m_uiEVar);

    // Rotate history so the current f(t_n, y_n) is saved for future steps.
    rotate_history();

    m_uiAppCtx->increment_ts_info();
    m_uiTimeInfo = m_uiAppCtx->get_ts_info();
    pMesh->waitAll();

    dendro::logger::debug(dendro::logger::Scope{"ETS_MSRK"},
                          "MSRK evolve step finished");
}

/*===================================================================
 * History rotation.
 *
 * After each MSRK step we need to save the "f(t_n, y_n)" evaluation
 * for use by subsequent steps.  For the 2-step methods this is
 * m_uiStVec[1] (= k1 = f(t_n, y_n)).  For the 3-step method this is
 * m_uiStVec[2] (= k2 = f(t_n, y_n)).
 *
 * In general, the "current time" evaluation is always at stage index
 * m_uiFirstFreshStage (the first fresh stage with c = 0).
 *
 * For 2-step (1 history slot):
 *   history[0] ← m_uiStVec[m_uiFirstFreshStage]
 *
 * For 3-step (2 history slots):
 *   history[0] ← history[1]   (shift older entry)
 *   history[1] ← m_uiStVec[m_uiFirstFreshStage]
 *=================================================================*/

template <typename T, typename Ctx>
void ETS_MSRK<T, Ctx>::rotate_history() {
    if (m_uiNumHistorySlots == 1) {
        m_uiHistVec[0].copy_data(m_uiStVec[m_uiFirstFreshStage]);

    } else if (m_uiNumHistorySlots == 2) {
        // Shift: oldest slot gets the previous "most recent".
        m_uiHistVec[0].copy_data(m_uiHistVec[1]);
        // Save current f(t_n, y_n).
        m_uiHistVec[1].copy_data(m_uiStVec[m_uiFirstFreshStage]);
    }
}

}  // end of namespace ts
