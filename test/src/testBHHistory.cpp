// BHHistory stores only the raw location stream and horizon QoIs: the puncture
// tracker's velocity lives in the solver, and every derived series is rebuilt
// from positions. Catches: extra state in the checkpoint layout, a restore
// that skips the replay, wrong relative position or velocity, wrong angular
// velocity, a missing post-merger freeze. Exit code 0 = pass.

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "b91_util.h"
#include "bh_history.h"

static int fails = 0;

static void check(bool ok, const char* what) {
    std::printf("  %-62s %s\n", what, ok ? "ok" : "FAIL");
    if (!ok) fails++;
}

static bool samePts(const std::vector<Point>& a, const std::vector<Point>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
        if (a[i].x() != b[i].x() || a[i].y() != b[i].y() ||
            a[i].z() != b[i].z())
            return false;
    return true;
}

int main() {
    // a tilted binary (z differs from x, so a z = x typo shows), then a merger
    dendro_bh::BHHistoryConfig cfg;
    cfg.merged_sep_tol = 0.5;
    dendro_bh::BHHistory h(cfg);
    const unsigned N = 12, nMerge = 9;  // samples, first merged sample
    for (unsigned n = 0; n < N; n++) {
        const double t = 0.25 * n, s = (n < nMerge) ? 1.0 : 0.1;
        h.append(Point(s * std::cos(t), s * std::sin(t), 0.3 * s * t),
                 Point(-s * std::cos(t), -s * std::sin(t), -0.2 * s), t);
    }
    dendro_aeh::HorizonQoI q;
    q.r_mean = 0.7;
    q.valid  = true;
    h.set_horizon_qoi(2, q);

    // the checkpoint layout: [B, N, Q], 3*B*N coordinates, N times, Q QoIs;
    // a tracker velocity stored here would add 3*B doubles per step or blob
    const std::string blob        = h.encode();
    const std::vector<double> raw = dendro_aeh::b91_decode<double>(blob);
    const size_t Q = 3, QP = dendro_aeh::HorizonQoI::NUM_PACKED;
    check(raw.size() == 3 + 3 * 2 * N + N + Q * QP && raw[0] == 2 &&
              raw[1] == N && raw[2] == Q,
          "checkpoint holds only positions, times and QoIs");

    // a blob carrying one extra point per body (a velocity) is refused
    std::vector<double> extra = raw;
    extra.insert(extra.end(), 6, 0.0);
    bool refused = false;
    try {
        dendro_bh::BHHistory r(cfg);
        r.restore(dendro_aeh::b91_encode(extra));
    } catch (const std::runtime_error&) {
        refused = true;
    }
    check(refused, "restore refuses a blob with extra per-body state");

    // round trip: raw stream and every derived series bit-identical
    dendro_bh::BHHistory r(cfg);
    r.restore(blob);
    bool same = r.times() == h.times() && r.size() == N;
    for (unsigned b = 0; b < 2; b++)
        same = same && samePts(r.body_locations()[b], h.body_locations()[b]);
    same =
        same &&
        samePts(r.relative_position_history(), h.relative_position_history()) &&
        samePts(r.relative_velocity_history(), h.relative_velocity_history()) &&
        samePts(r.angular_velocity_history(), h.angular_velocity_history()) &&
        r.angle_history() == h.angle_history() &&
        r.clean_wavelength_history() == h.clean_wavelength_history() &&
        r.has_horizon(2) && !r.has_horizon(0) &&
        r.mean_horizon_radius(2) == 0.7;
    check(same, "restore(encode()) rebuilds every series bit-identically");

    // relative position is body0 - body1; relative velocity is its backward
    // difference, so it is computed from stored positions only
    bool kin = h.relative_position_history().size() == N &&
               h.relative_velocity_history().size() == N - 1;
    for (unsigned n = 0; kin && n < N; n++) {
        const Point &a  = h.body_locations()[0][n],
                    &b  = h.body_locations()[1][n];
        const Point& r0 = h.relative_position_history()[n];
        kin             = r0.x() == a.x() - b.x() && r0.y() == a.y() - b.y() &&
                          r0.z() == a.z() - b.z();
        if (kin && n > 0) {
            const Point& rp = h.relative_position_history()[n - 1];
            const Point& v  = h.relative_velocity_history()[n - 1];
            const double dt = h.times()[n] - h.times()[n - 1];
            kin             = v.x() == (r0.x() - rp.x()) / dt &&
                              v.y() == (r0.y() - rp.y()) / dt &&
                              v.z() == (r0.z() - rp.z()) / dt;
        }
    }
    check(kin, "relative position and velocity come from the positions");

    // circular planar orbit, angle w*dt per step: r(n) x (r(n) - r(n-1))
    // / (dt |r|^2) = sin(w dt) / dt exactly; measured 2.8e-16, so 1e-12 is
    // rounding with margin and any formula error is far above it
    {
        dendro_bh::BHHistory c(cfg);
        const double w = 0.8, dt = 0.1, R = 2.0;
        for (unsigned n = 0; n < 6; n++) {
            const double a = w * dt * n;
            c.append(Point(R * std::cos(a), R * std::sin(a), 0.0),
                     Point(-R * std::cos(a), -R * std::sin(a), 0.0), dt * n);
        }
        const double want = std::sin(w * dt) / dt;
        double err        = 0.0;
        for (unsigned n = 1; n < 6; n++) {
            const Point& om = c.angular_velocity_history()[n];
            err = std::fmax(err, std::fabs(om.z() - want) / want +
                                     std::fabs(om.x()) + std::fabs(om.y()));
        }
        std::printf("  circular orbit: max relative error of omega %.2e\n",
                    err);
        check(err < 1e-12, "angular velocity on a circular orbit");
    }

    // after the merger omega is frozen at a pre-merger value, and the clean
    // wavelength never grows; either of the last two pre-merger samples is
    // accepted as the frozen value
    const std::vector<Point>& om = h.angular_velocity_history();
    auto eq                      = [&](unsigned i, unsigned j) {
        return om[i].x() == om[j].x() && om[i].y() == om[j].y() &&
               om[i].z() == om[j].z();
    };
    bool frozen = eq(nMerge, nMerge - 1) || eq(nMerge, nMerge - 2);
    for (unsigned n = nMerge + 1; n < N; n++) frozen = frozen && eq(n, nMerge);
    const std::vector<double>& lam = h.clean_wavelength_history();
    for (unsigned n = 1; n < N; n++) frozen = frozen && lam[n] <= lam[n - 1];
    check(frozen, "omega frozen after merger, wavelength non-increasing");

    // clean_wavelength_at: first value up to t_begin, last past the end, linear
    // between; a quarter step catches swapped weights (1e-14: roundoff floor)
    {
        const std::vector<double>& t = h.times();
        unsigned k                   = 1;
        while (k < N && lam[k] == lam[k - 1]) k++;
        bool ok = k < N && h.clean_wavelength_at(-1.0) == lam.front() &&
                  h.clean_wavelength_at(t.back() + 1.0) == lam.back() &&
                  h.clean_wavelength_at(t[k]) == lam[k];
        if (ok) {
            const double tq   = t[k - 1] + 0.25 * (t[k] - t[k - 1]);
            const double want = 0.75 * lam[k - 1] + 0.25 * lam[k];
            ok                = std::fabs(h.clean_wavelength_at(tq) - want) <=
                                1e-14 * std::fabs(want);
        }
        check(ok, "clean_wavelength_at ends and interpolation");
    }

    // separation: |body1 - body0| at the latest sample; a third body gives
    // current_separation_of a pair the binary series do not cover
    {
        const Point &a    = h.body_locations()[0].back(),
                    &b    = h.body_locations()[1].back();
        const double want = std::sqrt((a.x() - b.x()) * (a.x() - b.x()) +
                                      (a.y() - b.y()) * (a.y() - b.y()) +
                                      (a.z() - b.z()) * (a.z() - b.z()));
        dendro_bh::BHHistoryConfig c3 = cfg;
        c3.num_bodies                 = 3;
        dendro_bh::BHHistory t3(c3);
        t3.append({Point(9, 9, 9), Point(0, 0, 0), Point(9, 9, 9)}, 0.0);
        t3.append({Point(1, 2, 2), Point(0, 0, 0), Point(0, 3, 4)}, 1.0);
        const bool ok =
            h.separation() == want && h.current_separation_of(0, 1) == want &&
            h.current_separation_of(1, 0) == want &&
            h.current_separation_of(0, 2) == 0.0 &&
            t3.current_separation_of(0, 1) == 3.0 &&
            t3.current_separation_of(2, 1) == 5.0 && t3.separation() == 0.0;
        check(ok, "separation and current_separation_of");
    }

    // wrong body counts and malformed blobs throw instead of corrupting state
    {
        int thrown  = 0;
        auto throws = [&](auto f) {
            try {
                f();
            } catch (const std::runtime_error&) {
                thrown++;
            }
        };
        dendro_bh::BHHistory w(cfg);
        throws([&] { w.append({Point(0, 0, 0)}, 0.0); });
        throws([&] {
            w.append({Point(0, 0, 0), Point(1, 0, 0), Point(2, 0, 0)}, 0.0);
        });
        throws([&] { w.restore(""); });
        std::vector<double> cut = raw;  // one double short
        cut.pop_back();
        throws([&] { w.restore(dendro_aeh::b91_encode(cut)); });
        dendro_bh::BHHistoryConfig c3 = cfg;
        c3.num_bodies                 = 3;
        dendro_bh::BHHistory t3(c3);
        t3.append({Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1)}, 0.0);
        throws([&] { w.restore(t3.encode()); });  // 3 bodies into 2
        try {
            h.horizon_qoi(Q);  // past the stored QoIs
        } catch (const std::out_of_range&) {
            thrown++;
        }
        check(thrown == 6 && w.size() == 0,
              "bad counts, malformed blobs, QoI index throw");
    }

    std::printf("testBHHistory: %d check(s) failed\n", fails);
    return fails ? 1 : 0;
}
