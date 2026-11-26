#include "bh_history.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "b91_util.h"

namespace dendro_bh {

// (r x v) / |r|^2 ; zero vector if r is degenerate
static Point angular_velocity(const Point& r, const Point& v) {
    const double r2 = r.x() * r.x() + r.y() * r.y() + r.z() * r.z();
    if (r2 < std::numeric_limits<double>::epsilon()) {
        return Point(0.0, 0.0, 0.0);
    }
    return Point((r.y() * v.z() - r.z() * v.y()) / r2,
                 (r.z() * v.x() - r.x() * v.z()) / r2,
                 (r.x() * v.y() - r.y() * v.x()) / r2);
}

BHHistory::BHHistory(const BHHistoryConfig& cfg)
    : cfg_(cfg), body_track_(cfg.num_bodies) {}

void BHHistory::append(const std::vector<Point>& body_locations, double t) {
    if (body_locations.size() != cfg_.num_bodies) {
        throw std::runtime_error(
            "BHHistory::append given the wrong number of body locations");
    }
    for (unsigned int b = 0; b < cfg_.num_bodies; ++b) {
        body_track_[b].push_back(body_locations[b]);
    }
    t_.push_back(t);
    update_derived_(t_.size() - 1);
}

void BHHistory::append(const Point& body0, const Point& body1, double t) {
    append(std::vector<Point>{body0, body1}, t);
}

void BHHistory::update_derived_(size_t n) {
    // derived kinematics are only defined for a binary system
    if (cfg_.num_bodies != 2) return;

    const Point& b0 = body_track_[0][n];
    const Point& b1 = body_track_[1][n];

    // relative position (body0 - body1) and orbital-plane angle
    relpos_.emplace_back(b0.x() - b1.x(), b0.y() - b1.y(), b0.z() - b1.z());
    angle_.push_back(std::atan2(b0.y() - b1.y(), b0.x() - b1.x()));

    // relative velocity over interval [n-1, n]
    if (n >= 1) {
        const double dt      = t_[n] - t_[n - 1];
        const Point& rp      = relpos_[n];
        const Point& rp_prev = relpos_[n - 1];
        relvel_.emplace_back((rp.x() - rp_prev.x()) / dt,
                             (rp.y() - rp_prev.y()) / dt,
                             (rp.z() - rp_prev.z()) / dt);
    }

    // angular velocity: entry 0 from the seed, then frozen at the last
    // pre-merger value once separation drops below merged_sep_tol
    if (n == 0) {
        angvel_.push_back(angular_velocity(cfg_.initial_rel_position,
                                           cfg_.initial_rel_velocity));
    } else {
        Point r = relpos_[n];
        const double rabs =
            std::sqrt(r.x() * r.x() + r.y() * r.y() + r.z() * r.z());
        if (rabs > cfg_.merged_sep_tol) {
            angvel_.push_back(angular_velocity(r, relvel_[n - 1]));
            index_pre_merge_ = static_cast<int>(n - 1);
        } else {
            angvel_.push_back(angvel_[index_pre_merge_]);
        }
    }

    // clean wavelength lambda = 2*pi*c/|omega| (c = 1), forced monotone
    const Point& w         = angvel_[n];
    const double omega_mag = std::hypot(w.x(), w.y(), w.z());
    double lam = (omega_mag < std::numeric_limits<double>::epsilon())
                     ? std::numeric_limits<double>::max()
                     : (2.0 * M_PI * 1.0) / omega_mag;
    if (n >= 1) lam = std::min(lam, wavelength_[n - 1]);
    wavelength_.push_back(lam);
}

double BHHistory::clean_wavelength_at(double t_ret) const {
    if (wavelength_.empty()) return std::numeric_limits<double>::max();

    // before the start of the history -> first value
    if (t_ret <= cfg_.t_begin) return wavelength_.front();

    auto it            = std::lower_bound(t_.begin(), t_.end(), t_ret);
    const size_t index = std::distance(t_.begin(), it);

    // after the end of the history -> last value
    if (index == t_.size()) return wavelength_.back();
    // before the first stored time (guard against index-1 underflow)
    if (index == 0) return wavelength_.front();

    const double t0    = t_[index - 1];
    const double t1    = t_[index];
    const double alpha = (t_ret - t0) / (t1 - t0);

    return (1.0 - alpha) * wavelength_[index - 1] + alpha * wavelength_[index];
}

double BHHistory::separation() const {
    if (cfg_.num_bodies != 2 || relpos_.empty()) return 0.0;
    Point r = relpos_.back();
    return std::sqrt(r.x() * r.x() + r.y() * r.y() + r.z() * r.z());
}

double BHHistory::current_separation_of(unsigned int i, unsigned int j) const {
    if (i >= cfg_.num_bodies || j >= cfg_.num_bodies || t_.empty()) return 0.0;
    const Point& a  = body_track_[i].back();
    const Point& b  = body_track_[j].back();
    const double dx = a.x() - b.x(), dy = a.y() - b.y(), dz = a.z() - b.z();
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void BHHistory::set_horizon_qoi(unsigned int idx,
                                const dendro_aeh::HorizonQoI& q) {
    if (idx >= qoi_.size()) qoi_.resize(idx + 1);
    qoi_[idx] = q;
}

// ---- checkpoint ------------------------------------------------------------
// Flat double layout (self-describing header), base91-encoded:
//   [num_bodies, num_steps, num_qoi,
//    body0 coords (3*N), body1 coords (3*N), ... ,
//    times (N),
//    qoi (num_qoi * HorizonQoI::NUM_PACKED)]
// Only the raw stream + QoIs are stored; derived series are recomputed on
// restore (smaller checkpoints, no version skew).
std::string BHHistory::encode() const {
    const size_t N = t_.size();
    const size_t B = cfg_.num_bodies;
    const size_t Q = qoi_.size();

    std::vector<double> buf;
    buf.reserve(3 + B * N * 3 + N + Q * dendro_aeh::HorizonQoI::NUM_PACKED);
    buf.push_back(static_cast<double>(B));
    buf.push_back(static_cast<double>(N));
    buf.push_back(static_cast<double>(Q));

    for (size_t b = 0; b < B; ++b) {
        for (size_t s = 0; s < N; ++s) {
            const Point& p = body_track_[b][s];
            buf.push_back(p.x());
            buf.push_back(p.y());
            buf.push_back(p.z());
        }
    }
    for (size_t s = 0; s < N; ++s) buf.push_back(t_[s]);

    for (size_t q = 0; q < Q; ++q) {
        double packed[dendro_aeh::HorizonQoI::NUM_PACKED];
        qoi_[q].pack(packed);
        for (int k = 0; k < dendro_aeh::HorizonQoI::NUM_PACKED; ++k)
            buf.push_back(packed[k]);
    }

    return dendro_aeh::b91_encode(buf);
}

void BHHistory::restore(const std::string& blob) {
    std::vector<double> buf = dendro_aeh::b91_decode<double>(blob);
    if (buf.size() < 3) {
        throw std::runtime_error("BHHistory::restore got a malformed blob");
    }

    size_t off            = 0;
    const size_t B        = static_cast<size_t>(buf[off++]);
    const size_t N        = static_cast<size_t>(buf[off++]);
    const size_t Q        = static_cast<size_t>(buf[off++]);
    constexpr size_t QP   = dendro_aeh::HorizonQoI::NUM_PACKED;

    const size_t expected = 3 + B * N * 3 + N + Q * QP;
    if (B != cfg_.num_bodies || buf.size() != expected) {
        throw std::runtime_error(
            "BHHistory::restore size/num_bodies mismatch with config");
    }

    // reset everything
    body_track_.assign(B, {});
    t_.clear();
    angle_.clear();
    relpos_.clear();
    relvel_.clear();
    angvel_.clear();
    wavelength_.clear();
    index_pre_merge_ = 0;
    qoi_.clear();

    for (size_t b = 0; b < B; ++b) {
        body_track_[b].reserve(N);
        for (size_t s = 0; s < N; ++s) {
            const double x = buf[off++];
            const double y = buf[off++];
            const double z = buf[off++];
            body_track_[b].emplace_back(x, y, z);
        }
    }
    t_.reserve(N);
    for (size_t s = 0; s < N; ++s) t_.push_back(buf[off++]);

    qoi_.resize(Q);
    for (size_t q = 0; q < Q; ++q) {
        qoi_[q].unpack(&buf[off]);
        off += QP;
    }

    // rebuild derived series by replay
    for (size_t n = 0; n < N; ++n) update_derived_(n);
}

}  // namespace dendro_bh
