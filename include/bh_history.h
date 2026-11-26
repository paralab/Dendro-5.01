#pragma once

#include <string>
#include <vector>

#include "horizon_qoi.h"
#include "point.h"

namespace dendro_bh {

/**
 * @brief Inputs the derived kinematics need beyond the location stream itself.
 *
 * Supplied by the caller so the tracker carries no solver-specific globals.
 */
struct BHHistoryConfig {
    unsigned int num_bodies = 2;    // binary derived series require == 2
    double t_begin          = 0.0;  // start time (for wavelength sampling)
    // seed for angular-velocity entry 0; for a binary use body1 - body0
    Point initial_rel_position{0.0, 0.0, 0.0};
    Point initial_rel_velocity{0.0, 0.0, 0.0};
    double merged_sep_tol = 1.0;  // |r| below which omega is frozen
};

/**
 * @brief Physics-agnostic tracked-body history container.
 *
 * Owns the raw per-body location/time stream and incrementally maintains the
 * binary derived kinematics (relative position/velocity, angular velocity
 * frozen post-merger, inter-body angle, monotone "clean" orbital wavelength),
 * plus the latest per-horizon QoI snapshot from an AH finder. Maintaining the
 * series incrementally avoids recomputing it on every query.
 */
class BHHistory {
   public:
    explicit BHHistory(const BHHistoryConfig& cfg);

    // ---- ingestion (call once per evolution step) --------------------------
    void append(const std::vector<Point>& body_locations, double t);
    void append(const Point& body0, const Point& body1, double t);  // binary

    // ---- raw history -------------------------------------------------------
    const std::vector<std::vector<Point>>& body_locations() const {
        return body_track_;
    }
    const std::vector<double>& times() const { return t_; }
    size_t size() const { return t_.size(); }
    unsigned int num_bodies() const { return cfg_.num_bodies; }

    // ---- derived series (binary only; empty when num_bodies != 2) ----------
    const std::vector<double>& angle_history() const { return angle_; }
    const std::vector<Point>& relative_position_history() const {
        return relpos_;
    }
    const std::vector<Point>& relative_velocity_history() const {
        return relvel_;
    }
    const std::vector<Point>& angular_velocity_history() const {
        return angvel_;
    }
    const std::vector<double>& clean_wavelength_history() const {
        return wavelength_;
    }
    // interpolated clean wavelength at a (retarded) time
    double clean_wavelength_at(double t_ret) const;
    // |body1 - body0| at the latest sample (0 if no samples / not binary)
    double separation() const;
    // |body j - body i| at the latest sample (generic N-body)
    double current_separation_of(unsigned int i, unsigned int j) const;

    // ---- per-horizon QoI snapshots (fed from an AH finder) -----------------
    // idx follows the finder's convention (binary: 0,1 inspiral, 2 common).
    void set_horizon_qoi(unsigned int idx, const dendro_aeh::HorizonQoI& q);
    const dendro_aeh::HorizonQoI& horizon_qoi(unsigned int idx) const {
        return qoi_[idx];
    }
    bool has_horizon(unsigned int idx) const {
        return idx < qoi_.size() && qoi_[idx].valid;
    }
    double mean_horizon_radius(unsigned int idx) const {
        return has_horizon(idx) ? qoi_[idx].r_mean : -1.0;
    }

    // ---- checkpoint (base91; replaces encode_bh_locs/decode_bh_locs) -------
    std::string encode() const;
    void restore(const std::string& blob);

   private:
    // appends the derived-series entry for already-stored sample index n
    void update_derived_(size_t n);

    BHHistoryConfig cfg_;

    std::vector<std::vector<Point>> body_track_;  // [body][step]
    std::vector<double> t_;                       // [step]

    // cached derived (binary)
    std::vector<double> angle_;
    std::vector<Point> relpos_;
    std::vector<Point> relvel_;
    std::vector<Point> angvel_;
    std::vector<double> wavelength_;
    int index_pre_merge_ = 0;

    std::vector<dendro_aeh::HorizonQoI> qoi_;  // grows on demand
};

}  // namespace dendro_bh
