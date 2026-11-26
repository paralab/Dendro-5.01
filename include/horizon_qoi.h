#pragma once

namespace dendro_aeh {

/**
 * @brief Read-only QoI snapshot for one horizon from the latest BHaHAHA solve.
 *
 * Queryable by consumers (AMR refinement, BHHistory) and separate from the
 * finder's internal re-seeding history. Radii are coordinate radii w.r.t. the
 * centroid; spin uses -10.0 for "unavailable" and valid==false means "no
 * measurement yet". Lives in its own dependency-free header on purpose.
 */
struct HorizonQoI {
    double r_min = -1.0, r_mean = -1.0, r_max = -1.0, area = -1.0;
    double x_centroid = 0.0, y_centroid = 0.0, z_centroid = 0.0;
    double spin_x = -10.0, spin_y = -10.0, spin_z = -10.0;
    double t   = -1.0;   // simulation time of this measurement
    bool valid = false;  // false until the first successful find

    // flat double layout for MPI / checkpoint (bool -> 1.0/0.0)
    static constexpr int NUM_PACKED = 12;
    void pack(double* o) const {
        o[0]  = r_min;
        o[1]  = r_mean;
        o[2]  = r_max;
        o[3]  = area;
        o[4]  = x_centroid;
        o[5]  = y_centroid;
        o[6]  = z_centroid;
        o[7]  = spin_x;
        o[8]  = spin_y;
        o[9]  = spin_z;
        o[10] = t;
        o[11] = valid ? 1.0 : 0.0;
    }
    void unpack(const double* i) {
        r_min      = i[0];
        r_mean     = i[1];
        r_max      = i[2];
        area       = i[3];
        x_centroid = i[4];
        y_centroid = i[5];
        z_centroid = i[6];
        spin_x     = i[7];
        spin_y     = i[8];
        spin_z     = i[9];
        t          = i[10];
        valid      = (i[11] != 0.0);
    }
};

}  // namespace dendro_aeh
