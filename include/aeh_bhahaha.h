#pragma once

#include <math.h>

#include <cstdlib>
#include <functional>
#include <stdexcept>
#include <vector>

#include "BHaHAHA.h"
#include "base.h"
#include "daUtils.h"
#include "dendro.h"
#include "json.hpp"
#include "mesh.h"
#include "mpi.h"
#include "point.h"
#include "profiler.h"

namespace dendro_aeh {

template <typename T>
inline void print_vec(const std::vector<T>& vec, const std::string& prefix) {
    std::cout << prefix << ": ";
    for (const auto& val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

template <typename T>
inline std::string b91_encode(const std::vector<T>& vec) {
    return base<91>::encode(std::string(
        reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(T)));
}

template <typename T>
inline std::vector<T> b91_decode(const std::string& input) {
    std::string temp         = base<91>::decode(input);
    const size_t num_entries = temp.size() / sizeof(T);
    std::vector<T> output(num_entries);

    std::memcpy(output.data(), temp.data(), sizeof(T) * num_entries);

    return output;
}

#define PRINT_BAH_VAR(var) \
    std::cout << "    " << #var << " = " << var << std::endl;

inline void print_bha_param_data(bhahaha_params_and_data_struct* bha_param_data,
                                 unsigned int which_horizon) {
    // TODO: check to make sure input_metric_data isn't bad

    std::cout << "PARAMETER DUMP FOR PARAM DATA FOR HORIZON " << which_horizon
              << std::endl;
    PRINT_BAH_VAR(bha_param_data->time_external_input);
    PRINT_BAH_VAR(bha_param_data->iteration_external_input);
    PRINT_BAH_VAR(bha_param_data->Nr_external_input);
    PRINT_BAH_VAR(bha_param_data->r_min_external_input);
    PRINT_BAH_VAR(bha_param_data->dr_external_input);
    PRINT_BAH_VAR(bha_param_data->num_resolutions_multigrid);
    // TODO: ntheta array and nphi array
    PRINT_BAH_VAR(bha_param_data->use_fixed_radius_guess_on_full_sphere);
    PRINT_BAH_VAR(bha_param_data->cfl_factor);
    PRINT_BAH_VAR(bha_param_data->M_scale);
    PRINT_BAH_VAR(bha_param_data->eta_damping_times_M);
    PRINT_BAH_VAR(bha_param_data->KO_strength);
    PRINT_BAH_VAR(bha_param_data->max_iterations);
    PRINT_BAH_VAR(bha_param_data->Theta_Linf_times_M_tolerance);
    PRINT_BAH_VAR(bha_param_data->Theta_L2_times_M_tolerance);
    PRINT_BAH_VAR(bha_param_data->which_horizon);
    PRINT_BAH_VAR(bha_param_data->num_horizons);
    PRINT_BAH_VAR(bha_param_data->verbosity_level);
    PRINT_BAH_VAR(
        bha_param_data->enable_eta_varying_alg_for_precision_common_horizon);
    // TODO: t_m1, t_m2, and rmin, rmax, x, y, and z
}

template <typename T>
inline void restore_vector(std::vector<T>& original,
                           const std::vector<T>& restored,
                           const size_t expected_size) {
    if (restored.size() != expected_size) {
        throw std::runtime_error(
            "ERROR when restoring a vector for the AEH solver! The size is "
            "different from the expected size!");
    }

    // otherwise just copy straight in
    // NOTE: it MUST copy, due to pointer information for prev_horizon
    for (size_t i = 0; i < expected_size; ++i) {
        original[i] = restored[i];
    }
}

static inline double square(double x) { return x * x; }

static inline double dist(double x1, double x2, double y1, double y2, double z1,
                          double z2) {
    return sqrt(square(x1 - x2) + square(y1 - y2) + square(z1 - z2));
}

double inline idx3_spherical(const int ir, const int itheta, const int iphi,
                             const int n_theta, const int n_r) {
    return ir + n_r * (itheta + n_theta * iphi);
}

struct SimpleBlackHoleData {
    double x, y, z;
    double mass;

    SimpleBlackHoleData(double x, double y, double z, double mass)
        : x(x), y(y), z(z), mass(mass) {}
};

template <typename T>
void inline fill_vector_with_defaults(const std::vector<T>& input,
                                      std::vector<T>& output, T default_value,
                                      size_t size_requested) {
    output.resize(size_requested);
    size_t size_fill = std::min(input.size(), size_requested);

    // copy it all out
    std::copy(input.begin(), input.begin() + size_fill, output.begin());
    std::fill(output.begin() + size_fill, output.end(), default_value);
}

class AEH_BHaHAHA {
   private:
    profiler_t bhahaha_profiler_;
    bool initialized_ = false;

    unsigned int num_horizons_;
    bool is_bbh_;

    unsigned int file_output_freq_;

    // this is the most important vector, it's what is fed in
    std::vector<bhahaha_params_and_data_struct> bha_param_data_;
    std::vector<double> x_guess_, y_guess_, z_guess_;
    std::vector<double> r_min_guess_, r_max_guess_;

    std::vector<double> x_center_m1_, y_center_m1_, z_center_m1_;
    std::vector<double> x_center_m2_, y_center_m2_, z_center_m2_;
    std::vector<double> x_center_m3_, y_center_m3_, z_center_m3_;
    std::vector<double> t_m1_, t_m2_, t_m3_;
    std::vector<double> r_min_m1_, r_min_m2_, r_min_m3_;
    std::vector<double> r_max_m1_, r_max_m2_, r_max_m3_;
    std::vector<double> prev_horizon_m1_, prev_horizon_m2_, prev_horizon_m3_;
    std::vector<int> bah_horizon_active_;

    unsigned int num_resolutions_multigrid_;
    std::vector<int> ntheta_array_multigrid_;
    std::vector<int> nphi_array_multigrid_;

    int bah_verbosity_level_;
    int bah_enable_eta_varying_alg_for_precision_common_horizon_;

    std::vector<double> bah_m_scale_;
    std::vector<double> bah_cfl_factor_;
    std::vector<int> bah_max_iterations_;
    std::vector<double> bah_theta_l2_times_m_tolerance_;
    std::vector<double> bah_theta_linf_times_m_tolerance_;
    std::vector<double> bah_eta_damping_times_m_;
    std::vector<double> bah_ko_strength_;
    std::vector<int> bah_nr_interp_max_;
    std::vector<double> bah_max_search_radius_;
    // used for initial guess and filling out the initial data
    std::vector<double> bah_initial_grid_x_center_;
    std::vector<double> bah_initial_grid_y_center_;
    std::vector<double> bah_initial_grid_z_center;
    int max_ntheta_;
    int max_nphi_;
    std::string out_dir_;

    // not a parameter to set, this is used internally
    std::vector<int> bah_use_fixed_radius_guess_on_full_sphere_;

    Point grid_limits_[2];
    Point domain_limits_[2];

    // NOTE: this should be set if moving, but this way they're defined together
    // it wants them in gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
    // kzz
    std::vector<int> variable_indices_;
    std::function<std::vector<double>(const std::vector<double>&)> transform_;

    // internal flag to let us know if the first step has passed
    bool have_done_initial_test_ = false;

    std::vector<bool> failed_last_find_;
    std::vector<int> failed_last_find_int_;

   public:
    // add more input types for data...
    AEH_BHaHAHA(const unsigned int n_horizons, const bool is_binary_black_hole,
                const std::vector<double>& initial_x_center,
                const std::vector<double>& initial_y_center,
                const std::vector<double>& initial_z_center,
                const unsigned int n_resolutions_multigrid,
                const std::vector<double>& m_scale,
                const std::vector<double>& cfl_factor,
                const std::vector<int>& max_iterations,
                const std::vector<double>& theta_l2_m_tol,
                const std::vector<double>& theta_linf_m_tol,
                const std::vector<double>& eta_damp_m,
                const std::vector<double>& ko_strength,
                const std::vector<double>& max_search_radius,
                const std::vector<int>& nr_interp_max, const int ntheta_max,
                const int nphi_max, const std::string save_directory,
                const std::vector<SimpleBlackHoleData>& blackholes,
                const std::vector<int>& indices_extract,
                std::function<std::vector<double>(const std::vector<double>&)>
                    transform,
                const Point grid_limits[2], const Point domain_limits[2],
                const unsigned int file_output_freq,
                const int num_resolutions_after_find = 3,
                const std::vector<int>& ntheta_array = {8, 16, 32},
                const std::vector<int>& nphi_array   = {16, 32, 64},
                const int enable_eta_varying_alg     = 0,
                const int verbosity_level            = 1)
        : num_horizons_(n_horizons),
          is_bbh_(is_binary_black_hole),
          num_resolutions_multigrid_(n_resolutions_multigrid),
          bah_verbosity_level_(verbosity_level),
          bah_enable_eta_varying_alg_for_precision_common_horizon_(
              enable_eta_varying_alg),
          max_ntheta_(ntheta_max),
          max_nphi_(nphi_max),
          out_dir_(save_directory),
          variable_indices_(indices_extract),
          transform_(transform),
          file_output_freq_(file_output_freq) {
        allocate_data_structures();

        grid_limits_[0]   = grid_limits[0];
        grid_limits_[1]   = grid_limits[1];
        domain_limits_[0] = domain_limits[0];
        domain_limits_[1] = domain_limits[1];

        // then set parameters for each horizon of interest, fills in the data
        // as asked for by the user, and fills any un-filled values with good
        // defaults, note it just ignores any "additional" values past the
        // number of input horizons
        fill_vector_with_defaults(m_scale, bah_m_scale_, 0.0, num_horizons_);
        fill_vector_with_defaults(cfl_factor, bah_cfl_factor_, 1.0,
                                  num_horizons_);
        fill_vector_with_defaults(max_iterations, bah_max_iterations_, 1000000,
                                  num_horizons_);
        fill_vector_with_defaults(theta_l2_m_tol,
                                  bah_theta_l2_times_m_tolerance_, 1.e-2,
                                  num_horizons_);
        fill_vector_with_defaults(theta_linf_m_tol,
                                  bah_theta_linf_times_m_tolerance_, 1.e-5,
                                  num_horizons_);
        fill_vector_with_defaults(eta_damp_m, bah_eta_damping_times_m_, 7.0,
                                  num_horizons_);
        fill_vector_with_defaults(ko_strength, bah_ko_strength_, 0.0,
                                  num_horizons_);
        fill_vector_with_defaults(nr_interp_max, bah_nr_interp_max_, 48,
                                  num_horizons_);

        fill_vector_with_defaults(max_search_radius, bah_max_search_radius_,
                                  1.5, num_horizons_);

        // NOTE: ntheta_array and nphi_array *must* be the same size and must be
        // num_resolutions after_find
        num_resolutions_multigrid_ = num_resolutions_after_find;
        if (ntheta_array.size() != num_resolutions_multigrid_) {
            throw std::runtime_error(
                "ERROR: nTheta Array needs to be the same size as "
                "num_resolutions_after_find!");
        }
        if (nphi_array.size() != num_resolutions_multigrid_) {
            throw std::runtime_error(
                "ERROR: nPhi Array needs to be the same size as "
                "num_resolutions_after_find!");
        }
        // these are vector fields!
        ntheta_array_multigrid_ = ntheta_array;
        nphi_array_multigrid_   = nphi_array;

        // then initialize the data
        initialize_data(initial_x_center, initial_y_center, initial_z_center,
                        blackholes);

        initialized_ = true;
    }

    ~AEH_BHaHAHA() = default;

    void allocate_data_structures() {
        // initialize all internal structures that will store the persistent
        // data
        x_guess_              = std::vector<double>(num_horizons_, 0.0);
        y_guess_              = std::vector<double>(num_horizons_, 0.0);
        z_guess_              = std::vector<double>(num_horizons_, 0.0);
        r_min_guess_          = std::vector<double>(num_horizons_, 0.0);
        r_max_guess_          = std::vector<double>(num_horizons_, 0.0);

        x_center_m1_          = std::vector<double>(num_horizons_, 0.0);
        y_center_m1_          = std::vector<double>(num_horizons_, 0.0);
        z_center_m1_          = std::vector<double>(num_horizons_, 0.0);
        x_center_m2_          = std::vector<double>(num_horizons_, 0.0);
        y_center_m2_          = std::vector<double>(num_horizons_, 0.0);
        z_center_m2_          = std::vector<double>(num_horizons_, 0.0);
        x_center_m3_          = std::vector<double>(num_horizons_, 0.0);
        y_center_m3_          = std::vector<double>(num_horizons_, 0.0);
        z_center_m3_          = std::vector<double>(num_horizons_, 0.0);

        r_min_m1_             = std::vector<double>(num_horizons_, 0.0);
        r_min_m2_             = std::vector<double>(num_horizons_, 0.0);
        r_min_m3_             = std::vector<double>(num_horizons_, 0.0);
        r_max_m1_             = std::vector<double>(num_horizons_, 0.0);
        r_max_m2_             = std::vector<double>(num_horizons_, 0.0);
        r_max_m3_             = std::vector<double>(num_horizons_, 0.0);
        t_m1_                 = std::vector<double>(num_horizons_, 0.0);
        t_m2_                 = std::vector<double>(num_horizons_, 0.0);
        t_m3_                 = std::vector<double>(num_horizons_, 0.0);
        failed_last_find_     = std::vector<bool>(num_horizons_, false);
        failed_last_find_int_ = std::vector<int>(num_horizons_, 0);

        prev_horizon_m1_ =
            std::vector<double>(num_horizons_ * max_ntheta_ * max_nphi_);
        prev_horizon_m2_ =
            std::vector<double>(num_horizons_ * max_ntheta_ * max_nphi_);
        prev_horizon_m3_ =
            std::vector<double>(num_horizons_ * max_ntheta_ * max_nphi_);

        // turn on all horizons to start, bbh and shared will be deactivated
        // later
        bah_horizon_active_ = std::vector<int>(num_horizons_, 1);

        bah_use_fixed_radius_guess_on_full_sphere_ =
            std::vector<int>(num_horizons_, 1);
    }

    void initialize_data(const std::vector<double>& initial_x_center,
                         const std::vector<double>& initial_y_center,
                         const std::vector<double>& initial_z_center,
                         const std::vector<SimpleBlackHoleData>& blackholes) {
        for (unsigned int i = 0; i < num_horizons_; i++) {
            // initialize to initial guess
            x_center_m1_[i] = initial_x_center[i];
            y_center_m1_[i] = initial_y_center[i];
            z_center_m1_[i] = initial_z_center[i];
            t_m1_[i] = t_m2_[i] = t_m3_[i] = -1.0;

            // if bbh, force black hole to be set
            if (is_bbh_) {
                // i = 0 means we're on bh1
                if (i == 0) {
                    x_center_m1_[i] = blackholes[0].x;
                    y_center_m1_[i] = blackholes[0].y;
                    z_center_m1_[i] = blackholes[0].z;

                } else if (i == 1) {
                    x_center_m1_[i] = blackholes[1].x;
                    y_center_m1_[i] = blackholes[1].y;
                    z_center_m1_[i] = blackholes[1].z;
                }
                // assume third is zero, which is initalized
            }
            // TODO: otherwise we need to set it to the centers as parameters
        }

        // for binary black holes, activate two inspiral bh's, deactivate common
        // horizon
        if (is_bbh_) {
            bah_horizon_active_[0] = 1;
            bah_horizon_active_[1] = 1;

            // shared horizon inactive
            bah_horizon_active_[2] = 0;

            // set based on mass scale
            bah_m_scale_[0]        = blackholes[0].mass;
            bah_m_scale_[1]        = blackholes[1].mass;

            // NOTE: this is a temporary measure until we verify this is
            // correct, but the common horizon can just be a simple addition of
            // the masses
            bah_m_scale_[2]        = blackholes[0].mass + blackholes[1].mass;
        }

        // this is everything that is done the first time, we're not initialized
        // the creation of the param_structs and copying of data now happens in
        // find_horizons
        initialized_ = true;
    }

    void clear_bhahaha_param_structs() {
        // this fully clears and deletes everything in the bah_param_data_
        // vector
        bha_param_data_.clear();

        // the main find_horizons loop deletes allocated memory, everything else
        // is persistently stored
    }

    void create_bhahaha_param_structs() {
        // this allocates our vector to store the data we're using
        bha_param_data_ =
            std::vector<bhahaha_params_and_data_struct>(num_horizons_);

        for (int which_horizon = 0; which_horizon < num_horizons_;
             which_horizon++) {
            auto* bah_params_and_data = &bha_param_data_[which_horizon];

            // call the setting of parameters, which takes the data we're
            // storing here
            set_bhahaha_parameters(bah_params_and_data, which_horizon);
        }
    }

    void set_bhahaha_parameters(
        bhahaha_params_and_data_struct* bah_params_and_data,
        int which_horizon) {
        // initialize each data structure
        bah_poisoning_set_inputs(bah_params_and_data);
        // input metric data needs to be null!
        bah_params_and_data->input_metric_data        = NULL;

        // basic horizon metadata
        bah_params_and_data->which_horizon            = which_horizon + 1;
        bah_params_and_data->num_horizons             = num_horizons_;
        // this is the iteration number!
        bah_params_and_data->iteration_external_input = 0;
        bah_params_and_data->num_resolutions_multigrid =
            num_resolutions_multigrid_;

        for (int res = 0; res < num_resolutions_multigrid_; res++) {
            bah_params_and_data->Ntheta_array_multigrid[res] =
                ntheta_array_multigrid_[res];
            bah_params_and_data->Nphi_array_multigrid[res] =
                nphi_array_multigrid_[res];
        }

        bah_params_and_data->verbosity_level = bah_verbosity_level_;
        bah_params_and_data
            ->enable_eta_varying_alg_for_precision_common_horizon =
            bah_enable_eta_varying_alg_for_precision_common_horizon_;

        // per-horizon parameters
        bah_params_and_data->M_scale    = bah_m_scale_[which_horizon];
        bah_params_and_data->cfl_factor = bah_cfl_factor_[which_horizon];
        bah_params_and_data->use_fixed_radius_guess_on_full_sphere =
            bah_use_fixed_radius_guess_on_full_sphere_[which_horizon];
        bah_params_and_data->max_iterations =
            bah_max_iterations_[which_horizon];
        bah_params_and_data->Theta_L2_times_M_tolerance =
            bah_theta_l2_times_m_tolerance_[which_horizon];
        bah_params_and_data->Theta_Linf_times_M_tolerance =
            bah_theta_linf_times_m_tolerance_[which_horizon];
        bah_params_and_data->eta_damping_times_M =
            bah_eta_damping_times_m_[which_horizon];
        bah_params_and_data->KO_strength = bah_ko_strength_[which_horizon];

        // then make sure we get the data from persistent
        transfer_to_bhahaha_from_persistent(bah_params_and_data);
    }

    void transfer_to_bhahaha_from_persistent(
        bhahaha_params_and_data_struct* bhahaha_params_and_data) {
        const int which_horizon = bhahaha_params_and_data->which_horizon - 1;

        bhahaha_params_and_data->x_center_m1 = x_center_m1_[which_horizon];
        bhahaha_params_and_data->y_center_m1 = y_center_m1_[which_horizon];
        bhahaha_params_and_data->z_center_m1 = z_center_m1_[which_horizon];
        bhahaha_params_and_data->x_center_m2 = x_center_m2_[which_horizon];
        bhahaha_params_and_data->y_center_m2 = y_center_m2_[which_horizon];
        bhahaha_params_and_data->z_center_m2 = z_center_m2_[which_horizon];
        bhahaha_params_and_data->x_center_m3 = x_center_m3_[which_horizon];
        bhahaha_params_and_data->y_center_m3 = y_center_m3_[which_horizon];
        bhahaha_params_and_data->z_center_m3 = z_center_m3_[which_horizon];
        bhahaha_params_and_data->t_m1        = t_m1_[which_horizon];
        bhahaha_params_and_data->t_m2        = t_m2_[which_horizon];
        bhahaha_params_and_data->t_m3        = t_m3_[which_horizon];
        bhahaha_params_and_data->r_min_m1    = r_min_m1_[which_horizon];
        bhahaha_params_and_data->r_max_m1    = r_max_m1_[which_horizon];
        bhahaha_params_and_data->r_min_m2    = r_min_m2_[which_horizon];
        bhahaha_params_and_data->r_max_m2    = r_max_m2_[which_horizon];
        bhahaha_params_and_data->r_min_m3    = r_min_m3_[which_horizon];
        bhahaha_params_and_data->r_max_m3    = r_max_m3_[which_horizon];
        // NOTE: we're setting these as *pointers*, so they'll be updated every
        // time
        bhahaha_params_and_data->prev_horizon_m1 =
            &prev_horizon_m1_[which_horizon * max_ntheta_ * max_nphi_];
        bhahaha_params_and_data->prev_horizon_m2 =
            &prev_horizon_m2_[which_horizon * max_ntheta_ * max_nphi_];
        bhahaha_params_and_data->prev_horizon_m3 =
            &prev_horizon_m3_[which_horizon * max_ntheta_ * max_nphi_];
    }

    void transfer_to_persistent_from_bhahaha(
        bhahaha_params_and_data_struct* bhahaha_params_and_data) {
        const int which_horizon = bhahaha_params_and_data->which_horizon - 1;

        x_center_m1_[which_horizon] = bhahaha_params_and_data->x_center_m1;
        y_center_m1_[which_horizon] = bhahaha_params_and_data->y_center_m1;
        z_center_m1_[which_horizon] = bhahaha_params_and_data->z_center_m1;
        x_center_m2_[which_horizon] = bhahaha_params_and_data->x_center_m2;
        y_center_m2_[which_horizon] = bhahaha_params_and_data->y_center_m2;
        z_center_m2_[which_horizon] = bhahaha_params_and_data->z_center_m2;
        x_center_m3_[which_horizon] = bhahaha_params_and_data->x_center_m3;
        y_center_m3_[which_horizon] = bhahaha_params_and_data->y_center_m3;
        z_center_m3_[which_horizon] = bhahaha_params_and_data->z_center_m3;
        t_m1_[which_horizon]        = bhahaha_params_and_data->t_m1;
        t_m2_[which_horizon]        = bhahaha_params_and_data->t_m2;
        t_m3_[which_horizon]        = bhahaha_params_and_data->t_m3;
        r_min_m1_[which_horizon]    = bhahaha_params_and_data->r_min_m1;
        r_max_m1_[which_horizon]    = bhahaha_params_and_data->r_max_m1;
        r_min_m2_[which_horizon]    = bhahaha_params_and_data->r_min_m2;
        r_max_m2_[which_horizon]    = bhahaha_params_and_data->r_max_m2;
        r_min_m3_[which_horizon]    = bhahaha_params_and_data->r_min_m3;
        r_max_m3_[which_horizon]    = bhahaha_params_and_data->r_max_m3;
        // NOTE: due to the pointer nature above, we don't have to worry about
        // copying back the prev_horizon data
    }

    void bah_sum_shared_arrays(const ot::Mesh* mesh);

    void fill_domain_coords(const int which_horizon, const int n_r,
                            const int n_theta, const int n_phi,
                            const double* radii, const double x_center,
                            const double y_center, const double z_center,
                            std::vector<double>& domain_coords);

    void interpolate_metric_data(const ot::Mesh* mesh, const double** varData,
                                 const int which_horizon, const int which_rank,
                                 const int n_r, const int n_theta,
                                 const int n_phi, const double* radii,
                                 const double x_center, const double y_center,
                                 const double z_center,
                                 double* input_metric_data);

    void find_horizons(const ot::Mesh* mesh, const double** var,
                       const unsigned int current_step,
                       const double current_time,
                       const std::vector<Point> tracked_location_data = {});

    void synchronize_to_root(const ot::Mesh* mesh, const int targetProc = 0);

    void create_checkpoint(
        const ot::Mesh* mesh,
        const std::string& checkpoint_output = "aeh_bha_chkpt.json");

    void restore_checkpoint(
        const ot::Mesh* mesh,
        const std::string& checkpoint_file = "aeh_bha_chkpt.json");
};

}  // namespace dendro_aeh
