#include "BHaHAHA.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Time-symmetric Schwarzschild data in isotropic coordinates, M=1:
// gamma_ij = (1 + 1/(2r))^4 delta_ij, K_ij = 0.
// The exact horizon has coordinate radius 0.5 and area 16*pi.
int main(void) {
  const int max_nr = 96, nt = 32, np = 64;
  BHA_REAL radii[96];
  bhahaha_params_and_data_struct p = {0};
  bhahaha_diagnostics_struct d = {0};
  bah_poisoning_set_inputs(&p);
  p.time_external_input = 0;
  p.iteration_external_input = 0;
  bah_radial_grid_cell_centered_set_up(max_nr, 1.0, 0.0, 1.0,
      &p.Nr_external_input, &p.r_min_external_input, &p.dr_external_input, radii);
  const size_t points = (size_t)p.Nr_external_input * nt * np;
  p.input_metric_data = calloc(NUM_EXT_INPUT_CARTESIAN_GFS * points, sizeof(BHA_REAL));
  p.prev_horizon_m1 = calloc(nt * np, sizeof(BHA_REAL));
  p.prev_horizon_m2 = calloc(nt * np, sizeof(BHA_REAL));
  p.prev_horizon_m3 = calloc(nt * np, sizeof(BHA_REAL));
  if (!p.input_metric_data || !p.prev_horizon_m1 || !p.prev_horizon_m2 || !p.prev_horizon_m3)
    return 2;
  for (int k = 0; k < np; ++k)
    for (int j = 0; j < nt; ++j)
      for (int i = 0; i < p.Nr_external_input; ++i) {
        const size_t q = i + p.Nr_external_input * (j + nt * k);
        const BHA_REAL metric = pow(1.0 + 0.5 / radii[i], 4);
        p.input_metric_data[q + points * INTERP_GAMMADDXXGF] = metric;
        p.input_metric_data[q + points * INTERP_GAMMADDYYGF] = metric;
        p.input_metric_data[q + points * INTERP_GAMMADDZZGF] = metric;
      }
  p.num_resolutions_multigrid = 3;
  for (int i = 0; i < 3; ++i) {
    p.Ntheta_array_multigrid[i] = 8 << i;
    p.Nphi_array_multigrid[i] = 16 << i;
  }
  p.use_fixed_radius_guess_on_full_sphere = 1;
  p.cfl_factor = 1.05;
  p.M_scale = 1;
  p.eta_damping_times_M = 1.6;
  p.KO_strength = 0;
  p.max_iterations = 10000;
  p.Theta_Linf_times_M_tolerance = 1e-2;
  p.Theta_L2_times_M_tolerance = 2e-5;
  p.which_horizon = 0;
  p.num_horizons = 1;
  p.verbosity_level = 0;
  p.enable_eta_varying_alg_for_precision_common_horizon = 0;
  bah_poisoning_check_inputs(&p);
  const int status = bah_find_horizon(&p, &d);
  printf("status=%d area=%.12g radius=%.12g Theta_L2=%.12g\n",
      status, d.area, d.mean_coord_radius_wrt_centroid, d.Theta_L2_times_M);
  // Require 0.1% agreement with the analytic horizon and requested convergence.
  const int ok = status == BHAHAHA_SUCCESS && isfinite(d.area) &&
      fabs(d.area / (16 * M_PI) - 1) < 1e-3 &&
      fabs(d.mean_coord_radius_wrt_centroid / 0.5 - 1) < 1e-3 &&
      isfinite(d.Theta_L2_times_M) && d.Theta_L2_times_M < p.Theta_L2_times_M_tolerance;
  free(p.input_metric_data);
  free(p.prev_horizon_m1);
  free(p.prev_horizon_m2);
  free(p.prev_horizon_m3);
  return !ok;
}
