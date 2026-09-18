#include "BHaH_defines.h"
#include "interpolation_lagrange_uniform.h"

#define DEBUG
#define INTERP_ORDER (2 * NinterpGHOSTS + 1) // Interpolation order (== number of points in stencil, in each dimension).
#define DST_IDX4(g, i, j, k) ((i) + dst_Nxx_plus_2NGHOSTS0 * ((j) + dst_Nxx_plus_2NGHOSTS1 * ((k) + dst_Nxx_plus_2NGHOSTS2 * (g))))
#define SRC_IDX4(g, i, j, k) ((i) + src_Nxx_plus_2NGHOSTS0 * ((j) + src_Nxx_plus_2NGHOSTS1 * ((k) + src_Nxx_plus_2NGHOSTS2 * (g))))
// #define STANDALONE

#pragma GCC optimize("unroll-loops")

/**
 * Interpolates data from the source grid (external input) to the destination grid (interp_src) using 2D Lagrange interpolation.
 *
 * Perform 2D interpolation on a uniform grid in the theta and phi directions,
 * at all (perfectly overlapping) radial points on the src and dst grids.
 *
 * @param[in,out] commondata A structure containing grid parameters and data arrays for both source and destination grids.
 *
 * @return Returns an error code if any pointer is NULL or if the interpolation order exceeds grid boundaries. Returns 0 on success.
 *
 * Notes:
 * - Assumes that theta and phi values are uniform and consistent between source and destination grids.
 * - The stencil size for interpolation is defined by the INTERP_ORDER macro.
 * - The destination and src grid radial points overlap.
 */
int bah_interpolation_2d_external_input_to_interp_src_grid(commondata_struct *restrict commondata) {
  // UNPACK PARAMETERS:
  // src = external_input
  // dst = interp_src
  const BHA_REAL src_dxx1 = commondata->external_input_dxx1;
  const BHA_REAL src_dxx2 = commondata->external_input_dxx2;
  const BHA_REAL src_invdxx1 = commondata->external_input_invdxx1;
  const BHA_REAL src_invdxx2 = commondata->external_input_invdxx2;
  const int src_Nxx_plus_2NGHOSTS0 = commondata->external_input_Nxx_plus_2NGHOSTS0;
  const int src_Nxx_plus_2NGHOSTS1 = commondata->external_input_Nxx_plus_2NGHOSTS1;
  const int src_Nxx_plus_2NGHOSTS2 = commondata->external_input_Nxx_plus_2NGHOSTS2;
  const BHA_REAL *restrict src_gfs = commondata->external_input_gfs;

  const int dst_Nxx1 = commondata->interp_src_Nxx1;
  const int dst_Nxx2 = commondata->interp_src_Nxx2;
  const int dst_Nxx_plus_2NGHOSTS0 = commondata->interp_src_Nxx_plus_2NGHOSTS0;
  const int dst_Nxx_plus_2NGHOSTS1 = commondata->interp_src_Nxx_plus_2NGHOSTS1;
  const int dst_Nxx_plus_2NGHOSTS2 = commondata->interp_src_Nxx_plus_2NGHOSTS2;
  BHA_REAL *restrict dst_gfs = commondata->interp_src_gfs;

  // Compute (1/(src_dtheta * src_dphi))^(INTERP_ORDER-1) normalization factor, as pow() is expensive.
  // Orig code: params.inv_dxx0_dxx1_ORDERm1 = pow(1.0 / (dxx0 * dxx1), ORDER - 1);
  const BHA_REAL src_invdxx12_INTERP_ORDERm1 = pow(commondata->external_input_dxx1 * commondata->external_input_dxx2, -(INTERP_ORDER - 1));

  const BHA_REAL *restrict src_r_theta_phi[3] = {commondata->external_input_r_theta_phi[0], commondata->external_input_r_theta_phi[1],
                                             commondata->external_input_r_theta_phi[2]};

  // Perform debug checks if DEBUG is defined
  if (src_r_theta_phi[1] == NULL || src_r_theta_phi[2] == NULL || commondata->external_input_gfs == NULL ||
      commondata->interp_src_r_theta_phi[1] == NULL || commondata->interp_src_r_theta_phi[2] == NULL || commondata->interp_src_gfs == NULL)
    return INTERP2D_EXT_TO_INTERPSRC_NULL_PTRS;

  if (INTERP_ORDER > commondata->external_input_Nxx1 + 2 * NinterpGHOSTS || INTERP_ORDER > commondata->external_input_Nxx2 + 2 * NinterpGHOSTS)
    return INTERP2D_EXT_TO_INTERPSRC_INTERP_ORDER_GT_NXX_PLUS_2NINTERPGHOSTS12;

  // Precompute inverse denominators for Lagrange interpolation coefficients to optimize performance.
  BHA_REAL inv_denom[INTERP_ORDER];
  compute_inv_denom(INTERP_ORDER, inv_denom);

  // Perform interpolation for each destination angular point (theta, phi)
  const BHA_REAL xxmin_incl_ghosts1 = src_r_theta_phi[1][0];
  const BHA_REAL xxmin_incl_ghosts2 = src_r_theta_phi[2][0];
  int error_flag = BHAHAHA_SUCCESS;
  int i0_min_shift = 0;
  if (commondata->bhahaha_params_and_data->r_min_external_input == 0)
    i0_min_shift = NGHOSTS;
#pragma omp parallel for
  for (int ir = i0_min_shift; ir < dst_Nxx_plus_2NGHOSTS0; ir++) {
    for (int iphi = NGHOSTS; iphi < dst_Nxx2 + NGHOSTS; iphi++) { // Ignore ghost zones; these can be handled separately.
      const BHA_REAL phi_dst = commondata->interp_src_r_theta_phi[2][iphi];
      for (int itheta = NGHOSTS; itheta < dst_Nxx1 + NGHOSTS; itheta++) {
        const BHA_REAL theta_dst = commondata->interp_src_r_theta_phi[1][itheta];

        // printf("destination: (%e, %e, %e)\n", commondata->interp_src_r_theta_phi[0][ir], commondata->interp_src_r_theta_phi[1][itheta],
        //        commondata->interp_src_r_theta_phi[2][iphi]);

        // Calculate the central indices for the stencil in the theta & phi directions
        const int idx_center_th = (int)((theta_dst - xxmin_incl_ghosts1) * src_invdxx1 + 0.5);
        const int idx_center_ph = (int)((phi_dst - xxmin_incl_ghosts2) * src_invdxx2 + 0.5);

        {
          // Ensure the stencil is within valid grid bounds
          if ((idx_center_th - NinterpGHOSTS < 0) || (idx_center_th + NinterpGHOSTS >= src_Nxx_plus_2NGHOSTS1) ||
              (idx_center_ph - NinterpGHOSTS < 0) || (idx_center_ph + NinterpGHOSTS >= src_Nxx_plus_2NGHOSTS2)) {
#ifdef DEBUG
            fprintf(stderr, "ERROR: Interpolation stencil exceeds grid boundaries. %d %d %d %d\n", (idx_center_th - NinterpGHOSTS < 0),
                    (idx_center_th + NinterpGHOSTS >= src_Nxx_plus_2NGHOSTS1), (idx_center_ph - NinterpGHOSTS < 0),
                    (idx_center_ph + NinterpGHOSTS >= src_Nxx_plus_2NGHOSTS2));
            fprintf(stderr, "(th_dst,ph_dst) = (%.6f,%.6f). (itheta,iphi) = (%d,%d).\n", theta_dst, phi_dst, itheta, iphi);
            fprintf(stderr, "Grid bounds along theta direction: [0, %d], stencil indices: [%d, %d]\n", src_Nxx_plus_2NGHOSTS1 - 1,
                    idx_center_th - NinterpGHOSTS, idx_center_th + NinterpGHOSTS);
            fprintf(stderr, "Grid bounds along phi direction: [0, %d], stencil indices: [%d, %d]\n", src_Nxx_plus_2NGHOSTS2 - 1,
                    idx_center_ph - NinterpGHOSTS, idx_center_ph + NinterpGHOSTS);
            fprintf(stderr, "Ensure that the destination point is within grid bounds or adjust the interpolation stencil.\n");
#endif // DEBUG
#pragma omp critical
            {
              error_flag = INTERP2D_EXT_TO_INTERPSRC_HORIZON_OUT_OF_BOUNDS;
            }
            continue; // Skip further work for this iteration
          } // END IF: stencil in bounds

#ifdef DEBUG
          // Verify that the central index is the closest grid point to the destination radius
          const BHA_REAL TOLERANCE = 1e-13; // The TOLERANCE is needed, as we often interpolate to points exactly midway between src points.
          if (fabs(src_r_theta_phi[1][idx_center_th] - theta_dst) > src_dxx1 * (0.5 + TOLERANCE)) {
            fprintf(stderr, "ERROR: theta center index too far from destination point! %.15e > %.15e\n",
                    fabs(src_r_theta_phi[1][idx_center_th] - theta_dst), src_dxx1 * (0.5 + TOLERANCE));
          }
          if (fabs(src_r_theta_phi[2][idx_center_ph] - phi_dst) > src_dxx2 * (0.5 + TOLERANCE)) {
            fprintf(stderr, "ERROR: phi center index too far from destination point! %.15e > %.15e\n",
                    fabs(src_r_theta_phi[2][idx_center_ph] - phi_dst), src_dxx2 * (0.5 + TOLERANCE));
          } // END IF: central index is properly centered
#endif // DEBUG
        } // END BLOCK: theta/phi stencil bounds and center-index sanity checks

        const int base_idx_th = idx_center_th - NinterpGHOSTS;
        const int base_idx_ph = idx_center_ph - NinterpGHOSTS;

        // Step 1: Precompute all differences
        BHA_REAL diffs_th[INTERP_ORDER], diffs_ph[INTERP_ORDER];
        compute_diffs_xi(INTERP_ORDER, theta_dst, &src_r_theta_phi[1][base_idx_th], diffs_th);
        compute_diffs_xi(INTERP_ORDER, phi_dst, &src_r_theta_phi[2][base_idx_ph], diffs_ph);

        // Step 2: Precompute combined Lagrange coefficients to reduce computations
        BHA_REAL lagrange_basis_coeffs_th[INTERP_ORDER], lagrange_basis_coeffs_ph[INTERP_ORDER];
        compute_lagrange_basis_coeffs_xi(INTERP_ORDER, inv_denom, diffs_th, lagrange_basis_coeffs_th);
        compute_lagrange_basis_coeffs_xi(INTERP_ORDER, inv_denom, diffs_ph, lagrange_basis_coeffs_ph);
        BHA_REAL coeff_2d[INTERP_ORDER][INTERP_ORDER];
        for (int iph = 0; iph < INTERP_ORDER; iph++) {
          const BHA_REAL coeff_ph_i = lagrange_basis_coeffs_ph[iph];
          for (int ith = 0; ith < INTERP_ORDER; ith++) {
            coeff_2d[iph][ith] = coeff_ph_i * lagrange_basis_coeffs_th[ith];
          } // END LOOP: for ith over theta stencil
        } // END LOOP: for iph over phi stencil

        // Step 3: Perform the 2D Lagrange interpolation, optimizing memory accesses and enabling vectorization
        BHA_REAL sum[NUM_EXT_INPUT_CONFORMAL_GFS];
        for (int gf = 0; gf < NUM_EXT_INPUT_CONFORMAL_GFS; gf++) {
          sum[gf] = 0.0;
        } // END LOOP: for which_gf over grid functions
        for (int iph = 0; iph < INTERP_ORDER; iph++) {
          const int idx_ph = base_idx_ph + iph;
          for (int ith = 0; ith < INTERP_ORDER; ith++) {
            const int idx_th = base_idx_th + ith;
            const BHA_REAL coeff = coeff_2d[iph][ith];
#pragma omp simd
            for (int gf = 0; gf < NUM_EXT_INPUT_CONFORMAL_GFS; gf++) {
              sum[gf] += src_gfs[SRC_IDX4(gf, ir, idx_th, idx_ph)] * coeff;
            } // END LOOP: for which_gf over source gridfunctions
          } // END LOOP: for ith over theta stencil
        } // END LOOP: for iph over phi stencil
        for (int gf = 0; gf < NUM_EXT_INPUT_CONFORMAL_GFS; gf++) {
          dst_gfs[DST_IDX4(gf, ir, itheta, iphi)] = sum[gf] * src_invdxx12_INTERP_ORDERm1;
        } // END LOOP: for which_gf over destination gridfunctions
      } // END LOOP: for itheta over destination theta
    } // END LOOP: for iphi over destination phi
  } // END LOOP: for ir over radial shells
  return error_flag;
} // END FUNCTION: bah_interpolation_2d_external_input_to_interp_src_grid

#pragma GCC reset_options // Reset compiler optimizations after the function

#ifdef STANDALONE
#include <omp.h>

#define RADIAL_EXTENT 1.0
#define NUM_RESOLUTIONS 3

/**
 * Analytic function used for initializing the source grid function.
 *
 * Calculates the product of \( r^2 \sin(r) \sin(\theta) \cos(\phi) \).
 *
 * @param r       Radial coordinate.
 * @param theta   Angle in the theta direction.
 * @param phi     Angle in the phi direction.
 *
 * @return        The computed value of \( r^2 \sin(r) \sin(\theta) \cos(\phi) \).
 */
static inline BHA_REAL analytic_function(BHA_REAL r, BHA_REAL theta, BHA_REAL phi) { return r * r * sin(r) * sin(theta) * cos(phi); }

/**
 * Initializes the coordinate arrays for a uniform grid.
 *
 * Allocates and sets up the r, theta, and phi coordinate values based on the number of ghost zones and grid dimensions.
 * This version mirrors the 1D code's structure, returning an error code on failure.
 *
 * @param N_r                    Number of grid points in the radial direction.
 * @param N_theta                Number of grid points in the theta direction.
 * @param N_phi                  Number of grid points in the phi direction.
 * @param r_theta_phi            Arrays to store coordinate values for r, theta, and phi.
 * @param dxx                    Array to store grid spacings in r, theta, and phi directions.
 * @param Nxx_plus_2NGHOSTS      Array containing the total number of grid points in each direction, including ghost zones.
 *
 * @return                       0 on success, -1 on memory allocation failure.
 */
int initialize_coordinates(const int N_r, const int N_theta, const int N_phi, BHA_REAL *r_theta_phi[3], BHA_REAL dxx[3], const int Nxx_plus_2NGHOSTS[3]) {
  dxx[0] = RADIAL_EXTENT / N_r;
  dxx[1] = M_PI / N_theta;
  dxx[2] = 2.0 * M_PI / N_phi;

  r_theta_phi[0] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * Nxx_plus_2NGHOSTS[0]);
  r_theta_phi[1] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * Nxx_plus_2NGHOSTS[1]);
  r_theta_phi[2] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * Nxx_plus_2NGHOSTS[2]);

  if (!r_theta_phi[0] || !r_theta_phi[1] || !r_theta_phi[2]) {
    free(r_theta_phi[0]);
    free(r_theta_phi[1]);
    free(r_theta_phi[2]);
    return -1;
  } // END IF: coordinate-array allocation failed

  // Populate coordinate arrays. Note the r-direction is vertex-centered like the 1D example.
  for (int i = 0; i < Nxx_plus_2NGHOSTS[0]; i++) {
    r_theta_phi[0][i] = (i - NGHOSTS) * dxx[0];
  } // END LOOP: for i0 over radial direction
  for (int i = 0; i < Nxx_plus_2NGHOSTS[1]; i++) {
    r_theta_phi[1][i] = (i - NGHOSTS + 0.5) * dxx[1];
  } // END LOOP: for i1 over theta direction
  for (int i = 0; i < Nxx_plus_2NGHOSTS[2]; i++) {
    r_theta_phi[2][i] = -M_PI + (i - NGHOSTS + 0.5) * dxx[2];
  } // END LOOP: for i2 over phi direction
  return 0;
} // END FUNCTION: initialize_coordinates

/**
 * Initializes the source grid function using a pre-optimized analytic function evaluation.
 *
 * This version mirrors the 1D code's approach by pre-calculating r-dependent values
 * and using a nested loop structure optimized for memory access.
 *
 * @param src_Nxx_plus_2NGHOSTS  Array containing the total number of grid points in each direction, including ghost zones.
 * @param src_r_theta_phi        Array containing coordinate values for r, theta, and phi.
 * @param src_gf                 Array to store the initialized source grid function values.
 *
 * @return                       0 on success, -1 on memory allocation failure.
 */
int initialize_src_gf(const int src_Nxx_plus_2NGHOSTS[3], BHA_REAL *src_r_theta_phi[3], BHA_REAL *src_gf) {
  const int src_Nxx_plus_2NGHOSTS0 = src_Nxx_plus_2NGHOSTS[0];
  const int src_Nxx_plus_2NGHOSTS1 = src_Nxx_plus_2NGHOSTS[1];
  const int src_Nxx_plus_2NGHOSTS2 = src_Nxx_plus_2NGHOSTS[2];

  BHA_REAL *r_func_values = (BHA_REAL *)malloc(sizeof(BHA_REAL) * src_Nxx_plus_2NGHOSTS0);
  if (!r_func_values) {
    fprintf(stderr, "malloc failed for r_func_values.\n");
    return -1;
  } // END IF: r_func_values allocation failed

  for (int i = 0; i < src_Nxx_plus_2NGHOSTS0; i++) {
    const BHA_REAL r = src_r_theta_phi[0][i];
    r_func_values[i] = r * r * sin(r);
  } // END LOOP: for i over precomputed r-dependent values

#pragma omp parallel for collapse(2)
  for (int gf = 0; gf < NUM_EXT_INPUT_CONFORMAL_GFS; gf++) {
    for (int i2 = 0; i2 < src_Nxx_plus_2NGHOSTS2; i2++) {
      const BHA_REAL cosphi = cos(src_r_theta_phi[2][i2]);
      for (int i1 = 0; i1 < src_Nxx_plus_2NGHOSTS1; i1++) {
        const BHA_REAL sintheta = sin(src_r_theta_phi[1][i1]);
        for (int i0 = 0; i0 < src_Nxx_plus_2NGHOSTS0; i0++) {
          src_gf[SRC_IDX4(gf, i0, i1, i2)] = r_func_values[i0] * sintheta * cosphi;
        } // END LOOP: for i0 over radial grid
      } // END LOOP: for i1 over theta grid
    } // END LOOP: for i2 over phi grid
  } // END LOOP: for which_gf over grid functions
  free(r_func_values);
  return 0;
} // END FUNCTION: initialize_src_gf

/**
 * Main function to execute the standalone 2D interpolation and perform convergence validation tests.
 *
 * Adopts the structure of the 1D code, including error handling with a goto cleanup pattern,
 * designated initializers, and resource management within a loop.
 *
 * @return EXIT_SUCCESS on successful execution.
 *         EXIT_FAILURE or specific error codes if memory allocation or interpolation fails.
 */
int main() {
  int return_code = EXIT_SUCCESS;

  // Source grid data (re-allocated each loop)
  BHA_REAL *src_r_theta_phi[3] = {NULL, NULL, NULL};
  BHA_REAL *src_gf = NULL;

  // Destination grid data (allocated once)
  BHA_REAL *dst_r_theta_phi[3] = {NULL, NULL, NULL};
  BHA_REAL *dst_gfs = NULL;

  const int N_r = 64;
  const int Ntheta_src_arr[NUM_RESOLUTIONS] = {16, 32, 64};
  const int Nphi_src_arr[NUM_RESOLUTIONS] = {32, 64, 128};
  const int Ntheta_dst = 320;
  const int Nphi_dst = 640;

  BHA_REAL h_arr[NUM_RESOLUTIONS], error_L2_norm[NUM_RESOLUTIONS];

  // Set up destination grid once. Define scalar variables for the DST_IDX4 macro.
  const int dst_Nxx_plus_2NGHOSTS0 = N_r + 2 * NGHOSTS;
  const int dst_Nxx_plus_2NGHOSTS1 = Ntheta_dst + 2 * NGHOSTS;
  const int dst_Nxx_plus_2NGHOSTS2 = Nphi_dst + 2 * NGHOSTS;
  const int dst_Nxx_plus_2NGHOSTS[3] = {dst_Nxx_plus_2NGHOSTS0, dst_Nxx_plus_2NGHOSTS1, dst_Nxx_plus_2NGHOSTS2};
  BHA_REAL dst_dxx[3];
  if (initialize_coordinates(N_r, Ntheta_dst, Nphi_dst, dst_r_theta_phi, dst_dxx, dst_Nxx_plus_2NGHOSTS) != 0) {
    fprintf(stderr, "Memory allocation failed for destination coordinate arrays.\n");
    return_code = EXIT_FAILURE;
    goto cleanup;
  } // END IF: destination coordinate initialization failed

  const size_t total_dst_pts = (size_t)dst_Nxx_plus_2NGHOSTS0 * dst_Nxx_plus_2NGHOSTS1 * dst_Nxx_plus_2NGHOSTS2;
  dst_gfs = (BHA_REAL *)malloc(sizeof(BHA_REAL) * NUM_EXT_INPUT_CONFORMAL_GFS * total_dst_pts);
  if (dst_gfs == NULL) {
    fprintf(stderr, "Memory allocation failed for destination grid functions.\n");
    return_code = EXIT_FAILURE;
    goto cleanup;
  } // END IF: destination gridfunction allocation failed

  for (int res = 0; res < NUM_RESOLUTIONS; res++) {
    const int Ntheta_src = Ntheta_src_arr[res];
    const int Nphi_src = Nphi_src_arr[res];

    const int src_Nxx_plus_2NGHOSTS[3] = {N_r + 2 * NGHOSTS, Ntheta_src + 2 * NGHOSTS, Nphi_src + 2 * NGHOSTS};
    BHA_REAL src_dxx[3];
    if (initialize_coordinates(N_r, Ntheta_src, Nphi_src, src_r_theta_phi, src_dxx, src_Nxx_plus_2NGHOSTS) != 0) {
      fprintf(stderr, "Memory allocation failed for source coordinates at resolution %d.\n", res);
      return_code = EXIT_FAILURE;
      goto cleanup;
    } // END IF: source coordinate initialization failed

    const size_t total_src_pts = (size_t)src_Nxx_plus_2NGHOSTS[0] * src_Nxx_plus_2NGHOSTS[1] * src_Nxx_plus_2NGHOSTS[2];
    src_gf = (BHA_REAL *)malloc(sizeof(BHA_REAL) * NUM_EXT_INPUT_CONFORMAL_GFS * total_src_pts);
    if (!src_gf) {
      fprintf(stderr, "Memory allocation failed for source grid function at resolution %d.\n", res);
      return_code = EXIT_FAILURE;
      goto cleanup;
    } // END IF: source gridfunction allocation failed

    if (initialize_src_gf(src_Nxx_plus_2NGHOSTS, src_r_theta_phi, src_gf) != 0) {
      return_code = EXIT_FAILURE;
      goto cleanup;
    } // END IF: source gridfunction initialization failed

    // Use designated initializers for commondata_struct
    commondata_struct commondata = {
        .bhahaha_params_and_data = malloc(sizeof(bhahaha_params_and_data_struct)),
        .external_input_dxx1 = src_dxx[1],
        .external_input_dxx2 = src_dxx[2],
        .external_input_invdxx1 = 1.0 / src_dxx[1],
        .external_input_invdxx2 = 1.0 / src_dxx[2],
        .external_input_Nxx1 = Ntheta_src,
        .external_input_Nxx2 = Nphi_src,
        .external_input_Nxx_plus_2NGHOSTS0 = src_Nxx_plus_2NGHOSTS[0],
        .external_input_Nxx_plus_2NGHOSTS1 = src_Nxx_plus_2NGHOSTS[1],
        .external_input_Nxx_plus_2NGHOSTS2 = src_Nxx_plus_2NGHOSTS[2],
        .external_input_gfs = src_gf,
        .interp_src_Nxx1 = Ntheta_dst,
        .interp_src_Nxx2 = Nphi_dst,
        .interp_src_Nxx_plus_2NGHOSTS0 = dst_Nxx_plus_2NGHOSTS[0],
        .interp_src_Nxx_plus_2NGHOSTS1 = dst_Nxx_plus_2NGHOSTS[1],
        .interp_src_Nxx_plus_2NGHOSTS2 = dst_Nxx_plus_2NGHOSTS[2],
        .interp_src_gfs = dst_gfs,
    };
    commondata.bhahaha_params_and_data->r_min_external_input = 0.0;
    for (int i = 0; i < 3; i++) {
      commondata.external_input_r_theta_phi[i] = src_r_theta_phi[i];
      commondata.interp_src_r_theta_phi[i] = dst_r_theta_phi[i];
    } // END LOOP: for i over coordinate pointers in commondata

#ifdef _OPENMP
    double start_time = omp_get_wtime();
#endif
    int error_code = bah_interpolation_2d_external_input_to_interp_src_grid(&commondata);
#ifdef _OPENMP
    double end_time = omp_get_wtime();
#endif

    printf("--- Benchmarking for Resolution N_theta = %d, N_phi = %d ---\n", Ntheta_src, Nphi_src);
#ifdef _OPENMP
    // The number of points interpolated is the number of interior points on the destination grid.
    const long long num_interps = (long long)N_r * Ntheta_dst * Nphi_dst;
    printf("Interpolated %lld points in %.4f seconds.\n", num_interps, end_time - start_time);
    printf("Performance: %.4f million points per second.\n", (double)num_interps / (end_time - start_time) / 1e6);
#endif

    if (error_code != BHAHAHA_SUCCESS) {
      fprintf(stderr, "Interpolation failed with error code: %d at resolution %d.\n", error_code, res);
      free(commondata.bhahaha_params_and_data);
      continue; // Continue to next resolution
    } // END IF: interpolation failed
    free(commondata.bhahaha_params_and_data);

    h_arr[res] = src_dxx[1];
    BHA_REAL error_sum = 0.0;
#pragma omp parallel for reduction(+ : error_sum)
    for (int gf = 0; gf < NUM_EXT_INPUT_CONFORMAL_GFS; gf++) {
      for (int i0 = NGHOSTS; i0 < dst_Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        for (int i1 = NGHOSTS; i1 < dst_Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
          for (int i2 = NGHOSTS; i2 < dst_Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
            const BHA_REAL r = dst_r_theta_phi[0][i0];
            const BHA_REAL theta = dst_r_theta_phi[1][i1];
            const BHA_REAL phi = dst_r_theta_phi[2][i2];
            const BHA_REAL interpolated = dst_gfs[DST_IDX4(gf, i0, i1, i2)];
            const BHA_REAL exact = analytic_function(r, theta, phi);
            const BHA_REAL error = interpolated - exact;
            error_sum += error * error;
          } // END LOOP: for i2 over phi error calculation
        } // END LOOP: for i1 over theta error calculation
      } // END LOOP: for i0 over r error calculation
    } // END LOOP: for which_gf over grid functions for error calculation
    // Normalize L2 error by the number of interior destination points.
    error_L2_norm[res] = sqrt(error_sum / (N_r * Ntheta_dst * Nphi_dst * NUM_EXT_INPUT_CONFORMAL_GFS));
    printf("Resolution %d: N_theta = %d, N_phi = %d, h = %.5e, L2 error = %.5e\n", res, Ntheta_src, Nphi_src, h_arr[res], error_L2_norm[res]);

    // Free source-specific memory and nullify pointers
    for (int i = 0; i < 3; i++) {
      free(src_r_theta_phi[i]);
      src_r_theta_phi[i] = NULL;
    } // END LOOP: for i over source coordinate arrays
    free(src_gf);
    src_gf = NULL;
  } // END LOOP: for res over resolutions

  printf("\n--- Convergence Results ---\n");
  for (int res = 1; res < NUM_RESOLUTIONS; res++) {
    BHA_REAL observed_order = log(error_L2_norm[res - 1] / error_L2_norm[res]) / log(h_arr[res - 1] / h_arr[res]);
    printf("Observed order of convergence between resolutions %d and %d: %.2f\n", res - 1, res, observed_order);
  } // END LOOP: for res over resolutions for convergence
  printf("Expected order of convergence: %d\n", INTERP_ORDER);

cleanup:
  if (return_code != EXIT_SUCCESS)
    printf("\nAn error occurred. Cleaning up...\n");
  else
    printf("\nProgram finished successfully. Cleaning up...\n");

  // Free all potentially allocated memory
  for (int i = 0; i < 3; i++) {
    free(src_r_theta_phi[i]);
    free(dst_r_theta_phi[i]);
  } // END LOOP: for i over coordinate arrays to free
  free(src_gf);
  free(dst_gfs);

  return return_code;
} // END FUNCTION: main
#endif // STANDALONE
