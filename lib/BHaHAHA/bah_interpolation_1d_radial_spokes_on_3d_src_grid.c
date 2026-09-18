#include "BHaH_defines.h"
#include "interpolation_lagrange_uniform.h"

// #define STANDALONE
// Enables print statements:
// #define DEBUG
#define DST_IDX4(g, i, j, k) ((i) + dst_Nxx_plus_2NGHOSTS0 * ((j) + dst_Nxx_plus_2NGHOSTS1 * ((k) + dst_Nxx_plus_2NGHOSTS2 * (g))))
#define SRC_IDX4(g, i, j, k) ((i) + src_Nxx_plus_2NGHOSTS0 * ((j) + src_Nxx_plus_2NGHOSTS1 * ((k) + src_Nxx_plus_2NGHOSTS2 * (g))))
#define DST_IDX3(i, j, k) ((i) + dst_Nxx_plus_2NGHOSTS0 * ((j) + dst_Nxx_plus_2NGHOSTS1 * (k)))

#pragma GCC optimize("unroll-loops")

/**
 * Perform 1D radial Lagrange interpolation along the radial spokes of a 3D
 * spherical grid. It computes the interpolated values using Lagrange polynomials.
 *
 * @param[in] params           Pointer to destination grid parameters.
 * @param[in] commondata       Pointer to common data, including source grid parameters and source grid functions.
 * @param[in] dst_radii_aka_src_h_gf Array of destination radial points where interpolation is performed.
 * @param[out] dst_interp_gfs  Array to store the interpolated results.
 *
 * @return                     BHAHAHA_SUCCESS on success, or an appropriate error code on failure.
 *
 * @note
 * - Assumes that the source grid has uniform grid spacing in r, theta, and phi.
 * - The interpolation destination is radial index i = NGHOSTS.
 */
int bah_interpolation_1d_radial_spokes_on_3d_src_grid(const params_struct *restrict params, const commondata_struct *restrict commondata,
                                                      const BHA_REAL *restrict dst_radii_aka_src_h_gf, BHA_REAL *restrict dst_interp_gfs) {
  // Define the interpolation order based on the number of ghost zones.
  const int INTERP_ORDER = (2 * NinterpGHOSTS + 1); // Interpolation order corresponds to the number of points in the stencil per dimension.

  // UNPACK PARAMETERS
  // Source grid parameters from commondata.
  const int src_Nxx_plus_2NGHOSTS0 = commondata->interp_src_Nxx_plus_2NGHOSTS0;
  const int src_Nxx_plus_2NGHOSTS1 = commondata->interp_src_Nxx_plus_2NGHOSTS1;
  const int src_Nxx_plus_2NGHOSTS2 = commondata->interp_src_Nxx_plus_2NGHOSTS2;

  // Destination grid parameters from params.
  const int dst_Nxx1 = params->Nxx1;
  const int dst_Nxx2 = params->Nxx2;
  const int dst_Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int dst_Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int dst_Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  // Source grid coordinates for r, theta, and phi.
  const BHA_REAL *interp_src_r_theta_phi[3] = {commondata->interp_src_r_theta_phi[0], commondata->interp_src_r_theta_phi[1],
                                           commondata->interp_src_r_theta_phi[2]};
  const BHA_REAL xxmin_incl_ghosts0 = interp_src_r_theta_phi[0][0];

  // Precompute (1/(src_dr))^(INTERP_ORDER-1) normalization factor.
  const BHA_REAL src_invdxx0 = commondata->interp_src_invdxx0;
  const BHA_REAL src_invdxx0_INTERP_ORDERm1 = pow(commondata->interp_src_invdxx0, (INTERP_ORDER - 1));

  // Source grid functions.
  const BHA_REAL *interp_src_gfs = commondata->interp_src_gfs;

  // Perform sanity checks to ensure all required pointers are valid.
  if (interp_src_r_theta_phi[0] == NULL || interp_src_gfs == NULL || dst_radii_aka_src_h_gf == NULL || dst_interp_gfs == NULL)
    return INTERP1D_NULL_PTRS; // Return error if any pointer is NULL.

  // Ensure that the interpolation order does not exceed the source grid size.
  if (INTERP_ORDER > commondata->interp_src_Nxx0 + 2 * NinterpGHOSTS)
    return INTERP1D_INTERP_ORDER_GT_NXX_PLUS_2NINTERPGHOSTS0; // Return error if interpolation order is too high.

  // Initialize the error flag to track any interpolation issues.
  int error_flag = BHAHAHA_SUCCESS;

  // Precompute 1/denom coefficients for interpolation.
  BHA_REAL inv_denom[INTERP_ORDER];
  compute_inv_denom(INTERP_ORDER, inv_denom);

  // Parallelize the outer loops using OpenMP for better performance.
#pragma omp parallel for
  for (int iphi = NGHOSTS; iphi < dst_Nxx2 + NGHOSTS; iphi++) {         // Iterate over phi indices, ignoring ghost zones.
    for (int itheta = NGHOSTS; itheta < dst_Nxx1 + NGHOSTS; itheta++) { // Iterate over theta indices, ignoring ghost zones.
      // Perform interpolation only at radial index i = NGHOSTS.
      const BHA_REAL r_dst = dst_radii_aka_src_h_gf[DST_IDX3(NGHOSTS, itheta, iphi)];
      // Calculate the central index for the stencil in the radial direction.
      const int idx_center0 = (int)((r_dst - xxmin_incl_ghosts0) * src_invdxx0 + 0.5);

      // Determine the base index for the interpolation stencil.
      const int base_idx_x0 = idx_center0 - NinterpGHOSTS;

      {
        // Ensure the stencil is within valid grid bounds.
        if ((idx_center0 - NinterpGHOSTS < 0) || (idx_center0 + NinterpGHOSTS >= src_Nxx_plus_2NGHOSTS0)) {
#ifdef DEBUG
          // Provide detailed error messages in debug mode for easier troubleshooting.
          fprintf(stderr, "ERROR: Interpolation stencil exceeds grid boundaries for r_dst = %.6f. itheta, iphi = %d %d.\n", r_dst, itheta, iphi);
          fprintf(stderr, "Grid bounds along radial direction: [0, %d], idx_center = %d, stencil indices: [%d, %d]\n", src_Nxx_plus_2NGHOSTS0 - 1,
                  idx_center0, idx_center0 - NinterpGHOSTS, idx_center0 + NinterpGHOSTS);
          fprintf(stderr, "Ensure that the destination point is within grid bounds or adjust the interpolation stencil.\n");
#endif // DEBUG
#pragma omp critical
          {
            error_flag = INTERP1D_HORIZON_TOO_LARGE; // Set error flag if stencil is out of bounds.
            if (idx_center0 - NinterpGHOSTS < 0)
              error_flag = INTERP1D_HORIZON_TOO_SMALL; // Adjust error flag if stencil is too small.
          }
          continue; // Skip further work for this iteration.
        } // END IF: stencil in bounds

#ifdef DEBUG
        // Verify that the central index is the closest grid point to the destination radius.
        if (fabs(interp_src_r_theta_phi[0][idx_center0] - r_dst) > commondata->interp_src_dxx0 * 0.5) {
          fprintf(stderr, "ERROR: Radial center index too far from destination point!\n");
        } // END IF: central index is properly centered
#endif // DEBUG
      } // END BLOCK: sanity checks for stencil bounds and center index

      // Step 1: Precompute all differences between destination radius and source grid points within the stencil.
      BHA_REAL diffs_x0[INTERP_ORDER];
      compute_diffs_xi(INTERP_ORDER, r_dst, &interp_src_r_theta_phi[0][base_idx_x0], diffs_x0);

      // Step 2: Compute Lagrange basis coefficients for the radial direction.
      BHA_REAL lagrange_basis_coeffs_x0[INTERP_ORDER];
      compute_lagrange_basis_coeffs_xi(INTERP_ORDER, inv_denom, diffs_x0, lagrange_basis_coeffs_x0);

      // Step 3: Perform the 1D Lagrange interpolation along the radial direction.
      for (int gf = 0; gf < NUM_INTERP_SRC_GFS; gf++) {
        dst_interp_gfs[DST_IDX4(gf, NGHOSTS, itheta, iphi)] =
            sum_lagrange_x0_simd(INTERP_ORDER, &interp_src_gfs[SRC_IDX4(gf, base_idx_x0, itheta, iphi)], lagrange_basis_coeffs_x0) *
            src_invdxx0_INTERP_ORDERm1;
      } // END LOOP: for which_gf over grid functions
    } // END LOOP: for itheta over theta
  } // END LOOP: for iphi over phi

  return error_flag; // Return the status of the interpolation process.
} // END FUNCTION: bah_interpolation_1d_radial_spokes_on_3d_src_grid

#pragma GCC reset_options // Reset compiler optimizations after the function

#ifdef STANDALONE

#include <omp.h>

#define RADIAL_EXTENT 10.0
#define NUM_RESOLUTIONS 3

/**
 * Initializes the coordinate arrays for the source grid.
 *
 * Allocates and sets up the coordinate values for r, theta, and phi based on the number of ghost zones and grid dimensions.
 *
 * @param N_r                        Number of grid points in the radial direction.
 * @param N_theta                    Number of grid points in the theta direction.
 * @param N_phi                      Number of grid points in the phi direction.
 * @param src_r_theta_phi            Array to store coordinate values for r, theta, and phi.
 * @param src_dxx                    Array to store grid spacings in r, theta, and phi directions.
 * @param src_Nxx_plus_2NGHOSTS      Array containing the total number of grid points in each direction, including ghost zones.
 */
int initialize_coordinates(const int N_r, const int N_theta, const int N_phi, BHA_REAL *src_r_theta_phi[3], BHA_REAL src_dxx[3],
                           const int src_Nxx_plus_2NGHOSTS[3]) {
  src_dxx[0] = RADIAL_EXTENT / N_r;
  src_dxx[1] = M_PI / N_theta;
  src_dxx[2] = 2.0 * M_PI / N_phi;
  src_r_theta_phi[0] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * src_Nxx_plus_2NGHOSTS[0]);
  src_r_theta_phi[1] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * src_Nxx_plus_2NGHOSTS[1]);
  src_r_theta_phi[2] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * src_Nxx_plus_2NGHOSTS[2]);
  if (!src_r_theta_phi[0] || !src_r_theta_phi[1] || !src_r_theta_phi[2]) {
    free(src_r_theta_phi[0]);
    free(src_r_theta_phi[1]);
    free(src_r_theta_phi[2]);
    return -1;
  } // END IF: memory allocation issue
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS[0]; i++)
    src_r_theta_phi[0][i] = (i - NinterpGHOSTS) * src_dxx[0];
  for (int j = 0; j < src_Nxx_plus_2NGHOSTS[1]; j++)
    src_r_theta_phi[1][j] = 0.0 + ((BHA_REAL)(j - NGHOSTS) + 0.5) * src_dxx[1];
  for (int k = 0; k < src_Nxx_plus_2NGHOSTS[2]; k++)
    src_r_theta_phi[2][k] = -M_PI + ((BHA_REAL)(k - NGHOSTS) + 0.5) * src_dxx[2];
  return 0;
} // END FUNCTION: initialize_coordinates
/**
 * Initializes the source grid function using a provided analytic function.
 *
 * Evaluates the analytic function at each grid point and stores the result in the source grid function array.
 *
 * @param src_Nxx_plus_2NGHOSTS  Array containing the total number of grid points in each direction, including ghost zones.
 * @param src_r_theta_phi        Array containing coordinate values for r, theta, and phi.
 * @param src_gf                 Array to store the initialized source grid function values.
 */
int initialize_src_gf(const int src_Nxx_plus_2NGHOSTS[3], BHA_REAL *src_r_theta_phi[3], BHA_REAL *src_gf) {
  const int src_Nxx_plus_2NGHOSTS0 = src_Nxx_plus_2NGHOSTS[0], src_Nxx_plus_2NGHOSTS1 = src_Nxx_plus_2NGHOSTS[1],
            src_Nxx_plus_2NGHOSTS2 = src_Nxx_plus_2NGHOSTS[2];
  BHA_REAL *r_func_values = (BHA_REAL *)malloc(sizeof(BHA_REAL) * src_Nxx_plus_2NGHOSTS0);
  if (!r_func_values) {
    fprintf(stderr, "malloc failed for r_func_values.\n");
    return -1;
  } // END IF: memory allocation issue
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS0; i++) {
    r_func_values[i] = sin(src_r_theta_phi[0][i]) * src_r_theta_phi[0][i] * src_r_theta_phi[0][i];
  } // END LOOP: for i over precomputed radial function values
#pragma omp parallel for
  for (int gf = 0; gf < NUM_INTERP_SRC_GFS; gf++) {
    for (int k = 0; k < src_Nxx_plus_2NGHOSTS2; k++) {
      const BHA_REAL cosphi = cos(src_r_theta_phi[2][k]);
      for (int j = 0; j < src_Nxx_plus_2NGHOSTS1; j++) {
        const BHA_REAL sintheta = sin(src_r_theta_phi[1][j]);
        for (int i = 0; i < src_Nxx_plus_2NGHOSTS0; i++) {
          src_gf[SRC_IDX4(gf, i, j, k)] = r_func_values[i] * sintheta * cosphi;
        } // END LOOP: for i over radial grid
      } // END LOOP: for j over theta grid
    } // END LOOP: for k over phi grid
  } // END LOOP: for gf over source gridfunctions
  free(r_func_values);
  return 0;
} // END FUNCTION: initialize_src_gf

/**
 * Main function to execute the standalone 1D interpolation and perform convergence validation tests.
 *
 * Sets up different grid resolutions, performs interpolation,
 * computes the L2 norm of interpolation errors, and evaluates the observed order of convergence.
 *
 * @return EXIT_SUCCESS on successful execution.
 *         EXIT_FAILURE or specific error codes if memory allocation or interpolation fails.
 */
int main() {
  int return_code = EXIT_SUCCESS;
  BHA_REAL *src_r_theta_phi[3] = {NULL, NULL, NULL};
  BHA_REAL *src_gf = NULL, *dst_radii = NULL, *dst_interp_gfs = NULL;
  int N_r_arr[NUM_RESOLUTIONS] = {16, 32, 64};
  BHA_REAL h_arr[NUM_RESOLUTIONS], error_L2_norm[NUM_RESOLUTIONS];
  const int N_theta = 32, N_phi = 64;

  for (int res = 0; res < NUM_RESOLUTIONS; res++) {
    int N_r = N_r_arr[res];
    h_arr[res] = (N_r > 0) ? (RADIAL_EXTENT / N_r) : 0.0;
    int src_Nxx_plus_2NGHOSTS[3] = {N_r + 2 * NinterpGHOSTS, N_theta + 2 * NGHOSTS, N_phi + 2 * NGHOSTS};
    BHA_REAL src_dxx[3];
    if (initialize_coordinates(N_r, N_theta, N_phi, src_r_theta_phi, src_dxx, src_Nxx_plus_2NGHOSTS) != 0) {
      return_code = EXIT_FAILURE;
      goto cleanup;
    }

    params_struct params = {.Nxx0 = 1,
                            .Nxx1 = N_theta,
                            .Nxx2 = N_phi,
                            .Nxx_plus_2NGHOSTS0 = 1 + 2 * NGHOSTS,
                            .Nxx_plus_2NGHOSTS1 = N_theta + 2 * NGHOSTS,
                            .Nxx_plus_2NGHOSTS2 = N_phi + 2 * NGHOSTS};
    const int dst_Nxx_plus_2NGHOSTS0 = params.Nxx_plus_2NGHOSTS0;
    const int dst_Nxx_plus_2NGHOSTS1 = params.Nxx_plus_2NGHOSTS1;
    const int dst_Nxx_plus_2NGHOSTS2 = params.Nxx_plus_2NGHOSTS2;
    const size_t total_src_pts = (size_t)src_Nxx_plus_2NGHOSTS[0] * src_Nxx_plus_2NGHOSTS[1] * src_Nxx_plus_2NGHOSTS[2];
    src_gf = (BHA_REAL *)malloc(sizeof(BHA_REAL) * NUM_INTERP_SRC_GFS * total_src_pts);
    if (!src_gf) {
      fprintf(stderr, "malloc failed for src_gf.\n");
      return_code = EXIT_FAILURE;
      goto cleanup;
    }
    if (initialize_src_gf(src_Nxx_plus_2NGHOSTS, src_r_theta_phi, src_gf) != 0) {
      return_code = EXIT_FAILURE;
      goto cleanup;
    }

    // Use designated initializers for commondata_struct to avoid ambiguity and errors.
    commondata_struct commondata = {.interp_src_Nxx0 = N_r,
                                    .interp_src_Nxx_plus_2NGHOSTS0 = src_Nxx_plus_2NGHOSTS[0],
                                    .interp_src_Nxx_plus_2NGHOSTS1 = src_Nxx_plus_2NGHOSTS[1],
                                    .interp_src_Nxx_plus_2NGHOSTS2 = src_Nxx_plus_2NGHOSTS[2],
                                    .interp_src_invdxx0 = 1.0 / src_dxx[0],
                                    .interp_src_dxx0 = src_dxx[0],
                                    .interp_src_gfs = src_gf};
    for (int i = 0; i < 3; i++)
      commondata.interp_src_r_theta_phi[i] = src_r_theta_phi[i];

    const size_t total_dst_pts = (size_t)params.Nxx_plus_2NGHOSTS0 * params.Nxx_plus_2NGHOSTS1 * params.Nxx_plus_2NGHOSTS2;
    dst_radii = (BHA_REAL *)malloc(sizeof(BHA_REAL) * total_dst_pts);
    dst_interp_gfs = (BHA_REAL *)malloc(sizeof(BHA_REAL) * NUM_INTERP_SRC_GFS * total_dst_pts);
    if (!dst_radii || !dst_interp_gfs) {
      fprintf(stderr, "malloc failed for destination arrays.\n");
      return_code = EXIT_FAILURE;
      goto cleanup;
    }

    const BHA_REAL r_min_safe = src_r_theta_phi[0][NinterpGHOSTS] + 1e-6;
    const BHA_REAL r_max_safe = src_r_theta_phi[0][src_Nxx_plus_2NGHOSTS[0] - NinterpGHOSTS - 1] - 1e-6;
    srand(42 + res);
    for (int j = NGHOSTS; j < N_theta + NGHOSTS; j++) {
      for (int k = NGHOSTS; k < N_phi + NGHOSTS; k++) {
        // Generate a random value between r_min_safe and r_max_safe
        dst_radii[DST_IDX3(NGHOSTS, j, k)] = r_min_safe + ((BHA_REAL)rand() / RAND_MAX) * (r_max_safe - r_min_safe);
      }
    }

#ifdef _OPENMP
    double start_time = omp_get_wtime();
#endif
    int error_code = BHAHAHA_SUCCESS;
    const int NUM_TIMES = 90000;
    for (int n = 0; n < NUM_TIMES; n++) {
      error_code = bah_interpolation_1d_radial_spokes_on_3d_src_grid(&params, &commondata, dst_radii, dst_interp_gfs);
      if (error_code != BHAHAHA_SUCCESS)
        break;
    }
#ifdef _OPENMP
    double end_time = omp_get_wtime();
#endif

    printf("--- Benchmarking for Resolution N_r = %d ---\n", N_r);
#ifdef _OPENMP
    const long long num_interps = (long long)N_theta * N_phi * NUM_TIMES;
    printf("Interpolated %lld points in %.4f seconds.\n", num_interps, end_time - start_time);
    printf("Performance: %.4f million points per second.\n", (double)num_interps / (end_time - start_time) / 1e6);
#endif

    if (error_code != BHAHAHA_SUCCESS) {
      fprintf(stderr, "Interpolation failed with error code: %d\n", error_code);
      return_code = error_code;
      goto cleanup;
    }

    BHA_REAL error_sum = 0.0;
#pragma omp parallel for reduction(+ : error_sum)
    for (int k = NGHOSTS; k < N_phi + NGHOSTS; k++) {
      const BHA_REAL phi = src_r_theta_phi[2][k];
      const BHA_REAL cosphi = cos(phi);
      for (int j = NGHOSTS; j < N_theta + NGHOSTS; j++) {
        const BHA_REAL theta = src_r_theta_phi[1][j];
        const BHA_REAL sintheta = sin(theta);
        const BHA_REAL r_dst = dst_radii[DST_IDX3(NGHOSTS, j, k)];
        const BHA_REAL r_func = sin(r_dst) * r_dst * r_dst;

        BHA_REAL exact_value = r_func * sintheta * cosphi;
        // Only need to check the zeroth gridfunction; all gridfunctions are the same.
        BHA_REAL interp_value = dst_interp_gfs[DST_IDX4(0, NGHOSTS, j, k)];
        BHA_REAL error = interp_value - exact_value;
        error_sum += error * error;
      }
    }
    error_L2_norm[res] = sqrt(error_sum / (N_theta * N_phi));
    printf("Resolution %d: N_r = %d, h = %.5e, L2 error = %.5e\n", res, N_r, h_arr[res], error_L2_norm[res]);

    for (int i = 0; i < 3; i++) {
      free(src_r_theta_phi[i]);
      src_r_theta_phi[i] = NULL;
    }
    free(src_gf);
    src_gf = NULL;
    free(dst_radii);
    dst_radii = NULL;
    free(dst_interp_gfs);
    dst_interp_gfs = NULL;
  }

  const int INTERP_ORDER = (2 * NinterpGHOSTS + 1);
  printf("\n--- Convergence Results ---\n");
  for (int res = 1; res < NUM_RESOLUTIONS; res++) {
    BHA_REAL observed_order = log(error_L2_norm[res - 1] / error_L2_norm[res]) / log(h_arr[res - 1] / h_arr[res]);
    printf("Observed order of convergence between res %d and %d: %.2f\n", res - 1, res, observed_order);
  }
  printf("Expected order of convergence: %d\n", INTERP_ORDER);

cleanup:
  if (return_code != EXIT_SUCCESS)
    printf("\nAn error occurred. Cleaning up...\n");
  else
    printf("\nProgram finished successfully. Cleaning up...\n");
  for (int i = 0; i < 3; i++)
    free(src_r_theta_phi[i]);
  free(src_gf);
  free(dst_radii);
  free(dst_interp_gfs);
  return return_code;
} // END FUNCTION: main

#endif // STANDALONE
