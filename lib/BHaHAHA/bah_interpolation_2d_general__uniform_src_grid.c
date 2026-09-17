#include "interpolation_lagrange_uniform.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#ifndef BHA_REAL
#define BHA_REAL double
#endif
#define DEBUG

#ifdef STANDALONE
// Remove the bah_ prefix if compiling as a standalone code, as this function goes by other names in other codes.
#define bah_interpolation_2d_general__uniform_src_grid interpolation_2d_general__uniform_src_grid
#endif

// In case this code is compiled as C++:
#ifdef __cplusplus
#ifndef restrict
#define restrict __restrict__
#endif
#endif
//===============================================
// BHaHAHA Error Handling
//===============================================
// Error codes, set in error_message.py
typedef enum {
  BHAHAHA_SUCCESS,
  FIND_HORIZON_GETTIMEOFDAY_BROKEN,
  FIND_HORIZON_MAX_ITERATIONS_EXCEEDED,
  FIND_HORIZON_HORIZON_TOO_SMALL,
  BCSTRUCT_EIGENCOORD_FAILURE,
  BCSTRUCT_SET_PARITY_ERROR,
  INITIAL_DATA_MALLOC_ERROR,
  NUMGRID_EXTERN_MALLOC_ERROR_GFS,
  NUMGRID_EXTERN_MALLOC_ERROR_RTHETAPHI,
  NUMGRID_INTERP_MALLOC_ERROR_GFS,
  NUMGRID_INTERP_MALLOC_ERROR_RTHETAPHI,
  INTERP1D_NULL_PTRS,
  INTERP1D_INTERP_ORDER_GT_NXX_PLUS_2NINTERPGHOSTS0,
  INTERP1D_HORIZON_TOO_LARGE,
  INTERP1D_HORIZON_TOO_SMALL,
  INTERP2D_EXT_TO_INTERPSRC_NULL_PTRS,
  INTERP2D_EXT_TO_INTERPSRC_INTERP_ORDER_GT_NXX_PLUS_2NINTERPGHOSTS12,
  INTERP2D_EXT_TO_INTERPSRC_HORIZON_OUT_OF_BOUNDS,
  INTERP2D_GENERAL_NULL_PTRS,
  INTERP2D_GENERAL_INTERP_ORDER_GT_NXX_PLUS_2NGHOSTS12,
  INTERP2D_GENERAL_HORIZON_OUT_OF_BOUNDS,
  DIAG_PROPER_CIRCUM_MALLOC_ERROR,
} bhahaha_error_codes;
#pragma GCC optimize("unroll-loops")

/**
 * Performs 2D Lagrange interpolation on a uniform source grid.
 *
 * This function interpolates a scalar grid function from a uniform source grid defined in theta and phi
 * to a set of arbitrary destination points using Lagrange interpolation of a specified order.
 *
 * @param n_interp_ghosts          Number of ghost zones around each source point. Determines the interpolation order as 2 * n_interp_ghosts + 1.
 * @param src_dxx1                 Grid spacing in the theta direction on the source grid.
 * @param src_dxx2                 Grid spacing in the phi direction on the source grid.
 * @param src_Nxx_plus_2NGHOSTS1  Total number of grid points in the theta direction, including ghost zones.
 * @param src_Nxx_plus_2NGHOSTS2  Total number of grid points in the phi direction, including ghost zones.
 * @param[in,out] src_r_theta_phi  Arrays containing coordinate values for r, theta, and phi on the source grid.
 * @param[in] src_gf               Pointer to the source grid function data, organized as a flattened 2D array.
 * @param num_dst_pts              Number of destination points where interpolation is performed.
 * @param[in] dst_pts              Array of destination points' coordinates, each consisting of (theta, phi).
 * @param[out] dst_data            Output array to store interpolated values at each destination point.
 *
 * @return                         BHAHAHA_SUCCESS on successful interpolation.
 *                                  Appropriate error code if an error is encountered.
 *
 * @note
 * - Assumes that the source and destination grids are uniform in theta and phi directions.
 * - Ensures that destination points lie within the bounds of the source grid to prevent memory access violations.
 */
int bah_interpolation_2d_general__uniform_src_grid(const int n_interp_ghosts, const BHA_REAL src_dxx1, const BHA_REAL src_dxx2,
                                                   const int src_Nxx_plus_2NGHOSTS1, const int src_Nxx_plus_2NGHOSTS2,
                                                   BHA_REAL *restrict src_r_theta_phi[3], const BHA_REAL *restrict src_gf, const int num_dst_pts,
                                                   const BHA_REAL dst_pts[][2], BHA_REAL *restrict dst_data) {
  // Define the interpolation order based on the number of ghost zones.
  const int INTERP_ORDER = (2 * n_interp_ghosts + 1); // Interpolation order corresponds to the number of points in the stencil per dimension.

  // Calculate inverse grid spacings for efficient index calculations.
  const BHA_REAL src_invdxx1 = 1.0 / src_dxx1;
  const BHA_REAL src_invdxx2 = 1.0 / src_dxx2;

  // Compute normalization factor to scale the interpolated result appropriately.
  const BHA_REAL src_invdxx12_INTERP_ORDERm1 = pow(src_dxx1 * src_dxx2, -(INTERP_ORDER - 1));

  // Validate input pointers to prevent segmentation faults.
  if (src_r_theta_phi[1] == NULL || src_r_theta_phi[2] == NULL || src_gf == NULL || dst_data == NULL)
    return INTERP2D_GENERAL_NULL_PTRS; // Exit if any required pointer is NULL.

  // Ensure that the interpolation order does not exceed the grid dimensions in either direction.
  if (INTERP_ORDER > src_Nxx_plus_2NGHOSTS1 || INTERP_ORDER > src_Nxx_plus_2NGHOSTS2)
    return INTERP2D_GENERAL_INTERP_ORDER_GT_NXX_PLUS_2NGHOSTS12; // Exit if interpolation order is too high.

  // Precompute inverse denominators for Lagrange interpolation coefficients to reduce redundant calculations.
  BHA_REAL inv_denom[INTERP_ORDER];
  compute_inv_denom(INTERP_ORDER, inv_denom);

  // Define the minimum coordinate values including ghost zones for both theta and phi.
  const BHA_REAL xxmin_incl_ghosts1 = src_r_theta_phi[1][0];
  const BHA_REAL xxmin_incl_ghosts2 = src_r_theta_phi[2][0];

  // Initialize the error flag to track any interpolation issues.
  int error_flag = BHAHAHA_SUCCESS;

#pragma omp parallel for
  for (int dst_pt = 0; dst_pt < num_dst_pts; dst_pt++) {
    const BHA_REAL theta_dst = dst_pts[dst_pt][0]; // Destination point's theta coordinate.
    const BHA_REAL phi_dst = dst_pts[dst_pt][1];   // Destination point's phi coordinate.

    // Determine the central grid indices in theta and phi for the interpolation stencil.
    const int idx_center_th = (int)((theta_dst - xxmin_incl_ghosts1) * src_invdxx1 + 0.5);
    const int idx_center_ph = (int)((phi_dst - xxmin_incl_ghosts2) * src_invdxx2 + 0.5);

    {
      // Verify that the interpolation stencil does not exceed grid boundaries.
      if ((idx_center_th - n_interp_ghosts < 0) || (idx_center_th + n_interp_ghosts >= src_Nxx_plus_2NGHOSTS1) ||
          (idx_center_ph - n_interp_ghosts < 0) || (idx_center_ph + n_interp_ghosts >= src_Nxx_plus_2NGHOSTS2)) {
#ifdef DEBUG
        // Provide detailed error messages in debug mode for easier troubleshooting.
        fprintf(stderr, "ERROR: Interpolation stencil exceeds grid boundaries. %d %d %d %d\n", (idx_center_th - n_interp_ghosts < 0),
                (idx_center_th + n_interp_ghosts >= src_Nxx_plus_2NGHOSTS1), (idx_center_ph - n_interp_ghosts < 0),
                (idx_center_ph + n_interp_ghosts >= src_Nxx_plus_2NGHOSTS2));
        fprintf(stderr, "(th_dst, ph_dst) = (%.6f, %.6f). (dst_pt) = %d.\n", theta_dst, phi_dst, dst_pt);
        fprintf(stderr, "Grid bounds along theta direction: [0, %d], stencil indices: [%d, %d]\n", src_Nxx_plus_2NGHOSTS1 - 1,
                idx_center_th - n_interp_ghosts, idx_center_th + n_interp_ghosts);
        fprintf(stderr, "Grid bounds along phi direction: [0, %d], stencil indices: [%d, %d]\n", src_Nxx_plus_2NGHOSTS2 - 1,
                idx_center_ph - n_interp_ghosts, idx_center_ph + n_interp_ghosts);
        fprintf(stderr, "Ensure that the destination point is within grid bounds or adjust the interpolation stencil.\n");
#endif // DEBUG
#pragma omp critical
        {
          error_flag = INTERP2D_GENERAL_HORIZON_OUT_OF_BOUNDS; // Set error flag if stencil is out of bounds.
        }
        continue; // Skip interpolation for this destination point to prevent invalid memory access.
      } // END IF: theta/phi stencil exceeded source-grid bounds

      // Additional sanity checks to ensure central index is correctly positioned.
#ifdef DEBUG
      const BHA_REAL TOLERANCE = 1e-13; // Tolerance to account for floating-point precision.
      if (fabs(src_r_theta_phi[1][idx_center_th] - theta_dst) > src_dxx1 * (0.5 + TOLERANCE)) {
        fprintf(stderr, "ERROR: theta center index too far from destination point! %.15e > %.15e\n",
                fabs(src_r_theta_phi[1][idx_center_th] - theta_dst), src_dxx1 * (0.5 + TOLERANCE));
      }
      if (fabs(src_r_theta_phi[2][idx_center_ph] - phi_dst) > src_dxx2 * (0.5 + TOLERANCE)) {
        fprintf(stderr, "ERROR: phi center index too far from destination point! %.15e > %.15e\n", fabs(src_r_theta_phi[2][idx_center_ph] - phi_dst),
                src_dxx2 * (0.5 + TOLERANCE));
      } // END IF: Central index is properly centered
#endif // DEBUG
    } // END BLOCK: theta/phi stencil bounds and center-index sanity checks

    // Calculate the starting indices for the interpolation stencil in theta and phi directions.
    const int base_idx_th = idx_center_th - n_interp_ghosts;
    const int base_idx_ph = idx_center_ph - n_interp_ghosts;

    // Step 1: Precompute differences between destination theta and source grid theta points within the stencil.
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

    // Define a macro to calculate the flattened index for accessing the source grid function.
#define SRC_IDX2(j, k) ((j) + src_Nxx_plus_2NGHOSTS1 * (k))
      // Step 3: Perform the 1D Lagrange interpolation along the radial direction.
    BHA_REAL sum = 0.0;

    for (int iph = 0; iph < INTERP_ORDER; iph++) {
      const int idx_ph = base_idx_ph + iph;
      const int base_offset = base_idx_th + src_Nxx_plus_2NGHOSTS1 * idx_ph;

      sum += sum_lagrange_x0_simd(INTERP_ORDER, &src_gf[base_offset], &coeff_2d[iph][0]);
    } // END LOOP: for iph over phi direction

    // Store the interpolated value for this grid function and destination point.
    dst_data[dst_pt] = sum * src_invdxx12_INTERP_ORDERm1;

  } // END LOOP: for dst_pt over destination points

  return error_flag; // Return the status of the interpolation process.
} // END FUNCTION: bah_interpolation_2d_general__uniform_src_grid

#pragma GCC reset_options // Reset compiler optimizations after the function.

#ifdef STANDALONE

#include <omp.h>

#define NUM_INTERP_GFS 2
#define NUM_RESOLUTIONS 3
#define NUM_DST_PTS 40000000

// Analytic functions.
static inline BHA_REAL analytic_function1(BHA_REAL x0, BHA_REAL x1) { return sin(x0) * cos(x1); }
static inline BHA_REAL analytic_function2(BHA_REAL x0, BHA_REAL x1) { return cos(x0) * sin(x1); }

/**
 * Initializes the 1D coordinate arrays for a 2D uniform source grid.
 *
 * This function calculates the grid spacing (dx) for each dimension and allocates memory for
 * and populates the 1D coordinate arrays. The coordinate arrays include ghost zones. The first
 * coordinate pointer is unused and set to NULL to match the expected 3-pointer array format.
 *
 * @param n_interp_ghosts The number of ghost zones on each side for interpolation.
 * @param N_x0 The number of interior grid points in the x0-dimension.
 * @param N_x1 The number of interior grid points in the x1-dimension.
 * @param[out] src_x0x1 An array of 3 pointers. The function allocates memory for indices 1 and 2. Index 0 is unused.
 * @param[out] src_dxx0 Pointer to store the calculated grid spacing in the x0-dimension.
 * @param[out] src_dxx1 Pointer to store the calculated grid spacing in the x1-dimension.
 * @param src_Nxx_plus_2NGHOSTS0 The total number of points in the x0-dimension, including ghost zones.
 * @param src_Nxx_plus_2NGHOSTS1 The total number of points in the x1-dimension, including ghost zones.
 * @return 0 on success, -1 on memory allocation failure.
 */
int initialize_coordinates(const int n_interp_ghosts, const int N_x0, const int N_x1, BHA_REAL *src_x0x1[3], BHA_REAL *src_dxx0, BHA_REAL *src_dxx1,
                           const int src_Nxx_plus_2NGHOSTS0, const int src_Nxx_plus_2NGHOSTS1) {
  *src_dxx0 = (M_PI) / N_x0;
  *src_dxx1 = (2.0 * M_PI) / N_x1;
  // Index 0 is unused, mimicking the original 2D code's r-theta-phi structure.
  src_x0x1[0] = NULL;
  src_x0x1[1] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * src_Nxx_plus_2NGHOSTS0);
  src_x0x1[2] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * src_Nxx_plus_2NGHOSTS1);
  if (!src_x0x1[1] || !src_x0x1[2]) {
    free(src_x0x1[1]);
    free(src_x0x1[2]);
    src_x0x1[1] = src_x0x1[2] = NULL;
    return -1;
  } // END IF: coordinate-array allocation failed
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS0; i++)
    src_x0x1[1][i] = (i - n_interp_ghosts) * (*src_dxx0);
  // END LOOP: for i over x0 coordinates
  for (int i = 0; i < src_Nxx_plus_2NGHOSTS1; i++)
    src_x0x1[2][i] = (i - n_interp_ghosts) * (*src_dxx1);
  // END LOOP: for i over x1 coordinates
  return 0;
} // END FUNCTION: initialize_coordinates

/**
 * Populates a 2D source grid function with values from a given analytic function.
 *
 * This function iterates through all points of a 2D grid (including ghost zones)
 * and computes the value of the grid function at each point using the provided
 * analytic function pointer.
 *
 * @param src_Nxx_plus_2NGHOSTS0 The total number of points in the x0-dimension.
 * @param src_Nxx_plus_2NGHOSTS1 The total number of points in the x1-dimension.
 * @param src_x0x1 The pre-initialized 1D coordinate arrays (as an array of 3 pointers).
 * @param[out] src_gf The 2D source grid function data array to be populated.
 * @param func A function pointer to the analytic function used to compute the values.
 */
void initialize_src_gf(const int src_Nxx_plus_2NGHOSTS0, const int src_Nxx_plus_2NGHOSTS1, BHA_REAL *src_x0x1[3], BHA_REAL *src_gf,
                       BHA_REAL (*func)(BHA_REAL, BHA_REAL)) {
  for (int i1 = 0; i1 < src_Nxx_plus_2NGHOSTS1; i1++) {
    for (int i0 = 0; i0 < src_Nxx_plus_2NGHOSTS0; i0++) {
      const int idx = i0 + src_Nxx_plus_2NGHOSTS0 * i1;
      src_gf[idx] = func(src_x0x1[1][i0], src_x0x1[2][i1]);
    } // END LOOP: for i0 over x0 source-grid points
  } // END LOOP: for i1 over x1 source-grid points
} // END FUNCTION: initialize_src_gf

/**
 * Main driver for testing 2D Lagrange interpolation.
 *
 * This program tests a 2D Lagrange interpolation routine by performing the following steps:
 * 1. Sets up source grids at multiple resolutions.
 * 2. Populates the source grids with data from known analytic functions.
 * 3. Generates a set of random destination points within the source grid domain.
 * 4. Calculates the exact function values at these destination points.
 * 5. Calls the interpolation routine to compute interpolated values at the destination points.
 * 6. Measures the L2 norm of the error between the interpolated and exact values.
 * 7. Calculates and prints the observed order of convergence to verify the accuracy of the interpolator.
 * 8. Prints performance benchmarks.
 *
 * @return EXIT_SUCCESS on successful completion, EXIT_FAILURE otherwise.
 */
int main() {
  int return_code = EXIT_SUCCESS;
  BHA_REAL(*dst_pts)[2] = NULL;
  BHA_REAL *f_exact[NUM_INTERP_GFS] = {NULL};
  // Changed to an array of 3 pointers to match the interpolator's function signature.
  BHA_REAL *src_x0x1[3] = {NULL, NULL, NULL};
  BHA_REAL *src_gf[NUM_INTERP_GFS] = {NULL};
  BHA_REAL *dst_data[NUM_INTERP_GFS] = {NULL};

  const int n_interp_ghosts = 3;
  const int INTERP_ORDER = 2 * n_interp_ghosts + 1;

  int N_x0_arr[NUM_RESOLUTIONS] = {16, 32, 64};
  int N_x1_arr[NUM_RESOLUTIONS] = {16, 32, 64};
  BHA_REAL h_arr[NUM_RESOLUTIONS];
  BHA_REAL error_L2_norm[NUM_INTERP_GFS][NUM_RESOLUTIONS];

  dst_pts = (BHA_REAL(*)[2])malloc(sizeof(BHA_REAL) * NUM_DST_PTS * 2);
  if (!dst_pts) {
    fprintf(stderr, "malloc failed for dst_pts.\n");
    return_code = EXIT_FAILURE;
    goto cleanup;
  } // END IF: destination-point allocation failed

  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    f_exact[gf] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * NUM_DST_PTS);
    if (!f_exact[gf]) {
      fprintf(stderr, "malloc failed for f_exact.\n");
      return_code = EXIT_FAILURE;
      goto cleanup;
    } // END IF: exact-value allocation failed
  } // END LOOP: for gf over exact function value arrays

  for (int res = 0; res < NUM_RESOLUTIONS; res++) {
    // --- Main Loop Over Resolutions ---
    int N_x0 = N_x0_arr[res];
    int N_x1 = N_x1_arr[res];
    h_arr[res] = (N_x0 > 0) ? ((BHA_REAL)(M_PI) / N_x0) : 0.0;
    int src_Nxx_plus_2NGHOSTS0 = N_x0 + 2 * n_interp_ghosts;
    int src_Nxx_plus_2NGHOSTS1 = N_x1 + 2 * n_interp_ghosts;
    BHA_REAL src_dxx0_val, src_dxx1_val;
    if (initialize_coordinates(n_interp_ghosts, N_x0, N_x1, src_x0x1, &src_dxx0_val, &src_dxx1_val, src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1) !=
        0) {
      fprintf(stderr, "malloc failed for coordinates.\n");
      return_code = EXIT_FAILURE;
      goto cleanup;
    } // END IF: source-coordinate initialization failed
    BHA_REAL x0_min_safe = src_x0x1[1][n_interp_ghosts] + 1e-6;
    BHA_REAL x0_max_safe = src_x0x1[1][src_Nxx_plus_2NGHOSTS0 - n_interp_ghosts - 1] - 1e-6;
    BHA_REAL x1_min_safe = src_x0x1[2][n_interp_ghosts] + 1e-6;
    BHA_REAL x1_max_safe = src_x0x1[2][src_Nxx_plus_2NGHOSTS1 - n_interp_ghosts - 1] - 1e-6;
    srand(42 + res);
    for (int i = 0; i < NUM_DST_PTS; i++) {
      dst_pts[i][0] = x0_min_safe + ((BHA_REAL)rand() / RAND_MAX) * (x0_max_safe - x0_min_safe);
      dst_pts[i][1] = x1_min_safe + ((BHA_REAL)rand() / RAND_MAX) * (x1_max_safe - x1_min_safe);
      f_exact[0][i] = analytic_function1(dst_pts[i][0], dst_pts[i][1]);
      f_exact[1][i] = analytic_function2(dst_pts[i][0], dst_pts[i][1]);
    } // END LOOP: for dst_pt over destination points and exact values
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      const size_t size = (size_t)src_Nxx_plus_2NGHOSTS0 * src_Nxx_plus_2NGHOSTS1;
      src_gf[gf] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * size);
      if (!src_gf[gf]) {
        fprintf(stderr, "malloc failed for src_gf.\n");
        return_code = EXIT_FAILURE;
        goto cleanup;
      } // END IF: source-gridfunction allocation failed
      dst_data[gf] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * NUM_DST_PTS);
      if (!dst_data[gf]) {
        fprintf(stderr, "malloc failed for dst_data.\n");
        return_code = EXIT_FAILURE;
        goto cleanup;
      } // END IF: destination-data allocation failed
    } // END LOOP: for gf over source grid and destination data arrays
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_x0x1, src_gf[0], analytic_function1);
    initialize_src_gf(src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1, src_x0x1, src_gf[1], analytic_function2);

#ifdef _OPENMP
    double start_time = omp_get_wtime();
#endif
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      int error_code =
          bah_interpolation_2d_general__uniform_src_grid(n_interp_ghosts, src_dxx0_val, src_dxx1_val, src_Nxx_plus_2NGHOSTS0, src_Nxx_plus_2NGHOSTS1,
                                                         src_x0x1, src_gf[gf], NUM_DST_PTS, dst_pts, dst_data[gf]);
      if (error_code != BHAHAHA_SUCCESS) {
        fprintf(stderr, "Interpolation error code: %d for GF %d\n", error_code, gf + 1);
        return_code = error_code;
        goto cleanup;
      } // END IF: interpolation routine returned error
    } // END LOOP: for which_gf over grid functions for interpolation
#ifdef _OPENMP
    double elapsed_time = omp_get_wtime() - start_time;
#endif

    printf("\n--- Benchmarking for Resolution %d (%dx%d) ---\n", res, N_x0, N_x1);
#ifdef _OPENMP
    printf("Interpolated %d points for %d GFs in %.4f seconds.\n", NUM_DST_PTS, NUM_INTERP_GFS, elapsed_time);
    printf("Performance: %.4f million points per second.\n", (double)(NUM_DST_PTS * NUM_INTERP_GFS) / elapsed_time / 1e6);
#endif
    printf("--------------------------------------------------\n");

    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      BHA_REAL error_sum = 0.0;
      for (int i = 0; i < NUM_DST_PTS; i++) {
        BHA_REAL error = dst_data[gf][i] - f_exact[gf][i];
        error_sum += error * error;
      } // END LOOP: for dst_pt over squared error sum
      error_L2_norm[gf][res] = sqrt(error_sum / NUM_DST_PTS);
      printf("Resolution %d: N_x0=%d, h=%.5e, GF %d, L2 error=%.5e\n", res, N_x0, h_arr[res], gf + 1, error_L2_norm[gf][res]);
    } // END LOOP: for gf over L2 error norm
    for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
      free(src_gf[gf]);
      src_gf[gf] = NULL;
      free(dst_data[gf]);
      dst_data[gf] = NULL;
    } // END LOOP: for gf over data for current resolution
    // Changed loop to free coordinate arrays up to index 2
    for (int dim = 1; dim < 3; dim++) {
      free(src_x0x1[dim]);
      src_x0x1[dim] = NULL;
    } // END LOOP: for dim over coordinate arrays for current resolution
  } // END LOOP: for res over resolutions

  printf("\n--- Convergence Results ---\n");
  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    for (int res = 1; res < NUM_RESOLUTIONS; res++) {
      BHA_REAL observed_order = log(error_L2_norm[gf][res - 1] / error_L2_norm[gf][res]) / log(h_arr[res - 1] / h_arr[res]);
      printf("Observed order of convergence for GF %d between res %d and %d: %.2f\n", gf + 1, res - 1, res, observed_order);
    } // END LOOP: for res over observed convergence order
    printf("Expected order of convergence for GF %d: %d\n", gf + 1, INTERP_ORDER);
  } // END LOOP: for gf over convergence results

cleanup:
  if (return_code == EXIT_FAILURE)
    printf("\nAn error occurred. Cleaning up...\n");
  else
    printf("\nProgram finished successfully. Cleaning up...\n");
  // END IF: print final status message
  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    free(src_gf[gf]);
    free(dst_data[gf]);
  } // END LOOP: for gf over src_gf and dst_data cleanup
  // Changed loop to clean up all three pointers
  for (int dim = 0; dim < 3; dim++) {
    free(src_x0x1[dim]);
  } // END LOOP: for dim over src_x0x1 cleanup
  for (int gf = 0; gf < NUM_INTERP_GFS; gf++) {
    free(f_exact[gf]);
  } // END LOOP: for gf over f_exact cleanup
  free(dst_pts);
  return return_code;

} // END FUNCTION: main

#endif // STANDALONE
