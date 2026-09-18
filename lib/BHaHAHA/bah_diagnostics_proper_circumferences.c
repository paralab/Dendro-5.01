#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Complete elliptic integrals E(m) and K(m) via 8th-order midpoint quadrature.
 *
 * Computes the complete elliptic integrals in the **parameter** (Legendre) form
 *   E(m) = ∫₀^{π/2} sqrt(1 − m sin²θ) dθ,
 *   K(m) = ∫₀^{π/2} dθ / sqrt(1 − m sin²θ).
 *
 * The implementation uses a fixed 8th-order midpoint rule with periodic weights
 * over [0, π/2]. The sample count is fixed at 128 (power of two) for high accuracy
 * and good vectorization.
 *
 * @param m  Legendre **parameter** (a.k.a. k²). Real results require m ≤ 1.
 *           Negative m is supported (e.g. m ∈ [−1, 0] in this code path).
 *           @warning This is the parameter m, **not** the modulus k. If you have a
 *           modulus k, pass m = k*k.
 * @param[out] E  On return, E(m).
 * @param[out] K  On return, K(m) (diverges as m → 1⁻).
 *
 * @pre E and K are non-null.
 * @pre The internal weight generator expects the sample count to be compatible
 *      with the 8th-order periodic stencil (here fixed to 128).
 */
static void elliptic_E_and_K_integrals(const BHA_REAL k, BHA_REAL *restrict E, BHA_REAL *restrict K) {
  static const int N_sample_pts = 128; // Number of sample points for integration. Chosen for high precision.
  const BHA_REAL *restrict weights;        // Precomputed integration weights for accuracy in the midpoint method.
  int weight_stencil_size;             // Size of the weight stencil to correctly cycle through the weights.

  // Retrieve the integration weights based on the number of sample points (divisible by 8 -> 8th order).
  bah_diagnostics_integration_weights(N_sample_pts, N_sample_pts, &weights, &weight_stencil_size);

  const BHA_REAL a = 0.0;                            // Lower limit of integration (0 radians).
  const BHA_REAL b = M_PI / 2.0;                     // Upper limit of integration (pi/2 radians).
  const BHA_REAL h = (b - a) / ((BHA_REAL)N_sample_pts); // Step size for each subinterval based on sample points.

  BHA_REAL sum_E = 0.0; // Accumulator for the elliptic integral of the second kind E(k).
  BHA_REAL sum_K = 0.0; // Accumulator for the elliptic integral of the first kind K(k).

  // Parallelized loop to compute both integrals E(k) and K(k) using OpenMP.
  // #pragma omp parallel for reduction(+ : sum_E, sum_K) <- thread creation/destruction likely -> slower code here
  for (int i = 0; i < N_sample_pts; i++) {
    const BHA_REAL theta = a + ((BHA_REAL)i + 0.5) * h; // Compute the midpoint for the current subinterval.
    // Compute sin(theta). For optimization, we could use a lookup table if N_sample_pts is constant.
    const BHA_REAL sintheta = sin(theta);
    // Compute the integrands for the elliptic integrals of the second & first kinds (E(k) & K(k), respectively) at this midpoint.
    const BHA_REAL elliptic_E_integrand = sqrt(1.0 - k * sintheta * sintheta);
    const BHA_REAL elliptic_K_integrand = 1.0 / elliptic_E_integrand;
    // Update the running sums for E(k) & K(k), applying the corresponding weight.
    sum_E += weights[i % weight_stencil_size] * elliptic_E_integrand;
    sum_K += weights[i % weight_stencil_size] * elliptic_K_integrand;
  } // END LOOP: for i over sample points to compute both integrals

  // Multiply by the step size to complete the integration and store the results.
  *E = sum_E * h; // Elliptic integral of the second kind.
  *K = sum_K * h; // Elliptic integral of the first kind.
} // END FUNCTION: elliptic_E_and_K_integrals

/**
 * Estimates the spin parameter magnitude for equilibrium black holes based on the circumference ratio C_r.
 *
 * Rationale & safeguards:
 *   * Start from an analytic approximation (Alcubierre et al., Eq. 5.3) to land near the root basin.
 *   * Refine with Newton–Raphson using elliptic integrals (Eq. 5.2), but clamp any out-of-range iterates
 *     into [0,1] and terminate if we step negative; this avoids excursions where the modulus or square roots
 *     would be undefined or numerically fragile.
 *
 * @param C_r The circumference ratio parameter used to estimate the spin.
 * @return    The estimated spin parameter. Returns -10.0 if C_r is out of valid bounds or if convergence fails.
 *
 */
static BHA_REAL compute_spin(const BHA_REAL C_r) {
  // Validate the input parameter. Return an error code if C_r exceeds the valid range.
  if (C_r > 1)
    return -10.0;

  // Turns out, this is a conservative spin estimate.
  BHA_REAL spin = 0.9;

  // Refine the initial guess using an analytical approximation based on Eq. 5.3 of Alcubierre et al arXiv:gr-qc/0411149.
  const BHA_REAL spin_sq = 1 - (2.55 * C_r - 1.55) * (2.55 * C_r - 1.55);
  if (spin_sq >= 0 && spin_sq < 1)
    spin = sqrt(spin_sq);

  const BHA_REAL rel_tolerance = 1e-7; // Desired relative tolerance for convergence.
  BHA_REAL rel_diff = 1e10;            // Initialize relative difference to a large value.
  const int max_its = 20;          // Maximum number of iterations to prevent infinite loops.
  int it = 0;                      // Iteration counter.

  // Iteratively refine the spin estimate until the relative difference is within tolerance or max iterations are reached.
  while (rel_diff > rel_tolerance && it < max_its) {
    const BHA_REAL x = spin;
    BHA_REAL E, K;

    // Compute the elliptic integrals E and K based on the current spin estimate.
    elliptic_E_and_K_integrals(-((x * x) / pow(1 + sqrt(1 - (x * x)), 2)), &E, &K);

    // Next complete a Newton-Raphson iteration to improve the spin estimate.
    /*
     *  Original SymPy expression:
     *  "const BHA_REAL x_np1 = x - (-C_r + E*(sqrt(1 - x**2) + 1)/pi)/(-E*x/(pi*sqrt(1 - x**2)) - (E - K)*(-2*x**3/(sqrt(1 - x**2)*(sqrt(1 - x**2) +
     * 1)**3) - 2*x/(sqrt(1 - x**2) + 1)**2)*(sqrt(1 - x**2) + 1)**3/(2*pi*x**2))"
     */
    const BHA_REAL tmp2 = sqrt(1 - ((x) * (x)));
    const BHA_REAL tmp3 = tmp2 + 1;
    const BHA_REAL tmp4 = (1.0 / (tmp2));
    const BHA_REAL tmp5 = ((tmp3) * (tmp3) * (tmp3));
    const BHA_REAL x_np1 = x - (-C_r + E * tmp3 / M_PI) /
                               (-E * tmp4 * x / M_PI - 1.0 / 2.0 * tmp5 * (E - K) *
                                                           (-2 * tmp4 * ((x) * (x) * (x)) / tmp5 - 2 * x / ((tmp3) * (tmp3))) / (M_PI * ((x) * (x))));

    if (x_np1 > 1.0) {
      // Adjust the spin estimate to remain within valid bounds.
      spin = 2.0 - x_np1;
    } else if (x_np1 < 0) {
      // Terminate iteration if the spin magnitude estimate becomes negative.
      it = max_its;
      break;
    } else {
      // Calculate the relative difference and update the spin estimate.
      rel_diff = fabs(x_np1 - x) / x;
      spin = x_np1;
    } // END ELSE: spin adjustment to go back in-bounds
    it++;
  } // END WHILE: Refining spin estimate until convergence or maximum iterations

  // Assign spin=-10 if the Newton-Raphson did not converge within the allowed iterations.
  if (it >= max_its)
    spin = -10.0;

  return spin;
} // END FUNCTION: compute_spin

/**
 * Computes proper circumferences along the equator and polar directions for apparent horizon diagnostics.
 *
 * @param[in,out] commondata Pointer to common data structure containing shared parameters and settings.
 * @param[in,out] griddata Pointer to grid data structures for each grid, containing parameters and gridfunctions.
 * @return Status code indicating success or type of error (e.g., BHAHAHA_SUCCESS or INITIAL_DATA_MALLOC_ERROR).
 * @note This function uses OpenMP for parallel loops and performs interpolation and integration over grid data.
 */
int bah_diagnostics_proper_circumferences(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  const int NUM_DIAG_GFS = 2;
  const int grid = 0;
  // Extract grid dimensions, including ghost zones, for each coordinate direction. Needed for IDX4() macro.
  const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
  const int NUM_THETA = Nxx_plus_2NGHOSTS1; // Needed for IDX2() macro.

  BHA_REAL *restrict metric_data_gfs;
  BHAH_MALLOC(metric_data_gfs, Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_DIAG_GFS * sizeof(BHA_REAL));

  // Compute the line element in the phi direction (sqrt(q_{phi phi})) across the entire grid.
  {
    // Extract pointers to auxiliary and evolved gridfunctions, and coordinate arrays.
    const BHA_REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    const BHA_REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
    const BHA_REAL *restrict xx[3] = {griddata[grid].xx[0], griddata[grid].xx[1], griddata[grid].xx[2]};
    // Inverse grid spacings for theta and phi directions.
    const BHA_REAL invdxx1 = griddata[grid].params.invdxx1;
    const BHA_REAL invdxx2 = griddata[grid].params.invdxx2;
    const int i0 = NGHOSTS; // Fixed index for radial coordinate (r).

    // Loop over angular grid points (theta and phi) to compute:
    // 1. sqrt(q_{theta theta}), stored to metric_data_gfs[IDX4(0,...)], and
    // 2. sqrt(q_{phi phi}), stored to metric_data_gfs[IDX4(1,...)] at each point (theta, phi).
    // Notice we do clever indexing to ensure this 2D computation stays within the memory bounds of (3D) metric_data_gfs.
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      const MAYBE_UNUSED BHA_REAL xx2 = xx[2][i2]; // Phi coordinate at index i2.
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        const MAYBE_UNUSED BHA_REAL xx1 = xx[1][i1]; // Theta coordinate at index i1.
        static const BHA_REAL FDPart1_Rational_3_4 = 3.0 / 4.0;
        static const BHA_REAL FDPart1_Rational_3_20 = 3.0 / 20.0;
        static const BHA_REAL FDPart1_Rational_1_60 = 1.0 / 60.0;
        const BHA_REAL FDPart3tmp4 = sin(xx1);
        const BHA_REAL WW = auxevol_gfs[IDX4(WWGF, i0, i1, i2)];
        const BHA_REAL FDPart3tmp0 = (1.0 / ((WW) * (WW)));
        const BHA_REAL hDD00 = auxevol_gfs[IDX4(HDD00GF, i0, i1, i2)];
        const BHA_REAL FDPart3tmp2 = FDPart3tmp0 * (hDD00 + 1);
        const BHA_REAL hDD01 = auxevol_gfs[IDX4(HDD01GF, i0, i1, i2)];
        const BHA_REAL hDD02 = auxevol_gfs[IDX4(HDD02GF, i0, i1, i2)];
        const BHA_REAL hDD11 = auxevol_gfs[IDX4(HDD11GF, i0, i1, i2)];
        const BHA_REAL hDD22 = auxevol_gfs[IDX4(HDD22GF, i0, i1, i2)];
        const BHA_REAL hh_i2m3 = in_gfs[IDX4(HHGF, i0, i1, i2 - 3)];
        const BHA_REAL hh_i2m2 = in_gfs[IDX4(HHGF, i0, i1, i2 - 2)];
        const BHA_REAL hh_i2m1 = in_gfs[IDX4(HHGF, i0, i1, i2 - 1)];
        const BHA_REAL hh_i1m3 = in_gfs[IDX4(HHGF, i0, i1 - 3, i2)];
        const BHA_REAL hh_i1m2 = in_gfs[IDX4(HHGF, i0, i1 - 2, i2)];
        const BHA_REAL hh_i1m1 = in_gfs[IDX4(HHGF, i0, i1 - 1, i2)];
        const BHA_REAL hh = in_gfs[IDX4(HHGF, i0, i1, i2)];
        const BHA_REAL hh_i1p1 = in_gfs[IDX4(HHGF, i0, i1 + 1, i2)];
        const BHA_REAL hh_i1p2 = in_gfs[IDX4(HHGF, i0, i1 + 2, i2)];
        const BHA_REAL hh_i1p3 = in_gfs[IDX4(HHGF, i0, i1 + 3, i2)];
        const BHA_REAL hh_dD1 = invdxx1 * (FDPart1_Rational_1_60 * (-hh_i1m3 + hh_i1p3) + FDPart1_Rational_3_20 * (hh_i1m2 - hh_i1p2) +
                                       FDPart1_Rational_3_4 * (-hh_i1m1 + hh_i1p1));
        const BHA_REAL hh_i2p1 = in_gfs[IDX4(HHGF, i0, i1, i2 + 1)];
        const BHA_REAL hh_i2p2 = in_gfs[IDX4(HHGF, i0, i1, i2 + 2)];
        const BHA_REAL hh_i2p3 = in_gfs[IDX4(HHGF, i0, i1, i2 + 3)];
        const BHA_REAL hh_dD2 = invdxx2 * (FDPart1_Rational_1_60 * (-hh_i2m3 + hh_i2p3) + FDPart1_Rational_3_20 * (hh_i2m2 - hh_i2p2) +
                                       FDPart1_Rational_3_4 * (-hh_i2m1 + hh_i2p1));
        const BHA_REAL FDPart3tmp3 = ((hh) * (hh));
        const BHA_REAL FDPart3tmp1 = 2 * FDPart3tmp0 * hh;
        const BHA_REAL FDPart3tmp5 = FDPart3tmp3 * ((FDPart3tmp4) * (FDPart3tmp4));
        metric_data_gfs[IDX4pt(0, 0) + IDX2(i1, i2)] =
            sqrt(FDPart3tmp0 * (FDPart3tmp3 * hDD11 + FDPart3tmp3) + FDPart3tmp1 * hDD01 * hh_dD1 + FDPart3tmp2 * ((hh_dD1) * (hh_dD1)));
        metric_data_gfs[IDX4pt(1, 0) + IDX2(i1, i2)] = sqrt(FDPart3tmp0 * (FDPart3tmp5 * hDD22 + FDPart3tmp5) +
                                                            FDPart3tmp1 * FDPart3tmp4 * hDD02 * hh_dD2 + FDPart3tmp2 * ((hh_dD2) * (hh_dD2)));

      } // END LOOP: for i1 over theta points on the horizon surface
    } // END LOOP: for i2 over phi points on the horizon surface

    // Apply inner boundary conditions to the computed sqrt(q_{phi phi}) gridfunction.
    {
      bc_struct *restrict bcstruct = &griddata[grid].bcstruct; // Retrieve boundary condition structure for the grid.

      // Unpack bc_info from bcstruct.
      const bc_info_struct *bc_info = &bcstruct->bc_info;

      // Apply boundary conditions at inner boundary points for the selected gridfunctions.
#pragma omp parallel for collapse(2)
      for (int which_gf = 0; which_gf < 2; which_gf++) {
        for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
          const int dstpt = bcstruct->inner_bc_array[pt].dstpt; // Destination point index.
          const int srcpt = bcstruct->inner_bc_array[pt].srcpt; // Source point index for copying.

          // Extract the i0, i1, and i2 indices from dstpt and srcpt.
          //  -> idx3 = i + Nx0*(j + Nx1*k)
          //  -> i = mod(idx3, Nx0)
          //   tmp = (idx3-i)/Nx0
          //  -> j = mod(tmp, Nx1)
          //  -> k = (tmp-j)/Nx1

          const int dst_i0 = dstpt % Nxx_plus_2NGHOSTS0;
          const int dsttmp = (dstpt - dst_i0) / Nxx_plus_2NGHOSTS0;
          const int dst_i1 = dsttmp % Nxx_plus_2NGHOSTS1;
          const int dst_i2 = (dsttmp - dst_i1) / Nxx_plus_2NGHOSTS1;

          const int src_i0 = srcpt % Nxx_plus_2NGHOSTS0;
          const int srctmp = (srcpt - src_i0) / Nxx_plus_2NGHOSTS0;
          const int src_i1 = srctmp % Nxx_plus_2NGHOSTS1;
          const int src_i2 = (srctmp - src_i1) / Nxx_plus_2NGHOSTS1;

          // Apply boundary condition if at the radial interior point (i0 == NGHOSTS).
          if (dst_i0 == NGHOSTS) {
            metric_data_gfs[IDX4pt(which_gf, 0) + IDX2(dst_i1, dst_i2)] = metric_data_gfs[IDX4pt(which_gf, 0) + IDX2(src_i1, src_i2)];
          }
        } // END LOOP: for pt over inner boundary points
      } // END LOOP: for which_gf over gridfunctions
    } // END BLOCK: apply inner boundary conditions to circumference metric data
  } // END BLOCK: compute line-element gridfunctions on the horizon surface

  // Grid spacings in theta and phi directions.
  const BHA_REAL dxx1 = griddata[grid].params.dxx1;
  const BHA_REAL dxx2 = griddata[grid].params.dxx2;
  // Number of angular points to sample over 2 pi radians.
  const int N_angle = griddata[grid].params.Nxx2;
  // Allocate arrays for destination points (theta, phi) and circumference values.
  BHA_REAL(*dst_pts)[2] = malloc(N_angle * sizeof(*dst_pts));
  BHA_REAL *restrict circumference = malloc(N_angle * sizeof(BHA_REAL));
  // Check for successful memory allocation.
  if (dst_pts == NULL || circumference == NULL) {
    if (dst_pts != NULL)
      free(dst_pts);
    if (circumference != NULL)
      free(circumference);
    commondata->error_flag = DIAG_PROPER_CIRCUM_MALLOC_ERROR; // Return error code if allocation fails.
    return commondata->error_flag;
  }

  // Compute the angular increment for sampling over 2 pi radians, whether it be in theta (xz & yz planes) or phi (xy-plane).
  const BHA_REAL d_angle = (M_PI - (-M_PI)) / ((BHA_REAL)N_angle);

  // Equatorial (xy-plane) circumference first
  {
    // Initialize destination points along the equator (theta = pi/2) for interpolation.
#pragma omp parallel for
    for (int i2 = 0; i2 < N_angle; i2++) {
      dst_pts[i2][0] = M_PI / 2;                                   // Equator: theta = pi/2.
      dst_pts[i2][1] = -M_PI + ((BHA_REAL)i2 + (1.0 / 2.0)) * d_angle; // Equator: phi = [-pi, pi].
    } // END LOOP: for i2 over phi angles

    // Interpolate sqrt(q_{phi phi}) values onto the equator points to compute the circumference;
    //   note that sqrt(q_{phi phi}) is stored in metric_data_gfs[IDX4(1,...)]
    const int error =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, dxx1, dxx2, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, griddata[grid].xx,
                                                       &metric_data_gfs[IDX4pt(1, 0)], N_angle, dst_pts, circumference);
    if (error != BHAHAHA_SUCCESS)
      return error;

    // Retrieve integration weights for numerical integration over the sampled points.
    const BHA_REAL *restrict weights;
    int weight_stencil_size;
    bah_diagnostics_integration_weights(N_angle, N_angle, &weights, &weight_stencil_size);

    // Compute the total xy-plane circumference by integrating over the sampled points.
    BHA_REAL sum_circumference = 0.0;
#pragma omp parallel for reduction(+ : sum_circumference)
    for (int ic = 0; ic < N_angle; ic++) {
      const BHA_REAL weight = weights[ic % weight_stencil_size]; // Integration weight for this point.
      sum_circumference += circumference[ic] * weight;
    } // END LOOP: for ic over circumference samples
    // Multiply the sum by d[angle]
    commondata->bhahaha_diagnostics->xy_plane_circumference = sum_circumference * d_angle;
  } // END BLOCK: compute xy-plane proper circumference

  // Polar (xz-plane) circumference next
  {
    // Initialize destination points along the xz-plane for interpolation.
#pragma omp parallel for
    for (int i2 = 0; i2 < N_angle; i2++) {
      if (i2 < N_angle / 2) {
        // First half: Theta from 0 to pi, phi = 0
        dst_pts[i2][0] = ((BHA_REAL)i2 + 0.5) * (M_PI / ((BHA_REAL)(N_angle) / 2.0));
        dst_pts[i2][1] = 0.0;
      } else {
        // Second half: Theta from pi back to 0, phi = -pi
        dst_pts[i2][0] = ((BHA_REAL)(N_angle - i2) - 0.5) * (M_PI / ((BHA_REAL)(N_angle) / 2.0));
        dst_pts[i2][1] = -M_PI; // phi spans from [-pi, pi), so instead of interpolating at phi=pi, must interpolate at phi=-pi.
      } // END IF: theta is going from 0 to pi or vice-versa
    } // END LOOP: for i2 over angle samples

    // Interpolate sqrt(q_{theta theta}) values onto the polar (xz-plane) points to compute the circumference;
    //   note that sqrt(q_{theta theta}) is stored in metric_data_gfs[IDX4(0,...)]
    const int error =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, dxx1, dxx2, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, griddata[grid].xx,
                                                       &metric_data_gfs[IDX4pt(0, 0)], N_angle, dst_pts, circumference);
    if (error != BHAHAHA_SUCCESS)
      return error;

    // Retrieve integration weights for numerical integration over the sampled points.
    const BHA_REAL *restrict weights;
    int weight_stencil_size;
    bah_diagnostics_integration_weights(N_angle, N_angle, &weights, &weight_stencil_size);

    // Compute the total xz-plane circumference by integrating over the sampled points.
    BHA_REAL sum_circumference = 0.0;
#pragma omp parallel for reduction(+ : sum_circumference)
    for (int ic = 0; ic < N_angle; ic++) {
      const BHA_REAL weight = weights[ic % weight_stencil_size]; // Integration weight for this point.
      sum_circumference += circumference[ic] * weight;
    } // END LOOP: for ic over circumference samples
    // Multiply the sum by d[angle]
    commondata->bhahaha_diagnostics->xz_plane_circumference = sum_circumference * d_angle;
  } // END BLOCK: compute xz-plane proper circumference

  // Polar (yz-plane) circumference next
  {
    // Initialize destination points along the yz-plane for interpolation.
#pragma omp parallel for
    for (int i2 = 0; i2 < N_angle; i2++) {
      if (i2 < N_angle / 2) {
        // First half: Theta from 0 to pi, phi = pi/2
        dst_pts[i2][0] = ((BHA_REAL)i2 + 0.5) * (M_PI / ((BHA_REAL)(N_angle) / 2.0));
        dst_pts[i2][1] = M_PI / 2.0;
      } else {
        // Second half: Theta from pi back to 0, phi = -pi/2
        dst_pts[i2][0] = ((BHA_REAL)(N_angle - i2) - 0.5) * (M_PI / ((BHA_REAL)(N_angle) / 2.0));
        dst_pts[i2][1] = -M_PI / 2.0;
      } // END IF: theta is going from 0 to pi or vice-versa
    } // END LOOP: for i2 over angle samples

    // Interpolate sqrt(q_{theta theta}) values onto the polar (yz-plane) points to compute the circumference;
    //   note that sqrt(q_{theta theta}) is stored in metric_data_gfs[IDX4(0,...)]
    const int error =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, dxx1, dxx2, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, griddata[grid].xx,
                                                       &metric_data_gfs[IDX4pt(0, 0)], N_angle, dst_pts, circumference);
    if (error != BHAHAHA_SUCCESS)
      return error;

    // Retrieve integration weights for numerical integration over the sampled points.
    const BHA_REAL *restrict weights;
    int weight_stencil_size;
    bah_diagnostics_integration_weights(N_angle, N_angle, &weights, &weight_stencil_size);

    // Compute the total yz-plane circumference by integrating over the sampled points.
    BHA_REAL sum_circumference = 0.0;
#pragma omp parallel for reduction(+ : sum_circumference)
    for (int ic = 0; ic < N_angle; ic++) {
      const BHA_REAL weight = weights[ic % weight_stencil_size]; // Integration weight for this point.
      sum_circumference += circumference[ic] * weight;
    } // END LOOP: for ic over circumference samples
    // Multiply the sum by d[angle]
    commondata->bhahaha_diagnostics->yz_plane_circumference = sum_circumference * d_angle;
  } // END BLOCK: compute yz-plane proper circumference

  // Next estimate spin parameter magnitudes, valid for equilibrium BHs only.
  //   Based on Eq 5.2 of Alcubierre et al arXiv:gr-qc/0411149.
  {
    BHA_REAL C_xz_yz = commondata->bhahaha_diagnostics->xz_plane_circumference / commondata->bhahaha_diagnostics->yz_plane_circumference;
    BHA_REAL C_xy_yz = commondata->bhahaha_diagnostics->xy_plane_circumference / commondata->bhahaha_diagnostics->yz_plane_circumference;
    BHA_REAL C_yz_xz = commondata->bhahaha_diagnostics->yz_plane_circumference / commondata->bhahaha_diagnostics->xz_plane_circumference;
    BHA_REAL C_xy_xz = commondata->bhahaha_diagnostics->xy_plane_circumference / commondata->bhahaha_diagnostics->xz_plane_circumference;
    BHA_REAL C_xz_xy = commondata->bhahaha_diagnostics->xz_plane_circumference / commondata->bhahaha_diagnostics->xy_plane_circumference;
    BHA_REAL C_yz_xy = commondata->bhahaha_diagnostics->yz_plane_circumference / commondata->bhahaha_diagnostics->xy_plane_circumference;

    // Compute spin estimates for each direction based on different circumference ratios
    commondata->bhahaha_diagnostics->spin_a_x_from_xz_over_yz_prop_circumfs = compute_spin(C_xz_yz);
    commondata->bhahaha_diagnostics->spin_a_x_from_xy_over_yz_prop_circumfs = compute_spin(C_xy_yz);
    commondata->bhahaha_diagnostics->spin_a_y_from_yz_over_xz_prop_circumfs = compute_spin(C_yz_xz);
    commondata->bhahaha_diagnostics->spin_a_y_from_xy_over_xz_prop_circumfs = compute_spin(C_xy_xz);
    commondata->bhahaha_diagnostics->spin_a_z_from_xz_over_xy_prop_circumfs = compute_spin(C_xz_xy);
    commondata->bhahaha_diagnostics->spin_a_z_from_yz_over_xy_prop_circumfs = compute_spin(C_yz_xy);
  }
  // Free allocated memory for destination points, circumference values, and metric_data_gfs.
  free(dst_pts);
  free(circumference);
  free(metric_data_gfs);

  return BHAHAHA_SUCCESS; // Return success status code.
} // END FUNCTION: bah_diagnostics_proper_circumferences
