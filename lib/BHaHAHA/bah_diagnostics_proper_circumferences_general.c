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

// Apply inner BCs for a selection of gridfunctions
// Note: Nxx_plus_2NGHOSTS2 is needed for IDX4pt()
static void apply_inner_bc_for_selected_gfs(bc_struct *restrict bcstruct, BHA_REAL *restrict metric_data_gfs, const int Nxx_plus_2NGHOSTS0,
                                            const int Nxx_plus_2NGHOSTS1, const int Nxx_plus_2NGHOSTS2, const int *which_gfs, const int num_gfs) {

  const bc_info_struct *bc_info = &bcstruct->bc_info;
  const int NUM_THETA = Nxx_plus_2NGHOSTS1; // Needed for IDX2

  // Apply boundary conditions at inner boundary points for the selected gridfunctions.
#pragma omp parallel for collapse(2)
  for (int gf_idx = 0; gf_idx < num_gfs; gf_idx++) {
    for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
      const int which_gf = which_gfs[gf_idx];
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

      if (dst_i0 == NGHOSTS) {
        metric_data_gfs[IDX4pt(which_gf, 0) + IDX2(dst_i1, dst_i2)] = metric_data_gfs[IDX4pt(which_gf, 0) + IDX2(src_i1, src_i2)];
      }
    } // END LOOP: for pt over inner boundary points
  } // END LOOP: for which_gf over gridfunctions
} // END FUNCTION: apply_inner_bc_for_selected_gfs

// ---------- small vector & math helpers (pure C) ----------
// These are intentionally tiny & inlined: they live in tight OpenMP loops,
// and we want predictable codegen and good autovectorization.
static inline BHA_REAL dot3(const BHA_REAL a[3], const BHA_REAL b[3]) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
static inline void cross3(const BHA_REAL a[3], const BHA_REAL b[3], BHA_REAL c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
} // END FUNCTION: cross3
static inline BHA_REAL norm3(const BHA_REAL a[3]) { return sqrt(dot3(a, a)); }
static inline void normalize3(BHA_REAL a[3]) {
  const BHA_REAL n = norm3(a);
  if (n > 0.0) {
    a[0] /= n;
    a[1] /= n;
    a[2] /= n;
  }
}
// Build an orthonormal basis {e1, e2} orthogonal to s; avoid near-collinearity for numerical stability.
static inline void build_basis_from_s(const BHA_REAL s[3], BHA_REAL e1[3], BHA_REAL e2[3]) {
  BHA_REAL a[3] = {1.0, 0.0, 0.0};
  if (fabs(dot3(a, s)) > 0.9) {
    a[0] = 0.0;
    a[1] = 1.0;
    a[2] = 0.0;
  }
  cross3(s, a, e1);
  normalize3(e1);
  cross3(s, e1, e2);
  normalize3(e2);
} // END FUNCTION: build_basis_from_s
// Cartesian unit vector -> spherical angles (theta, phi); inputs here are unit by construction.
static inline void cart_to_sph(const BHA_REAL r[3], BHA_REAL *theta, BHA_REAL *phi) {
  const BHA_REAL x = r[0], y = r[1], z = r[2];
  BHA_REAL zz = z;
  if (zz < -1.0)
    zz = -1.0;
  if (zz > 1.0)
    zz = 1.0;
  *theta = acos(zz);
  *phi = atan2(y, x);
} // END FUNCTION: cart_to_sph

// Midpoint integrate over alpha with precomputed 8th-order weights (N_angle == Nxx2 is typically divisible by 8 -> 8th order).
static inline BHA_REAL integrate_over_alpha(const BHA_REAL *vals, int N_angle, BHA_REAL d_alpha) {
  const BHA_REAL *restrict weights;
  int weight_stencil_size;
  bah_diagnostics_integration_weights(N_angle, N_angle, &weights, &weight_stencil_size);
  BHA_REAL sum = 0.0;
  // #pragma omp parallel for reduction(+ : sum) // <- N_angle ~64 => thread creation/destruction makes OMP SLOWER
  for (int i = 0; i < N_angle; i++)
    sum += vals[i] * weights[i % weight_stencil_size];
  return sum * d_alpha;
} // END FUNCTION: integrate_over_alpha
/* ---------- end helpers ---------- */

/**
 * Computes proper circumferences along the equator and polar directions with respect to the provided spin axis,
 * using the induced 2-metric q_{AB} on the surface r = h(theta,phi).
 *
 * Line element integrated (alpha parametrizes the great circle):
 *     ds = sqrt( q_tt * (dtheta/dalpha)^2 + 2*q_tp*(dtheta/dalpha)*(dphi/dalpha) + q_pp * (dphi/dalpha)^2 ) dalpha,
 * with q_tt = (sqrt_qtt)^2, q_pp = (sqrt_qpp)^2, and q_tp as-is. Storing sqrt(diagonals) improves interpolation stability
 * and preserves positivity; keeping q_tp raw preserves the sign required by the cross term.
 *
 * @param[in,out] commondata Pointer to common data structure containing shared parameters and settings.
 * @param[in,out] griddata Pointer to grid data structures for each grid, containing parameters and gridfunctions.
 * @return Status code indicating success or type of error (e.g., BHAHAHA_SUCCESS or INITIAL_DATA_MALLOC_ERROR).
 * @note Uses OpenMP for parallel loops; performs interpolation and midpoint integration over angular samples.
 *
 * Precomputation strategy on the (theta,phi) grid at fixed i0=NGHOSTS:
 *   sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi} are generated via NRPy's SymPy expressions
 *   (nrpy/equations/general_relativity/bhahaha/area.py::circumference_metric_roots). In addition, we precompute scalar integrands
 *   f_eq(theta,phi) and f_pol(theta,phi) using smooth great-circle tangent fields, so the integration phase
 *   only interpolates a single scalar per sample. This avoids finite differencing of angles entirely
 *   and aligns the derivative accuracy with the high-order midpoint quadrature.
 */
int bah_diagnostics_proper_circumferences_general(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  // which_gf in {0: sqrt(q_tt), 1: sqrt(q_pp), 2: q_tp, 3: f_eq(θ,φ), 4: f_pol(θ,φ)};
  // diagonals as sqrt for stability; f_eq/f_pol are scalar line-element integrands.
  const int NUM_DIAG_GFS = 5;
  const int grid = 0;
  // Extract grid dimensions, including ghost zones, for each coordinate direction. Needed for IDX4() macro.
  const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
  const int NUM_THETA = Nxx_plus_2NGHOSTS1; // Needed for IDX2() macro.

  BHA_REAL *restrict metric_data_gfs;
  // Single heap allocation sized for a 2D (theta,phi) slab at i0 = NGHOSTS.
  BHAH_MALLOC(metric_data_gfs, Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_DIAG_GFS * sizeof(BHA_REAL));

  // Compute sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi} across the entire (theta,phi) grid at i0=NGHOSTS.
  {
    // Extract pointers to auxiliary and evolved gridfunctions, and coordinate arrays.
    const BHA_REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    const BHA_REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
    const BHA_REAL *restrict xx[3] = {griddata[grid].xx[0], griddata[grid].xx[1], griddata[grid].xx[2]};
    // Inverse grid spacings for theta and phi directions.
    const BHA_REAL invdxx1 = griddata[grid].params.invdxx1;
    const BHA_REAL invdxx2 = griddata[grid].params.invdxx2;
    const int i0 = NGHOSTS; // Fixed index for radial coordinate (r) where r=h(theta,phi) lives.

    // Loop over angular grid points (theta and phi) to compute:
    // 1. sqrt(q_{theta theta}), stored to metric_data_gfs[IDX4(0,...)], and
    // 2. sqrt(q_{phi phi}),    stored to metric_data_gfs[IDX4(1,...)] at each point (theta, phi).
    // 3. q_{theta phi},       stored to metric_data_gfs[IDX4(2,...)]. Keep sign; no sqrt.
    // Clever IDX math ensures this 2D computation stays within the memory bounds of the 3D allocation.
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      const MAYBE_UNUSED BHA_REAL xx2 = xx[2][i2]; // Phi coordinate at index i2.
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        const MAYBE_UNUSED BHA_REAL xx1 = xx[1][i1]; // Theta coordinate at index i1.
        static const BHA_REAL FDPart1_Rational_3_4 = 3.0 / 4.0;
        static const BHA_REAL FDPart1_Rational_3_20 = 3.0 / 20.0;
        static const BHA_REAL FDPart1_Rational_1_60 = 1.0 / 60.0;
        const BHA_REAL FDPart3tmp5 = sin(xx1);
        const BHA_REAL WW = auxevol_gfs[IDX4(WWGF, i0, i1, i2)];
        const BHA_REAL FDPart3tmp0 = (1.0 / ((WW) * (WW)));
        const BHA_REAL hDD00 = auxevol_gfs[IDX4(HDD00GF, i0, i1, i2)];
        const BHA_REAL FDPart3tmp3 = FDPart3tmp0 * (hDD00 + 1);
        const BHA_REAL hDD01 = auxevol_gfs[IDX4(HDD01GF, i0, i1, i2)];
        const BHA_REAL hDD02 = auxevol_gfs[IDX4(HDD02GF, i0, i1, i2)];
        const BHA_REAL hDD11 = auxevol_gfs[IDX4(HDD11GF, i0, i1, i2)];
        const BHA_REAL hDD12 = auxevol_gfs[IDX4(HDD12GF, i0, i1, i2)];
        const BHA_REAL hDD22 = auxevol_gfs[IDX4(HDD22GF, i0, i1, i2)];
        const BHA_REAL hh_i2m3 = in_gfs[IDX4(HHGF, i0, i1, i2 - 3)];
        const BHA_REAL hh_i2m2 = in_gfs[IDX4(HHGF, i0, i1, i2 - 2)];
        const BHA_REAL hh_i2m1 = in_gfs[IDX4(HHGF, i0, i1, i2 - 1)];
        const BHA_REAL hh_i1m3 = in_gfs[IDX4(HHGF, i0, i1 - 3, i2)];
        const BHA_REAL hh_i1m2 = in_gfs[IDX4(HHGF, i0, i1 - 2, i2)];
        const BHA_REAL hh_i1m1 = in_gfs[IDX4(HHGF, i0, i1 - 1, i2)];
        const BHA_REAL hh = in_gfs[IDX4(HHGF, i0, i1, i2)];
        const BHA_REAL FDPart3tmp2 = FDPart3tmp0 * hDD01 * hh;
        const BHA_REAL FDPart3tmp6 = FDPart3tmp0 * FDPart3tmp5 * hDD02 * hh;
        const BHA_REAL FDPart3tmp4 = ((hh) * (hh));
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
        const BHA_REAL FDPart3tmp7 = FDPart3tmp4 * ((FDPart3tmp5) * (FDPart3tmp5));
        metric_data_gfs[IDX4pt(0, 0) + IDX2(i1, i2)] =
            sqrt(FDPart3tmp0 * (FDPart3tmp4 * hDD11 + FDPart3tmp4) + 2 * FDPart3tmp2 * hh_dD1 + FDPart3tmp3 * ((hh_dD1) * (hh_dD1)));
        metric_data_gfs[IDX4pt(1, 0) + IDX2(i1, i2)] =
            sqrt(FDPart3tmp0 * (FDPart3tmp7 * hDD22 + FDPart3tmp7) + FDPart3tmp3 * ((hh_dD2) * (hh_dD2)) + 2 * FDPart3tmp6 * hh_dD2);
        metric_data_gfs[IDX4pt(2, 0) + IDX2(i1, i2)] =
            FDPart3tmp0 * FDPart3tmp4 * FDPart3tmp5 * hDD12 + FDPart3tmp2 * hh_dD2 + FDPart3tmp3 * hh_dD1 * hh_dD2 + FDPart3tmp6 * hh_dD1;

      } // END LOOP: for i1 over theta points on the horizon surface
    } // END LOOP: for i2 over phi points on the horizon surface

    // Apply inner boundary conditions to q-metric gridfunctions sqrt(qtt), sqrt(qpp), and qtp:
    {
      const int which_gfs_1[3] = {0, 1, 2};
      apply_inner_bc_for_selected_gfs(&griddata[grid].bcstruct, metric_data_gfs, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
                                      which_gfs_1, 3);
    } // END BLOCK: apply inner boundary conditions to q-metric root gridfunctions
  } // END BLOCK: compute q-metric root gridfunctions on the horizon surface

  // Number of angular points to sample over 2 pi radians (controls resolution of great-circle integrals).
  const int N_angle = griddata[grid].params.Nxx2;
  // Uniform alpha step for midpoint samples in [-pi, pi).
  const BHA_REAL d_alpha = (M_PI - (-M_PI)) / ((BHA_REAL)N_angle);

  BHA_REAL dst_pts[N_angle][2];
  BHA_REAL theta[N_angle];
  BHA_REAL phi[N_angle];
  BHA_REAL integrand[N_angle];

  // Normalize spin axis; if zero-length is provided, fall back to z-axis for determinism.
  BHA_REAL s[3] = {
      commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_X, // unit spin direction, x-component
      commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_Y, // unit spin direction, y-component
      commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_Z  // unit spin direction, z-component
  };
  BHA_REAL s_norm = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
  if (s_norm > 0.0) {
    s[0] /= s_norm;
    s[1] /= s_norm;
    s[2] /= s_norm;
  } else {
    s[0] = 0.0;
    s[1] = 0.0;
    s[2] = 1.0;
  } // END ELSE: fallback to z-axis when spin axis norm vanishes
  commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_X = s[0];
  commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_Y = s[1];
  commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_Z = s[2];

  // Build an orthonormal basis {e1, e2} orthogonal to s, to define great circles.
  BHA_REAL e1[3], e2[3];
  build_basis_from_s(s, e1, e2);
  // Normal of the polar great-circle plane (spanned by s and e1):
  BHA_REAL nvec[3];
  cross3(s, e1, nvec);

  // ================================================================
  // Precompute integrand scalars f_eq(θ,φ) and f_pol(θ,φ) on the i0=NGHOSTS slab.
  // Each uses a smooth tangent field on the unit sphere:
  //   equator: t = normalize(s × r),   polar: t = normalize((s × e1) × r)
  // and the spherical basis (θ̂, φ̂) to build dθ/dalpha = t·θ̂, dφ/dalpha = (t·φ̂)/sinθ.
  // This avoids any finite differencing of angles and eliminates branch cuts.
  // ================================================================
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    const BHA_REAL phi_c = griddata[grid].xx[2][i2];
    const BHA_REAL sinph = sin(phi_c), cosph = cos(phi_c);
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      const BHA_REAL theta_c = griddata[grid].xx[1][i1];
      const BHA_REAL sinth = sin(theta_c), costh = cos(theta_c);
      // unit position and spherical basis
      const BHA_REAL rx = sinth * cosph, ry = sinth * sinph, rz = costh;
      const BHA_REAL thx = costh * cosph, thy = costh * sinph, thz = -sinth;
      const BHA_REAL phx = -sinph, phy = cosph, phz = 0.0;

      // load q-metric pieces (recall diagonals are stored as square-roots)
      const BHA_REAL sqtt = metric_data_gfs[IDX4pt(0, 0) + IDX2(i1, i2)];
      const BHA_REAL sqpp = metric_data_gfs[IDX4pt(1, 0) + IDX2(i1, i2)];
      const BHA_REAL qtp = metric_data_gfs[IDX4pt(2, 0) + IDX2(i1, i2)];
      const BHA_REAL qtt = sqtt * sqtt;
      const BHA_REAL qpp = sqpp * sqpp;

      // ---------- equator tangent: t = normalize(s × r) ----------
      BHA_REAL tex = s[1] * rz - s[2] * ry;
      BHA_REAL tey = s[2] * rx - s[0] * rz;
      BHA_REAL tez = s[0] * ry - s[1] * rx;
      BHA_REAL tnorm = sqrt(tex * tex + tey * tey + tez * tez);
      if (tnorm > 1e-14) {
        tex /= tnorm;
        tey /= tnorm;
        tez /= tnorm;
      } else {
        tex = tey = tez = 0.0;
      }
      const BHA_REAL dth_eq = tex * thx + tey * thy + tez * thz;
      const BHA_REAL dph_eq = (tex * phx + tey * phy + tez * phz) / fmax(1e-14, sinth);
      const BHA_REAL feq = sqrt(qtt * dth_eq * dth_eq + 2.0 * qtp * dth_eq * dph_eq + qpp * dph_eq * dph_eq);
      metric_data_gfs[IDX4pt(3, 0) + IDX2(i1, i2)] = feq;

      // ---------- polar tangent: plane normal n = s × e1; t = normalize(n × r) ----------
      BHA_REAL tpx = nvec[1] * rz - nvec[2] * ry;
      BHA_REAL tpy = nvec[2] * rx - nvec[0] * rz;
      BHA_REAL tpz = nvec[0] * ry - nvec[1] * rx;
      tnorm = sqrt(tpx * tpx + tpy * tpy + tpz * tpz);
      if (tnorm > 1e-14) {
        tpx /= tnorm;
        tpy /= tnorm;
        tpz /= tnorm;
      } else {
        tpx = tpy = tpz = 0.0;
      }
      const BHA_REAL dth_pol = tpx * thx + tpy * thy + tpz * thz;
      const BHA_REAL dph_pol = (tpx * phx + tpy * phy + tpz * phz) / fmax(1e-14, sinth);
      const BHA_REAL fpol = sqrt(qtt * dth_pol * dth_pol + 2.0 * qtp * dth_pol * dph_pol + qpp * dph_pol * dph_pol);
      metric_data_gfs[IDX4pt(4, 0) + IDX2(i1, i2)] = fpol;
    } // END LOOP: for i1 over theta
  } // END LOOP: for i2 over phi

  // Apply inner boundary conditions to the newly computed scalar integrands f_eq (which_gf=3)
  // and f_pol (which_gf=4), so their ghost zones are valid prior to interpolation.
  {
    const int which_gfs_2[2] = {3, 4};
    apply_inner_bc_for_selected_gfs(&griddata[grid].bcstruct, metric_data_gfs, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
                                    which_gfs_2, 2);
  } // END BLOCK: apply inner boundary conditions to equatorial and polar integrands

  // Convenience base pointers for scalar integrands:
  const BHA_REAL *restrict src_feq = &metric_data_gfs[IDX4pt(3, 0)];
  const BHA_REAL *restrict src_fpol = &metric_data_gfs[IDX4pt(4, 0)];

  // Some compilers warn on constness; build an explicit coords pointer matching prototype.
  BHA_REAL *restrict(*coords)[3] = (BHA_REAL *restrict(*)[3]) & griddata[grid].xx;

  // ================================================================
  // Equatorial great circle: orthogonal to s.
  //   Parameterization: r(alpha) = cos(alpha) e1 + sin(alpha) e2, alpha in [-pi, pi).
  //   Convert to (theta, phi), interpolate precomputed scalar integrand, and integrate.
  // ================================================================
  // #pragma omp parallel for // <- N_angle ~64 => thread creation/destruction makes OMP SLOWER
  for (int i = 0; i < N_angle; i++) {
    const BHA_REAL alpha = -M_PI + ((BHA_REAL)i + 0.5) * d_alpha;
    BHA_REAL rvec[3] = {cos(alpha) * e1[0] + sin(alpha) * e2[0], cos(alpha) * e1[1] + sin(alpha) * e2[1], cos(alpha) * e1[2] + sin(alpha) * e2[2]};
    cart_to_sph(rvec, &theta[i], &phi[i]);
    dst_pts[i][0] = theta[i];
    dst_pts[i][1] = phi[i];
  } // END LOOP: for i over alpha: cell-centered sampling 2 pi across N_angle points

  {
    // Interpolate precomputed equator integrand f_eq onto the alpha-midpoints.
    int err =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, griddata[grid].params.dxx1, griddata[grid].params.dxx2, Nxx_plus_2NGHOSTS1,
                                                       Nxx_plus_2NGHOSTS2, (BHA_REAL *restrict *)(*coords), src_feq, N_angle, dst_pts, integrand);
    if (err != BHAHAHA_SUCCESS) {
      free(metric_data_gfs);
      return err;
    } // END IF: equatorial-integrand interpolation failed
  } // END BLOCK: interpolate equatorial integrand onto great-circle samples

  const BHA_REAL C_equator = integrate_over_alpha(integrand, N_angle, d_alpha);

  // ================================================================
  // Polar great circle: contains s.
  //   Parameterization: r(alpha) = cos(alpha) s + sin(alpha) e1, alpha in [-pi, pi).
  //   Proceed as above.
  // ================================================================
  // #pragma omp parallel for reduction(+ : sum) // <- N_angle ~64 => thread creation/destruction makes OMP SLOWER
  for (int i = 0; i < N_angle; i++) {
    const BHA_REAL alpha = -M_PI + ((BHA_REAL)i + 0.5) * d_alpha;
    BHA_REAL rvec[3] = {cos(alpha) * s[0] + sin(alpha) * e1[0], cos(alpha) * s[1] + sin(alpha) * e1[1], cos(alpha) * s[2] + sin(alpha) * e1[2]};
    cart_to_sph(rvec, &theta[i], &phi[i]);
    dst_pts[i][0] = theta[i];
    dst_pts[i][1] = phi[i];
  } // END LOOP: for i over alpha: cell-centered sampling 2 pi across N_angle points

  {
    // Interpolate precomputed polar integrand f_pol onto the alpha-midpoints.
    int err =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, griddata[grid].params.dxx1, griddata[grid].params.dxx2, Nxx_plus_2NGHOSTS1,
                                                       Nxx_plus_2NGHOSTS2, (BHA_REAL *restrict *)(*coords), src_fpol, N_angle, dst_pts, integrand);
    if (err != BHAHAHA_SUCCESS) {
      free(metric_data_gfs);
      return err;
    } // END IF: polar-integrand interpolation failed
  } // END BLOCK: interpolate polar integrand onto great-circle samples

  const BHA_REAL C_polar = integrate_over_alpha(integrand, N_angle, d_alpha);

  // Compute ratio and spin estimate; write results to macro-backed fields.
  const BHA_REAL ratio = C_polar / C_equator; // C_polar / C_equator wrt spin axis
  const BHA_REAL a_est = compute_spin(ratio); // compute_spin(circumf_ratio_polar_over_equator)

  commondata->bhahaha_diagnostics->BHAHAHA_CIRC_GENERAL_POLAR = C_polar;     // polar_circ
  commondata->bhahaha_diagnostics->BHAHAHA_CIRC_GENERAL_EQUATOR = C_equator; // equator_circ
  commondata->bhahaha_diagnostics->BHAHAHA_CIRC_GENERAL_SPIN = a_est;        // spin_a_along_spin_axis_from_circumf

  // Free the only heap allocation from this scope.
  free(metric_data_gfs);

  return BHAHAHA_SUCCESS; // Return success status code.
} // END FUNCTION: bah_diagnostics_proper_circumferences_general
