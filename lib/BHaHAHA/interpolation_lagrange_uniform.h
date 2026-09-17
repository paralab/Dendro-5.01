/**
 * @file interpolation_lagrange_uniform.h
 * @brief Lagrange interpolation helpers on a uniform source grid.
 */

#ifndef INTERPOLATION_LAGRANGE_UNIFORM_H_
#define INTERPOLATION_LAGRANGE_UNIFORM_H_

#include "intrinsics/simd_intrinsics.h"

/**
 * @brief Define BHA_REAL data type. Defaults to double.
 */
#ifndef BHA_REAL
#define BHA_REAL double
#endif

#ifndef RESTRICT
#if defined(_MSC_VER)
/* MSVC */
#define RESTRICT __restrict
#elif defined(__cplusplus)
/* GNU/Clang C++ */
#define RESTRICT __restrict__
#else
/* C99 */
#define RESTRICT restrict
#endif
#endif

/**
 * @brief Precompute inverse denominators for Lagrange interpolation coefficients.
 *
 * Precomputes inverse denominators to optimize performance by avoiding repeated divisions during interpolation.  This significantly speeds up the
 * computation of Lagrange basis coefficients.
 *
 * @param INTERP_ORDER The order of the interpolation.
 * @param inv_denom Array to store the precomputed inverse denominators.  Must be of size INTERP_ORDER.
 */
static inline void compute_inv_denom(const int INTERP_ORDER, BHA_REAL *RESTRICT inv_denom) {
  for (int i = 0; i < INTERP_ORDER; i++) {
    BHA_REAL denom = 1.0;
    for (int j = 0; j < i; j++)
      denom *= (BHA_REAL)(i - j);
    for (int j = i + 1; j < INTERP_ORDER; j++)
      denom *= (BHA_REAL)(i - j);
    inv_denom[i] = 1.0 / denom;
  } // END LOOP: Precompute inverse denominators.
} // END FUNCTION: compute_inv_denom()

/**
 * @brief Compute differences between a destination point and source stencil points.
 *
 * Calculates the differences between a destination coordinate and each coordinate in the source stencil.  These differences are crucial for computing
 * Lagrange basis coefficients.
 *
 * @param INTERP_ORDER The order of the interpolation.
 * @param dst_xi The destination coordinate.
 * @param src_xi_stencil Pointer to the array of source stencil coordinates.
 * @param diffs Array to store the computed differences. Must be of size INTERP_ORDER.
 */
static inline void compute_diffs_xi(const int INTERP_ORDER, const BHA_REAL dst_xi, const BHA_REAL *RESTRICT src_xi_stencil, BHA_REAL *RESTRICT diffs) {
#pragma omp simd
  for (int j = 0; j < INTERP_ORDER; j++) {
    diffs[j] = dst_xi - src_xi_stencil[j];
  } // END LOOP over j.
} // END FUNCTION: compute_diffs_xi()

/**
 * @brief Compute Lagrange basis coefficients.
 *
 * Computes the Lagrange basis coefficients using precomputed inverse denominators and the differences between the destination and source stencil
 * coordinates.
 *
 * @param INTERP_ORDER The order of the interpolation.
 * @param inv_denom Array of precomputed inverse denominators.
 * @param diffs Array of differences between destination and source stencil coordinates.
 * @param lagrange_basis_coeffs_xi Array to store the computed Lagrange basis coefficients. Must be of size INTERP_ORDER.
 */
static inline void compute_lagrange_basis_coeffs_xi(const int INTERP_ORDER, const BHA_REAL *RESTRICT inv_denom, const BHA_REAL *RESTRICT diffs,
                                                    BHA_REAL *RESTRICT lagrange_basis_coeffs_xi) {
#pragma omp simd
  for (int i = 0; i < INTERP_ORDER; i++) {
    BHA_REAL numer_i = 1.0;
    // Compute product for j < i (scalar loop).
    for (int j = 0; j < i; j++) {
      numer_i *= diffs[j];
    } // END LOOP over j < i.
    // Compute product for j > i (scalar loop).
    for (int j = i + 1; j < INTERP_ORDER; j++) {
      numer_i *= diffs[j];
    } // END LOOP over j > i.
    lagrange_basis_coeffs_xi[i] = numer_i * inv_denom[i];
  } // END LOOP over i.
} // END FUNCTION: compute_lagrange_basis_coeffs_xi()

/**
 * @brief Compute the weighted sum of Lagrange basis coefficients and source data using SIMD instructions where possible.
 *
 * This function computes the weighted sum of the Lagrange basis coefficients and the corresponding source data values.  It leverages SIMD
 * instructions for improved performance when the interpolation order allows.  The code handles cases where the interpolation order is not a multiple
 * of the SIMD vector width.
 *
 * @param INTERP_ORDER The order of the interpolation.
 * @param src_gf_base_idx Pointer to the array of source data values.
 * @param lagrange_basis_coeffs_x0_base_idx Pointer to the array of Lagrange basis coefficients.
 * @return The weighted sum of the Lagrange interpolation.
 */
static inline BHA_REAL sum_lagrange_x0_simd(const int INTERP_ORDER, const BHA_REAL *RESTRICT src_gf_base_idx,
                                        const BHA_REAL *RESTRICT lagrange_basis_coeffs_x0_base_idx) {
  BHA_REAL sum = 0;

  // Vectorized loop over ix0 using SIMD with FMA, if available
  int ix0 = 0;
  REAL_SIMD_ARRAY vec_sum = SetZeroSIMD;
  for (; ix0 <= INTERP_ORDER - SIMD_WIDTH; ix0 += SIMD_WIDTH) {
    vec_sum = FusedMulAddSIMD(ReadSIMD(&src_gf_base_idx[ix0]), ReadSIMD(&lagrange_basis_coeffs_x0_base_idx[ix0]), vec_sum);
  } // END LOOP x0 direction up to integer number of SIMD widths
  // Accumulate SIMD result
  sum += HorizAddSIMD(vec_sum);

  // Handle remaining elements that don't fit into a full SIMD register
  for (; ix0 < INTERP_ORDER; ix0++) {
    sum += src_gf_base_idx[ix0] * lagrange_basis_coeffs_x0_base_idx[ix0];
  } // END LOOP remainder of x0 direction
  return sum;
} // END FUNCTION: lagrange_sum_x0()

#endif // INTERPOLATION_LAGRANGE_UNIFORM_H_
