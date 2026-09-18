#include "BHaH_defines.h"

/**
 * Finite difference function for operator dKOD1, with FD accuracy order 4.
 */
static BHA_REAL fd_function_dKOD1_fdorder4(const BHA_REAL FDPROTO, const BHA_REAL FDPROTO_i1m1, const BHA_REAL FDPROTO_i1m2, const BHA_REAL FDPROTO_i1m3,
                                       const BHA_REAL FDPROTO_i1p1, const BHA_REAL FDPROTO_i1p2, const BHA_REAL FDPROTO_i1p3, const BHA_REAL invdxx1) {
  static const BHA_REAL FDPart1_Rational_5_16 = 5.0 / 16.0;
  static const BHA_REAL FDPart1_Rational_3_32 = 3.0 / 32.0;
  static const BHA_REAL FDPart1_Rational_1_64 = 1.0 / 64.0;
  static const BHA_REAL FDPart1_Rational_15_64 = 15.0 / 64.0;
  const BHA_REAL FD_result = invdxx1 * (-FDPROTO * FDPart1_Rational_5_16 + FDPart1_Rational_15_64 * (FDPROTO_i1m1 + FDPROTO_i1p1) +
                                    FDPart1_Rational_1_64 * (FDPROTO_i1m3 + FDPROTO_i1p3) + FDPart1_Rational_3_32 * (-FDPROTO_i1m2 - FDPROTO_i1p2));

  return FD_result;
} // END FUNCTION: fd_function_dKOD1_fdorder4
/**
 * Finite difference function for operator dKOD2, with FD accuracy order 4.
 */
static BHA_REAL fd_function_dKOD2_fdorder4(const BHA_REAL FDPROTO, const BHA_REAL FDPROTO_i2m1, const BHA_REAL FDPROTO_i2m2, const BHA_REAL FDPROTO_i2m3,
                                       const BHA_REAL FDPROTO_i2p1, const BHA_REAL FDPROTO_i2p2, const BHA_REAL FDPROTO_i2p3, const BHA_REAL invdxx2) {
  static const BHA_REAL FDPart1_Rational_5_16 = 5.0 / 16.0;
  static const BHA_REAL FDPart1_Rational_3_32 = 3.0 / 32.0;
  static const BHA_REAL FDPart1_Rational_1_64 = 1.0 / 64.0;
  static const BHA_REAL FDPart1_Rational_15_64 = 15.0 / 64.0;
  const BHA_REAL FD_result = invdxx2 * (-FDPROTO * FDPart1_Rational_5_16 + FDPart1_Rational_15_64 * (FDPROTO_i2m1 + FDPROTO_i2p1) +
                                    FDPart1_Rational_1_64 * (FDPROTO_i2m3 + FDPROTO_i2p3) + FDPart1_Rational_3_32 * (-FDPROTO_i2m2 - FDPROTO_i2p2));

  return FD_result;
} // END FUNCTION: fd_function_dKOD2_fdorder4

/**
 * Apply KO
 */
void bah_KO_apply(const commondata_struct *restrict commondata, const params_struct *restrict params, const rfm_struct *restrict rfmstruct,
                  const BHA_REAL *restrict auxevol_gfs, const BHA_REAL *restrict in_gfs, BHA_REAL *restrict rhs_gfs) {
#include "set_CodeParameters.h"
  if (commondata->KO_diss_strength == 0.0)
    return;
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      MAYBE_UNUSED const BHA_REAL f1_of_xx1 = rfmstruct->f1_of_xx1[i1];
      MAYBE_UNUSED const BHA_REAL f1_of_xx1__D1 = rfmstruct->f1_of_xx1__D1[i1];
      MAYBE_UNUSED const BHA_REAL f1_of_xx1__DD11 = rfmstruct->f1_of_xx1__DD11[i1];

      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        MAYBE_UNUSED const BHA_REAL f0_of_xx0 = rfmstruct->f0_of_xx0[i0];

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
        const BHA_REAL hh_dKOD1 = fd_function_dKOD1_fdorder4(hh, hh_i1m1, hh_i1m2, hh_i1m3, hh_i1p1, hh_i1p2, hh_i1p3, invdxx1);
        const BHA_REAL hh_i2p1 = in_gfs[IDX4(HHGF, i0, i1, i2 + 1)];
        const BHA_REAL hh_i2p2 = in_gfs[IDX4(HHGF, i0, i1, i2 + 2)];
        const BHA_REAL hh_i2p3 = in_gfs[IDX4(HHGF, i0, i1, i2 + 3)];
        const BHA_REAL hh_dKOD2 = fd_function_dKOD2_fdorder4(hh, hh_i2m1, hh_i2m2, hh_i2m3, hh_i2p1, hh_i2p2, hh_i2p3, invdxx2);
        const BHA_REAL FDPart3tmp0 = KO_diss_strength / hh;
        const BHA_REAL vv_i2m3 = in_gfs[IDX4(VVGF, i0, i1, i2 - 3)];
        const BHA_REAL vv_i2m2 = in_gfs[IDX4(VVGF, i0, i1, i2 - 2)];
        const BHA_REAL vv_i2m1 = in_gfs[IDX4(VVGF, i0, i1, i2 - 1)];
        const BHA_REAL vv_i1m3 = in_gfs[IDX4(VVGF, i0, i1 - 3, i2)];
        const BHA_REAL vv_i1m2 = in_gfs[IDX4(VVGF, i0, i1 - 2, i2)];
        const BHA_REAL vv_i1m1 = in_gfs[IDX4(VVGF, i0, i1 - 1, i2)];
        const BHA_REAL vv = in_gfs[IDX4(VVGF, i0, i1, i2)];
        const BHA_REAL vv_i1p1 = in_gfs[IDX4(VVGF, i0, i1 + 1, i2)];
        const BHA_REAL vv_i1p2 = in_gfs[IDX4(VVGF, i0, i1 + 2, i2)];
        const BHA_REAL vv_i1p3 = in_gfs[IDX4(VVGF, i0, i1 + 3, i2)];
        const BHA_REAL vv_dKOD1 = fd_function_dKOD1_fdorder4(vv, vv_i1m1, vv_i1m2, vv_i1m3, vv_i1p1, vv_i1p2, vv_i1p3, invdxx1);
        const BHA_REAL vv_i2p1 = in_gfs[IDX4(VVGF, i0, i1, i2 + 1)];
        const BHA_REAL vv_i2p2 = in_gfs[IDX4(VVGF, i0, i1, i2 + 2)];
        const BHA_REAL vv_i2p3 = in_gfs[IDX4(VVGF, i0, i1, i2 + 3)];
        const BHA_REAL vv_dKOD2 = fd_function_dKOD2_fdorder4(vv, vv_i2m1, vv_i2m2, vv_i2m3, vv_i2p1, vv_i2p2, vv_i2p3, invdxx2);
        const BHA_REAL FDPart3tmp1 = FDPart3tmp0 / f1_of_xx1;
        rhs_gfs[IDX4(HHGF, i0, i1, i2)] = FDPart3tmp0 * hh_dKOD1 + FDPart3tmp1 * hh_dKOD2 + rhs_gfs[IDX4(HHGF, i0, i1, i2)];
        rhs_gfs[IDX4(VVGF, i0, i1, i2)] = FDPart3tmp0 * vv_dKOD1 + FDPart3tmp1 * vv_dKOD2 + rhs_gfs[IDX4(VVGF, i0, i1, i2)];

      } // END LOOP: for i0 over [NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS)
    } // END LOOP: for i1 over [NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS)
  } // END LOOP: for i2 over [NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS)
} // END FUNCTION: bah_KO_apply
