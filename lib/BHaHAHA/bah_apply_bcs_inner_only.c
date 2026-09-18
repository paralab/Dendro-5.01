#include "BHaH_defines.h"

static inline BHA_REAL apply_parity_branchless(const BHA_REAL v, const int8_t p) {
#ifdef __cplusplus
  static_assert(sizeof(BHA_REAL) == sizeof(uint32_t) || sizeof(BHA_REAL) == sizeof(uint64_t), "BHA_REAL must be float or double");
#else
  _Static_assert(sizeof(BHA_REAL) == sizeof(uint32_t) || sizeof(BHA_REAL) == sizeof(uint64_t), "BHA_REAL must be float or double");
#endif

  if (sizeof(BHA_REAL) == sizeof(uint64_t)) {
    uint64_t bits;
    BHA_REAL out;
    memcpy(&bits, &v, sizeof(bits));
    bits ^= (uint64_t)(p < 0) << 63;
    memcpy(&out, &bits, sizeof(out));
    return out;
  } else {
    uint32_t bits;
    BHA_REAL out;
    memcpy(&bits, &v, sizeof(bits));
    bits ^= (uint32_t)(p < 0) << 31;
    memcpy(&out, &bits, sizeof(out));
    return out;
  } // END IF 64 bits vs 32 bits
} // END FUNCTION: apply_parity_branchless

/**
 * Apply BCs to inner boundary points only,
 * using data stored in bcstruct->inner_bc_array.
 * These structs are set in bcstruct_set_up().
 * Inner boundary points map to either the grid
 * interior ("pure inner") or to pure outer
 * boundary points ("inner maps to outer").
 */
void bah_apply_bcs_inner_only(const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct,
                              BHA_REAL *restrict gfs) {
#include "set_CodeParameters.h"
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

#pragma omp parallel for schedule(static)
  for (int which_gf = 0; which_gf < NUM_EVOL_GFS; ++which_gf) {
    const int parity_idx = evol_gf_parity[which_gf];
    BHA_REAL *restrict gf = &gfs[IDX4pt(which_gf, 0)];

    for (int pt = 0; pt < bc_info->num_inner_boundary_points; ++pt) {
      const innerpt_bc_struct *restrict bc = &bcstruct->inner_bc_array[pt];
      const int dstpt = bc->dstpt;
      //  -> idx3 = i0 + Nx0*(i1 + Nx1*i2)
      //  -> i0 = mod(idx3, Nx0)
      // Only apply boundary condition if at the radial interior point (i0 == NGHOSTS).
      if (dstpt % Nxx_plus_2NGHOSTS0 != NGHOSTS)
        continue;

      const BHA_REAL v = gf[bc->srcpt];
      const int8_t p = bc->parity[parity_idx];
      gf[dstpt] = apply_parity_branchless(v, p);
    } // END LOOP: for pt over inner boundary points
  } // END LOOP: for which_gf over evolution gridfunctions
} // END FUNCTION: bah_apply_bcs_inner_only
