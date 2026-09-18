#include "BHaH_defines.h"

/**
 * Free intermediate-level (k_i) storage for the "SSPRK33" Method of Lines (MoL) scheme.
 *
 * This routine is registered as "MoL_free_intermediate_stage_gfs".
 * It frees intermediate-level gridfunction storage.
 *
 * @param commondata Pointer to common, read-only runtime data. Included for a uniform
 *                   function signature across BHaH routines; it is not modified here.
 * @param params     Pointer to grid parameter struct providing Nxx_plus_2NGHOSTS0/1/2 and
 *                   related metadata needed to compute allocation sizes.
 * @param gridfuncs  Pointer to the MoL gridfunctions struct whose intermediate-level
 *                   arrays will be allocated by this routine.
 *
 * @return void
 */
void bah_MoL_free_intermediate_stage_gfs(MoL_gridfunctions_struct *restrict gridfuncs) {
  BHAH_FREE(gridfuncs->next_y_input_gfs);
  BHAH_FREE(gridfuncs->k1_gfs);
  BHAH_FREE(gridfuncs->k2_gfs);
  BHAH_FREE(gridfuncs->k3_gfs);
} // END FUNCTION: bah_MoL_free_intermediate_stage_gfs
