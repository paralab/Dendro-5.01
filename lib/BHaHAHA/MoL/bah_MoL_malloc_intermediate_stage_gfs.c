#include "BHaH_defines.h"

/**
 * Allocate intermediate-level (k_i) storage for the "SSPRK33" Method of Lines (MoL) scheme.
 *
 * This routine is registered as "MoL_malloc_intermediate_stage_gfs". At runtime it computes
 * the total number of grid points including ghost zones as:
 *     Nxx_plus_2NGHOSTS_tot = params->Nxx_plus_2NGHOSTS0
 *                            * params->Nxx_plus_2NGHOSTS1
 *                            * params->Nxx_plus_2NGHOSTS2;
 * For each required intermediate-level gridfunction array (with "auxevol_gfs" excluded),
 * the routine allocates:
 *     sizeof(BHA_REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot
 * The allocation macro is chosen by the build/runtime configuration:
 *     BHAH_MALLOC_DEVICE for CUDA builds, otherwise BHAH_MALLOC.
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
void bah_MoL_malloc_intermediate_stage_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                           MoL_gridfunctions_struct *restrict gridfuncs) {
  const int Nxx_plus_2NGHOSTS_tot = params->Nxx_plus_2NGHOSTS0 * params->Nxx_plus_2NGHOSTS1 * params->Nxx_plus_2NGHOSTS2;
  BHAH_MALLOC(gridfuncs->next_y_input_gfs, sizeof(BHA_REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  BHAH_MALLOC(gridfuncs->k1_gfs, sizeof(BHA_REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  BHAH_MALLOC(gridfuncs->k2_gfs, sizeof(BHA_REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  BHAH_MALLOC(gridfuncs->k3_gfs, sizeof(BHA_REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
} // END FUNCTION: bah_MoL_malloc_intermediate_stage_gfs
