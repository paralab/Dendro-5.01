#include "BHaH_defines.h"

/**
 * @file rfm_precompute_free.c
 * Frees memory for precomputed reference metric quantities.
 *
 * - params->is_host==true: free HOST memory.
 * - params->is_host==false: free DEVICE memory.
 *
 * @param[in]  commondata  Global simulation metadata (unused).
 * @param[in]  params      Grid parameters.
 * @param[out] rfmstruct   Struct whose members will be freed.
 *
 */
void bah_rfm_precompute_free(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct) {
  // If params->is_host==true: rfmstruct, xx[] are HOST pointers; CPU path executes.
  // If params->is_host==false: rfmstruct, xx[] are DEVICE or UVM pointers; GPU path executes.

  // rfm_precompute_free: free rfmstruct arrays from host or device
  if (params->is_host) {
    BHAH_FREE__PtrMember(rfmstruct, f0_of_xx0);
    BHAH_FREE__PtrMember(rfmstruct, f1_of_xx1);
    BHAH_FREE__PtrMember(rfmstruct, f1_of_xx1__D1);
    BHAH_FREE__PtrMember(rfmstruct, f1_of_xx1__DD11);
  } else {
    IFCUDARUN({
      BHAH_FREE_DEVICE__PtrMember(rfmstruct, f0_of_xx0);
      BHAH_FREE_DEVICE__PtrMember(rfmstruct, f1_of_xx1);
      BHAH_FREE_DEVICE__PtrMember(rfmstruct, f1_of_xx1__D1);
      BHAH_FREE_DEVICE__PtrMember(rfmstruct, f1_of_xx1__DD11);
    });
  } // END IF params->is_host
} // END FUNCTION: bah_rfm_precompute_free
