#include "BHaH_defines.h"

/**
 * @file rfm_precompute_malloc.c
 * Allocates memory for precomputed reference metric quantities.
 *
 * - When params->is_host==true, perform HOST allocation.
 * - When params->is_host==false, perform DEVICE allocation.
 *
 * @param[in]  commondata  Global simulation metadata (unused).
 * @param[in]  params      Grid parameters (is_host, dims).
 * @param[out] rfmstruct   Destination struct for allocated pointers.
 *
 */
void bah_rfm_precompute_malloc(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct) {
  // If params->is_host==true: rfmstruct, xx[] are HOST pointers; CPU path executes.
  // If params->is_host==false: rfmstruct, xx[] are DEVICE or UVM pointers; GPU path executes.

  // rfm_precompute_malloc: allocate rfmstruct arrays on host or device
  if (params->is_host) {
    BHAH_MALLOC__PtrMember(rfmstruct, f0_of_xx0, sizeof(BHA_REAL) * ((size_t)params->Nxx_plus_2NGHOSTS0));
    BHAH_MALLOC__PtrMember(rfmstruct, f1_of_xx1, sizeof(BHA_REAL) * ((size_t)params->Nxx_plus_2NGHOSTS1));
    BHAH_MALLOC__PtrMember(rfmstruct, f1_of_xx1__D1, sizeof(BHA_REAL) * ((size_t)params->Nxx_plus_2NGHOSTS1));
    BHAH_MALLOC__PtrMember(rfmstruct, f1_of_xx1__DD11, sizeof(BHA_REAL) * ((size_t)params->Nxx_plus_2NGHOSTS1));
  } else {
    IFCUDARUN({
      BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f0_of_xx0, sizeof(BHA_REAL) * ((size_t)params->Nxx_plus_2NGHOSTS0));
      BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f1_of_xx1, sizeof(BHA_REAL) * ((size_t)params->Nxx_plus_2NGHOSTS1));
      BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f1_of_xx1__D1, sizeof(BHA_REAL) * ((size_t)params->Nxx_plus_2NGHOSTS1));
      BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f1_of_xx1__DD11, sizeof(BHA_REAL) * ((size_t)params->Nxx_plus_2NGHOSTS1));
    });
  } // END IF params->is_host
} // END FUNCTION: bah_rfm_precompute_malloc
