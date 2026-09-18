#include "BHaH_defines.h"

#ifdef __CUDACC__
__global__ static void rfm_precompute_defines__f0_of_xx0(const size_t N, rfm_struct *restrict d_rfm, const BHA_REAL *restrict dx0) {
  size_t thread_linear_index = (size_t)threadIdx.x + (size_t)blockIdx.x * (size_t)blockDim.x;
  size_t thread_linear_stride = (size_t)blockDim.x * (size_t)gridDim.x;
  for (size_t i0 = thread_linear_index; i0 < N; i0 += thread_linear_stride) {
    const BHA_REAL xx0 = dx0[i0];
    d_rfm->f0_of_xx0[i0] = xx0;
  }
}
#endif // __CUDACC__

#ifdef __CUDACC__
__global__ static void rfm_precompute_defines__f1_of_xx1(const size_t N, rfm_struct *restrict d_rfm, const BHA_REAL *restrict dx1) {
  size_t thread_linear_index = (size_t)threadIdx.x + (size_t)blockIdx.x * (size_t)blockDim.x;
  size_t thread_linear_stride = (size_t)blockDim.x * (size_t)gridDim.x;
  for (size_t i1 = thread_linear_index; i1 < N; i1 += thread_linear_stride) {
    const BHA_REAL xx1 = dx1[i1];
    d_rfm->f1_of_xx1[i1] = sin(xx1);
  }
}
#endif // __CUDACC__

#ifdef __CUDACC__
__global__ static void rfm_precompute_defines__f1_of_xx1__D1(const size_t N, rfm_struct *restrict d_rfm, const BHA_REAL *restrict dx1) {
  size_t thread_linear_index = (size_t)threadIdx.x + (size_t)blockIdx.x * (size_t)blockDim.x;
  size_t thread_linear_stride = (size_t)blockDim.x * (size_t)gridDim.x;
  for (size_t i1 = thread_linear_index; i1 < N; i1 += thread_linear_stride) {
    const BHA_REAL xx1 = dx1[i1];
    d_rfm->f1_of_xx1__D1[i1] = cos(xx1);
  }
}
#endif // __CUDACC__

#ifdef __CUDACC__
__global__ static void rfm_precompute_defines__f1_of_xx1__DD11(const size_t N, rfm_struct *restrict d_rfm, const BHA_REAL *restrict dx1) {
  size_t thread_linear_index = (size_t)threadIdx.x + (size_t)blockIdx.x * (size_t)blockDim.x;
  size_t thread_linear_stride = (size_t)blockDim.x * (size_t)gridDim.x;
  for (size_t i1 = thread_linear_index; i1 < N; i1 += thread_linear_stride) {
    const BHA_REAL xx1 = dx1[i1];
    d_rfm->f1_of_xx1__DD11[i1] = -sin(xx1);
  }
}
#endif // __CUDACC__

/**
 * @file rfm_precompute_defines.c
 * Computes and populates arrays with precomputed reference metric quantities.
 *
 * - params->is_host==true: rfmstruct & xx[] are HOST; run CPU loops.
 * - params->is_host==false: rfmstruct & xx[] are DEVICE/UVM; launch CUDA kernels.
 *
 * @param[in]  commondata  Global simulation metadata (unused).
 * @param[in]  params      Grid parameters and dimensions.
 * @param[out] rfmstruct   Struct containing arrays to populate.
 * @param[in]  xx          Pointers to 1D coordinate arrays.
 *
 */
void bah_rfm_precompute_defines(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                BHA_REAL *restrict xx[3]) {
  // If params->is_host==true: rfmstruct, xx[] are HOST pointers; CPU path executes.
  // If params->is_host==false: rfmstruct, xx[] are DEVICE or UVM pointers; GPU path executes.
  MAYBE_UNUSED const BHA_REAL *restrict x0 = xx[0];
  MAYBE_UNUSED const BHA_REAL *restrict x1 = xx[1];
  MAYBE_UNUSED const BHA_REAL *restrict x2 = xx[2];

  /* f0_of_xx0: 1D precompute */
  if (params->is_host) {
    {
      const size_t N = (size_t)params->Nxx_plus_2NGHOSTS0;
      for (size_t i0 = 0; i0 < N; i0++) {
        const BHA_REAL xx0 = x0[i0];
        rfmstruct->f0_of_xx0[i0] = xx0;
      }
    }
  } else {
    IFCUDARUN({
      const size_t N = (size_t)params->Nxx_plus_2NGHOSTS0;

      const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_DEFAULT;
      dim3 threads_per_block(threads_in_x_dir, 1, 1); // 1-D block to avoid duplicate writes
      dim3 blocks_per_grid(((N) + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
      size_t sm = 0;
      const size_t streamid = params->grid_idx % NUM_STREAMS;
      rfm_precompute_defines__f0_of_xx0<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(N, rfmstruct, x0);
      cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f0_of_xx0 failure");
    });
  }

  /* f1_of_xx1: 1D precompute */
  if (params->is_host) {
    {
      const size_t N = (size_t)params->Nxx_plus_2NGHOSTS1;
      for (size_t i1 = 0; i1 < N; i1++) {
        const BHA_REAL xx1 = x1[i1];
        rfmstruct->f1_of_xx1[i1] = sin(xx1);
      }
    }
  } else {
    IFCUDARUN({
      const size_t N = (size_t)params->Nxx_plus_2NGHOSTS1;

      const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_DEFAULT;
      dim3 threads_per_block(threads_in_x_dir, 1, 1); // 1-D block to avoid duplicate writes
      dim3 blocks_per_grid(((N) + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
      size_t sm = 0;
      const size_t streamid = params->grid_idx % NUM_STREAMS;
      rfm_precompute_defines__f1_of_xx1<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(N, rfmstruct, x1);
      cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f1_of_xx1 failure");
    });
  }

  /* f1_of_xx1__D1: 1D precompute */
  if (params->is_host) {
    {
      const size_t N = (size_t)params->Nxx_plus_2NGHOSTS1;
      for (size_t i1 = 0; i1 < N; i1++) {
        const BHA_REAL xx1 = x1[i1];
        rfmstruct->f1_of_xx1__D1[i1] = cos(xx1);
      }
    }
  } else {
    IFCUDARUN({
      const size_t N = (size_t)params->Nxx_plus_2NGHOSTS1;

      const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_DEFAULT;
      dim3 threads_per_block(threads_in_x_dir, 1, 1); // 1-D block to avoid duplicate writes
      dim3 blocks_per_grid(((N) + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
      size_t sm = 0;
      const size_t streamid = params->grid_idx % NUM_STREAMS;
      rfm_precompute_defines__f1_of_xx1__D1<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(N, rfmstruct, x1);
      cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f1_of_xx1__D1 failure");
    });
  }

  /* f1_of_xx1__DD11: 1D precompute */
  if (params->is_host) {
    {
      const size_t N = (size_t)params->Nxx_plus_2NGHOSTS1;
      for (size_t i1 = 0; i1 < N; i1++) {
        const BHA_REAL xx1 = x1[i1];
        rfmstruct->f1_of_xx1__DD11[i1] = -sin(xx1);
      }
    }
  } else {
    IFCUDARUN({
      const size_t N = (size_t)params->Nxx_plus_2NGHOSTS1;

      const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_DEFAULT;
      dim3 threads_per_block(threads_in_x_dir, 1, 1); // 1-D block to avoid duplicate writes
      dim3 blocks_per_grid(((N) + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
      size_t sm = 0;
      const size_t streamid = params->grid_idx % NUM_STREAMS;
      rfm_precompute_defines__f1_of_xx1__DD11<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(N, rfmstruct, x1);
      cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f1_of_xx1__DD11 failure");
    });
  }
} // END FUNCTION: bah_rfm_precompute_defines
