#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

// Indexing macros
#define EX_IDX4(g, i, j, k)                                                                                                                          \
  ((i) + commondata->external_input_Nxx_plus_2NGHOSTS0 *                                                                                             \
             ((j) + commondata->external_input_Nxx_plus_2NGHOSTS1 * ((k) + commondata->external_input_Nxx_plus_2NGHOSTS2 * (g))))
#define EX_NOGZ_IDX4(g, i, j, k)                                                                                                                     \
  ((i) + bhahaha_params_and_data->Nr_external_input * ((j) + commondata->external_input_Nxx1 * ((k) + commondata->external_input_Nxx2 * (g))))

/**
 * Initializes and processes the external input metric data (gamma_{ij}, K_{ij}).
 *
 * This function performs the following steps:
 * 1. Unpacks input parameters from the common data structure.
 * 2. Adds ghost zones to the external input arrays to facilitate boundary condition application.
 * 3. Allocates memory for the external input gridfunctions with ghost zones and assigns it to the common data structure.
 * 4. Transfers metric data from arrays without ghost zones into the newly allocated arrays with ghost zones.
 * 5. Sets up coordinate arrays for a uniform, cell-centered spherical grid.
 * 6. Transforms the metric components (gamma_{ij}, K_{ij}) from Cartesian to spherical coordinates, including necessary rescaling.
 * 7. Sets up boundary condition structures and applies inner boundary conditions, including parity corrections for all gridfunctions.
 *
 * @param commondata - Pointer to the common data structure containing simulation parameters and data.
 * @param n_resolutions - Number of angular resolutions.
 * @param Ntheta - Array containing the number of theta points for each resolution.
 * @param Nphi - Array containing the number of phi points for each resolution.
 *
 * @return BHAHAHA_SUCCESS on successful setup, or an error code indicating the failure reason.
 *
 */
int bah_numgrid__external_input_set_up(commondata_struct *restrict commondata, const int n_resolutions, const int *restrict Ntheta,
                                       const int *restrict Nphi) {

  // Step 1: Unpack input parameters from the common data structure.
  const bhahaha_params_and_data_struct *restrict bhahaha_params_and_data = commondata->bhahaha_params_and_data;

  // Calculate the number of interior (non-ghost) radial points by subtracting ghost zones.
  // Nr from external includes r ~ r_max NGHOSTS.
  commondata->external_input_Nxx0 = bhahaha_params_and_data->Nr_external_input - NGHOSTS;
  if (bhahaha_params_and_data->r_min_external_input > 0) {
    commondata->external_input_Nxx0 = bhahaha_params_and_data->Nr_external_input - 2 * NGHOSTS;
  }
  int i0_min_shift = 0;
  if (bhahaha_params_and_data->r_min_external_input == 0)
    i0_min_shift = NGHOSTS;

  // Set fixed angular resolutions for theta and phi directions.
  {
    const int max_resolution_i = bhahaha_params_and_data->num_resolutions_multigrid - 1;
    commondata->external_input_Nxx1 = bhahaha_params_and_data->Ntheta_array_multigrid[max_resolution_i];
    commondata->external_input_Nxx2 = bhahaha_params_and_data->Nphi_array_multigrid[max_resolution_i];
  }

  // Step 1.a: Calculate grid spacing in each coordinate direction based on the simulation domain and resolution.
  // x_i = min_i + (j + 0.5) * dx_i, where dx_i = (max_i - min_i) / N_i

  commondata->external_input_dxx0 = bhahaha_params_and_data->dr_external_input;
  commondata->external_input_dxx1 = M_PI / ((BHA_REAL)commondata->external_input_Nxx1);
  commondata->external_input_dxx2 = 2 * M_PI / ((BHA_REAL)commondata->external_input_Nxx2);

  // Precompute inverse grid spacings for performance optimization in calculations.
  commondata->external_input_invdxx0 = 1.0 / commondata->external_input_dxx0;
  commondata->external_input_invdxx1 = 1.0 / commondata->external_input_dxx1;
  commondata->external_input_invdxx2 = 1.0 / commondata->external_input_dxx2;

  // Step 2: Add ghost zones to the external input data.
  // Ghost zones are added to so that inner boundary conditions may be applied; 2 * NGHOSTS in each angular direction and NGHOSTS in the radial
  // direction.
  commondata->external_input_Nxx_plus_2NGHOSTS0 = commondata->external_input_Nxx0 + 2 * NGHOSTS;
  commondata->external_input_Nxx_plus_2NGHOSTS1 = commondata->external_input_Nxx1 + 2 * NGHOSTS;
  commondata->external_input_Nxx_plus_2NGHOSTS2 = commondata->external_input_Nxx2 + 2 * NGHOSTS;

  // Calculate the total number of grid points including ghost zones.
  const int total_elements_incl_gzs =
      commondata->external_input_Nxx_plus_2NGHOSTS0 * commondata->external_input_Nxx_plus_2NGHOSTS1 * commondata->external_input_Nxx_plus_2NGHOSTS2;

  // Pointers to the external input gridfunctions without ghost zones.
  BHA_REAL *restrict external_input_gfs_no_gzs = commondata->external_input_gfs_Cart_basis_no_gzs;

  // Allocate memory for the external input gridfunctions with ghost zones.
  BHA_REAL *restrict external_input_gfs = (BHA_REAL *)malloc(NUM_EXT_INPUT_CONFORMAL_GFS * total_elements_incl_gzs * sizeof(BHA_REAL));
  if (external_input_gfs == NULL) {
    return NUMGRID_EXTERN_MALLOC_ERROR_GFS;
  } // END IF memory allocation for external_input_gfs failed

  // Step 3: Assign the allocated array to commondata for use outside this function.
  commondata->external_input_gfs = external_input_gfs;

  // Step 4: Transfer data from the no-ghost zones array to the array with ghost zones.
  // This involves copying metric data (gamma_{ij} and K_{ij}) from the Cartesian basis into the newly allocated arrays with ghost zones.
  LOOP_OMP("omp parallel for",                                //
           i0, 0, bhahaha_params_and_data->Nr_external_input, //
           i1, 0, commondata->external_input_Nxx1,            //
           i2, 0, commondata->external_input_Nxx2) {
    for (int gf = 0; gf < NUM_EXT_INPUT_CARTESIAN_GFS; gf++) {
      external_input_gfs[EX_IDX4(gf, i0 + i0_min_shift, i1 + NGHOSTS, i2 + NGHOSTS)] = external_input_gfs_no_gzs[EX_NOGZ_IDX4(gf, i0, i1, i2)];
    }
  } // END LOOP: iterating through the external input grid points

  // Step 5: Set up coordinate arrays for a uniform, cell-centered spherical grid.
  {
    const int Nxx_plus_2NGHOSTS0 = commondata->external_input_Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = commondata->external_input_Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = commondata->external_input_Nxx_plus_2NGHOSTS2;

    // Step 5.a: Allocate memory for coordinate arrays in radial, theta, and phi directions.
    commondata->external_input_r_theta_phi[0] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * Nxx_plus_2NGHOSTS0);
    commondata->external_input_r_theta_phi[1] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * Nxx_plus_2NGHOSTS1);
    commondata->external_input_r_theta_phi[2] = (BHA_REAL *)malloc(sizeof(BHA_REAL) * Nxx_plus_2NGHOSTS2);
    if (commondata->external_input_r_theta_phi[0] == NULL || commondata->external_input_r_theta_phi[1] == NULL ||
        commondata->external_input_r_theta_phi[2] == NULL) {
      free(external_input_gfs);
      return NUMGRID_EXTERN_MALLOC_ERROR_RTHETAPHI;
    } // END IF memory allocation for external_input_r_theta_phi arrays failed

    // Step 5.b: Initialize coordinate arrays for a uniform, cell-centered spherical grid.
    // The coordinates are centered within each cell by adding 0.5 to the index before scaling.
    const BHA_REAL xxmin0 = bhahaha_params_and_data->r_min_external_input;
    const BHA_REAL xxmin1 = 0.0;
    const BHA_REAL xxmin2 = -M_PI;

    for (int j = 0; j < Nxx_plus_2NGHOSTS0; j++)
      commondata->external_input_r_theta_phi[0][j] = xxmin0 + ((BHA_REAL)(j - NGHOSTS) + (1.0 / 2.0)) * commondata->external_input_dxx0;
    for (int j = 0; j < Nxx_plus_2NGHOSTS1; j++)
      commondata->external_input_r_theta_phi[1][j] = xxmin1 + ((BHA_REAL)(j - NGHOSTS) + (1.0 / 2.0)) * commondata->external_input_dxx1;
    for (int j = 0; j < Nxx_plus_2NGHOSTS2; j++)
      commondata->external_input_r_theta_phi[2][j] = xxmin2 + ((BHA_REAL)(j - NGHOSTS) + (1.0 / 2.0)) * commondata->external_input_dxx2;
  } // END BLOCK: setting up coordinate arrays

  // Step 6: Transform the metric components (gamma_{ij}, K_{ij}) from Cartesian to spherical coordinates,
  // including necessary rescaling.
  {
    const int Nxx_plus_2NGHOSTS0 = commondata->external_input_Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = commondata->external_input_Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = commondata->external_input_Nxx_plus_2NGHOSTS2;

    // Step 6.a: Extract coordinate arrays for easier access during transformation.
    BHA_REAL *restrict external_input_r_theta_phi[3];
    for (int ww = 0; ww < 3; ww++)
      external_input_r_theta_phi[ww] = commondata->external_input_r_theta_phi[ww];

    // Step 6.b: Metric components: basis transform from Cartesian to Spherical & convert ADM->rescaled BSSN.
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < commondata->external_input_Nxx2 + NGHOSTS; i2++) {
      const BHA_REAL xx2 = external_input_r_theta_phi[2][i2];
      for (int i1 = NGHOSTS; i1 < commondata->external_input_Nxx1 + NGHOSTS; i1++) {
        const BHA_REAL xx1 = external_input_r_theta_phi[1][i1];
        // Include all valid points, including those near r ~ r_max.
        for (int i0 = i0_min_shift; i0 < commondata->external_input_Nxx_plus_2NGHOSTS0; i0++) {
          const BHA_REAL xx0 = external_input_r_theta_phi[0][i0];

          // Read Cartesian metric components at the current grid point.
          const BHA_REAL Cart_gammaDD00 = external_input_gfs[IDX4(INTERP_GAMMADDXXGF, i0, i1, i2)];
          const BHA_REAL Cart_gammaDD01 = external_input_gfs[IDX4(INTERP_GAMMADDXYGF, i0, i1, i2)];
          const BHA_REAL Cart_gammaDD02 = external_input_gfs[IDX4(INTERP_GAMMADDXZGF, i0, i1, i2)];
          const BHA_REAL Cart_gammaDD11 = external_input_gfs[IDX4(INTERP_GAMMADDYYGF, i0, i1, i2)];
          const BHA_REAL Cart_gammaDD12 = external_input_gfs[IDX4(INTERP_GAMMADDYZGF, i0, i1, i2)];
          const BHA_REAL Cart_gammaDD22 = external_input_gfs[IDX4(INTERP_GAMMADDZZGF, i0, i1, i2)];
          const BHA_REAL Cart_KDD00 = external_input_gfs[IDX4(INTERP_KDDXXGF, i0, i1, i2)];
          const BHA_REAL Cart_KDD01 = external_input_gfs[IDX4(INTERP_KDDXYGF, i0, i1, i2)];
          const BHA_REAL Cart_KDD02 = external_input_gfs[IDX4(INTERP_KDDXZGF, i0, i1, i2)];
          const BHA_REAL Cart_KDD11 = external_input_gfs[IDX4(INTERP_KDDYYGF, i0, i1, i2)];
          const BHA_REAL Cart_KDD12 = external_input_gfs[IDX4(INTERP_KDDYZGF, i0, i1, i2)];
          const BHA_REAL Cart_KDD22 = external_input_gfs[IDX4(INTERP_KDDZZGF, i0, i1, i2)];
          const BHA_REAL tmp0 = ((xx0) * (xx0) * (xx0) * (xx0));
          const BHA_REAL tmp1 = sin(xx1);
          const BHA_REAL tmp4 = sin(xx2);
          const BHA_REAL tmp5 = cos(xx1);
          const BHA_REAL tmp10 = cos(xx2);
          const BHA_REAL tmp13 = Cart_gammaDD01 * xx0;
          const BHA_REAL tmp22 = ((xx0) * (xx0));
          const BHA_REAL tmp35 = 2 * Cart_gammaDD01;
          const BHA_REAL tmp75 = 2 * Cart_KDD01;
          const BHA_REAL tmp98 = (2.0 / 3.0) * Cart_gammaDD01;
          const BHA_REAL tmp2 = ((tmp1) * (tmp1));
          const BHA_REAL tmp6 = tmp1 * tmp5;
          const BHA_REAL tmp11 = ((tmp4) * (tmp4));
          const BHA_REAL tmp16 = tmp10 * tmp4;
          const BHA_REAL tmp19 = ((tmp10) * (tmp10));
          const BHA_REAL tmp25 = Cart_gammaDD02 * tmp10;
          const BHA_REAL tmp28 = Cart_gammaDD12 * tmp4;
          const BHA_REAL tmp29 = ((tmp5) * (tmp5));
          const BHA_REAL tmp72 = Cart_KDD02 * tmp10;
          const BHA_REAL tmp73 = Cart_KDD12 * tmp4;
          const BHA_REAL tmp91 = (1.0 / (tmp1));
          const BHA_REAL tmp93 = (1.0 / (tmp22));
          const BHA_REAL tmp3 = (1.0 / (tmp2));
          const BHA_REAL tmp7 = tmp6 * xx0;
          const BHA_REAL tmp12 = tmp11 * tmp2;
          const BHA_REAL tmp15 = tmp2 * xx0;
          const BHA_REAL tmp20 = tmp19 * tmp2;
          const BHA_REAL tmp23 = tmp2 * tmp22;
          const BHA_REAL tmp26 = tmp22 * tmp6;
          const BHA_REAL tmp30 = tmp22 * tmp29;
          const BHA_REAL tmp39 = 2 * tmp6;
          const BHA_REAL tmp42 = tmp16 * tmp2;
          const BHA_REAL tmp53 = tmp29 * xx0;
          const BHA_REAL tmp60 = tmp13 * tmp16;
          const BHA_REAL tmp97 = (2.0 / 3.0) * tmp6;
          const BHA_REAL tmp17 = tmp15 * tmp16;
          const BHA_REAL tmp27 = 2 * tmp26;
          const BHA_REAL tmp32 = Cart_gammaDD00 * tmp19 * tmp30;
          const BHA_REAL tmp36 = tmp16 * tmp30;
          const BHA_REAL tmp44 = tmp10 * tmp23;
          const BHA_REAL tmp48 = tmp16 * tmp26;
          const BHA_REAL tmp51 = Cart_gammaDD22 * tmp7;
          const BHA_REAL tmp52 = tmp15 * tmp25;
          const BHA_REAL tmp54 = tmp25 * tmp53;
          const BHA_REAL tmp55 = tmp15 * tmp28;
          const BHA_REAL tmp56 = tmp28 * tmp53;
          const BHA_REAL tmp57 = Cart_gammaDD00 * tmp19 * tmp7;
          const BHA_REAL tmp82 = Cart_KDD00 * tmp20 + Cart_KDD11 * tmp12 + Cart_KDD22 * tmp29 + tmp39 * tmp72 + tmp39 * tmp73 + tmp42 * tmp75;
          const BHA_REAL tmp9 = Cart_gammaDD02 * tmp4 * tmp7;
          const BHA_REAL tmp34 = Cart_gammaDD11 * tmp11 * tmp30;
          const BHA_REAL tmp43 = Cart_gammaDD00 * tmp20 + Cart_gammaDD11 * tmp12 + Cart_gammaDD22 * tmp29 + tmp25 * tmp39 + tmp28 * tmp39 + tmp35 * tmp42;
          const BHA_REAL tmp47 = Cart_gammaDD01 * tmp11 * tmp26;
          const BHA_REAL tmp59 = Cart_gammaDD11 * tmp11 * tmp7;
          const BHA_REAL tmp63 = Cart_gammaDD00 * tmp11 * tmp23;
          const BHA_REAL tmp65 = Cart_gammaDD11 * tmp19 * tmp23;
          const BHA_REAL tmp66 = tmp4 * tmp44;
          const BHA_REAL tmp76 =
              Cart_KDD00 * tmp19 * tmp30 + Cart_KDD11 * tmp11 * tmp30 + Cart_KDD22 * tmp23 - tmp27 * tmp72 - tmp27 * tmp73 + tmp36 * tmp75;
          const BHA_REAL tmp84 = Cart_KDD00 * tmp19 * tmp7 + Cart_KDD11 * tmp11 * tmp7 - Cart_KDD22 * tmp7 - tmp15 * tmp72 - tmp15 * tmp73 +
                             tmp16 * tmp7 * tmp75 + tmp53 * tmp72 + tmp53 * tmp73;
          const BHA_REAL tmp86 = Cart_KDD00 * tmp48 - Cart_KDD01 * tmp1 * tmp19 * tmp22 * tmp5 + Cart_KDD01 * tmp11 * tmp26 -
                             Cart_KDD02 * tmp2 * tmp22 * tmp4 - Cart_KDD11 * tmp1 * tmp10 * tmp22 * tmp4 * tmp5 + Cart_KDD12 * tmp44;
          const BHA_REAL tmp87 = Cart_KDD00 * tmp17 + Cart_KDD01 * tmp12 * xx0 - Cart_KDD01 * tmp19 * tmp2 * xx0 + Cart_KDD02 * tmp4 * tmp7 -
                             Cart_KDD11 * tmp10 * tmp2 * tmp4 * xx0 - Cart_KDD12 * tmp1 * tmp10 * tmp5 * xx0;
          const BHA_REAL tmp21 = Cart_gammaDD00 * tmp17 - Cart_gammaDD11 * tmp17 - Cart_gammaDD12 * tmp10 * tmp7 + tmp12 * tmp13 - tmp13 * tmp20 + tmp9;
          const BHA_REAL tmp37 = Cart_gammaDD22 * tmp23 - tmp25 * tmp27 - tmp27 * tmp28 + tmp32 + tmp34 + tmp35 * tmp36;
          const BHA_REAL tmp50 = Cart_gammaDD00 * tmp48 - Cart_gammaDD01 * tmp19 * tmp26 - Cart_gammaDD02 * tmp23 * tmp4 - Cart_gammaDD11 * tmp48 +
                             Cart_gammaDD12 * tmp44 + tmp47;
          const BHA_REAL tmp61 = -2 * tmp51 - 2 * tmp52 + 2 * tmp54 - 2 * tmp55 + 2 * tmp56 + 2 * tmp57 + 2 * tmp59 + 4 * tmp6 * tmp60;
          const BHA_REAL tmp67 = -tmp35 * tmp66 + tmp63 + tmp65;
          const BHA_REAL tmp68 = tmp39 * tmp60 - tmp51 - tmp52 + tmp54 - tmp55 + tmp56 + tmp57 + tmp59;
          const BHA_REAL tmp83 = Cart_KDD00 * tmp11 * tmp23 + Cart_KDD11 * tmp19 * tmp23 - tmp66 * tmp75;
          const BHA_REAL tmp71 = -tmp37 * tmp43 * tmp67 + tmp67 * ((tmp68) * (tmp68));
          const BHA_REAL tmp81 = (1.0 / (-((tmp21) * (tmp21)) * tmp37 + tmp21 * tmp50 * tmp61 - tmp43 * ((tmp50) * (tmp50)) - tmp71));
          const BHA_REAL tmp85 = 2 * tmp81;
          const BHA_REAL tmp89 = cbrt(tmp0 * tmp2 * tmp81);
          const BHA_REAL tmp88 = tmp76 * tmp81 * (-((tmp21) * (tmp21)) + tmp43 * tmp67) + tmp81 * tmp82 * (tmp37 * tmp67 - ((tmp50) * (tmp50))) +
                             tmp81 * tmp83 * (tmp37 * tmp43 - ((tmp68) * (tmp68))) + tmp84 * tmp85 * (tmp21 * tmp50 - tmp67 * tmp68) -
                             tmp85 * tmp86 * (-tmp21 * tmp68 + tmp43 * tmp50) - tmp85 * tmp87 * (tmp21 * tmp37 - tmp50 * tmp68);
          const BHA_REAL tmp90 = tmp89 / xx0;
          const BHA_REAL tmp95 = tmp89 * tmp91 * tmp93;
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_WWGF, i0, i1, i2)] =
              (1.0 / (sqrt(cbrt(fabs(tmp3 * (((tmp21) * (tmp21)) * tmp37 - tmp21 * tmp50 * tmp61 + tmp43 * ((tmp50) * (tmp50)) + tmp71) / tmp0)))));
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_TRKGF, i0, i1, i2)] = tmp88;
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_HDD00GF, i0, i1, i2)] = tmp43 * tmp89 - 1;
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_HDD01GF, i0, i1, i2)] = tmp68 * tmp90;
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_HDD02GF, i0, i1, i2)] = -tmp21 * tmp90 * tmp91;
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_HDD11GF, i0, i1, i2)] = tmp93 * (-tmp22 + tmp37 * tmp89);
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_HDD12GF, i0, i1, i2)] = -tmp50 * tmp95;
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_HDD22GF, i0, i1, i2)] = tmp3 * tmp93 * (-tmp23 + tmp67 * tmp89);
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_ADD00GF, i0, i1, i2)] =
              tmp89 * (tmp82 - tmp88 * ((1.0 / 3.0) * Cart_gammaDD00 * tmp20 + (1.0 / 3.0) * Cart_gammaDD11 * tmp12 +
                                        (1.0 / 3.0) * Cart_gammaDD22 * tmp29 + tmp25 * tmp97 + tmp28 * tmp97 + tmp42 * tmp98));
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_ADD01GF, i0, i1, i2)] =
              tmp90 * (tmp84 - tmp88 * (-1.0 / 3.0 * tmp51 - 1.0 / 3.0 * tmp52 + (1.0 / 3.0) * tmp54 - 1.0 / 3.0 * tmp55 + (1.0 / 3.0) * tmp56 +
                                        (1.0 / 3.0) * tmp57 + (1.0 / 3.0) * tmp59 + tmp60 * tmp97));
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_ADD02GF, i0, i1, i2)] =
              tmp90 * tmp91 *
              (-tmp87 - tmp88 * (-1.0 / 3.0 * Cart_gammaDD00 * tmp17 + (1.0 / 3.0) * Cart_gammaDD01 * tmp19 * tmp2 * xx0 +
                                 (1.0 / 3.0) * Cart_gammaDD11 * tmp10 * tmp2 * tmp4 * xx0 + (1.0 / 3.0) * Cart_gammaDD12 * tmp1 * tmp10 * tmp5 * xx0 -
                                 1.0 / 3.0 * tmp12 * tmp13 - 1.0 / 3.0 * tmp9));
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_ADD11GF, i0, i1, i2)] =
              tmp89 * tmp93 *
              (tmp76 - tmp88 * ((1.0 / 3.0) * Cart_gammaDD22 * tmp23 - 2.0 / 3.0 * tmp25 * tmp26 - 2.0 / 3.0 * tmp26 * tmp28 + (1.0 / 3.0) * tmp32 +
                                (1.0 / 3.0) * tmp34 + tmp36 * tmp98));
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_ADD12GF, i0, i1, i2)] =
              tmp95 * (-tmp86 - tmp88 * (-1.0 / 3.0 * Cart_gammaDD00 * tmp48 + (1.0 / 3.0) * Cart_gammaDD01 * tmp1 * tmp19 * tmp22 * tmp5 +
                                         (1.0 / 3.0) * Cart_gammaDD02 * tmp2 * tmp22 * tmp4 +
                                         (1.0 / 3.0) * Cart_gammaDD11 * tmp1 * tmp10 * tmp22 * tmp4 * tmp5 - 1.0 / 3.0 * Cart_gammaDD12 * tmp44 -
                                         1.0 / 3.0 * tmp47));
          external_input_gfs[IDX4(EXTERNAL_SPHERICAL_ADD22GF, i0, i1, i2)] =
              tmp3 * tmp89 * tmp93 * (tmp83 - tmp88 * ((1.0 / 3.0) * tmp63 + (1.0 / 3.0) * tmp65 - tmp66 * tmp98));

        } // END LOOP over i0
      } // END LOOP over i1
    } // END LOOP over i2
  } // END BLOCK: transformation and rescaling

  // Step 7: Set up boundary condition structures and apply inner boundary conditions.
  {
    const int Nxx_plus_2NGHOSTS0 = commondata->external_input_Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = commondata->external_input_Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = commondata->external_input_Nxx_plus_2NGHOSTS2;

    // Assign grid spacings and sizes to the boundary condition structure within commondata.
    commondata->bcstruct_dxx0 = commondata->external_input_dxx0;
    commondata->bcstruct_dxx1 = commondata->external_input_dxx1;
    commondata->bcstruct_dxx2 = commondata->external_input_dxx2;

    commondata->bcstruct_Nxx_plus_2NGHOSTS0 = Nxx_plus_2NGHOSTS0;
    commondata->bcstruct_Nxx_plus_2NGHOSTS1 = Nxx_plus_2NGHOSTS1;
    commondata->bcstruct_Nxx_plus_2NGHOSTS2 = Nxx_plus_2NGHOSTS2;

    // Initialize the boundary condition structure for external input data.
    bc_struct external_input_bcstruct;
    bah_bcstruct_set_up(commondata, commondata->external_input_r_theta_phi, &external_input_bcstruct);

    // Step 7.a: Unpack boundary condition information from the boundary condition structure.
    const bc_info_struct *restrict bc_info = &external_input_bcstruct.bc_info;

    // Step 7.b: Apply inner boundary conditions to all gridfunctions.
    // This involves copying values from source points to destination points with parity corrections.
#pragma omp parallel
    for (int which_gf = 0; which_gf < NUM_EXT_INPUT_CONFORMAL_GFS; which_gf++) {
#pragma omp for
      for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
        const int dstpt = external_input_bcstruct.inner_bc_array[pt].dstpt;
        const int srcpt = external_input_bcstruct.inner_bc_array[pt].srcpt;
        // Apply the boundary condition by copying values from the source point to the destination point,
        // applying the appropriate parity correction for the gridfunction.
        commondata->external_input_gfs[IDX4pt(which_gf, dstpt)] =
            external_input_bcstruct.inner_bc_array[pt].parity[external_input_gf_parity[which_gf]] *
            commondata->external_input_gfs[IDX4pt(which_gf, srcpt)];
      } // END LOOP over inner boundary points
    } // END LOOP over gridfunctions

    // Step 7.c: Free allocated memory for boundary condition structures to prevent memory leaks.
    free(external_input_bcstruct.inner_bc_array);
    for (int ng = 0; ng < NGHOSTS * 3; ng++)
      free(external_input_bcstruct.pure_outer_bc_array[ng]);
  } // END BLOCK: applying boundary conditions

  return BHAHAHA_SUCCESS;
} // END FUNCTION bah_numgrid__external_input_set_up
