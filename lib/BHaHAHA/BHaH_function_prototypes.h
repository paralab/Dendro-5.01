void bah_apply_bcs_inner_only(const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct,
                              BHA_REAL *restrict gfs);
void bah_apply_bcs_r_maxmin_partial_r_hDD_upwinding(const commondata_struct *restrict commondata, BHA_REAL *restrict xx[3], BHA_REAL *restrict gfs,
                                                    const bool fill_r_min_ghosts);
int bah_bcstruct_set_up(const commondata_struct *restrict commondata, BHA_REAL *restrict xx[3], bc_struct *restrict bcstruct);
void bah_cfl_limited_timestep_based_on_h_equals_r(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
void bah_commondata_struct_set_to_default(commondata_struct *restrict commondata);
void bah_diagnostics(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
void bah_diagnostics_area_centroid_and_Theta_norms(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
void bah_diagnostics_file_output(const bhahaha_diagnostics_struct *diags, const bhahaha_params_and_data_struct *bhahaha_params_and_data,
                                 int N_horizons, const BHA_REAL x_center_input, const BHA_REAL y_center_input, const BHA_REAL z_center_input,
                                 const char *output_directory);
void bah_diagnostics_integration_weights(int Nxx1, int Nxx2, const BHA_REAL *restrict *weights, int *restrict weight_stencil_size);
void bah_diagnostics_min_max_mean_radii_wrt_centroid(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
int bah_diagnostics_proper_circumferences(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
const char *bah_error_message(const bhahaha_error_codes error_code);
int bah_find_horizon(bhahaha_params_and_data_struct *restrict bhahaha_params_and_data, bhahaha_diagnostics_struct *restrict bhahaha_diags);
void bah_hDD_dD_and_W_dD_in_interp_src_grid_interior(commondata_struct *restrict commondata);
int bah_initial_data(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
int bah_interpolation_1d_radial_spokes_on_3d_src_grid(const params_struct *restrict params, const commondata_struct *restrict commondata,
                                                      const BHA_REAL *restrict dst_radii_aka_src_h_gf, BHA_REAL *restrict dst_interp_gfs);
int bah_interpolation_2d_external_input_to_interp_src_grid(commondata_struct *restrict commondata);
int bah_interpolation_2d_general__uniform_src_grid(const int N_interp_GHOSTS, const BHA_REAL src_dxx1, const BHA_REAL src_dxx2,
                                                   const int src_Nxx_plus_2NGHOSTS1, const int src_Nxx_plus_2NGHOSTS2,
                                                   BHA_REAL *restrict src_r_theta_phi[3], const BHA_REAL *restrict src_gf, const int num_dst_pts,
                                                   const BHA_REAL dst_pts[][2], BHA_REAL *restrict dst_data);
void bah_KO_apply(const commondata_struct *restrict commondata, const params_struct *restrict params, const rfm_struct *restrict rfmstruct,
                  const BHA_REAL *restrict auxevol_gfs, const BHA_REAL *restrict in_gfs, BHA_REAL *restrict rhs_gfs);
void bah_MoL_free_memory_non_y_n_gfs(MoL_gridfunctions_struct *restrict gridfuncs);
void bah_MoL_free_memory_y_n_gfs(MoL_gridfunctions_struct *restrict gridfuncs);
void bah_MoL_malloc_non_y_n_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                MoL_gridfunctions_struct *restrict gridfuncs);
void bah_MoL_malloc_y_n_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                            MoL_gridfunctions_struct *restrict gridfuncs);
void bah_MoL_step_forward_in_time(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
void bah_numgrid__evol_set_up(commondata_struct *restrict commondata, griddata_struct *restrict griddata, const int Nx_evol_grid[3]);
int bah_numgrid__external_input_set_up(commondata_struct *restrict commondata, const int n_resolutions, const int *restrict Ntheta,
                                       const int *restrict Nphi);
int bah_numgrid__interp_src_set_up(commondata_struct *restrict commondata, const int Nx_evol_grid[3]);
void bah_over_relaxation(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
void bah_params_struct_set_to_default(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
void bah_poisoning_check_inputs(const bhahaha_params_and_data_struct *restrict params);
void bah_poisoning_set_inputs(bhahaha_params_and_data_struct *restrict params);
BHA_REAL bah_quadratic_extrapolation(const BHA_REAL times[3], const BHA_REAL y_tm1, const BHA_REAL y_tm2, const BHA_REAL y_tm3, const BHA_REAL dst_time);
void bah_radial_grid_cell_centered_set_up(const int Nr_interp_max, const BHA_REAL max_search_radius, const BHA_REAL input_r_min, const BHA_REAL input_r_max,
                                          int *restrict output_Nr_interp, BHA_REAL *restrict output_r_min_interior, BHA_REAL *restrict output_dr,
                                          BHA_REAL radii[Nr_interp_max]);
void bah_rfm_precompute_defines(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                BHA_REAL *restrict xx[3]);
void bah_rfm_precompute_free(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct);
void bah_rfm_precompute_malloc(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct);
void bah_rhs_eval(const commondata_struct *restrict commondata, const params_struct *restrict params, const rfm_struct *restrict rfmstruct,
                  const BHA_REAL *restrict auxevol_gfs, const BHA_REAL *restrict in_gfs, BHA_REAL *restrict rhs_gfs);
void bah_xx_to_Cart(const params_struct *restrict params, BHA_REAL xx[3], BHA_REAL xCart[3]);
void bah_xyz_center_r_minmax(const bhahaha_params_and_data_struct *restrict pars, BHA_REAL *restrict x_center, BHA_REAL *restrict y_center,
                             BHA_REAL *restrict z_center, BHA_REAL *restrict r_min, BHA_REAL *restrict r_max);