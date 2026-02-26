MAYBE_UNUSED const BHA_REAL Cart_originx = params->Cart_originx; // nrpy.grid::Cart_originx
MAYBE_UNUSED const BHA_REAL Cart_originy = params->Cart_originy; // nrpy.grid::Cart_originy
MAYBE_UNUSED const BHA_REAL Cart_originz = params->Cart_originz; // nrpy.grid::Cart_originz
MAYBE_UNUSED const BHA_REAL CFL_FACTOR = commondata->CFL_FACTOR; // nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all::CFL_FACTOR
char CoordSystemName[100];                                   // nrpy.reference_metric::CoordSystemName
{
  // Safely copy string with snprintf, which guarantees null termination
  snprintf(CoordSystemName, sizeof(CoordSystemName), "%s", params->CoordSystemName);
}
MAYBE_UNUSED const BHA_REAL dt = commondata->dt;                   // nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all::dt
MAYBE_UNUSED const BHA_REAL dxx0 = params->dxx0;                   // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::dxx0
MAYBE_UNUSED const BHA_REAL dxx1 = params->dxx1;                   // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::dxx1
MAYBE_UNUSED const BHA_REAL dxx2 = params->dxx2;                   // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::dxx2
MAYBE_UNUSED const BHA_REAL eta_damping = commondata->eta_damping; // nrpydev.infrastructures.BHaH.BHaHAHA.rhs_eval_KO_apply::eta_damping
MAYBE_UNUSED const BHA_REAL external_input_dxx0 =
    commondata->external_input_dxx0; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_dxx0
MAYBE_UNUSED const BHA_REAL external_input_dxx1 =
    commondata->external_input_dxx1; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_dxx1
MAYBE_UNUSED const BHA_REAL external_input_dxx2 =
    commondata->external_input_dxx2; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_dxx2
MAYBE_UNUSED const BHA_REAL external_input_invdxx0 =
    commondata->external_input_invdxx0; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_invdxx0
MAYBE_UNUSED const BHA_REAL external_input_invdxx1 =
    commondata->external_input_invdxx1; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_invdxx1
MAYBE_UNUSED const BHA_REAL external_input_invdxx2 =
    commondata->external_input_invdxx2; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_invdxx2
MAYBE_UNUSED const int external_input_Nxx0 =
    commondata->external_input_Nxx0; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_Nxx0
MAYBE_UNUSED const int external_input_Nxx1 =
    commondata->external_input_Nxx1; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_Nxx1
MAYBE_UNUSED const int external_input_Nxx2 =
    commondata->external_input_Nxx2; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_Nxx2
MAYBE_UNUSED const int external_input_Nxx_plus_2NGHOSTS0 =
    commondata
        ->external_input_Nxx_plus_2NGHOSTS0; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_Nxx_plus_2NGHOSTS0
MAYBE_UNUSED const int external_input_Nxx_plus_2NGHOSTS1 =
    commondata
        ->external_input_Nxx_plus_2NGHOSTS1; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_Nxx_plus_2NGHOSTS1
MAYBE_UNUSED const int external_input_Nxx_plus_2NGHOSTS2 =
    commondata
        ->external_input_Nxx_plus_2NGHOSTS2; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__external_input_set_up::external_input_Nxx_plus_2NGHOSTS2
MAYBE_UNUSED const BHA_REAL grid_hole_radius = params->grid_hole_radius;     // nrpy.reference_metric::grid_hole_radius
MAYBE_UNUSED const BHA_REAL grid_physical_size = params->grid_physical_size; // nrpy.reference_metric::grid_physical_size
MAYBE_UNUSED const bool grid_rotates = params->grid_rotates;             // nrpy.grid::grid_rotates
MAYBE_UNUSED const BHA_REAL interp_src_dxx0 =
    commondata->interp_src_dxx0; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_dxx0
MAYBE_UNUSED const BHA_REAL interp_src_dxx1 =
    commondata->interp_src_dxx1; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_dxx1
MAYBE_UNUSED const BHA_REAL interp_src_dxx2 =
    commondata->interp_src_dxx2; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_dxx2
MAYBE_UNUSED const BHA_REAL interp_src_invdxx0 =
    commondata->interp_src_invdxx0; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_invdxx0
MAYBE_UNUSED const BHA_REAL interp_src_invdxx1 =
    commondata->interp_src_invdxx1; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_invdxx1
MAYBE_UNUSED const BHA_REAL interp_src_invdxx2 =
    commondata->interp_src_invdxx2; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_invdxx2
MAYBE_UNUSED const int interp_src_Nxx0 =
    commondata->interp_src_Nxx0; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_Nxx0
MAYBE_UNUSED const int interp_src_Nxx1 =
    commondata->interp_src_Nxx1; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_Nxx1
MAYBE_UNUSED const int interp_src_Nxx2 =
    commondata->interp_src_Nxx2; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_Nxx2
MAYBE_UNUSED const int interp_src_Nxx_plus_2NGHOSTS0 =
    commondata->interp_src_Nxx_plus_2NGHOSTS0; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_Nxx_plus_2NGHOSTS0
MAYBE_UNUSED const int interp_src_Nxx_plus_2NGHOSTS1 =
    commondata->interp_src_Nxx_plus_2NGHOSTS1; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_Nxx_plus_2NGHOSTS1
MAYBE_UNUSED const int interp_src_Nxx_plus_2NGHOSTS2 =
    commondata->interp_src_Nxx_plus_2NGHOSTS2;     // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__interp_src_set_up::interp_src_Nxx_plus_2NGHOSTS2
MAYBE_UNUSED const BHA_REAL invdxx0 = params->invdxx0; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::invdxx0
MAYBE_UNUSED const BHA_REAL invdxx1 = params->invdxx1; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::invdxx1
MAYBE_UNUSED const BHA_REAL invdxx2 = params->invdxx2; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::invdxx2
MAYBE_UNUSED const BHA_REAL KO_diss_strength = commondata->KO_diss_strength; // nrpydev.infrastructures.BHaH.BHaHAHA.rhs_eval_KO_apply::KO_diss_strength
MAYBE_UNUSED const BHA_REAL max_radius_wrt_grid_center =
    commondata
        ->max_radius_wrt_grid_center; // nrpydev.infrastructures.BHaH.BHaHAHA.diagnostics_area_centroid_and_Theta_norms::max_radius_wrt_grid_center
MAYBE_UNUSED const BHA_REAL min_radius_wrt_grid_center =
    commondata
        ->min_radius_wrt_grid_center; // nrpydev.infrastructures.BHaH.BHaHAHA.diagnostics_area_centroid_and_Theta_norms::min_radius_wrt_grid_center
MAYBE_UNUSED const int nn = commondata->nn;             // nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all::nn
MAYBE_UNUSED const int nn_0 = commondata->nn_0;         // nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all::nn_0
MAYBE_UNUSED const int NUMGRIDS = commondata->NUMGRIDS; // nrpy.grid::NUMGRIDS
MAYBE_UNUSED const int Nxx0 = params->Nxx0;             // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::Nxx0
MAYBE_UNUSED const int Nxx1 = params->Nxx1;             // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::Nxx1
MAYBE_UNUSED const int Nxx2 = params->Nxx2;             // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::Nxx2
MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 =
    params->Nxx_plus_2NGHOSTS0; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::Nxx_plus_2NGHOSTS0
MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 =
    params->Nxx_plus_2NGHOSTS1; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::Nxx_plus_2NGHOSTS1
MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 =
    params->Nxx_plus_2NGHOSTS2; // nrpydev.infrastructures.BHaH.BHaHAHA.numgrid__evol_set_up::Nxx_plus_2NGHOSTS2
MAYBE_UNUSED const int output_diagnostics_every_nn =
    commondata->output_diagnostics_every_nn;           // nrpydev.infrastructures.BHaH.BHaHAHA.diagnostics::output_diagnostics_every_nn
MAYBE_UNUSED const BHA_REAL PI = params->PI;               // nrpy.reference_metric::PI
MAYBE_UNUSED const BHA_REAL RMAX = params->RMAX;           // nrpy.reference_metric::RMAX
MAYBE_UNUSED const BHA_REAL t_0 = commondata->t_0;         // nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all::t_0
MAYBE_UNUSED const BHA_REAL t_final = commondata->t_final; // nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all::t_final
MAYBE_UNUSED const BHA_REAL time = commondata->time;       // nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all::time
