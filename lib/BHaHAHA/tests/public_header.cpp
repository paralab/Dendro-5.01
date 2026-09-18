// The embedding application may define REAL independently of BHaHAHA.
#define REAL float
struct bhahaha_params_and_data_struct;
struct bhahaha_diagnostics_struct;
#include "BHaHAHA.h"
#include <cmath>
#include <type_traits>
static_assert(std::is_same<BHA_REAL, double>::value, "BHaHAHA must use its own real type");
static_assert(std::is_same<REAL, float>::value, "The application's REAL must remain unchanged");
int main() {
  bhahaha_params_and_data_struct params{};
  bah_poisoning_set_inputs(&params);
  BHA_REAL radii[48], r_min, dr;
  int nr;
  bah_radial_grid_cell_centered_set_up(48, 1.0, 0.0, 1.0, &nr, &r_min, &dr, radii);
  return !(std::isnan(params.cfl_factor) && nr > 0 && nr <= 48 && dr > 0 &&
           bah_error_message(BHAHAHA_SUCCESS) != nullptr);
}
