.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

# Replace GNU Make's built-in cc, while honoring environment and command-line overrides.
ifeq ($(origin CC),default)
CC := gcc
endif

PROJECT_CPPFLAGS := -I.
ALL_CPPFLAGS = $(CPPFLAGS) $(PROJECT_CPPFLAGS)
CFLAGS = -std=gnu99 -O2 -march=native -g -Wall
CXXFLAGS = -O2 -g -Wall -Wno-unknown-pragmas -march=native
VALGRIND_CFLAGS = -std=gnu99 -O2 -g -Wall -Wno-unknown-pragmas
PIC_CFLAGS := 
LDFLAGS ?=
PROJECT_LDLIBS := -lm
ALL_LDLIBS = $(LDLIBS) $(PROJECT_LDLIBS)
DEPDIR := .deps
DEPFLAGS = -MMD -MP -MF $(DEPDIR)/$(@:.o=.d) -MT $@

SOURCES :=

# Add or remove sources only in the records below.
# Supported source suffixes: .c, .cc, .cpp, .cxx
# ADD_SOURCE(source, direct project headers...)
define ADD_SOURCE
$(eval $(basename $(strip $(1))).o: $(strip $(1)) $(strip $(2)))$(strip $(1))
endef

SOURCES += $(call ADD_SOURCE,bah_apply_bcs_inner_only.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_apply_bcs_r_maxmin_partial_r_hDD_upwinding.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_bcstruct_set_up.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_cfl_limited_timestep_based_on_h_equals_r.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_commondata_struct_set_to_default.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_diagnostics.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,bah_diagnostics_area_centroid_and_Theta_norms.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,bah_diagnostics_file_output.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_diagnostics_integration_weights.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_diagnostics_min_max_mean_radii_wrt_centroid.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,bah_diagnostics_proper_circumferences.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,bah_diagnostics_proper_circumferences_general.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,bah_error_message.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_find_horizon.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,bah_hDD_dD_and_W_dD_in_interp_src_grid_interior.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_initial_data.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,bah_interpolation_1d_radial_spokes_on_3d_src_grid.c,BHaH_defines.h interpolation_lagrange_uniform.h)
SOURCES += $(call ADD_SOURCE,bah_interpolation_2d_external_input_to_interp_src_grid.c,BHaH_defines.h interpolation_lagrange_uniform.h)
SOURCES += $(call ADD_SOURCE,bah_interpolation_2d_general__uniform_src_grid.c,interpolation_lagrange_uniform.h)
SOURCES += $(call ADD_SOURCE,bah_KO_apply.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_numgrid__evol_set_up.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,bah_numgrid__external_input_set_up.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,bah_numgrid__interp_src_set_up.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,bah_over_relaxation.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,bah_params_struct_set_to_default.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_poisoning_check_inputs.c)
SOURCES += $(call ADD_SOURCE,bah_poisoning_set_inputs.c)
SOURCES += $(call ADD_SOURCE,bah_quadratic_extrapolation.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_radial_grid_cell_centered_set_up.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_rfm_precompute_defines.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_rfm_precompute_free.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_rfm_precompute_malloc.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_rhs_eval.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_xx_to_Cart.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,bah_xyz_center_r_minmax.c,BHaH_defines.h BHaH_function_prototypes.h)
SOURCES += $(call ADD_SOURCE,MoL/bah_MoL_free_intermediate_stage_gfs.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,MoL/bah_MoL_malloc_intermediate_stage_gfs.c,BHaH_defines.h)
SOURCES += $(call ADD_SOURCE,MoL/bah_MoL_step_forward_in_time.c,BHaH_defines.h BHaH_function_prototypes.h)
OBJECTS := $(addsuffix .o,$(basename $(SOURCES)))
DEPFILES := $(addprefix $(DEPDIR)/,$(OBJECTS:.o=.d))

C_SOURCES := $(filter %.c,$(SOURCES))
CXX_SOURCES := $(filter %.cc %.cpp %.cxx,$(SOURCES))
ifeq ($(strip $(CXX_SOURCES)),)
LINKER := $(CC)
else
LINKER := $(CXX)
endif

OPENMP ?= 1

OPENMP_SUPPORTED := 0
ifeq ($(OPENMP),1)
hash := $(shell printf '\043')
OPENMP_TEST = printf '%s\n' '$(hash)include <omp.h>' '$(hash)ifndef _OPENMP' '$(hash)error OpenMP is unavailable' '$(hash)endif' 'int main(void) { return omp_get_max_threads() < 1; }'
OPENMP_C_SUPPORTED := 1
OPENMP_CXX_SUPPORTED := 1
ifneq ($(strip $(C_SOURCES)),)
OPENMP_C_SUPPORTED := $(shell $(OPENMP_TEST) | $(CC) $(ALL_CPPFLAGS) $(CFLAGS) $(PIC_CFLAGS) -fopenmp -x c - $(LDFLAGS) -o /dev/null >/dev/null 2>&1 && printf 1 || printf 0)
endif
ifneq ($(strip $(CXX_SOURCES)),)
OPENMP_CXX_SUPPORTED := $(shell $(OPENMP_TEST) | $(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) $(PIC_CFLAGS) -fopenmp -x c++ - $(LDFLAGS) -o /dev/null >/dev/null 2>&1 && printf 1 || printf 0)
endif
ifneq ($(strip $(C_SOURCES)$(CXX_SOURCES)),)
ifeq ($(OPENMP_C_SUPPORTED)$(OPENMP_CXX_SUPPORTED),11)
OPENMP_SUPPORTED := 1
endif
endif
endif

ifeq ($(OPENMP_SUPPORTED),1)
OMP_CFLAGS := -fopenmp
OMP_CXXFLAGS := -fopenmp
OMP_LDFLAGS := -fopenmp
else
OMP_CFLAGS := -Wno-unknown-pragmas
OMP_CXXFLAGS := -Wno-unknown-pragmas
OMP_LDFLAGS :=
endif

ALL_CFLAGS = $(CFLAGS) $(OMP_CFLAGS) $(PIC_CFLAGS)
ALL_CXXFLAGS = $(CXXFLAGS) $(OMP_CXXFLAGS) $(PIC_CFLAGS)
ALL_LDFLAGS = $(LDFLAGS) $(OMP_LDFLAGS)

.PHONY: all valgrind clean

all: libBHaHAHA.a

%.o: %.c
	@mkdir -p $(dir $(DEPDIR)/$(@:.o=.d))
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(DEPFLAGS) -c $< -o $@

%.o: %.cc
	@mkdir -p $(dir $(DEPDIR)/$(@:.o=.d))
	$(CXX) $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) $(DEPFLAGS) -c $< -o $@

%.o: %.cpp
	@mkdir -p $(dir $(DEPDIR)/$(@:.o=.d))
	$(CXX) $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) $(DEPFLAGS) -c $< -o $@

%.o: %.cxx
	@mkdir -p $(dir $(DEPDIR)/$(@:.o=.d))
	$(CXX) $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJECTS): Makefile
libBHaHAHA.a: $(OBJECTS)
	$(RM) $@
	$(AR) rcs $@ $(OBJECTS)

valgrind:
	+$(MAKE) clean
	+$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" OPENMP=0 all

# Remove nested build products and root runtime output.
clean:
	$(RM) -r libBHaHAHA.a $(DEPDIR) *.o */*.o */*/*.o *.d */*.d */*/*.d *.txt *.gp *.dat *.out *.log *.avi *.png *.bin

-include $(DEPFILES)
