#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "common_testing_functions.h"
#include "source_components_common.h"
#include <mwa_hyperbeam.h>
#include "hyperbeam_error.h"
// #include "source_component_common_common.h"

//External GPU code we're linking in
extern void copy_outputs_source_component_common_gpu(int num_of_each_flux_type,
           source_t *d_chunked_source,
           d_beam_gains_t *d_beam_gains,
           woden_settings_t *woden_settings,
           beam_settings_t *beam_settings,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           user_precision_t *extrap_flux_I, user_precision_t *extrap_flux_Q,
           user_precision_t *extrap_flux_U, user_precision_t *extrap_flux_V,
           double *ls, double *ms, double *ns,
           e_component_type comptype);

extern double * malloc_freqs_gpu(int num_extrap_freqs, double *extrap_freqs);
extern void free_freqs_gpu(double *d_extrap_freqs);

extern void free_d_components(source_t *d_chunked_source,
                                  e_component_type comptype);

extern void free_extrapolated_flux_arrays(components_t *d_components);

extern void free_beam_gains_gpu(d_beam_gains_t *d_beam_gains, e_beamtype beamtype);

extern source_t * copy_chunked_source_to_GPU(source_t *chunked_source);

extern void malloc_extrapolated_flux_arrays_gpu(components_t *d_components, int num_comps,
                                        int num_freqs);

extern void malloc_beam_gains_gpu(d_beam_gains_t *d_component_beam_gains,
                                     int beamtype, int num_gains);

extern void calc_lmn_for_components_gpu(components_t *d_components,
                                        int num_components,
                                        woden_settings_t *woden_settings);


//Actual tests
void test_source_component_common_ConstantDecFEEBeam_multiant(e_component_type comptype,
                                                              int do_gpu);
void test_source_component_common_ConstantDecFEEBeamInterp_multiant(e_component_type comptype,
                                                              int do_gpu);