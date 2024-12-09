#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
// #include "gpucomplex.h"

// #include "calculate_visibilities_gpu.h"
// #include "fundamental_coords_gpu.h"
#include "constants.h"
// #include "calculate_visibilities_gpu.h"
#include "source_components_common.h"
// #include "primary_beam_gpu.h"
#include "hyperbeam_error.h"
#include "visibility_set.h"

//On link the external GPU functions if we're not compiling for CPU
//Having this here means
// #ifndef !defined(__NVCC__) && !defined(__HIPCC__)



//External GPU code we're linking in--------------------------------------------
extern void copy_CPU_beam_gains_to_GPU(components_t *components,
  beam_gains_t *d_beam_gains, int num_gains);

extern calc_visi_inouts_t * create_calc_visi_inouts_gpu(array_layout_t *array_layout,
        visibility_set_t *visibility_set, visibility_set_t *d_visibility_set,
        user_precision_t *sbf, woden_settings_t *woden_settings,
        int num_shapelets, int use_twobeams);

extern void fill_ant_to_baseline_mapping_gpu(int num_ants, int *d_ant1_to_baseline_map,
                                                 int *d_ant2_to_baseline_map);

extern void set_visibilities_to_zero_gpu(visibility_set_t *d_visibility_set,
                              visibility_set_t *chunk_visibility_set, int num_visis);

extern source_t * copy_chunked_source_to_GPU(source_t *chunked_source);

extern void calc_uvw_gpu(double *d_X_diff, double *d_Y_diff, double *d_Z_diff,
                             user_precision_t *d_u_metres,
                             user_precision_t *d_v_metres, user_precision_t *d_w_metres,
                             user_precision_t *d_us, user_precision_t *d_vs,
                             user_precision_t *d_ws, user_precision_t *d_allsteps_wavelengths,
                             double *d_allsteps_cha0s, double *d_allsteps_sha0s,
                             woden_settings_t *woden_settings);

extern void calc_uv_shapelet_gpu(user_precision_t *d_u_shapes, user_precision_t *d_v_shapes,
                                    int num_shapes,
                                    double *d_X_diff, double *d_Y_diff, double *d_Z_diff,
                                    double *d_ras, double *d_decs,
                                    woden_settings_t *woden_settings);
                                    

extern void set_auto_uvw_to_zero_gpu(int num_cross, int num_autos,
                              user_precision_t *d_u, user_precision_t *d_v,
                              user_precision_t *d_w);

extern void free_extrapolated_flux_arrays(components_t *d_components);

extern void free_d_components(source_t *d_chunked_source, e_component_type comptype);

extern void free_beam_gains_gpu(beam_gains_t *d_beam_gains, e_beamtype beamtype);

extern void calc_visi_point_or_gauss_gpu(components_t d_components,
                                        beam_gains_t d_component_beam_gains,
                                        calc_visi_inouts_t *d_calc_visi_inouts,
                                        visibility_set_t *d_visibility_set,
                                        int num_components, e_beamtype beamtype,
                                        e_component_type comptype,
                                        woden_settings_t *woden_settings);

extern void calc_visi_shapelets_gpu(components_t d_components,
                                        beam_gains_t d_component_beam_gains,
                                        calc_visi_inouts_t *d_calc_visi_inouts,
                                        visibility_set_t *d_visibility_set,
                                        int num_shapes, int num_shape_coeffs,
                                        e_beamtype beamtype,
                                        woden_settings_t *woden_settings);

extern void copy_gpu_visi_set_to_host(visibility_set_t *d_visibility_set,
                                      calc_visi_inouts_t *d_calc_visi_inouts,
                                      visibility_set_t *chunk_visibility_set,
                                      int num_visis);

extern void free_calc_visi_inouts_gpu(calc_visi_inouts_t *d_calc_visi_inouts,
                                          visibility_set_t *d_visibility_set,
                                          int num_shapelets, int use_twobeams);

//Function defined in calculate_visibilities_common.c---------------------------
void calculate_visibilities(array_layout_t *array_layout,
  source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
  woden_settings_t *woden_settings, visibility_set_t *visibility_set,
  user_precision_t *sbf);


// void calculate_component_visis(e_component_type comptype,
//                                calc_visi_inouts_t *d_calc_visi_inouts,
//                                woden_settings_t *woden_settings,
//                                beam_settings_t *beam_settings,
//                                source_t *source, source_t *d_chunked_source,
//                                int num_components, int num_shapelet_coeffs,
//                                visibility_set_t *d_visibility_set,
//                                int num_beams, int use_twobeams,
//                                int do_gpu);