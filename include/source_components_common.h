#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"

//Put all the GPU code we are linking in here

void malloc_extrapolated_flux_arrays_gpu(components_t *d_components, int num_comps,
                                        int num_freqs);

extern void extrap_power_laws_stokesI_gpu(components_t d_components,
                                   int n_powers, double *d_extrap_freqs,
                                   int num_extrap_freqs);

extern void extrap_curved_power_laws_stokesI_gpu(components_t d_components,
                                   int n_curves, double *d_extrap_freqs,
                                   int num_extrap_freqs);

extern void extrap_list_fluxes_stokesI_gpu(components_t d_components,
                                   int n_lists, double *d_extrap_freqs,
                                   int num_extrap_freqs);

extern void extrap_power_laws_stokesV_gpu(components_t d_components,
                                              double *d_extrap_freqs,
                                              int num_extrap_freqs);

extern void extrap_curved_power_laws_stokesV_gpu(components_t d_components,
                                                     double *d_extrap_freqs,
                                                     int num_extrap_freqs);

extern void polarisation_fraction_stokesV_gpu(components_t d_components,
                                                     double *d_extrap_freqs,
                                                     int num_extrap_freqs);

extern void extrap_list_fluxes_stokesV_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);

extern void extrap_power_laws_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);

extern void extrap_curved_power_laws_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);

extern void polarisation_fraction_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);

extern void extrap_list_fluxes_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);

extern void extrap_p_list_fluxes_linpol_gpu(components_t d_components,
                                                double *d_extrap_freqs,
                                                int num_extrap_freqs);

extern void apply_rotation_measure_gpu(components_t d_components,
                                                double *d_extrap_freqs,
                                                int num_extrap_freqs);