#ifndef TELESCOPE_C_H
#define TELESCOPE_C_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// There is a fancy EveryBeam::telescope::Telescope class in C++ that we
// want to pass around in our C code. We can't do that directly, so we
// define a struct here that will hold a pointer to the C++ class.
typedef struct Telescope Telescope;

// Function to check the telescope type of a Measurement Set
char* check_ms_telescope_type(const char *ms_path);

// Function to create a new Telescope instance
Telescope* load_everybeam_telescope(int * status, const char *ms_path,
                                    const char *element_response_model,
                                    bool use_differential_beam,
                                    bool use_channel_frequency,
                                    const char *coeff_path,
                                    bool use_local_mwa);

// Function to delete a Telescope instance
void destroy_everybeam_telescope(Telescope* telescope);

void run_lofar_beam(Telescope *telescope,
                               int num_stations, int *station_idxs,
                               int num_dirs,
                               double ra0, double dec0,
                               double *ras, double *decs,
                               int num_times, double *mjd_sec_times,
                               int num_freqs, double *freqs,
                               bool apply_beam_norms, bool rotate,
                               bool element_only,
                               double _Complex * jones);

int load_and_run_lofar_beam(const char *ms_path,
                            const char *element_response_model,
                            const char *coeff_path,
                            int num_stations, int *station_idxs,
                            int num_dirs,
                            double ra0, double dec0,
                            double *ras, double *decs,
                            int num_times, double *mjd_sec_times,
                            int num_freqs, double *freqs,
                            bool apply_beam_norms, bool rotate,
                            bool element_only,
                            double _Complex * jones);



int load_and_run_mwa_beam(const char *ms_path,
                          const char *element_response_model,
                          const char *coeff_path,
                          int num_stations, int *station_idxs,
                          int num_dirs,
                          double ra0, double dec0,
                          double *azs, double *zas,
                          double *para_angles,
                          int num_times, double *mjd_sec_times,
                          int num_freqs, double *freqs,
                          bool apply_beam_norms, bool rotate,
                          bool element_only, bool ia_order,
                          double _Complex * jones);

#ifdef __cplusplus
}
#endif

#endif // TELESCOPE_C_H