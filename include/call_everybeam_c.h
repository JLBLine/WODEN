/*! \file
Functions to call the EveryBeam C++ code, from C. All functionality is
defined in `call_everybeam.cc`, and this header file is used to define
the C interface to that code. To avoid double documentation, I haven't actually
included the function definitions in `call_everybeam.h`, as `WODEN` only ever
calls these functions from `C`.

*/

#ifndef TELESCOPE_C_H
#define TELESCOPE_C_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @typedef Telescope
 * 
 * @brief A structure representing an EveryBeam telescope.
 * 
 * @details There is a fancy `EveryBeam::telescope::Telescope class` in C++ that we
 * want to pass around in our C code. We can't do that directly, so we
 * define a struct here that will hold a pointer to the C++ class. This is
 * classic compiled language trickery, which is probably bad practice, but I don't
 * know how else to do it, and it works.
 * 
 */
typedef struct Telescope Telescope;


/**
 * @brief Checks the telescope type of a measurement set (MS).
 *
 * @details This function takes the path to a measurement set (MS) file and determines
 * the type of telescope used to create the MS. As this is a wrapper around
 * EveryBeam code designed to work with `WODEN`, only recognises LOFAR, OSKAR, and MWA
 * as these are the only supported EVeryBeam telescopes in `WODEN`. Anything else will return
 * "UNKNOWN".
 *
 * @param ms_path The file path to the measurement set (MS).
 * @return A string representing the telescope type. The caller is responsible
 *         for freeing the returned string.
 */
char* check_ms_telescope_type(const char *ms_path);


/**
 * @brief Loads the EveryBeam Telescope object. This Telescope class can then
 * be used to calculate beam responses for a given set of directions.
 * 
 * @details This function loads the EveryBeam Telescope object from a Measurement Set (MS)
 * It runs `check_ms_telescope_type` to determine the type of telescope in the MS,
 * and errors if not LOFAR, OSKAR, or MWA. It also checks the `element_response_model`
 * (allowed values listed below). Note that two of the LOFAR models, `hamaker_lba` and `lobes`,
 * should in theory work, in practice I've never had any success with them. So I would
 * recommend sticking to `hamaker` for LOFAR.
 * 
 *
 * @param status Pointer to an integer that will hold the status of the operation.
 *               A value of 0 indicates success, while any other value indicates an error.
 * @param ms_path Path to the measurement set (MS) file.
 * @param element_response_model The type of element response. Allowed values are:
 * 'hamaker', 'hamaker_lba', 'lobes', 'skala40_wave', 'MWA'.
 * @param use_differential_beam Boolean flag indicating whether to use differential beam a.k.a normalisation out of EveryBeam
 * @param use_channel_frequency Boolean flag indicating whether to use channel frequency (not entirely sure what this does, always set to True in defualt `WODEN`).
 * @param coeff_path Path to the coefficients file, which should point to an HDF5 file in the case of MWA; not needed for LOFAR or OSKAR.
 * @param use_local_mwa Boolean flag indicating whether to use the local MWA (Murchison Widefield Array). This means loading an EveryBeam MWA telescope that uses `az,za` as input coords instead of `ra,dec`.
 *
 * @return A pointer to the initialized Telescope structure, or NULL if an error occurred.
 */
Telescope* load_everybeam_telescope(int * status, const char *ms_path,
                                    const char *element_response_model,
                                    bool use_differential_beam,
                                    bool use_channel_frequency,
                                    const char *coeff_path,
                                    bool use_local_mwa);


/**
 * @brief Destroys the given EveryBeam telescope object.
 *
 * @details Frees the memory associated with the given Telescope object; as 
 * the Telescope is actually a pointer to a C++ object, this function is
 * a wrapper around the C++ destructor.
 *
 * @param telescope A pointer to the Telescope object to be destroyed.
 */
void destroy_everybeam_telescope(Telescope* telescope);


/**
 * @brief Runs the phased array beamforming process for a given telescope.
 * 
 * @details This function is used for the LOFAR and OSKAR EveryBeam models. 
 * It calculates the beam response for a given set of directions, times, and frequencies.
 * Outputs are stored in the `jones` array, which should be pre-allocated to
 * `num_stations*num_times*num_freqs*num_dirs*4` (as there are 4 polarisations).
 * 
 * Outputs are ordered by station (slowest changing), time, frequency,
 * direction, and polarisation (fasting changing).
 *
 * @param telescope Pointer to the EveryBeam Telescope structure.
 * @param num_stations Number of stations.
 * @param station_idxs Array of station indices; must match the number of stations.
 * @param num_dirs Number of directions to process.
 * @param ra0 Right Ascension of the reference direction (radians).
 * @param dec0 Declination of the reference direction (radians).
 * @param ras Array of Right Ascensions for the directions (radians).
 * @param decs Array of Declinations for the directions (radians).
 * @param num_times Number of time samples.
 * @param mjd_sec_times Array of Modified Julian Dates (seconds).
 * @param num_freqs Number of frequency channels.
 * @param freqs Array of frequency channels (Hz).
 * @param apply_beam_norms Boolean flag to apply beam normalisation.
 * @param parallactic_rotate Boolean flag to apply parallactic angle rotation.
 * @param element_only Boolean flag to only calculate a single element (e.g. a single dipole instead of a beam-formed station).
 * @param iau_order Boolean flag to use IAU order for outputs; in IAU, X = North-South.
 * @param jones Output array of Jones matrices (complex values).
 */
void run_phased_array_beam(Telescope *telescope,
                    int num_stations, int *station_idxs,
                    int num_dirs,
                    double ra0, double dec0,
                    double *ras, double *decs,
                    int num_times, double *mjd_sec_times,
                    int num_freqs, double *freqs,
                    bool apply_beam_norms, bool parallactic_rotate,
                    bool element_only, bool iau_order,
                    double _Complex * jones);


/**
 * @brief Loads a LOFAR telescope from a Measurement Set (MS) file and runs
 * the phased array beamforming process for a given set of directions.
 * 
 * @details Runs `load_everybeam_telescope` with some default values, 
 * then runs `run_phased_array_beam` to calculate the beam response, and finally
 * destroys the telescope object. Currently identical to `load_and_run_oskar_beam`,
 * but kept separate in case we need to add LOFAR-specific functionality in the future.
 * 
 * Outputs are stored in the `jones` array, which should be pre-allocated to
 * `num_stations*num_times*num_freqs*num_dirs*4` (as there are 4 polarisations).
 * Outputs are ordered by station (slowest changing), time, frequency,
 * direction, and polarisation (fasting changing).
 *
 * @param ms_path Path to the measurement set (MS) file.
 * @param element_response_model The type of element response. Allowed values are:
 * 'hamaker', 'hamaker_lba', 'lobes'.
 * @param num_stations Number of stations.
 * @param station_idxs Array of station indices; must match the number of stations.
 * @param num_dirs Number of directions to process.
 * @param ra0 Right Ascension of the reference direction (radians).
 * @param dec0 Declination of the reference direction (radians).
 * @param ras Array of Right Ascensions for the directions (radians).
 * @param decs Array of Declinations for the directions (radians).
 * @param num_times Number of time samples.
 * @param mjd_sec_times Array of Modified Julian Dates (seconds).
 * @param num_freqs Number of frequency channels.
 * @param freqs Array of frequency channels (Hz).
 * @param apply_beam_norms Boolean flag to apply beam normalisation.
 * @param parallactic_rotate Boolean flag to apply parallactic angle rotation.
 * @param element_only Boolean flag to only calculate a single element (e.g. a single dipole instead of a beam-formed station).
 * @param iau_order Boolean flag to use IAU order for outputs; in IAU, X = North-South.
 * @param jones Output array of Jones matrices (complex values).
 */
int load_and_run_lofar_beam(const char *ms_path,
                            const char *element_response_model,
                            int num_stations, int *station_idxs,
                            int num_dirs,
                            double ra0, double dec0,
                            double *ras, double *decs,
                            int num_times, double *mjd_sec_times,
                            int num_freqs, double *freqs,
                            bool apply_beam_norms, bool parallactic_rotate,
                            bool element_only, bool iau_order,
                            double _Complex * jones);

                            /**
 * @brief Loads an OSKAR telescope from a Measurement Set (MS) file and runs
 * the phased array beamforming process for a given set of directions.
 * 
 * @details Runs `load_everybeam_telescope` with some default values, 
 * then runs `run_phased_array_beam` to calculate the beam response, and finally
 * destroys the telescope object. Currently identical to `load_and_run_lofar_beam`,
 * but kept separate in case we need to add OSKAR-specific functionality in the future.
 * 
 * Outputs are stored in the `jones` array, which should be pre-allocated to
 * `num_stations*num_times*num_freqs*num_dirs*4` (as there are 4 polarisations).
 * Outputs are ordered by station (slowest changing), time, frequency,
 * direction, and polarisation (fasting changing).
 *
 * @param ms_path Path to the measurement set (MS) file.
 * @param element_response_model The type of element response. Allowed values are:
 * 'skala40_wave'.
 * @param num_stations Number of stations.
 * @param station_idxs Array of station indices; must match the number of stations.
 * @param num_dirs Number of directions to process.
 * @param ra0 Right Ascension of the reference direction (radians).
 * @param dec0 Declination of the reference direction (radians).
 * @param ras Array of Right Ascensions for the directions (radians).
 * @param decs Array of Declinations for the directions (radians).
 * @param num_times Number of time samples.
 * @param mjd_sec_times Array of Modified Julian Dates (seconds).
 * @param num_freqs Number of frequency channels.
 * @param freqs Array of frequency channels (Hz).
 * @param apply_beam_norms Boolean flag to apply beam normalisation.
 * @param parallactic_rotate Boolean flag to apply parallactic angle rotation.
 * @param element_only Boolean flag to only calculate a single element (e.g. a single dipole instead of a beam-formed station).
 * @param iau_order Boolean flag to use IAU order for outputs; in IAU, X = North-South.
 * @param jones Output array of Jones matrices (complex values).
 */
int load_and_run_oskar_beam(const char *ms_path,
                            const char *element_response_model,
                            int num_stations, int *station_idxs,
                            int num_dirs,
                            double ra0, double dec0,
                            double *ras, double *decs,
                            int num_times, double *mjd_sec_times,
                            int num_freqs, double *freqs,
                            bool apply_beam_norms, bool parallactic_rotate,
                            bool element_only, bool iau_order,
                            double _Complex * jones);

/**
 * @brief Runs the MWA primary beam for a given MWALocal telescope.
 * 
 * @details Calculates the beam response for a given set of directions, times, and frequencies.
 * Outputs are stored in the `jones` array, which should be pre-allocated to
 * `num_stations*num_times*num_freqs*num_dirs*4` (as there are 4 polarisations).
 * Outputs are ordered by station (slowest changing), time, frequency,
 * direction, and polarisation (fasting changing).
 * 
 * Specifically, this function takes in azimuth and zenith angle (az,za) coordinates,
 * allowing us to control the precession of directions in the beam response.
 * 
 * EveryBeam MWA also doesn't seem to do a parallactic rotation, so do that rotation
 * in this wrapper function. Rather than link to erfa or try and use casacore,
 * just take the parallactic angle as an input. Very easy to calculate using
 * Python erfa via erfa.hd2pa, which is how `run_woden.py` does it. 
 * 
 * Note that although you can ask for multiple stations, this wrapper can
 * only calculate identical beams for each station (no way to pass bespoke
 * dipole ampltiudes for example). I've kept the machinery for future expansion,
 * but you should only ask for one station, otherwise you're just wasting
 * compute.
 *
 * @param telescope Pointer to the EveryBeam Telescope structure.
 * @param num_stations Number of stations.
 * @param station_idxs Array of station indices; must match the number of stations.
 * @param num_dirs Number of directions to process.
 * @param azs Azimuth for all directions (radians).
 * @param zas Zenith angles for all direcitons (radians).
 * @param para_angles Parallactic angle for all directions (radians).
 * @param num_times Number of time samples.
 * @param mjd_sec_times Array of Modified Julian Dates (seconds).
 * @param num_freqs Number of frequency channels.
 * @param freqs Array of frequency channels (Hz).
 * @param apply_beam_norms Boolean flag to apply beam normalisation.
 * @param parallactic_rotate Boolean flag to apply parallactic angle rotation.
 * @param element_only Boolean flag to only calculate a single element (e.g. a single dipole instead of a beam-formed station).
 * @param iau_order Boolean flag to use IAU order for outputs; in IAU, X = North-South.
 * @param jones Output array of Jones matrices (complex values).
 */
void run_mwa_beam(Telescope *telescope,
                  int num_stations, int *station_idxs,
                  int num_dirs,
                  double *azs, double *zas,
                  double *para_angles,
                  int num_times, double *mjd_sec_times,
                  int num_freqs, double *freqs,
                  bool apply_beam_norms, bool parallactic_rotate,
                  bool element_only, bool iau_order,
                  double _Complex * jones);


/**
 * @brief Loads and runs the MWA primary beam model.
 * 
 * Runs `load_everybeam_telescope` with some default values, 
 * then runs `run_mwa_beam` to calculate the beam response, and finally
 * destroys the telescope object.
 * 
 * Outputs are stored in the `jones` array, which should be pre-allocated to
 * `num_stations*num_times*num_freqs*num_dirs*4` (as there are 4 polarisations).
 * Outputs are ordered by station (slowest changing), time, frequency,
 * direction, and polarisation (fasting changing).
 * 
 * Note that although you can ask for multiple stations, this wrapper can
 * only calculate identical beams for each station (no way to pass bespoke
 * dipole ampltiudes for example). I've kept the machinery for future expansion,
 * but you should only ask for one station, otherwise you're just wasting
 * compute.
 *
 * @param ms_path Path to the measurement set.
 * @param element_response_model Path to the element response model.
 * @param coeff_path Path to the coefficient file.
 * @param num_stations Number of stations.
 * @param station_idxs Array of station indices.
 * @param num_dirs Number of directions.
 * @param azs Array of azimuth angles (in radians).
 * @param zas Array of zenith angles (in radians).
 * @param para_angles Array of parallactic angles (in radians).
 * @param num_times Number of time samples.
 * @param mjd_sec_times Array of times in Modified Julian Date seconds.
 * @param num_freqs Number of frequency samples.
 * @param freqs Array of frequencies (in Hz).
 * @param apply_beam_norms Flag to apply beam normalization.
 * @param parallactic_rotate Flag to apply parallactic rotation.
 * @param element_only Flag to use element-only model.
 * @param iau_order Flag to use IAU order.
 * @param jones Output array for Jones matrices.
 * 
 * @return 0 on success, non-zero on failure.
 */
int load_and_run_mwa_beam(const char *ms_path,
                          const char *element_response_model,
                          const char *coeff_path,
                          int num_stations, int *station_idxs,
                          int num_dirs,
                          double *azs, double *zas,
                          double *para_angles,
                          int num_times, double *mjd_sec_times,
                          int num_freqs, double *freqs,
                          bool apply_beam_norms, bool parallactic_rotate,
                          bool element_only, bool iau_order,
                          double _Complex * jones);

#ifdef __cplusplus
}
#endif

#endif // TELESCOPE_C_H