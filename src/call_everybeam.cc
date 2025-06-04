
/*! \file
Note that all documentation is in the header file `call_everybeam_c.h`, as this
is the C interface to the C++ code. We only ever bother calling this code from
C so made sense to chuck the doc-strings there.

*/

#include "call_everybeam.h"


#include <EveryBeam/telescope/lofar.h>
#include <EveryBeam/coords/itrfdirection.h>
#include <EveryBeam/coords/itrfconverter.h>

using casacore::MeasurementSet;
using everybeam::BeamNormalisationMode;
using everybeam::GetTelescopeType;
using everybeam::Load;
using everybeam::telescope::Telescope;
using everybeam::telescope::PhasedArray;
using everybeam::Station;
using everybeam::vector3r_t;
// using aocommon::CoordinateSystem;
using everybeam::BeamMode;
using everybeam::BeamNormalisationMode;
using everybeam::pointresponse::PointResponse;
using everybeam::pointresponse::PhasedArrayPoint;
using everybeam::real_t;
using everybeam::telescope::LOFAR;



#include <mutex>

using everybeam::coords::ItrfConverter;

std::mutex mtx;

extern "C" char* check_ms_telescope_type(const char *ms_path) {
  MeasurementSet ms(ms_path);

  const everybeam::TelescopeType telescope_type =
      everybeam::GetTelescopeType(ms);

  const char* telescope;

  switch (telescope_type) {
    case everybeam::TelescopeType::kLofarTelescope:
      telescope = "LOFAR";
      break;
    case everybeam::TelescopeType::kOSKARTelescope:
      telescope = "OSKAR";
      break;
    case everybeam::TelescopeType::kMWATelescope:
      telescope = "MWA";
      break;
    default:
      telescope = "UNKNOWN";
  }

  char *telescope_return = new char[strlen(telescope) + 1];
  strcpy(telescope_return, telescope);

  return telescope_return;
}


extern "C" Telescope* load_everybeam_telescope(int * status,  const char *ms_path,
                                               const char *element_response_model,
                                               bool use_differential_beam,
                                               bool use_channel_frequency,
                                               const char *coeff_path,
                                               bool use_local_mwa) {
  char message[256];
  std::string data_column = "DATA";

  BeamNormalisationMode beam_normalisation_mode =
            use_differential_beam ? BeamNormalisationMode::kPreApplied
                                  : BeamNormalisationMode::kNone;
  

  MeasurementSet ms(ms_path);

  const everybeam::TelescopeType telescope_type =
      everybeam::GetTelescopeType(ms);

  switch (telescope_type) {
    // case everybeam::TelescopeType::kAARTFAAC:
    case everybeam::TelescopeType::kLofarTelescope:
    case everybeam::TelescopeType::kOSKARTelescope:
    // case everybeam::TelescopeType::kSkaMidTelescope:
    case everybeam::TelescopeType::kMWATelescope:
      break;
    default:
      log_message("Currently WODEN only supports LOFAR, MWA and OSKAR EveryBeam models. The telescope in the linked MeasurementSet is one of these");
      * status = 1;
      return nullptr;
  }

  // Fill everybeam options
  everybeam::Options options;

  // log_message(element_response_model);

  std::string element_response_tmp = element_response_model;
  std::for_each(element_response_tmp.begin(), element_response_tmp.end(),
                [](char& c) { c = ::toupper(c); });

  // std::printf("Element response model: %s\n", element_response_tmp.c_str());

  if (element_response_tmp == "HAMAKER") {
    options.element_response_model = everybeam::ElementResponseModel::kHamaker;
  }
  else if (element_response_tmp == "HAMAKER_LBA") {
    options.element_response_model =
        everybeam::ElementResponseModel::kHamakerLba;
  }
  else if (element_response_tmp == "LOBES") {
    options.element_response_model = everybeam::ElementResponseModel::kLOBES;
  }
  // else if (element_response_tmp == "LWA") {
  //   options.element_response_model = everybeam::ElementResponseModel::kLwa;
  // }
  else if (element_response_tmp == "OSKAR_DIPOLE") {
    options.element_response_model =
        everybeam::ElementResponseModel::kOSKARDipole;
  }
  else if (element_response_tmp == "SKALA40_WAVE") {
    options.element_response_model =
        everybeam::ElementResponseModel::kOSKARSphericalWave;
  }
  // else if (element_response_tmp == "SKAMID_ANALYTICAL") {
  //   options.element_response_model =
  //       everybeam::ElementResponseModel::kSkaMidAnalytical;
  //  }
  else if (element_response_tmp == "MWA") {
    options.element_response_model = everybeam::ElementResponseModel::kDefault;
  }
  else {
    * status = 1;

    sprintf(message, "ERROR: Requested EveryBeam element response model '%s' not recognised. ", element_response_model);
    log_message(message);

    return nullptr;

  }

  // sprintf(message, "Attempting to read telescope from '%s'. ", ms_path);
  // log_message(message);

  options.data_column_name = "DATA";
  options.beam_normalisation_mode = beam_normalisation_mode;
  options.use_channel_frequency = use_channel_frequency;
  if (element_response_tmp == "MWA") {
    options.coeff_path = coeff_path;
  } else {
    options.coeff_path = "";
  }
  options.use_local_mwa = use_local_mwa;


  std::unique_ptr<Telescope> telescope = Load(ms_path, options);
  return telescope.release();
}

extern "C" void destroy_everybeam_telescope(Telescope* telescope) {
  delete telescope;
}

void print_arr(const char* name, double _Complex *jones, int offset, int num_times, int num_freqs,
               int num_dirs, int num_ants, int num_baselines, int *ant_to_baseline_map) {

  double _Complex arr[num_times*num_baselines];

  for (int ant = 0; ant < num_ants; ant++) {
    for (int time = 0; time < num_times; time++) {
      int ind = num_ants*time + ant;
      int jones_index = 4*(ant*num_times*num_freqs*num_dirs + time*num_freqs*num_dirs ); 

      // std::printf("THIS %d %d %d %d %d\n", baseline, time, si, ind, jones_index);

      arr[ind] = jones[jones_index + offset];
      // std::printf("jones_index: %d, offset: %d, ind: %d, %.5f + %.5f*I\n", jones_index, offset, ind, __real__ arr[ind], __imag__ arr[ind]);
      // std::printf("arr: %d %.5f + %.5f*I\n", ind, __real__ arr[ind], __imag__ arr[ind]);
    }
  }

  std::printf("double _Complex %s[NUM_ANTS*NUM_TIME_STEPS] = {", name);
  for (int gain = 0; gain < num_times*num_ants; gain++) {
    std::printf("%.10f + %.10f*I, ", __real__ arr[gain], __imag__ arr[gain]);
  }
  std::printf("};\n");
}

void print_gains_from_jones(int num_times, int num_freqs,
                            int num_dirs, int num_ants,
                            double _Complex * jones){

  int num_baselines = ((num_ants - 1)*num_ants) / 2;

  int *ant1_to_baseline_map = NULL;
  int *ant2_to_baseline_map = NULL;

  ant1_to_baseline_map = (int *)malloc(num_baselines*sizeof(int));
  ant2_to_baseline_map = (int *)malloc(num_baselines*sizeof(int));

  //These functions only do cross correlations, so create all combos of antennas
  //that make up all the crosses
  int cross_index = 0;
  for (int ant1 = 0; ant1 < num_ants-1; ant1++)
  {
    for (int ant2 = ant1 + 1; ant2 < num_ants; ant2++)
    {
      ant1_to_baseline_map[cross_index] = ant1;
      ant2_to_baseline_map[cross_index] = ant2;

      cross_index += 1;
    }
  }
  
  print_arr("gxs", jones, 0, num_times, num_freqs,
            num_dirs, num_ants, num_baselines, ant1_to_baseline_map);
  print_arr("Dxs", jones, 1, num_times, num_freqs,
            num_dirs, num_ants, num_baselines, ant1_to_baseline_map);
  print_arr("Dys", jones, 2, num_times, num_freqs,
            num_dirs, num_ants, num_baselines, ant1_to_baseline_map);
  print_arr("gys", jones, 3, num_times, num_freqs,
            num_dirs, num_ants, num_baselines, ant1_to_baseline_map);
  
}


extern "C" void run_phased_array_beam(Telescope *telescope,
                               int num_stations, int *station_idxs,
                               int num_dirs,
                               double ra0, double dec0,
                               double *ras, double *decs,
                               int num_times, double *mjd_sec_times,
                               int num_freqs, double *freqs,
                               bool apply_beam_norms, bool parallactic_rotate,
                               bool element_only, bool iau_order,
                               double _Complex * jones) {

  std::vector<vector3r_t> direction_itrfs(num_dirs);
  vector3r_t direction, phase_itrf; //, phase_itrf, station_itrf;

  aocommon::MC2x2 norm, normed, response; //, response;

  double freq, mjd_time;

  BeamMode beammode = element_only ? BeamMode::kElement : BeamMode::kFull;

  for (int ti = 0; ti < num_times; ti++) {
    mjd_time = mjd_sec_times[ti];

    std::unique_ptr<PointResponse> point_response = telescope->GetPointResponse(mjd_time);

    PhasedArrayPoint& phased_array_point =
          static_cast<PhasedArrayPoint&>(*point_response);
    phased_array_point.SetParalacticRotation(parallactic_rotate);
    // phased_array_point.SetParalacticRotation(false);

    everybeam::coords::ItrfConverter itrf_converter(mjd_time);

    // phase_itrf = itrf_converter.RaDecToItrf(ra0, dec0);
    // station_itrf = itrf_converter.RaDecToItrf(ra0, dec0);

    for (int ci = 0; ci < num_dirs; ci++) {
      direction_itrfs[ci] = itrf_converter.RaDecToItrf(ras[ci], decs[ci]);
      // std::printf("direction_itrfs ci: %d %.5f %.5f  %.8f %.8f %.8f\n", ci, ras[ci], decs[ci], direction_itrfs[ci][0], direction_itrfs[ci][1], direction_itrfs[ci][2]);
    }

    for (int fi = 0; fi < num_freqs; fi++) {
      freq = freqs[fi];
      for (int si = 0; si < num_stations; si++) {
        int station_idx = station_idxs[si];

        if (apply_beam_norms) {
          phase_itrf = itrf_converter.RaDecToItrf(ra0, dec0);
          norm = phased_array_point.Response(
                  beammode, station_idx, freq, phase_itrf, &mtx);

          // std::printf("norm: %.5f %.5f, %.5f %.5f, %.5f %.5f, %.5f %.5f\n",
          //                                         norm[0].real(), norm[0].imag(),
          //                                         norm[1].real(), norm[1].imag(),
          //                                         norm[2].real(), norm[2].imag(),
          //                                         norm[3].real(), norm[3].imag());

          aocommon::Matrix2x2::Invert(reinterpret_cast<std::complex<double>*>(&norm));

          // std::printf("norm inv: %.5f %.5f, %.5f %.5f, %.5f %.5f, %.5f %.5f\n",
          //                                         norm[0].real(), norm[0].imag(),
          //                                         norm[1].real(), norm[1].imag(),
          //                                         norm[2].real(), norm[2].imag(),
          //                                         norm[3].real(), norm[3].imag());
        }

        for (int ci = 0; ci < num_dirs; ci++) {
          direction = direction_itrfs[ci];
          //The mtx argument is a thread lock method. Lock is necessary to
          // prevent Response from being used simultaneously by different 
          //threads. This is because casacore::Direction is not thread-safe,
          //and it's laced all throughout EveryeBeam
          response = phased_array_point.Response(
                   beammode, station_idx, freq, direction, &mtx);

          int jones_index = 4*(si*num_times*num_freqs*num_dirs + ti*num_freqs*num_dirs + fi*num_dirs + ci);

          // std::printf("%d %d %d %d: %.5f %.5f, %.5f %.5f, %.5f %.5f, %.5f %.5f\n",
          //                               si, ti, fi, ci,
          //                               response[0].real(), response[0].imag(),
          //                               response[1].real(), response[1].imag(),
          //                               response[2].real(), response[2].imag(),
          //                               response[3].real(), response[3].imag());

          if (apply_beam_norms) {
            // std::printf("response[0] before: %.8f\n", response[0].real());
            aocommon::MC2x2::ATimesB(normed, norm, response);
            response = normed; 
            // std::printf("response[0] after: %.8f\n", response[0].real());
          }

          // std::printf("normed %d %d %d %d: %.5f %.5f, %.5f %.5f, %.5f %.5f, %.5f %.5f\n",
          //                               si, ti, fi, ci,
          //                               response[0].real(), response[0].imag(),
          //                               response[1].real(), response[1].imag(),
          //                               response[2].real(), response[2].imag(),
          //                               response[3].real(), response[3].imag());

          if (iau_order){
            jones[jones_index + 0] = {response[0].real(), response[0].imag()};
            jones[jones_index + 1] = {response[1].real(), response[1].imag()};
            jones[jones_index + 2] = {response[2].real(), response[2].imag()};
            jones[jones_index + 3] = {response[3].real(), response[3].imag()};
          } else {
            jones[jones_index + 0] = {response[3].real(), response[3].imag()};
            jones[jones_index + 1] = {response[2].real(), response[2].imag()};
            jones[jones_index + 2] = {response[1].real(), response[1].imag()};
            jones[jones_index + 3] = {response[0].real(), response[0].imag()};
          }
          // std::printf("s%d t%d f%d c%d: %.10f + %.10f*I, %.10f + %.10f*I, %.10f + %.10f*I, %.10f + %.10f*I\n",
          //   si, ti, fi, ci,
          //   __real__ jones[jones_index + 0], __imag__ jones[jones_index + 0],
          //   __real__ jones[jones_index + 1], __imag__ jones[jones_index + 1],
          //   __real__ jones[jones_index + 2], __imag__ jones[jones_index + 2],
          //   __real__ jones[jones_index + 3], __imag__ jones[jones_index + 3]);
        }
      }
    }
  }
  // print_gains_from_jones(num_times, num_freqs,
  //                        num_dirs, num_stations,
  //                        jones);

}


extern "C" int load_and_run_lofar_beam(const char *ms_path,
                                       const char *element_response_model,
                                       int num_stations, int *station_idxs,
                                       int num_dirs,
                                       double ra0, double dec0,
                                       double *ras, double *decs,
                                       int num_times, double *mjd_sec_times,
                                       int num_freqs, double *freqs,
                                       bool apply_beam_norms, bool parallactic_rotate,
                                       bool element_only,  bool iau_order,
                                       double _Complex * jones) {

  int status = 0;
  bool use_differential_beam = false;
  bool use_channel_frequency = true;
  bool use_local_mwa = true;
  const char *coeff_path = "";
  
  Telescope *telescope = load_everybeam_telescope(&status, ms_path, element_response_model,
                                                  use_differential_beam, use_channel_frequency,
                                                  coeff_path, use_local_mwa);

  run_phased_array_beam(telescope, num_stations, station_idxs,
                 num_dirs, ra0, dec0, ras, decs,
                 num_times, mjd_sec_times, num_freqs, freqs,
                 apply_beam_norms, parallactic_rotate, element_only, iau_order,
                 jones);

  destroy_everybeam_telescope(telescope);

  return status;
}

extern "C" int load_and_run_oskar_beam(const char *ms_path,
                                       const char *element_response_model,
                                       int num_stations, int *station_idxs,
                                       int num_dirs,
                                       double ra0, double dec0,
                                       double *ras, double *decs,
                                       int num_times, double *mjd_sec_times,
                                       int num_freqs, double *freqs,
                                       bool apply_beam_norms, bool parallactic_rotate,
                                       bool element_only,  bool iau_order,
                                       double _Complex * jones) {

  int status = 0;
  bool use_differential_beam = false;
  bool use_channel_frequency = true;
  bool use_local_mwa = true;
  const char *coeff_path = "";

  Telescope *telescope = load_everybeam_telescope(&status, ms_path, element_response_model,
                                                  use_differential_beam, use_channel_frequency,
                                                  coeff_path, use_local_mwa);

  run_phased_array_beam(telescope, num_stations, station_idxs,
                 num_dirs, ra0, dec0, ras, decs,
                 num_times, mjd_sec_times, num_freqs, freqs,
                 apply_beam_norms, parallactic_rotate, element_only, iau_order,
                 jones);

  destroy_everybeam_telescope(telescope);

  return status;
}



extern "C" void run_mwa_beam(Telescope *telescope,
                               int num_stations, int *station_idxs,
                               int num_dirs,
                               double *azs, double *zas,
                               double *para_angles,
                               int num_times, double *mjd_sec_times,
                               int num_freqs, double *freqs,
                               bool apply_beam_norms, bool parallactic_rotate,
                               bool element_only, bool iau_order,
                               double _Complex * jones) {

  double freq, mjd_time, az, za, para_angle;

  BeamMode beammode = element_only ? BeamMode::kElement : BeamMode::kFull;

  std::complex<float>* buffer = new std::complex<float>[4];
  std::complex<float>* rot_mat = new std::complex<float>[4];
  std::complex<float>* rotated = new std::complex<float>[4];

  for (int ti = 0; ti < num_times; ti++) {
    mjd_time = mjd_sec_times[ti];

    std::unique_ptr<PointResponse> point_response = telescope->GetPointResponse(mjd_time);

    // everybeam::coords::ItrfConverter itrf_converter(mjd_time);

    for (int fi = 0; fi < num_freqs; fi++) {
      freq = freqs[fi];
      // std::printf("DOING freq: %.3e\n", freq);
      for (int si = 0; si < num_stations; si++) {
        int station_idx = station_idxs[si];

        for (int ci = 0; ci < num_dirs; ci++) {

          az = azs[ ci*num_times + ti];
          za = zas[ ci*num_times + ti];
          para_angle = para_angles[ ci*num_times + ti];

          // std::printf("s%d t%d f%d c%d this %d: az, za, para %.8f, %.8f, %.8f\n",
          //           si, ti, fi, ci, ci*num_times + ti,
          //           az, za, para_angle);
          // std::printf("s%d t%d f%d c%d : time %.4f, freq %.1e, station %d\n",
          //           si, ti, fi, ci, mjd_time, freq, station_idx);

          //wtf is the fieled number
          int field = 0;
          point_response->Response(beammode, buffer,
                                  za, az, freq,
                                  station_idx, field);

          if (parallactic_rotate) {

            rot_mat[0] = sin(-para_angle);
            rot_mat[1] = -cos(-para_angle);
            rot_mat[2] = -cos(-para_angle);
            rot_mat[3] = -sin(-para_angle);

            aocommon::Matrix2x2::ATimesB(rotated, buffer, rot_mat);

            buffer[0] = rotated[0]; 
            buffer[1] = rotated[1];
            buffer[2] = rotated[2];
            buffer[3] = rotated[3];
          }

          int jones_index = 4*(si*num_times*num_freqs*num_dirs + ti*num_freqs*num_dirs + fi*num_dirs + ci);

          if (iau_order){
            jones[jones_index + 0] = {buffer[3].real(), buffer[3].imag()};
            jones[jones_index + 1] = {buffer[2].real(), buffer[2].imag()};
            jones[jones_index + 2] = {buffer[1].real(), buffer[1].imag()};
            jones[jones_index + 3] = {buffer[0].real(), buffer[0].imag()};
          } else {
            jones[jones_index + 0] = {buffer[0].real(), buffer[0].imag()};
            jones[jones_index + 1] = {buffer[1].real(), buffer[1].imag()};
            jones[jones_index + 2] = {buffer[2].real(), buffer[2].imag()};
            jones[jones_index + 3] = {buffer[3].real(), buffer[3].imag()};
          }

          // std::printf("s%d t%d f%d c%d: %.10f %.10f, %.10f %.10f, %.10f %.10f, %.10f %.10f\n",
          //           si, ti, fi, ci,
          //           __real__ jones[jones_index + 0], __imag__ jones[jones_index + 0],
          //           __real__ jones[jones_index + 1], __imag__ jones[jones_index + 1],
          //           __real__ jones[jones_index + 2], __imag__ jones[jones_index + 2],
          //           __real__ jones[jones_index + 3], __imag__ jones[jones_index + 3]);
        }
      }
    }
  }
}

extern "C" int load_and_run_mwa_beam(const char *ms_path,
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
                                     double _Complex * jones) {

  // std::printf("Loading MWA beam\n");

  int status = 0;
  bool use_differential_beam = false;
  bool use_channel_frequency = true;
  bool use_local_mwa = true;
  
  Telescope *telescope = load_everybeam_telescope(&status, ms_path, element_response_model,
                                                  use_differential_beam, use_channel_frequency,
                                                  coeff_path, use_local_mwa);

  run_mwa_beam(telescope, num_stations, station_idxs,
                 num_dirs, azs, zas, para_angles,
                 num_times, mjd_sec_times, num_freqs, freqs,
                 apply_beam_norms, parallactic_rotate, element_only, iau_order,
                 jones);

  destroy_everybeam_telescope(telescope);

  return status;
}
