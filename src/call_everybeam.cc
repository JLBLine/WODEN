
// #include <EveryBeam/pointresponse/pointresponse.h

#include "call_everybeam.h"

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

#include <EveryBeam/coords/itrfdirection.h>
#include <EveryBeam/coords/itrfconverter.h>

using everybeam::coords::ItrfConverter;


extern "C" Telescope* load_everybeam_telescope(int * status,  const char *ms_path,
                                               const char *element_response_model,
                                               bool use_differential_beam,
                                               bool use_channel_frequency,
                                               const char *coeff_path,
                                               bool use_local_mwa) {

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
  // else if (element_response_tmp == "LOBES") {
  //   options.element_response_model = everybeam::ElementResponseModel::kLOBES;
  // }
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

    char message[128];

    sprintf(message, "ERROR: Requested EveryBeam element response model '%s' not recognised. ", element_response_model);
    log_message(message);

    return nullptr;

  }

  options.data_column_name = "DATA";
  options.beam_normalisation_mode = beam_normalisation_mode;
  options.use_channel_frequency = use_channel_frequency;
  options.coeff_path = coeff_path;
  options.use_local_mwa = use_local_mwa;


  std::unique_ptr<Telescope> telescope = Load(ms_path, options);

  return telescope.release();
}

extern "C" void destroy_everybeam_telescope(Telescope* telescope) {
  delete telescope;
}

extern "C" void run_lofar_beam(Telescope *telescope,
                               int num_stations, int *station_idxs,
                               int num_dirs,
                               double ra0, double dec0,
                               double *ras, double *decs,
                               int num_times, double *mjd_sec_times,
                               int num_freqs, double *freqs,
                               bool apply_beam_norms, bool rotate,
                               double _Complex * jones) {

  std::printf("This is the apply_beam_norms: %d\n", apply_beam_norms);
  
  PhasedArray& lofar_phased_array = static_cast<PhasedArray&>(*telescope);

  vector3r_t direction, phase_itrf;// = {0.0, 0.0, 0.0};
  double freq, time;

  std::vector<vector3r_t> direction_itrfs(num_dirs);

  // lock, since casacore::Direction is not thread-safe
  // The lock prevents different PhasedArrayPoints to calculate the
  // the station response simultaneously
  // std::unique_lock<std::mutex> lock(mutex_);

  aocommon::MC2x2 norm, normed, response;

  for (int ti = 0; ti < num_times; ti++) {
    time = mjd_sec_times[ti];
    std::unique_ptr<PointResponse> point_response =
      lofar_phased_array.GetPointResponse(time);
    PhasedArrayPoint& phased_array_point =
        static_cast<PhasedArrayPoint&>(*point_response);
    phased_array_point.SetParalacticRotation(rotate);

    everybeam::coords::ItrfConverter itrf_converter(time);

    phase_itrf = itrf_converter.RaDecToItrf(ra0, dec0);

    // std::printf("phase_itrf: %.12f %.12f %.12f\n", phase_itrf[0], phase_itrf[1], phase_itrf[2]);

    for (int ci = 0; ci < num_dirs; ci++) {
      direction_itrfs[ci] = itrf_converter.RaDecToItrf(ras[ci], decs[ci]);
    }

    for (int fi = 0; fi < num_freqs; fi++) {
      freq = freqs[fi];
      for (int si = 0; si < num_stations; si++) {
        int station_idx = station_idxs[si];

        if (apply_beam_norms) {
          norm = phased_array_point.UnnormalisedResponse(
          BeamMode::kFull, station_idx, freq, phase_itrf, phase_itrf, phase_itrf);

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

          response = phased_array_point.UnnormalisedResponse(
          BeamMode::kFull, station_idx, freq, direction, phase_itrf, phase_itrf);

          // std::printf("Response: %.5f %.5f, %.5f %.5f, %.5f %.5f, %.5f %.5f\n",
          //                                       real(response[0]), imag(response[0]),
          //                                       real(response[1]), imag(response[1]),
          //                                       real(response[2]), imag(response[2]),
          //                                       real(response[3]), imag(response[3]));

          int jones_index = 4*(si*num_times*num_freqs*num_dirs + ti*num_freqs*num_dirs + fi*num_dirs + ci);

          if (apply_beam_norms) {
            // std::printf("response[0] before: %.8f\n", response[0].real());
            aocommon::MC2x2::ATimesB(normed, norm, response);
            response = normed; 
            // std::printf("response[0] after: %.8f\n", response[0].real());
          }

          // std::printf("jones_index: %d\n", jones_index);
          jones[jones_index + 0] = {response[0].real(), response[0].imag()};
          jones[jones_index + 1] = {response[1].real(), response[1].imag()};
          jones[jones_index + 2] = {response[2].real(), response[2].imag()};
          jones[jones_index + 3] = {response[3].real(), response[3].imag()};
        }
      }
    }
  }
}


extern "C" int load_and_run_lofar_beam(const char *ms_path,
                                       const char *element_response_model,
                                       const char *coeff_path,
                                       int num_stations, int *station_idxs,
                                       int num_dirs,
                                       double ra0, double dec0,
                                       double *ras, double *decs,
                                       int num_times, double *mjd_sec_times,
                                       int num_freqs, double *freqs,
                                       bool apply_beam_norms, bool rotate,
                                       double _Complex * jones) {

  int status = 0;

  bool use_differential_beam = false;
  bool use_channel_frequency = true;
  bool use_local_mwa = true;
  
  Telescope *telescope = load_everybeam_telescope(&status, ms_path, element_response_model,
                                                  use_differential_beam, use_channel_frequency,
                                                  coeff_path, use_local_mwa);


  run_lofar_beam(telescope, num_stations, station_idxs,
                 num_dirs, ra0, dec0, ras, decs,
                 num_times, mjd_sec_times, num_freqs, freqs,
                 apply_beam_norms, rotate, jones);


  return status;
}