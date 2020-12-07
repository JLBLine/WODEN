#include "calculate_visibilities.h"

__device__ void extrap_flux(float *d_wavelengths, float *d_freqs,
           float *d_fluxes, int iComponent, int iBaseline,
           float * extrap_flux);

__global__ void kern_extrap_stokes(int num_visis, int num_components,
           float *d_wavelengths, float *d_ref_freqs, float *d_SIs,
           float *d_ref_stokesI, float *d_ref_stokesQ,
           float *d_ref_stokesU, float *d_ref_stokesV,
           float *d_flux_I, float *d_flux_Q,
           float *d_flux_U, float *d_flux_V );

extern "C" void test_extrap_flux(catsource_t *catsource,
           const int num_visis, const int num_components,
           float *wavelengths, float *flux_I, float *flux_Q,
           float *flux_U, float *flux_V );

__device__ cuFloatComplex calc_measurement_equation(float *d_us,
           float *d_vs, float *d_ws, float *d_ls, float *d_ms, float *d_ns,
           int iBaseline, int iComponent);

__device__ void apply_beam_gains(cuFloatComplex g1xx, cuFloatComplex g1xy,
          cuFloatComplex g1yx, cuFloatComplex g1yy,
          cuFloatComplex g2xx, cuFloatComplex g2xy,
          cuFloatComplex g2yx, cuFloatComplex g2yy,
          float flux_I, float flux_Q,
          float flux_U, float flux_V,
          cuFloatComplex visi,
          cuFloatComplex * visi_XX, cuFloatComplex * visi_XY,
          cuFloatComplex * visi_YX, cuFloatComplex * visi_YY,
          int beamtype );

__device__ void get_beam_gains(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
           cuFloatComplex * g1xx, cuFloatComplex * g1xy,
           cuFloatComplex * g1yx, cuFloatComplex * g1yy,
           cuFloatComplex * g2xx, cuFloatComplex * g2xy,
           cuFloatComplex * g2yx, cuFloatComplex * g2yy);

__device__ void update_sum_visis(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
           cuFloatComplex visi,
           float flux_I, float flux_Q, float flux_U, float flux_V,
           float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
           float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
           float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
           float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag);

void source_component_common(int num_components, int num_beam_values,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
           float *d_freqs, float *d_ls, float *d_ms, float *d_ns,
           float *d_ras, float *d_decs, float *azs, float *zas,
           float *sin_para_angs, float *cos_para_angs,
           float fwhm_lm, float cos_theta, float sin_theta, float sin_2theta,
           float *beam_has, float *beam_decs,
           woden_settings_t *woden_settings,
           beam_settings_t beam_settings,
           copy_primary_beam_t *FEE_beam);

// void setup_point_sources(float *d_point_ras, float *d_point_decs,
//     float *d_point_freqs, float *d_point_stokesI, float *d_point_stokesQ,
//     float *d_point_stokesU, float *d_point_stokesV, float *d_point_SIs,
//     float *d_ls, float *d_ms, float *d_ns, float *d_freqs,
//     catsource_t catsource, beam_settings_t beam_settings,
//     woden_settings_t *woden_settings, copy_primary_beam_t *FEE_beam,
//     float fwhm_lm, float cos_theta, float sin_theta, float sin_2theta,
//     cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
//     cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11);

__global__ void kern_calc_visi_point(float *d_point_ras, float *d_point_decs,
      float *d_point_freqs, float *d_point_stokesI, float *d_point_stokesQ,
      float *d_point_stokesU, float *d_point_stokesV, float *d_point_SIs,
      float *d_us, float *d_vs, float *d_ws,
      float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
      float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
      float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
      float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
      float *d_wavelengths,
      float *d_ls, float *d_ms, float *d_ns,
      int num_points, int num_baselines, int num_freqs, int num_visis,
      int num_times, int beamtype,
      cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
      cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11);

__global__ void kern_calc_visi_gaussian(float *d_gauss_ras, float *d_gauss_decs,
      float *d_gauss_freqs, float *d_gauss_stokesI, float *d_gauss_stokesQ,
      float *d_gauss_stokesU, float *d_gauss_stokesV, float *d_gauss_SIs,
      float *d_us, float *d_vs, float *d_ws,
      float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
      float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
      float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
      float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
      float *d_wavelengths,
      float *d_ls, float *d_ms, float *d_ns,
      float *d_gauss_pas, float *d_gauss_majors, float *d_gauss_minors,
      int num_gauss, int num_baselines, int num_freqs, int num_visis,
      int num_times, int beamtype,
      cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
      cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11);

__global__ void kern_calc_visi_shapelets(float *d_shape_ras,float *d_shape_decs,
      float *d_shape_freqs, float *d_shape_stokesI, float *d_shape_stokesQ,
      float *d_shape_stokesU, float *d_shape_stokesV, float *d_shape_SIs,
      float *d_us, float *d_vs, float *d_ws,
      float *d_wavelengths,
      float *d_u_s_metres, float *d_v_s_metres, float *d_w_s_metres,
      float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
      float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
      float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
      float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
      float *d_shape_pas, float *d_shape_majors, float *d_shape_minors,
      float *d_shape_n1s, float *d_shape_n2s, float *d_shape_coeffs,
      float *d_shape_param_indexes,
      float *d_shape_ls, float *d_shape_ms, float *d_shape_ns,
      float *d_sbf,
      int num_shapes, int num_baselines, int num_freqs, int num_visis,
      const int num_coeffs, int num_times, int beamtype,
      cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
      cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11);
