
__device__ void extrap_flux(float *d_wavelengths,
                float *d_freqs, float *d_fluxes,
                int iComponent, int iBaseline,
                float * extrap_flux );


__global__ void calc_visi_point(float *d_point_ras, float *d_point_decs, float *d_point_fluxes, float *d_point_freqs,
      float *d_u_metres, float *d_v_metres, float *d_w_metres,
      float *d_u, float *d_v, float *d_w,
      float *d_sum_visi_real, float *d_sum_visi_imag,
      float *d_angles_array, float *d_wavelengths,
      float *d_ls, float *d_ms, float *d_ns,
      int num_points, int num_visis);
