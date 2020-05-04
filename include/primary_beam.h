__device__ void twoD_Gaussian(float x, float y, float xo, float yo,
           float sigma_x, float sigma_y, float cos_theta, float sin_theta,
           float * d_beam_real, float * d_beam_imag);

__global__ void kern_gaussian_beam(float *d_beam_ls, float *d_beam_ms,
           float *ref_freq, float *d_freqs,
           float fwhm_lm, float cos_theta, float sin_theta, float sin_2theta,
           int num_freqs, int num_times, int num_components,
           float *d_beam_reals, float *d_beam_imags);

void calculate_gaussian_beam(int num_points, int num_time_steps, int num_freqs,
     float fwhm_lm, float cos_theta, float sin_theta, float sin_2theta,
     float *d_beam_ref_freq, float *d_freqs, float *d_beam_angles_array,
     float *beam_point_has, float *beam_point_decs,
     float *d_beam_reals, float *d_beam_imags);


extern "C" void testing_gaussian_beam( float *beam_has, float *beam_decs,
           float *beam_angles_array, float *beam_freqs, float *ref_freq_array,
           float *beam_ls, float *beam_ms,
           int num_components, int num_times, int num_freqs,
           float *beam_reals, float *beam_imags);
