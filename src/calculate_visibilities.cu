#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "calculate_visibilities.h"
#include "shapelet_basis.h"
#include "cudacomplex.h"
#include "fundamental_coords.h"
#include "constants.h"
#include "source_components.h"


extern "C" void calculate_visibilities(float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                    catsource_t catsource, float *angles_array,
                    const int num_baselines, const int num_visis,
                    visibility_set_t *visibility_set,
                    float *sbf2) {

  /*START We should be able to do all this outside of this function and transfer in--------------*/
  //==========================================================================================
  float *d_X_diff = NULL;
  float *d_Y_diff = NULL;
  float *d_Z_diff = NULL;

  cudaMalloc( (void**)&d_X_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_X_diff, X_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_Y_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_Y_diff, Y_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_Z_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_Z_diff, Z_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );

  float *d_angles_array;
  cudaMalloc( (void**)&d_angles_array, 3*sizeof(float) );
  cudaMemcpy( d_angles_array, angles_array, 3*sizeof(float), cudaMemcpyHostToDevice );

  float *d_sha0s = NULL;
  float *d_cha0s = NULL;
  float *d_wavelengths = NULL;
  cudaMalloc( (void**)&d_sha0s, num_visis*sizeof(float) );
  cudaMemcpy( d_sha0s, visibility_set->sha0s, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_cha0s, num_visis*sizeof(float) );
  cudaMemcpy( d_cha0s, visibility_set->cha0s, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_wavelengths, num_visis*sizeof(float) );
  cudaMemcpy( d_wavelengths, visibility_set->wavelengths, num_visis*sizeof(float), cudaMemcpyHostToDevice );

  //
  /* END We should be able to do all this outside of this function and transfer in--------------*/

  float *d_sum_visi_real = NULL;
  float *d_sum_visi_imag = NULL;
  float *d_u_metres = NULL;
  float *d_v_metres = NULL;
  float *d_w_metres = NULL;
  float *d_us = NULL;
  float *d_vs = NULL;
  float *d_ws = NULL;

  cudaMalloc( (void**)&d_sum_visi_real, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_sum_visi_imag, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_u_metres, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_v_metres, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_w_metres, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_us, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_vs, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_ws, num_visis*sizeof(float) );

  dim3 grid, threads;

  threads.x = 128;
  threads.y = 1;
  threads.z = 1;
  grid.x = (int)ceil( (float)num_visis / (float)threads.x );
  grid.y = 1;
  grid.z = 1;

  kern_calc_uvw<<< grid, threads >>>( d_X_diff,
          d_Y_diff, d_Z_diff,
          d_u_metres, d_v_metres, d_w_metres,
          d_us, d_vs, d_ws, d_wavelengths,
          d_angles_array, d_cha0s, d_sha0s,
          num_visis, num_baselines);

  int num_points = catsource.n_points;
  int num_gauss = catsource.n_gauss;
  int num_shapes = catsource.n_shapes;

  if (num_points > 0) {

    float *d_point_ras=NULL;
    float *d_point_decs=NULL;
    float *d_point_fluxes=NULL;
    float *d_point_freqs=NULL;

    cudaMalloc( (void**)&(d_point_ras), num_points*sizeof(float) );
    cudaMemcpy( d_point_ras, catsource.point_ras, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_decs), num_points*sizeof(float) );
    cudaMemcpy( d_point_decs, catsource.point_decs, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_fluxes), num_points*sizeof(float) );
    cudaMemcpy( d_point_fluxes, catsource.point_fluxes, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_freqs), num_points*sizeof(float) );
    cudaMemcpy( d_point_freqs, catsource.point_freqs, num_points*sizeof(float), cudaMemcpyHostToDevice );

    float *d_ls=NULL;
    float *d_ms=NULL;
    float *d_ns=NULL;

    cudaMalloc( (void**)&d_ls, num_points*sizeof(float) );
    cudaMalloc( (void**)&d_ms, num_points*sizeof(float) );
    cudaMalloc( (void**)&d_ns, num_points*sizeof(float) );

    threads.x = 128;
    threads.y = 1;
    threads.z = 1;
    grid.x = (int)ceil( (float)num_points / (float)threads.x );
    grid.y = 1;
    grid.z = 1;

    kern_calc_lmn<<< grid, threads >>>(d_angles_array, d_point_ras, d_point_decs,
                               d_ls, d_ms, d_ns, num_points);

    if (num_points == 1) {
      threads.x = 128;
      threads.y = 1;
      grid.x = grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;
    }
    else {
      threads.x = 64;
      threads.y = 4;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = (int)ceil( ((float)num_points) / ((float)threads.y) );
    }

    kern_calc_visi_point<<<grid , threads>>>(d_point_ras,
            d_point_decs, d_point_fluxes, d_point_freqs,
            d_us, d_vs, d_ws,
            d_sum_visi_real, d_sum_visi_imag,
            d_angles_array, d_wavelengths,
            d_ls, d_ms, d_ns,
            num_points, num_visis);

    cudaFree( d_ns);
    cudaFree( d_ms);
    cudaFree( d_ls);
    cudaFree( d_point_freqs );
    cudaFree( d_point_fluxes );
    cudaFree( d_point_decs);
    cudaFree( d_point_ras);

  }//if point sources

  if (num_gauss > 0) {

    float *d_gauss_ras=NULL;
    float *d_gauss_decs=NULL;
    float *d_gauss_fluxes=NULL;
    float *d_gauss_freqs=NULL;
    float *d_gauss_pas=NULL;
    float *d_gauss_majors=NULL;
    float *d_gauss_minors=NULL;

    cudaMalloc( (void**)&(d_gauss_ras), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_ras, catsource.gauss_ras, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_decs), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_decs, catsource.gauss_decs, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_fluxes), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_fluxes, catsource.gauss_fluxes, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_pas), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_pas, catsource.gauss_pas, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_majors), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_majors, catsource.gauss_majors, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_minors), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_minors, catsource.gauss_minors, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_freqs), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_freqs, catsource.gauss_freqs, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    float *d_ls=NULL;
    float *d_ms=NULL;
    float *d_ns=NULL;

    cudaMalloc( (void**)&d_ls, num_gauss*sizeof(float) );
    cudaMalloc( (void**)&d_ms, num_gauss*sizeof(float) );
    cudaMalloc( (void**)&d_ns, num_gauss*sizeof(float) );

    threads.x = 128;
    threads.y = 1;
    threads.z = 1;
    grid.x = (int)ceil( (float)num_gauss / (float)threads.x );
    grid.y = 1;
    grid.z = 1;

    kern_calc_lmn<<< grid , threads >>>(d_angles_array, d_gauss_ras, d_gauss_decs,
                               d_ls, d_ms, d_ns, num_gauss);

    if (num_gauss == 1) {
     threads.x = 128;
     threads.y = 1;
     grid.x = grid.x = (int)ceil( (float)num_visis / (float)threads.x );
     grid.y = 1;
    }
    else {
     threads.x = 64;
     threads.y = 4;
     grid.x = (int)ceil( (float)num_visis / (float)threads.x );
     grid.y = (int)ceil( ((float)num_gauss) / ((float)threads.y) );
    }

    kern_calc_visi_gaussian<<< grid , threads>>>(d_gauss_ras,
            d_gauss_decs, d_gauss_fluxes, d_gauss_freqs,
            d_us, d_vs, d_ws,
            d_sum_visi_real, d_sum_visi_imag,
            d_angles_array, d_wavelengths,
            d_ls, d_ms, d_ns,
            d_gauss_pas, d_gauss_majors, d_gauss_minors,
            num_gauss, num_visis);

    cudaFree( d_ns);
    cudaFree( d_ms);
    cudaFree( d_ls);
    cudaFree( d_gauss_minors );
    cudaFree( d_gauss_majors);
    cudaFree( d_gauss_pas);
    cudaFree( d_gauss_freqs );
    cudaFree( d_gauss_fluxes );
    cudaFree( d_gauss_decs);
    cudaFree( d_gauss_ras);

  }//if gauss sources

  if (num_shapes > 0) {
    float *d_shape_ras=NULL;
    float *d_shape_decs=NULL;
    float *d_shape_fluxes=NULL;
    float *d_shape_freqs=NULL;
    float *d_shape_pas=NULL;
    float *d_shape_majors=NULL;
    float *d_shape_minors=NULL;

    float *d_shape_coeffs=NULL;
    float *d_shape_n1s=NULL;
    float *d_shape_n2s=NULL;
    float *d_shape_param_indexes=NULL;

    float *d_sbf=NULL;
    float *d_lsts=NULL;

    float *d_u_s_metres = NULL;
    float *d_v_s_metres = NULL;
    float *d_w_s_metres = NULL;

    float *d_shape_ls=NULL;
    float *d_shape_ms=NULL;
    float *d_shape_ns=NULL;

    //Who likes cudaMalloc cudaMallocs? We like cudaMalloc cudaMallocs
    cudaMalloc( (void**)&(d_shape_ras), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_ras, catsource.shape_ras, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_decs), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_decs, catsource.shape_decs, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_fluxes), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_fluxes, catsource.shape_fluxes, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_pas), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_pas, catsource.shape_pas, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_majors), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_majors, catsource.shape_majors, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_minors), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_minors, catsource.shape_minors, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_freqs), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_freqs, catsource.shape_freqs, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_coeffs), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_coeffs, catsource.shape_coeffs, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_n1s), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_n1s, catsource.shape_n1s, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_n2s), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_n2s, catsource.shape_n2s, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_param_indexes), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_param_indexes, catsource.shape_param_indexes, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_sbf), sbf_N*sbf_L*sizeof(float) );
    cudaMemcpy( d_sbf, sbf2, sbf_N*sbf_L*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_ls), num_shapes*sizeof(float) );
    cudaMalloc( (void**)&(d_shape_ms), num_shapes*sizeof(float) );
    cudaMalloc( (void**)&(d_shape_ns), num_shapes*sizeof(float) );

    cudaMalloc( (void**)&(d_lsts), num_visis*sizeof(float) );
    cudaMemcpy( d_lsts, visibility_set->lsts, num_visis*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&d_u_s_metres, num_shapes*num_visis*sizeof(float) );
    cudaMalloc( (void**)&d_v_s_metres, num_shapes*num_visis*sizeof(float) );
    cudaMalloc( (void**)&d_w_s_metres, num_shapes*num_visis*sizeof(float) );


    threads.x = 128;
    threads.y = 1;
    grid.x = (int)ceil( ((float)num_shapes / (float)threads.x) );
    grid.y = 1;

    kern_calc_lmn<<< grid , threads, 0 >>>(d_angles_array, d_shape_ras, d_shape_decs,
                             d_shape_ls, d_shape_ms, d_shape_ns, num_shapes);


    if (num_shapes == 1) {
      threads.x = 128;
      threads.y = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;
    }
    else {
      threads.x = 64;
      threads.y = 2;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = (int)ceil( ((float)num_shapes) / ((float)threads.y) );
    }

    kern_calc_uvw_shapelet<<< grid, threads >>>(d_X_diff, d_Y_diff, d_Z_diff,
          d_u_s_metres, d_v_s_metres, d_w_s_metres,
          d_lsts, d_shape_ras, d_shape_decs,
          num_baselines, num_visis, num_shapes);

    if (catsource.n_shape_coeffs == 1) {
      threads.x = 64;
      threads.y = 1;
      threads.z = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;
    }
    else {
      threads.x = 64;
      threads.y = 2;
      threads.z = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = (int)ceil( ((float)catsource.n_shape_coeffs) / ((float)threads.y) );
    }

    kern_calc_visi_shapelets<<< grid , threads, 0 >>>(d_shape_ras,
            d_shape_decs, d_shape_fluxes, d_shape_freqs,
            d_us, d_vs, d_ws, d_wavelengths,
            d_u_s_metres, d_v_s_metres, d_w_s_metres,
            d_sum_visi_real, d_sum_visi_imag,
            d_angles_array, d_shape_pas, d_shape_majors, d_shape_minors,
            d_shape_n1s, d_shape_n2s, d_shape_coeffs, d_shape_param_indexes,
            d_shape_ls, d_shape_ms, d_shape_ns,
            d_sbf,
            num_baselines, num_visis, catsource.n_shape_coeffs);

    cudaFree( d_shape_ns );
    cudaFree( d_shape_ms );
    cudaFree( d_shape_ls );
    cudaFree(d_w_s_metres);
    cudaFree(d_v_s_metres);
    cudaFree(d_u_s_metres);
    cudaFree(d_lsts);
    cudaFree(d_sbf);
    cudaFree(d_shape_param_indexes);
    cudaFree(d_shape_n2s);
    cudaFree(d_shape_n1s);
    cudaFree(d_shape_coeffs);
    cudaFree( d_shape_minors );
    cudaFree( d_shape_majors);
    cudaFree( d_shape_pas);
    cudaFree( d_shape_freqs );
    cudaFree( d_shape_fluxes );
    cudaFree( d_shape_decs);
    cudaFree( d_shape_ras);

  }//if shapelet

  //Get the results into host memory
  cudaMemcpy(visibility_set->sum_visi_real,d_sum_visi_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->sum_visi_imag,d_sum_visi_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  cudaMemcpy(visibility_set->us_metres,d_u_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->vs_metres,d_v_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->ws_metres,d_w_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  //Free up the GPU memory
  cudaFree( d_ws );
  cudaFree( d_vs );
  cudaFree( d_us );
  cudaFree( d_w_metres );
  cudaFree( d_v_metres );
  cudaFree( d_u_metres );
  cudaFree( d_sum_visi_imag );
  cudaFree( d_sum_visi_real );
  cudaFree( d_wavelengths );
  cudaFree( d_cha0s );
  cudaFree( d_sha0s );
  cudaFree( d_angles_array );
  cudaFree( d_Z_diff );
  cudaFree( d_Y_diff );
  cudaFree( d_X_diff );

}

extern "C" void copy_XYZ_to_GPU(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
                                float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                                const int num_baselines) {
  // printf("What the fuck chico man\n");

  cudaError_t err;
  err = cudaMalloc( (void**)&d_X_diff, num_baselines*sizeof(float) );
  if( cudaSuccess != err ) { \
    fprintf( stderr, "CUDA Error in Malloc\n");	\
  }

  if( cudaSuccess != err ) { \
    fprintf( stderr, "CUDA Error in Memcpy\n");	\
  }
  err = cudaMemcpy( d_X_diff, X_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );

  cudaMalloc( (void**)&d_Y_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_Y_diff, Y_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_Z_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_Z_diff, Z_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );

}