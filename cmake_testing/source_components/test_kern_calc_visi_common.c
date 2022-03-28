#include "test_kern_calc_visi_common.h"

void malloc_args_for_testing(args_for_testing_t *args_ft,
                             components_t *components,
                             int num_baselines,  int num_times,
                             int num_freqs, int num_components,
                             int num_coeffs,
                             e_component_type component_type) {

  int num_visis = num_baselines*num_times*num_freqs;
  int num_beam_values = num_freqs*num_times*num_components;

  args_ft->num_baselines = num_baselines;
  args_ft->num_times = num_times;
  args_ft->num_freqs = num_freqs;
  args_ft->num_visis = num_visis;
  args_ft->num_beam_values = num_beam_values;

  args_ft->us = malloc(num_visis*sizeof(user_precision_t));
  args_ft->vs = malloc(num_visis*sizeof(user_precision_t));
  args_ft->ws = malloc(num_visis*sizeof(user_precision_t));
  args_ft->allsteps_wavelengths = malloc(num_visis*sizeof(user_precision_t));

  args_ft->primay_beam_J00 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  args_ft->primay_beam_J01 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  args_ft->primay_beam_J10 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  args_ft->primay_beam_J11 = malloc(num_beam_values*sizeof(user_precision_complex_t));

  //Make sure arrays to hold summed visis are initialised to zero
  args_ft->sum_visi_XX_real = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_XX_imag = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_XY_real = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_XY_imag = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_YX_real = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_YX_imag = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_YY_real = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_YY_imag = calloc(num_visis, sizeof(user_precision_t));

  //Setup the component information

  // components->ras = malloc(num_components*sizeof(float));
  // components->decs = malloc(num_components*sizeof(float));
  components->ref_stokesI = malloc(num_components*sizeof(user_precision_t));
  components->ref_stokesQ = calloc(num_components, sizeof(user_precision_t));
  components->ref_stokesU = calloc(num_components, sizeof(user_precision_t));
  components->ref_stokesV = calloc(num_components, sizeof(user_precision_t));
  components->SIs = malloc(num_components*sizeof(user_precision_t));
  components->ref_freqs = malloc(num_components*sizeof(double));
  //GAUSS/SHAPELET STUFF
  components->pas = malloc(num_components*sizeof(user_precision_t));
  components->majors = malloc(num_components*sizeof(user_precision_t));
  components->minors = malloc(num_components*sizeof(user_precision_t));
  //
  components->ls = malloc(num_components*sizeof(double));
  components->ms = malloc(num_components*sizeof(double));
  components->ns = malloc(num_components*sizeof(double));

  //
  if (component_type == SHAPELET) {
    args_ft->u_shapes = malloc(num_components*num_visis*sizeof(user_precision_t));
    args_ft->v_shapes = malloc(num_components*num_visis*sizeof(user_precision_t));
    args_ft->w_shapes = malloc(num_components*num_visis*sizeof(user_precision_t));
    args_ft->num_coeffs = num_coeffs;
    args_ft->sbf = malloc( sbf_N * sbf_L * sizeof(user_precision_t) );
    args_ft->sbf = create_sbf(args_ft->sbf);


    // printf("%.5f\n",args_ft->sbf[5001] );

    components->n1s = malloc(num_coeffs*sizeof(user_precision_t));
    components->n2s = malloc(num_coeffs*sizeof(user_precision_t));
    components->shape_coeffs = malloc(num_coeffs*sizeof(user_precision_t));
    components->param_indexes = malloc(num_coeffs*sizeof(user_precision_t));
  }
}

void free_args_for_testing(args_for_testing_t *args_ft,
                           components_t components,
                           e_component_type component_type) {

  free( args_ft->us );
  free( args_ft->vs );
  free( args_ft->ws );
  free( args_ft->allsteps_wavelengths );

  free( args_ft->sum_visi_XX_real );
  free( args_ft->sum_visi_XX_imag );
  free( args_ft->sum_visi_XY_real );
  free( args_ft->sum_visi_XY_imag );
  free( args_ft->sum_visi_YX_real );
  free( args_ft->sum_visi_YX_imag );
  free( args_ft->sum_visi_YY_real );
  free( args_ft->sum_visi_YY_imag );

  free( args_ft->primay_beam_J00 );
  free( args_ft->primay_beam_J01 );
  free( args_ft->primay_beam_J10 );
  free( args_ft->primay_beam_J11 );

  // free( components.ras );
  // free( components.decs );
  free( components.ref_stokesI );
  free( components.ref_stokesQ );
  free( components.ref_stokesU );
  free( components.ref_stokesV );
  free( components.SIs );
  free( components.ref_freqs );
  free( components.pas );
  free( components.majors );
  free( components.minors );

  free( components.ls );
  free( components.ms );
  free( components.ns );

  if (component_type == SHAPELET) {
    free( args_ft->sbf );
    free( args_ft->u_shapes );
    free( args_ft->v_shapes );
    free( args_ft->w_shapes );
    free( components.n1s );
    free( components.n2s );
    free( components.shape_coeffs );
    free( components.param_indexes );
  }

  free( args_ft );

}

void create_lmn(components_t components) {
  int count = 0;
  for (int l_ind = 0; l_ind < 5; l_ind++) {
    for (int m_ind = 0; m_ind < 5; m_ind++) {
      components.ls[count] = -0.5 + 0.25*l_ind;
      components.ms[count] = -0.5 + 0.25*m_ind;
      components.ns[count] = sqrt(1 - components.ls[count]*components.ls[count] - components.ms[count]*components.ms[count]);

      count ++;
    }
  }
}

double _Complex visi_env_shape(int comp, int visi, int coeff,
                              args_for_testing_t *args_ft,
                              components_t components  ) {
  double pa = (double)components.pas[comp];
  double sinpa = sin(pa);
  double cospa = cos(pa);

  double u_shape = args_ft->u_shapes[comp*args_ft->num_visis + visi];
  double v_shape = args_ft->v_shapes[comp*args_ft->num_visis + visi];

  double x = (cospa*v_shape + sinpa*u_shape); // major axis
  double y = (-sinpa*v_shape + cospa*u_shape); // minor axis

  //Scales the FWHM to std to match basis functions, and account for the
  //basis functions being stored with beta = 1.0
  //Basis functions have been stored in such a way that x is in the same
  //direction as on sky, but y is opposite, so include negative here
  double const_x = (components.majors[comp]*SQRT_M_PI_2_2_LN_2)/sbf_dx;
  double const_y = -(components.minors[comp]*SQRT_M_PI_2_2_LN_2)/sbf_dx;

  // printf("%d %.5f %.5f\n", visi, const_x, const_y  );

  double _Complex Ipow_lookup[] = { 1.0 + I*0.0,
                           0.0 + I*1.0,
                          -1.0 + I*0.0,
                           0.0 + I*-1.0 };

  double xlow, xhigh, ylow, yhigh, u_value, v_value, f_hat;
  user_precision_t *sbf_n;

  // find the indices in the basis functions for u*beta_u and v*beta_v

  double xpos = x*const_x + sbf_c;
  double ypos = y*const_y + sbf_c;

  int xindex = (int)floor(xpos);
  int yindex = (int)floor(ypos);
  //
  int n1 = (int)components.n1s[coeff];
  int n2 = (int)components.n2s[coeff];

  // if ( n1 < 0 || n2 < 0 || n1 >= sbf_N || n2 >= sbf_N ) continue;

  f_hat = components.shape_coeffs[coeff];
  //
  sbf_n = &args_ft->sbf[n1*sbf_L];
  xlow  = (double)sbf_n[xindex];
  xhigh = (double)sbf_n[xindex+1];
  u_value = xlow + (xhigh-xlow)*(xpos-xindex);

  sbf_n = &args_ft->sbf[n2*sbf_L];
  ylow  = (double)sbf_n[yindex];
  yhigh = (double)sbf_n[yindex+1];
  v_value = ylow + (yhigh-ylow)*(ypos-yindex);

  // accumulate the intensity model for baseline pair (u,v)
  double _Complex visi_env = 0.0 + I*0.0;
  visi_env = visi_env + Ipow_lookup[(n1+n2) % 4] * f_hat * u_value*v_value;

  // printf("INSIDE THE THING %.5f %.5f\n",creal(visi_env), cimag(visi_env) );

  return visi_env;
}

double _Complex visi_env_gauss(int comp, int visi,
                              args_for_testing_t *args_ft,
                              components_t components  ) {
  double pa = components.pas[comp];
  double u = args_ft->us[visi];
  double v = args_ft->vs[visi];
  double sinpa = sin(pa);
  double cospa = cos(pa);

  double x =  cospa*v + sinpa*u; // major axis
  double y = -sinpa*v + cospa*u; // minor axis
  double invsig_x = components.majors[comp];
  double invsig_y = components.minors[comp];

  double _Complex visi_env = exp( -0.5 * ( x*x*invsig_x*invsig_x*M_PI_2_2_LN_2 + y*y*invsig_y*invsig_y*M_PI_2_2_LN_2 ) ) + I*0.0;

  return visi_env;
}

/*
Basic implementation of the measurement equation to get expected visibilities
Loops over components, gets expected flux and beam gain and sum. DO everything
in double precision to test against
*/
void get_expected(int visi, int num_components, int num_baselines,
                  int num_freqs, int beamtype,
                  args_for_testing_t *args_ft,
                  components_t components,
                  e_component_type component_type,
                  double * expec_re, double * expec_im) {
  double expec_re_inc, expec_im_inc;
  double temp, visi_freq;
  double flux_ratio, flux_extrap, xx_gain;
  int time_ind, freq_ind, beam_ind;
  * expec_re = 0.0;
  * expec_im = 0.0;
  for (int comp = 0; comp < num_components; comp++) {
    // printf("%.5f %.5f %.5f\n", components.ls[comp], components.ms[comp], components.ns[comp] );
    temp = 2*M_PI*( args_ft->us[visi]*components.ls[comp] + args_ft->vs[visi]*components.ms[comp] + args_ft->ws[visi]*(components.ns[comp]-1) );
    // sincos(temp, &(expec_im_inc), &(expec_re_inc));

    expec_im_inc = sin(temp);
    expec_re_inc = cos(temp);

    visi_freq = VELC / args_ft->allsteps_wavelengths[visi];
    flux_ratio = pow(visi_freq / components.ref_freqs[comp], components.SIs[comp]);
    flux_extrap = flux_ratio*components.ref_stokesI[comp];

    //Do some indexing so we can call up correct instrumental gain
    time_ind = (int)floorf( (float)visi / ((float)num_baselines * (float)num_freqs));
    freq_ind = (int)floorf( ((float)visi - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
    beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + comp;

    if (beamtype == NO_BEAM) {
      xx_gain = 1.0;
    } else {
      //In tests below have set gains purely real, and all pols equal, so
      //only need to call one pol and square it
      xx_gain = creal(args_ft->primay_beam_J00[beam_ind]);
      xx_gain = xx_gain*xx_gain;
    }

    expec_re_inc = expec_re_inc*flux_extrap*xx_gain;
    expec_im_inc = expec_im_inc*flux_extrap*xx_gain;

    double _Complex visi_env = 0.0 + I*0.0;

    if (component_type == POINT) {
      visi_env = 1.0 + I*0.0;
    }
    else if (component_type == GAUSSIAN) {
      visi_env = visi_env_gauss(comp, visi, args_ft, components);
      // visi_env = 1.0 + I*0.0;
    }
    else if (component_type == SHAPELET) {
      for (int coeff = 0; coeff < args_ft->num_coeffs; coeff++) {
        if (components.param_indexes[coeff] == comp) {
          // visi_env = visi_env_gauss(comp, visi, args_ft, components);
          visi_env += visi_env_shape(comp, visi, coeff, args_ft, components);

        }
      }
    }
    * expec_re += (expec_re_inc*creal(visi_env) - expec_im_inc*cimag(visi_env));
    * expec_im += (expec_re_inc*cimag(visi_env) + expec_im_inc*creal(visi_env));
  }
}

//Take input parameters and test whether GPU outputs match expectations
void test_visi_outputs(int num_visis, int num_components,
                       int num_baselines, int num_freqs,
                       e_beamtype beamtype,  args_for_testing_t *args_ft,
                       components_t components,
                       e_component_type component_type) {

  double frac_tol;
  if (component_type == POINT || component_type == GAUSSIAN) {
    #ifdef DOUBLE_PRECISION
      frac_tol = 1e-13;
    #else
      frac_tol = 7e-5;
    #endif
  }
  else {
    #ifdef DOUBLE_PRECISION
      frac_tol = 1e-12;
    #else
      frac_tol = 5e-3;
    #endif
  }


  double expec_re, expec_im;
  for (int visi = 0; visi < num_visis; visi++) {

      get_expected(visi, num_components, num_baselines, num_freqs,
                  beamtype, args_ft, components, component_type,
                  &expec_re, &expec_im);

    // // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
    // //       args_ft->sum_visi_XX_real[visi], args_ft->sum_visi_XX_imag[visi],
    // //       args_ft->sum_visi_XY_real[visi], args_ft->sum_visi_XY_imag[visi],
    // //       args_ft->sum_visi_YX_real[visi], args_ft->sum_visi_YX_imag[visi],
    // //       args_ft->sum_visi_YY_real[visi], args_ft->sum_visi_YY_imag[visi]);
    //
    if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY ) {
      //MWA beam has cross pols which double everything when using just Stokes I
      //and setting the cross pols to 1.0 as well as the gains
      //Also means cross-pols are non-zero

      // printf("%.8f %.8f\n",2*frac_tol*expec_im, 2.0*expec_im - args_ft->sum_visi_XX_imag[visi]);

      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_re, 2.0*expec_re, args_ft->sum_visi_XX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_im, 2.0*expec_im, args_ft->sum_visi_XX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_re, 2.0*expec_re, args_ft->sum_visi_XY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_im, 2.0*expec_im, args_ft->sum_visi_XY_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_re, 2.0*expec_re, args_ft->sum_visi_YX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_im, 2.0*expec_im, args_ft->sum_visi_YX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_re, 2.0*expec_re, args_ft->sum_visi_YY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_im, 2.0*expec_im, args_ft->sum_visi_YY_imag[visi]);
    }
    else {

      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re, expec_re, args_ft->sum_visi_XX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im, expec_im, args_ft->sum_visi_XX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol, 0.0, args_ft->sum_visi_XY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol, 0.0, args_ft->sum_visi_XY_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol, 0.0, args_ft->sum_visi_YX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol, 0.0, args_ft->sum_visi_YX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re, expec_re, args_ft->sum_visi_YY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im, expec_im, args_ft->sum_visi_YY_imag[visi]);

    }
  }
}

/*
Test the __global__ code that calculates visibilities for different type
of COMPONENTs
Vary the l,m,n coords but keep all other variables constant
*/
void test_kern_calc_visi_Varylmn(e_beamtype beamtype, e_component_type comptype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Component information
  components_t components;
  //Allocate memory
  malloc_args_for_testing(args_ft, &components, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, comptype);



  //Setup l,m args that span a decent chunk of sky
  create_lmn(components);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (size_t visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (size_t comp = 0; comp < num_components; comp++) {
    components.ref_stokesI[comp] = 1.0;
    components.SIs[comp] = 0.0;
    components.ref_freqs[comp] = 150e+6;

    //Set major,minor to 3 arcmins
    components.pas[comp] = 0.0;
    components.majors[comp] = 3.0*(DD2R / 60.0);
    components.minors[comp] = 3.0*(DD2R / 60.0);

  }

  //Make up some u,v,w values and scale by wavelength in correct order
  int count = 0;
  float freq_base = 150e+6;
  float freq_inc = 25e+6;
  float wavelength, frequency;
  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
      frequency = freq_base + freq_step*freq_inc;
      wavelength = VELC / frequency;
      for (int baseline = 0; baseline < num_baselines; baseline++) {
        args_ft->us[count] = ((baseline + 1)*10) / wavelength;
        args_ft->vs[count] = ((baseline + 1)*10) / wavelength;
        args_ft->ws[count] = ((baseline + 1)*10) / wavelength;

        args_ft->allsteps_wavelengths[count] = wavelength;

        count ++;
      }
    }
  }

  if (comptype == SHAPELET) {
    //Set the shapelet u,v,w same as the measurement equation one (this is not
    //true in reality but works fine for testing)

    for (size_t comp = 0; comp < num_components; comp++) {
      for (size_t visi = 0; visi < num_visis; visi++) {
        args_ft->u_shapes[comp*num_visis + visi] = args_ft->us[visi];
        args_ft->v_shapes[comp*num_visis + visi] = args_ft->vs[visi];
        args_ft->w_shapes[comp*num_visis + visi] = args_ft->ws[visi];
      }
    }

    //This means we have a single basis function for every component
    for (size_t coeff = 0; coeff < num_coeffs; coeff++) {
      components.n1s[coeff] = 0.0;
      components.n2s[coeff] = 0.0;
      components.shape_coeffs[coeff] = 1.0;
      components.param_indexes[coeff] = coeff;
    }
  }

  test_kern_calc_visi_all(num_components, num_baselines, num_coeffs,
          num_freqs, num_visis, num_times, beamtype, comptype,
          components, args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->u_shapes, args_ft->v_shapes, args_ft->w_shapes,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths, args_ft->sbf,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);

  //Check all results are within 0.1% of expected value
  // double frac_tol = 1e-3;
  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    beamtype, args_ft,components, comptype);

  free_args_for_testing( args_ft, components, comptype );
}

/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the component visibilities and beam gains constant and vary the fluxes
Test works for all primary beam types
*/
void test_kern_calc_visi_VarylmnVaryFlux(e_beamtype beamtype,
                                         e_component_type comptype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Component information
  components_t components;
  //Allocate memory
  malloc_args_for_testing(args_ft, &components, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, comptype);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(components);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (size_t visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (size_t comp = 0; comp < num_components; comp++) {
    components.ref_stokesI[comp] = comp;
    components.SIs[comp] = -0.8;
    components.ref_freqs[comp] = 150e+6;

    //Set major,minor to 3 arcmins
    components.pas[comp] = 0.0;
    components.majors[comp] = 3.0*(DD2R / 60.0);
    components.minors[comp] = 3.0*(DD2R / 60.0);
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  int count = 0;
  float freq_base = 150e+6;
  float freq_inc = 25e+6;
  float wavelength, frequency;


  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
      frequency = freq_base + freq_step*freq_inc;
      wavelength = VELC / frequency;

      for (int baseline = 0; baseline < num_baselines; baseline++) {
        args_ft->us[count] = ((baseline + 1)*10) / wavelength;
        args_ft->vs[count] = ((baseline + 1)*10) / wavelength;
        args_ft->ws[count] = ((baseline + 1)*10) / wavelength;

        args_ft->allsteps_wavelengths[count] = wavelength;

        count ++;
      }
    }
  }

  if (comptype == SHAPELET) {
  //Set the shapelet u,v,w same as the measurement equation one (this is not
  //true in reality but works fine for testing)

    for (size_t comp = 0; comp < num_components; comp++) {
      for (size_t visi = 0; visi < num_visis; visi++) {
        args_ft->u_shapes[comp*num_visis + visi] = args_ft->us[visi];
        args_ft->v_shapes[comp*num_visis + visi] = args_ft->vs[visi];
        args_ft->w_shapes[comp*num_visis + visi] = args_ft->ws[visi];
      }
    }

    //This means we have a single basis function for every component
    for (size_t coeff = 0; coeff < num_coeffs; coeff++) {
      components.n1s[coeff] = 0.0;
      components.n2s[coeff] = 0.0;
      components.shape_coeffs[coeff] = 1.0;
      components.param_indexes[coeff] = coeff;
    }
  }

  if (comptype == SHAPELET) {
    //Set the shapelet u,v,w same as the measurement equation one (this is not
    //true in reality but works fine for testing)

    for (size_t comp = 0; comp < num_components; comp++) {
      for (size_t visi = 0; visi < num_visis; visi++) {
        args_ft->u_shapes[comp*num_visis + visi] = args_ft->us[visi];
        args_ft->v_shapes[comp*num_visis + visi] = args_ft->vs[visi];
        args_ft->w_shapes[comp*num_visis + visi] = args_ft->ws[visi];
      }
    }

    //This means we have a single basis function for every component
    for (size_t coeff = 0; coeff < num_coeffs; coeff++) {
      components.n1s[coeff] = 0.0;
      components.n2s[coeff] = 0.0;
      components.shape_coeffs[coeff] = 1.0;
      components.param_indexes[coeff] = coeff;
    }
  }

  //Run the CUDA code
  test_kern_calc_visi_all(num_components, num_baselines, num_coeffs,
          num_freqs, num_visis, num_times, beamtype, comptype,
          components, args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->u_shapes, args_ft->v_shapes, args_ft->w_shapes,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths, args_ft->sbf,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);

  //Check all results are within 0.1% of expected value
  // double frac_tol = 5e-3;

  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    beamtype, args_ft, components, comptype);

  free_args_for_testing( args_ft, components, comptype );
}

/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the component visibilities and fluxes constant and vary the beam gains
Test works for all primary beam types
*/
void test_kern_calc_visi_VarylmnVaryBeam(e_beamtype beamtype,
                                         e_component_type comptype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Component information
  components_t components;
  //Allocate memory
  malloc_args_for_testing(args_ft, &components, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, comptype);


  //Setup l,m args that span a decent chunk of sky
  create_lmn(components);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (size_t beam = 0; beam < num_beam_values; beam++) {
    args_ft->primay_beam_J00[beam] = beam + I*0.0;
    args_ft->primay_beam_J01[beam] = beam + I*0.0;
    args_ft->primay_beam_J10[beam] = beam + I*0.0;
    args_ft->primay_beam_J11[beam] = beam + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (size_t comp = 0; comp < num_components; comp++) {
    components.ref_stokesI[comp] = 1.0;
    components.SIs[comp] = 0.0;
    components.ref_freqs[comp] = 150e+6;

    //Set major,minor to 3 arcmins
    components.pas[comp] = 0.0;
    components.majors[comp] = 3.0*(DD2R / 60.0);
    components.minors[comp] = 3.0*(DD2R / 60.0);
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  int count = 0;
  float freq_base = 150e+6;
  float freq_inc = 25e+6;
  float wavelength, frequency;

  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
      frequency = freq_base + freq_step*freq_inc;
      wavelength = VELC / frequency;

      for (int baseline = 0; baseline < num_baselines; baseline++) {
        args_ft->us[count] = ((baseline + 1)*10) / wavelength;
        args_ft->vs[count] = ((baseline + 1)*10) / wavelength;
        args_ft->ws[count] = ((baseline + 1)*10) / wavelength;

        args_ft->allsteps_wavelengths[count] = wavelength;

        count ++;
      }
    }
  }

  if (comptype == SHAPELET) {
    //Set the shapelet u,v,w same as the measurement equation one (this is not
    //true in reality but works fine for testing)

    for (size_t comp = 0; comp < num_components; comp++) {
      for (size_t visi = 0; visi < num_visis; visi++) {
        args_ft->u_shapes[comp*num_visis + visi] = args_ft->us[visi];
        args_ft->v_shapes[comp*num_visis + visi] = args_ft->vs[visi];
        args_ft->w_shapes[comp*num_visis + visi] = args_ft->ws[visi];
      }
    }

    //This means we have a single basis function for every component
    for (size_t coeff = 0; coeff < num_coeffs; coeff++) {
      components.n1s[coeff] = 0.0;
      components.n2s[coeff] = 0.0;
      components.shape_coeffs[coeff] = 1.0;
      components.param_indexes[coeff] = coeff;
    }
  }

  //Run the CUDA code
  test_kern_calc_visi_all(num_components, num_baselines, num_coeffs,
          num_freqs, num_visis, num_times, beamtype, comptype,
          components, args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->u_shapes, args_ft->v_shapes, args_ft->w_shapes,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths, args_ft->sbf,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);

  //Check all results are within 1% of expected value
  // double frac_tol = 1e-2;
  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    beamtype, args_ft,components, comptype);

  free_args_for_testing( args_ft, components, comptype );
}




void test_kern_calc_visi_VarylmnVaryPAMajMin(e_beamtype beamtype,
                                             e_component_type comptype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Component information
  components_t components;
  //Allocate memory
  malloc_args_for_testing(args_ft, &components, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, comptype);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(components);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (int visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (size_t comp = 0; comp < num_components; comp++) {
    components.ref_stokesI[comp] = 1.0;
    components.SIs[comp] = 0.0;
    components.ref_freqs[comp] = 150e+6;

    //Set major,minor to 3 arcmins
    components.pas[comp] = (comp + 1)*DD2R;
    components.majors[comp] = (comp + 1)*(DD2R / 60.0);
    components.minors[comp] = (comp + 2)*(DD2R / 60.0);
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  int count = 0;
  double freq_base = 150e+6;
  double freq_inc = 25e+6;
  user_precision_t wavelength, frequency;
  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
      frequency = freq_base + freq_step*freq_inc;
      wavelength = VELC / frequency;
      for (int baseline = 0; baseline < num_baselines; baseline++) {
        args_ft->us[count] = ((baseline + 1)*100) / wavelength;
        args_ft->vs[count] = ((baseline + 1)*100) / wavelength;
        //ws are usually smaller than u,v
        args_ft->ws[count] = ((baseline + 2)*10) / wavelength;

        args_ft->allsteps_wavelengths[count] = wavelength;

        count ++;
      }
    }
  }

  if (comptype == SHAPELET) {
    //Set the shapelet u,v,w same as the measurement equation one (this is not
    //true in reality but works fine for testing)

    for (size_t comp = 0; comp < num_components; comp++) {
      for (size_t visi = 0; visi < num_visis; visi++) {
        args_ft->u_shapes[comp*num_visis + visi] = args_ft->us[visi];
        args_ft->v_shapes[comp*num_visis + visi] = args_ft->vs[visi];
        args_ft->w_shapes[comp*num_visis + visi] = args_ft->ws[visi];
      }
    }

    //This means we have a single basis function for every component
    for (size_t coeff = 0; coeff < num_coeffs; coeff++) {
      components.n1s[coeff] = 0.0;
      components.n2s[coeff] = 0.0;
      components.shape_coeffs[coeff] = 1.0;
      components.param_indexes[coeff] = coeff;
    }
  }

  //Run the CUDA code
  test_kern_calc_visi_all(num_components, num_baselines, num_coeffs,
          num_freqs, num_visis, num_times, beamtype, comptype,
          components, args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->u_shapes, args_ft->v_shapes, args_ft->w_shapes,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths, args_ft->sbf,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);

  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    beamtype, args_ft,components, comptype);

  free_args_for_testing( args_ft, components, comptype );
}
