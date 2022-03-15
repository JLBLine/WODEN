#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <assert.h>

#include "FEE_primary_beam.h"

/****************************
Find the nearest freq in the element pattern hdf file
****************************/
herr_t RTS_op_func (hid_t loc_id, const char *name, const H5L_info_t *info,
            void *operator_data) {
    herr_t          status = 0;
    H5O_info_t      infobuf;
    struct opdata   *od = (struct opdata *) operator_data;
                                /* Type conversion */
    char *ret;

    float freq_entry,freq_diff;

    status = H5Oget_info_by_name (loc_id, name, &infobuf, H5P_DEFAULT);

    ret = strstr(name,"X1_");
    if(ret != NULL){
      sscanf(name,"X1_%f",&freq_entry);
      freq_diff = fabs(od->freq_in - freq_entry);
      if(od->least_diff < 0 || freq_diff < od->least_diff){
	       od->least_diff = freq_diff;
	        od->freq_out = freq_entry;
      }
    }
    return status;
}

int RTS_MWAFEEInit(const char *h5filename, double freq_Hz,
                   RTS_MWA_FEE_beam_t *pb, user_precision_t *FEE_delays){

  hid_t       file, group;         /* handles */
  herr_t          status;
  H5O_info_t      infobuf;
  struct opdata   od;

  static float last_freq_in=0.0;
  static float last_freq_out=0.0;

  file = H5Fopen(h5filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(freq_Hz != last_freq_in) {
    /***************************************************************/
    /* Get nearest frequency by iterating over dataset in HDF file */

    //file = H5Fopen(H5FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
    status = H5Oget_info(file, &infobuf);
    od.freq_in = freq_Hz;
    od.freq_out = -1.0;
    od.least_diff = -1.0;

    group = H5Gopen(file, "/", H5P_DEFAULT);

    H5Literate(group, H5_INDEX_NAME, H5_ITER_INC, NULL, RTS_op_func, (void *) &od);

    H5Gclose(group);
      /***************************************************************/

    last_freq_out = od.freq_out;
    last_freq_in = freq_Hz;

  } else {
    od.freq_out = last_freq_out;
  }

  /* Set complex excitation voltages */

  user_precision_t lam;
  user_precision_complex_t Vcplx[N_COPOL][NUM_DIPOLES];

  lam = VELC/(od.freq_out);

  for (int pol=0; pol < N_COPOL; pol++){
    for (int i=0; i<NUM_DIPOLES; i++) {
      user_precision_t phase;
      phase = (-2.0*M_PI*DQ/lam)*FEE_delays[i];

      //TODO implement dipole flagging and/or dipole amplitudes here
      //eg Vcplx[pol][i] = amps[pol][i] * cexp(1.0j * phase);
      Vcplx[pol][i] = cexp(1.0j * phase);
    } //end for dipole
  } // end for pol

  /* Get max length of spherical wave table */

  hid_t  dataset, dataspace;
  hsize_t dims_out[2];

  int n_ant = 16;
  int max_length=0;
  char table_name[80];
  int n_pols = 2;

  /* Determine the maximum required size for the beam mode data */
  /* Check both the X and Y pols*/

  for(int i = 1; i<n_ant+1;i++){
    sprintf(table_name,"X%d_%d",i,(int)od.freq_out);
    dataset = H5Dopen(file, table_name, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    status  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if(dims_out[1]/2 > max_length){
       max_length = dims_out[1]/2; // Complex number
    }
    H5Dclose(dataset);
    H5Sclose(dataspace);

    sprintf(table_name,"Y%d_%d",i,(int)od.freq_out);
    dataset = H5Dopen(file, table_name, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    status  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if(dims_out[1]/2 > max_length){
       max_length = dims_out[1]/2; // Complex number
    }
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }

  double _Complex *Q1, *Q2;

  /* Allocate space for two polarizations */

  pb->Q1 = (double _Complex **)malloc(n_pols*sizeof(double _Complex *));
  pb->Q2 = (double _Complex **)malloc(n_pols*sizeof(double _Complex *));
  pb->M = (double **)malloc(n_pols*sizeof(double *));
  pb->N = (double **)malloc(n_pols*sizeof(double *));

  for(int i=0;i<n_pols;i++){
    pb->Q1[i]= (double _Complex *)malloc(sizeof(double _Complex)*max_length);
    pb->Q2[i]= (double _Complex *)malloc(sizeof(double _Complex)*max_length);
    pb->M[i] = NULL;
    pb->N[i] = NULL;
    for(int ii=0;ii<max_length;ii++){
      pb->Q1[i][ii] = 0.0 + 0.0j;
      pb->Q2[i][ii] = 0.0 + 0.0j;
    }
  }

  Q1 = (double _Complex *)malloc(sizeof(double _Complex) * max_length);
  Q2 = (double _Complex *)malloc(sizeof(double _Complex) * max_length);

  // Read in Q_modes_all
  // double Q_modes_all[3][2046];

  dataset = H5Dopen(file,"modes", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);
  status  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

  int q_modes_all_dim1 = (int)dims_out[1];

  // printf("About to malloc %dx%d\n",dims_out[0], dims_out[1] );
  double *Q_modes_all = malloc(dims_out[0]*dims_out[1]*sizeof(double));

  status = H5Dread(dataset, H5T_NATIVE_DOUBLE,
	     H5S_ALL, H5S_ALL, H5P_DEFAULT, Q_modes_all);

  H5Dclose(dataset);

  // Read in modes for each antenna
  double *Q_all, *Q_modes;

  int Nmax[n_pols];

  for(int pol=0;pol<n_pols;pol++){
    Nmax[pol]=0;
    for(int ant_i=1; ant_i<n_ant+1;ant_i++){

      // re-initialise accumulation variables
      for(int i=0;i<max_length;i++){
        Q1[i] = 0.0 + 0.0j;
        Q2[i] = 0.0 + 0.0j;
      }

      if(pol==0){
        sprintf(table_name, "X%d_%d", ant_i, (int)od.freq_out);
      } else {
        sprintf(table_name, "Y%d_%d", ant_i, (int)od.freq_out);
      }

      dataset = H5Dopen(file, table_name, H5P_DEFAULT);
      dataspace = H5Dget_space(dataset);
      status  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
      // Q_all = malloc(sizeof(double) * dims_out[0] * dims_out[1]);
      Q_all = malloc(sizeof(double) * dims_out[0] * dims_out[1]);
      status = H5Dread(dataset, H5T_NATIVE_DOUBLE,
      	     H5S_ALL, H5S_ALL, H5P_DEFAULT, Q_all);
      H5Dclose(dataset);

      int my_len;
      my_len = (int)dims_out[1]; // max(Q_all.shape)

      // Get Q_modes for this antenna

      // Q_modes = Q_modes_all[0:my_len,:]
      // 3 is first dimension of Q_modes_all (why??)

      Q_modes = malloc(sizeof(double) * my_len * 3);

      for(int q_i =0; q_i < my_len; q_i++){
        for(int q_j = 0; q_j < 3; q_j++){
          //	  Q_modes[q_i + my_len*q_j] = Q_modes_all[q_i][q_j];
          // Q_modes[q_i + my_len*q_j] = (double)Q_modes_all[q_j][q_i];
          Q_modes[q_i + my_len*q_j] = Q_modes_all[q_i + q_modes_all_dim1*q_j];
        }
      }

      // convert Qall to M, N, Q1, Q2 vectors for processing

      // find s=1 and s=2 indices
      // only find s1 and s2 for this antenna

      int *s1, *s2;
      double *M_vec, *N_vec;
      int nM, nN;

      s1 = malloc(sizeof(int) * my_len);
      s2 = malloc(sizeof(int) * my_len);
      nM = nN = 0;

      for(int s = 0; s<my_len;s++){
        if(Q_modes[s] <=1){
          s1[s] = 1;
          nM++;
        } else {
          s1[s] = 0;
        }
        if(Q_modes[s] > 1){
          s2[s] = 1;
          nN++;
        } else {
          s2[s] = 0;
        }
      }

      M_vec = malloc(sizeof(double) * nM);
      N_vec = malloc(sizeof(double) * nN);

      int nC, mC, max_n;

      nC = mC = max_n = 0;

      for(int s = 0; s < my_len; s++){
        if(s1[s]){
          M_vec[mC] = Q_modes[s + my_len];
          mC++;
        }
        if(s2[s]){
          N_vec[nC] = Q_modes[s + 2 * my_len];
          if(N_vec[nC] > max_n) max_n = N_vec[nC];
          nC++;
        }
      }

      if(max_n > Nmax[pol]){
        if(pb->M[pol] != NULL){
          free(pb->M[pol]);
          free(pb->N[pol]);
        }
        pb->M[pol] = (double *)malloc(sizeof(double) * nM);
        pb->N[pol] = (double *)malloc(sizeof(double) * nN);

        memcpy(pb->M[pol],M_vec,sizeof(double) * nM);
        memcpy(pb->N[pol],N_vec,sizeof(double) * nN);

        Nmax[pol] = max_n;
        pb->nmax = (int)Nmax[pol];
        assert(nM==nN);
        pb->nMN = nM;
        if (pol==0) {
          pb->nMN0 = nM;
        }
        else if (pol==1) {
          pb->nMN1 = nM;
        }
      }

      nC = mC = 0;

      for(int s = 0; s < my_len; s++){
        if(s1[s]){
          Q1[mC] = Q_all[s] * cexp(1.0j * Q_all[s + my_len] * M_PI/180.0);
          mC++;
        }
        if(s2[s]){
          Q2[nC] = Q_all[s] * cexp(1.0j * Q_all[s + my_len] * M_PI/180.0);
          nC++;
        }
      }

      for(int i = 0; i < max_length; i++){
        pb->Q1[pol][i] += Q1[i] * Vcplx[pol][ant_i-1];
        pb->Q2[pol][i] += Q2[i] * Vcplx[pol][ant_i-1];
      }

      free(M_vec);
      free(N_vec);
      free(s1);
      free(s2);
      free(Q_modes);
      free(Q_all);
    // printf("END OF DIPOLE LOOP %d %d %d %d\n",nM, nN, max_n, pb->nMN );
    } // loop over ant_i
  } // loop over XY pols

  //Unfortunately, for 119 MHz, the Y pol has a higher order of spherical
  //harmonic that the X pol. This means the arrays in pb->M[1], pb->N[1]
  //are larger than in pb->M[0], pb->N[0], which hurts the GPU later down the
  //road. So do some realloc and chuck some extra zeros on the end
  if( Nmax[1] > Nmax[0]){
    printf("\nNmax[1] > Nmax[0], fiddling the length of pb->M[0], pb->N[0] \n");

    //Add extra memory
    pb->M[0] = realloc(pb->M[0],sizeof(double)*pb->nMN1);
    pb->N[0] = realloc(pb->N[0],sizeof(double)*pb->nMN1);
    //Set new contents to zero. The coeffs in Q1, Q2 should be zero also,
    //meaning these should contribute zero power to the overall beam.
    //Just stops the GPU kernel wigging out from different size arrays
    for (int i = pb->nMN0; i < pb->nMN1; i++) {
      pb->M[0][i] = 0;
      pb->N[0][i] = 0;
      pb->Q1[0][i] = 0.0 + 0.0j;
      pb->Q2[0][i] = 0.0 + 0.0j;
    }

  }

  free(Q1);
  free(Q2);
  free(Q_modes_all);

  H5Fclose(file);
  H5Sclose(dataspace);

  pb->m_range = malloc((2*pb->nmax + 1)*sizeof(user_precision_t) );

  for (int m = 0; m < 2*pb->nmax + 1; m++) {
    pb->m_range[m] = -(user_precision_t)pb->nmax + (user_precision_t)m;
  }

  return status;

}

void RTS_freeHDFBeam(RTS_MWA_FEE_beam_t *pb){

  if(pb==NULL) return;

  int n_pols = 2;

  for(int n=0; n < n_pols; n++){
    free(pb->Q1[n]);
    free(pb->Q2[n]);
    free(pb->M[n]);
    free(pb->N[n]);
  }

  free(pb->Q1);
  free(pb->Q2);
  free(pb->M);
  free(pb->N);

}


int multifreq_RTS_MWAFEEInit(beam_settings_t *beam_settings,
                             woden_settings_t *woden_settings,
                             double *beam_freqs) {

  //Find out which frequency in the stored MWA FEE coeffs is closest to
  //our lowest frequency
  herr_t status;
  hid_t file, group;         /* handles */
  H5O_info_t infobuf;
  struct opdata od;

  file = H5Fopen(woden_settings->hdf5_beam_path, H5F_ACC_RDONLY, H5P_DEFAULT);
  status = H5Oget_info(file, &infobuf);

  beam_settings->num_MWAFEE = (int)woden_settings->num_freqs;
  beam_settings->MWAFEE_freqs = malloc(beam_settings->num_MWAFEE*sizeof(double));

  //Setup arrays to hold the beam and zenith beams for all frequencies
  beam_settings->FEE_beams = malloc(beam_settings->num_MWAFEE*sizeof(RTS_MWA_FEE_beam_t));
  beam_settings->FEE_beam_zeniths = malloc(beam_settings->num_MWAFEE*sizeof(RTS_MWA_FEE_beam_t));

  printf("Initialising MWA FEE beams...");

  for (int freq_ind = 0; freq_ind < beam_settings->num_MWAFEE; freq_ind++) {
    od.freq_in = (float)beam_freqs[freq_ind];
    od.freq_out = -1.0;
    od.least_diff = -1.0;

    group = H5Gopen(file, "/", H5P_DEFAULT);
    H5Literate(group, H5_INDEX_NAME, H5_ITER_INC, NULL, RTS_op_func, (void *) &od);

    // printf("In freq, Out freq, %.1f %.1f\n",beam_freqs[freq_ind], od.freq_out );

    user_precision_t float_zenith_delays[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    beam_settings->MWAFEE_freqs[freq_ind] = (double)od.freq_out;

    //Setup zenith beams
    status = RTS_MWAFEEInit(woden_settings->hdf5_beam_path,
                            (double)od.freq_out,
                            &beam_settings->FEE_beam_zeniths[freq_ind],
                            float_zenith_delays);

    //Setup beam at desired pointing
    status = RTS_MWAFEEInit(woden_settings->hdf5_beam_path,
                            (double)od.freq_out,
                            &beam_settings->FEE_beams[freq_ind],
                            woden_settings->FEE_ideal_delays);

    H5Gclose(group);
  }

  H5Fclose(file);

  printf(" done.\n");

  return status;
}
