#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <assert.h>

#include "FEE_primary_beam.h"
#include "legendre_polynomial.h"


static float _Complex *norm_fac=NULL;

static const double factorials[] = {
  1.0,
  1.0,
  2.0,
  6.0,
  24.0,
  120.0,
  720.0,
  5040.0,
  40320.0,
  362880.0,
  3628800.0,
  39916800.0,
  479001600.0,
  6227020800.0,
  87178291200.0,
  1.307674368e+12,
  2.0922789888e+13,
  3.55687428096e+14,
  6.40237370573e+15,
  1.21645100409e+17,
  2.43290200818e+18,
  5.10909421717e+19,
  1.12400072778e+21,
  2.58520167389e+22,
  6.20448401733e+23,
  1.55112100433e+25,
  4.03291461127e+26,
  1.08888694504e+28,
  3.04888344612e+29,
  8.84176199374e+30,
  2.65252859812e+32,
  8.22283865418e+33,
  2.63130836934e+35,
  8.68331761881e+36,
  2.9523279904e+38,
  1.03331479664e+40,
  3.7199332679e+41,
  1.37637530912e+43,
  5.23022617467e+44,
  2.03978820812e+46,
  8.15915283248e+47,
  3.34525266132e+49,
  1.40500611775e+51,
  6.04152630634e+52,
  2.65827157479e+54,
  1.19622220865e+56,
  5.50262215981e+57,
  2.58623241511e+59,
  1.24139155925e+61,
  6.08281864034e+62,
  3.04140932017e+64,
  1.55111875329e+66,
  8.06581751709e+67,
  4.27488328406e+69,
  2.30843697339e+71,
  1.26964033537e+73,
  7.10998587805e+74,
  4.05269195049e+76,
  2.35056133128e+78,
  1.38683118546e+80,
  8.32098711274e+81,
  5.07580213877e+83,
  3.14699732604e+85,
  1.9826083154e+87,
  1.26886932186e+89,
  8.24765059208e+90,
  5.44344939077e+92,
  3.64711109182e+94,
  2.48003554244e+96,
  1.71122452428e+98,
  1.197857167e+100,
  8.50478588568e+101,
  6.12344583769e+103,
  4.47011546151e+105,
  3.30788544152e+107,
  2.48091408114e+109,
  1.88549470167e+111,
  1.45183092028e+113,
  1.13242811782e+115,
  8.94618213078e+116
};

// /* function prototypes for LAPACK */
// void cgetrf_( int *m, int *n, float _Complex *a, int *lda, int *ipiv, int *info );
// void cgetri_( int *n, float _Complex *a, int *lda, int *ipiv, float _Complex *work, int *lwork, int *info );

static int count=0;
static int count1=0;
static int count2=0;
static int count3=0;


// /****************************
// reorder the delays from RTS order into engineer order
// (see comments in dipole_files.c)
// ****************************/
// void RTS_reorderDelays(float in[NUM_DIPOLES], float out[NUM_DIPOLES]) {
//     int i;
//     for (i=0; i< 4; i++) {
//         out[0+i*4] = in[3-i];
//         out[1+i*4] = in[7-i];
//         out[2+i*4] = in[11-i];
//         out[3+i*4] = in[15-i];
//     }
// }
//
//
// void RTS_reorderDelays2RTS(float in[NUM_DIPOLES], float out[NUM_DIPOLES]) {
//     int i;
//     for (i=0; i< 4; i++) {
//       out[3-i] = in[0+i*4];
//       out[7-i] = in[1+i*4];
//       out[11-i] = in[2+i*4];
//       out[15-i] = in[3+i*4];
//     }
// }

herr_t RTS_op_func (hid_t loc_id, const char *name, const H5L_info_t *info,
            void *operator_data) {
    herr_t          status, return_val = 0;
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
    return return_val;
}



/****************************
Find the nearest freq in the element pattern hdf file
****************************/
//#define H5FILE_NAME "/home/bpindor/code/MWAtools_pb/MWA_embedded_element_pattern_V02.h5"
//#define H5FILE_NAME "/home/bpindor/temp/MWAtools_pb/MWA_embedded_element_pattern_V02.h5"

int RTS_HDFBeamInit(char *h5filename, float freq_Hz, copy_primary_beam_t *pb, float *FEE_delays, int stn){

  hid_t       file, group;         /* handles */
  herr_t          status;
  H5O_info_t      infobuf;
  struct opdata   od;

  static float last_freq_in=0.0;
  static float last_freq_out=0.0;

  file = H5Fopen(h5filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  // if(norm_fac != NULL){
  //   memcpy(pb->norm_fac,norm_fac,MAX_POLS*sizeof(float _Complex));
  // }

  if(freq_Hz != last_freq_in) {
    /***************************************************************/
    /* Get nearest frequency by iterating over dataset in HDF file */

    //file = H5Fopen(H5FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
    status = H5Oget_info (file, &infobuf);
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

  float lam;
  float _Complex Vcplx[N_COPOL][NUM_DIPOLES];

  lam = VEL_LIGHT/(od.freq_out); // Should this be freq_out??

  for (int pol=0; pol < N_COPOL; pol++){
    for (int i=0; i<NUM_DIPOLES; i++) {
      float phase;
      phase = (-2.0*M_PI*DQ/lam)*FEE_delays[i];
      // phase = (-2.0*M_PI/lam)*FEE_delays[i]; // rts delays are already *DQ
      // if (debug) fprintf(fpd,"delay: %g, phase: %g\n",FEE_delays[i],phase);
      // Vcplx=amps[pol]*np.exp(1.0j*phases) #complex excitation col voltage
      // No dipole amplitudes implemented in RTS but cant flag dipoles

      // ant_i == -1 used to override dipole_flagging when setting reference_beam

      // if(options->disable_dipole_flags || stn==-1){
      //    Vcplx[pol][i] = cexp(1.0j * phase);
      // } else {
      //   // Dipole flags are written into dead_dip array in read_metafits.c:read_dipoles
      //   // Flagged dipoles have dead_dip==1, so set their amplitudes to abs(dead_dip-1)
      //
      //    Vcplx[pol][i] = fabs(options->dead_dip[(stn+pol)*NUM_DIPOLES + i] - 1.0) * cexp(1.0j * phase);
      // }

      //TODO implement dipole flagging and/or dipole amplitudes here
      //eg Vcplx[pol][i] = amps[pol][i] * cexp(1.0j * phase);

      Vcplx[pol][i] = cexp(1.0j * phase);
    } //end for dipole
  } // end for pol

    /* Get max length of spherical wave table */

  hid_t    dataset, dataspace;
  hsize_t dims_out[2];

  int n_ant = 16;
  int max_length=0;
  char table_name[80];
  int n_pols = 2;

  /* Determine the maximum required size for the beam mode data */
  /* Should check Y pol as well? */

  for(int i = 1; i<n_ant+1;i++){
    sprintf(table_name,"X%d_%d",i,(int)od.freq_out);
    dataset = H5Dopen(file,table_name, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    status  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if(dims_out[1]/2 > max_length){
       max_length = dims_out[1]/2; // Complex number
    }
    H5Dclose(dataset);
  }

  double _Complex *Q1, *Q2;

  /* Allocate space for two polarizations */

  pb->Q1 = (double _Complex **)malloc(n_pols*sizeof(double _Complex *));
  pb->Q2 = (double _Complex **)malloc(n_pols*sizeof(double _Complex *));
  pb->p_P = (double _Complex **)malloc(n_pols*sizeof(double _Complex *));
  pb->p_T = (double _Complex **)malloc(n_pols*sizeof(double _Complex *));
  pb->M = (double **)malloc(n_pols*sizeof(double *));
  pb->N = (double **)malloc(n_pols*sizeof(double *));

  for(int i=0;i<n_pols;i++){
    pb->Q1[i]= (double _Complex *)malloc(sizeof(double _Complex)*max_length);
    pb->Q2[i]= (double _Complex *)malloc(sizeof(double _Complex)*max_length);
    pb->p_P[i]= (double _Complex *)malloc(sizeof(double _Complex)*max_length);
    pb->p_T[i]= (double _Complex *)malloc(sizeof(double _Complex)*max_length);
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

  //    double Q_modes_all[2046][3];
  double Q_modes_all[3][2046];

  dataset = H5Dopen(file,"modes", H5P_DEFAULT);

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
        sprintf(table_name,"X%d_%d",ant_i,(int)od.freq_out);
      } else {
        sprintf(table_name,"Y%d_%d",ant_i,(int)od.freq_out);
      }
      dataset = H5Dopen(file,table_name, H5P_DEFAULT);
      dataspace = H5Dget_space(dataset);
      status  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
      Q_all = malloc(sizeof(double) * dims_out[0] * dims_out[1]);
      status = H5Dread(dataset, H5T_NATIVE_DOUBLE,
      	     H5S_ALL, H5S_ALL, H5P_DEFAULT, Q_all);
      H5Dclose(dataset);

      int my_len;

      my_len = dims_out[1]; // max(Q_all.shape)

      // Get Q_modes for this antenna

      // Q_modes = Q_modes_all[0:my_len,:]
      // 3 is first dimension of Q_modes_all (why??)

      Q_modes = malloc(sizeof(double) * my_len * 3);

      for(int q_i =0; q_i < my_len; q_i++){
        for(int q_j = 0; q_j < 3; q_j++){
          //	  Q_modes[q_i + my_len*q_j] = Q_modes_all[q_i][q_j];
          Q_modes[q_i + my_len*q_j] = Q_modes_all[q_j][q_i];
        }
      }

      // convert Qall to M, N, Q1, Q2 vectors for processing

      // find s=1 and s=2 indices
      // only find s1 and s2 for this antenna

      int *s1, *s2;
      double *M_vec, *N_vec;
      int nM,nN;

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
        pb->nmax = Nmax[pol];
        assert(nM==nN);
        pb->nMN = nM;
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

    } // loop over ant_i

    for(int i = 0; i < pb->nMN; i++){
      pb->p_T[pol][i] = cpow(1.0j,pb->N[pol][i]);
      pb->p_P[pol][i] = cpow(1.0j,(pb->N[pol][i]+1.0));
    }

  } // loop over XY pols

  if(Nmax[0] != Nmax[1]){
    printf("Beam Nmax is not equal for XY pols\n");
    exit(1);
  }

  free(Q1);
  free(Q2);

  H5Fclose(file);

  pb->m_range = malloc((2*pb->nmax + 1)*sizeof(float) );

  // printf("Managed to do m_range the Malloc\n");

  for (int m = 0; m < 2*pb->nmax + 1; m++) {
    pb->m_range[m] = -(float)pb->nmax + (float)m;
    // printf("%d %d %.2f\n", m, pb->nmax, pb->m_range[m] );
  }

  return 0;

}

void RTS_P1SIN(int nmax, float theta, double _Complex **P_sin, double _Complex **P1){



  *P_sin =  (double _Complex *)malloc(sizeof(double _Complex)*(nmax*nmax + 2*nmax));
  *P1 =  (double _Complex *)malloc(sizeof(double _Complex)*(nmax*nmax + 2*nmax));

  double u = cos(theta); // array??
  double sin_th = sin(theta);
  double delu = 1.0e-6;
  double pm_in[1];

  /* Test Legendre Polynomial values */

  double *pm_vals;

  /* Optimization testing */

  static int last_nmax=0;
  static int count4 = 0;

  /***** Syntax Testing *********/
  /*
  for(int orders=0;orders<=3;orders++){
    pm_vals = pm_polynomial_value (1, 3, orders, u); // as python lpmv
    printf("P3m(%3.2lf) %3.2lf m=%d\n",u[0],pm_vals[3],orders);
    free(pm_vals);
  }
  */
  /*******************************/

  double *P,*Pm1;
  double _Complex *Pm_sin;

  if(nmax != last_nmax){
    count4++;
    last_nmax = nmax;
  }

  int n_vals;
  double *all_vals;

  // Total number of Pmn values
  // P00 not used but just included to make indexing simpler
  n_vals = ((nmax+1) * (nmax+2) / 2);

  all_vals = malloc(n_vals * sizeof(double));
  int m_incr=0; // for incrementing array index

  pm_in[0] = u;

  for(int m=0;m<=nmax;m++){
    pm_vals = pm_polynomial_value(1,nmax,m,pm_in);
    for(int i=m;i<=nmax;i++){
      if(!(i==0 && m==0)){
	all_vals[(i-m) + m_incr] = pm_vals[i];
      }
    }
    m_incr+=nmax-m+1;
  }

  m_incr=0;
  int index;

  for(int n=1;n<=nmax;n++){
    P = malloc(sizeof(double )*(n+1));
    Pm1 = malloc(sizeof(double )*(n+1));
    Pm_sin = malloc(sizeof(double _Complex)*(n+1));
    pm_in[0] = u;

    m_incr=0;
    for(int order=0;order<=n;order++){
      index = n + m_incr;
      //      pm_vals = pm_polynomial_value(1,n,order,pm_in);
      //  P[order] = pm_vals[n];
      P[order] = all_vals[index];
      if(order>0){
	//	Pm1[order-1]=pm_vals[n];
	      Pm1[order-1] = all_vals[index];
      }
      if(order == n){
	      Pm1[order] = 0.0;
      }
      //free(pm_vals);
      Pm_sin[order] = 0.0;
      m_incr += nmax-order;
    }

    if(u==1.0){
      pm_in[0] = u-delu;
      pm_vals = pm_polynomial_value(1,n,0,pm_in);
      Pm_sin[1] = -(P[0] - pm_vals[n])/delu;
      free(pm_vals);
    } else if (u == -1.0){
      pm_in[0] = u-delu;
      pm_vals = pm_polynomial_value(1,n,0,pm_in);
      Pm_sin[1] = -(pm_vals[n] - P[0])/delu;
      free(pm_vals);
    } else {
      for(int order=0;order<=n;order++){
	       Pm_sin[order] = P[order]/sin_th;
      }
    }

    int ind_start,ind_stop;

    ind_start = (n-1)*(n-1) + 2*(n-1);
    ind_stop = n*n + 2*n;

    for(int i=ind_start;i<ind_stop;i++){
      int index = i-ind_start;
      if(index<n){
	(*P_sin)[i] = Pm_sin[n-index];
	(*P1)[i] = Pm1[n-index];
      } else {
	(*P_sin)[i] = Pm_sin[index-n];
	(*P1)[i] = Pm1[index-n];
      }
    }

    free(P);
    free(Pm1);
    free(Pm_sin);

  } // loop over n

  free(all_vals);

}
//
// // Version of get_FF which uses cached values to avoid recalculating values
// // when repeatedly for same phi, theta values
//
int RTS_get_FF2(float phi, float theta, copy_primary_beam_t *pb, float _Complex result[4], int scaling, int clean_up)
{

  int n_pols = 2;
  static double _Complex Jones[2][2];
  static double _Complex prevJones[2][2];
  int nmax;
  int beam_reset=0;
  static int nMN=0;
  static float last_phi=-100.0;
  static float last_theta=-100.0;
  static float prev_phi=-100.0;
  static float prev_theta=-100.0;

  static double **C_MN=NULL;
  static double **M_absM=NULL;
  static double _Complex **emn_T=NULL;
  static double _Complex **emn_P=NULL;
  double factor1, factor2;
  int index;
  static double _Complex **phi_comp;


  if(clean_up){
    // release cached values
    for(int pol=0;pol<n_pols;pol++){
      free(C_MN[pol]);
      free(M_absM[pol]);
      free(phi_comp[pol]);
      free(emn_T[pol]);
      free(emn_P[pol]);
    }
    free(C_MN);
    free(M_absM);
    free(phi_comp);
    free(emn_T);
    free(emn_P);
    last_phi=-1.0;
    last_theta=-1.0;
    prev_phi=-1.0;
    prev_theta=-1.0;
    nMN = 0;
    /*
    printf("Count0 %d\n",count);
    printf("Count1 %d\n",count1);
    printf("Count2 %d\n",count2);
    printf("Count3 %d\n",count3);
    */

    return 0;
  }

  index=0;
  phi = (M_PI/2.0 - phi);
  if(phi < 0) phi += 2.0*M_PI;

  if(theta == prev_theta && phi == prev_phi && nMN == pb->nMN){
    result[3] = (float _Complex) prevJones[0][0];
    result[2] = (float _Complex) prevJones[0][1];
    result[1] = (float _Complex) prevJones[1][0];
    result[0] = (float _Complex) prevJones[1][1];

    prevJones[0][0] = Jones[0][0];
    prevJones[0][1] = Jones[0][1];
    prevJones[1][0] = Jones[1][0];
    prevJones[1][1] = Jones[1][1];

    Jones[0][0] = result[3];
    Jones[0][1] = result[2];
    Jones[1][0] = result[1];
    Jones[0][0] = result[0];

    float tmpt,tmpp;
    tmpt = last_theta;
    last_theta = prev_theta;
    prev_theta = tmpt;

    tmpp = last_phi;
    last_phi = prev_phi;
    prev_phi = tmpp;

    return 0;
  }


  if(theta == last_theta && phi == last_phi && nMN == pb->nMN){
    result[3] = (float _Complex) Jones[0][0];
    result[2] = (float _Complex) Jones[0][1];
    result[1] = (float _Complex) Jones[1][0];
    result[0] = (float _Complex) Jones[1][1];

    return 0;
  } else {
    prevJones[0][0] = Jones[0][0];
    prevJones[0][1] = Jones[0][1];
    prevJones[1][0] = Jones[1][0];
    prevJones[1][1] = Jones[1][1];

    prev_theta = last_theta;
    prev_phi = last_phi;
  }


  nmax = pb->nmax;

  /* C_MN and M_absM need to be recalculated only if beam model changes */

  if(nMN==0 || nMN != pb->nMN){
    beam_reset=1;
    count1++;
    if(nMN == 0){
      C_MN = malloc(n_pols*sizeof(double *));
      M_absM = malloc(n_pols*sizeof(double *));
      phi_comp = malloc(n_pols*sizeof(double _Complex *));
      emn_T = malloc(n_pols*sizeof(double _Complex *));
      emn_P = malloc(n_pols*sizeof(double _Complex *));
    } else {
      for(int pol=0;pol<n_pols;pol++){
	       free(C_MN[pol]);
	       free(M_absM[pol]);
	       free(phi_comp[pol]);
	       free(emn_T[pol]);
	       free(emn_P[pol]);
      }
    }
    nMN = pb->nMN; // Need to consider other static values
    for(int pol=0;pol<n_pols;pol++){
      C_MN[pol] = malloc(sizeof(double) * pb->nMN);
      M_absM[pol] = malloc(sizeof(double) * pb->nMN);
      phi_comp[pol] = malloc(sizeof(double _Complex) * pb->nMN);
      emn_T[pol] = malloc(sizeof(double _Complex) * pb->nMN);
      emn_P[pol] = malloc(sizeof(double _Complex) * pb->nMN);
    }


    for(int pol=0;pol<n_pols;pol++){

      for(int i=0;i<nMN;i++){

	index = (int) (pb->N[pol][i] - fabs(pb->M[pol][i]));
	if(index >=80){
	  printf("Maximum factorial exceeded\n");
	  exit(1);
	}
	factor1 = factorials[index];
	index = (int) (pb->N[pol][i] + fabs(pb->M[pol][i]));
	if(index >=80){
	  printf("Maximum factorial exceeded\n");
	  exit(1);
	}
	factor2 = factorials[index];

	C_MN[pol][i] = sqrt(0.5 * (2 * pb->N[pol][i] + 1) * factor1 / factor2);
	if(pb->M[pol][i] == 0){
	  M_absM[pol][i] = 1;
	} else {
	  M_absM[pol][i] = -(pb->M[pol][i] / fabs(pb->M[pol][i]));
	}
	M_absM[pol][i] = pow(M_absM[pol][i],pb->M[pol][i]);
      }
    }

  } // End Calculation of C_MN and M_absM


  /* phi_comp needs to be recalculated if phi has changed or if beam model changes */


  if(phi != last_phi || beam_reset){
    count2++;
    for(int pol=0;pol<n_pols;pol++){
      for(int i=0;i<nMN;i++){
	phi_comp[pol][i] = cexp(1.0j * pb->M[pol][i] * phi) * C_MN[pol][i] * M_absM[pol][i] / sqrt(pb->N[pol][i] * (pb->N[pol][i] + 1.0));
      }
    }
  }

  /* Now phi_comp has been recomputed or stored value can be used from previous call */

  double _Complex *P1, *P_sin;

  if(theta != last_theta || beam_reset){


    count3++;
    RTS_P1SIN(nmax,theta,&P_sin,&P1);

    for(int pol=0;pol<n_pols;pol++){

      double u = cos(theta);

      for(int i=0;i<nMN;i++){
	       emn_T[pol][i] = (pb->p_T[pol][i])*(P_sin[i]*(fabs(pb->M[pol][i])*pb->Q2[pol][i]*u - pb->M[pol][i] * pb->Q1[pol][i]) + pb->Q2[pol][i] * P1[i]);
	       emn_P[pol][i] = (pb->p_P[pol][i])*(P_sin[i]*(pb->M[pol][i]*pb->Q2[pol][i] - fabs(pb->M[pol][i])*pb->Q1[pol][i]*u) - pb->Q1[pol][i] * P1[i]);

      }

    }
    free(P1);
    free(P_sin);

  } // END (theta != last_theta)

  // If either phi or theta have changed, recalculate Jones elements

  double _Complex Sigma_P, Sigma_T;

  if(theta != last_theta || phi != last_phi || beam_reset){

    for(int pol=0;pol<n_pols;pol++){

      Sigma_P = Sigma_T = 0.0 + 0.0j;

      for(int i=0;i<nMN;i++){
	       Sigma_P += emn_P[pol][i] * phi_comp[pol][i];
	       Sigma_T += emn_T[pol][i] * phi_comp[pol][i];
      }

      Jones[pol][0] = Sigma_T;
      Jones[pol][1] = Sigma_P;

    }
  }

  last_phi = phi;
  last_theta = theta;


  // RTS polarizations are (of course) opposite orientation to FEKO

  result[3] = (float _Complex) Jones[0][0];
  result[2] = (float _Complex) Jones[0][1];
  result[1] = (float _Complex) Jones[1][0];
  result[0] = (float _Complex) Jones[1][1];

  /*
  prevJones[0][0] = Jones[0][0];
  prevJones[0][1] = Jones[0][1];
  prevJones[1][0] = Jones[1][0];
  prevJones[1][1] = Jones[1][1];
  */

  if(scaling==1) {

      if(pb->norm_fac == NULL){
          printf("ERROR: Beam normalization must be set with getHDFBeamNormalization before scaling\n");
	        exit(1);
	    }
      // Zenith normalized

      Jones[0][0] /= norm_fac[3];
      Jones[0][1] /= norm_fac[2];
      Jones[1][0] /= norm_fac[1];
      Jones[1][1] /= norm_fac[0];

  }

  result[3] = (float _Complex) Jones[0][0];
  result[2] = (float _Complex) Jones[0][1];
  result[1] = (float _Complex) Jones[1][0];
  result[0] = (float _Complex) Jones[1][1];

  return 0;

}
//
// /*****************************************************************************
//  * This function is part of the MWA-LFD realtime software library
//  !----------------------------------------------------------------------------
//  !
//  ! NAME:        getJonesSphHarm
//  ! PURPOSE:     Calc the Jones matrix response for a simulated dipole
//  ! ARGUMENTS:   freq_Hz reference frequency in Hz
//  !              az: azimuth angle (radian)
//  !              za: zenith angle (radian)
//  !              pb: primary beam data structure with dipole height set
//  !              result is an array of dimension 2x2 for complex array X and Y mapping unit vecs
//  !              in AZ and ZA onto unit vectors in East and North.
//  !              order of output products is:
//  !                  0: ZA onto E
//  !                  1: AZ onto E
//  !                  2: ZA onto W
//  !                  3: AZ onto W
//  !
//  *****************************************************************************/
int RTS_getJonesSphHarm(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4],int scaling) {
    int res=0;

    // Straight port of get_FF from python version
    //get_FF(az,za,pb,result,scaling);

    // Version which caches previous vaules
    RTS_get_FF2(az,za,pb,result,scaling,0);

    res = 0;

    // Counter for hand profiling
    count++;

    return res;
}
//
int RTS_getJonesDipole(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4],int scaling) {
    // int res=0;

    return RTS_getJonesSphHarm(freq_Hz, az, za, pb, result,scaling);

    // switch (pb->type) {
    //     case ANT_CROSSED_DIPOLES_ON_GROUNDPLANE :
    //         //if(debug) fprintf(fpd,"%s: Type %d.\n",__func__,pb->type);
    //         return getJonesShortDipole(freq_Hz, az, za, pb, result) ;
    //         break;
    //     case ANT_SIMULATED :
    //
	  // break;
    //
    //     default:
	  //        fprintf(stderr,"%s: Unknown dipole type %d\n",__func__,pb->type);
	  //        res=1;
    //     }
    // return res;
}
//
// /*****************************************************************************
//  * This function is part of the MWA-LFD realtime software library
//  !----------------------------------------------------------------------------
//  !
//  ! NAME:        getTileResponse
//  ! PURPOSE:     Calc the Jones matrix response for a tile
//  ! ARGUMENTS:   freq_Hz reference frequency in Hz
//  !              az: azimuth angle (radian)
//  !              za: zenith angle (radian)
//  !              pb: primary beam data structure with dipole height set
//  !              result is an array of dimension 2x2 for complex array X and Y mapping unit vecs
//  !              in AZ and ZA onto unit vectors in East and North.
//  !              order of output products is:
//  !                  0: ZA onto E
//  !                  1: AZ onto E
//  !                  2: ZA onto W
//  !                  3: AZ onto W
//  !
// from RTS in InstrumentCalibration.c:
//     might need to take into account local zenith effects
//     ground_plane = 2.0*sin(2.0*pi*dpl_hgt/lambda*cos(slaDsep(beam->coord_pnt.ha,beam->coord_pnt.dec,ha,dec)));
//     BP 2019: This is the actual beam function called inside the RTS. The rotation term
//              takes the beam from the alt-az coordinates to ra-dec coordinates.
//  *****************************************************************************/
int RTS_getTileResponse(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4], int scaling, float rotation) {

    int res=0;
    float _Complex prerot[4];


    // BP 2019: The FEE model section

    res = RTS_getJonesDipole(freq_Hz, az, za, pb, prerot, scaling);

    // The FEE beam is given relative to local az,za but the RTS wants J
    // relative to local ha,dec. So we need to rotate through the paralactic
    // angle (actually para + pi/2)

    if(rotation!=0.0){

      float crot,srot;

      crot = cos(rotation);
      srot = sin(rotation);

      //cgemm2x2??
      // The pols have already been flipped in get_FF, so need to flip them
      // back when calculating rotation (ugh)
      result[3] = prerot[3]*crot + prerot[2]*srot;
      result[2] = -prerot[3]*srot + prerot[2]*crot;
      result[1] = prerot[1]*crot + prerot[0]*srot;
      result[0] = -prerot[1]*srot + prerot[0]*crot;

    } else {

      for(int i=0;i<MAX_POLS;i++){
        result[i] = prerot[i];
      }
    }

    return res;
}
//
void RTS_freeHDFBeam(copy_primary_beam_t *pb){

  if(pb==NULL) return;

  int n_pols = 2;

  for(int n=0; n < n_pols; n++){
    free(pb->Q1[n]);
    free(pb->Q2[n]);
    free(pb->p_P[n]);
    free(pb->p_T[n]);
    free(pb->M[n]);
    free(pb->N[n]);
  }

  free(pb->Q1);
  free(pb->Q2);
  free(pb->p_P);
  free(pb->p_T);
  free(pb->M);
  free(pb->N);

  if(pb->parameters!=NULL) free(pb->parameters);


}
//
// int RTS_getHDFBeamNormalization(char *h5filename, float freq, copy_primary_beam_t *pb){
//
//   /* Get zenith response of zenith pointing */
//   /* This factor remains the same for a given frequency, so only need to calculate once for
//      any and all beams in a given rts instance */
//   /* From primary_beam_full_EE.py:
//
//   Calculate normalisation factors for the Jones vector for this
//         ApertureArray object. For MWA, these are at the zenith of a zenith pointed beam,
//         which is the maximum for all beam pointings.
//         The FEKO simulations include all ph angles at za=0. These are not redundant,
//         and the ph value determines the unit vector directions of both axes.
//         For the E-W dipoles, the projection of the theta unit vec will be max when
//         pointing east, i.e. when ph_EtN=0 (ph_NtE=90). For the phi unit vec,
//         this will be when ph_EtN=-90 or 90 (ph_NtE=180 or 0: we use 180)
//         For the N-S dipoles, projection of ZA onto N-S is max az ph_EtN=90 (ph_NtE=0) and
//         proj of ph onto N-S is max when ph_EtN=0 (ph_NtE=90
//
//   BP: Not sure this is okay for sky x y as defined by the RTS though
//   */
//
//   int res=0;
//
//   // We dont want any dipole flagging when calculating the zenith values
//   // Creating a options structure will set the dipole flagging to FALSE
//
//   copy_primary_beam_t *primary_beam;
//   primary_beam = malloc(sizeof(copy_primary_beam_t));
//
//   if(norm_fac == NULL){
//
//     float zenith_delays[NUM_DIPOLES] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//
//     float _Complex j_norm[MAX_POLS];
//     float delays_copy[NUM_DIPOLES];
//
//     norm_fac = malloc(MAX_POLS*sizeof(float _Complex));
//
//     // for(unsigned int k=0; k<NUM_DIPOLES; k++){
//     //   delays_copy[k] = ->delay[k];
//     //   ->delay[k] = 0.0;
//     // }
//
//     int scaling = 0; // We want the 'true' value at zenith
//     res = RTS_HDFBeamInit(h5filename, freq, primary_beam, zenith_delays, 0);
//
//     // BP: Not quite sure about the pol ordering here RTS v FEKO
//     //     although in pratical term the differences are very small
//
//     float rotation=0.0;
//
//     res = RTS_getTileResponse(freq, (M_PI/2.0), 0.0, primary_beam, j_norm ,scaling ,rotation);
//     pb->norm_fac[0]=j_norm[0];
//     res = RTS_getTileResponse(freq, 0.0, 0.0, primary_beam, j_norm ,scaling ,rotation);
//     pb->norm_fac[1]=j_norm[1];
//     res = RTS_getTileResponse(freq, M_PI, 0.0, primary_beam, j_norm ,scaling ,rotation);
//     pb->norm_fac[2]=j_norm[2];
//     res = RTS_getTileResponse(freq, (M_PI/2.0), 0.0, primary_beam, j_norm ,scaling ,rotation);
//     pb->norm_fac[3]=j_norm[3];
//
//     RTS_freeHDFBeam(primary_beam);
//     RTS_HDFBeamCleanUp();
//
//   }
//
//   return res;
// //
// }
//
void RTS_HDFBeamCleanUp(){

  // Makes special call to get_FF2 which release cached arrays
  RTS_get_FF2(0,0,NULL,NULL,0,1);

}
