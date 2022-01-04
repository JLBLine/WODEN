/*******************************************************************************
*  Methods to read read/generate array layouts, and move the array back to J2000.
*  @author J.L.B. Line
*
*  Please see documentation in ../include/array_layout.h or online at
*  https://woden.readthedocs.io/en/latest/index.html
*******************************************************************************/
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <json.h>
#include <pal.h>

void RTS_ENH2XYZ_local(double E, double N, double H,
                       double lat,
                       double *X, double *Y, double *Z) {
  double sl, cl;
  sl = sin(lat);
  cl = cos(lat);

  *X = -N*sl + H*cl;
  *Y = E;
  *Z = N*cl + H*sl;
}

void RTS_precXYZ(double rmat[3][3], double x, double y, double z, double lmst,
         double *xp, double *yp, double *zp, double lmst2000) {

  double sep, cep, s2000, c2000;
  double xpr, ypr, zpr, xpr2, ypr2, zpr2;

  sep = sin(lmst);
  cep = cos(lmst);
  s2000 = sin(lmst2000);
  c2000 = cos(lmst2000);

  /* rotate to frame with x axis at zero RA */
  xpr = cep*x - sep*y;
  ypr = sep*x + cep*y;
  zpr = z;

  /* apply the rotation matrix to account for precession/nutation */
  xpr2 = (rmat[0][0])*xpr + (rmat[0][1])*ypr + (rmat[0][2])*zpr;
  ypr2 = (rmat[1][0])*xpr + (rmat[1][1])*ypr + (rmat[1][2])*zpr;
  zpr2 = (rmat[2][0])*xpr + (rmat[2][1])*ypr + (rmat[2][2])*zpr;

  /* rotate back to frame with xp pointing out at lmst2000 */
  *xp = c2000*xpr2 + s2000*ypr2;
  *yp = -s2000*xpr2 + c2000*ypr2;
  *zp = zpr2;

}

/**************************
//RTS matrix transpose function
***************************/
void RTS_mat_transpose(double rmat1[3][3], double rmat2[3][3]) {

  int i, j;
  for(i=0;i<3;++i) {
    for(j=0;j<3;++j) {
      rmat2[j][i] = rmat1[i][j];
    }
  }
}

/**
As said in the include file/help, this is an RTS function to rotate the array
back to J2000. There are a number of lines commented out, which I don't think
I need, but am leaving here in case I do one day
*/
void RTS_PrecessXYZtoJ2000( array_layout_t *array_layout,
                       woden_settings_t *woden_settings) {

  //The LST and latitude at the current observation time
  double lst_current = (double)woden_settings->lst_base;
  double latitude_current = woden_settings->latitude;

  int n_tile = array_layout->num_tiles;

  double rmatpn[3][3];
  double mjd = woden_settings->jd_date - 2400000.5;
  double J2000_transformation[3][3];

  //Calculate a rotation matrix that accounts for precession and nutation
  //between the current modified julian date (mjd) and J2000
  // palPrenut calls:
  //  - palPrec( 2000.0, palEpj(mjd), rmatp ); // form precession matrix: v_mean(mjd epoch) = rmatp * v_mean(J2000)
  //  - palNut( mjd, rmatn );                  // form nutation matrix: v_true(mjd epoch) = rmatn * v_mean(mjd epoch)
  //  - palDmxm( rmatn, rmatp, rmatpn );       // Combine the matrices:  pn = n x p

  palPrenut(2000.0, mjd, rmatpn);
  RTS_mat_transpose( rmatpn, J2000_transformation );

  // printf("Rotation matrix\n");
  // printf("[%.8f %.8f %.8f\n",J2000_transformation[0][0],J2000_transformation[0][1], J2000_transformation[0][2] );
  // printf("[%.8f %.8f %.8f\n",J2000_transformation[1][0],J2000_transformation[1][1], J2000_transformation[1][2] );
  // printf("[%.8f %.8f %.8f]\n",J2000_transformation[2][0],J2000_transformation[2][1], J2000_transformation[2][2] );

  //Hold some things for me please
  double X_epoch, Y_epoch, Z_epoch, X_prec, Y_prec, Z_prec;
  double v1[3], v2[3], lst_J2000, latitude_J2000;

  /**
  ****************************************************************************
  * Change the various coordinates to the J2000 mean system
  ****************************************************************************
  * palDcs2c   - convert the apparent direction to direction cosines
  * palDmxv    - perform the 3-d forward unitary transformation: v2 = tmatpn * v1
  * palDcc2s   - convert cartesian coordinates back to spherical coordinates (i.e. zenith in the J2000 mean system).
  * palDranrm  - normalize into range 0-2 pi.
  */

  // Change the coordinates of the initial zenith
  palDcs2c(lst_current, latitude_current, v1);
  palDmxv(J2000_transformation, v1, v2);
  palDcc2s(v2, &lst_J2000, &latitude_J2000);
  lst_J2000 = palDranrm(lst_J2000);

  //Change things to be in J2000
  woden_settings->lst_base = lst_J2000;
  woden_settings->latitude = latitude_J2000;

  for (int st=0; st < n_tile; st++ ) {

    // Calculate antennas positions in the J2000 frame
    X_epoch = array_layout->ant_X[st];
    Y_epoch = array_layout->ant_Y[st];
    Z_epoch = array_layout->ant_Z[st];
    RTS_precXYZ( J2000_transformation, X_epoch, Y_epoch, Z_epoch, lst_current, &X_prec, &Y_prec, &Z_prec, lst_J2000 );

    // Update stored coordinates
    array_layout->ant_X[st] = X_prec;
    array_layout->ant_Y[st] = Y_prec;
    array_layout->ant_Z[st] = Z_prec;
  }
} // RTS_PrecessXYZtoJ2000

array_layout_t * calc_XYZ_diffs(woden_settings_t *woden_settings,
                                int do_precession){

  array_layout_t * array_layout;
  array_layout = malloc( sizeof(array_layout_t) );

  int num_tiles = 0;

  FILE *fp=NULL;
  char line[BUFSIZ];

  if ((fp=fopen(woden_settings->array_layout_file_path,"r"))==NULL) {
    printf("Reading of array_layout file:\n %s\nhas failed", woden_settings->array_layout_file_path);
    exit(1);
  }

  //gcc 7.5.0 on my desktop will not perform realloc later in the code
  //unless I do an initial malloc here
  array_layout->ant_east = malloc(sizeof(double));
  array_layout->ant_north = malloc(sizeof(double));
  array_layout->ant_height = malloc(sizeof(double));

  while(fgets(line,BUFSIZ,fp) != NULL) {

    num_tiles += 1;

    array_layout->ant_east = realloc(array_layout->ant_east,
                                    sizeof(double)*num_tiles);
    array_layout->ant_north = realloc(array_layout->ant_north,
                                    sizeof(double)*num_tiles);
    array_layout->ant_height = realloc(array_layout->ant_height,
                                    sizeof(double)*num_tiles);

    sscanf( line, "%lf %lf %lf", &array_layout->ant_east[num_tiles-1],
                                 &array_layout->ant_north[num_tiles-1],
                                 &array_layout->ant_height[num_tiles-1] );
  }

  array_layout->num_tiles = num_tiles;

  //malloc some arrays for holding array coords

  array_layout->latitude = woden_settings->latitude;
  array_layout->num_baselines = (array_layout->num_tiles*(array_layout->num_tiles-1)) / 2;
  woden_settings->num_baselines = array_layout->num_baselines;

  array_layout->ant_X = malloc( array_layout->num_tiles * sizeof(double) );
  array_layout->ant_Y = malloc( array_layout->num_tiles * sizeof(double) );
  array_layout->ant_Z = malloc( array_layout->num_tiles * sizeof(double) );

  for (int i = 0; i < array_layout->num_tiles; i++) {
    //Convert to local X,Y,Z
    RTS_ENH2XYZ_local(array_layout->ant_east[i], array_layout->ant_north[i],
                      array_layout->ant_height[i],
                      array_layout->latitude,
                      &(array_layout->ant_X[i]), &(array_layout->ant_Y[i]),
                      &(array_layout->ant_Z[i]));
  }

  if (do_precession == 1) {
    RTS_PrecessXYZtoJ2000(array_layout, woden_settings);
    printf("After precession, latitude of array is set to %.16f deg\n",
              woden_settings->latitude/DD2R);
  } else {
    printf("No precession, latitude of array is set to %.16f  deg\n",
              woden_settings->latitude/DD2R);
  }

  array_layout->X_diff_metres = malloc( array_layout->num_baselines * sizeof(double) );
  array_layout->Y_diff_metres = malloc( array_layout->num_baselines * sizeof(double) );
  array_layout->Z_diff_metres = malloc( array_layout->num_baselines * sizeof(double) );

  int baseline_ind = 0;
  for (int ant1 = 0; ant1 < array_layout->num_tiles - 1; ant1++) {
    for (int ant2 = ant1 + 1; ant2 < array_layout->num_tiles; ant2++) {
      array_layout->X_diff_metres[baseline_ind] = array_layout->ant_X[ant1] - array_layout->ant_X[ant2];
      array_layout->Y_diff_metres[baseline_ind] = array_layout->ant_Y[ant1] - array_layout->ant_Y[ant2];
      array_layout->Z_diff_metres[baseline_ind] = array_layout->ant_Z[ant1] - array_layout->ant_Z[ant2];
      baseline_ind++;
    }
  }

  return array_layout;
}
