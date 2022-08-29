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



void RTS_Precess_LST_Lat_to_J2000(double lst_current, double latitude_current,
                              double mjd,
                              double * lst_J2000, double * latitude_J2000){

    double rmatpn[3][3];
    double J2000_transformation[3][3];

    //Calculate a rotation matrix that accounts for precession and nutation
    //between the current modified julian date (mjd) and J2000
    // palPrenut calls:
    //  - palPrec( 2000.0, palEpj(mjd), rmatp ); // form precession matrix: v_mean(mjd epoch) = rmatp * v_mean(J2000)
    //  - palNut( mjd, rmatn );                  // form nutation matrix: v_true(mjd epoch) = rmatn * v_mean(mjd epoch)
    //  - palDmxm( rmatn, rmatp, rmatpn );       // Combine the matrices:  pn = n x p

    palPrenut(2000.0, mjd, rmatpn);
    RTS_mat_transpose( rmatpn, J2000_transformation );
    /**
    ****************************************************************************
    * Change the various coordinates to the J2000 mean system
    ****************************************************************************
    * palDcs2c   - convert the apparent direction to direction cosines
    * palDmxv    - perform the 3-d forward unitary transformation: v2 = tmatpn * v1
    * palDcc2s   - convert cartesian coordinates back to spherical coordinates (i.e. zenith in the J2000 mean system).
    * palDranrm  - normalize into range 0-2 pi.
    */

    double v1[3], v2[3];
    // Change the coordinates of the initial zenith
    palDcs2c(lst_current, latitude_current, v1);
    palDmxv(J2000_transformation, v1, v2);
    palDcc2s(v2, lst_J2000, latitude_J2000);
    * lst_J2000 = palDranrm(* lst_J2000);

    // printf("INSIDE Prec function %f %f %f\n",lst_current, latitude_current,
    //                                                           * lst_J2000 );
}

/**
As said in the include file/help, this is an RTS function to rotate the array
back to J2000. There are a number of lines commented out, which I don't think
I need, but am leaving here in case I do one day
*/
void RTS_PrecessXYZtoJ2000( array_layout_t *array_layout,
                            woden_settings_t *woden_settings) {

  // //The LST and latitude at the current observation time
  // double lst_current = (double)woden_settings->lst_base;
  // double latitude_current = woden_settings->latitude;

  int n_tile = array_layout->num_tiles;
  double J2000_transformation[3][3];
  double rmatpn[3][3];
  // double mjd = woden_settings->jd_date - 2400000.5;
  double lst_J2000;
  double X_epoch, Y_epoch, Z_epoch, X_prec, Y_prec, Z_prec;
  double lst_current, mjd_current;

  //Rotate the array positions for each time step - they have different
  //mjd dates and current epoch lsts so yield different XYZ over time
  for (int time_step = 0; time_step < woden_settings->num_time_steps; time_step++) {


    lst_J2000 = woden_settings->lsts[time_step];
    mjd_current = woden_settings->mjds[time_step];

    //Add on the angle accrued by current time step to the base LST
    lst_current = woden_settings->lst_obs_epoch_base + time_step*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;
    //Add half a time_res so we are sampling centre of each time step
    lst_current += 0.5*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;


    // double mjd = woden_settings->jd_date - 2400000.5;
    // double mjd_current = mjd + ((time_step + 0.5)*woden_settings->time_res)/(24.0*60.0*60.0);
    // printf("mjd %.16f\n",mjd_current );
    //double latitude_J2000;
    //
    // RTS_Precess_LST_Lat_to_J2000(lst_current, woden_settings->latitude,
    //                               mjd_current,
    //                               &lst_J2000, &latitude_J2000);
    //
    // printf("coords %.16f %.16f\n",lst_current, lst_J2000 );

    //Calculate a rotation matrix that accounts for precession and nutation
    //between the current modified julian date (mjd) and J2000
    // palPrenut calls:
    //  - palPrec( 2000.0, palEpj(mjd), rmatp ); // form precession matrix: v_mean(mjd epoch) = rmatp * v_mean(J2000)
    //  - palNut( mjd, rmatn );                  // form nutation matrix: v_true(mjd epoch) = rmatn * v_mean(mjd epoch)
    //  - palDmxm( rmatn, rmatp, rmatpn );       // Combine the matrices:  pn = n x p

    palPrenut(2000.0, mjd_current, rmatpn);
    RTS_mat_transpose( rmatpn, J2000_transformation );

    // printf("Rotation matrix\n");
    // printf("[%.8f %.8f %.8f\n",J2000_transformation[0][0],J2000_transformation[0][1], J2000_transformation[0][2] );
    // printf("[%.8f %.8f %.8f\n",J2000_transformation[1][0],J2000_transformation[1][1], J2000_transformation[1][2] );
    // printf("[%.8f %.8f %.8f]\n",J2000_transformation[2][0],J2000_transformation[2][1], J2000_transformation[2][2] );

    //Hold some things for me please

    //precess the current LST and latitude back to equivalent ones in the J2000 frame

    // RTS_Precess_LST_Lat_to_J2000(lst_current, woden_settings->latitude,
    //                              mjd_current,
    //                              &lst_J2000, &latitude_J2000);



    //Change things to be in J2000
    // woden_settings->lst_base = lst_J2000;
    // woden_settings->latitude = latitude_J2000;

    for (int st=0; st < n_tile; st++ ) {

      //Where do this XYZ for this time step start
      int st_offset = time_step*n_tile;

      // Calculate antennas positions in the J2000 frame
      X_epoch = array_layout->ant_X[st_offset + st];
      Y_epoch = array_layout->ant_Y[st_offset + st];
      Z_epoch = array_layout->ant_Z[st_offset + st];
      RTS_precXYZ( J2000_transformation, X_epoch, Y_epoch, Z_epoch,
                   lst_current, &X_prec, &Y_prec, &Z_prec, lst_J2000 );

      // Update stored coordinates
      array_layout->ant_X[st_offset + st] = X_prec;
      array_layout->ant_Y[st_offset + st] = Y_prec;
      array_layout->ant_Z[st_offset + st] = Z_prec;
    }
  } //time step loop
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
  woden_settings->num_ants = num_tiles;

  //malloc some arrays for holding array coords

  array_layout->latitude = woden_settings->latitude;
  array_layout->num_baselines = (array_layout->num_tiles*(array_layout->num_tiles-1)) / 2;
  woden_settings->num_baselines = array_layout->num_baselines;

  const int num_cross = woden_settings->num_baselines * woden_settings->num_time_steps * woden_settings->num_freqs;

  woden_settings->num_cross = num_cross;

  int num_autos;
  //If user asked for auto-correlations, set this number to total autos to work out
  //Otherwise stick it to zero
  if (woden_settings->do_autos){
    num_autos = woden_settings->num_ants * woden_settings->num_time_steps * woden_settings->num_freqs;
  } else {
    num_autos = 0;
  }

  woden_settings->num_autos = num_autos;
  woden_settings->num_visis = num_cross + num_autos;

  array_layout->ant_X = malloc( woden_settings->num_time_steps * array_layout->num_tiles * sizeof(double) );
  array_layout->ant_Y = malloc( woden_settings->num_time_steps * array_layout->num_tiles * sizeof(double) );
  array_layout->ant_Z = malloc( woden_settings->num_time_steps * array_layout->num_tiles * sizeof(double) );


  //When precessing the array back to J2000, we apply a different precession
  //based on mjd of each time step. So store enough copies of XYZ for
  //all time steps
  for (int time_step = 0; time_step < woden_settings->num_time_steps; time_step++) {
    for (int i = 0; i < array_layout->num_tiles; i++) {

      //Where do this XYZ for this time step start
      int st_offset = time_step*array_layout->num_tiles;

      //Convert to local X,Y,Z
      RTS_ENH2XYZ_local(array_layout->ant_east[i],
                        array_layout->ant_north[i],
                        array_layout->ant_height[i],
                        woden_settings->latitude_obs_epoch_base,
                        &(array_layout->ant_X[st_offset + i]),
                        &(array_layout->ant_Y[st_offset + i]),
                        &(array_layout->ant_Z[st_offset + i]));

    }
  }

  if (do_precession == 1) {
    RTS_PrecessXYZtoJ2000(array_layout, woden_settings);
  }

  array_layout->X_diff_metres = malloc( woden_settings->num_time_steps * array_layout->num_baselines * sizeof(double) );
  array_layout->Y_diff_metres = malloc( woden_settings->num_time_steps * array_layout->num_baselines * sizeof(double) );
  array_layout->Z_diff_metres = malloc( woden_settings->num_time_steps * array_layout->num_baselines * sizeof(double) );

  for (int time_step = 0; time_step < woden_settings->num_time_steps; time_step++) {

    //Where do these XYZ for this time step start
    int ant_off = array_layout->num_tiles*time_step;

    //Where do these baselines for this time step start
    int base_off = array_layout->num_baselines*time_step;

    int baseline_ind = 0;
    for (int ant1 = 0; ant1 < array_layout->num_tiles - 1; ant1++) {
      for (int ant2 = ant1 + 1; ant2 < array_layout->num_tiles; ant2++) {
        array_layout->X_diff_metres[base_off + baseline_ind] = array_layout->ant_X[ant_off + ant1] - array_layout->ant_X[ant_off + ant2];
        array_layout->Y_diff_metres[base_off + baseline_ind] = array_layout->ant_Y[ant_off + ant1] - array_layout->ant_Y[ant_off + ant2];
        array_layout->Z_diff_metres[base_off + baseline_ind] = array_layout->ant_Z[ant_off + ant1] - array_layout->ant_Z[ant_off + ant2];
        // printf("%.3f %.3f %.3f\n",array_layout->X_diff_metres[base_off + baseline_ind],
        //                           array_layout->Y_diff_metres[base_off + baseline_ind],
        //                           array_layout->Z_diff_metres[base_off + baseline_ind] );
        // printf("%.3f %.3f %.3f %.3f %.3f %.3f\n",array_layout->ant_X[ant_off + ant1], array_layout->ant_X[ant_off + ant2],
        //                                          array_layout->ant_Y[ant_off + ant1], array_layout->ant_Y[ant_off + ant2],
        //                                          array_layout->ant_Z[ant_off + ant1], array_layout->ant_Z[ant_off + ant2]);

        baseline_ind++;


      }
    }
  }

  return array_layout;
}
