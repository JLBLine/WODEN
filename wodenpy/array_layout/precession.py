import numpy as np
import erfa
import palpy

# void RTS_mat_transpose(double rmat1[3][3], double rmat2[3][3]) {

#   int i, j;
#   for(i=0;i<3;++i) {
#     for(j=0;j<3;++j) {
#       rmat2[j][i] = rmat1[i][j];
#     }
#   }
# }

def RTS_mat_transpose(rmat1):

    rmat2 = np.empty(rmat1.shape)

    for i in range(3):
        for j in range(3):
            rmat2[j][i] = rmat1[i][j]
  
    return rmat2
  
def RTS_Precess_LST_Lat_to_J2000(lst_current, latitude_current,
                                 mjd):
    

    # Calculate a rotation matrix that accounts for precession and nutation
    # between the current modified julian date (mjd) and J2000
    #  palPrenut calls:
    #   - palPrec( 2000.0, palEpj(mjd), rmatp ); // form precession matrix: v_mean(mjd epoch) = rmatp * v_mean(J2000)
    #   - palNut( mjd, rmatn );                  // form nutation matrix: v_true(mjd epoch) = rmatn * v_mean(mjd epoch)
    #   - palDmxm( rmatn, rmatp, rmatpn );       // Combine the matrices:  pn = n x p
    
    rmatpn = palpy.prenut(2000, mjd)
    J2000_transformation = RTS_mat_transpose(rmatpn)
    
    # /**
    # ****************************************************************************
    # * Change the various coordinates to the J2000 mean system
    # ****************************************************************************
    # * palDcs2c   - convert the apparent direction to direction cosines
    # * palDmxv    - perform the 3-d forward unitary transformation: v2 = tmatpn * v1
    # * palDcc2s   - convert cartesian coordinates back to spherical coordinates (i.e. zenith in the J2000 mean system).
    # * palDranrm  - normalize into range 0-2 pi.
    # */

    # // Change the coordinates of the initial zenith
    v1 = palpy.dcs2c(lst_current, latitude_current)
    v2 = palpy.dmxv(J2000_transformation, v1)
    lst_J2000, latitude_J2000 = palpy.dcc2s(v2)
    lst_J2000 = palpy.dranrm(lst_J2000)

    return lst_J2000, latitude_J2000

if __name__ == "__main__":
    
    mjd = 59945
    
    lst_current = 0.0
    latitude_current = -30*(np.pi/180.0)
    
    lst_J2000, latitude_J2000 = RTS_Precess_LST_Lat_to_J2000(lst_current,
                                                             latitude_current,
                                                             mjd)
    
    print(lst_J2000, latitude_J2000)