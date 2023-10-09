import numpy as np
import palpy as pal
from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation
from astropy import units as u

DD2R = np.pi/180.0

VELC  = 299792458.0
MWA_LAT = -26.703319
MWA_LONG = 116.67081524

def calc_b(phi_simple, num_mult):
    """Given the target angle `phi_simple` (radians), and
    an integer multiplier `num_mult`, find a baseline
    length `b` that results in the same sine/cosine
    angle"""

    if phi_simple == 0:

        b = 2*np.pi*num_mult

    else:
        b = (phi_simple + 2*np.pi*num_mult) / phi_simple

    return b

def enh2xyz(east, north, height, latitiude):
    """Convert local e,n,h to local x,y,z coords"""
    sl = np.sin(latitiude)
    cl = np.cos(latitiude)
    X = -north*sl + height*cl
    Y = east
    Z = north*cl + height*sl

    return X, Y, Z

def get_uvw(xdiff, ydiff, zdiff, dec0, ha0):
    """Get the u,v,w given a baseline length in local X, Y, Z coords,
    and a phase centre ha0, dec0"""
    u = np.sin(ha0)*xdiff + np.cos(ha0)*ydiff
    v = -np.sin(dec0)*np.cos(ha0)*xdiff + np.sin(dec0)*np.sin(ha0)*ydiff + np.cos(dec0)*zdiff
    w = np.cos(dec0)*np.cos(ha0)*xdiff - np.cos(dec0)*np.sin(ha0)*ydiff + np.sin(dec0)*zdiff

    return u,v,w

def rotation_applied_by_woden(current_lst, mjd, latitude, ra0 = 0.0, dec0 = 0.0):

    current_time_rotation = np.zeros((3,3))
    current_time_rotation[0,0] = np.cos(current_lst)
    current_time_rotation[0,1] = -np.sin(current_lst)
    current_time_rotation[1,0] = np.sin(current_lst)
    current_time_rotation[1,1] = np.cos(current_lst)
    current_time_rotation[2,2] = 1


    # print("INPUTS INTO PALPRENUT", 2000.0, mjd)

    pal_rotation_matrix = pal.prenut(2000.0, mjd)
    pal_rotation_matrix = np.transpose(pal_rotation_matrix)

    # print(pal_rotation_matrix)

    v1 = pal.dcs2c(current_lst, latitude)
    v2 = pal.dmxv(pal_rotation_matrix, v1)
    J2000_lst, J2000_lat = pal.dcc2s(v2)

    J2000_lst = pal.dranrm(J2000_lst)


    print("Current_LST, Current_Latitude", current_lst, latitude)

    if J2000_lst > np.pi:
        print("New_LST, New_Latitude", 2*np.pi - J2000_lst, J2000_lat)
    else:
        print("New_LST, New_Latitude", J2000_lst, J2000_lat)

    return np.matmul(pal_rotation_matrix, current_time_rotation), J2000_lst, J2000_lat


def rotate_enh_to_account_for_nonzero_lst(enh, latitude):

    lat = latitude

    hen = np.array([enh[2], enh[0], enh[1]])

    account_for_latitude_rot_matrix = np.array([[np.cos(lat), 0, np.sin(lat)],
                                                [0, 1, 0],
                                                [-np.sin(lat), 0, np.cos(lat)]])

    hen = np.matmul(account_for_latitude_rot_matrix, hen)

    return account_for_latitude_rot_matrix, np.array([hen[1], hen[2], hen[0]])

def calculate_uvw_from_enh_in_J2000(enh_current, woden_rotation_matrix,
                                    current_lat, J2000_lst, J2000_lat):
    x, y, z = enh2xyz(enh_current[0], enh_current[1], enh_current[2], current_lat)

    xyz_current = np.array([x,y,z])
    xyz_J2000 = np.matmul(woden_rotation_matrix, xyz_current)

    u,v,w = get_uvw(xyz_J2000[0], xyz_J2000[1], xyz_J2000[2], J2000_lat, J2000_lst)

    return u,v,w

if __name__ == '__main__':

    ##observation date
    date = "2020-01-01T12:00:00.0"

    #Lat/long of the array at the observation date
    latitude = 0.10950738359636049*(np.pi/180.0)
    longitude = 79.638150061479*(np.pi/180.0)

    ##Setup location
    observing_location = EarthLocation(lat=latitude*u.rad, lon=longitude*u.rad, height=0.0)
    ##Setup time at that locatoin
    observing_time = Time(date, scale='utc', location=observing_location)
    ##Grab the LST
    LST = observing_time.sidereal_time('mean')
    LST_deg = LST.value*15.0

    ##This is the LST at the observing_time
    current_lst = LST_deg*(np.pi/180.0)

    ##Modified julian date
    mjd = observing_time.jd - 2400000.5

    ##This is what WODEN applies as a rotation to the x,y,z coords
    woden_rotation_matrix, J2000_lst, J2000_lat = rotation_applied_by_woden(current_lst,
                                                mjd, latitude)

    ##Get the opposite rotation to apply to our desired e,n,h
    inv_woden_rotation_matrix = np.transpose(woden_rotation_matrix)


    ##Known list of angles that have predictable sin/cos outputs
    phi_simples = [0.0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, 2*np.pi/3,
                   3*np.pi/4, 5*np.pi/6, np.pi, 7*np.pi/6, 5*np.pi/4]

    ##For known angles, calculate l,m,n coords and check they are legit
    # for ind, phi_simple in enumerate(phi_simples):
    for ind, phi_simple in enumerate(phi_simples):


        for num_mult in [1, 10, 100, 1000, 10000]:
        # for num_mult in [1000]:

            b = calc_b(phi_simple, num_mult)
            enh = np.array([b,b,b])

            # print("WANT", b)

            ##This rotation means this set of enh should give the desired
            ##x,y,z = np.array([b,b,b]) in the current epoch
            account_for_latitude_rot_matrix, enh_corr = rotate_enh_to_account_for_nonzero_lst(enh, latitude)

            ##This rotation should undo what the rotation applied by WODEN does

            undo_all_rotations = np.matmul(account_for_latitude_rot_matrix, inv_woden_rotation_matrix)

            ##Reorder for correct rotations
            hen = np.array([enh[2], enh[0], enh[1]])
            hen = np.matmul(undo_all_rotations, hen)

            ##Put back into e,n,h order
            enh = np.array([hen[1], hen[2], hen[0]])

            ##Check if we get out what we want (this is what should happen
            ##inside WODEN)
            u,v,w = calculate_uvw_from_enh_in_J2000(enh, woden_rotation_matrix,
                                            latitude, J2000_lst, J2000_lat)

            ##This prints out how accurate the final u,v,w calculations
            ##inside WODEN should be
            print("Offsets to desired u,v,w", u-b, v-b, w-b)

            ##Write the array layouts
            with open(f"array_layouts/array_layout_ang{ind:02d}_n{num_mult:05d}.txt", 'w') as outfile:
                outfile.write(f"{enh[0]:.16f} {enh[1]:.16f} {enh[2]:.16f}\n")
                outfile.write("0.0 0.0 0.0")
