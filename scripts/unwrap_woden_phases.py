import numpy as np
import os
import warnings
import sys
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation
from astropy import units as u
import matplotlib.pyplot as plt
import palpy as pal

D2R = np.pi / 180

MWA_LAT = -26.7033194444
MWA_LAT_RAD = MWA_LAT*D2R

MWA_LONG = 116.670813889
MWA_HEIGHT = 377.

VELC = 299792458.0

DS2R = 7.2722052166430399038487115353692196393452995355905e-5
SOLAR2SIDEREAL = 1.00274

def decode_baseline(baselines):
    b2 = baselines % 256
    b1 = (baselines - b2) / 256
    return b1, b2

def calc_lmn(ra,ra0,dec,dec0):
    '''Calculate l,m for a given phase centre ra0,dec0 and sky point ra,dec
    Enter angles in radians'''

    ##RTS way of doing it
    cdec0 = np.cos(dec0)
    sdec0 = np.sin(dec0)
    cdec = np.cos(dec)
    sdec = np.sin(dec)
    cdra = np.cos(ra-ra0)
    sdra = np.sin(ra-ra0)
    l = cdec*sdra
    m = sdec*cdec0 - cdec*sdec0*cdra
    n = sdec*sdec0 + cdec*cdec0*cdra

    return l,m,n

def get_parser():
    """
    Runs the argument parser to get command line inputs - used by sphinx and
    argparse extension to unpack the help below into the online readthedocs
    documentation.

    Returns
    -------
    parser : `argparse.ArgumentParser`
        The populated argument parser used by `run_woden.py`

    """
    import argparse
    from argparse import RawTextHelpFormatter

    class SmartFormatter(argparse.HelpFormatter):
        """Argparse by default ignores all \n and \t formatters. If you start
        a help class with R| the formatters will be respected."""
        def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()
            # this is the RawTextHelpFormatter._split_lines
            return argparse.HelpFormatter._split_lines(self, text, width)

    parser = argparse.ArgumentParser(description="Use this to unwrap the phase "
                "tracking applied in WODEN simulation so the uvfits can be "
                "read directly into the RTS. This script does NOT change the "
                "u,v,w coords as the RTS ignores those, so be warned.", formatter_class=SmartFormatter)

    parser.add_argument('--input_uvfits_prepend', default=False,
        help='The uvfits files to be converted should all have the same prepend '
             'and end in a band number. For example, if you have a uvfits'
             'file "data/uvfits_name_band01.uvfits" you should enter'
             '"--uvfits_prepend=data/uvfits_name_band". The script will then'
             'use --band_nums to add as many band numbers as required')
    parser.add_argument('--band_nums', default='all',
        help='Defaults to running 24 coarse bands. Alternatively, enter required'
             ' numbers delineated by commas, e.g. --band_nums=1,7,9')
    parser.add_argument('--single_uvfits', default=False,
        help='Alternatively, just give a single uvfits file to transform.'
             'Overrides --band_num and --uvfits_prepend if give.')
    parser.add_argument('--output_uvfits_prepend', default=False,
        help='What to name the outputs with.')
    parser.add_argument('--metafits', default=False,
        help='The metafits associated with the observation (if not avaible, '
             'you must provide the --lst of the observation')
    parser.add_argument('--lst', default=False,
        help='RTS needs the "phase centre" to point to zenith, so needs the '
             'LST. Can be read directly from a metafits file using --metafits '
             'or provided directly here.')

    return parser


def calc_uvw(X, Y, Z, d, H):

    u = np.sin(H)*X + np.cos(H)*Y
    v = -np.sin(d)*np.cos(H)*X + np.sin(d)*np.sin(H)*Y + np.cos(d)*Z
    w = np.cos(d)*np.cos(H)*X - np.cos(d)*np.sin(H)*Y + np.sin(d)*Z

    return u,v,w

def rotation_applied_by_rts(current_lst, mjd, latitude):

    current_time_rotation = np.zeros((3,3))
    current_time_rotation[0,0] = np.cos(current_lst)
    current_time_rotation[0,1] = -np.sin(current_lst)
    current_time_rotation[1,0] = np.sin(current_lst)
    current_time_rotation[1,1] = np.cos(current_lst)
    current_time_rotation[2,2] = 1

    pal_rotation_matrix = pal.prenut(2000.0, mjd)
    pal_rotation_matrix = np.transpose(pal_rotation_matrix)

    v1 = pal.dcs2c(current_lst, latitude)
    v2 = pal.dmxv(pal_rotation_matrix, v1)
    J2000_lst, J2000_lat = pal.dcc2s(v2)

    J2000_lst = pal.dranrm(J2000_lst)

    # print("Current_LST, Current_Latitude", current_lst, latitude)

    # if J2000_lst > np.pi:
    #     print("New_LST, New_Latitude", 2*np.pi - J2000_lst, J2000_lat)
    # else:
    #     print("New_LST, New_Latitude", J2000_lst, J2000_lat)

    return np.matmul(pal_rotation_matrix, current_time_rotation), J2000_lst, J2000_lat

def rotate_xyz_yo_J2000(Xs, Ys, Zs, rts_rotation_matrix, J2000_lst):

    xyz_current = np.array([Xs,Ys,Zs])
    xyz_J2000 = np.matmul(rts_rotation_matrix, xyz_current)


    J2000_time_rotation = np.zeros((3,3))
    J2000_time_rotation[0,0] = np.cos(J2000_lst)
    J2000_time_rotation[0,1] = np.sin(J2000_lst)
    J2000_time_rotation[1,0] = -np.sin(J2000_lst)
    J2000_time_rotation[1,1] = np.cos(J2000_lst)
    J2000_time_rotation[2,2] = 1

    X_J2000 = np.cos(J2000_lst)*xyz_J2000[0] + np.sin(J2000_lst)*xyz_J2000[1]
    Y_J2000 = -np.sin(J2000_lst)*xyz_J2000[0] + np.cos(J2000_lst)*xyz_J2000[1]
    Z_J2000 = xyz_J2000[2]

    return X_J2000, Y_J2000, Z_J2000


def remove_phase_tracking_from_v_container(frequencies=None, wws_seconds=None,
                          num_time_steps=None, v_container=None,
                          num_baselines=None, lst=None, time_res=None,
                          Xs=None, Ys=None, Zs=None,
                          ra_phase=None, dec_phase=None, mjd=None):
    """
    WARNING - currently does not change the :math:`u,v,w` coordinates, so they
    are still defined via the original phase centre. This function really is
    just to feed uvfits into the RTS (which generates it's own u,v,w using the
    antenna table)

    Undoes phase tracking applied by WODEN - to phase track, a phase was applied
    to counter the delay term caused by :math:`w` term of baseline - so just
    apply the opposite effect of the w term, i.e.

    .. math::
        V^\\prime = V \\exp(2\pi i w)

    where :math:`V` is the phase tracked visibility and :math:`V^\\prime` is
    the visibility after removing phase tracking.

    Parameters
    ----------

    frequencies : float array
        Frequencies of all fine channels (Hz)
    wws_seconds : float array
        The :math:`w` coordinates (seconds)
    num_baselines : int
        Number of baselines
    v_container : float array
        Complex visibility data out of WODEN with phase tracking, with
        `shape=(num_time_steps*num_baselines,1,1,num_freq_channels,4,3))`

    Returns
    -------
    v_container : float array
        Same visibility data as before, with phase tracking returned.

    """

    sign = 1
    PhaseConst = 2j * np.pi * sign

    num_freqs = len(frequencies)

    rts_rotation_matrix, J2000_lst, J2000_lat = rotation_applied_by_rts(lst, mjd, MWA_LAT_RAD)

    ##Calculate the u,v,w at zenith in the J2000 frame
    u_metres, v_metres, w_metres = calc_uvw(Xs, Ys, Zs, J2000_lat, 0.0000000000)

    for time_ind in np.arange(num_time_steps):

        # these_wws_secs = wws_seconds[time_ind*num_baselines:(time_ind + 1)*num_baselines]

        lst = J2000_lst
        prec_lat = J2000_lat

        ra_zenith = J2000_lst + (time_ind + 0.5) * time_res *SOLAR2SIDEREAL*DS2R;

        l,m,n = calc_lmn(ra_phase, ra_zenith, dec_phase, J2000_lat)


        ha_phase = ra_zenith - ra_phase

        for freq_ind, freq in enumerate(frequencies):

            xx_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0]
            xx_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1]

            yy_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0]
            yy_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1]

            xy_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0]
            xy_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1]

            yx_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0]
            yx_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1]

            xx_comp = xx_re + 1j*xx_im
            yy_comp = yy_re + 1j*yy_im
            xy_comp = xy_re + 1j*xy_im
            yx_comp = yx_re + 1j*yx_im

            ##theory - so normal phase delay is caused by path difference across
            ##a base line, which is u*l + v*m + w*n
            ##To phase track, you insert a phase to make sure there is no w contribution at
            ##phase centre; this is when n = 1, so you insert a phase thus:
            ##a base line, which is u*l + v*m + w*(n - 1)
            ##So we just need to remove the effect of the -w term

            # if time_ind == 0 and freq_ind == 0:
            #     print("xx_comp[:10]", xx_comp[:10])

            # wws = these_wws_secs * (freq)
            # phase_rotate = np.exp( PhaseConst * wws)
            # print(wws[:5])

            wavelength = VELC / freq
            u_lamb = u_metres / wavelength
            v_lamb = v_metres / wavelength
            w_lamb = w_metres / wavelength

            phase_rotate = np.exp( PhaseConst * (u_lamb*l + v_lamb*m + w_lamb*n))

            xx_comp = xx_comp * phase_rotate
            yy_comp = yy_comp * phase_rotate
            xy_comp = xy_comp * phase_rotate
            yx_comp = yx_comp * phase_rotate

            # if time_ind == 0 and freq_ind == 0:
            #     print("xx_comp[:10]", xx_comp[:10])

            # v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0] = np.real(xx_comp)
            # v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1] = np.imag(xx_comp)
            # v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0] = np.real(yy_comp)
            # v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1] = np.imag(yy_comp)
            # v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0] = np.real(xy_comp)
            # v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1] = np.imag(xy_comp)
            # v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0] = np.real(yx_comp)
            # v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1] = np.imag(yx_comp)

            ##Swap the polarisation ordering because fuck you that's why
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0] = np.real(yy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1] = np.imag(yy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0] = np.real(xx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1] = np.imag(xx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0] = np.real(yx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1] = np.imag(yx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0] = np.real(xy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1] = np.imag(xy_comp)

    return v_container


def unwrap_phases_from_uvfits(filename, band_num, args):

    with fits.open(filename) as hdu:



        ##Setup location
        observing_location = EarthLocation(lat=MWA_LAT*u.deg, lon=MWA_LONG*u.deg,
                                           height=MWA_HEIGHT)
        ##Setup time at that locatoin

        # date = hdu[1].header['RDATE']

        v_container = hdu[0].data.data

        num_freq_channels = hdu[0].header['NAXIS4']
        freq_res = hdu[0].header['CDELT4']
        freq_cent_pix = hdu[0].header['CRPIX4'] - 1


        cent_freq = hdu[0].header['CRVAL4']

        ##FUCKING AROUND HERE NOW
        hdu[0].header['CRVAL4'] = cent_freq + freq_res / 2

        frequencies = cent_freq + np.arange(-freq_cent_pix, num_freq_channels-freq_cent_pix)*freq_res

        num_antennas = hdu[1].header['NAXIS2']

        ##ARGH this assumes no auto-correlations
        num_baselines = int((num_antennas*(num_antennas - 1)) / 2)

        num_time_steps = int(hdu[0].header['GCOUNT'] / num_baselines)

        print(f"Found {int(num_baselines)} baselines, {int(num_time_steps)} time steps")

        wws_seconds = hdu[0].data['WW']

        intial_jd_date = hdu[0].data['DATE'][0]
        observing_time = Time(intial_jd_date, format='jd', location=observing_location)

        # time_res = hdu[0].data['DATE'][num_baselines] - hdu[0].data['DATE'][0]
        # print(2 / (time_res*24*60*60))

        # print(intial_jd_date)
        time_res = 2.0

        # print(date)

        # ##Grab the LST
        LST = observing_time.sidereal_time('apparent')
        LST_deg = LST.value*15.0
        lst_rad = LST_deg*D2R

        # lst_rad = 0.0077478248

        print(lst_rad)


        b1s, b2s = decode_baseline(hdu[0].data['BASELINE'])

        b1s = b1s[:num_baselines]
        b2s = b2s[:num_baselines]

        all_Xs = []
        all_Ys = []
        all_Zs = []

        rts_Xs, rts_Ys, rts_Zs = np.loadtxt("rts_XYZ.txt", delimiter=',', unpack=True)

        X_locs = hdu[1].data['STABXYZ'][:, 0]
        Y_locs = hdu[1].data['STABXYZ'][:, 1]
        Z_locs = hdu[1].data['STABXYZ'][:, 2]

        rts_rotation_matrix, J2000_lst, J2000_lat = rotation_applied_by_rts(lst_rad, observing_time.mjd, MWA_LAT_RAD)


        # J2000_lst = 0.00425996

        # print("before", X_locs[0])

        X_locs, Y_locs, Z_locs = rotate_xyz_yo_J2000(X_locs, Y_locs, Z_locs,
                                                    rts_rotation_matrix, J2000_lst)

        # print("after", X_locs[0])


        # print(rts_Xs[:5])
        # print(X_locs[:5])

        # X_locs = rts_Xs
        # Y_locs = rts_Ys
        # Z_locs = rts_Zs

        for b1, b2 in zip(b1s, b2s):
            ant1 = int(b1 - 1)
            ant2 = int(b2 - 1)

            all_Xs.append(X_locs[ant1] - X_locs[ant2])
            all_Ys.append(Y_locs[ant1] - Y_locs[ant2])
            all_Zs.append(Z_locs[ant1] - Z_locs[ant2])

        all_Xs = np.array(all_Xs)
        all_Ys = np.array(all_Ys)
        all_Zs = np.array(all_Zs)

        ra_phase = hdu[0].header['CRVAL5']*D2R
        dec_phase = hdu[0].header['CRVAL6']*D2R

        v_container = remove_phase_tracking_from_v_container(frequencies=frequencies,
                                  wws_seconds=wws_seconds,
                                  num_time_steps=int(num_time_steps),
                                  v_container=v_container,
                                  num_baselines=int(num_baselines),
                                  lst=lst_rad, time_res=time_res,
                                  Xs = all_Xs, Ys = all_Ys, Zs = all_Zs,
                                  ra_phase=ra_phase, dec_phase=dec_phase,
                                  mjd=observing_time.mjd)

        # hdu[0].header['CRVAL5'] = 0.0077478248 / D2R #LST_deg
        # # # hdu[0].header['CRVAL6'] = -0.4675872326 / D2R #MWA_LAT
        # hdu[0].header['CRVAL6'] = -26.703319

        hdu[0].header['CRVAL5'] = LST_deg
        hdu[0].header['CRVAL6'] = MWA_LAT

        hdu.writeto(f"{args.output_uvfits_prepend}{band_num:02d}.uvfits", overwrite=True)


def main(args):

    uvfits_filenames = []


    if args.single_uvfits:
        uvfits_filenames.append(args.single_uvfits)

    else:

        if args.band_nums == 'all':
            args.band_nums = range(1,25)
        else:
            try:
                args.band_nums = list(np.array(args.band_nums.split(','),dtype=int))
            except:
                message = ("ERROR - failed to convert --band_nums into a list of ints"
                           " correctly. You entered:\n"
                           "    --band_nums={:s}\n"
                           "Exiting now.")
                exit(message)

        for band_num in args.band_nums:
            uvfits_filenames.append(f"{args.input_uvfits_prepend}{band_num:02d}.uvfits")

    for band_num, filename in zip(args.band_nums, uvfits_filenames):

        if os.path.exists(filename):

            unwrap_phases_from_uvfits(filename, band_num, args)
        else:
            print(f"{filename} does not exist, skipping conversion")

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    main(args)
