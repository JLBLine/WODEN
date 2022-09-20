#!/usr/bin/env python3
from astropy.io import fits
import numpy as np
from sys import exit
from astropy.time import Time
# import matplotlib.pyplot as plt

def RTS_decode_baseline(blcode):
    """The ancient aips/miriad extended way of decoding a baseline. Takes
    the baseline code from the 'BASELINE' array of a uvfits, and returns the
    index of antennas 1 and 2 that form the baseline.

    Parameters
    ----------
    blcode : int
        Baseline code from a uvfits file encoded the two antennas

    Returns
    -------
    b1 : int
        Index of first antenna
    b2 : int
        Index of second antenna

    """
    blcode = int(blcode)

    if blcode > 65536:
        blcode -= 65536
        b2 = int(blcode % 2048)
        b1 = int((blcode - b2) / 2048)
    else:
        b2 = int(blcode % 256)
        b1 = int((blcode - b2) / 256)

    return b1,b2

# def add_bandpass_error(uvfits_name_in, uvfits_name_out, bandpass_error):
#     """Reads in the visibilities from `uvfits_name_in`, and applies the
#     frequency-dependent bandpass gain error `bandpass_error`. Results are
#     saved into `uvfits_name_out`.

#     `bandpass_error` must have the same length as the frequency axis of the
#     uvfits file. In `WODEN`, this is the 3rd axis (when 0 indexing).

#     Parameters
#     -----------
#     uvfits_name_in : string
#         Name of uvfits file to modify
#     uvfits_name_out : string
#         Name of uvfits file to write modified data to
#     bandpass_error : array
#         Frequency dependent

#     """

#     return

# def add_bandpass_error_uvfits(uvfits_prepend=False, uvfits_name=False,
#                               num_bands=24, bandpass_error):
#     """Applies the frequency-dependent bandpass gain error `bandpass_error` to
#     uvfits files as specified by either a combination of `uvfits_prepend` and
#     `num_bands`

#     `bandpass_error` must have the same length as the frequency axis of the
#     uvfits file. In `WODEN`, this is the 3rd axis (when 0 indexing).

#     Parameters
#     -----------
#     uvfits_name_in : string
#         Name of uvfits file to modify
#     uvfits_name_out : string
#         Name of uvfits file to write modified data to
#     bandpass_error : array
#         Frequency dependent
#
#     """

def make_single_polarsiation_jones_element(num_antennas, num_freqs, 
        amp_err=0.05, phase_err=10, freq_dep_jones_entry=False):
    """Make some jones matrix errors
    
    Returns a (num_antennas, num_freqs, 2, 2) shape complex array

    """

    ##If making jones_entry with a frequency dependence, do this
    if freq_dep_jones_entry:

        ##First up, make the real scalar gain error - one per antenna
        jones_entry = 1 + np.random.uniform(-amp_err, amp_err, (num_antennas, num_freqs))
        
        ##Make things complex
        jones_entry = jones_entry + 1j*np.zeros((num_antennas, num_freqs))

    ##If making a single flat gain per antenna, do this
    else:

        ##First up, make the real scalar gain error - one per antenna
        jones_entry = 1 + np.random.uniform(-amp_err, amp_err, num_antennas)
        
        ##Make things complex
        jones_entry = jones_entry + 1j*np.zeros(num_antennas)

        ##Make a frequency axis
        jones_entry = np.repeat(jones_entry, num_freqs)
        jones_entry.shape = (num_antennas, num_freqs)

    for ant in range(num_antennas):

        phase =  np.random.uniform(-phase_err*(np.pi/180.0), phase_err*(np.pi/180.0))

        phase_grad = np.linspace(-phase, phase, num_freqs)

        jones_entry[ant] = jones_entry[ant]*np.exp(1j*phase_grad)

    ##First gain is reference gain
    jones_entry[0] = 1.0 + 0.0j


    ##TODO - write these out to a `hyperdrive` style gain FITS?
    # np.savetxt("applied_jones_entry.txt", jones_entry, fmt='%.10e')

    return jones_entry


def make_antenna_jones_matrices(num_antennas, num_freqs, 
        gain_amp_err=0.05, gain_phase_err=10, freq_dep_gains=False,
        leakage_amp_err=0.05, leakage_phase_err=10, freq_dep_leakages=False,
        equal_gains=False, equal_leakages=True, zero_leakage=True):
    
    
    antenna_jones_matrices = np.zeros((num_antennas, num_freqs, 2, 2),
                                      dtype=complex)
    
    if equal_gains:
    
        gx = make_single_polarsiation_jones_element(num_antennas, num_freqs, 
                amp_err=gain_amp_err, phase_err=gain_phase_err,
                freq_dep_jones_entry=freq_dep_gains)
        gy = gx
        
    else:
        gx = make_single_polarsiation_jones_element(num_antennas, num_freqs, 
                amp_err=gain_amp_err, phase_err=gain_phase_err,
                freq_dep_jones_entry=freq_dep_gains)
        gy = make_single_polarsiation_jones_element(num_antennas, num_freqs, 
                amp_err=gain_amp_err, phase_err=gain_phase_err,
                freq_dep_jones_entry=freq_dep_gains)
        
    antenna_jones_matrices[:, :, 0, 0] = gx
    antenna_jones_matrices[:, :, 1, 1] = gy
    # antenna_jones_matrices[:, :, 0, 0] = complex(1, 0)
    # antenna_jones_matrices[:, :, 1, 1] = complex(1, 0)
    
        
    if zero_leakage:
        Dx = Dy = np.zeros(num_antennas)
    else:
        if equal_leakages:
            Dx = make_single_polarsiation_jones_element(num_antennas,
                                    num_freqs, amp_err=leakage_amp_err,
                                    phase_err=leakage_phase_err,
                                    freq_dep_jones_entry=freq_dep_gains)
            Dy = Dx
        else:
            Dx = make_single_polarsiation_jones_element(num_antennas,
                                    num_freqs, amp_err=leakage_amp_err,
                                    phase_err=leakage_phase_err,
                                    freq_dep_jones_entry=freq_dep_gains)
            Dy = make_single_polarsiation_jones_element(num_antennas,
                                    num_freqs, amp_err=leakage_amp_err,
                                    phase_err=leakage_phase_err,
                                    freq_dep_jones_entry=freq_dep_gains)
        antenna_jones_matrices[:, :, 0, 1] = Dx        
        antenna_jones_matrices[:, :, 1, 0] = Dy
        
        
    antenna_jones_matrices
    
    np.savez_compressed("gains_applied_woden.npz",
                        gx=gx, Dx=Dx, Dy=Dy, gy=gy,
                        antenna_jones_matrices=antenna_jones_matrices)
        
    
    return antenna_jones_matrices

def apply_antenna_jones_matrices(num_antennas, visibilities,
                                 antenna_jones_matrices, b1s, b2s, 
                                 frequency_dependent=False):

    # print(visibilities.shape)
    # print(len(b1s))
    # print(len(gains))

    jones_b1s = antenna_jones_matrices[b1s]
    jones_b2s = antenna_jones_matrices[b2s]
    
    print(jones_b1s.shape)

    reshape_visi = np.empty((visibilities.shape[0], visibilities.shape[1],
                             2, 2), dtype=complex)
    
    reshape_visi[:, :, 0, 0] = visibilities[:, :, 0]
    reshape_visi[:, :, 1, 1] = visibilities[:, :, 1]
    reshape_visi[:, :, 0, 1] = visibilities[:, :, 2]
    reshape_visi[:, :, 1, 0] = visibilities[:, :, 3]
    
    reshape_visi = np.matmul(np.matmul(jones_b1s, reshape_visi), np.conjugate((jones_b2s)))
    
    visibilities[:, :, 0] = reshape_visi[:, :, 0, 0]
    visibilities[:, :, 1] = reshape_visi[:, :, 1, 1]
    visibilities[:, :, 2] = reshape_visi[:, :, 0, 1]
    visibilities[:, :, 3] = reshape_visi[:, :, 1, 0]
    
    return visibilities

def visibility_noise_stddev(freq_vec, time_res, freq_res,
                          Trec=50, Aeff=20.35):
    """
    For an input set of frequencies, and time in seconds, calculate
    the visibility noise using the radiometer equation. default MWA
    values are given.

    See Equation 6.50 in TMS Third Edition for details of radiometer eq.
    Parameters
        ----------
        freq_vec : numpy array, float
            Vector of fine channels in Hz.
        time_res : float
            Observation time in seconds
        freq_res : float
            Fine channel width in Hz, default is 80kHz.
        Trec : float
            Reciever temperature in Kelvin.
        Aeff : float
            Effective MWA tile area, default is 20.35 for chan 136. Chan 136\n
            is the default middle channel for high band obs.
        Returns
        -------
            Vector of sigma values for each fine channel. These are used to generate\n
            random Gaussian noise.
    """
    # Boltzmans constant
    kb = 1380.648 #[Jy K^-1 m^2]
    freq_temp_vec = freq_vec/1e+6 # [MHz]
    # calculate sky temperature.
    Tsky_vec = 228*(freq_temp_vec/150)**(-2.53)

    # print(Tsky_vec)

    # Standard deviation term for the noise:
    sigma_vec = (np.sqrt(2)*kb*(Tsky_vec + Trec)) / (Aeff*np.sqrt(freq_res*time_res)) #[Jy]

    return sigma_vec

def add_visi_noise(visibilities, all_freqs, freq_res, time_res):
    """
    Add visibilities noise based on the radiometer equation.
    Shape of `visibilities` should be (num_visis, num_freqs, num_pols)
    
    """

    noise_stddev = visibility_noise_stddev(all_freqs, time_res, freq_res)

    print(f'First freq std dev {noise_stddev[0]:.2e}')

    num_visi = visibilities.shape[0]
    num_pols = visibilities.shape[2]

    for freq_ind, stddev in enumerate(noise_stddev):

        ##I think noise on real and imag are uncorrelated??
        ##The noise is calculated for complex values so divide by root two
        ##when calcing real or imag
        real_noise = np.random.normal(0, stddev, (num_visi, num_pols))
        imag_noise = np.random.normal(0, stddev, (num_visi, num_pols))

        # if freq_ind == 0:
        #     plt.hist(real_noise, histtype='step')
        #     plt.show()

        visibilities[:, freq_ind, :] += real_noise + 1j*imag_noise

    return visibilities

def get_parser():
    """
    Runs the argument parser to get command line inputs - used by sphinx and
    argparse extension to unpack the help below into the online readthedocs
    documentation.

    Returns
    -------
    parser : `argparse.ArgumentParser`
        The populated argument parser used by `add_instrumental_effects_woden.py`

    """
    import argparse

    parser = argparse.ArgumentParser(description="Do daa ")

    band_group = parser.add_argument_group('ADDING EFFECTS TO COARSE BAND UVFITS')
    band_group.add_argument('--uvfits_prepend', default=False,
        help='Prepend for a set of uvfits files e.g. ./data/filename_band')
    band_group.add_argument('--num_bands', default=24, type=int,
        help='How many pairs of files to combine - default is 24')
    band_group.add_argument('--output_name_prepend', default="instrumental_band",
        help='Name for start of output uvfits file, default: instrumental_band')

    sing_group = parser.add_argument_group('ADDING EFFECTS TO SINGLE UVFITS')
    sing_group.add_argument('--uvfits', default=False,
        help='Name of the uvfits file to add instrumental effects to e.g. filename.uvfits')
    sing_group.add_argument('--output_name', default="instrumental.uvfits",
        help='Name for output uvfits file, default: instrumental.uvfits')

    eff_group = parser.add_argument_group('POSSIBLE INSTRUMENTAL EFFECTS')
    eff_group.add_argument('--antenna_gain_error', default=0, type=float,
        help='Add a single multiplicative gain error per antenna, of 1 +/- '
             'the value given (e.g. gains between 0.95 and 1.05 if '
             '--antenna_gain_error=0.05. If no value given, defaults to 0.05')
    eff_group.add_argument('--antenna_phase_error', default=0, type=float,
        help='Add a phase error (degrees) per antenna, which will make a '
             'frequency dependent phase gradient between value to value, '
             'e.g if --antenna_phase_error=10, a phase error of up to '
             '-10 deg  will be added to the lowest frequency, a phase'
             ' error of up to +10 deg will be added to the highest '
             'frequency, with a smooth graident over phases added for '
             'all other frequencies.')
    eff_group.add_argument('--antenna_gain_error_freq_random',
        default=False, action='store_true',
        help='By default, a single gain error is added per antenna, with '
              'frequency dependence. Add this switch to add a random '
              'gain for all frequencies')

    eff_group.add_argument('--add_visi_noise',
        default=False, action='store_true',
        help='Add visibility noise via the radiometer equation. '
             'Defaults to MWA-like parameters for reciever '
             'temperature and effective tile area')



    return parser

# def check_args(args):
#     """Check that the args returned by """
#
#     if not args.uvfits_prepend1 and not args.uvfits_prepend2 and not args.uvfits1 and not args.uvfits2:
#        exit("You must set either a combination of --uvfits_prepend1 and "
#             "--uvfits_prepend2 or --uvfits1 and --uvfits2 otherwise I don't "
#             "know what to sum")
#
#     if args.uvfits_prepend1 and not args.uvfits_prepend2:
#         exit("If you supply --uvfits_prepend1 you must supply --uvfits_prepend2")
#     if args.uvfits_prepend2 and not args.uvfits_prepend1:
#         exit("If you supply --uvfits_prepend2 you must supply --uvfits_prepend1")
#
#     if args.uvfits1 and not args.uvfits2:
#         exit("If you supply --uvfits1 you must supply --uvfits2")
#     if args.uvfits2 and not args.uvfits1:
#         exit("If you supply --uvfits2 you must supply --uvfits1")

def add_all_errors_single_uvfits(args):
    """Add all the errors to a single uvfits file"""

    with fits.open(args.uvfits) as hdu:


        ##Leave the weights alone
        visibilities = hdu[0].data.data[:,0,0,:,:,:2]
        num_antennas = int(hdu[1].header['NAXIS2'])
        b1s, b2s = [], []

        for blcode in hdu[0].data['BASELINE']:
            b1, b2 = RTS_decode_baseline(blcode)
            b1s.append(b1)
            b2s.append(b2)

        ##BLCODE is one indexed, python is zero indexed
        b1s = np.array(b1s) - 1
        b2s = np.array(b2s) - 1

        num_freqs = hdu[0].header['NAXIS4']
        cent_freq = hdu[0].header['CRVAL4']
        ##subtract one because this is one indexed not zero
        cent_pix = hdu[0].header['CRPIX4'] - 1
        freq_res = hdu[0].header['CDELT4']

        all_freqs = cent_freq + (np.arange(num_freqs) - cent_pix)*freq_res

        ##Look to see how many antennas (tiles) there are
        num_ants = hdu[1].data['STABXYZ'].shape[0]

        ##This is total number of visibilities (for all time steps)
        num_visis = hdu[0].header['GCOUNT']

        ##Number of cross-correlations and auto-correlations
        num_cross = int((num_ants * (num_ants - 1)) / 2)
        num_autos = num_ants

        ##Work out if there are auto-correlations or not
        if num_visis % (num_cross + num_autos) == 0:
            num_visi_per_time = num_cross + num_autos
        else:
            num_visi_per_time = num_cross

        num_times = int(num_visis / num_visi_per_time)

        # print(num_visi_per_time)

        # time_res = hdu[0].data['DATE'][num_visi_per_time] - hdu[0].data['DATE'][0]

        # time_res *= (24*60*60.0)

        time1 = Time(hdu[0].data['DATE'][num_visi_per_time], format='jd')
        time0 = Time(hdu[0].data['DATE'][0], format='jd')

    time_res = time1 - time0
    time_res = time_res.to_value('s')

    print("Found the following in the `uvfits`:")
    print(f"\tNum visi per time step: {num_visi_per_time}")
    print(f"\tNum time steps: {num_times}")
    print(f"\tTime res: {time_res:.2f} s")
    print(f"\tFreq res: {freq_res:.5f} Hz")


    ##Make them complex so we can do complex maths easier
    visibilities = visibilities[:,:,:,0] + 1j*visibilities[:,:,:,1]

    print("NUM ANTENNAS",num_antennas)
    if args.antenna_gain_error or args.antenna_phase_error:
        print("Adding antenna gains... ")
        antenna_jones_matrices = make_antenna_jones_matrices(num_antennas, num_freqs,
                                   gain_amp_err=args.antenna_gain_error,
                                   gain_phase_err=args.antenna_phase_error,
                                   freq_dep_gains=args.antenna_gain_error_freq_random)

        visibilities = apply_antenna_jones_matrices(num_antennas,
                            visibilities, antenna_jones_matrices, b1s, b2s)
        print("Finished adding antenna gains.")

    if args.add_visi_noise:
        print("Adding visibility noise... ")
        visibilities = add_visi_noise(visibilities, all_freqs, 
                                      freq_res, time_res)
        print("Finished adding visibility noise")


    with fits.open(args.uvfits) as hdu:
        ##Leave the weights alone
        hdu[0].data.data[:,0,0,:,:,0] = np.real(visibilities)
        hdu[0].data.data[:,0,0,:,:,1] = np.imag(visibilities)

        hdu.writeto(args.output_name, overwrite=True)

def main():
    """Runs all the things"""
    parser = get_parser()
    args = parser.parse_args()

    # check_args(args)

    if args.uvfits:
        add_all_errors_single_uvfits(args)


if __name__ == '__main__':
    np.random.seed(87234)
    main()

    # with fits.open('error_added_test.uvfits') as hdu:
    #     data = hdu[0].data.data[:,0,0,:,:,:]

    # # print(data[:8128, 0, 0, 0])
    # # print(data[345, 0, 0, :])
    # # print(data[9873, 0, 0, :])
