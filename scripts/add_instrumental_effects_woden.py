#!/usr/bin/env python3
from astropy.io import fits
import numpy as np
from sys import exit
from astropy.time import Time
from wodenpy.uvfits.wodenpy_uvfits import RTS_decode_baseline
from argparse import Namespace
from astropy.constants import c as speed_of_light
from astropy.constants import k_B
import os
from typing import Tuple
import importlib_resources
import wodenpy
from scipy.interpolate import CubicSpline

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

    parser = argparse.ArgumentParser(description="Adds instrumental effects to a WODEN uvfits, given command line inputs")

    sing_group = parser.add_argument_group('INPUT/OUTPUT OPTIONS')
    sing_group.add_argument('--uvfits', default=False,
        help='Name of the uvfits file to add instrumental effects to e.g. filename.uvfits')
    sing_group.add_argument('--output_name', default="instrumental.uvfits",
        help='Name for output uvfits file, default: instrumental.uvfits')
    sing_group.add_argument('--numpy_seed', default=0, type=int,
        help='A specific np.random.seed to use for reproducibility. Otherwise numpy is left to seed itself.')
    
    noise_group = parser.add_argument_group('NOISE EFFECTS')
    noise_group.add_argument('--add_visi_noise',
        default=False, action='store_true',
        help='Add visibility noise via the radiometer equation. '
             'Defaults to MWA-like parameters for reciever '
             'temperature and effective tile area')
    
    noise_group.add_argument('--visi_noise_int_time',
        default=False, type=float,
        help='Use a different integration time (seconds) to what is actually in the '
             'data to calculate the noise, e.g. even if you data is at 2s '
             'resolution, --visi_noise_int_time=60 will add the noise for a '
             'one minute integration.')
    noise_group.add_argument('--visi_noise_freq_reso',
        default=False, type=float,
        help='Use a different frequency channel width (Hz) to what is actually in the '
             'data to calculate the noise, e.g. even if you data is at 40kHz '
             'resolution, --visi_noise_freq_reso=1e+6 will add the noise for a '
             'channel width of 1MHz instead.')
    
    reflec_group = parser.add_argument_group('CABLE REFLECTION GROUP')
    reflec_group.add_argument('--cable_reflection_from_metafits', default=False,
        help='Given a metafits file with cable length information, add in '
             'cable reflections to the antenna gains. This will use the '
             '--cable_reflection_coeff_amp to set the magnitude.')
    
    reflec_group.add_argument('--cable_reflection_coeff_amp', default=0.05,
        type=float, help='Amplitude of the cable will be drawn from a uniform '
            'distribution between 0 and this value. Default is 0.05.')

    gain_group = parser.add_argument_group('ANTENNA (tile) GAIN EFFECTS')
    gain_group.add_argument('--ant_gain_amp_error', default=0, type=float,
        help='Add a single multiplicative gain error per antenna, of 1 +/- '
             'the value given (e.g. gains between 0.95 and 1.05 if '
             '--ant_gain_amp_error=0.05.')
    gain_group.add_argument('--ant_gain_phase_error', default=0, type=float,
        help='Add a phase error (degrees) per antenna, which will make a '
             'frequency dependent phase gradient between value to value, '
             'e.g if --ant_gain_phase_error=10, a phase error of up to '
             '-10 deg  will be added to the lowest frequency, a phase'
             ' error of up to +10 deg will be added to the highest '
             'frequency, with a smooth graident over phases added for '
             'all other frequencies.')
    gain_group.add_argument('--ant_leak_errs', default=[0,0], type=float, nargs = '*',
        help='Use as  `--ant_leak_errs psi_err chi_err` (degrees). Adds an engineering '
             'tolerance error for linear dipole alignment (see TMS '
             'eqn A4.5) to add leakage terms to the antenna Jones matrix. '
             'This leakage is based off two angles where `Dx = psi_err - 1j*chi_err` '
             'and `Dy = -psi_err + 1j*chi_err`. '
             'A random angle between 0 and the given value will be added for each angle. ')




    return parser

class UVFITS(object):
    """
    A class for reading and storing data from UVFITS files.

    Parameters:
    -----------
    filename : str
        The path to the UVFITS file to be read.

    Attributes:
    -----------
    visibilities : ndarray
        The complex visibilities of the UVFITS file.
    num_antennas : int
        The number of antennas in the UVFITS file.
    b1s : ndarray
        The first antenna indices of the baselines in the UVFITS file.
    b2s : ndarray
        The second antenna indices of the baselines in the UVFITS file.
    num_freqs : int
        The number of frequency channels in the UVFITS file.
    cent_freq : float
        The central frequency of the UVFITS file.
    cent_pix : float
        The central pixel of the UVFITS file.
    freq_res : float
        The frequency resolution of the UVFITS file.
    all_freqs : ndarray
        The frequencies of all the channels in the UVFITS file.
    num_ants : int
        The number of antennas (tiles) in the UVFITS file.
    num_visis : int
        The total number of visibilities (for all time steps) in the UVFITS file.
    num_cross : int
        The number of cross-correlations in the UVFITS file.
    num_autos : int
        The number of auto-correlations in the UVFITS file.
    has_autos : bool
        Whether auto-correlations are present in the UVFITS file.
    time_res : float
        The time resolution of the UVFITS file.
    """
class UVFITS(object):
    def __init__(self, filename):
        with fits.open(filename) as hdu:
            
            ##Leave the weights alone
            visibilities = hdu[0].data.data[:,0,0,:,:,:2]
            self.visibilities = visibilities[:,:,:,0] + 1j*visibilities[:,:,:,1]
            
            self.num_antennas = int(hdu[1].header['NAXIS2'])
            b1s, b2s = [], []

            for blcode in hdu[0].data['BASELINE']:
                b1, b2 = RTS_decode_baseline(blcode)
                b1s.append(b1)
                b2s.append(b2)

            ##BLCODE is one indexed, python is zero indexed
            b1s = np.array(b1s) - 1
            b2s = np.array(b2s) - 1
            
            self.b1s = b1s
            self.b2s = b2s

            self.num_freqs = hdu[0].header['NAXIS4']
            self.cent_freq = hdu[0].header['CRVAL4']
            ##subtract one because this is one indexed not zero
            self.cent_pix = hdu[0].header['CRPIX4'] - 1
            self.freq_res = hdu[0].header['CDELT4']

            self.all_freqs = self.cent_freq + (np.arange(self.num_freqs) - self.cent_pix)*self.freq_res

            ##Look to see how many antennas (tiles) there are
            self.num_ants = hdu[1].data['STABXYZ'].shape[0]

            ##This is total number of visibilities (for all time steps)
            self.num_visis = hdu[0].header['GCOUNT']

            ##Number of cross-correlations and auto-correlations
            self.num_cross = int((self.num_ants * (self.num_ants - 1)) / 2)
            self.num_autos = self.num_ants

            ##Work out if there are auto-correlations or not
            if self.num_visis % (self.num_cross + self.num_autos) == 0:
                num_visi_per_time = self.num_cross + self.num_autos
                self.has_autos = True
            else:
                num_visi_per_time = self.num_cross
                self.has_autos = False

            num_times = int(self.num_visis / num_visi_per_time)

            time1 = Time(hdu[0].data['DATE'][num_visi_per_time], format='jd')
            time0 = Time(hdu[0].data['DATE'][0], format='jd')
            
            self.uus = hdu[0].data['UU']

        self.time_res = time1 - time0
        self.time_res = self.time_res.to_value('s')

        print(f"Found the following in the {filename}:")
        if self.has_autos:
            print(f"\tAutocorrelations are present")
        else:
            print(f"\tAutocorrelations are not present")
        print(f"\tNum visi per time step: {num_visi_per_time}")
        print(f"\tNum time steps: {num_times}")
        print(f"\tTime res: {self.time_res:.2f} s")
        print(f"\tFreq res: {self.freq_res:.5f} Hz")
        
        

def make_single_polarsiation_jones_gain(num_antennas, num_freqs, 
        amp_err=0.05, phase_err=10):
    """
    Generate a Jones gain matrix for a single polarisation, with random amplitude and phase errors.

    Parameters
    ----------
    num_antennas : int
        Number of antennas in the array.
    num_freqs : int
        Number of frequency channels.
    amp_err : float, optional
        Maximum amplitude error (as a fraction of the true gain) to apply to each antenna. Default is 0.05.
    phase_err : float, optional
        Maximum phase error (in degrees) to apply to each antenna. Default is 10.

    Returns
    -------
    jones_entry : ndarray
        Complex Jones gain matrix of shape (num_antennas, num_freqs), with the first row set to (1+0j).
    """

    ##First up, make the real scalar gain error - one per antenna
    jones_entry = 1.0 + np.random.uniform(-amp_err, amp_err, num_antennas)
    
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
    jones_entry[0, :] = 1.0 + 0.0j


    ##TODO - write these out to a `hyperdrive` style gain FITS?

    return jones_entry

def make_jones_leakage(num_antennas, num_freqs, leak_psi_err=0.0, leak_chi_err=0.0):
    """
    Generate Jones leakage terms for a set of antennas and frequencies.

    Parameters
    -----------
    num_antennas : int
        Number of antennas.
    num_freqs : int
        Number of frequencies.
    leak_psi_err : float, optional
        Small engineering error on the linear polarization alignment of the dipoles, represented by psi_err.
    leak_chi_err : float, optional
        Small engineering error on the linear polarization alignment of the dipoles, represented by chi_err.

    Returns
    --------
    Dx : numpy.ndarray
        Jones leakage term for the x-polarization.
    Dy : numpy.ndarray
        Jones leakage term for the y-polarization.
    """
    Dx = np.repeat(np.random.uniform(0, leak_psi_err, num_antennas), num_freqs) - 1j*np.repeat(np.random.uniform(0, leak_chi_err, num_antennas), num_freqs)
    Dx.shape = (num_antennas, num_freqs)
    
    Dy = -np.repeat(np.random.uniform(0, leak_psi_err, num_antennas), num_freqs) + 1j*np.repeat(np.random.uniform(0, leak_chi_err, num_antennas), num_freqs)
    Dy.shape = (num_antennas, num_freqs)
    
    return Dx, Dy

def make_antenna_jones_matrices(num_antennas : int, num_freqs : int, 
        gain_amp_err=0.0, gain_phase_err=0.0,
        leak_psi_err=0.0, leak_chi_err=0.0):
    """Given input error parameters, create a set of instrumental Jones matrices for each antenna and frequency.

    Parameters
    ----------
    num_antennas : int
        Number of antennas.
    num_freqs : int
        Number of frequencies.
    gain_amp_err : float, optional
        Maximum amplitude error (as a fraction of the true gain) to apply to each antenna, by default 0.0
    gain_phase_err : float, optional
        Maximum phase error (in degrees) to apply to each antenna.
    leak_psi_err : float, optional
        Small engineering error on the linear polarization alignment of the dipoles, represented by psi_err, by default 0.0
    leak_chi_err : float, optional
        Small engineering error on the linear polarization alignment of the dipoles, represented by chi_err., by default 0.0

    Returns
    -------
    antenna_jones_matrices : np.ndarray
        Complex antenna jones matrices for all antennas, of shape (num_antennas, num_freqs, 2, 2)
    """
    
    antenna_jones_matrices = np.zeros((num_antennas, num_freqs, 2, 2),
                                      dtype=complex)
    
    gx = make_single_polarsiation_jones_gain(num_antennas, num_freqs, 
            amp_err=gain_amp_err, phase_err=gain_phase_err)
    gy = make_single_polarsiation_jones_gain(num_antennas, num_freqs, 
            amp_err=gain_amp_err, phase_err=gain_phase_err)
    
    Dx, Dy = make_jones_leakage(num_antennas, num_freqs,
                                leak_psi_err=leak_psi_err,
                                leak_chi_err=leak_chi_err)
    
    antenna_jones_matrices[:, :, 0, 0] = gx
    antenna_jones_matrices[:, :, 1, 1] = gy
    antenna_jones_matrices[:, :, 0, 1] = Dx        
    antenna_jones_matrices[:, :, 1, 0] = Dy
    
    np.savez_compressed("gains_applied_woden.npz",
                        gx=gx, Dx=Dx, Dy=Dy, gy=gy,
                        antenna_jones_matrices=antenna_jones_matrices)
        
    
    return antenna_jones_matrices

def apply_antenna_jones_matrices(visibilities, antenna_jones_matrices,
                                 b1s, b2s):
    """
    Apply the antenna Jones matrices to the visibilities for a set of antenna
    pairs defined by b1s, b2s

    Parameters
    ----------
    visibilities : numpy.ndarray
        The visibilities to which the Jones matrices will be applied. The
        shape of the array should be (Nvis, Nfreqs, 4), where Nvis is the
        number of visibilities, Nfreqs is the number of frequency channels, and
        the last dimension contains the complex visibilities for all 4 instrumental
        stokes polarisations.
    antenna_jones_matrices : numpy.ndarray
        The Jones matrices for each antenna in the array. The shape of the array
        should be (Nants, Nfreqs, 2, 2), where Nants is the number of antennas and the
        last two dimensions contain the complex Jones matrix for each antenna.
    b1s : numpy.ndarray
        The indices of the first antenna in each baseline.
    b2s : numpy.ndarray
        The indices of the second antenna in each baseline.

    Returns
    -------
    numpy.ndarray
        The visibilities with the Jones matrices applied. The shape of the array
        is the same as the input `visibilities`.

    """

    jones_b1s = antenna_jones_matrices[b1s]
    jones_b2s = antenna_jones_matrices[b2s]
    
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

def get_cable_delay(length, velocity_factor=0.81):
    """
    Calculate the delay imparted by a cable of specific length and velocity factor.

    Parameters
    ----------
    length : float
        Length of the cable in meters.
    velocity_factor : float
        Velocity factor of the cable (0.81 from Beardsley et al. 2016)

    Returns
    -------
    float
        The delay of the cable in seconds.
    """
    return (2*length) / (speed_of_light.value * velocity_factor)

def create_single_pol_reflections(freqs : np.ndarray,
                                  cable_reflection_coeff_amp : float,
                                  delays : np.ndarray) -> np.ndarray:
    """Calculate the reflections caused by a set of cables with a
    given delays, assigns a random reflection coefficient amplitude
    for the given frequency vector. Returns a 2D array of shape (len(delays), len(freqs)).

    Parameters
    ----------
    freqs : np.ndarray
        Frequencies to calculate the reflections at (Hz)
    cable_reflection_coeff_amp : float
        Maximum amplitude of the cable reflection coefficient.
    delays : np.ndarray
        Delays of the cable reflections (seconds)

    Returns
    -------
    reflections : np.ndarray
        A 2D array of reflections of shape (len(delays), len(freqs)) for
        all delays and frequencies.
    """
    
    reflection_coeffs = np.random.uniform(0, cable_reflection_coeff_amp, len(delays)) + 1j*np.random.uniform(0, cable_reflection_coeff_amp, len(delays))
    
    ftimesd = np.outer(delays, freqs)
    reflections = reflection_coeffs[:, np.newaxis]*np.exp(-2j * np.pi * ftimesd)
    
    return reflections

def create_cable_reflections(args : Namespace, uvfits : UVFITS) -> Tuple[np.ndarray, np.ndarray]:
    """Given the input arguments and a UVFITS object, create the cable reflections for the x and y pols.

    Parameters
    ----------
    args : Namespace
        Namespace object containing command line arguments.
    uvfits : UVFITS
        UVFITS object containing visibilities to which reflections will be added.

    Returns
    -------
    reflections_x, reflections_y : np.ndarray, np.ndarray
        2D arrays of reflections to add to the x and y pols, of shape (num_antennas, num_freqs).
    """
    
    
    if not os.path.isfile(args.cable_reflection_from_metafits):
            exit('Could not open metafits specified by user as:\n'
                 '\t--metafits_filename={:s}.\n'
                 'Cannot get required observation settings, exiting now'.format(args.metafits_filename))

    with fits.open(args.cable_reflection_from_metafits) as f:
        antenna_order = np.argsort(f[1].data['Tile'])
        
        selection = np.arange(0,len(antenna_order),2)
        ##this has both XX and YY pol, assume they have the same cable lengths
        antenna_order = antenna_order[selection]
        
        ##the length causing reflections is in the flavour column
        flavours = f[1].data['Flavors'][antenna_order]
        cable_lengths = np.array([float(flavour.split('_')[1]) for flavour in flavours])
        
    delays = get_cable_delay(cable_lengths)
    
    reflections_x = create_single_pol_reflections(uvfits.all_freqs,
                                  args.cable_reflection_coeff_amp, delays)
    
    reflections_y = create_single_pol_reflections(uvfits.all_freqs,
                                  args.cable_reflection_coeff_amp, delays)
    
    return reflections_x, reflections_y

def add_complex_ant_gains(args : Namespace, uvfits : UVFITS):
    """
    Add complex antenna gains to the UVFITS object.

    Parameters
    ------------
    args : Namespace
        Namespace object containing command line arguments.
    uvfits : UVFITS
        UVFITS object containing visibilities to which antenna gains will be added.

    Returns
    ---------
    UVFITS:
        UVFITS object with antenna gains added.
    """
    
    antenna_jones_matrices = make_antenna_jones_matrices(uvfits.num_antennas,
                                   uvfits.num_freqs,
                                   gain_amp_err=args.ant_gain_amp_error,
                                   gain_phase_err=args.ant_gain_phase_error,
                                   leak_psi_err=args.ant_leak_errs[0]*np.pi/180.0,
                                   leak_chi_err=args.ant_leak_errs[1]*np.pi/180.0)
    
    if args.cable_reflection_from_metafits:
        reflections_x, reflections_y = create_cable_reflections(args, uvfits)
        
        antenna_jones_matrices[:, :, 0, 0] += reflections_x
        antenna_jones_matrices[:, :, 1, 1] += reflections_y

    uvfits.visibilities = apply_antenna_jones_matrices(uvfits.visibilities,
                                                       antenna_jones_matrices,
                                                       uvfits.b1s, uvfits.b2s)
    
    return uvfits

def visibility_noise_stddev(freq_vec, time_res, freq_res,
                          Trec=50, Aeff=20.35):
    """
    For an input set of frequencies, and time in seconds, calculate
    the visibility noise using the radiometer equation. default MWA
    values are given.

    See Equation 6.50 in TMS Third Edition for details of radiometer eq.
    
    Parameters
    -----------
    freq_vec : numpy array, float
        Vector of fine channels in Hz.
    time_res : float
        Observation time in seconds
    freq_res : float
        Fine channel width in Hz, default is 80kHz.
    Trec : float
        Reciever temperature in Kelvin.
    Aeff : float
        Effective MWA tile area, default is 20.35 for chan 136. Chan 136
        is the default middle channel for high band obs.
    Returns
    -------
    sigma : numpy array, float
        $\sigma$ values for each given frequency
    """
    # Boltzmans constant
    kb = k_B.value*1e+26 #[Jy K^-1 m^2]
    freq_temp_vec = freq_vec/1e+6 # [MHz]
    # calculate sky temperature.
    Tsky_vec = 230*(freq_temp_vec/150)**(-2.53)
    
    # Tsky_vec = np.ones(len(freq_vec))*230.0

    # Standard deviation term for the noise:
    sigma = (np.sqrt(2)*kb*(Tsky_vec + Trec)) / (Aeff*np.sqrt(freq_res*time_res)) #[Jy]
    
    return sigma

def add_visi_noise(args : Namespace, uvfits : UVFITS):
    """
    Adds instrumental thermal noise to the visibilities of a UVFITS object,
    given the input `args`. The noise is calculated using the radiometer equation.
    Everything is currently hardcoded to MWA values.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments parsed by argparse.
    uvfits : UVFITS
        UVFITS object to which noise will be added.

    Returns
    -------
    UVFITS
        The input UVFITS object with noise added to its visibilities.
    """
    if args.visi_noise_int_time:
        print(f"--visi_noise_int_time was set to {args.visi_noise_int_time:.1e} "
                  "seconds, using instead of time resolution inside uvfits")
            
        time_res = args.visi_noise_int_time
    else:
        time_res = uvfits.time_res
        
    if args.visi_noise_freq_reso:
        print(f"--visi_noise_freq_reso was set to {args.visi_noise_freq_reso:.3e} "
                "Hz, using instead of freq resolution inside uvfits")
        
        freq_res = args.visi_noise_freq_reso
    else:
        freq_res = uvfits.freq_res
        
    noise_stddev = visibility_noise_stddev(uvfits.all_freqs, time_res, freq_res)

    print(f'First freq std dev {noise_stddev[0]:.2e}')

    baseline_cross = np.where(uvfits.b1s != uvfits.b2s)[0]
    baseline_auto = np.where(uvfits.b1s == uvfits.b2s)[0]

    for freq_ind, stddev in enumerate(noise_stddev):
        
        ##Add noise to all the XX and YY cross-correlations
        for pol_ind in range(4):

            ##Add noise to the cross-correlations
            real_noise = np.random.normal(0, stddev, len(baseline_cross))
            imag_noise = np.random.normal(0, stddev, len(baseline_cross))
            uvfits.visibilities[baseline_cross, freq_ind, pol_ind] += real_noise + 1j*imag_noise
            
            ##Add noise to the auto-correlations
            real_noise = np.random.normal(0, np.sqrt(2)*stddev, len(baseline_auto))
            imag_noise = np.random.normal(0, np.sqrt(2)*stddev, len(baseline_auto))
            uvfits.visibilities[baseline_auto, freq_ind, pol_ind] += real_noise + 1j*imag_noise
    
    return uvfits

def calculate_bandpass(freq_res : float, num_bands : int = 24):
    
    bandpass = importlib_resources.files(wodenpy).joinpath("bandpass_1kHz.txt")
    bandpass = np.loadtxt(bandpass)
    
    cs = CubicSpline(np.arange(5e+2, 1.28e+6, 1e+3), bandpass)
    
    single_bandpass = cs(np.arange(freq_res/2, 1.28e+6, freq_res))
    
    full_bandpass = np.tile(single_bandpass, num_bands)
    
    return single_bandpass, full_bandpass


def add_bandpass(args : Namespace, uvfits : UVFITS):
    
    print(uvfits.freq_res)
    
    single_bandpass, full_bandpass = calculate_bandpass(uvfits.freq_res)
    
    return


def main(argv=None):
    """Adds instrumental effects to a WODEN uvfits, given command line inputs

    Parameters
    ----------
    argv : _type_, optional
        Will be parsed from the command line if run in main script. Can be
        passed in explicitly for easy testing
    """
    parser = get_parser()
    args = parser.parse_args(argv)
    
    if args.numpy_seed:
        np.random.seed(args.numpy_seed)

    uvfits = UVFITS(args.uvfits)
    
    if args.add_visi_noise or args.visi_noise_int_time or args.visi_noise_freq_reso:
        print("Adding visibility noise... ")
        uvfits = add_visi_noise(args, uvfits)
        print("Finished adding visibility noise")
    
    if args.ant_gain_amp_error or args.ant_gain_phase_error or args.ant_leak_errs != [0,0] or args.cable_reflection_from_metafits:
        print("Adding antenna gains... ")
        uvfits = add_complex_ant_gains(args, uvfits)
        print("Finished adding antenna gains.")
        
    with fits.open(args.uvfits) as hdu:
        ##Leave the weights alone
        hdu[0].data.data[:,0,0,:,:,0] = np.real(uvfits.visibilities)
        hdu[0].data.data[:,0,0,:,:,1] = np.imag(uvfits.visibilities)

        hdu.writeto(args.output_name, overwrite=True)


if __name__ == '__main__':
    ##Here we goooo
    main()