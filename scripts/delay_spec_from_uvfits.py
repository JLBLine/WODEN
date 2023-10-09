import numpy as np
from astropy.io import fits
import matplotlib
##useful when using a super cluster to specify Agg
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import blackmanharris
from astropy.cosmology import LambdaCDM

##Speed of light in m/s
SPEED_LIGHT = 299792458.0

##21cm wavelength in m
WAVELENGTH_21CM = 0.21

##Boltzmann constant
BOLTZMANN = 1.38064852e-23

def add_colourbar(fig=None,ax=None,im=None,label=False,top=False):
    """
    Adds a colourbar (colorbar, fine) in a nice way to a subplot

    Parameters
    ----------
    fig : matplotlib.pyplot.figure instance
        The figure that the plot lives on
    ax : matplotlib.pyplot.figure.add_subplot instance
        The axis to append a colorbar to
    im : ax.imshow output
        The output of imshow to base the colourbar on
    label : string
        Optional - add a label to the colorbar
    top : Bool
        Optional - put the colorbar above the axis instead of to the right
    """

    divider = make_axes_locatable(ax)
    if top == True:
        cax = divider.append_axes("top", size="5%", pad=0.05)
        cbar = fig.colorbar(im, cax = cax,orientation='horizontal')
        cax.xaxis.set_ticks_position('top')
        cax.xaxis.set_label_position('top')
    else:
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, cax = cax)
    if label:
        cbar.set_label(label)


def calc_conversion_factor(freqs, args):
    """Calculate the conversion from Jy**2 to P(k)"""

    central_freq = freqs[int(len(freqs)/2)]

    ##21cm radiation frequency in m/s
    f21 = SPEED_LIGHT / WAVELENGTH_21CM

    ##Redshift based on central frequency
    z = f21 / central_freq - 1.

    ##Create the cosmology and report the proprties
    cosmology = LambdaCDM(H0=args.hubble,
                          Om0=args.omega_matter,
                          Ode0=args.omega_lambda,
                          Ob0=args.omega_baryon)

    nu_21 = (SPEED_LIGHT)/(WAVELENGTH_21CM) #[Hz]

    # Cosmological scaling parameter:
    h = cosmology.H(0).value/100 # Hubble parameter.

    E_z = cosmology.efunc(z) ## Scaling function, see (Hogg 2000)

    # Cosmological distances:
    Dm = cosmology.comoving_distance(z).value*h #[Mpc/h] Transverse co-moving distance.
    DH = 3000 # [Mpc/h] Hubble distance.
    # ##New way of doing it===================================================
    cent_wavelength = SPEED_LIGHT / central_freq
    beam_area_steradian = 0.07597

    ##Frequency bandwidth of the data
    bandwidth = freqs[-1] - freqs[0]

    norm_numer = cent_wavelength**4 * Dm**2*DH * bandwidth * (1 + z)**2 * 1.e6
    norm_denom = len(freqs)*(2*BOLTZMANN*1e26)**2*beam_area_steradian*f21*E_z

    # self.normalisation = norm_numer / norm_denom

    return norm_numer / norm_denom


def calc_delay_spec_multi_baselines(visi, freqs, args):
    """The actual delay spec calculation - multiple visibilities by
    a blackman harris, FTT, multiply by complex conj and convert into P(b)"""

    spec_taper = blackmanharris(len(freqs))

    spec_taper_visi = visi*spec_taper[np.newaxis, :]

    ft_visi = np.fft.fftshift(np.fft.fft(spec_taper_visi), axes=1)

    conversion = calc_conversion_factor(freqs, args)

    return np.real(ft_visi*np.conj(ft_visi))*conversion


def do_delay_spec(bin_edges, bin_indexes, all_visi_data_by_bin, freqs, args):
    """Calculate the delay spectrum for every baseline, and average within
    each bin"""

    freq_res = freqs[1] - freqs[0]

    delays = np.fft.fftshift(np.fft.fftfreq(len(freqs), freq_res))
    all_delay_spec_by_bin = np.empty((len(delays), len(bin_indexes)))

    for bin_ind, bins in enumerate(bin_indexes):

        delay_slice_xx = calc_delay_spec_multi_baselines(all_visi_data_by_bin[bin_ind],
                                                         freqs, args)
        all_delay_spec_by_bin[:,bin_ind] = np.mean(delay_slice_xx, axis=0)

    return delays, all_delay_spec_by_bin


def plot_delay_spec(delays, bin_edges, all_delay_spec_by_bin, name,
                    vmin=False, vmax=False):
    """Plots the 2D delay spectra in `all_delay_spec_by_bin`, using `delays`
    and `bin_edges` to get the extent correct,  and outputs to `name`"""

    fig, ax = plt.subplots(1, 1, figsize=(6,6))

    if vmin and vmax:

        im = ax.imshow(np.log10(all_delay_spec_by_bin), origin='lower',
                        aspect="auto", cmap='Spectral_r',
                        extent=[bin_edges[0], bin_edges[-1],
                                delays[0]*1e+9, delays[-1]*1e+9],
                        vmin=np.log10(vmin), vmax=np.log10(vmax))
    else:
        im = ax.imshow(np.log10(all_delay_spec_by_bin), origin='lower',
                        aspect="auto", cmap='Spectral_r',
                        extent=[bin_edges[0], bin_edges[-1],
                                delays[0]*1e+9, delays[-1]*1e+9])

    add_colourbar(ax=ax, im=im, fig=fig, label=r'log10 P(b) mK$^2$ $h^{-3}$ Mpc$^3$')

    ax.set_xlabel("Baseline length (m)")
    ax.set_ylabel("Delay (ns)")

    ax.set_xscale("log")
    # ax.set_yscale("symlog")

    plt.tight_layout()
    fig.savefig(name,bbox_inches='tight')
    plt.close()



def get_uvfits_freqs(hdu):
    cent_freq = hdu[0].header['CRVAL4']
    ##subtract one because this is one indexed not zero
    cent_pix = hdu[0].header['CRPIX4'] - 1
    freq_res = hdu[0].header['CDELT4']

    num_freqs = hdu[0].data.data.shape[3]

    freqs = cent_freq + (np.arange(num_freqs) - cent_pix)*freq_res

    return freqs

def plot_delay_spec_multi_uvfits(args):
    """Gather data from multiple coarse band uvfits (ala RTS or WODEN)
    and plot a delay spectrum"""

    uvfits_prepend = args.uvfits_prepend
    num_bands = args.num_bands

    reso = 5
    # bin_edges = np.arange(reso, 1000 + reso, reso)
    bin_edges = np.arange(4, 500 + reso, reso)

    bin_indexes = []

    with fits.open(f"{uvfits_prepend}01.uvfits") as hdu:
        uu = hdu[0].data['UU']*SPEED_LIGHT
        vv = hdu[0].data['VV']*SPEED_LIGHT

        num_fine_chan = hdu[0].data.data.shape[3]

        lengths = np.sqrt(uu**2 + vv**2)

        for low_ind in range(len(bin_edges) - 1):
            bin_low = bin_edges[low_ind]
            bin_high = bin_edges[low_ind + 1]

            indexes = np.where((bin_low < lengths) & (lengths < bin_high))

            bin_indexes.append(indexes)

    freqs = np.empty(num_bands*num_fine_chan)

    all_visi_data_by_bin = []

    for bins in bin_indexes:

        delay_spec_bins = np.empty((len(bins[0]), num_bands*num_fine_chan), dtype=complex)

        all_visi_data_by_bin.append(delay_spec_bins)

    for band in range(num_bands):
        # print(f"Doing band {band+1}")

        with fits.open(f"{uvfits_prepend}{band+1:02d}.uvfits") as hdu:
            all_data = hdu[0].data.data[:,0,0,:,:,:]

            these_freqs = get_uvfits_freqs(hdu)

            low = band*num_fine_chan
            high = (band+1)*num_fine_chan

            freqs[low:high] = these_freqs

            for bin_ind, bins in enumerate(bin_indexes):
                data_slice_xx_re = all_data[bins[0], :, 0, 0]
                data_slice_xx_im = all_data[bins[0], :, 0, 1]

                data_slice_xx = data_slice_xx_re + 1j*data_slice_xx_im

                zeros = np.where(all_data[bins[0], :, 0, 2] == 0.0)
                data_slice_xx[zeros] = 0.0

                all_visi_data_by_bin[bin_ind][:, low:high] = data_slice_xx

    ##Sometimes the data out of the RTS are in descending freq order,
    ##so make sure re reorder #lesigh
    freq_order = np.argsort(freqs)

    for visi_data in all_visi_data_by_bin:
        visi_data = visi_data[:, freq_order]

    freqs = freqs[freq_order]

    delays, all_delay_spec_by_bin = do_delay_spec(bin_edges, bin_indexes,
                                                  all_visi_data_by_bin,
                                                  freqs, args)

    # print(all_delay_spec_by_bin.min(), all_delay_spec_by_bin.max())

    if args.output_name:
        name = args.output_name
    else:
        name = uvfits_prepend.split('/')[-1] + ".png"

    plot_delay_spec(delays, bin_edges, all_delay_spec_by_bin, name,
                    args.vmin, args.vmax)

    return all_delay_spec_by_bin


def plot_delay_spec_single_uvfits(args):
    """Gather visibilities from a single uvfits (ala hyperdrive) and
    plot a delay spectrum"""
    uvfits_name = args.uvfits

    reso = 10
    bin_edges = np.arange(reso, 1000 + reso, reso)

    bin_indexes = []

    with fits.open(uvfits_name) as hdu:
        uu = hdu[0].data['UU']*SPEED_LIGHT
        vv = hdu[0].data['VV']*SPEED_LIGHT

        all_data = hdu[0].data.data[:,0,0,:,:,:]

        num_fine_chan = hdu[0].data.data.shape[3]

        lengths = np.sqrt(uu**2 + vv**2)

        freqs = get_uvfits_freqs(hdu)

        for low_ind in range(len(bin_edges) - 1):
            bin_low = bin_edges[low_ind]
            bin_high = bin_edges[low_ind + 1]

            indexes = np.where((bin_low < lengths) & (lengths < bin_high))

            bin_indexes.append(indexes)

        all_visi_data_by_bin = []

        for bin_ind, bins in enumerate(bin_indexes):
            data_slice_xx_re = all_data[bins[0], :, 0, 0]
            data_slice_xx_im = all_data[bins[0], :, 0, 1]

            data_slice_xx = data_slice_xx_re + 1j*data_slice_xx_im

            zeros = np.where(all_data[bins[0], :, 0, 2] == 0.0)
            data_slice_xx[zeros] = 0.0

            all_visi_data_by_bin.append(data_slice_xx)

    delays, all_delay_spec_by_bin = do_delay_spec(bin_edges, bin_indexes,
                                                  all_visi_data_by_bin,
                                                  freqs, args)

    delay_cut = np.where(np.abs(delays) < 3000e-9)[0]
    
    print(delays.shape)
    print(all_delay_spec_by_bin.shape)

    delays = delays[delay_cut]
    all_delay_spec_by_bin = all_delay_spec_by_bin[delay_cut, :]

    if args.output_name:
        name = args.output_name
    else:
        name = uvfits_name.split('/')[-1][:-7] + ".png"

    plot_delay_spec(delays, bin_edges, all_delay_spec_by_bin, name,
                    args.vmin, args.vmax)

    return all_delay_spec_by_bin


def get_parser():
    """
    Runs the argument parser to get command line inputs - used by sphinx and
    argparse extension to unpack the help below into the online readthedocs
    documentation.

    Returns
    -------
    parser : `argparse.ArgumentParser`
        The populated argument parser used by `delay_spec_from_uvfits.py`

    """
    import argparse

    parser = argparse.ArgumentParser(description="Gimme some delay spectra.")


    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    parser.add_argument('--output_name', default=False,
        help='Name for output uvfits file, defaults to using the input uvfits name')
    parser.add_argument('--vmax', default=False, type=float,
        help='Maximum value for colourscale (e.g. 1e+16)')
    parser.add_argument('--vmin', default=False, type=float,
        help='Minimum value for colourscale (e.g. 1e+6)')

    band_group = parser.add_argument_group('DELAY FROM COARSE BAND UVFITS')
    band_group.add_argument('--uvfits_prepend', default=False,
        help='Prepend for a set of uvfits files e.g. ./data/uvdump_')
    band_group.add_argument('--num_bands', default=24, type=int,
        help='How many coarse bands to add into the delay spectra')

    sing_group = parser.add_argument_group('DELAY FROM SINGLE UVFITS')
    sing_group.add_argument('--uvfits', default=False,
        help='Single uvfits to make delay spectra from e.g. data/calibrated.uvfits'
             '(alternative to using --uvfits_prepend)')

    cosmo_group = parser.add_argument_group('COSMOLOGICAL CONSTANTS')
    cosmo_group.add_argument("--omega_matter", default=0.272, type=float,
        help="Matter density parameter. Default = 0.272")
    cosmo_group.add_argument("--omega_baryon", default=0.046, type=float,
        help="Baryon density parameter. Default = 0.046")
    cosmo_group.add_argument("--omega_lambda", default=0.7, type=float,
        help="Dark energy dentisty parameter. Default = 0.7")
    cosmo_group.add_argument("--hubble" ,default=70.4, type=float,
        help="Hubble param in km/s/Mpc, default=70.4")

    return parser

def main():
    """DO THE THINGS"""

    parser = get_parser()

    args = parser.parse_args()

    print('Cosmology being used has the following parameters:')
    print(f"\tH0 = {args.hubble:.2f} km / (Mpc s)")
    print(f"\tOmega Matter = {args.omega_matter:.4f}")
    print(f"\tOmega Lambda = {args.omega_lambda:.4f}")
    print(f"\tOmega Baryon = {args.omega_baryon:.4f}")

    if args.uvfits_prepend:
        plot_delay_spec_multi_uvfits(args)
    else:
        plot_delay_spec_single_uvfits(args)

if __name__ == '__main__':

    main()
