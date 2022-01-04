import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from subprocess import call
import os
import argparse

DD2R = np.pi/180.0
VELC  = 299792458.0

MWA_LAT = -26.703319
MWA_LONG = 116.67081524

def read_uvfits_data(filename):
    """Reads in data from the uvfits files. As we run with no primary beam,
    only need to look at one polarisation so just grab the xx"""
    with fits.open(filename) as hdu:
        data = hdu[0].data.data[0,0,0,0,:,:]

        xx_re = data[0,0]
        xx_im = data[0,1]

        uu = hdu[0].data['UU']*VELC
        vv = hdu[0].data['VV']*VELC
        ww = hdu[0].data['WW']*VELC

    return xx_re, xx_im, uu

def make_command(ang_ind, num_mult, precision):
    """Create the WODEN command with the """

    ##This will NOT work for anyone else but Jack - needed it to make
    ##the comparison to the old fully float code
    if precision == "fully_float":
        os.environ["WODEN_DIR"] = "/home/jline/software/WODEN_master/build/"
        cmd = "/home/jline/software/WODEN_master/build/run_woden.py "
    else:
        cmd = "run_woden.py "
        cmd += f"--precision={precision} "
    cmd += f"--band_nums=1 --lowest_channel_freq={VELC:.1f} "
    cmd += "--coarse_band_width=1 --freq_res=1 "
    cmd += "--num_time_steps=1 --time_res=1e-9 "
    cmd += "--ra0=0.0 --dec0=0.0 --date=2020-01-01T12:00:00.0 "
    cmd += f"--latitude=0.1095073835963605 --longitude=79.6423588359480874 "
    cmd += f"--array_layout=array_layouts/array_layout_ang{ang_ind:02d}_n{num_mult:05d}.txt "
    cmd += "--primary_beam=none "
    cmd += f"--cat_filename=sky_models/srclist_ang{ang_ind:02d}.txt "
    cmd += f"--output_uvfits_prepend=uvfits/single_point_ang{ang_ind:02d}_n{num_mult:05d}_{precision} "
    cmd += "--band_nums=1"

    return cmd

phi_simples = [0.0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, 2*np.pi/3,
               3*np.pi/4, 5*np.pi/6, np.pi, 7*np.pi/6, 5*np.pi/4]

##Known sin outputs for input phi_simples
known_sine_angles = [0.0, 0.5, np.sqrt(2)/2, np.sqrt(3)/2, 1.0, np.sqrt(3)/2,
                     np.sqrt(2)/2, 0.5, 0.0, -0.5, -np.sqrt(2)/2]

##Known sin outputs for input phi_simples
known_cosine_angles = [1.0, np.sqrt(3)/2, np.sqrt(2)/2, 0.5, 0.0, -0.5,
                      -np.sqrt(2)/2, -np.sqrt(3)/2, -1.0, -np.sqrt(3)/2,
                      -np.sqrt(2)/2]

num_mults = [1, 10, 100, 1000, 10000]
# num_mults = [10000]

def run_simulation(precision):
    """Runs the simulation to the given precision"""

    call("mkdir -p uvfits", shell=True)

    u_lengths = []

    abs_diff_re = []
    abs_diff_im = []

    all_data = np.empty((len(phi_simples)*len(num_mults), 5))

    for ang_ind in range(len(phi_simples)):
    # for ang_ind in range(10,11):
        for num_ind, num_mult in enumerate(num_mults):

            cmd = make_command(ang_ind, num_mult, precision)

            call(cmd, shell=True)

            uvfits_name = f"uvfits/single_point_ang{ang_ind:02d}_n{num_mult:05d}_{precision}_band01.uvfits"

            xx_re, xx_im, uu = read_uvfits_data(uvfits_name)

            u_lengths.append(np.sqrt(3*uu**2))

            expec_re = known_cosine_angles[ang_ind]
            expec_im = known_sine_angles[ang_ind]

            if expec_re == 0:
                diff_re = np.abs((xx_re - expec_re))
            else:
                diff_re = np.abs((xx_re - expec_re)/ expec_re)

            if expec_im == 0:
                diff_im = np.abs((xx_im - expec_im))
            else:
                diff_im = np.abs((xx_im - expec_im)/ expec_im)

            abs_diff_re.append(diff_re*100.0)
            abs_diff_im.append(diff_im*100.0)

            # print(type(ang_ind*len(phi_simples) + num_ind))

            all_data[ang_ind*len(num_mults) + num_ind, 0] = np.sqrt(3*uu**2)
            all_data[ang_ind*len(num_mults) + num_ind, 1] = diff_re*100.0
            all_data[ang_ind*len(num_mults) + num_ind, 2] = diff_im*100.0
            all_data[ang_ind*len(num_mults) + num_ind, 3] = xx_re
            all_data[ang_ind*len(num_mults) + num_ind, 4] = xx_im

            ##tidy up the data, no need to store the uvfits
            call(f"rm {uvfits_name}", shell=True)

    np.savez_compressed(f"accuracy_test_outputs_{precision}.npz",
                        u_lengths=u_lengths,
                        abs_diff_re=abs_diff_re,
                        abs_diff_im=abs_diff_im,
                        all_data=all_data)


def plot_individual_precision(fig, axs, label, filename):

    data = np.load(filename)

    if label == 'float (mix)':
        marker = 'C1s'
        ms = 7
    elif label == 'double':
        marker = 'C0x'
        ms = 4
    else:
        marker = 'C2^'
        ms = 7

    # print(len(data['u_lengths']))

    axs[0].plot(data['u_lengths'], data['abs_diff_re'], marker, mfc='none', ms=ms,
                label=label)
    axs[1].plot(data['u_lengths'], data['abs_diff_im'], marker, mfc='none', ms=ms,
                label=label)

def plot_accuracy(float_filename, double_filename, float_full_filename):
    fig, axs = plt.subplots(1,2, figsize=(8,4.0))

    plot_individual_precision(fig, axs, "float (full)", float_full_filename)
    plot_individual_precision(fig, axs, "float (mix)", float_filename)
    plot_individual_precision(fig, axs, "double", double_filename)

    for ax in axs:
        ax.axhline(0.2, color='k', alpha=0.4, label='0.2% error', linestyle='--')
        small_err = 0.000002
        ax.axhline(small_err, color='k', alpha=0.4, label=f'{small_err:.0e}% error',
                   linestyle=':')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Baseline length (metres)')
        ax.legend(ncol=2)

        ax.set_ylim(5e-7, 1e3)

        ax.set_xlim(1, 1e6)

    axs[0].set_ylabel('Percent fractional difference from expected')

    axs[0].set_title('Real component')
    axs[1].set_title('Imaginary component')

    plt.tight_layout()
    fig.savefig('quantify_woden_accuracy.png',bbox_inches='tight')
    plt.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Run the WODEN simulator and profit. ")

    parser.add_argument('--only_plot', default=False, action='store_true',
        help='Whether to run the simulations or not. If added, just look for '
             'previous outputs and make the plot')

    args = parser.parse_args()

    if args.only_plot:
        pass
    else:

        run_simulation('float')
        run_simulation('double')

        ##This next line won't work for you, I'm just leaving this here to show how I got the
        ##old fully float results, using an old version of master
        # run_simulation('fully_float')

    plot_accuracy("accuracy_test_outputs_float.npz", "accuracy_test_outputs_double.npz",
                  "accuracy_test_outputs_fully_float.npz")
