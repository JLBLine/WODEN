import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from subprocess import call

VELC  = 299792458.0

def read_uvfits_data(filename):
    with fits.open(filename) as hdu:
        data = hdu[0].data.data[0,0,0,0,:,:]
        # print(data.shape)

        xx_re = data[0,0]
        xx_im = data[0,1]

        # print(xx_re)
        # print(xx_im)

        uu = hdu[0].data['UU']*VELC
        vv = hdu[0].data['VV']*VELC
        ww = hdu[0].data['WW']*VELC

        # print(uu, vv, ww)

    return xx_re, xx_im, uu

def make_command(ang_ind, num_mult):
    """Create the WODEN command"""

    ##want a frequency to make a wavelength of 1m, so set freq to speed of light
    ##--longitude=260.03426640562714 gives an LST=0.0 on UTC 2000-01-01T00:00:00

    cmd = "run_woden.py "
    cmd += f"--band_nums=1 --lowest_channel_freq={VELC:.1f} "
    cmd += "--coarse_band_width=1 --freq_res=1 "
    cmd += "--num_time_steps=1 --time_res=1e-9 "
    cmd += "--ra0=0.0 --dec0=0.0 --date=2000-01-01T00:00:00 "
    cmd += "--latitude=0.0 --longitude=260.03426640562714 --no_precession "
    cmd += f"--array_layout=array_layouts/array_layout_ang{ang_ind:02d}_n{num_mult:05d}.txt "
    cmd += "--primary_beam=none "
    cmd += f"--cat_filename=sky_models/srclist_ang{ang_ind:02d}.txt "
    cmd += f"--output_uvfits_prepend=uvfits/single_point_ang{ang_ind:02d}_n{num_mult:05d}"

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

# num_mults = [1]

if __name__ == '__main__':

    run_cuda = False

    u_lengths = []
    # expec_re = []
    # expec_im = []
    # calc_re = []
    # calc_im = []

    abs_diff_re = []
    abs_diff_im = []

    all_data = np.empty((len(phi_simples)*len(num_mults), len(num_mults)))
    print(all_data.shape)

    for ang_ind in range(len(phi_simples)):
    # for ang_ind in range(1):

        for num_ind, num_mult in enumerate(num_mults):

            cmd = make_command(ang_ind, num_mult)

            if run_cuda: call(cmd, shell=True)

            xx_re, xx_im, uu = read_uvfits_data(f"uvfits/single_point_ang{ang_ind:02d}_n{num_mult:05d}_band01.uvfits")

            u_lengths.append(np.sqrt(3*uu**2))
            # expec_re.append(known_cosine_angles[ang_ind])
            # expec_im.append(known_sine_angles[ang_ind])
            # calc_re.append(xx_re)
            # calc_im.append(xx_im)

            expec_re = known_cosine_angles[ang_ind]
            expec_im = known_sine_angles[ang_ind]

            # print(xx_re, expec_re)

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

    fig, axs = plt.subplots(1,2)


    axs[0].plot(u_lengths, abs_diff_re, 'C0o', mfc='none' )
    axs[1].plot(u_lengths, abs_diff_im, 'C0o', mfc='none' )

    axs[0].set_title("Real visibilities")
    axs[1].set_title("Imag visibilities")

    for ax in axs:
        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_xlabel('Baseline length (metres)')

    axs[0].set_ylabel('Percent fractional difference from expected')

    plt.tight_layout()
    fig.savefig('quantify_woden_accuracy.png',bbox_inches='tight')
    plt.close()

    np.savez_compressed("accuracy_test_outputs.npz",
                        u_lengths=u_lengths,
                        abs_diff_re=abs_diff_re,
                        abs_diff_im=abs_diff_im,
                        all_data=all_data)
