import matplotlib.pyplot as plt
import numpy as np

NUM_COORDS = 5

def plot_coord(axs, freqs, gx_re, gx_im, gy_re, gy_im,
               style, label):

    for coord_ind in np.arange(NUM_COORDS):

        coord = np.arange(coord_ind, len(freqs)*NUM_COORDS, NUM_COORDS)

        axs[coord_ind, 0].plot(freqs, gx_re[coord], style, mfc='none', label=label)
        axs[coord_ind, 1].plot(freqs, gx_im[coord], style, mfc='none', label=label)
        axs[coord_ind, 2].plot(freqs, gy_re[coord], style, mfc='none', label=label)
        axs[coord_ind, 3].plot(freqs, gy_im[coord], style, mfc='none', label=label)

def plot_hyperdrive(fig, axs, freqs, filename):

    data = np.loadtxt(filename)

    gx_re = data[:,2]
    gx_im = data[:,3]

    gy_re = data[:,8]
    gy_im = data[:,9]

    style = 'C0-'
    label = 'Hyperbeam'

    plot_coord(axs, freqs, gx_re, gx_im, gy_re, gy_im, style, label)

def plot_woden(fig, axs, freqs, filename):

    data = np.loadtxt(filename)

    gx_re = data[:,0]
    gx_im = data[:,1]

    gy_re = data[:,6]
    gy_im = data[:,7]

    style = 'C1--'
    label = 'WODEN'

    plot_coord(axs, freqs, gx_re, gx_im, gy_re, gy_im, style, label)

def plot_pointing_and_freq(hyperfile, woden_file, freqs, plotname):

    fig, axs = plt.subplots(5,4, figsize=(10, 12))

    plot_hyperdrive(fig, axs, freqs / 1e+6, hyperfile)
    plot_woden(fig, axs, freqs / 1e+6, woden_file)

    axs[0,0].set_title('Real gx')
    axs[0,1].set_title('Imag gx')
    axs[0,2].set_title('Real gy')
    axs[0,3].set_title('Imag gy')


    for ax in axs[:,0]: ax.set_ylabel('Gain (normed)')
    for ax in axs[4,:]: ax.set_xlabel('Frequency (MHz)')

    for ax in axs.flatten(): ax.legend()

    plt.tight_layout()
    fig.savefig(plotname, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    freqs1 = np.arange(167e+6, 167e+6 + 32*80e+3, 80e+3)
    freqs2 = np.arange(167e+6, 167e+6 + 12*1.28e+6, 1.28e+6)
    freqs3 = np.arange(190e+6, 190e+6 + 12*320e+3, 320e+3)

    plot_pointing_and_freq('hyperbeam_multi_zenith_freqs1_rot.txt',
        "MWA_FEE_multifreq_gains_delays1_freqs1.txt",
        freqs1, "multi_zenith_freqs1.png")

    plot_pointing_and_freq('hyperbeam_multi_offzen1_freqs2_rot.txt',
        "MWA_FEE_multifreq_gains_delays2_freqs2.txt",
        freqs2, "multi_offzen1_freqs2.png")

    plot_pointing_and_freq('hyperbeam_multi_offzen2_freqs3_rot.txt',
        "MWA_FEE_multifreq_gains_delays3_freqs3.txt",
        freqs3, "multi_offzen2_freqs3.png")
