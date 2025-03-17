import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

known_angles = [0.0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, 2*np.pi/3,
                3*np.pi/4, 5*np.pi/6, np.pi,  7*np.pi/6, 5*np.pi/4]

known_angles_strings = ["$0.0$", "$\pi/6$", "$\pi/4$", "$\pi/3$",
                        "$\pi/2$", "$2\pi/3$", "$3\pi/4$", "$5\pi/6$",
                        "$\pi$",  "$7\pi/6$", "$5\pi/4$"]


markers = ["o", "v", "^", "<", ">", "8",
           "s", "p", "P", "*", "h", "H", "+"]


num_baselines = 5
num_angles = len(known_angles)

def plot_offsets(data, label, fig, axs):
    """Take the output data, cut them up into different phi_simple
    angles and plot the fractional offset from expectation as a function
    of baseline length"""

    all_re_diffs = np.empty(num_baselines*num_angles)
    all_im_diffs = np.empty(num_baselines*num_angles)

    num_neg = []
    num_pos = []

    for angle_ind, angle in enumerate(known_angles):

        slice_low = angle_ind*num_baselines
        slice_high = (angle_ind + 1)*num_baselines

        u_lens = data[slice_low:slice_high, 0]
        expec_re = data[slice_low:slice_high, 1]
        calc_re = data[slice_low:slice_high, 2]
        expec_im = data[slice_low:slice_high, 3]
        calc_im = data[slice_low:slice_high, 4]

        u_lens = np.sqrt(3*u_lens**2)

        expec_re_div = deepcopy(expec_re)
        expec_re_div[expec_re_div == 0] = 1.0

        expec_im_div = deepcopy(expec_im)
        expec_im_div[expec_im_div == 0] = 1.0

        abs_diff_re = np.abs((expec_re - calc_re) / expec_re_div)*100.0
        abs_diff_im = np.abs((expec_im - calc_im) / expec_im_div)*100.0

        for diff in expec_re - calc_re:
            if diff < 0:
                num_neg.append(diff)
            else:
                num_pos.append(diff)

        for diff in expec_im - calc_im:
            if diff < 0:
                num_neg.append(diff)
            else:
                num_pos.append(diff)

        all_re_diffs[slice_low:slice_high] = abs_diff_re
        all_im_diffs[slice_low:slice_high] = abs_diff_im

        if label == 'float':
            colour = 'C0'
        else:
            colour = 'C1'

        axs[0].plot(u_lens, abs_diff_re, color=colour,linestyle='none',
                    marker=markers[angle_ind],mfc='none',
                    label=label+ ' ' +known_angles_strings[angle_ind])
        axs[1].plot(u_lens, abs_diff_im, color=colour,linestyle='none',
                    marker=markers[angle_ind],mfc='none',
                    label=label+ ' ' +known_angles_strings[angle_ind])

    print(f"{label} Max fractional offset {all_re_diffs.max():.2E}%", )

    print(f"{label} Num positive offsets {len(num_pos)}")
    print(f"{label} Num negative offsets {len(num_neg)}")

if __name__ == '__main__':

    ##Load up the data, assuming you built everthing in ../../build
    data_float = np.loadtxt('../../../build/cmake_testing/GPU_or_C_code/source_components/measurement_eq_outcomes_float.txt')
    data_double = np.loadtxt('../../../build/cmake_testing/GPU_or_C_code/source_components/measurement_eq_outcomes_double.txt')

    ##Make the fig, axs
    fig, axs  = plt.subplots(1,2, figsize=(8.5,5.5))

    ##Do the plots
    plot_offsets(data_float, 'float', fig, axs)
    plot_offsets(data_double, 'double', fig, axs)

    ##Tidy up labels and make it a log/log plot
    for ax in axs.flatten():

        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_xlabel('$|\mathrm{u}| \,(\lambda)$')

    axs[1].legend(bbox_to_anchor=[1.3, 0.5], loc='center')

    axs[0].set_ylabel('Percentage difference from expected')

    axs[0].set_title("Real")
    axs[1].set_title("Imaginary")

    plt.tight_layout()
    fig.savefig('measure_eq_results.png',bbox_inches='tight',dpi=300)
    plt.close()
