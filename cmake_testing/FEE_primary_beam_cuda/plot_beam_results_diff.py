import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

def make_ax_do_diff_plot(fig, az, za, plot_ind, data1, data2, numx=2, numy=2,
                    s=15, vmax=False, alpha=0.5 ):
    """Do a real simple polar scatter plot to show the beam. Crude but quick"""

    labels = ['gx real', 'gx real', 'Dx real', 'Dx imag',
              'Dy real', 'Dy imag', 'gy real', 'gy imag']
    ax = fig.add_subplot(numy, numx, plot_ind, projection='polar')

    data = data1 - data2

    print(f"Min / Max diff for {labels[ind]}: {data.min():.1e}, {data.max():.1e}")

    sc = ax.scatter(az, za, c=data, s=s, alpha=alpha, cmap='gnuplot')
    plt.colorbar(sc)

    # add_colourbar(ax=ax, fig=fig, im=sc)
    ax.set_ylim(0,np.pi/2)
    ax.set_yticklabels([])

    ax.set_title(labels[ind])

def make_ax_do_plot(fig, az, za, plot_ind, data, numx=2, numy=2,
                    s=15, vmax=False, alpha=0.5 ):
    """Do a real simple polar scatter plot to show the beam. Crude but quick"""

    labels = ['gx real', 'gx imag', 'Dx real', 'Dx imag',
              'Dy real', 'Dy imag', 'gy real', 'gy imag']
    ax = fig.add_subplot(numy, numx, plot_ind, projection='polar')

    sc = ax.scatter(az, za, c=data, s=s, alpha=alpha)
    plt.colorbar(sc)

    ax.set_ylim(0,np.pi/2)
    ax.set_yticklabels([])

    ax.set_title(labels[ind])

if __name__ == '__main__':

    ##Couldn't be bothered to google argparse so being lazy and using argv
    if len(sys.argv) != 3:
        sys.exit('Usage: python plot_beam_results.py [beam_values1.txt] [beam_values2.txt]')

    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
        sys.exit('Usage: python plot_beam_results.py [beam_values1.txt] [beam_values2.txt]')

    filename1 = sys.argv[1]
    filename2 = sys.argv[2]

    print("Doing these files", filename1, filename2)

    ##Read in data
    az_1, za_1, g1r_1, g1i_1, D1r_1, D1i_1, D2r_1, D2i_1, g2r_1, g2i_1 = np.loadtxt(filename1, unpack=True)
    az_2, za_2, g1r_2, g1i_2, D1r_2, D1i_2, D2r_2, D2i_2, g2r_2, g2i_2 = np.loadtxt(filename2, unpack=True)

    ##Plotting times
    fig = plt.figure(figsize=(8,16))

    data_sets1 = [g1r_1, g1i_1, D1r_1, D1i_1, D2r_1, D2i_1, g2r_1, g2i_1]
    data_sets2 = [g1r_2, g1i_2, D1r_2, D1i_2, D2r_2, D2i_2, g2r_2, g2i_2]

    for ind, data_set1 in enumerate(data_sets1):
        make_ax_do_diff_plot(fig, az_1, za_1, ind + 1, data_set1, data_sets2[ind],
                        numx=2, numy=4,
                        s=30, vmax=0.3, alpha=0.7)

    plt.tight_layout()

    plotname = "beam_plots/" + filename1.split('.')[0] + filename2.split('.')[0] + '_diff.png'

    fig.savefig(plotname,bbox_inches='tight')
    plt.close()


    ##Data set one
    fig = plt.figure(figsize=(8,16))

    for ind, data_set1 in enumerate(data_sets1):
        make_ax_do_plot(fig, az_1, za_1, ind + 1, data_set1,
                        numx=2, numy=4,
                        s=30, vmax=0.3, alpha=0.7)

    plt.tight_layout()

    plotname = "beam_plots/" + filename1.split('.')[0] + '_jones.png'

    fig.savefig(plotname,bbox_inches='tight')
    plt.close()

    ##Data set two
    fig = plt.figure(figsize=(8,16))

    for ind, data_set2 in enumerate(data_sets2):
        make_ax_do_plot(fig, az_1, za_1, ind + 1, data_set2,
                        numx=2, numy=4,
                        s=30, vmax=0.3, alpha=0.7)

    plt.tight_layout()

    plotname = "beam_plots/" + filename2.split('.')[0] + '_jones.png'

    fig.savefig(plotname,bbox_inches='tight')
    plt.close()
