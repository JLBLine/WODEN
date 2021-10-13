
def make_ax_do_plot(fig, az, za, plot_ind, data, numx=2, numy=4,
                    s=15, vmax=False, alpha=0.5 ):
    """Do a real simple polar scatter plot to show the beam. Crude but quick"""

    labels = ['XX real', 'XX imag', 'YY real', 'YY imag']
    ax = fig.add_subplot(numy, numx, plot_ind, projection='polar')

    ax.scatter(az, za, c=np.log10(data), s=s, alpha=alpha)

    ax.set_ylim(0,np.pi/2)
    ax.set_yticklabels([])

    ax.set_title(labels[ind])


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import numpy as np

    import sys

    ##Couldn't be bothered to google argparse so being lazy and using argv
    if len(sys.argv) != 2:
        sys.exit('Usage: python plot_beam_results.py [beam_values.txt]')

    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
        sys.exit('Usage: python plot_beam_results.py [beam_values.txt]')

    filename = sys.argv[1]

    ##Read in data
    az, za, g1r, g1i, D1r, D1i, D2r, D2i, g2r, g2i = np.loadtxt(filename, unpack=True)

    ##make complex
    g1 = g1r + 1j*g1i
    D1 = D1r + 1j*D1i
    D2 = D2r + 1j*D2i
    g2 = g2r + 1j*g2i

    ##Covert to XX and YY
    XX = g1*np.conjugate(g1) + D1*np.conjugate(D1)
    YY = g2*np.conjugate(g2) + D2*np.conjugate(D2)

    ##Plotting times
    fig = plt.figure(figsize=(10,10))

    data_sets = [np.real(XX), np.imag(XX), np.real(YY), np.imag(YY)]

    for ind, data in enumerate(data_sets):
        make_ax_do_plot(fig, az, za, ind + 1, data, numx=2, numy=2,
                        s=30, vmax=0.3, alpha=0.7)

    plt.tight_layout()

    plotname = filename.split('.')[0] + '.png'

    fig.savefig(plotname,bbox_inches='tight')
    plt.close()
