"""Plots outputs of the tests run by cmake (ctest). You gotta run that code
for this to woooooooork"""

from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from copy import deepcopy
import erfa
from subprocess import call
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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


##Read in the outputs of running `ctest`
nside = 101

used_az, used_za, gx_re, gx_im, Dx_re, Dx_im, Dy_re, Dy_im, gy_re, gy_im, freqs = np.loadtxt(f"../../../build/cmake_testing/primary_beam_cuda/hyperbeam_zenith_200_rot_double.txt",unpack=True)

print(used_za.max()*(180.0/np.pi))

def reshape_and_plot(data, ax, label, fig, vmin=False, vmax=False):

    square_plot = data[:nside*nside]

    square_plot.shape = (nside, nside)

    if vmin is not False and vmax is not False:
        im = ax.imshow(square_plot, origin='lower', vmin=vmin, vmax=vmax)
    else:
        im = ax.imshow(square_plot, origin='lower')

    add_colourbar(ax=ax, fig=fig, im=im)

    ax.set_title(label)

    ax.set_xticks([])
    ax.set_yticks([])


fig, axs = plt.subplots(4, 2, figsize=(6,10))

reshape_and_plot(gx_re, axs[0,0], 'Real $g_x$', fig)
reshape_and_plot(gx_im, axs[0,1], 'Imag $g_x$', fig)
reshape_and_plot(Dx_re, axs[1,0], 'Real $D_x$', fig)
reshape_and_plot(Dx_im, axs[1,1], 'Imag $D_x$', fig)
reshape_and_plot(Dy_re, axs[2,0], 'Real $D_y$', fig)
reshape_and_plot(Dy_im, axs[2,1], 'Imag $D_y$', fig)
reshape_and_plot(gy_re, axs[3,0], 'Real $g_y$', fig)
reshape_and_plot(gy_im, axs[3,1], 'Imag $g_y$', fig)

plt.tight_layout()
fig.savefig('hyperbeam_jones.png',bbox_inches='tight')
plt.close()

gx = gx_re + 1j*gx_im
Dx = Dx_re + 1j*Dx_im
Dy = Dy_re + 1j*Dy_im
gy = gy_re + 1j*gy_im

gx_conj = np.conjugate(gx)
Dx_conj = np.conjugate(Dx)
Dy_conj = np.conjugate(Dy)
gy_conj = np.conjugate(gy)

sI = complex(1,0)
sQ = complex(0,0)
sU = complex(0,0)
sV = complex(0,0)

XX = (gx*gx_conj + Dx*Dx_conj)*sI
XX += (gx*gx_conj - Dx*Dx_conj)*sQ
XX += (gx*Dx_conj + Dx*gx_conj)*sU
XX += 1j*(gx*Dx_conj - Dx*gx_conj)*sV

XY = (gx*Dy_conj + Dx*gy_conj)*sI
XY += (gx*Dy_conj - Dx*gy_conj)*sQ
XY += (gx*gy_conj + Dx*Dy_conj)*sU
XY += 1j*(gx*gy_conj - Dx*Dy_conj)*sV

YX = (Dy*gx_conj + gy*Dx_conj)*sI
YX += (Dy*gx_conj - gy*Dx_conj)*sQ
YX += (Dy*Dx_conj + gy*gx_conj)*sU
YX += 1j*(Dy*Dx_conj - gy*gx_conj)*sV

YY = (Dy*Dy_conj + gy*gy_conj)*sI
YY += (Dy*Dy_conj - gy*gy_conj)*sQ
YY += (Dy*gy_conj + gy*Dy_conj)*sU
YY += 1j*(Dy*gy_conj - gy*Dy_conj)*sV

fig, axs = plt.subplots(4, 2, figsize=(6,10))

reshape_and_plot(np.real(XX), axs[0,0], 'Real XX', fig, vmin=0, vmax=0.1)
reshape_and_plot(np.imag(XX), axs[0,1], 'Imag XX', fig)
reshape_and_plot(np.real(XY), axs[1,0], 'Real XY', fig)
reshape_and_plot(np.imag(XY), axs[1,1], 'Imag XY', fig)
reshape_and_plot(np.real(YX), axs[2,0], 'Real YX', fig)
reshape_and_plot(np.imag(YX), axs[2,1], 'Imag YX', fig)
reshape_and_plot(np.real(YY), axs[3,0], 'Real YY', fig, vmin=0, vmax=0.1)
reshape_and_plot(np.imag(YY), axs[3,1], 'Imag YY', fig)

plt.tight_layout()
fig.savefig('hyperbeam_linear_pols.png',bbox_inches='tight')
plt.close()



gx_re, gx_im, Dx_re, Dx_im, Dy_re, Dy_im, gy_re, gy_im, freqs = np.loadtxt(f"../../../build/cmake_testing/primary_beam_cuda/hyperbeam_interp_delays1_freqs1.txt",unpack=True)


coord_ind = 2
NUM_COORDS = 5
# freqs2 = np.arange(167e+6, 167e+6 + 12*1.28e+6, 1.28e+6)
# freqs = np.arange(190e+6, 190e+6 + 12*320e+3, 320e+3)

freqs = np.arange(167e+6, 167e+6 + 32*80e+3, 80e+3)

coord_slice = np.arange(coord_ind, len(freqs)*NUM_COORDS, NUM_COORDS)



gx = gx_re[coord_slice] + 1j*gx_im[coord_slice]
Dx = Dx_re[coord_slice] + 1j*Dx_im[coord_slice]
Dy = Dy_re[coord_slice] + 1j*Dy_im[coord_slice]
gy = gy_re[coord_slice] + 1j*gy_im[coord_slice]

gx_conj = np.conjugate(gx)
Dx_conj = np.conjugate(Dx)
Dy_conj = np.conjugate(Dy)
gy_conj = np.conjugate(gy)

XX = (gx*gx_conj + Dx*Dx_conj)
YY = (Dy*Dy_conj + gy*gy_conj)

fig, ax = plt.subplots(1,1, figsize=(5,4))

ax.plot(freqs / 1e6, np.real(XX), '-o', mfc='none', label='XX (real)')
ax.plot(freqs / 1e6, np.real(YY), '-o', mfc='none', label='YY (real)')

ax.legend()

ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel('Beam gain')

plt.tight_layout()
fig.savefig('hyperbeam_vs_freq.svg',bbox_inches='tight')
plt.close()
