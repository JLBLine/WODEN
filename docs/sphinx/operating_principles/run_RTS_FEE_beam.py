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

D2R = np.pi / 180.0

nside = 301

lst_deg = 0.0

##Setup a dummy FITS header with appropriate settings
header = fits.Header()

##Give it 301 pixel for each axis
nside = 301

##This resolution seems to cover the full sky nicely
cpix = int(nside // 2)
cdelt = 0.25
cdelt = 125 / nside

header['NAXIS']   = 2
header['NAXIS1']  = nside
header['NAXIS2']  = nside
header['CTYPE1']  = 'RA---SIN'
header['CRPIX1']  = cpix
header['CRVAL1']  = lst_deg
header['CDELT1']  = cdelt
header['CUNIT1']  = 'deg     '
header['CTYPE2']  = 'DEC--SIN'
header['CRPIX2']  = cpix
header['CRVAL2']  = -26.7
header['CDELT2']  = cdelt
header['CUNIT2']  = 'deg     '

##Make a world coord system
wcs = WCS(header)

##Set up x/y pixels that cover the whole image
x_mesh, y_mesh = np.meshgrid(np.arange(nside), np.arange(nside))

x_pixels = x_mesh.flatten()
y_pixels = y_mesh.flatten()

##convert to ra, dec
ras, decs = wcs.all_pix2world(x_pixels, y_pixels, 0.0)

##Then use erfa to convert these values into azs, els
has = lst_deg - ras

##use this erfa function to convert to azimuth and elevation
##were using degrees, but erfa uses rads, so convert here
az_grid, els = erfa.hd2ae(has*D2R, decs*D2R, -26.7*D2R);

##convert elevation to zenith angle
za_grid = np.pi/2 - els


##Only feed az/za above the horizon to save on CUDA memory
##write out the az/za to feed into the C/CUDA code
with open('azza_values.txt','w') as outfile:
    for az, za in zip(az_grid[za_grid <= np.pi/2], za_grid[za_grid < np.pi/2]):
        outfile.write("{:.8f} {:.8f}\n".format(az, za))


# ##Mask below horizon for better plots
# below_horizon = np.where(za_grid > np.pi/2)
# az_grid[below_horizon] = np.nan
# za_grid[below_horizon] = np.nan


# ##compile the C/CUDA code
# call('source make_run_RTS_FEE_beam.sh', shell=True)
#
# ##run the C/CUDA code
# call('./run_RTS_FEE_beam', shell=True)


##read in outputs from C/CUDA code

used_az, used_za, gx_re, gx_im, Dx_re, Dx_im, Dy_re, Dy_im, gy_re, gy_im = np.loadtxt('MWAFEE_beamvalues_180MHz.txt',unpack=True)

def reshape_and_plot(data, ax, label, fig, vmin=False, vmax=False):

    square_plot = np.zeros(nside*nside)*np.nan

    square_plot[za_grid <= np.pi/2] = data

    square_plot.shape = (nside, nside)

    if vmin is not False and vmax is not False:
        im = ax.imshow(square_plot, origin='lower', vmin=vmin, vmax=vmax)
    else:
        im = ax.imshow(square_plot, origin='lower')

    add_colourbar(ax=ax, fig=fig, im=im)

    ax.set_title(label)

    ax.set_xticks([])
    ax.set_yticks([])


# za_grid.shape = (nside,nside)
# az_grid.shape = (nside,nside)

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
fig.savefig('MWAFEE_jones.png',bbox_inches='tight')
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
fig.savefig('MWAFEE_instrumental_pols.png',bbox_inches='tight')
plt.close()



freqs, gx_re, gx_im, Dx_re, Dx_im, Dy_re, Dy_im, gy_re, gy_im = np.loadtxt('MWAFEE_beamvalues-vs-freq.txt',unpack=True)

gx = gx_re + 1j*gx_im
Dx = Dx_re + 1j*Dx_im
Dy = Dy_re + 1j*Dy_im
gy = gy_re + 1j*gy_im

gx_conj = np.conjugate(gx)
Dx_conj = np.conjugate(Dx)
Dy_conj = np.conjugate(Dy)
gy_conj = np.conjugate(gy)

XX = (gx*gx_conj + Dx*Dx_conj)
YY = (Dy*Dy_conj + gy*gy_conj)

fig, ax = plt.subplots(1,1, figsize=(4,3))

ax.plot(freqs / 1e6, np.abs(XX),label='XX')
ax.plot(freqs / 1e6, np.abs(YY),label='YY')

ax.legend()

ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel('Beam gain')

plt.tight_layout()
fig.savefig('MWAFEE_beam_vs_freq.svg',bbox_inches='tight')
plt.close()
