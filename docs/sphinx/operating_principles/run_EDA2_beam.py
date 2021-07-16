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


def reshape_and_plot(data, ax, label, fig, vmin=False, vmax=False):

    square_plot = np.zeros(nside*nside)*np.nan

    square_plot[za_grid <= np.pi/2] = data

    square_plot.shape = (nside, nside)

    if vmin is not False and vmax is not False:
        im = ax.imshow(square_plot, origin='lower', vmin=vmin, vmax=vmax)
    else:
        im = ax.imshow(square_plot, origin='lower',cmap='gnuplot')

    add_colourbar(ax=ax, fig=fig, im=im)

    ax.set_title(label)

    ax.set_xticks([])
    ax.set_yticks([])

gx_re, gy_re = np.loadtxt('EDA2_beam_70MHz.txt',unpack=True)

fig, axs = plt.subplots(1, 2, figsize=(6,3))

reshape_and_plot(gx_re, axs[0], 'Real $g_x$', fig)
reshape_and_plot(gy_re, axs[1], 'Real $g_y$', fig)

plt.tight_layout()
fig.savefig('EDA2_jones.png',bbox_inches='tight')
plt.close()
