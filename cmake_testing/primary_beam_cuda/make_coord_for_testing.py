from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from copy import deepcopy
import erfa
import matplotlib.pyplot as plt
from shamfi.shamfi_plotting import add_colourbar

# from rts_analytic_beam import RTS_analytic_beam_array

D2R = np.pi / 180.0
VELC = 299792458.0
MWA_LAT = -26.7033194444
MWA_LAT_RAD = MWA_LAT * D2R

lst_rad = 0.0
lst_deg = lst_rad / D2R

##Setup a dummy FITS header with appropriate settings
header = fits.Header()

##Give it 301 pixel for each axis
nside = 201

##This resolution seems to cover the full sky nicely
cpix = int(nside // 2) + 1
cdelt = 0.25
cdelt = 80 / nside

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
header['CRVAL2']  = MWA_LAT
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
az_grid, els = erfa.hd2ae(has*D2R, decs*D2R, MWA_LAT_RAD)



##convert elevation to zenith angle
za_grid = np.pi/2 - els


print(az_grid[-1], za_grid[-1])

za_grid.shape = (nside,nside)
az_grid.shape = (nside,nside)


##Mask below horizon for better plots
below_horizon = np.where(za_grid > np.pi/2)
az_grid[below_horizon] = np.nan
za_grid[below_horizon] = np.nan

az_flat = az_grid.flatten()
za_flat = za_grid.flatten()


with open('azza_radec_nside{:03d}.h'.format(nside),'w') as outfile:

    az_arr = 'user_precision_t nside{:03d}_azs[] = {{'.format(nside)
    za_arr = 'user_precision_t nside{:03d}_zas[] = {{'.format(nside)
    ra_arr = 'user_precision_t nside{:03d}_ras[] = {{'.format(nside)
    ha_arr = 'user_precision_t nside{:03d}_has[] = {{'.format(nside)
    dec_arr = 'user_precision_t nside{:03d}_decs[] = {{'.format(nside)

    for coord in range(len(az_flat)-1):
        az_arr += '{:.7f},'.format(az_flat[coord])
        za_arr += '{:.7f},'.format(za_flat[coord])
        ra_arr += '{:.7f},'.format(ras[coord]*D2R)
        ha_arr += '{:.7f},'.format(has[coord]*D2R)
        dec_arr += '{:.7f},'.format(decs[coord]*D2R)

    az_arr += '{:.7f}}};\n'.format(az_flat[-1])
    za_arr += '{:.7f}}};\n'.format(za_flat[-1])
    ra_arr += '{:.7f}}};\n'.format(ras[-1]*D2R)
    ha_arr += '{:.7f}}};\n'.format(has[-1]*D2R)
    dec_arr += '{:.7f}}};\n'.format(decs[-1]*D2R)

    outfile.write(az_arr)
    outfile.write(za_arr)
    outfile.write(ra_arr)
    outfile.write(ha_arr)
    outfile.write(dec_arr)

##Do some fiddling for plotting reasons
has_plot = deepcopy(has)
decs_plot = deepcopy(decs)
azs_plot = az_grid.flatten() / D2R
zas_plot = za_grid.flatten() / D2R

fig, axs = plt.subplots(2,2,figsize=(8,8))

arr_plots = [has_plot, decs_plot, azs_plot, zas_plot]
labels = ['HA (deg)', 'DEC (deg)', 'AZ (deg)', 'ZA (deg)']

for arr, ax, label in zip(arr_plots, axs.flatten(),labels):
    arr.shape = (nside,nside)
    im = ax.imshow(arr,origin='lower')

    ax.set_title(label)

    add_colourbar(ax=ax, im=im, fig=fig)

plt.tight_layout()

fig.savefig('coords_used_nside{:d}.png'.format(nside), bbox_inches='tight')


# delays = [0]*16
# delays = [6, 4, 2, 0, 8, 6, 4, 2, 10, 8, 6, 4, 12, 10, 8, 6]

# rts_jones = RTS_analytic_beam_array(az_flat, za_flat, delays, 100e+6, norm=True)
#
# xx = np.real(rts_jones[:,0])**2 + np.real(rts_jones[:,1])**2
# xy = np.real(rts_jones[:,0])*np.real(rts_jones[:,2]) + np.real(rts_jones[:,1])*np.real(rts_jones[:,3])
# yx = np.real(rts_jones[:,2])*np.real(rts_jones[:,0]) + np.real(rts_jones[:,3])*np.real(rts_jones[:,1])
# yy = np.real(rts_jones[:,3])**2 + np.real(rts_jones[:,2])**2
#
# xx.shape = (nside, nside)
#
#
# fig, axs = plt.subplots(2,2,figsize=(8,8))
#
# arr_plots = [xx, xy, yx, yy]
# labels = ['XX','XY','YX','YY']
#
# for arr, ax, label in zip(arr_plots, axs.flatten(), labels):
#     arr.shape = (nside,nside)
#     im = ax.imshow(np.log10(arr),origin='lower') #,vmin=0,vmax=0.3)
#
#     # print(arr.max())
#
#     ax.set_title(label)
#
#     add_colourbar(ax=ax, im=im, fig=fig)
#
# plt.tight_layout()
#
# fig.savefig('python_MWA_analytic_nside{:d}.png'.format(nside), bbox_inches='tight')
