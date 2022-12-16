#!/usr/bin/env python3
from astropy.io import fits
import numpy as np

def check_uvfits_freq_order(uvfits_prepend, num_bands):
    """By default, WODEN should output the band numbers in ascending
    frequency, so check that this is true otherwise things are probably
    going to go wrong"""

    freqs = []

    for band in range(1, num_bands+1):

        with fits.open(f"{uvfits_prepend}{band:02d}.uvfits") as hdu:

            freqs.append(hdu[0].header['CRVAL4'])

    bands = np.arange(num_bands, dtype=int)
    freq_order = np.argsort(freqs)

    if np.array_equal(bands, freq_order):
        pass
    else:
        print("WARNING - the order of the frequencies in the uvfits files "
              "isn't the same as the ordering of the bands. Reordering "
              "them to be correct, but something is probably wrong and your "
              "outputs will probably suck.")

    return freq_order


def concat_uvfits(uvfits_prepend, bands, output_name, reverse_pols=False):
    """For band in `band_nums`, make a list of uvfits files with name
    `{uvfits_prepend}_band{band}.uvfits`. Then concanenate them by frequency,
    assuming that the set of uvfits are contiguous in frequency. Save output
    uvfits to `output_name`."""

    num_bands = bands[-1] - bands[0] + 1

    with fits.open(f"{uvfits_prepend}{bands[0]:02d}.uvfits") as hdu:
        data1 = hdu[0].data.data
        shape1 = data1.shape
        orig_hdu = hdu

        uu = hdu[0].data['UU']
        vv = hdu[0].data['VV']
        ww = hdu[0].data['WW']
        date_array = hdu[0].data['DATE'] - hdu[0].header['PZERO4']
        baselines_array = hdu[0].data['BASELINE']
        ant1_array = hdu[0].data['ANTENNA1']
        ant2_array = hdu[0].data['ANTENNA2']
        subs_array = hdu[0].data['SUBARRAY']

        ant_hdu = hdu[1]

        num_chans = data1.shape[3]

        all_data = np.empty((data1.shape[0], data1.shape[1], data1.shape[2],
                             num_bands*num_chans, data1.shape[4], data1.shape[5]))

        ##If swapping XX and YY (which really means swapping E-W with N-S)
        if reverse_pols:
            all_data[:,:,:,0:num_chans, 0, :] = data1[:,:,:,:,1,:]
            all_data[:,:,:,0:num_chans, 1, :] = data1[:,:,:,:,0,:]
            all_data[:,:,:,0:num_chans, 2, :] = data1[:,:,:,:,3,:]
            all_data[:,:,:,0:num_chans, 3, :] = data1[:,:,:,:,2,:]
        else:
            all_data[:,:,:,0:num_chans, :, :] = data1

        for band in bands[1:]:
            with fits.open(f"{uvfits_prepend}{band:02d}.uvfits") as this_hdu:
                this_data = this_hdu[0].data.data
                base_channel = num_chans*(band - 1)
                ##If swapping XX and YY (which really means swapping E-W with N-S)
                if reverse_pols:
                    all_data[:,:,:,base_channel:base_channel+num_chans,0,:] = this_data[:,:,:,:,1,:]
                    all_data[:,:,:,base_channel:base_channel+num_chans,1,:] = this_data[:,:,:,:,0,:]
                    all_data[:,:,:,base_channel:base_channel+num_chans,2,:] = this_data[:,:,:,:,3,:]
                    all_data[:,:,:,base_channel:base_channel+num_chans,3,:] = this_data[:,:,:,:,2,:]
                else:

                
                    all_data[:,:,:,base_channel:base_channel+num_chans,:,:] = this_data

        uvparnames = ['UU','VV','WW','DATE','BASELINE', 'ANTENNA1', 'ANTENNA2', 'SUBARRAY']

        parvals = [uu, vv, ww, date_array, baselines_array, ant1_array,
                   ant2_array, subs_array]

        ##This sets up the new size uvfits
        uvhdu = fits.GroupData(all_data, parnames=uvparnames, pardata=parvals, bitpix=-32)
        uvhdu = fits.GroupsHDU(uvhdu)

        uvhdu.header['CTYPE4'] = orig_hdu[0].header['CTYPE4']
        uvhdu.header['CRVAL4'] = orig_hdu[0].header['CRVAL4']  ##Middle pixel value in Hz
        uvhdu.header['CRPIX4'] = orig_hdu[0].header['CRPIX4'] ##Middle pixel number
        uvhdu.header['CDELT4'] = orig_hdu[0].header['CDELT4']

        ##Now just copy across all the old header info to the new header
        uvhdu.header['PSCAL1'] = orig_hdu[0].header['PSCAL1']
        uvhdu.header['PZERO1'] = orig_hdu[0].header['PZERO1']
        uvhdu.header['PSCAL2'] = orig_hdu[0].header['PSCAL2']
        uvhdu.header['PZERO2'] = orig_hdu[0].header['PZERO2']
        uvhdu.header['PSCAL3'] = orig_hdu[0].header['PSCAL3']
        uvhdu.header['PZERO3'] = orig_hdu[0].header['PZERO3']
        uvhdu.header['PSCAL4'] = orig_hdu[0].header['PSCAL4']
        uvhdu.header['PZERO4'] = orig_hdu[0].header['PZERO4']
        uvhdu.header['PSCAL5'] = orig_hdu[0].header['PSCAL5']
        uvhdu.header['PZERO5'] = orig_hdu[0].header['PZERO5']

        ###uvfits standards
        uvhdu.header['CTYPE2'] = orig_hdu[0].header['CTYPE2']
        uvhdu.header['CRVAL2'] = orig_hdu[0].header['CRVAL2']
        uvhdu.header['CRPIX2'] = orig_hdu[0].header['CRPIX2']
        uvhdu.header['CDELT2'] = orig_hdu[0].header['CDELT2']

        ##This means it's linearly polarised
        uvhdu.header['CTYPE3'] = orig_hdu[0].header['CTYPE3']
        uvhdu.header['CRVAL3'] = orig_hdu[0].header['CRVAL3']
        uvhdu.header['CRPIX3'] = orig_hdu[0].header['CRPIX3']
        uvhdu.header['CDELT3'] = orig_hdu[0].header['CDELT3']


        uvhdu.header['CTYPE5'] = orig_hdu[0].header['CTYPE5']
        uvhdu.header['CRVAL5'] = orig_hdu[0].header['CRVAL5']
        uvhdu.header['CRPIX5'] = orig_hdu[0].header['CRPIX5']
        uvhdu.header['CDELT5'] = orig_hdu[0].header['CDELT5']

        uvhdu.header['CTYPE6'] = orig_hdu[0].header['CTYPE6']
        uvhdu.header['CRVAL6'] = orig_hdu[0].header['CRVAL6']
        uvhdu.header['CRPIX6'] = orig_hdu[0].header['CRPIX6']
        uvhdu.header['CDELT6'] = orig_hdu[0].header['CDELT6']

        ##We're outputting into J2000
        uvhdu.header['EPOCH'] = orig_hdu[0].header['EPOCH']

        ##Old observation parameters that were/are needed in CHIPS
        uvhdu.header['OBJECT']  = orig_hdu[0].header['OBJECT']
        uvhdu.header['OBSRA']   = orig_hdu[0].header['OBSRA']
        uvhdu.header['OBSDEC']  = orig_hdu[0].header['OBSDEC']

        uvhdu.header['TELESCOP'] = orig_hdu[0].header['TELESCOP']
        uvhdu.header['INSTRUME'] = orig_hdu[0].header['INSTRUME']
        # uvhdu.header['SOFTWARE'] = orig_hdu[0].header['SOFTWARE']
        uvhdu.header['GITLABEL'] = orig_hdu[0].header['GITLABEL']

        uvhdu.header['DATEOBS'] = orig_hdu[1].header['RDATE']

        uvhdu.header['LAT'] = orig_hdu[0].header['LAT']
        uvhdu.header['LON'] = orig_hdu[0].header['LON']
        uvhdu.header['ALT'] = orig_hdu[0].header['ALT']
        uvhdu.header['INSTRUME'] = orig_hdu[0].header['INSTRUME']

        ##Before WODEN version 1.4.0, this header value wasn't
        ##used, so try and copy but don't crash if it isn't a key
        try:
            uvhdu.header['IAUORDER'] = orig_hdu[0].header['IAUORDER']
        
        except KeyError:
            pass

        ## Create hdulist and write out file
        hdulist = fits.HDUList(hdus=[uvhdu, orig_hdu[1]])
        hdulist.writeto(output_name, overwrite=True)
        hdulist.close()

def get_parser():
    """
    Runs the argument parser to get command line inputs - used by sphinx and
    argparse extension to unpack the help below into the online readthedocs
    documentation.

    Returns
    -------
    parser : `argparse.ArgumentParser`
        The populated argument parser used by `concat_woden_uvfits.py`

    """
    import argparse

    parser = argparse.ArgumentParser(description="Concatenate a number "
     "of uvfits files by frequency.")

    parser.add_argument('--num_bands', default=24, type=int,
        help='How many files to concatenate')
    parser.add_argument('--uvfits_prepend', default=False,
        help=r'Prepend for the uvfits files e.g. ./data/uvdump_')
    parser.add_argument('--output_name', default="concanenated.uvfits",
        help='Name for output concatenated uvfits file, default: concanenated.uvfits')
    parser.add_argument('--swap_pols', default=False, action='store_true',
        help='Reverse the order of polarisations')

    return parser


def main():
    """Runs all the things"""

    parser = get_parser()
    args = parser.parse_args()

    freq_order = check_uvfits_freq_order(args.uvfits_prepend, args.num_bands)

    bands = np.arange(1, args.num_bands+1)

    bands = bands[freq_order]

    concat_uvfits(args.uvfits_prepend, bands, args.output_name,
                  args.swap_pols)

if __name__ == '__main__':

    main()










#
