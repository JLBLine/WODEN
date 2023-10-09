#!/usr/bin/env python3
from astropy.io import fits
import numpy as np
from sys import exit


def check_same_dimensions(hdu1, hdu2):
    """Bear minimum check that the uvfits data sections are the same size"""

    shape1 = hdu1[0].data.data.shape
    shape2 = hdu2[0].data.data.shape

    if shape1 != shape2:
        exit("Shape of the data in the uvfits are not the same. Data "
             "might have different number of times steps, frequencies, "
             "or baselines.")

def combine_uvfits(uvfits1, uvfits2, outname):
    """Opens uvfits1, uvfits2, reads in the visibilities from both, and
    sums them. Saves the result using the hdu of uvfits1, and leaves the
    weights of uvfits1 in the outputs. So output data have the array layout
    of uvfits1, u,v,w of uvfits1 etc, only visibilities are changed"""
    with fits.open(uvfits1) as hdu1, fits.open(uvfits2) as hdu2:

        check_same_dimensions(hdu1, hdu2)

        data1 = hdu1[0].data.data[:,0,0,:,:,:2]
        data2 = hdu2[0].data.data[:,0,0,:,:,:2]

        new_data = data1 + data2

        hdu1[0].data.data[:,0,0,:,:,:2] = new_data

        hdu1.writeto(outname, overwrite=True)

def combine_uvfits_bands(uvfits_prepend1, uvfits_prepend2, num_bands, output_name):
    """Runs `combine_uvfits` on multiple pairs of uvfits, with names generated
    by adding `{band}.uvfits` onto the end of `uvfits_prepend1` and
    `uvfits_prepend2`"""
    for band in range(1, num_bands + 1):

        combine_uvfits(f"{uvfits_prepend1}{band:02d}.uvfits",
                       f"{uvfits_prepend2}{band:02d}.uvfits",
                       f"{output_name}{band:02d}.uvfits")

def get_parser():
    """
    Runs the argument parser to get command line inputs - used by sphinx and
    argparse extension to unpack the help below into the online readthedocs
    documentation.

    Returns
    -------
    parser : `argparse.ArgumentParser`
        The populated argument parser used by `add_woden_uvfits.py`

    """
    import argparse

    parser = argparse.ArgumentParser(description="Sum the visibilities in a "
        "set of uvfits. There are two modes to add uvfits as decribed here. "
        "Mode 1: Assumes that the uvfits still end in '*{num}.uvfits', and so "
        "iterates over a given number of bands, so if you have 8 bands "
        "to combine, you end up with 8 summed uvfits files."
        "Mode 1 uses the `--num_bands`, `--uvfits_prepend1`, `--uvfits_prepend2`, "
        "`--output_name_prepend` arguments. Mode 2: This just adds two specific uvfits "
        "files together via `--uvfits1`, `--uvfits2` and `--output_name` "
        "DISCLAIMER - only very minimal checks are made that the two sets of "
        "visibilities make sense to be summed. Add uvfits at your own discretion.")

    band_group = parser.add_argument_group('ADDING COARSE BAND UVFITS')
    band_group.add_argument('--uvfits_prepend1', default=False,
        help='Prepend for the first set of uvfits files e.g. ./data/uvdump_')
    band_group.add_argument('--uvfits_prepend2', default=False,
        help='Prepend for the second set of uvfits files e.g. ./data2/uvdump_')
    band_group.add_argument('--output_name_prepend', default="combined_band",
        help='Name for output summed uvfits file, default: combined_band')
    band_group.add_argument('--num_bands', default=24, type=int,
        help='How many pairs of files to combine - default is 24')

    sing_group = parser.add_argument_group('ADDING SINGLE PAIR OF UVFITS')
    sing_group.add_argument('--uvfits1', default=False,
        help='Name of the first uvfits file to sum e.g. uvfits1.uvfits')
    sing_group.add_argument('--uvfits2', default=False,
        help='Name of the second uvfits file to sum e.g. uvfits2.uvfits')
    sing_group.add_argument('--output_name', default="combined.uvfits",
        help='Name for output summed uvfits file, default: combined.uvfits')

    return parser

def check_args(args):
    """Check that the args returned by """

    if not args.uvfits_prepend1 and not args.uvfits_prepend2 and not args.uvfits1 and not args.uvfits2:
       exit("You must set either a combination of --uvfits_prepend1 and "
            "--uvfits_prepend2 or --uvfits1 and --uvfits2 otherwise I don't "
            "know what to sum")

    if args.uvfits_prepend1 and not args.uvfits_prepend2:
        exit("If you supply --uvfits_prepend1 you must supply --uvfits_prepend2")
    if args.uvfits_prepend2 and not args.uvfits_prepend1:
        exit("If you supply --uvfits_prepend2 you must supply --uvfits_prepend1")

    if args.uvfits1 and not args.uvfits2:
        exit("If you supply --uvfits1 you must supply --uvfits2")
    if args.uvfits2 and not args.uvfits1:
        exit("If you supply --uvfits2 you must supply --uvfits1")


def main():
    """Runs all the things"""
    parser = get_parser()
    args = parser.parse_args()

    check_args(args)

    if args.uvfits_prepend1 and args.uvfits_prepend2:

        combine_uvfits_bands(args.uvfits_prepend1, args.uvfits_prepend2,
                      args.num_bands, args.output_name_prepend)
    else:
        combine_uvfits(args.uvfits1, args.uvfits2, args.output_name)


if __name__ == '__main__':
    main()
