#!/usr/bin/env python3
"""Script to feed to casa to convert uvfits files to measurement sets.
Run via `casa -c uv2ms.py --options`
"""
from subprocess import call
import numpy as np
import os
import argparse
import sys

def make_ms(uvfits_file, no_delete=False):
    """Takes a path to a uvfits file, checks it exists, and converts it to
    a measurement set of the same name at the same location"""

    ##If file exists, attempt to tranform it
    if os.path.isfile(uvfits_file):
        ##Get rid of '.uvfits' off the end of the file
        name = uvfits_file[:-7]
        ##Delete old measurement set if requested
        if args.no_delete:
            pass
        else:
            call("rm -r %s.ms" %name,shell=True)

        ##Do the conversion to measurement set
        importuvfits(fitsfile=uvfits_file, vis="{:s}.ms".format(name))
    else:
        ##print warning that uvfits file doesn't exist
        print("Could not find the uvfits specified by the path: "
              "\t{:s}. Skipping this tranformation".format(uvfits_file))

def get_parser():
    """
    Runs the argument parser to get command line inputs - used by sphinx and
    argparse extension to unpack the help below into the online readthedocs
    documentation.

    Returns
    -------
    parser : `argparse.ArgumentParser`
        The populated argument parser used by `uv2ms.py`

    """
    parser = argparse.ArgumentParser(description="Script to transform uvfits files"
                    "into measurement sets. Should be run via `casa -c uv2ms.py`."
                    "Uses same name of uvfits for the output measurement set")

    parser.add_argument('--single_uvfits', default=False,
        help="Convert this single uvfits file into a measurement set "
            "e.g. --single_uvfits=./data/nice_uvfits_band01.uvfits")
    parser.add_argument('--uvfits_prepend', default=False,
        help="Use in conjunction with `--band_nums` to process multiple uvfits "
            "files, where '{:02d}.uvfits'.format(band_num) is tacked onto "
            "the end of --uvfits_prepend. For example,  setting:\n"
            "\t--uvfits_prepend=example_band\n"
            "\t--band_nums=1,2\n"
            "will tranform both example_band01.uvfits and example_band02.uvfits")
    parser.add_argument('--band_nums', default='all',
        help="Add these numbers to the end of `--uvfits_prepend` to convert "
            "multiple uvfits files. Alternatively, enter required numbers delineated"
            " by commas, e.g. --band_nums=1,7,9")
    parser.add_argument('--no_delete', default=False, action='store_true',
        help="Defaults to deleting a measurement set with the destination name "
            "of the measurement set to be written. "
            " by commas, e.g. --band_nums=1,7,9")

    return parser

if __name__ == '__main__':
    ##Get command line arguments
    parser = get_parser()
    args = parser.parse_args()
    ##if single uvfits, convert single uvfits
    if args.single_uvfits:
        make_ms(args.single_uvfits, no_delete=args.no_delete)

    ##if multiple uvfits, work out which numbers to convert then do it
    if args.uvfits_prepend:
        if args.band_nums == 'all':
            args.band_nums = range(1,25)
        else:
            try:
                args.band_nums = list(np.array(args.band_nums.split(','),dtype=int))
            except:
                message = ("ERROR - failed to convert --band_nums into a list of ints"
                           " correctly. You entered:\n"
                           "    --band_nums={:s}\n"
                           "Exiting now.")
                exit(message)

        for band in args.band_nums:
            make_ms("{:s}{:02d}.uvfits".format(args.uvfits_prepend, band), no_delete=args.no_delete)
