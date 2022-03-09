# from astropy.io import fits
import numpy as np
# from astropy.wcs import WCS
# from copy import deepcopy
# import erfa
# import matplotlib.pyplot as plt
# from shamfi.shamfi_plotting import add_colourbar

nside = 201

azs, zas, gx_re, gx_im, Dx_re, Dx_im, Dy_re, Dy_im, gy_re, gy_im, freqs = np.loadtxt('../../build/cmake_testing/primary_beam_cuda/MWA_analy_gains_azza40401.txt', unpack=True)

def make_C_array_string(data, name, precision='user_precision_t'):
    """Takes a numpy array `data` and creates a string to represent a C array
    of type `precision` with name `name`"""

    Carray = f'{precision} {name}[] = {{'

    for d in data[:-1]:
        Carray +=  f'{d:.8f},'

    Carray += f'{data[-1]:.8f}}};\n'

    return Carray

with open(f'MWA_analy_expected_nside{nside:03d}.h','w') as outfile:

    gx_re_array = make_C_array_string(gx_re, "gx_re_expec")
    Dx_re_array = make_C_array_string(Dx_re, "Dx_re_expec")
    Dy_re_array = make_C_array_string(Dy_re, "Dy_re_expec")
    gy_re_array = make_C_array_string(gy_re, "gy_re_expec")

    outfile.write(gx_re_array)
    outfile.write(Dx_re_array)
    outfile.write(Dy_re_array)
    outfile.write(gy_re_array)
