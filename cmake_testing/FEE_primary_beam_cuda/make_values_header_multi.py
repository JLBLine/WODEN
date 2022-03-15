"""Make the header values used to test the gains in test_RTS_FEE_beam.c"""

import matplotlib.pyplot as plt
import numpy as np


with open('test_multifreq_RTS_FEE_beam.h', 'w') as outfile:

    outfile.write('#include "woden_precision_defs.h"\n\n')

    filenames = ['hyperbeam_multi_zenith_freqs1_rot.txt',
             'hyperbeam_multi_offzen1_freqs2_rot.txt',
             'hyperbeam_multi_offzen2_freqs3_rot.txt']

    array_names = ['hyper_f1d1', 'hyper_f2d2', 'hyper_f3d3']

    for filename, array_name in zip(filenames, array_names):
        az, za, g1r, g1i, D1r, D1i, D2r, D2i, g2r, g2i = np.loadtxt(filename, unpack=True)
        outfile.write('double {:s}[] = {{ {:.18f}, {:.18f}, {:.18f}, {:.18f},\n    {:.18f}, {:.18f}, {:.18f}, {:.18f},\n'.format(array_name,
                                            g1r[0], g1i[0], D1r[0], D1i[0], D2r[0], D2i[0], g2r[0], g2i[0]))

        for ind in range(1,len(az)-1):
            outfile.write('    {:.18f}, {:.18f}, {:.18f}, {:.18f},\n    {:.18f}, {:.18f}, {:.18f}, {:.18f},\n'.format(g1r[ind],
                                         g1i[ind], D1r[ind], D1i[ind], D2r[ind], D2i[ind], g2r[ind], g2i[ind]))

        outfile.write('    {:.18f}, {:.18f}, {:.18f}, {:.18f},\n    {:.18f}, {:.18f}, {:.18f}, {:.18f} }};\n'.format(g1r[-1],
                                     g1i[-1], D1r[-1], D1i[-1], D2r[-1], D2i[-1], g2r[-1], g2i[-1]))

        outfile.write('\n\n')
