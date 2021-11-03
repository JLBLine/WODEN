"""Make the header values used to test the gains in test_RTS_FEE_beam.c"""

import matplotlib.pyplot as plt
import numpy as np


with open('test_RTS_FEE_beam.h', 'w') as outfile:

    outfile.write('#include "woden_precision_defs.h"\n\n')

    names = ['zenith_100', 'zenith_150', 'zenith_200',
             'offzen1_100', 'offzen1_150', 'offzen1_200',
             'offzen2_100', 'offzen2_150', 'offzen2_200',
             'zenith_100_rot', 'zenith_150_rot', 'zenith_200_rot',
             'offzen1_100_rot', 'offzen1_150_rot', 'offzen1_200_rot',
             'offzen2_100_rot', 'offzen2_150_rot', 'offzen2_200_rot']

    for name in names:
        az, za, g1r, g1i, D1r, D1i, D2r, D2i, g2r, g2i = np.loadtxt('hyperbeam_{:s}.txt'.format(name), unpack=True)
        outfile.write('user_precision_t {:s}[] = {{ {:.18f}, {:.18f}, {:.18f}, {:.18f},\n    {:.18f}, {:.18f}, {:.18f}, {:.18f},\n'.format(name,
                                            g1r[0], g1i[0], D1r[0], D1i[0], D2r[0], D2i[0], g2r[0], g2i[0]))

        for ind in range(1,len(az)-1):
            outfile.write('    {:.18f}, {:.18f}, {:.18f}, {:.18f},\n    {:.18f}, {:.18f}, {:.18f}, {:.18f},\n'.format(g1r[ind],
                                         g1i[ind], D1r[ind], D1i[ind], D2r[ind], D2i[ind], g2r[ind], g2i[ind]))

        outfile.write('    {:.18f}, {:.18f}, {:.18f}, {:.18f},\n    {:.18f}, {:.18f}, {:.18f}, {:.18f} }};\n'.format(g1r[-1],
                                     g1i[-1], D1r[-1], D1i[-1], D2r[-1], D2i[-1], g2r[-1], g2i[-1]))

        outfile.write('\n\n')
