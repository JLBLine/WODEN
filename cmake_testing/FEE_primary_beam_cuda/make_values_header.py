import matplotlib.pyplot as plt
import numpy as np


with open('test_RTS_FEE_beam.h', 'w') as outfile:

    # names = ['zenith_100']

    names = ['zenith_100', 'zenith_150', 'zenith_200',
             'offzen1_100', 'offzen1_150', 'offzen1_200',
             'offzen2_100', 'offzen2_150', 'offzen2_200']

    for name in names:
        az, za, g1r, g1i, D1r, D1i, D2r, D2i, g2r, g2i = np.loadtxt('{:s}.txt'.format(name), unpack=True)
        outfile.write('float {:s}[] = {{ {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f},\n'.format(name,
                                            g1r[0], g1i[0], D1r[0], D1i[0], D2r[0], D2i[0], g2r[0], g2i[0]))

        for ind in range(1,len(az)-1):
            outfile.write('    {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f},\n'.format(g1r[ind],
                                         g1i[ind], D1r[ind], D1i[ind], D2r[ind], D2i[ind], g2r[ind], g2i[ind]))

        outfile.write('    {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f} }};\n'.format(g1r[-1],
                                     g1i[-1], D1r[-1], D1i[-1], D2r[-1], D2i[-1], g2r[-1], g2i[-1]))

        outfile.write('\n\n')
