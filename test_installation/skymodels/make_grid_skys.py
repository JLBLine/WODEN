"""Make three different sky models, each in a grid, of point source, Gaussian
source, and shapelet sources. """

import numpy as np

np.random.seed(98743)

##ras are input as hours
ras = np.arange(58, 63) / 15.0
##dec are input as degrees
decs = np.arange(-42,-37)

source_ind = 0

##Iterate through all source positions and add in a new source at that position
##for all point, gaussian, and shapelet grid srclists
with open('srclist_point_source_grid.txt', 'w') as outpoint, \
     open('srclist_gauss_source_grid.txt', 'w') as outgauss, \
     open('srclist_shapelet_source_grid.txt', 'w') as outshape:

    for ra in ras:
        for dec in decs:
            outpoint.write('SOURCE {:02d} P 1 G 0 S 0 0\n'.format(source_ind))
            outpoint.write('COMPONENT POINT {:.8f} {:.7f}\n'.format(ra, dec))
            ##Set a source with a flat spectrum (SI = 0)
            outpoint.write('LINEAR 160e+6 1.0 0.0 0.0 0.0 0.0\n')
            outpoint.write('ENDCOMPONENT\n')
            outpoint.write('ENDSOURCE\n')

            outgauss.write('SOURCE {:02d} P 0 G 1 S 0 0\n'.format(source_ind))
            outgauss.write('COMPONENT GAUSSIAN {:.8f} {:.7f}\n'.format(ra, dec))
            ##Set a source with a flat spectrum (SI = 0)
            outgauss.write('LINEAR 160e+6 1.0 0.0 0.0 0.0 0.0\n')
            ##Give them a random PA to look cooler
            outgauss.write('GPARAMS {:.5f} 6.0 3.0\n'.format(np.random.uniform(0, 360)))
            outgauss.write('ENDCOMPONENT\n')
            outgauss.write('ENDSOURCE\n')

            outshape.write('SOURCE {:02d} P 0 G 0 S 1 22\n'.format(source_ind))
            outshape.write('COMPONENT SHAPELET {:.8f} {:.7f}\n'.format(ra, dec))
            ##Set a source with a flat spectrum (SI = 0)
            outshape.write('LINEAR 160e+6 1.0 0.0 0.0 0.0 0.0\n')
            ##These settings make a low-res PicA model at each location
            ##Give them a random PA to look cooler
            outshape.write('SPARAMS {:.5f} 2.0 1.5\n'.format(np.random.uniform(0, 360)))
            outshape.write('SCOEFF 12.0 0.0 1.582860954004\n')
            outshape.write('SCOEFF 8.0 0.0 1.499375585180\n')
            outshape.write('SCOEFF 4.0 6.0 -0.849974208075\n')
            outshape.write('SCOEFF 0.0 0.0 0.193440189215\n')
            outshape.write('SCOEFF 10.0 0.0 -1.987810132787\n')
            outshape.write('SCOEFF 5.0 5.0 0.152588248959\n')
            outshape.write('SCOEFF 4.0 4.0 0.604291976946\n')
            outshape.write('SCOEFF 6.0 4.0 -0.243182264246\n')
            outshape.write('SCOEFF 0.0 9.0 -0.065661839036\n')
            outshape.write('SCOEFF 3.0 7.0 0.024018473907\n')
            outshape.write('SCOEFF 4.0 0.0 0.111487284767\n')
            outshape.write('SCOEFF 3.0 5.0 -0.036255583012\n')
            outshape.write('SCOEFF 6.0 6.0 0.218560803924\n')
            outshape.write('SCOEFF 2.0 8.0 -0.097033457549\n')
            outshape.write('SCOEFF 5.0 7.0 -0.065892120142\n')
            outshape.write('SCOEFF 2.0 6.0 0.206179169091\n')
            outshape.write('SCOEFF 4.0 8.0 0.297133083972\n')
            outshape.write('SCOEFF 0.0 11.0 0.029920212310\n')
            outshape.write('SCOEFF 1.0 8.0 0.002164620831\n')
            outshape.write('SCOEFF 5.0 3.0 -0.104120035659\n')
            outshape.write('SCOEFF 2.0 0.0 0.172110209202\n')
            outshape.write('SCOEFF 0.0 7.0 0.023341102846\n')
            outshape.write('ENDCOMPONENT\n')
            outshape.write('ENDSOURCE\n')

            source_ind += 1
