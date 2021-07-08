"""Make three different sky models, each in a grid, of point source, Gaussian
source, and shapelet sources. """

import numpy as np

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
            outgauss.write('GPARAMS 45.0000000000 6.0 3.0\n')
            outgauss.write('ENDCOMPONENT\n')
            outgauss.write('ENDSOURCE\n')

            outshape.write('SOURCE {:02d} P 0 G 0 S 1 1\n'.format(source_ind))
            outshape.write('COMPONENT SHAPELET {:.8f} {:.7f}\n'.format(ra, dec))
            ##Set a source with a flat spectrum (SI = 0)
            outshape.write('LINEAR 160e+6 1.0 0.0 0.0 0.0 0.0\n')
            ##These settings mean the shapelet should come out the same as
            ##the Gaussian
            outshape.write('SPARAMS 45.0000000000 6.0 3.0\n')
            outshape.write('SCOEFF 0 0 1.0\n')
            outshape.write('ENDCOMPONENT\n')
            outshape.write('ENDSOURCE\n')

            source_ind += 1
