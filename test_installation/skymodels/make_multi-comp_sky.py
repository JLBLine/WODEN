"""Make a big patch of sky with a component on every degree. Use a mix of
component types to make sure all components work together."""

import numpy as np

##centre grid around a certain point and extend up to an edge
ra0 = 60.0 / 15.0
dec0 = -40.0
edge = 20.0

##ras are input as hours
ras = ra0 + np.arange(-edge, edge + 1)/ 15.0
decs = dec0 + np.arange(-edge, edge + 1)


source_ind = 0
with open('srclist_multi-comp_grid.txt', 'w') as outfile:

    for ra in ras:
        for dec in decs:
            ##Make the top of the grid as point sources
            if dec > dec0:
                outfile.write('SOURCE {:02d} P 1 G 0 S 0 0\n'.format(source_ind))
                outfile.write('COMPONENT POINT {:.8f} {:.7f}\n'.format(ra, dec))
                ##Set a source with a flat spectrum (SI = 0)
                outfile.write('LINEAR 160e+6 1.0 0.0 0.0 0.0 0.0\n')
                outfile.write('ENDCOMPONENT\n')
                outfile.write('ENDSOURCE\n')

                source_ind += 1

            ##make the bottom two corners gaussian or shapelet
            else:

                if ra < ra0 + 15/15.0:
                    outfile.write('SOURCE {:02d} P 0 G 1 S 0 0\n'.format(source_ind))
                    outfile.write('COMPONENT GAUSSIAN {:.8f} {:.7f}\n'.format(ra, dec))
                    ##Set a source with a flat spectrum (SI = 0)
                    outfile.write('LINEAR 160e+6 1.0 0.0 0.0 0.0 0.0\n')
                    outfile.write('GPARAMS 45.0000000000 6.0 3.0\n')
                    outfile.write('ENDCOMPONENT\n')
                    outfile.write('ENDSOURCE\n')

                    source_ind += 1

                else:
                    outfile.write('SOURCE {:02d} P 0 G 0 S 1 6\n'.format(source_ind))
                    outfile.write('COMPONENT SHAPELET {:.8f} {:.7f}\n'.format(ra, dec))
                    ##Set a source with a flat spectrum (SI = 0)
                    outfile.write('LINEAR 160e+6 1.0 0.0 0.0 0.0 0.0\n')
                    ##Whack in a number of random components to make an odd
                    ##source
                    outfile.write('SPARAMS 45.0000000000 6.0 3.0\n')
                    outfile.write('SCOEFF 2 54 -7.15974974e-05\n')
                    outfile.write('SCOEFF 18 33 9.02112356e-05\n')
                    outfile.write('SCOEFF 56 12 -6.93747818e-05\n')
                    outfile.write('SCOEFF 27 41 1.97941884e-05\n')
                    outfile.write('SCOEFF 53 17 -2.39217077e-05\n')
                    outfile.write('SCOEFF 69 2 7.93187138e-06\n')
                    outfile.write('ENDCOMPONENT\n')
                    outfile.write('ENDSOURCE\n')

                    source_ind += 1

print(source_ind)
