"""Make a big patch of sky with a component on every degree. Use a mix of
component types to make sure all components work together."""

import numpy as np

np.random.seed(98743)

##centre grid around a certain point and extend up to an edge
ra0 = 60.0 / 15.0
dec0 = -40.0
edge = 10.0

##ras are input as hours
ras = ra0 + np.arange(-edge, edge + 1)/ 15.0
decs = dec0 + np.arange(-edge, edge + 1)


source_ind = 0
with open('srclist_multi-comp_grid.txt', 'w') as outfile:

    for ra in ras:
        for dec in decs:
            ##Make the top of the grid as point sources
            if dec > dec0:
                pass
                outfile.write('SOURCE {:02d} P 1 G 0 S 0 0\n'.format(source_ind))
                outfile.write('COMPONENT POINT {:.8f} {:.7f}\n'.format(ra, dec))
                ##Set a source with a flat spectrum (SI = 0)
                outfile.write('LINEAR 160e+6 1.0 0.0 0.0 0.0 0.0\n')
                outfile.write('ENDCOMPONENT\n')
                outfile.write('ENDSOURCE\n')

                source_ind += 1

            ##make the bottom two corners gaussian or shapelet
            else:

                if ra < ra0:
                    pass
                    outfile.write('SOURCE {:02d} P 0 G 1 S 0 0\n'.format(source_ind))
                    outfile.write('COMPONENT GAUSSIAN {:.8f} {:.7f}\n'.format(ra, dec))
                    ##Set a source with a flat spectrum (SI = 0)
                    outfile.write('LINEAR 160e+6 3.0 0.0 0.0 0.0 0.0\n')
                    outfile.write('GPARAMS {:.5f} 8.0 4.0\n'.format(np.random.uniform(0, 360)))
                    outfile.write('ENDCOMPONENT\n')
                    outfile.write('ENDSOURCE\n')

                    source_ind += 1

                else:
                    outfile.write('SOURCE {:02d} P 0 G 0 S 1 22\n'.format(source_ind))
                    outfile.write('COMPONENT SHAPELET {:.8f} {:.7f}\n'.format(ra, dec))
                    ##Set a source with a flat spectrum (SI = 0)
                    outfile.write('LINEAR 160e+6 3.0 0.0 0.0 0.0 0.0\n')
                    ##Whack in a number of random components to make an odd
                    ##source
                    outfile.write('SPARAMS {:.5f} 2.0 1.5\n'.format(np.random.uniform(0, 360)))
                    outfile.write('SCOEFF 12.0 0.0 1.582860954004\n')
                    outfile.write('SCOEFF 8.0 0.0 1.499375585180\n')
                    outfile.write('SCOEFF 4.0 6.0 -0.849974208075\n')
                    outfile.write('SCOEFF 0.0 0.0 0.193440189215\n')
                    outfile.write('SCOEFF 10.0 0.0 -1.987810132787\n')
                    outfile.write('SCOEFF 5.0 5.0 0.152588248959\n')
                    outfile.write('SCOEFF 4.0 4.0 0.604291976946\n')
                    outfile.write('SCOEFF 6.0 4.0 -0.243182264246\n')
                    outfile.write('SCOEFF 0.0 9.0 -0.065661839036\n')
                    outfile.write('SCOEFF 3.0 7.0 0.024018473907\n')
                    outfile.write('SCOEFF 4.0 0.0 0.111487284767\n')
                    outfile.write('SCOEFF 3.0 5.0 -0.036255583012\n')
                    outfile.write('SCOEFF 6.0 6.0 0.218560803924\n')
                    outfile.write('SCOEFF 2.0 8.0 -0.097033457549\n')
                    outfile.write('SCOEFF 5.0 7.0 -0.065892120142\n')
                    outfile.write('SCOEFF 2.0 6.0 0.206179169091\n')
                    outfile.write('SCOEFF 4.0 8.0 0.297133083972\n')
                    outfile.write('SCOEFF 0.0 11.0 0.029920212310\n')
                    outfile.write('SCOEFF 1.0 8.0 0.002164620831\n')
                    outfile.write('SCOEFF 5.0 3.0 -0.104120035659\n')
                    outfile.write('SCOEFF 2.0 0.0 0.172110209202\n')
                    outfile.write('SCOEFF 0.0 7.0 0.023341102846\n')

                    outfile.write('ENDCOMPONENT\n')
                    outfile.write('ENDSOURCE\n')

                    source_ind += 1

print(source_ind)
