#!/usr/bin/env python
from __future__ import print_function,division,absolute_import
from numpy import *
import argparse
import os
from subprocess import check_output
from astropy.coordinates import SkyCoord
from astropy import units as u

##Find out where the git repo is, cd in and grab the git label
##TODO do this in a better way
fileloc = os.path.realpath(__file__)
cwd = os.getcwd()
os.chdir(('/').join(fileloc.split('/')[:-1]))
gitlabel = check_output(["git", "describe", "--always"],universal_newlines=True).strip()
##Get back to where we were before
os.chdir(cwd)

parser = argparse.ArgumentParser(description="Converts a multiscale component model output by WSClean into a WODEN style srclist")
parser.add_argument('--file', help='WSClean source list to convert')
parser.add_argument('--outname',default=False,help='Name for output srclist - defaults to using the current name of the file in output name')
args = parser.parse_args()

##Format of WSClean source list at time of writing of this script
##Format = Name, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='188790000', MajorAxis, MinorAxis, Orientation

##Load data and name
data = genfromtxt(args.file, usecols=(1,2,3,4,7,8,9,10),skip_header=1,delimiter=',',
                  names=['type','ra','dec','flux','frequency','major','minor','pa'],
                  dtype="S16,S16,S16,f16,f16,f16,f16,f16")

ras = data['ra']
decs = data['dec']
##Reformat decs into something astropy can read
decs = array(["%02d:%02d:%06.4f" %(int(dec.split('.')[0]),int(dec.split('.')[1]),float('.'.join(dec.split('.')[2:]))) for dec in decs])

##Set up astropy coords to convert coords into degrees
coords = SkyCoord(ras,decs,unit=(u.hourangle, u.deg))
##Converts into degress
ras = coords.ra.value
decs = coords.dec.value

##Find the number of each type of component
num_gauss = len(where(data['type'] == 'GAUSSIAN')[0])
num_point = len(where(data['type'] == 'POINT')[0])

print('There are %d point sources and %d gaussians' %(num_point,num_gauss))

##If a name given for srclist then use it
##Otherwise use the current name and modify
if args.outname:
    outname = args.outname
else:
    name_strip = args.file.split('/')[-1]
    outname = 'srclist-woden_wsclean%s' %name_strip

##Write the WODEN srclist
with open(outname,'w+') as outfile:
    outfile.write('SOURCE msclean_source P %d G %d S 0 0\n' %(num_point,num_gauss))
    for component,type in enumerate(data['type']):
        if type == 'POINT':
            outfile.write('COMPONENT POINT %.6f %.5f\n' %(ras[component]/15.0,decs[component]))
            outfile.write('FREQ %.5e %.10f 0 0 0\n' %(data['frequency'][component],data['flux'][component]))
            outfile.write('ENDCOMPONENT\n')
        elif type == 'GAUSSIAN':
            outfile.write('COMPONENT GAUSSIAN %.6f %.5f\n' %(ras[component]/15.0,decs[component]))
            outfile.write('FREQ %.5e %.10f 0 0 0\n' %(data['frequency'][component],data['flux'][component]))
            outfile.write('GPARAMS %.5f %.3f %.3f\n' %(data['pa'][component],data['major'][component] / 60.0,data['minor'][component] / 60.0))
            outfile.write('ENDCOMPONENT\n')

    outfile.write('ENDSOURCE')
