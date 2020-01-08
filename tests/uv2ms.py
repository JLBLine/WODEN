from subprocess import call
from numpy import *
from sys import argv

def make_ms(uvfits_file):
    name = uvfits_file[:-7]
    call("rm -r %s.ms" %name,shell=True)
    importuvfits(fitsfile=uvfits_file, vis="%s.ms" %name)

for band in arange(1,4):
    make_ms("%s_band%02d.uvfits" %(argv[3],band))
