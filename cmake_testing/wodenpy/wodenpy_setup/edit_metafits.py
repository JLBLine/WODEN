from astropy.io import fits
import numpy as np

if __name__ == "__main__":
    with fits.open('1202815152_metafits_ppds.fits') as hdu:
        
        flags = hdu[1].data['Delays']
        
    with fits.open('1088285600_DipAmps.metafits') as hdu:
        
        hdu[1].data['Delays'] = flags
        
        hdu.writeto('1088285600_DipAmps_withflags.metafits', overwrite=True)
        
    
        