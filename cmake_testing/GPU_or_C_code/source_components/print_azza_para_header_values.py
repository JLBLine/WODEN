import numpy as np
from astropy.coordinates import EarthLocation
import astropy.units as u
import erfa

MWA_LAT=-0.46606083776035967
RA0, DEC0 = 0.0, 0.0

if __name__ == "__main__":
    ras = np.array([(3*np.pi)/2, (5*np.pi)/3, (7*np.pi)/4, (11*np.pi)/6,
                   0.0, np.pi/6, np.pi/4, np.pi/3, np.pi/2])
    num_coords = len(ras)
    
    decs = np.zeros(num_coords)
    
    num_times = 3
    
    # azs = np.zeros(num_times*num_coords)
    # zas = np.zeros(num_times*num_coords)
    # para_angs = np.zeros(num_times*num_coords)
   
     


    lsts = ras
    
    azs, els = erfa.hd2ae(lsts - RA0, DEC0, 0.0)
    para_angles = erfa.hd2pa(lsts - RA0, DEC0, np.radians(MWA_LAT))
    
    zas = np.pi/2 - els
    
    azs = np.repeat(azs[::-1], num_times)
    zas = np.repeat(zas[::-1], num_times)
    para_angles = np.repeat(para_angles[::-1], num_times)
    
    arrays = [azs, zas, para_angles]
    array_names = ["azs", "zas", "para_angles"]
    
    for array, array_name in zip(arrays, array_names):

        arr_string = f"user_precision_t {array_name}[] = \u007b"
        for ind in range(len(array)-1):
            arr_string += f'{array[ind]:.8f}, '

        arr_string += f'{array[-1]:.8f}}};\n'

        print(arr_string)