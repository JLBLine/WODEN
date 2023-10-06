from wodenpy.array_layout.precession import RTS_Precess_LST_Lat_to_J2000
from astropy.time import Time, TimeDelta
import numpy as np
from astropy.coordinates import EarthLocation
from astropy import units as u
import scipy.optimize as opt

def get_lat_lon(inputs, lat, lon):
    
    # lat, lon = inputs
    
    height = inputs
    date = "2020-01-01T12:00:00.0"
    
    
    
    ##Setup location
    observing_location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height)
    ##Setup time at that location

    lst_type = 'mean'

    observing_time = Time(date, scale='utc', location=observing_location)
    ##Grab the LST
    LST = observing_time.sidereal_time(lst_type)
    LST_deg = LST.value*15.0
    
    lst_current = LST_deg*(np.pi/180.0)
    latitude_current = lat*(np.pi/180.0)
    
    t = Time(date)
    mjd = t.jd - 2400000.5

    lst_J2000, latitude_J2000 = RTS_Precess_LST_Lat_to_J2000(lst_current, latitude_current,
                                 mjd)
    
    return lst_J2000, latitude_J2000
    
##Do an optimisation to grab a good lat,lon to put into the simulation
if __name__ == "__main__":
    
    intitial_guess = [0.1095073835963605, 79.6423588359480874]
    
    desried_lat_lon = [0.0, 0.0]
    
    height=0.0
    date="2020-01-01T12:00:00.0"
    
    inputs = height
    
    popt, pcov = opt.curve_fit(get_lat_lon, inputs, desried_lat_lon, p0=intitial_guess)
    
    print(popt[0])
    print(popt[1])