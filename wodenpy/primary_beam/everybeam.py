import everybeam as eb
import numpy as np
from astropy.coordinates import ITRS, SkyCoord #AltAz, EarthLocation, 
from astropy.time import Time, TimeDelta
import astropy.units as u
import argparse
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes

##This call is so we can use it as a type annotation
woden_struct_classes = Woden_Struct_Classes()
Source_Catalogue = woden_struct_classes.Source_Catalogue
Woden_Settings = woden_struct_classes.Woden_Settings


# https://everybeam.readthedocs.io/en/latest/tree/demos/lofar-array-factor.html

def radec_to_xyz(ra : float, dec : float, time : Time):
    """
    Convert RA and Dec ICRS coordinates to ITRS cartesian coordinates.

    Args:
        ra (astropy.coordinates.Angle): Right ascension
        dec (astropy.coordinates.Angle): Declination
        time (float): astropy time instance

    Returns:
        pointing_xyz (ndarray): NumPy array containing the ITRS X, Y and Z coordinates
    """
    
    
    coord = SkyCoord(ra*u.rad, dec*u.rad, frame='icrs')
    coord_itrs = coord.transform_to(ITRS(obstime=time))
    
    return np.asarray(coord_itrs.cartesian.xyz.transpose())

def load_OSKAR_telescope(ms_path):

    ##TODO what does this mean, what is differential beam?
    use_differential_beam = False

    # Set element response model
    element_response_model = "skala40_wave"

    # Load the telescope
    telescope = eb.load_telescope(
        ms_path,
        use_differential_beam=use_differential_beam,
        element_response_model=element_response_model,
    )
    
    # assert type(telescope) == eb.OSKAR
    if type(telescope) != eb.OSKAR:
        print(f'WARNING: Telescope specified in {ms_path} is not an OSKAR telescope. Proceeding, but you might get nonsense results.')
    
    return telescope


def get_everybeam_norm(ra0, dec0, time, freq, telescope, station_id = 0):
    
    phase_itrf = radec_to_xyz(ra0, dec0, time)
    
    # Full beam for station 0
    response = telescope.station_response(time.mjd*3600*24, station_id,
                                          freq, phase_itrf, phase_itrf)
    
    norm_x = np.abs(response[0,0])
    norm_y = np.abs(response[1,1])
    
    return norm_x, norm_y

def run_everybeam(ra, dec, ra0, dec0, time, freq, telescope,
                  station_id = 0, beam_norms = np.ones(2)):
    
    # Convert RA and Dec to ITRS coordinates
    dir_itrf = radec_to_xyz(ra, dec, time)
    phase_itrf = radec_to_xyz(ra0, dec0, time)
    
    response = telescope.station_response(time.mjd*3600*24, station_id, freq,
                                          dir_itrf, phase_itrf)
    
    response[0,:] *= beam_norms[0] 
    response[1,:] *= beam_norms[1]
    
    return response
            
    
            

# def run_everybeam(source_catalogue : Source_Catalogue, args : argparse.Namespace,
#                   woden_settings : Woden_Settings):
    
    
#     everybeam_models = ["everybeam_OSKAR"]
    
#     ##Only need to run everybeam if we are using an everybeam model
#     if args.primary_beam in everybeam_models:
        
#         if args.primary_beam == "everybeam_OSKAR":
#             telescope = load_OSKAR_telescope(args.beam_ms_path)
        
    
#         obs_time = Time(args.date, scale='utc')
        
#         for time_step in range(woden_settings.num_time_steps):
        
#             time_current = obs_time + TimeDelta((time_step + 0.5)*woden_settings.time_res, format='sec')
            
#             for source_ind in range(source_catalogue.num_sources):
#                 source = source_catalogue.sources[source_ind]
                
#                 print(source.n_points)
                
#                 # if source.n_points > 0:
#                 #     print(source.n_points)
                
                
        
#     return