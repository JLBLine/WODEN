import numpy as np
from casacore.tables import table
from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation
import astropy.units as u
import erfa
from wodenpy.primary_beam.use_everybeam import run_lofar_beam, run_mwa_beam, run_everybeam

LATITUDE_LOFAR=52.905329712
LONGITUDE_LOFAR=6.867996528

LONGITUDE_MWA=116.67081523611111
LATITUDE_MWA=-26.703319405555554


def print_variables(ms_path, prepend, latitude, longitude):
    
    with table(ms_path + "::FIELD") as field_subtable:
        phase_dir = field_subtable.getcol('PHASE_DIR')
    
    
    
    
    RA0 = np.squeeze(phase_dir)[0]
    if RA0 < 0:
        RA0 += 2*np.pi
    DEC0 = np.squeeze(phase_dir)[1]


    with table(ms_path) as ms:
        time_col = ms.getcol("TIME")
        mjd_secs = time_col[0]
    
    date = Time(mjd_secs*u.s, format='mjd')
    first_date = date.datetime.strftime('%Y-%m-%dT%H:%M:%S')
    # print(first_date, date.mjd*24*3600)
    observing_location = EarthLocation(lat=latitude*u.deg,
                                    lon=LONGITUDE_LOFAR*u.deg, height=0)


    time_inc = 3*3600

    times = [Time(date, scale='utc', location=observing_location) + TimeDelta(time_inc*i, format='sec') for i in range(2)]

    lsts = [time.sidereal_time('mean').rad for time in times]
    mjds = [time.mjd for time in times]
    
    
    if "MWA" in ms_path:
        RA0 = lsts[0]
        DEC0 = latitude

    azs, els = erfa.hd2ae(lsts - RA0, DEC0, latitude)
    
    zas = np.pi/2 - els
    para_angles = erfa.hd2pa(lsts - RA0, DEC0, latitude)
    
    print(f"#define {prepend}_RA0 {RA0:.8f}")
    print(f"#define {prepend}_DEC0 {DEC0:.8f}")
    print(f"double {prepend}_mjds[] ={{{mjds[0]:.8f}, {mjds[1]:.8f}}};")
    print(f"double {prepend}_lsts[] ={{{lsts[0]:.8f}, {lsts[1]:.8f}}};")
    print(f"double {prepend}_azs[] ={{{azs[0]:.8f}, {azs[1]:.8f}}};")
    print(f"double {prepend}_zas[] ={{{zas[0]:.8f}, {zas[1]:.8f}}};")
    print(f"double {prepend}_para_angles[] ={{{para_angles[0]:.8f}, {para_angles[1]:.8f}}};")
    
    
    
    
    element_response_model = 'MWA'
    coeff_path = '/home/jack-line/software/mwa_beam_files/mwa_full_embedded_element_pattern.h5'
    apply_beam_norms = False
    rotate = True
    element_only = False
    
    mjd_sec_times  = np.array([time.mjd*24*3600 for time in times])
    
    freqs = [1.2e+8]
    station_ids = [0]
    
    new_jones = run_mwa_beam(ms_path, element_response_model,
                               coeff_path, 
                               station_ids,
                               [RA0], [DEC0],
                               mjd_sec_times,
                               lsts, [latitude]*2,
                               freqs,
                               apply_beam_norms=apply_beam_norms,
                               parallactic_rotate=rotate,
                               iau_order=True, element_only=element_only)
    
    # print(new_jones)
    
    


if __name__ == "__main__":
    ms_path = "../../../test_installation/everybeam/MWA-single-timeslot.ms"
    print_variables(ms_path, "MWA", np.radians(LATITUDE_MWA), np.radians(LONGITUDE_MWA))