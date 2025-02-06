from wodenpy.primary_beam import use_everybeam
from astropy.time import Time, TimeDelta
import numpy as np
from astropy.coordinates import EarthLocation
from astropy import units as u
import everybeam as eb

LATITUDE_LOFAR=52.905329712
LONGITUDE_LOFAR=6.867996528

NUM_COORDS = 10
NUM_DIRS = NUM_COORDS * NUM_COORDS
NUM_TIMES = 2
NUM_FREQS = 2
NUM_STATIONS = 1
D2R = np.pi / 180.0

RA_WIDTH_LOFAR = np.pi*2
DEC_WIDTH_LOFAR = 40.0 * D2R
RA0_LOFAR = 0.0
DEC0_LOFAR = 89.9*D2R
LOW_RA_LOFAR = 0.0
LOW_DEC_LOFAR = 49.5*D2R-DEC_WIDTH_LOFAR

LOW_FREQ_HBA=160e+6
FREQ_INC_HBA=1e+6
TIME_RES=10.0

def make_radec(low_ra, low_dec, ra_width, dec_width, ra0, dec0):
    ras = np.zeros(NUM_DIRS)
    decs = np.zeros(NUM_DIRS)

    ra_inc = (2 * np.pi) / (NUM_COORDS - 1)
    dec_inc = dec_width / (NUM_COORDS - 1)

    for rai in range(NUM_COORDS):
        for deci in range(NUM_COORDS):
            ra = low_ra + rai * ra_inc
            dec = low_dec + deci * dec_inc

            ras[rai * NUM_COORDS + deci] = ra
            decs[rai * NUM_COORDS + deci] = dec

            # print(f"ra: {ra / D2R:.2f}, dec: {dec / D2R:.2f}")

    return ras, decs


def write_C_array(f, name, jones):
    f.write(f"const double {name}[] = {{\n")
        
    for station in range(jones.shape[0]):
        for time in range(jones.shape[1]):
            for freq in range(jones.shape[2]):
                for coord in range(jones.shape[3]):
                    for pol1 in range(2):
                        for pol2 in range(2):
                            f.write(f"    {jones[station, time, freq, coord, pol1, pol2].real}, {jones[station, time, freq, coord, pol1, pol2].imag},\n")
    f.write("};\n\n")

def make_hba_values():
    
    ras, decs = make_radec(LOW_RA_LOFAR, LOW_DEC_LOFAR,
                           RA_WIDTH_LOFAR, DEC_WIDTH_LOFAR,
                           RA0_LOFAR, DEC0_LOFAR)
    
    observing_location = EarthLocation(lat=LATITUDE_LOFAR*u.deg,
                                       lon=LONGITUDE_LOFAR*u.deg, height=0)
    lst_type = 'mean'

    date = "2013-11-18T15:59:59"
    observing_time = Time(date, scale='utc', location=observing_location)
    
    lst = observing_time.sidereal_time(lst_type).value*15.0
    
    times = np.array([observing_time + TimeDelta(i*TIME_RES, format='sec') for i in range(NUM_TIMES)])
    mjd_sec_times = np.array([time.mjd * 86400.0 for time in times])
    
    ms_path = "../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms"
    telescope = use_everybeam.load_LOFAR_telescope(ms_path)
    
    j2000_latitudes = np.array([LATITUDE_LOFAR] * NUM_TIMES)
    j2000_lsts = np.array([lst] * NUM_TIMES)
    
    freqs = LOW_FREQ_HBA + np.arange(NUM_FREQS) * FREQ_INC_HBA
    station_ids = np.array([0] * NUM_STATIONS)
    
    full_accuracy = True
    parallactic_rotate = False
    
    jones = use_everybeam.run_everybeam(ras, decs,
                  RA0_LOFAR, DEC0_LOFAR, j2000_latitudes, j2000_lsts,
                  LATITUDE_LOFAR, LONGITUDE_LOFAR,
                  times, freqs,
                  telescope, station_ids,
                  full_accuracy=full_accuracy,
                  parallactic_rotate=parallactic_rotate,
                  apply_beam_norms=False,
                  reorder_jones=False,
                  element_only=False,
                  eb_rotate=False)
    
    jones_rotate = use_everybeam.run_everybeam(ras, decs,
                  RA0_LOFAR, DEC0_LOFAR, j2000_latitudes, j2000_lsts,
                  LATITUDE_LOFAR, LONGITUDE_LOFAR,
                  times, freqs,
                  telescope, station_ids,
                  full_accuracy=full_accuracy,
                  parallactic_rotate=parallactic_rotate,
                  apply_beam_norms=False,
                  reorder_jones=False,
                  element_only=False,
                  eb_rotate=True)
    
    jones_normed = use_everybeam.run_everybeam(ras, decs,
                  RA0_LOFAR, DEC0_LOFAR, j2000_latitudes, j2000_lsts,
                  LATITUDE_LOFAR, LONGITUDE_LOFAR,
                  times, freqs,
                  telescope, station_ids,
                  full_accuracy=full_accuracy,
                  parallactic_rotate=parallactic_rotate,
                  apply_beam_norms=True,
                  reorder_jones=False,
                  element_only=False,
                  eb_rotate=False)
    
    jones_rotate_normed = use_everybeam.run_everybeam(ras, decs,
                  RA0_LOFAR, DEC0_LOFAR, j2000_latitudes, j2000_lsts,
                  LATITUDE_LOFAR, LONGITUDE_LOFAR,
                  times, freqs,
                  telescope, station_ids,
                  full_accuracy=full_accuracy,
                  parallactic_rotate=parallactic_rotate,
                  apply_beam_norms=True,
                  reorder_jones=False,
                  element_only=False,
                  eb_rotate=True)
    
    
    with open("lofar_jones_values.h", "w") as f:
        f.write("#ifndef LOFAR_JONES_VALUES_H\n")
        f.write("#define LOFAR_JONES_VALUES_H\n\n")
        
        
        write_C_array(f, "mjd_sec_times", jones)
        
        write_C_array(f, "lofar_hba_jones", jones)
        write_C_array(f, "lofar_hba_jones_rotate", jones_rotate)
        write_C_array(f, "lofar_hba_jones_normed", jones_normed)
        write_C_array(f, "lofar_hba_jones_rotate_normed", jones_rotate_normed)
        
        
        f.write("#endif // LOFAR_JONES_VALUES_H\n")

if __name__ == "__main__":
    
    make_hba_values()