import numpy as np
import ctypes
from wodenpy.primary_beam import use_everybeam
from astropy.time import Time, TimeDelta
import numpy as np
from astropy.coordinates import EarthLocation
from astropy import units as u
from wodenpy.use_libwoden.skymodel_structs import c_double_complex
import numpy.testing as npt
from time import time
from line_profiler import LineProfiler, profile

LATITUDE_LOFAR=52.905329712
LONGITUDE_LOFAR=6.867996528

NUM_COORDS = 200
NUM_DIRS = NUM_COORDS * NUM_COORDS
NUM_TIMES = 10
NUM_FREQS = 10
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

@profile
def run_lofar_beam(ms_path : str, element_response_model : bool,
                   coeff_path : str,
                   num_stations : int, station_idxs : np.ndarray,
                   num_dirs : int, ra0 : float, dec0 : float,
                   ras : np.ndarray, decs : np.ndarray,
                   num_times : int, mjd_sec_times : np.ndarray,
                   num_freqs : int, freqs : np.ndarray,
                   apply_beam_norms : bool, rotate : bool):
    
    # lib_path = importlib_resources.files(wodenpy).joinpath(f"libwoden_{args.precision}.so")
    
    woden_lib = ctypes.cdll.LoadLibrary("/home/jline/software/WODEN_dev/build/libuse_everybeam.so")
    
    load_and_run_lofar_beam = woden_lib.load_and_run_lofar_beam
    
    station_idxs_ctypes = station_idxs.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    ras_ctypes = ras.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    decs_ctypes = decs.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    mjd_sec_times_ctypes = mjd_sec_times.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    freqs_ctypes = freqs.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    
    ms_path_ctypes = ctypes.c_char_p(ms_path.encode('utf-8'))
    element_response_model_ctypes = ctypes.c_char_p(element_response_model.encode('utf-8'))
    coeff_path_ctypes = ctypes.c_char_p(coeff_path.encode('utf-8'))
    
    jones = ((num_stations*num_times*num_freqs*num_dirs*4)*c_double_complex)()
    
    load_and_run_lofar_beam.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
                                        ctypes.c_int, ctypes.POINTER(ctypes.c_int),
                                        ctypes.c_int, ctypes.c_double, ctypes.c_double,
                                        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                        ctypes.c_int, ctypes.POINTER(ctypes.c_double),
                                        ctypes.c_int, ctypes.POINTER(ctypes.c_double),
                                        ctypes.c_bool, ctypes.c_bool,
                                        ctypes.POINTER(c_double_complex)]
    
    load_and_run_lofar_beam(ms_path_ctypes,
                            element_response_model_ctypes,
                            coeff_path_ctypes,
                            num_stations, station_idxs_ctypes,
                            num_dirs,
                            ra0, dec0,
                            ras_ctypes, decs_ctypes,
                            num_times, mjd_sec_times_ctypes,
                            num_freqs, freqs_ctypes,
                            apply_beam_norms, rotate,
                            jones)
    
    # print(jones)
    
    jones_py = np.ctypeslib.as_array(jones, shape=(num_stations*num_times*num_freqs*num_dirs*4))
    jones_py = jones_py['real'] + 1j*jones_py['imag']
    
    
    return jones_py
    
@profile
def main():
    ras, decs = make_radec(LOW_RA_LOFAR, LOW_DEC_LOFAR,
                           RA_WIDTH_LOFAR, DEC_WIDTH_LOFAR,
                           RA0_LOFAR, DEC0_LOFAR)
    
    observing_location = EarthLocation(lat=LATITUDE_LOFAR*u.deg,
                                       lon=LONGITUDE_LOFAR*u.deg, height=0)
    lst_type = 'mean'

    date = "2013-11-18T15:59:59"
    observing_time = Time(date, scale='utc', location=observing_location)
    
    lst = observing_time.sidereal_time(lst_type).value*15.0
    
    ms_path = "../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms"
    telescope = use_everybeam.load_LOFAR_telescope(ms_path)
    
    j2000_latitudes = np.array([LATITUDE_LOFAR] * NUM_TIMES)
    j2000_lsts = np.array([lst] * NUM_TIMES)
    
    
    
    times = np.array([observing_time + TimeDelta(i*TIME_RES, format='sec') for i in range(NUM_TIMES)])
    
    freqs = LOW_FREQ_HBA + np.arange(NUM_FREQS) * FREQ_INC_HBA
    station_ids = np.array([0] * NUM_STATIONS)
    
    # full_accuracy = True
    
    element_response_model = 'hamaker'
    coeff_path = ''
    apply_beam_norms = False
    # apply_beam_norms = False
    # rotate = False
    rotate = True
    
    mjd_sec_times = np.array([time.mjd * 86400.0 for time in times])
    
    # print(mjd_sec_times)
    
    start = time()
    
    new_jones = run_lofar_beam(ms_path, element_response_model,
                               coeff_path, 
                               NUM_STATIONS, station_ids,
                               NUM_DIRS, RA0_LOFAR, DEC0_LOFAR,
                               ras, decs,
                               NUM_TIMES, mjd_sec_times,
                               NUM_FREQS, freqs,
                               apply_beam_norms, rotate)
    
    print(f"New time taken: {(time() - start):1f} secs")
    
    start = time()
    
    telescope = use_everybeam.load_LOFAR_telescope(ms_path)
    
    old_jones = use_everybeam.run_everybeam(ras, decs,
                  RA0_LOFAR, DEC0_LOFAR, j2000_latitudes, j2000_lsts,
                  LATITUDE_LOFAR, LONGITUDE_LOFAR,
                  times, freqs,
                  telescope, station_ids,
                  full_accuracy=True,
                  parallactic_rotate=False,
                  apply_beam_norms=apply_beam_norms,
                  reorder_jones=False,
                  element_only=False,
                  eb_rotate=rotate)
    
    old_jones = old_jones.flatten()
    
    print(f"Old time taken: {(time() - start):1f} secs")
    
    npt.assert_allclose(np.imag(new_jones), np.imag(old_jones), atol=5e-3)
    npt.assert_allclose(np.real(new_jones), np.real(old_jones), atol=5e-3)
    
    
if __name__ == "__main__":
    main()
    
    