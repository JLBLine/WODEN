import numpy as np
import ctypes
from wodenpy.primary_beam import use_everybeam
from astropy.time import Time, TimeDelta
import numpy as np
from astropy.coordinates import EarthLocation
from astropy import units as u
import wodenpy
from wodenpy.use_libwoden.skymodel_structs import c_double_complex
import numpy.testing as npt
from time import time
from line_profiler import LineProfiler, profile
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from casacore.tables import table, taql
from astropy.coordinates import SkyCoord, FK5
from astropy.time import Time
import importlib_resources
from wodenpy.primary_beam.use_everybeam import run_lofar_beam, run_mwa_beam, run_everybeam
import erfa
import logging
import matplotlib

logging.getLogger("matplotlib").setLevel(logging.WARNING)

D2R = np.pi / 180.0

# LATITUDE_LOFAR=52.905329712
# LONGITUDE_LOFAR=6.867996528

# RA0_LOFAR = 120.0*D2R
# DEC0_LOFAR = 45.0*D2R

LATITUDE_LOFAR=-26.703319405555554
LONGITUDE_LOFAR=116.67081523611111

RA0_LOFAR = 153.1889998*D2R
DEC0_LOFAR = -26.703319405555554*D2R

NUM_COORDS = 50
NUM_DIRS = NUM_COORDS * NUM_COORDS
NUM_TIMES = 1
NUM_FREQS = 1
NUM_STATIONS = 1



# RA_WIDTH_LOFAR = np.pi*2
# DEC_WIDTH_LOFAR = 40.0 * D2R
# RA0_LOFAR = 123.4002825*D2R
# DEC0_LOFAR = 48.21738361*D2R

# RA0_LOFAR = 0.0*D2R
# DEC0_LOFAR = 90.0*D2R



# LOW_RA_LOFAR = 0.0
# LOW_DEC_LOFAR = 49.5*D2R-DEC_WIDTH_LOFAR

LOW_FREQ_HBA=160e+6
FREQ_INC_HBA=50e+6
TIME_RES=3600.0*1

def plot_beams(jones, wcs, title, observing_time,
               num_dirs=NUM_DIRS, num_stations=NUM_STATIONS,
               num_times=NUM_TIMES, num_freqs=NUM_FREQS):
    
    # jones = jones.reshape(num_stations, num_times, num_freqs, num_dirs, 2, 2)
    
    fig, axs = plt.subplots(num_times, num_freqs, figsize=(10, 10),
                            subplot_kw={'projection': wcs})
    
    
    # for station in range(num_stations):
    si = 0
    for ti in range(num_times):
        for fi in range(num_freqs):
            # jones_low = 4*(si*num_times*num_freqs*num_dirs + ti*num_freqs*num_dirs + fi*num_dirs)
            # xx_samps = np.arange(0, 4*NUM_DIRS, 4)
            # subset = jones[jones_low:jones_low+4*NUM_DIRS][xx_samps]
            
            # subset = jones[si, ti, fi, :, 0,0]
            subset = jones[si, ti, fi, :, 1,1]
            
            subset.shape = (NUM_COORDS, NUM_COORDS)
            
            if num_times == 1 and num_freqs == 1:
                ax = axs
            elif num_times == 1:
                ax = axs[fi]
            elif num_freqs == 1:
                ax = axs[ti]
            else:
                ax = axs[ti, fi]
            im = ax.imshow(np.abs(subset), origin='lower')
            
            ax.set_title(f"Time: {ti}, Freq: {fi}")
            
            ax.grid()
            
            plt.colorbar(im, ax=ax)
            
            x_beam, y_beam = wcs.all_world2pix([RA0_LOFAR/D2R], [DEC0_LOFAR/D2R], 0)
            ax.plot(x_beam, y_beam, 'C1o', mfc='none')
            
            # x_beam, y_beam = wcs.all_world2pix(coord_precess.ra.value, coord_precess.dec.value, 0)
            # ax.plot(x_beam, y_beam, 'cs', mfc='none')
    plt.tight_layout()
    
    fig.savefig(title, bbox_inches='tight')
    
def create_filtered_ms(ms_path, new_ms_path, ra0, dec0):
    # Step 1: Open original MS and apply selection (first time step & first channel)
    with table(ms_path, readonly=True) as ms:
        time_col = ms.getcol("TIME")
        ddid_col = ms.getcol("DATA_DESC_ID")

        first_time = time_col[0]  # First timestamp
        first_ddid = ddid_col[0]  # First frequency channel

        # Use TaQL (Table Query Language) to select the subset efficiently
        query = f"SELECT * FROM {ms_path} WHERE TIME = {first_time} AND DATA_DESC_ID = {first_ddid}"
        filtered_ms = taql(query)
        
        filtered_ms.copy(new_ms_path, deep=True)
        
        filtered_ms.close()
        
    with table(new_ms_path+'::FIELD', readonly=False) as field_table:
        
        delay_dir = field_table.getcol('DELAY_DIR')
        print(delay_dir.shape)
        
        field_table.putcol('DELAY_DIR', np.array([[[ra0, dec0]]]))
        # field_table.putcol('LOFAR_TILE_BEAM_DIR', np.array([[[ra0, dec0]]]))
        # field_table.putcol('REFERENCE_DIR', np.array([[[ra0, dec0]]]))
    
@profile
def main():
    
    ra_inc = 100 / (NUM_COORDS - 1)
    dec_inc = 100 / (NUM_COORDS - 1)
    
    
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [NUM_COORDS // 2, NUM_COORDS // 2]
    wcs.wcs.cdelt = np.array([-ra_inc, dec_inc])
    wcs.wcs.crval = [RA0_LOFAR / D2R, DEC0_LOFAR / D2R]
    wcs.wcs.ctype = ["RA---SIN", "DEC--SIN"]
    
    xrange, yrange = np.meshgrid(np.arange(NUM_COORDS), np.arange(NUM_COORDS))
    xrange = xrange.flatten()
    yrange = yrange.flatten()
    
    ras, decs = wcs.all_pix2world(xrange, yrange, 0)
    ras *= D2R
    decs *= D2R

    
    observing_location = EarthLocation(lat=LATITUDE_LOFAR*u.deg,
                                       lon=LONGITUDE_LOFAR*u.deg, height=0)
    lst_type = 'mean'
    
    # ms_path = "/home/jack-line/software/WODEN_dev/test_installation/everybeam/LOFAR_HBA_MOCK.ms"
    ms_path = '/home/jack-line/software/WODEN_dev/test_installation/everybeam/MWA-single-timeslot.ms'
    
    
    new_ms_path = 'filtered_measurementset.ms'
    # create_filtered_ms(ms_path, new_ms_path, RA0_LOFAR, DEC0_LOFAR)
    
    # with table(ms_path) as ms:
            
    #         first_time_mjd = ms.getcol("TIME_CENTROID")[0]
    #         time_res = ms.getcol("INTERVAL")[0]
            
    #         times = ms.getcol("TIME")
    #         num_time_steps = np.unique(times).size
            
    #         date = Time((first_time_mjd - time_res/2.0)*u.s, format='mjd')
    #         date = date.datetime.strftime('%Y-%m-%dT%H:%M:%S')
            
    # print(date)
    
    # date = "2024-07-21T20:13:00"
    date = "2013-05-16T10:48:48"
    observing_time = Time(date, scale='utc', location=observing_location)
    # print('DIS', observing_time.mjd*86400.0)
    
    lst = observing_time.sidereal_time(lst_type).value*15.0
    
    
    
    
    j2000_latitudes = np.array([LATITUDE_LOFAR*D2R] * NUM_TIMES)
    # j2000_lsts = np.array([lst] * NUM_TIMES)
    
    
    
    times = np.array([observing_time + TimeDelta(i*TIME_RES, format='sec') for i in range(NUM_TIMES)])
    
    j2000_lsts = np.radians([time.sidereal_time(lst_type).value*15.0 for time in times])
    
    freqs = LOW_FREQ_HBA + np.arange(NUM_FREQS) * FREQ_INC_HBA
    station_ids = np.array([0] * NUM_STATIONS)
    
    # full_accuracy = True
    
    # element_response_model = 'hamaker'
    # coeff_path = ''
    # apply_beam_norms = False
    # rotate = True
    # element_only = False
    
    element_response_model = 'MWA'
    coeff_path = '/home/jack-line/software/mwa_beam_files/mwa_full_embedded_element_pattern.h5'
    apply_beam_norms = False
    rotate = False
    element_only = False
    
    mjd_sec_times = np.array([time.mjd * 86400.0 for time in times])
    
    # print(mjd_sec_times)
    
    # start = time()
    
    # # new_jones = run_lofar_beam(new_ms_path, element_response_model,
    # #                            coeff_path, 
    # #                            NUM_STATIONS, station_ids,
    # #                            NUM_DIRS, RA0_LOFAR, DEC0_LOFAR,
    # #                            ras, decs,
    # #                            NUM_TIMES, mjd_sec_times,
    # #                            NUM_FREQS, freqs,
    # #                            apply_beam_norms, rotate, element_only)
    
    # # new_jones = run_mwa_beam(ms_path, element_response_model,
    # #                            coeff_path, 
    # #                            station_ids,
    # #                            RA0_LOFAR, DEC0_LOFAR,
    # #                            ras, decs,
    # #                            mjd_sec_times,
    # #                            j2000_lsts, j2000_latitudes,
    # #                            freqs,
    # #                            apply_beam_norms, rotate, element_only)
    
    # new_jones = run_everybeam(ms_path, coeff_path,
    #               ras, decs,
    #               RA0_LOFAR, DEC0_LOFAR, j2000_latitudes, j2000_lsts,
    #               LATITUDE_LOFAR*D2R, LONGITUDE_LOFAR*D2R,
    #               times, freqs,
    #               station_ids,
    #               element_response_model=element_response_model,
    #               apply_beam_norms=apply_beam_norms,
    #               parallactic_rotate=True,
    #               reorder_jones=False,
    #               element_only=element_only,
    #               eb_rotate=rotate)
    
    # print(f"New time taken: {(time() - start):1f} secs")
    
    # plot_beams(new_jones, wcs, "new_beam_on_sky.png", observing_time)
    
    start = time()
    
    # telescope = use_everybeam.load_LOFAR_telescope(ms_path)
    
    # old_jones = use_everybeam.run_everybeam(ras, decs,
    #               RA0_LOFAR, DEC0_LOFAR, j2000_latitudes, j2000_lsts,
    #               LATITUDE_LOFAR*D2R, LONGITUDE_LOFAR*D2R,
    #               times, freqs,
    #               telescope, station_ids,
    #               full_accuracy=True,
    #               parallactic_rotate=False,
    #               apply_beam_norms=apply_beam_norms,
    #               reorder_jones=False,
    #               element_only=element_only,
    #               eb_rotate=rotate)
    
    num_threads=1
    old_jones = use_everybeam.run_everybeam_over_threads(num_threads,
                  ms_path, coeff_path,
                  ras, decs,
                  RA0_LOFAR, DEC0_LOFAR, j2000_latitudes, j2000_lsts,
                  LATITUDE_LOFAR*D2R, LONGITUDE_LOFAR*D2R,
                  times, freqs,
                  station_ids,
                  element_response_model=element_response_model,
                  full_accuracy=True,
                  parallactic_rotate=False,
                  apply_beam_norms=apply_beam_norms,
                  reorder_jones=True,
                  element_only=element_only,
                  eb_rotate=rotate)
    
    # old_jones = old_jones.flatten()
    
    print(f"New threaded time taken: {(time() - start):1f} secs")
    
    all_para_angles = np.empty(NUM_DIRS*NUM_TIMES)
    
    for ind, lst in enumerate(j2000_lsts):
        has = lst - ras
        para_angles = erfa.hd2pa(has, decs, LATITUDE_LOFAR*D2R)
        all_para_angles[ind*NUM_DIRS:(ind+1)*NUM_DIRS] = para_angles
    
    
    rot_matrix = np.empty((NUM_DIRS*NUM_TIMES, 2,2))
    
    rot_matrix[:,0,0] = -np.sin(-all_para_angles)
    rot_matrix[:,0,1] = np.cos(-all_para_angles)
    rot_matrix[:,1,0] = np.cos(-all_para_angles)
    rot_matrix[:,1,1] = np.sin(-all_para_angles)
    old_jones_rotated = np.einsum('klm,kmn->kln', old_jones[0, 0, 0, :, :, :], rot_matrix)
    
    
    buffer = old_jones[0, 0, 0, 1250, :, :]
    rot_mat = rot_matrix[1250, :, :]
    rotated = old_jones_rotated[1250, :, :]
    
    print("THIS THING", buffer[0,0] * rot_mat[0,1] + buffer[0,1] * rot_mat[1,1])
    
    print(f"buffer: {np.real(buffer[0,0]):.5f} {np.imag(buffer[0,0]):.5f}, {np.real(buffer[0,1]):.5f} {np.imag(buffer[0,1]):.5f}, "
          f"{np.real(buffer[1,0]):.5f} {np.imag(buffer[1,0]):.5f}, {np.real(buffer[1,1]):.5f} {np.imag(buffer[1,1]):.5f}")
    print(f"rot_mat: {np.real(rot_mat[0,0]):.5f} {np.imag(rot_mat[0,0]):.5f}, {np.real(rot_mat[0,1]):.5f} {np.imag(rot_mat[0,1]):.5f}, "
            f"{np.real(rot_mat[1,0]):.5f} {np.imag(rot_mat[1,0]):.5f}, {np.real(rot_mat[1,1]):.5f} {np.imag(rot_mat[1,1]):.5f}")
    print(f"rotated: {np.real(rotated[0,0]):.5f} {np.imag(rotated[0,0]):.5f}, {np.real(rotated[0,1]):.5f} {np.imag(rotated[0,1]):.5f}, "
            f"{np.real(rotated[1,0]):.5f} {np.imag(rotated[1,0]):.5f}, {np.real(rotated[1,1]):.5f} {np.imag(rotated[1,1]):.5f}")
    
    
    old_jones[0, 0, 0, :, :, :] = old_jones_rotated
    
    
    # rot_matrix = np.array([[-np.sin(-all_para_angles), np.cos(-all_para_angles)],
    #                    [np.cos(-all_para_angles), np.sin(-all_para_angles)]])
    # semi_jones = np.array([[old_jones[0, 0, 0, :, 0, 0], old_jones[0, 0, 0, :, 0, 1]],
    #                        [old_jones[0, 0, 0, :, 1, 0], old_jones[0, 0, 0, :, 1, 1]]])
    # old_jones_rotated = -np.einsum('ijk,jlk->ilk', semi_jones, rot_matrix)

    
    # old_jones[0, 0, 0, :, 0, 0] = old_jones_rotated[0,0,:]
    # old_jones[0, 0, 0, :, 0, 1] = old_jones_rotated[0,1,:]
    # old_jones[0, 0, 0, :, 1, 0] = old_jones_rotated[1,0,:]
    # old_jones[0, 0, 0, :, 1, 1] = old_jones_rotated[1,1,:]
    
    # # plot_beams(old_jones, wcs, "old_beam_on_sky.png", observing_time)
    plot_beams(old_jones, wcs, "new_beam_threaded.png", observing_time)
    
    # # npt.assert_allclose(np.imag(new_jones), np.imag(old_jones), atol=5e-3)
    # # npt.assert_allclose(np.real(new_jones), np.real(old_jones), atol=5e-3)
    
    
if __name__ == "__main__":
    main()
    
    