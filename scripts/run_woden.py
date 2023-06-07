#!/usr/bin/env python3
"""Wrapper script to run the GPU WODEN code. Author: J.L.B. Line
"""
import numpy as np
import sys
import os
from time import time
from multiprocessing import Pool
import argparse
from queue import Queue
from threading import Thread

from ctypes import POINTER, c_double, c_float
from scipy.special import factorial,eval_hermite

##If we are performing a ctest, this check means we use the code we are
##testing and NOT what has been pip or conda installed
try:
    testdir = os.environ['CMAKE_CURRENT_SOURCE_DIR']
    sys.path.append('{:s}/../../../wodenpy'.format(testdir))
    
    from use_libwoden.woden_settings import create_woden_settings, setup_lsts_and_phase_centre
    from use_libwoden.visibility_set import setup_visi_set, load_visibility_set
    from wodenpy_setup.run_setup import get_parser, check_args, get_code_version
    from use_libwoden.use_libwoden import load_in_woden_library
    from observational.calc_obs import get_uvfits_date_and_position_constants, calc_jdcal
    from skymodel.woden_skymodel import Component_Type_Counter, CompTypes, crop_below_horizon
    from skymodel.read_yaml_skymodel import read_yaml_radec_count_components, read_yaml_skymodel_chunks, convert_source_catalogue_to_ctypes, read_yaml_skymodel_chunks_noctype
    from skymodel.chunk_sky_model import create_skymodel_chunk_map
    from array_layout.create_array_layout import calc_XYZ_diffs, enh2xyz
    from uvfits.wodenpy_uvfits import make_antenna_table, make_baseline_date_arrays, create_uvfits
    from phase_rotate.remove_phase_track import remove_phase_tracking

except KeyError:


    from wodenpy.use_libwoden.woden_settings import create_woden_settings, setup_lsts_and_phase_centre
    from wodenpy.use_libwoden.visibility_set import setup_visi_set, load_visibility_set
    from wodenpy.wodenpy_setup.run_setup import get_parser, check_args, get_code_version
    from wodenpy.use_libwoden.use_libwoden import load_in_woden_library
    from wodenpy.observational.calc_obs import get_uvfits_date_and_position_constants, calc_jdcal
    from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes, crop_below_horizon
    from wodenpy.skymodel.read_yaml_skymodel import read_yaml_radec_count_components, read_yaml_skymodel_chunks, convert_source_catalogue_to_ctypes, read_yaml_skymodel_chunks_noctype
    from wodenpy.skymodel.chunk_sky_model import create_skymodel_chunk_map
    from wodenpy.array_layout.create_array_layout import calc_XYZ_diffs, enh2xyz
    from wodenpy.uvfits.wodenpy_uvfits import make_antenna_table, make_baseline_date_arrays, create_uvfits
    from wodenpy.phase_rotate.remove_phase_track import remove_phase_tracking

##Constants
R2D = 180.0 / np.pi
D2R = np.pi / 180.0
MWA_LAT = -26.703319405555554
MWA_LONG = 116.67081523611111
MWA_HEIGHT = 377.827
VELC = 299792458.0
SOLAR2SIDEREAL = 1.00274

sbf_N = 101
sbf_L = 10001
sbf_c = 5000
sbf_dx = 0.01

def calc_basis_func_1D(x, n, beta=1):
    '''explicitly set beta=1'''
    
    norm = np.sqrt(1)*np.sqrt((2**n*factorial(n)))
    hermite = eval_hermite(n, x)
    gauss = np.exp(-0.5*((x*beta)**2))
    
    return (hermite*gauss) / norm

def create_sbf(precision = "double", sbf_N = 101, sbf_c = 5000, sbf_dx = 0.01):
    
    x_range = np.arange(-sbf_c*sbf_dx, sbf_c*sbf_dx + sbf_dx, sbf_dx)
    sbf_L = len(x_range)
    
    if precision == 'float':
        sbf = (c_float*(sbf_L*sbf_N))()
    else:
        sbf = (c_double*(sbf_L*sbf_N))()
    
    for n in range(sbf_N):
        low_ind = n*sbf_L
        basis = calc_basis_func_1D(x_range, n)
        
        for basis_ind, sbf_ind in enumerate(range(low_ind, low_ind + sbf_L)):
            sbf[sbf_ind] = basis[basis_ind]
    
    return sbf


def woden_thread(the_queue, run_woden, woden_settings, visibility_set, array_layout,
                 sbf):
    while True:
        source_catalogue = the_queue.get(block=True)
        
        if source_catalogue is None:
            return
        
        run_woden(woden_settings, visibility_set, source_catalogue, array_layout,
                  sbf)
        
def read_skymodel_thread(the_queue, chunked_skymodel_maps: list,
                                    max_num_chunks : int,
                                    lsts : np.ndarray,
                                    args : argparse.Namespace,
                                    beamtype : int):
    
    for lower_chunk_iter in range(0, len(chunked_skymodel_maps), max_num_chunks):
        
        chunk_map_subset = chunked_skymodel_maps[lower_chunk_iter:lower_chunk_iter+max_num_chunks]
        
        print(f"Reading chunks skymodel chunks {lower_chunk_iter}:{lower_chunk_iter+max_num_chunks}")
        t_before = time()
    
        chunked_sky_models = read_yaml_skymodel_chunks_noctype(args.cat_filename, chunk_map_subset,
                                                        args.num_freq_channels,
                                                        args.num_time_steps,
                                                        beamtype,
                                                        lsts, args.latitude*D2R)
        
        source_catalogue = convert_source_catalogue_to_ctypes(chunked_sky_models,
                                chunk_map_subset,
                                args.num_freq_channels,
                                args.num_time_steps,
                                beamtype,
                                precision=args.precision)
        the_queue.put(source_catalogue, block=True)
        
        t_after = time()
        print(f"Reading chunks {lower_chunk_iter}:{lower_chunk_iter+max_num_chunks} took {(t_after - t_before)/60.0:.1f} mins")
        # print(f"Skymodel chunks {lower_chunk_iter}:{lower_chunk_iter+max_num_chunks} have been read in")

    ##finally, stick a None into the_queue, this means `c_thead` knows to
    ##stop
    the_queue.put(None)


def main():
    """Run the things you fool"""

    ##Grab the parser and parse some args
    parser = get_parser()
    args = parser.parse_args()

    ##Check that the input arguments make sense
    args = check_args(args)

    ##Find out what git/release version we are using, and where the code lives
    gitlabel = get_code_version()

    lst_deg, gst0_deg, degpdy, ut1utc = get_uvfits_date_and_position_constants(latitude=args.latitude, longitude=args.longitude,
                        height=args.array_height, date=args.date)

    int_jd, float_jd = calc_jdcal(args.date)
    jd_date = int_jd + float_jd

    ##Check the uvfits prepend to make sure we end in .uvfits
    output_uvfits_prepend = args.output_uvfits_prepend
    if output_uvfits_prepend[-7:] == '.uvfits': output_uvfits_prepend = output_uvfits_prepend[-7:]

    if args.dry_run:
        ##User only wants to check whether the arguments would have worked or not
        pass
        
    else:
        ##Depending on what precision was selected by the user, load in the
        ##C/CUDA library and return the `run_woden` function

        run_woden = load_in_woden_library(args.precision)
        
        num_baselines = int(((args.num_antennas - 1)*args.num_antennas) / 2)

        if args.do_autos:
            num_visis = args.num_time_steps*args.num_freq_channels*(num_baselines + args.num_antennas)
        else:
            num_visis = args.num_time_steps*args.num_freq_channels*num_baselines

        ##This essentially does a malloc so we can shove this straight into
        ##the C library I think
        visibility_set = setup_visi_set(num_visis, args.precision)
        
        ##populates a ctype equivalent of woden_settings struct to pass
        ##to the C library
        woden_settings = create_woden_settings(args, jd_date, lst_deg)

        ##fill the lst and mjds fields, precessing if necessary
        lsts = setup_lsts_and_phase_centre(woden_settings)

        ##calculate the array layout
        array_layout = calc_XYZ_diffs(woden_settings, args)
        
        # ##read in and chunk the sky model=======================================
        print("Doing the initial reading/mapping of sky model into chunks")
        t_before = time()
        comp_counter = read_yaml_radec_count_components(args.cat_filename)
        t_after = time()
        print(f"Mapping took {(t_after - t_before)/60.0:.1f} mins")
    
        crop_by_component = True
        if args.sky_crop_sources:
            crop_by_component = False
            
        comp_counter = crop_below_horizon(lst_deg*D2R, args.latitude*D2R, comp_counter,
                                            crop_by_component=crop_by_component)
        
        chunked_skymodel_maps = create_skymodel_chunk_map(comp_counter,
                                            args.chunking_size, woden_settings.num_baselines,
                                            args.num_freq_channels,
                                            args.num_time_steps)
        
        
        sbf = create_sbf(precision=args.precision)
        
        ###---------------------------------------------------------------------
        ### heavy lifting area - here we setup running the sky model reading
        ### and running GPU code at the same time. Means we can limit the
        ### absolute amount of RAM used, and save time doing IO at the same
        ### time as compute
        
        ##TODO this needs to be calculated from actual avaible RAM
        max_num_chunks = 50
        
        ##setup a queue system to run out GPU code at the same time as the
        ##model reading
        the_queue = Queue(maxsize=1)
        
        t1 = Thread(target = read_skymodel_thread,
                    args =(the_queue, chunked_skymodel_maps,
                                    max_num_chunks, lsts, args,
                                    woden_settings.beamtype), daemon=True)
        
        t2 = Thread(target = woden_thread, args =(the_queue, run_woden,
                                                  woden_settings, visibility_set,
                                                  array_layout, sbf), daemon=True)
        
        t1.start()
        t2.start()
        t1.join()
        t2.join()
        
        ### we've now calculated all the visibilities
        ###---------------------------------------------------------------------
        

        ##I think we want to X,Y,Z to be in the current frame for writing
        ##out to the uvfits, so calculate again
        X,Y,Z = enh2xyz(args.east, args.north, args.height, args.latitude*D2R)
        ##X,Y,Z are stored in a 2D array in units of seconds in the uvfits file
        XYZ_array = np.empty((args.num_antennas,3))
        XYZ_array[:,0] = X
        XYZ_array[:,1] = Y
        XYZ_array[:,2] = Z

        ##If we were to grab them out of the code used in the simulation,
        ##need to do this
        # expec_num = woden_settings.num_time_steps*args.num_antennas
        # XYZ_array[:,0] = np.ctypeslib.as_array(array_layout.ant_X, shape=(expec_num, ))[:args.num_antennas]
        # XYZ_array[:,1] = np.ctypeslib.as_array(array_layout.ant_Y, shape=(expec_num, ))[:args.num_antennas]
        # XYZ_array[:,2] = np.ctypeslib.as_array(array_layout.ant_Z, shape=(expec_num, ))[:args.num_antennas]

        ##Get the central frequency channels, used in the uvfits header
        central_freq_chan = int(np.floor(args.num_freq_channels / 2.0))
        ##Useful number
        num_baselines = int(((args.num_antennas - 1)*args.num_antennas) / 2)

        ##Loop over coarse frequency band
        for band in args.band_nums:
            print('Converting binary to uvfits band',band)

            output_uvfits_name = output_uvfits_prepend + '_band%02d.uvfits' %band

            band_low_freq = args.lowest_channel_freq + (band - 1)*args.coarse_band_width
            central_freq_chan_value = band_low_freq + central_freq_chan*args.freq_res


            uus,vvs,wws,v_container = load_visibility_set(visibility_set=visibility_set,
                                                num_baselines=num_baselines,
                                                num_freq_channels=args.num_freq_channels,
                                                num_time_steps=args.num_time_steps,
                                                precision=args.precision,
                                                do_autos=args.do_autos,
                                                num_ants=args.num_antennas)

            if args.remove_phase_tracking:
                frequencies = band_low_freq + np.arange(args.num_freq_channels)*args.freq_res

                v_container = remove_phase_tracking(frequencies=frequencies,
                                          wws_seconds=wws,
                                          num_time_steps=args.num_time_steps,
                                          v_container=v_container,
                                          num_baselines=num_baselines)

            hdu_ant = make_antenna_table(XYZ_array=XYZ_array,telescope_name=args.telescope_name,
                          num_antennas=args.num_antennas, freq_cent=central_freq_chan_value,
                          date=args.date, gst0_deg=gst0_deg, degpdy=degpdy,
                          ut1utc=ut1utc, longitude=args.longitude, latitude=args.latitude,
                          array_height=args.array_height,
                          ant_names=args.ant_names)

            baselines_array, date_array = make_baseline_date_arrays(args.num_antennas,
                                          args.date, args.num_time_steps, args.time_res,
                                          do_autos=args.do_autos)

            create_uvfits(v_container=v_container, freq_cent=central_freq_chan_value,
                          ra_point=args.ra0, dec_point=args.dec0,
                          output_uvfits_name=output_uvfits_name,
                          uu=uus, vv=vvs, ww=wws, baselines_array=baselines_array,
                          date_array=date_array,
                          central_freq_chan=central_freq_chan,
                          ch_width=args.freq_res, int_jd=int_jd,
                          hdu_ant=hdu_ant, gitlabel=gitlabel,
                          longitude=args.longitude, latitude=args.latitude,
                          array_height=args.array_height,
                          telescope_name=args.telescope_name)


if __name__ == "__main__":
    main()
    
