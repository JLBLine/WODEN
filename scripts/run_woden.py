#!/usr/bin/env python3
"""Wrapper script to generate json input files for, and to run,
the GPU WODEN code. Author: J.L.B. Line
"""
from wodenpy.woden_lib import *
from wodenpy.wodenpy_setup.run_setup import *

##Constants
R2D = 180.0 / np.pi
D2R = np.pi / 180.0
MWA_LAT = -26.703319405555554
MWA_LONG = 116.67081523611111
MWA_HEIGHT = 377.827
VELC = 299792458.0
SOLAR2SIDEREAL = 1.00274
    
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

    ##Write json file
    json_band_str =  '-'.join(map(str, args.band_nums))
    json_name = 'run_woden_{:s}.json'.format(json_band_str)

    write_json(json_name=json_name, jd_date=jd_date, lst=lst_deg, args=args)

    ##Check the uvfits prepend to make sure we end in .uvfits
    output_uvfits_prepend = args.output_uvfits_prepend
    if output_uvfits_prepend[-7:] == '.uvfits': output_uvfits_prepend = output_uvfits_prepend[-7:]

    if args.dry_run:
        ##User only wants to check whether the arguments would have worked or not
        ##TODO stick all the arguments that would have been set into the log?
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
        visibility_set = setup_visi_set(num_visis)
        
        ##This calls the WODEN main function, and populates visibility set
        ##with the outputs that we want
        run_woden(json_name.encode('utf-8'), visibility_set)

        # command(f'{WODEN_DIR}/{woden_lib} {json_name}')

        X,Y,Z = enh2xyz(args.east, args.north, args.height, args.latitude*D2R)
        ##X,Y,Z are stored in a 2D array in units of seconds in the uvfits file
        XYZ_array = np.empty((args.num_antennas,3))
        XYZ_array[:,0] = X
        XYZ_array[:,1] = Y
        XYZ_array[:,2] = Z

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


            uus,vvs,wws,v_container = load_data(visibility_set=visibility_set,
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
                          telescope_name=args.telescope_name,)


        ##Tidy up or not
        if args.no_tidy:
            pass
        else:
            command("rm {:s}".format(json_name))
            ##if we generated a text file containing the array layout
            ##from the metafits, delete it now
            # if args.array_layout == 'from_the_metafits':
            command("rm WODEN_array_layout_band{:s}.txt".format(json_band_str))
            #     command("rm {:s}".format("WODEN_array_layout.txt"))

    

if __name__ == "__main__":
    main()
    