#!/usr/bin/env python3
"""Wrapper script to profile WODEN code.
"""
import sys
import yappi

sys.path.append('../../scripts')
import run_woden as rw


if __name__ == "__main__":
    num_times = 4
    num_freqs = 4
    
    ra0 = 0.0
    dec0 = -26.7
    date = "2024-07-21T20:13:00"
    profile_cat = "profiling_source.fits"

    args = []
    args.append("--band_nums=1")
    args.append("--lowest_channel_freq=160e+6")
    args.append(f"--num_freq_channels={num_freqs}")
    args.append("--freq_res=80e+3")
    args.append(f"--num_time_steps={num_times}")
    args.append("--time_res=2")
    args.append(f"--ra0={ra0}")
    args.append(f"--dec0={dec0}")
    args.append(f"--date={date}")
    args.append(f"--cat_filename={profile_cat}")
    args.append("--output_uvfits_prepend=profile_run_woden")
    args.append(f'--beam_ms_path=create_OSKAR-MWA_ms/OSKAR-MWA-flags-layout.ms')
    args.append("--primary_beam=everybeam_OSKAR")
    
    yappi.set_clock_type("wall") # Use set_clock_type("cpu") for CPU time
    # yappi.start()


    with yappi.run():
        ##Run the actual WODEN code
        rw.main(args)

    yappi.stop()


    print('\nWODENPY MODULE STATS------------------------------------------------------------')
    stats = yappi.get_func_stats(filter_callback=lambda x: 'wodenpy' in x.module
                                ).print_all()
    print('--------------------------------------------------------------------------------\n')


    print('\nGPU CALL THREAD-----------------------------------------------------------------')
    threads = yappi.get_thread_stats()
    ##Second thread calls the GPU code. Search within that for the function called
    ##`_run_run_woden`
    stats = yappi.get_func_stats(ctx_id=2,
                                    filter_callback=lambda x: '_run_run_woden' in x.name
                                ).print_all()
    
    print('--------------------------------------------------------------------------------\n')