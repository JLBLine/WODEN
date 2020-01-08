#!/usr/bin/env python
'''Wrapper script to generate json input files for, and to run,
the GPU WODEN code.
'''
from __future__ import print_function
from astropy.io import fits
from numpy import *
from struct import unpack
from subprocess import call, check_output
from jdcal import gcal2jd
from os import environ

gitlabel = check_output(["git", "describe", "--always"]).strip()
# print(check_output(["git", "describe", "--always"]))

R2D = 180.0 / pi
D2R = pi / 180.0
MWA_LAT = -26.7033194444
VELC = 299792458.0
SOLAR2SIDEREAL = 1.00274

WODEN_DIR = environ['WODEN_DIR']

def command(cmd):
    call(cmd,shell=True)
    # print(cmd)

def calc_jdcal(date):
    '''Takes a string format date and returns julian date in
    the two chunks a uvfits file likes'''

    dmy, hms = date.split()

    day,month,year = map(int,dmy.split('-'))
    hour,mins,secs = map(float,hms.split(':'))

    ##For some reason jdcal gives you the date in two pieces
    ##Gives you the time up until midnight of the day
    jd1,jd2 = gcal2jd(year,month,day)
    jd3 = (hour + (mins / 60.0) + (secs / 3600.0)) / 24.0

    jd = jd1 + jd2 + jd3

    jd_day = jd1 + floor(jd2)
    jd_fraction = (jd2 - floor(jd2)) + jd3

    ##The header of the uvdata file takes the integer, and
    ##then the fraction goes into the data array for PTYPE5
    return jd_day, jd_fraction

def create_uvfits(v_container=None,freq_cent=None, ra_point=None, dec_point=None,
                  output_uvfits_name=None,uu=None,vv=None,ww=None,
                  baselines_array=None,date_array=None,date=None,
                  central_freq_chan=None,ch_width=None,template_uvfits=None,
                  int_jd=None):
    '''Takes visibility data and writes out a uvfits format file'''
    ##TODO replace all of this with an interface with pyuvdata

    uvparnames = ['UU','VV','WW','BASELINE','DATE']
    parvals = [uu,vv,ww,baselines_array,date_array]

    uvhdu = fits.GroupData(v_container,parnames=uvparnames,pardata=parvals,bitpix=-32)
    uvhdu = fits.GroupsHDU(uvhdu)

    ###uvfits standards
    uvhdu.header['CTYPE2'] = 'COMPLEX '
    uvhdu.header['CRVAL2'] = 1.0
    uvhdu.header['CRPIX2'] = 1.0
    uvhdu.header['CDELT2'] = 1.0

    ##This means it's linearly polarised
    uvhdu.header['CTYPE3'] = 'STOKES '
    uvhdu.header['CRVAL3'] = -5.0
    uvhdu.header['CRPIX3'] =  1.0
    uvhdu.header['CDELT3'] = -1.0

    uvhdu.header['CTYPE4'] = 'FREQ'
    uvhdu.header['CRVAL4'] = freq_cent  ##Middle pixel value in Hz
    uvhdu.header['CRPIX4'] = int(central_freq_chan) + 1 ##Middle pixel number
    uvhdu.header['CDELT4'] = ch_width

    uvhdu.header['CTYPE5'] = template_uvfits[0].header['CTYPE5']
    uvhdu.header['CRVAL5'] = template_uvfits[0].header['CRVAL5']
    uvhdu.header['CRPIX5'] = template_uvfits[0].header['CRPIX5']
    uvhdu.header['CDELT5'] = template_uvfits[0].header['CDELT5']

    uvhdu.header['CTYPE6'] = template_uvfits[0].header['CTYPE6']
    uvhdu.header['CRVAL6'] = ra_point
    uvhdu.header['CRPIX6'] = template_uvfits[0].header['CRPIX6']
    uvhdu.header['CDELT6'] = template_uvfits[0].header['CDELT6']

    uvhdu.header['CTYPE7'] = template_uvfits[0].header['CTYPE7']
    uvhdu.header['CRVAL7'] = dec_point
    uvhdu.header['CRPIX7'] = template_uvfits[0].header['CRPIX7']
    uvhdu.header['CDELT7'] = template_uvfits[0].header['CDELT7']

    ## Write the parameters scaling explictly because they are omitted if default 1/0
    uvhdu.header['PSCAL1'] = 1.0
    uvhdu.header['PZERO1'] = 0.0
    uvhdu.header['PSCAL2'] = 1.0
    uvhdu.header['PZERO2'] = 0.0
    uvhdu.header['PSCAL3'] = 1.0
    uvhdu.header['PZERO3'] = 0.0
    uvhdu.header['PSCAL4'] = 1.0
    uvhdu.header['PZERO4'] = 0.0
    uvhdu.header['PSCAL5'] = 1.0

    uvhdu.header['PZERO5'] = float(int_jd)

    ##Old observation parameters that were/are needed in CHIPS
    uvhdu.header['OBJECT']  = 'Undefined'
    uvhdu.header['OBSRA']   = ra_point
    uvhdu.header['OBSDEC']  = dec_point

    ##ANTENNA TABLE MODS======================================================================

    template_uvfits[1].header['FREQ'] = freq_cent
    template_uvfits[1].header['ARRNAM'] = 'MWA'

    ##MAJICK uses this date to set the LST
    dmy, hms = date.split()
    day,month,year = map(int,dmy.split('-'))
    hour,mins,secs = map(float,hms.split(':'))

    rdate = "%d-%02d-%02dT%02d:%02d:%.2f" %(year,month,day,hour,mins,secs)

    template_uvfits[1].header['RDATE'] = rdate

    ## Create hdulist and write out file
    hdulist = fits.HDUList(hdus=[uvhdu,template_uvfits[1]])
    hdulist.writeto(output_uvfits_name,overwrite=True)
    hdulist.close()

def enh2xyz(east,north,height,latitude=MWA_LAT*D2R):
    '''Takes east, north, height coords out of the metafits file,
    and returns local X,Y,Z coords to put in the uvfits file'''
    sl = sin(latitude)
    cl = cos(latitude)
    X = -north*sl + height*cl
    Y = east
    Z = north*cl + height*sl
    return X,Y,Z

def load_data(filename=None,num_baselines=None,num_freq_channels=None,num_time_steps=None):
    '''Read the WODEN binary output and shove into a numpy arrays, ready to be put
    into a uvfits file'''

    with open(filename,'rb') as f:
        read_data = f.read()
    f.close()

    data = frombuffer(read_data,dtype=float32)

    n_data = num_time_steps * num_baselines
    v_container = zeros((n_data,1,1,1,num_freq_channels,4,3))
    uus = zeros(n_data)
    vvs = zeros(n_data)
    wws = zeros(n_data)

    num_visi = num_time_steps * num_freq_channels * num_baselines

    u_base = 0
    v_base = num_visi
    w_base = 2*num_visi
    re_base = 3*num_visi
    im_base = 4*num_visi

    num_cols = 5
    for time_ind in arange(num_time_steps):

        time_step = num_baselines * time_ind * num_freq_channels
        u_ind = u_base + time_step
        v_ind = v_base + time_step
        w_ind = w_base + time_step

        uus[time_ind*num_baselines:(time_ind + 1)*num_baselines] = data[u_ind:u_ind+num_baselines] / VELC
        vvs[time_ind*num_baselines:(time_ind + 1)*num_baselines] = data[v_ind:v_ind+num_baselines] / VELC
        wws[time_ind*num_baselines:(time_ind + 1)*num_baselines] = data[w_ind:w_ind+num_baselines] / VELC

        for freq_ind in arange(num_freq_channels):

            freq_step = num_baselines * (time_ind * num_freq_channels + freq_ind)

            real_ind = re_base + freq_step
            imag_ind = im_base + freq_step

            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,0,freq_ind,0,0] = data[real_ind:real_ind+num_baselines]
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,0,freq_ind,0,1] = data[imag_ind:imag_ind+num_baselines]
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,0,freq_ind,1,0] = data[real_ind:real_ind+num_baselines]
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,0,freq_ind,1,1] = data[imag_ind:imag_ind+num_baselines]

            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,0,freq_ind,0,2] = ones(num_baselines)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,0,freq_ind,1,2] = ones(num_baselines)

    return uus, vvs, wws, v_container

def write_json(ra0=None,dec0=None,num_freqs=None,num_time_steps=None,
               cat_filename=None,metafits_filename=None,band_nums=None,
               json_name=None,freq_res=None,time_res=None):
    '''Populate a json parameter file used to run WODEN'''

    outfile = open(json_name,'w+')

    outfile.write('{\n')
    outfile.write('  "ra0": %.10f,\n' %ra0)
    outfile.write('  "dec0": %.10f,\n' %dec0)
    outfile.write('  "num_freqs": %d,\n' %num_freqs)
    outfile.write('  "num_time_steps": %d,\n' %num_time_steps)
    outfile.write('  "cat_filename": "%s",\n' %cat_filename)
    outfile.write('  "metafits_filename": "%s",\n' %metafits_filename)
    outfile.write('  "time_res": %.5f,\n' %time_res)
    outfile.write('  "frequency_resolution": %.3f,\n' %freq_res)

    if len(band_nums) == 1:
        band_str = '[%d]' %band_nums[0]
    else:

        band_str = '[%d' %band_nums[0]
        for band in band_nums[1:-1]:
            band_str += ',%d' %band
        band_str += ',%d]' %band_nums[-1]

    outfile.write('  "band_nums": %s\n' %band_str)
    outfile.write('}\n')
    outfile.close()

if __name__ == "__main__":
    import argparse

    print("You are using WODEN commit %s" %gitlabel)

    parser = argparse.ArgumentParser(description='Run the woden simulator')
    parser.add_argument('--ra0', type=float,
        help='RA of the desired phase centre (deg)')
    parser.add_argument('--dec0', type=float,
        help='Dec of the desired phase centre (deg)')
    parser.add_argument('--cat_filename',
        help='RTS-v2-like srclist to simulate')
    parser.add_argument('--metafits_filename',
        help='Metafits file to base the simulation on')
    parser.add_argument('--template_uvfits', default="%s/template_MWA_128T.uvfits" %WODEN_DIR,
        help='Template uvfits to base outputs on - defaults to template_MWA_128T.uvfits')
    parser.add_argument('--output_uvfits_prepend',
        help='Prepend name for uvfits - will append band%%02d.uvfits %%band_num at the end')
    parser.add_argument('--num_freq_channels', type=int,
        help='Number of fine frequency channels to simulate')
    parser.add_argument('--num_time_steps', type=int,
        help='The number of time steps to simualte')
    parser.add_argument('--band_nums', default='all',
        help='Defaults to running all 24 course bands. Alternatively, enter required numbers delineated by commas, e.g. --band_nums=1,7,9')
    parser.add_argument('--no_tidy', default=False, action='store_true',
        help='Defaults to deleting output binary files from woden and json files. Add this flag to not delete those files')
    parser.add_argument('--nvprof', default=False, action='store_true',
        help='Add to switch on the nvidia profiler when running woden')
    parser.add_argument('--freq_res', type=float,default=False,
        help='Fine channel frequnecy resolution (Hz) - will default to what is in the metafits')
    parser.add_argument('--time_res', type=float,default=False,
        help='Time resolution (s) - will default to what is in the metafits')

    args = parser.parse_args()

    # def false_test(option,name):
    #     if option == False:
    #         print '-----------------------------------------------------'
    #         print '%s not entered but is required. Exiting now!!' %name
    #         print '-----------------------------------------------------'
    #         exit()
    #     else:
    #         return option
    #
    # ##Get inputs
    # time = false_test(options.time,'"time"')
    # metafits = false_test(options.metafits,'"metafits"')
    # output_dir = false_test(options.output_dir,'"output_dir"')
    # #srclist = false_test(options.srclist,'"srclist"')
    # time_int = false_test(options.time_int,'"time_int"')
    # oskar_uvfits_tag = false_test(options.oskar_uvfits_tag,'"oskar_uvfits_tag"')

    if args.band_nums == 'all':
        band_nums = range(1,25)
    else:
        try:
            band_nums = list(array(args.band_nums.split(','),dtype=int))
        except:
            print('-----------------------------------------------------')
            print('Failed to convert --band_nums into something sensible. Exiting now!!')
            print('-----------------------------------------------------')
            exit()

    ##Setup simulation parameters
    num_time_steps = args.num_time_steps
    num_freq_channels = args.num_freq_channels

    with fits.open(args.metafits_filename) as f:
        initial_date = f[0].header['DATE-OBS']
        ##Change in to oskar date format
        date,time = initial_date.split('T')
        year,month,day = date.split('-')
        oskar_date = "%s-%s-%s %s" %(day,month,year,time)

        ##Get the east, north, height antenna positions from the metafits
        east = f[1].data['East']
        north = f[1].data['North']
        height = f[1].data['Height']

        ##Read observation parameters from the metafits file
        time_res = float(f[0].header['INTTIME'])
        ch_width = float(f[0].header['FINECHAN'])*1e+3
        freqcent = float(f[0].header['FREQCENT'])*1e+6
        b_width = float(f[0].header['BANDWDTH'])*1e+6
        base_low_freq = freqcent - (b_width/2) - (ch_width/2)

        ##Replace default obs settings if specified in arguments
        if args.time_res:
            time_res = args.time_res

        if args.freq_res:
            ch_width = args.freq_res

        f.close()

    ##Write json file
    json_name = 'run_woden_%s.json' %args.band_nums
    write_json(ra0=args.ra0,dec0=args.dec0,num_freqs=num_freq_channels,num_time_steps=num_time_steps,
                   cat_filename=args.cat_filename,metafits_filename=args.metafits_filename,band_nums=band_nums,
                   json_name=json_name,freq_res=ch_width,time_res=time_res)

    ##Check the uvfits prepend to make sure we end in .uvfits
    output_uvfits_prepend = args.output_uvfits_prepend
    if output_uvfits_prepend[-7:] == '.uvfits': output_uvfits_prepend = output_uvfits_prepend[-7:]

    if args.nvprof:
        command('nvprof --log-file %s_%02d_nvprof.txt  %s/woden %s' %(output_uvfits_prepend,band_nums[0],WODEN_DIR,json_name))
    else:
        command('%s/woden %s' %(WODEN_DIR,json_name))

    ##Prepare the uvfits information
    ##Create and fill a layout array
    array_layout = zeros((int(len(east)/2),3))
    ##Tiles are listed as YY,XX,YY,XX,etc so only use half positions
    selection = arange(0,len(east),2)
    array_layout[:,0] = east[selection]
    array_layout[:,1] = north[selection]
    array_layout[:,2] = height[selection]

    X,Y,Z = enh2xyz(east[selection], north[selection],height[selection],MWA_LAT*D2R)

    ##Get the central frequency channels, used in the uvfits header
    central_freq_chan = int(floor(num_freq_channels / 2.0))

    ##Loop over coarse frequency band
    for band in band_nums:
        print('Converting binary to uvfits band',band)

        output_uvfits_name = output_uvfits_prepend + '_band%02d.uvfits' %band

        band_low_freq = base_low_freq + (band - 1)*1.28e+6
        central_freq_chan_value = band_low_freq + (1.28e+6 / 2.0)

        with fits.open(args.template_uvfits) as template_uvfits:
            template_data = template_uvfits[0].data
            template_baselines = template_uvfits[0].data['BASELINE'].copy()
            num_baselines = len(template_data)

            template_uvfits[1].data['STABXYZ'][:,0] = X
            template_uvfits[1].data['STABXYZ'][:,1] = Y
            template_uvfits[1].data['STABXYZ'][:,2] = Z

            int_jd, float_jd = calc_jdcal(oskar_date)

            ##Need an array the length of number of baselines worth of the fractional jd date
            float_jd_array = ones(num_baselines)*float_jd

            n_data = num_time_steps * num_baselines

            ##Create empty data structures for final uvfits file
            baselines_array = zeros(n_data)
            date_array = zeros(n_data)

            filename = "output_visi_band%02d.dat" %band
            uus,vvs,wws,v_container = load_data(filename=filename,num_baselines=num_baselines,
                                                num_freq_channels=num_freq_channels,num_time_steps=num_time_steps)

            for time_ind,time in enumerate(arange(0,num_time_steps*time_res,time_res)):
                time_ind_lower = time_ind*num_baselines
                baselines_array[time_ind_lower:time_ind_lower+num_baselines] = template_baselines

                ##Fill in the fractional julian date, after adding on the appropriate amount of
                ##time - /(24*60*60) because julian number is a fraction of a whole day
                adjust_float_jd_array = float_jd_array + (float(time) / (24.0*60.0*60.0))
                date_array[time_ind_lower:time_ind_lower+num_baselines] = adjust_float_jd_array

            create_uvfits(v_container=v_container,freq_cent=central_freq_chan_value,ra_point=args.ra0,dec_point=args.dec0,
                        output_uvfits_name=output_uvfits_name,uu=uus,vv=vvs,ww=wws,baselines_array=baselines_array,
                        date_array=date_array,date=oskar_date,central_freq_chan=central_freq_chan,ch_width=ch_width,
                        template_uvfits=template_uvfits,int_jd=int_jd)


        ##Tidy up or not
        if args.no_tidy:
            pass
        else:
            command("rm %s" %filename)

    ##Tidy up or not
    if args.no_tidy:
        pass
    else:
        command("rm %s" %json_name)
