"""Functions for setting up the logger for wodenpy."""

import logging
import logging.handlers
from ctypes import CDLL, CFUNCTYPE, c_char_p
import ctypes
import importlib_resources
import wodenpy
from multiprocessing import Process
from logging.handlers import QueueHandler, QueueListener
import sys
from wodenpy.use_libwoden.beam_settings import BeamTypes, BeamGroups
import numpy as np
import argparse

class MultiLineFormatter(logging.Formatter):
    """Multi-line formatter for logging messages. Means that the indentation
    comes out inline with end of time and level message, if you happen submit
    a long multi-line string to a single log call.
    
    Taken from https://stackoverflow.com/questions/58590731/how-to-indent-multiline-message-printed-by-python-logger"""
    def get_header_length(self, record):
        """Get the header length of a given record."""
        return len(super().format(logging.LogRecord(
            name=record.name,
            level=record.levelno,
            pathname=record.pathname,
            lineno=record.lineno,
            msg='', args=(), exc_info=None
        )))

    def format(self, record):
        """Format a record with added indentation."""
        indent = ' ' * self.get_header_length(record)
        head, *trailing = super().format(record).splitlines(True)
        return head + ''.join(indent + line for line in trailing)

    
def set_logger_header(logger : logging.Logger, gitlabel : str):
    """
    Writes a snazzy header to the logger, including the version or git hash label.

    Parameters
    ----------
    logger : logging.Logger
        The logger instance to which the header will be added.
    gitlabel : str
        The version or git hash label to be included in the header.

    Returns
    -------
    """

    
    logo_string = rf"""
                 )  (              )  
     (  (     ( /(  )\ )        ( /(  
     )\))(   ')\())(()/(   (    )\()) 
    ((_)()\ )((_)\  /(_))  )\  ((_)\  
    _(())\_)() ((_)(_))_  ((_)  _((_) 
    \ \((_)/ // _ \ |   \ | __|| \| | 
     \ \/\/ /| (_) || |) || _| | .` | 
      \_/\_/  \___/ |___/ |___||_|\_| 
      
    You are using wodenpy version/git hash: {gitlabel}
    """
    
    logger.info(logo_string)
    
def summarise_input_args(logger : logging.Logger, args : argparse.Namespace):
    """
    Summarises and logs the input arguments after parsing.
    
    Parameters
    ----------
    logger : logging.Logger
        The logger instance used to log the summarized input arguments.
    args : argparse.Namespace
        The parsed input arguments containing various configuration parameters.
    Attributes of args
    ------------------
    ra0 : float
        Right ascension of the phase center in degrees.
    dec0 : float
        Declination of the phase center in degrees.
    latitude : float
        Central latitude of the array in degrees.
    longitude : float
        Central longitude of the array in degrees.
    array_height : float
        Height of the array center in meters.
    lowest_channel_freq : float
        Lowest channel frequency in Hz.
    freq_res : float
        Frequency resolution in Hz.
    coarse_band_width : float
        Coarse band bandwidth in Hz.
    num_freq_channels : int
        Number of frequency channels per coarse band.
    date : str
        Start date of the observation.
    time_res : float
        Time resolution in seconds.
    num_time_steps : int
        Number of time steps.
    array_layout : str
        Source of the array layout information.
    num_antennas : int
        Number of antennas.
    metafits_filename : str
        Filename of the metafits file (if applicable).
    beam_ms_path : str
        Path to the measurement set (if applicable).
    output_dir : str
        Directory where outputs will be written.
    """
    
    
    input_str = "Input arguments after parsing:"
    
    input_str += f"\n\tPhase centre: {args.ra0:.5f}, {args.dec0:.5f} deg"
    input_str += f"\n\tArray central latitude: {args.latitude:.3f} deg"
    input_str += f"\n\tArray central longitude: {args.longitude:.3f} deg"
    input_str += f"\n\tArray central height: {args.array_height:.3f} m"
    input_str += f"\n\tLowest channel frequency: {args.lowest_channel_freq/1e+6:.3f} MHz"
    input_str += f"\n\tChannel frequency resolution: {args.freq_res/1e+3:.3f} kHz"
    input_str += f"\n\tCoarse band bandwidth: {args.coarse_band_width/1e+6:.3f} MHz"
    input_str += f"\n\tNum channels per coarse band: {args.num_freq_channels}"
    input_str += f"\n\tStart date: {args.date}"
    input_str += f"\n\tTime resolution: {args.time_res} (s)"
    input_str += f"\n\tNum time steps: {args.num_time_steps}"
    
    if args.array_layout == "from_the_metafits":
        input_str += f"\n\tHave read {args.num_antennas} antenna positions from metafits file: {args.metafits_filename}"
    elif args.array_layout == "from_ms":
        input_str += f"\n\tHave read {args.num_antennas} antenna positions from measurement set: {args.beam_ms_path}"
    else:
        input_str += f"\n\tHave read {args.num_antennas} antenna positions from array layout file: {args.array_layout}"
        
    input_str += f"\n\tWill write outputs to: {args.output_dir}"
    
    logger.info(input_str)
    
def get_log_callback(logger: logging.Logger, logging_level: int = logging.DEBUG):
    """
    Creates a C callback function that logs messages using the provided Python logger.
    Parameters
    ----------
    logger : logging.Logger
        The Python logger instance to use for logging messages.
    logging_level : int, optional
        The logging level to use (default is logging.DEBUG).
    Returns
    -------
    CFUNCTYPE
        A C callback function that logs messages using the provided Python logger.
    """
    
    LogCallbackType = CFUNCTYPE(None, c_char_p)

    # Define the Python logging function for libwoden C library
    def log_callback(message):
        if logging_level == logging.DEBUG:
            logger.debug(f"libwoden: {message.decode('utf-8')}")
        else:
            logger.info(f"libwoden: {message.decode('utf-8')}")

    # Wrap the Python function as a C callback
    c_log_callback = LogCallbackType(log_callback)
    
    return c_log_callback

def set_woden_logger(logging_level: int = logging.DEBUG, log_file: str = False) -> logging.Logger:
    """
    Sets up a logger for WODEN with a specified logging level and optional log file.
    
    Parameters
    ----------
    logging_level : int, optional
        The logging level to set for the logger (default is logging.DEBUG).
        Accepted values are logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL.
    log_file : str, optional
        The file path to log to. If False, logging to a file is disabled (default is False).
        
    Returns
    -------
    logger : logging.Logger
        Configured logger instance.
    """
    
    
    formatter = MultiLineFormatter("%(asctime)s - %(levelname)s - %(message)s",
                                   '%Y-%m-%d %H:%M:%S')
    
    logger = logging.getLogger(__name__)
    logger.setLevel(logging_level)
    
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging_level)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    
    if log_file:
        file_handler = logging.FileHandler(log_file, 'w+')
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging_level)
        logger.addHandler(file_handler)
    
    return logger

def simple_logger(logging_level: int = logging.DEBUG):
    """
    Use this to set a default logger for wodenpy functions when the user doesn't provide one.
    Basically useful for unit tests.
    
    Parameters
    ----------
    logging_level : int, optional
        The logging level to use for the logger. Default is logging.DEBUG.
        Accepted values are logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL.
        
    Returns
    -------
    logger : logging.Logger
        Configured logger instance.
    """
    
    logger = logging.getLogger("WODEN")
    
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s',
                                           '%Y-%m-%d %H:%M:%S'))
    logger.addHandler(handler)
    
    return logger


def _mwa_beam_settings_string(logger: logging.Logger, woden_settings: object,
                              args: argparse.Namespace, package = 'hyperbeam') -> str:
    """
    Generate a string describing the MWA primary beam settings.
    
    Parameters
    ----------
    logger : logging.Logger
        The logger instance to use for logging messages.
    woden_settings : object
        An object containing WODEN settings, including the path to the hdf5 beam file and FEE ideal delays.
    args : argparse.Namespace
        Command-line arguments containing flags and filenames for MWA dipole amplitudes and flags.
        
    Returns
    -------
    str
        A formatted string describing the MWA primary beam settings.
    """
    
    out_string = f"Using MWA primary beam via {package} with the following parameters:\n" \
                    f"\thdf5 file: { woden_settings.hdf5_beam_path}\n" \
                    f"\tdelays: {woden_settings.FEE_ideal_delays[:16]}"
    if args.use_MWA_dipamps:
        out_string += f"\n\twill use dipole amplitudes from given metafits file"
    else:
        out_string += "\n\tsetting all dipole amplitudes to 1.0"
    
    if args.use_MWA_dipflags:
        out_string += f"\n\twill use dipole flags from given metafits file"
    else:
        out_string += "\n\twill not flag any dipoles"
        
    if args.use_MWA_dipamps or args.use_MWA_dipflags:
        out_string += f"\n\tmetafits file: {args.metafits_filename}"
            
    return out_string

def _everybeam_settings_string(logger: logging.Logger, woden_settings: object, args: argparse.Namespace) -> str:
    """
    Generates a string describing the settings for the EveryBeam primary beam.
    
    Parameters
    ----------
    logger : logging.Logger
        The logger instance to use for logging messages.
    woden_settings : object
        An object containing the WODEN settings, including the beam type and paths.
    args : argparse.Namespace
        The command-line arguments, including paths and pointing information.
        
    Returns
    -------
    str
        A formatted string describing the EveryBeam primary beam settings.
    """
    
    if woden_settings.beamtype == BeamTypes.EB_OSKAR.value:
        beam = "OSKAR"
    elif woden_settings.beamtype == BeamTypes.EB_LOFAR.value:
        beam = "LOFAR"
    # elif woden_settings.beamtype == BeamTypes.EB_MWA.value:
    #     beam = "MWA"
    
    out_string = f"Will run with EveryBeam {beam} primary beam, based on this measurement set:"
    out_string +=  f"\n\t{args.beam_ms_path}"
    
    # if beam == "MWA":
    #     out_string += f"\nUsing the following hdf5 file:"
    #     out_string += f"\n\t{woden_settings.hdf5_beam_path}"
        
    if args.pointed_ms_file_name:
        out_string += f"\nCreated the following minimal MS to point the beam:"
        out_string += f"\n\t{args.pointed_ms_file_name}"
        
    # if beam == "LOFAR" or beam == "OSKAR":
    out_string += f"\nPrimary beam is pointed at RA,Dec = {args.eb_ra_point}, {args.eb_dec_point} deg"
        
    return out_string


def _uvbeam_settings_string(logger: logging.Logger, args: argparse.Namespace,
                              beam : str) -> str:
    """
    Generates a string describing the settings for the EveryBeam primary beam.
    
    Parameters
    ----------
    logger : logging.Logger
        The logger instance to use for logging messages.
    woden_settings : object
        An object containing the WODEN settings, including the beam type and paths.
    beam : str
        Name of the beam model (e.g., 'HERA').
        
    Returns
    -------
    str
        A formatted string describing the EveryBeam primary beam settings.
    """
    
    if hasattr(args, 'cst_freqs'):
    
        out_string = f"Will create a {beam} beam from CST files listed in this file:"
        out_string +=  f"\n\t{args.cst_file_list}"
        out_string += f"\nHave read the following frequencies from this file:"
        out_string += f"\n\t{args.cst_freqs}"

    else:
        out_string = f"Will create a {beam} beam from the following file:"
        out_string +=  f"\n\t{args.uvbeam_file_path}"
        
    return out_string

def log_chosen_beamtype(logger: logging.Logger, woden_settings: object,
                        args : argparse.Namespace):
    """
    Logs information about the chosen beam type based on the provided settings.
    
    Parameters
    ----------
    logger : logging.Logger
        The logger instance used to log the information.
    woden_settings : wodenpy.use_libwoden.woden_settings.Woden_Settings_Python
        An object containing the settings for WODEN, including the beam type and related parameters.
    args : object
        The args as parsed by argparse, containing command-line arguments and options.
    
    """
    
    if woden_settings.beamtype == BeamTypes.NO_BEAM.value:
        logger.info("No primary beam selected, no beam attenuation will be applied.")
        
    elif woden_settings.beamtype == BeamTypes.GAUSS_BEAM.value:
        logger.info("Using Gaussian primary beam with the following parameters:\n"
                    f"\tLocked to pointing: HA {np.degrees(woden_settings.lst_base - woden_settings.gauss_ra_point):.1f} deg, "
                    f"Dec {np.degrees(woden_settings.gauss_dec_point)} deg\n"
                    f"\tFWHM: {woden_settings.gauss_beam_FWHM} at "
                    f"reference frequency: {woden_settings.gauss_beam_ref_freq/1e+6} MHz")

    elif woden_settings.beamtype == BeamTypes.FEE_BEAM_INTERP.value or woden_settings.beamtype == BeamTypes.FEE_BEAM.value:
        mwa_beam_string = _mwa_beam_settings_string(logger, woden_settings, args)
        logger.info(mwa_beam_string)
        
    elif woden_settings.beamtype == BeamTypes.ANALY_DIPOLE.value:
        logger.info("Using an analytical dipole primary beam (a.k.a each element is an MWA dipole e.g. EDA2 array).")

    elif woden_settings.beamtype == BeamTypes.MWA_ANALY.value:
        logger.info("Using MWA analytic primary beam with:\n" \
                    f"\tdelays: {woden_settings.FEE_ideal_delays}")

    elif woden_settings.beamtype in BeamGroups.eb_ms_beam_values:
        eb_beam_string = _everybeam_settings_string(logger, woden_settings, args)
        logger.info(eb_beam_string)
        
    elif woden_settings.beamtype == BeamTypes.EB_MWA.value:
        mwa_beam_string = _mwa_beam_settings_string(logger, woden_settings, args, 'EveryBeam MWA')
        logger.info(mwa_beam_string)
        
    elif woden_settings.beamtype == BeamTypes.UVB_MWA.value:
        mwa_beam_string = _mwa_beam_settings_string(logger, woden_settings, args, 'pyuvdata.UVBeam')
        logger.info(mwa_beam_string)
        
    elif woden_settings.beamtype == BeamTypes.UVB_HERA.value:
        hera_beam_string = _uvbeam_settings_string(logger, args, 'HERA')
        logger.info(hera_beam_string)

    else:
        logger.error("Primary beam type not recognised. This shouldn't be possible "
                     "if you've used wodenpy.woden_setup.run_setup.check_args().")

