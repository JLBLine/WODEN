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

    
def set_logger_header(logger, gitlabel):
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
    
def summarise_input_args(logger, args):
    # logger.info("Input arguments:")
    # for arg in vars(args):
    #     logger.info(f"\t{arg}: {getattr(args, arg)}")
    
    input_str = "Input arguments after parsing:"
    
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
    
    logger.info(input_str)
    
def get_log_callback(logger, logging_level = logging.DEBUG):
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

def set_woden_logger(logging_level = logging.DEBUG, log_file = False):
    """WODEN main logger"""
    
    formatter = MultiLineFormatter("%(asctime)s - %(levelname)s - %(message)s",
                                   '%Y-%m-%d %H:%M:%S')
    
    logger = logging.getLogger(__name__)
    logger.setLevel(logging_level)
    
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging_level)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    
    if log_file:
        file_handler = logging.FileHandler(log_file, 'w+')
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging_level)
        logger.addHandler(file_handler)
    
    return logger

def simple_logger(logging_level = logging.DEBUG):
    """Use this to set a default logger for wodenpy functions when
    the user doesn't provide one. Basically useful for unit tests."""
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging_level,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    return logger


def _mwa_beam_settings_string(logger, woden_settings, args):
    
    out_string = "Using MWA primary beam via hyperbeam with the following parameters:\n" \
                    f"\thdf5 file: { woden_settings.hdf5_beam_path}\n" \
                    f"\tdelays: {woden_settings.FEE_ideal_delays}"
    if args.use_MWA_dipamps:
        out_string += f"\n\twill use dipole amplitudes from given metafits file"
    else:
        out_string += "\n\tsetting all dipole amplitudes to 1.0"
    
    if args.use_MWA_dipflags:
        out_string += f"\n\twill use dipole flags from given metafits file"
    else:
        out_string += "\n\twill not flag any dipoles"
        
    if args.use_MWA_dipamps or args.use_MWA_dipflags:
        out_string += f"\n\tmetafits file: {args.metafits_path}"
            
    return out_string

def log_chosen_beamtype(logger, woden_settings, args):
    
    if woden_settings.beamtype == BeamTypes.NO_BEAM.value:
        logger.info("No primary beam selected, no beam attenuation will be applied.")
        
    elif woden_settings.beamtype == BeamTypes.GAUSS_BEAM.value:
        logger.info("Using Gaussian primary beam with the following parameters:\n"
                    f"\tLocked to pointing: HA {np.degrees(woden_settings.lst_base - woden_settings.gauss_ra_point):.1f} deg, "
                    f"Dec {np.degrees(woden_settings.gauss_dec_point)} deg\n"
                    f"\tFWHM: {woden_settings.gauss_beam_FWHM} at "
                    f"reference frequency: {woden_settings.gauss_beam_ref_freq/1e+6} MHz")
    elif woden_settings.beamtype == BeamTypes.FEE_BEAM.value:
       mwa_beam_string = _mwa_beam_settings_string(logger, woden_settings, args)
       logger.info(mwa_beam_string)
       
    elif woden_settings.beamtype == BeamTypes.FEE_BEAM_INTERP.value:
        mwa_beam_string = _mwa_beam_settings_string(logger, woden_settings, args)
        logger.info(mwa_beam_string)
        
    elif woden_settings.beamtype == BeamTypes.ANALY_DIPOLE.value:
        logger.info("Using an analytical dipole primary beam (a.k.a each element is an MWA dipole e.g. EDA2 array).")

    elif woden_settings.beamtype == BeamTypes.MWA_ANALY.value:
        logger.info("Using MWA analytic primary beam with:\n" \
                    f"\tdelays: {woden_settings.FEE_ideal_delays}")

    elif woden_settings.beamtype == BeamTypes.EB_OSKAR.value:
        logger.info("Will run with EveryBeam OSKAR primary beam, based on this measurement set:\n"
                    f"\t{args.beam_ms_path}")

    elif woden_settings.beamtype == BeamTypes.EB_LOFAR.value:
        logger.info("Will run with EveryBeam LOFAR primary beam, based on this measurement set:\n"
                    f"\t{args.beam_ms_path}")

    elif woden_settings.beamtype == BeamTypes.EB_MWA.value:
        logger.info("Will run with EveryBeam MWA primary beam, based on this measurement set:\n"
                    f"\t{args.beam_ms_path}")
    
    else:
        logger.error("Primary beam type not recognised. This shouldn't be possible "
                     "if you've used wodenpy.woden_setup.run_setup.check_args().")

