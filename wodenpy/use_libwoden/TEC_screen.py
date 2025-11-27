import numpy as np
from ctypes import POINTER, c_double, c_float, c_int, Structure
from wodenpy.use_libwoden.woden_settings import Woden_Settings_Python

class TEC_Screen_Python(object):
    """A class structured equivalently to a `TEC_screen_t` struct, used by
    the C/C++/GPU code in libwoden_float.so or libwoden_double.so. Retain a
    copy of the settings in Python for easy access and manipulation."""
    def __init__(self):
        self.screen = None
        self.resolution = None
        self.screen_size = None
        self.height = None

def create_TEC_screen_struct(precision : str = "double"):
    """Creates a `TEC_Screen` class structured equivalently to a `TEC_screen_t`
    struct in the C/CUDA code. Created dynamically based on the `precision`,
    to match the compile time precision flag `-DUSE_DOUBLE` in the C code.

    Parameters
    ----------
    precision : str, optional
        Either "float" or "double:, by default "double"

    Returns
    -------
    TEC_Screen
        The TEC_Screen class structured equivalently to a `TEC_screen_t` struct
    """
    
    if precision == "float":
        c_user_precision = c_float
    else:
        c_user_precision = c_double

    class TEC_Screen(Structure):
        """A class structured equivalently to a `TEC_screen_t` struct, used by 
        the C and CUDA code in libwoden_float.so or libwoden_double.so.
        
        Created by the function `create_TEC_screen_struct`, which sets
        `user_precision_t` to either `c_float` or `c_double`.
        
        :cvar POINTER(c_user_precision) screen: pixel data of TEC screen
        :cvar c_int resolution: resolution of TEC screen
        :cvar c_user_precision height: height of TEC screen (meters)
        """
        
        _fields_ = [("screen", POINTER(c_user_precision)),
                ("resolution", c_int),
                ("screen_size", c_user_precision),
                ("height", c_user_precision),]
        
    return TEC_Screen

def setup_TEC_screen_python(woden_settings_python : Woden_Settings_Python) -> TEC_Screen_Python:
    """Given the populated Woden_Settings_Python class, set up the `Array_Layout_Python` class, 
    and fill it with the correct values.

    Parameters
    ----------
    woden_settings_python : Woden_Settings_Python
        Populated Woden_Settings_Python class.

    Returns
    -------
    TEC_screen : TEC_Screen_Python
        Initialised TEC_Screen_Python class.
    """

    TEC_screen = TEC_Screen_Python()

    if (woden_settings_python.do_ionosphere == 0):
        TEC_screen.resolution = 0
        TEC_screen.height = 0
        TEC_screen.screen_size = 0
        TEC_screen.screen = np.empty(0)
        return TEC_screen

    TEC_grad_x = woden_settings_python.TEC_grad_x
    TEC_grad_y = woden_settings_python.TEC_grad_y

    TEC_screen.resolution = 10000
    TEC_screen.height = 200000
    TEC_screen.screen_size = 200000

    xaxis = np.linspace(-TEC_screen.screen_size / 2, TEC_screen.screen_size / 2, TEC_screen.resolution)
    yaxis = np.linspace(-TEC_screen.screen_size / 2, TEC_screen.screen_size / 2, TEC_screen.resolution)
    x, y = np.meshgrid(xaxis, yaxis)
    TEC_screen_2d = TEC_grad_x * x + TEC_grad_y * y

    TEC_screen.screen = TEC_screen_2d.flatten()
    
    return TEC_screen

TEC_Screen = create_TEC_screen_struct()

def convert_TEC_screen_to_ctypes(TEC_screen_python : TEC_Screen_Python,
                                   TEC_screen_ctypes : TEC_Screen, #type: ignore
                                   precision='double') -> TEC_Screen: #type: ignore
    
    if precision == "float":
        c_user_precision = c_float
    else:
        c_user_precision = c_double
    
    TEC_screen_ctypes.screen = TEC_screen_python.screen.ctypes.data_as(POINTER(c_user_precision))
    TEC_screen_ctypes.resolution = TEC_screen_python.resolution
    TEC_screen_ctypes.screen_size = TEC_screen_python.screen_size
    TEC_screen_ctypes.height = TEC_screen_python.height
    
    return TEC_screen_ctypes