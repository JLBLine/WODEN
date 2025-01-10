from ctypes import POINTER, c_double, c_int, Structure

class Array_Layout_Ctypes(Structure):
    """
    A ctypes structure representing the layout of an array of antennas.

    Attributes:
        ant_X (POINTER(c_double)): Pointer to an array of antenna X positions.
        ant_Y (POINTER(c_double)): Pointer to an array of antenna Y positions.
        ant_Z (POINTER(c_double)): Pointer to an array of antenna Z positions.
        X_diff_metres (POINTER(c_double)): Pointer to an array of X position differences in metres.
        Y_diff_metres (POINTER(c_double)): Pointer to an array of Y position differences in metres.
        Z_diff_metres (POINTER(c_double)): Pointer to an array of Z position differences in metres.
        ant_east (POINTER(c_double)): Pointer to an array of antenna east positions.
        ant_north (POINTER(c_double)): Pointer to an array of antenna north positions.
        ant_height (POINTER(c_double)): Pointer to an array of antenna height positions.
        latitude (c_double): The latitude of the array.
        num_baselines (c_int): The number of baselines in the array.
        num_tiles (c_int): The number of tiles in the array.
        lst_base (c_double): The local sidereal time of the array.
    """
    _fields_ = [("ant_X", POINTER(c_double)),
                ("ant_Y", POINTER(c_double)),
                ("ant_Z", POINTER(c_double)),
                ("X_diff_metres", POINTER(c_double)),
                ("Y_diff_metres", POINTER(c_double)),
                ("Z_diff_metres", POINTER(c_double)),
                ("ant_east", POINTER(c_double)),
                ("ant_north", POINTER(c_double)),
                ("ant_height", POINTER(c_double)),
                ("latitude", c_double),
                ("num_baselines", c_int),
                ("num_tiles", c_int),
                ("lst_base", c_double)]