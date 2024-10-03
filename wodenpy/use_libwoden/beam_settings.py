import ctypes 
# import importlib_resources
import numpy as np
from typing import Union
from ctypes import POINTER, c_float, c_double, c_int, c_char
import os
import sys
from enum import Enum, auto

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes, Component_Info
from wodenpy.skymodel.chunk_sky_model import Skymodel_Chunk_Map


class BeamTypes(Enum):
    """
    Enumeration of beam types used in WODEN.

    This class allows us to label the WODEN primary beam models with a unique
    name/value, but as it's an enum each label only takes 8 bytes of memory, so we can
    stack loads of them into an array. We can also do numpy operations on them like np.where.

    Attributes:
        NO_BEAM (int): No beam type.
        GAUSS_BEAM (int): Gaussian beam type.
        FEE_BEAM (int): FEE beam type.
        ANALY_DIPOLE (int): Analytical dipole beam type.
        FEE_BEAM_INTERP (int): Interpolated FEE beam type.
        MWA_ANALY (int): MWA analytical beam type.
        EB_OSKAR (int): EveryBeam OSKAR beam type.
        EB_LOFAR (int): EveryBeam LOFAR beam type.
        EB_MWA (int): EveryBeam MWA beam type.
    """
    NO_BEAM = 0
    GAUSS_BEAM = 1
    FEE_BEAM = 2
    ANALY_DIPOLE = 3
    FEE_BEAM_INTERP = 4
    MWA_ANALY = 5
    EB_OSKAR = 6
    EB_LOFAR = 7
    EB_MWA = 8