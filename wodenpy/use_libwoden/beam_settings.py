"""Classes to match the BeamTypes enum in libwoden, and to group beam types for WODEN.
Beam group are useful to calculate associated values for specific beam types,
e.g. azimuth/altitude, hour angle/declination, etc.
"""
from enum import Enum, auto

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
        UVB_MWA (int): pyuvdata UVData MWA beam type.
        UVB_MWA (int): pyuvdata UVData HERA beam type.
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
    UVB_MWA = 9
    UVB_HERA = 10
    
class BeamGroups:
    """A class to represent different groups of beam types.
    Each group contains a list of `BeamTypes` values that are used in WODEN.
    
    Attributes
    ----------
    eb_beam_values : list
        All the EveryBeam beam types.
    azza_beam_values : list
        All beams that need azimuth/altitude calculated before sending to CPU/GPU.
    hadec_beam_values : list
        All beams that need Hour Angle/Declination calculated/filled before sending to CPU/GPU.
    off_cardinal_beam_values : list
        All beams that have off-cardinal dipole alignments, e.g. not N/S or E/W at 0,90 degrees,
        but rather 45, 135 degrees. Currently empty pending further development (something in
        LOFAR is definitely off-cardinal but needs more research).
    python_calc_beams : list
        All beams that are calculated in Python before sending to CPU/GPU (e.g.
        any `UVBeam` beam types).
    uvbeam_beams : list
        All the UVBeam beam types.
    """
    
    eb_beam_values = [BeamTypes.EB_OSKAR.value, BeamTypes.EB_LOFAR.value, BeamTypes.EB_MWA.value]
    ##these are beam types that need azimuth/altitude calculated before sending to CPU/GPU via ctypes
    azza_beam_values = [BeamTypes.MWA_ANALY.value, BeamTypes.FEE_BEAM.value,
                        BeamTypes.FEE_BEAM_INTERP.value, BeamTypes.ANALY_DIPOLE.value,
                        BeamTypes.EB_MWA.value]
    
    needs_MWA_delays = [BeamTypes.MWA_ANALY.value, BeamTypes.FEE_BEAM.value,
                        BeamTypes.FEE_BEAM_INTERP.value, BeamTypes.EB_MWA.value]
    
    needs_MWA_hdf5_path = [BeamTypes.FEE_BEAM.value, BeamTypes.FEE_BEAM_INTERP.value,
                           BeamTypes.EB_MWA.value]
    
    hadec_beam_values = [BeamTypes.GAUSS_BEAM.value, BeamTypes.MWA_ANALY.value]
    # needs_fee_hdf5 = [BeamTypes.FEE_BEAM.value, BeamTypes.FEE_BEAM_INTERP, BeamTypes.EB_MWA.value]
    # needs_mwa_delays = [BeamTypes.FEE_BEAM.value, BeamTypes.FEE_BEAM_INTERP, BeamTypes.MWA_ANALY.value]
    ##TODO investigate what beam types are off-cardinal, and add to this list below
    off_cardinal_beam_values = []
    python_calc_beams = [BeamTypes.UVB_MWA.value, BeamTypes.UVB_HERA.value]
    uvbeam_beams = [BeamTypes.UVB_MWA.value, BeamTypes.UVB_HERA.value]
        
        