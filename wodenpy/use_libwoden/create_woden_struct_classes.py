from wodenpy.use_libwoden.skymodel_structs import create_components_struct, create_source_struct, create_source_catalogue_struct, setup_source_catalogue
from wodenpy.use_libwoden.woden_settings import create_woden_settings_struct
from wodenpy.use_libwoden.visibility_set import create_visi_set_struct


class Woden_Struct_Classes:
    """Create and store the classes for the WODEN data structures.
    
    :cvar str precision: Either "float" or "double" to determine the precision of the data types.
    :cvar Components Components: Equivalent to the C struct `components_t`.
    :cvar Source Source: Equivalent to the C struct `source_t`.
    :cvar Source_Catalogue Source_Catalogue: Equivalent to the C struct `source_catalogue_t`.
    :cvar Woden_Settings Woden_Settings: Equivalent to the C struct `woden_settings_t`.
    :cvar Visi_Set Visi_Set: Equivalent to the C struct `visibility_set_t`.
    
    """
    def __init__(self, precision="double"):
        """Creates all necessary classes for the WODEN data structures. Uses
        the supplied `precision` to determine the precision of the data types;
        these are used to match the `-DUSE_DOUBLE` flag in the C/GPU code.
        
        

        Parameters
        ----------
        precision : str, optional
            "float" to use float precision, "double" to use double. Default "double".
        """
        self.precision = precision
        
        ##Sky model bits
        self.Components = create_components_struct(precision)
        self.Source = create_source_struct(self.Components)
        self.Source_Catalogue = create_source_catalogue_struct(self.Source)
        
        ##Woden settings
        self.Woden_Settings = create_woden_settings_struct(precision)
        
        ##visibility set
        self.Visi_Set = create_visi_set_struct(precision)