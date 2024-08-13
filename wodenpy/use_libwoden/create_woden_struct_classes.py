from wodenpy.use_libwoden.skymodel_structs import create_components_struct, create_source_struct, create_source_catalogue_struct, setup_source_catalogue



class Woden_Struct_Classes:
    """Create and store the classes for the WODEN data structures.
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