from sys import path
import os
import unittest
import numpy as np
import ctypes

# ##Code we are testing
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
from wodenpy.wodenpy_setup import run_setup
from wodenpy.use_libwoden.skymodel_structs import setup_chunked_source, setup_source_catalogue
import numpy.testing as npt

##annoying path hack to find where the C library is
test_dir = os.environ['CMAKE_CURRENT_SOURCE_DIR'] + "/../../../build/cmake_testing/wodenpy/use_libwoden/"

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

D2R = np.pi/180.0

##Vehicle for running tests
class Test(unittest.TestCase):
    """"""
        
    def read_in_C_functions(self, precision="double"):
        """Read in the C library function that reads the ctypes source_catalogue
        structure and writes the content to a text file"""
        
        ## Read in the C library for float version
        libwoden = ctypes.cdll.LoadLibrary(f"{test_dir}/libread_source_catalogue_{precision}.so")
        self.read_source_catalogue = libwoden.read_source_catalogue
        self.read_source_catalogue.argtypes = [ctypes.POINTER(self.woden_struct_classes.Source_Catalogue)]
        
        # ## Read in the C library for double version
        # libwoden_double = ctypes.cdll.LoadLibrary(f"{test_dir}/libread_source_catalogue_double.so")
        # self.read_source_catalogue_double = libwoden_double.read_source_catalogue
        # self.read_source_catalogue_double.argtypes = [ctypes.POINTER(ws.source_catalogue_Double)]

    
    def read_C_output_textfile(self):
        
        output_dict = {}
        
        with open('cropped_sky_models.txt', 'r') as infile:
            for line in infile.read().split('\n'):
                if line == '':
                    pass
                else:
                    output_dict[line.split()[0]] = float(line.split()[1])
            
        return output_dict

    def test_access_ra(self):
        """Super basic test. Just create a source catalogue with 2 sources, and
        fill in the RA values for the gauss components. Read in the source
        catalogue using C, and check we get back those RA. """
        
        
        num_shapelets = 0
        num_sources = 2
        num_comps = 3
        
        for precision in ["float", "double"]:
            
            print(precision)
        
            self.woden_struct_classes = Woden_Struct_Classes(precision)
            
            self.read_in_C_functions(precision)
            
            cropped_sky_models = setup_source_catalogue(self.woden_struct_classes.Source_Ctypes,
                                                        self.woden_struct_classes.Source_Catalogue, num_sources, num_shapelets,
                                                        precision)
            
            for i in range(num_sources):
                
                cropped_sky_models.sources[i] = self.woden_struct_classes.Source_Ctypes()
                
                components = cropped_sky_models.sources[i].gauss_components
                components.ras = (num_comps*ctypes.c_double)()
                
                for j in range(num_comps):
                    components.ras[j] = i + j
                    
        self.read_source_catalogue(cropped_sky_models)
        data = self.read_C_output_textfile()
        
        for i in range(num_sources):
            for j in range(num_comps):
                self.assertEqual(data[f'cropped_sky_models->sources[{i}].gauss_components.ras[{j}]'], i+j)
        

##Run the test
if __name__ == '__main__':
   unittest.main()