"""
Test `wodenpy.skymodel.crop_below_horizon`, which should setup flags of one
above the horizon, zero below the horizon. Also can crop by entire SOURCE
or just COMPONENT (crop_type). Test using a dummy model of two SOURCEs, each with three
COMPONENTs, with four different LSTs. Various combinations of crop_type and
LSTs should yield specific horizon flags
"""

from sys import path
import os
import unittest
import numpy as np

# ##Code we are testing
from wodenpy.skymodel import read_yaml_skymodel
# import wodenpy
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, crop_below_horizon


D2R = np.pi/180.0

LST0 = 0.0
LST1 = np.pi / 2
LST2 = np.pi
LST3 = 230.0*D2R
LATITUDE = -0.46606083776035967

TEST_RAS = np.array([0.0, D2R, 2*D2R, 106*D2R, 108*D2R, 110*D2R])
TEST_DECS = np.array([-0.46606083776035967, -0.46606083776035967,
                      -0.46606083776035967, -30*D2R, -30*D2R, -30*D2R])

##Vehicle for running tests
class Test(unittest.TestCase):
    """Test cropping everything below the horizon works"""
    
    def make_comp_counter(self):
        """Here fill a comp counter with two SOURCEs, each with 3 components"""
        
        comp_counter = Component_Type_Counter()
        
        comp_counter.num_sources = 2
        comp_counter.total_comps = 6
        comp_counter.source_indexes = [0, 0, 0, 1, 1, 1]
        
        comp_counter.comp_ras = TEST_RAS
        comp_counter.comp_decs = TEST_DECS
        
        ##just stick these equal to the index for easy testig
        comp_counter.comp_types = np.arange(comp_counter.total_comps)
        comp_counter.num_list_fluxes = np.arange(comp_counter.total_comps)
        comp_counter.num_shape_coeffs = np.arange(comp_counter.total_comps)
        comp_counter.file_line_nums = np.arange(comp_counter.total_comps)
        
        return comp_counter
    
    def check_results(self, expec_orig_comp_indexes : np.ndarray,
                      comp_counter : Component_Type_Counter):
        """Given the expected indexes of the components that survive cropping,
        check the outputs are correct"""
        
        ##indexes of components above horizon are correct
        self.assertTrue(np.array_equal(expec_orig_comp_indexes,
                                       comp_counter.orig_comp_indexes))
        
        ##correct RA/Dec recovered
        self.assertTrue(np.array_equal(TEST_RAS[expec_orig_comp_indexes],
                                        comp_counter.comp_ras))
        self.assertTrue(np.array_equal(TEST_DECS[expec_orig_comp_indexes],
                                        comp_counter.comp_decs))
        
        ##We stuck these array equal to index, so should just equal expec indexes
        self.assertTrue(np.array_equal(expec_orig_comp_indexes,
                                       comp_counter.comp_types))
        self.assertTrue(np.array_equal(expec_orig_comp_indexes,
                                       comp_counter.num_list_fluxes))
        self.assertTrue(np.array_equal(expec_orig_comp_indexes,
                                       comp_counter.num_shape_coeffs))
        self.assertTrue(np.array_equal(expec_orig_comp_indexes,
                                       comp_counter.file_line_nums))
        
    ##For LST0, all COMPONENTs of the first SOURCE should have been retained.
    ##Only the first COMPONENT of the second SOURCE is above the horizon, so
    ##if cropping by SOURCE there should only be three COMPONENTs retained.
    ##If cropping by COMPONENT, there should be a 4th component
    
    def test_LST0_crop_by_source(self):
        comp_counter = self.make_comp_counter()
        
        crop_below_horizon(LST0, LATITUDE, comp_counter,
                           crop_by_component=False)
        
        ## expected indexes of components above the horizon
        expec_orig_comp_indexes = np.arange(3)
        
        self.check_results(expec_orig_comp_indexes, comp_counter)
        
        
    def test_LST0_crop_by_component(self):
        comp_counter = self.make_comp_counter()
        
        crop_below_horizon(LST0, LATITUDE, comp_counter,
                           crop_by_component=True)
        
        ## expected indexes of components above the horizon
        expec_orig_comp_indexes = np.arange(4)
        self.check_results(expec_orig_comp_indexes, comp_counter)
        
    ##For LST1, both SOURCEs should be retained
    def test_LST1_crop_by_source(self):
        comp_counter = self.make_comp_counter()
        
        crop_below_horizon(LST1, LATITUDE, comp_counter,
                           crop_by_component=False)
        
        ## expected indexes of components above the horizon
        expec_orig_comp_indexes = np.arange(6)
        self.check_results(expec_orig_comp_indexes, comp_counter)
        
    def test_LST1_crop_by_component(self):
        comp_counter = self.make_comp_counter()
        
        crop_below_horizon(LST1, LATITUDE, comp_counter,
                           crop_by_component=True)
        
        ## expected indexes of components above the horizon
        expec_orig_comp_indexes = np.arange(6)
        self.check_results(expec_orig_comp_indexes, comp_counter)
    
    ##For LST2, first source should have been cropped, second retained
    def test_LST2_crop_by_source(self):
        comp_counter = self.make_comp_counter()
        
        crop_below_horizon(LST2, LATITUDE, comp_counter,
                           crop_by_component=False)
        
        ## expected indexes of components above the horizon
        expec_orig_comp_indexes = np.arange(3, 6)
        self.check_results(expec_orig_comp_indexes, comp_counter)
        
    def test_LST2_crop_by_component(self):
        comp_counter = self.make_comp_counter()
        
        crop_below_horizon(LST2, LATITUDE, comp_counter,
                           crop_by_component=True)
        
        ## expected indexes of components above the horizon
        expec_orig_comp_indexes = np.arange(3, 6)
        self.check_results(expec_orig_comp_indexes, comp_counter)
    
    ##For LST3, everything should be below the horizon
    def test_LST3_crop_by_source(self):
        comp_counter = self.make_comp_counter()
        
        crop_below_horizon(LST3, LATITUDE, comp_counter,
                           crop_by_component=False)
        
        ## expected indexes of components above the horizon
        expec_orig_comp_indexes = np.arange(0)
        self.check_results(expec_orig_comp_indexes, comp_counter)
        
    def test_LST3_crop_by_component(self):
        comp_counter = self.make_comp_counter()
        
        crop_below_horizon(LST3, LATITUDE, comp_counter,
                           crop_by_component=True)
        
        ## expected indexes of components above the horizon
        expec_orig_comp_indexes = np.arange(0)
        self.check_results(expec_orig_comp_indexes, comp_counter)
    


##Run the test
if __name__ == '__main__':
    unittest.main()
    