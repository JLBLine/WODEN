"""Test `wodenpy.skymodel.read_fits_skymodel.read_fits_radec_count_components`,
which should read in and index different types of components from a FITS
sky model. First set of tests here use existing FITS sky models to test
against expected outcomes. Second set generate FITS sky models on the fly,
and then compares them to similarly generated expected outcomes."""

from sys import path
import os
import unittest
import numpy as np

##This is where our code lives
# code_dir = os.environ['CMAKE_CURRENT_SOURCE_DIR']

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

##Append the location of run_woden.py to the sys.path to import it
path.append('{:s}/../../../wodenpy/skymodel'.format(code_dir))



# ##Code we are testing
from wodenpy.skymodel import read_fits_skymodel
import fits_skymodel_common
from read_skymodel_common import Skymodel_Settings, make_expected_comp_counter, check_comp_counter
# import wodenpy
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes



D2R = np.pi/180.0

##Expected values based on the test sky models
POINT0_RA = 15*D2R
POINT0_DEC = -30.0*D2R
POINT0_NUM_LIST_ENTRIES = 4

GAUSS0_RA = 15*2.0*D2R
GAUSS0_DEC = -30.0*D2R
GAUSS0_NUM_LIST_ENTRIES = 3

SHAPE0_RA = 15*3.0*D2R
SHAPE0_DEC = 20.0*D2R
SHAPE0_NUM_COEFFS = 4
SHAPE0_NUM_LIST_ENTRIES = 5


class Expected_Total_Nums(object):
    """
    Something to hold all the summed component type expected values. Initialise
    with everything equal to zero, unless explicitly set
    """
    def __init__(self, num_sources = 0, total_comps = 0, total_point_comps = 0,
                       total_gauss_comps = 0, total_shape_comps = 0, 
                       total_shape_basis = 0, total_point_list_fluxes = 0,
                       total_gauss_list_fluxes = 0, total_shape_list_fluxes = 0,
                       num_point_flux_powers = 0, num_point_flux_curves = 0,
                       num_point_flux_lists = 0, num_gauss_flux_powers = 0,
                       num_gauss_flux_curves = 0, num_gauss_flux_lists = 0,
                       num_shape_flux_powers = 0, num_shape_flux_curves = 0,
                       num_shape_flux_lists = 0) -> None:

        self.num_sources = num_sources
        self.total_comps = total_comps
        self.total_point_comps = total_point_comps
        self.total_gauss_comps = total_gauss_comps
        self.total_shape_comps = total_shape_comps
        self.total_shape_basis = total_shape_basis
        # self.total_point_list_fluxes = total_point_list_fluxes
        # self.total_gauss_list_fluxes = total_gauss_list_fluxes
        # self.total_shape_list_fluxes = total_shape_list_fluxes

        self.num_point_flux_powers = num_point_flux_powers
        self.num_point_flux_curves = num_point_flux_curves
        self.num_point_flux_lists = num_point_flux_lists
        self.num_gauss_flux_powers = num_gauss_flux_powers
        self.num_gauss_flux_curves = num_gauss_flux_curves
        self.num_gauss_flux_lists = num_gauss_flux_lists
        self.num_shape_flux_powers = num_shape_flux_powers
        self.num_shape_flux_curves = num_shape_flux_curves
        self.num_shape_flux_lists = num_shape_flux_lists
        
class Expected_Component_Nums(object):
    """
    Something to hold the expected values for a single component. Initialise
    `comp_type`, and everything else equal to zero, unless explicitly set
    """
    def __init__(self, comp_type : int, num_list_flux=0,
                 num_shape_basis=0) -> None:
        
        self.comp_type = comp_type
        self.comp_ra = None
        self.comp_dec = None
        self.num_list_flux = num_list_flux
        self.num_shape_basis = num_shape_basis

##Vehicle for running tests
class Test(unittest.TestCase):
    """"""

    def check_total_numbers(self, comp_counter : Component_Type_Counter,
                            expec_nums : Expected_Total_Nums):
        """Check that the overall counts of component types are as expected"""
        # print("num_sources", comp_counter.num_sources)
        # print("total_comps", comp_counter.total_comps)
        # print("total_point_comps", comp_counter.total_point_comps)
        # print("total_gauss_comps", comp_counter.total_gauss_comps)
        # print("total_shape_comps", comp_counter.total_shape_comps)
        # print("num_shape_basis", comp_counter.num_shape_basis)
        
        # print("num_point_list_fluxes", comp_counter.num_point_list_fluxes)
        # print("num_gauss_list_fluxes", comp_counter.num_gauss_list_fluxes)
        # print("num_shape_list_fluxes", comp_counter.num_shape_list_fluxes)
        
        self.assertEqual(expec_nums.num_sources,
                         comp_counter.num_sources)
        self.assertEqual(expec_nums.total_comps,
                         comp_counter.total_comps)
        self.assertEqual(expec_nums.total_point_comps,
                         comp_counter.total_point_comps)
        self.assertEqual(expec_nums.total_gauss_comps,
                         comp_counter.total_gauss_comps)
        self.assertEqual(expec_nums.total_shape_comps,
                         comp_counter.total_shape_comps)

        self.assertEqual(expec_nums.total_shape_basis,
                         comp_counter.total_shape_basis)

        self.assertEqual(expec_nums.num_point_flux_powers,
                         comp_counter.num_point_flux_powers)
        self.assertEqual(expec_nums.num_point_flux_curves,
                         comp_counter.num_point_flux_curves)
        self.assertEqual(expec_nums.num_point_flux_lists,
                         comp_counter.num_point_flux_lists)
        self.assertEqual(expec_nums.num_gauss_flux_powers,
                         comp_counter.num_gauss_flux_powers)
        self.assertEqual(expec_nums.num_gauss_flux_curves,
                         comp_counter.num_gauss_flux_curves)
        self.assertEqual(expec_nums.num_gauss_flux_lists,
                         comp_counter.num_gauss_flux_lists)
        self.assertEqual(expec_nums.num_shape_flux_powers,
                         comp_counter.num_shape_flux_powers)
        self.assertEqual(expec_nums.num_shape_flux_curves,
                         comp_counter.num_shape_flux_curves)
        self.assertEqual(expec_nums.num_shape_flux_lists,
                         comp_counter.num_shape_flux_lists)
        
    def check_single_component(self, comp_counter : Component_Type_Counter,
                               expec_comp : Expected_Component_Nums,
                               comp_ind : int):
        
        point_vals = [CompTypes.POINT_POWER.value, CompTypes.POINT_CURVE.value,
                      CompTypes.POINT_LIST.value]
        gauss_vals = [CompTypes.GAUSS_POWER.value, CompTypes.GAUSS_CURVE.value,
                      CompTypes.GAUSS_LIST.value]
        shape_vals = [CompTypes.SHAPE_POWER.value, CompTypes.SHAPE_CURVE.value,
                      CompTypes.SHAPE_LIST.value]
    
        if expec_comp.comp_type in point_vals:
            expec_comp.comp_ra = POINT0_RA
            expec_comp.comp_dec = POINT0_DEC
            
        elif expec_comp.comp_type in gauss_vals:
            expec_comp.comp_ra = GAUSS0_RA
            expec_comp.comp_dec = GAUSS0_DEC
            
        elif expec_comp.comp_type in shape_vals:
            expec_comp.comp_ra = SHAPE0_RA
            expec_comp.comp_dec = SHAPE0_DEC
            
        # print("comp_ras", comp_counter.comp_ras)
        # print("comp_decs", comp_counter.comp_decs)
            
        self.assertEqual(expec_comp.comp_ra,
                         comp_counter.comp_ras[comp_ind])
        self.assertEqual(expec_comp.comp_dec,
                         comp_counter.comp_decs[comp_ind])
        self.assertEqual(expec_comp.comp_type,
                         comp_counter.comp_types[comp_ind])
        self.assertEqual(expec_comp.num_list_flux,
                         comp_counter.num_list_fluxes[comp_ind])
        self.assertEqual(expec_comp.num_shape_basis,
                         comp_counter.num_shape_coeffs[comp_ind])
        
        # print(expec_comp.num_list_flux, comp_counter.num_list_fluxes[comp_ind])
        
    def run_and_check_tots_read_fits_radec_count_components(self,
                                    skymodel : str,
                                    expec_nums : Expected_Total_Nums,
                                    expec_source_inds : np.ndarray,
                                    expec_file_nums : np.ndarray):
        """Runs the code we are testing, and checks the total overall number
        counts, and that the component to source mapping index list is correct.
        Returns the comp counter so individual components can be checked"""
        
        ##code we are testing
        comp_counter = read_fits_skymodel.read_fits_radec_count_components(skymodel)
        
        ##Check the total numbers read in are correct
        self.check_total_numbers(comp_counter, expec_nums)
        
        ##only needed for text-style sky models
        ##check the source indexes (that map components to a source)
        ##comp_counter.source_indexes is initialised to a certain size and
        ##only filled up to number of components, so only test up to that point
        # print(expec_source_inds, comp_counter.source_indexes[:comp_counter.total_comps])
        # self.assertTrue((expec_source_inds == comp_counter.source_indexes[:comp_counter.total_comps]).all())
        
        # ##similar test for what line number each component starts at in the file
        # self.assertTrue((expec_file_nums == comp_counter.file_line_nums[:comp_counter.total_comps]).all())
        
        return comp_counter

    ##POWER LAW TESTS-----------------------------------------------------------------
    ##--------------------------------------------------------------------------------
        
    def test_SinglePointPower(self):
        """srclist that contains a single point source with power law"""
        
        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_singlepoint_power.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1, total_comps = 1,
                                         total_point_comps = 1,
                                         num_point_flux_powers = 1)

        expec_source_inds = np.arange(1)
        expec_file_nums = np.arange(1,2)
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)

        ##check the individual components are correct
        expec_comp = Expected_Component_Nums(CompTypes.POINT_POWER.value)
        self.check_single_component(comp_counter, expec_comp, 0)
        
        
        
    def test_SingleGaussianPower(self):
        """srclist that contains a single gaussian with a power law"""
        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_singlegauss_power.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1, total_comps = 1,
                                         total_gauss_comps = 1,
                                         num_gauss_flux_powers = 1)
        expec_source_inds = np.arange(1)
        expec_file_nums = np.arange(1,2)
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        ##check the individual components are correct
        expec_comp = Expected_Component_Nums(CompTypes.GAUSS_POWER.value)
        self.check_single_component(comp_counter, expec_comp, 0)
        
    def test_SingleShapeletPower(self):
        """srclist that contains a single shapelet with a power law"""
        
        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_singleshape_power.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1, total_comps = 1,
                                         total_shape_comps = 1,
                                         total_shape_basis = SHAPE0_NUM_COEFFS,
                                         num_shape_flux_powers = 1)
        expec_source_inds = np.arange(1)
        expec_file_nums = np.arange(1,2)
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        ##check the individual components are correct
        expec_comp = Expected_Component_Nums(CompTypes.SHAPE_POWER.value,
                                             num_shape_basis = SHAPE0_NUM_COEFFS)
        
        self.check_single_component(comp_counter, expec_comp, 0)
        
        
    def test_ThreeSourcesPower(self):
        """srclist contains three SOURCES, first a point, second a gaussian,
        third a shapelet, all power laws"""

        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_threesources_power.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 3,
                                         total_comps = 3,
                                         total_point_comps = 1,
                                         num_point_flux_powers = 1,
                                         total_gauss_comps = 1,
                                         num_gauss_flux_powers = 1,
                                         total_shape_comps = 1,
                                         num_shape_flux_powers = 1,
                                         total_shape_basis = SHAPE0_NUM_COEFFS)

        expec_source_inds = np.arange(3)
        expec_file_nums = np.array([1, 16, 33])
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        ##All the individual expected components in a list
        expec_comps = [Expected_Component_Nums(CompTypes.POINT_POWER.value),
                       Expected_Component_Nums(CompTypes.GAUSS_POWER.value),
                       Expected_Component_Nums(CompTypes.SHAPE_POWER.value,
                                       num_shape_basis = SHAPE0_NUM_COEFFS)]

        ##iterate over component details and check they are correct
        for comp_ind, expec_comp in enumerate(expec_comps):
            self.check_single_component(comp_counter, expec_comp, comp_ind)

    def test_ThreeComponentsPower(self):
        """srclist contains one SOURCE with three COMPONENTS, first a point,
        second a gaussian, third a shapelet, all power laws"""

        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_threecomponents_power.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1,
                                         total_comps = 3,
                                         total_point_comps = 1,
                                         num_point_flux_powers = 1,
                                         total_gauss_comps = 1,
                                         num_gauss_flux_powers = 1,
                                         total_shape_comps = 1,
                                         num_shape_flux_powers = 1,
                                         total_shape_basis = SHAPE0_NUM_COEFFS)

        expec_source_inds = np.zeros(3)
        expec_file_nums = np.array([1, 15, 31])
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        ##All the individual expected components in a list
        expec_comps = [Expected_Component_Nums(CompTypes.POINT_POWER.value),
                       Expected_Component_Nums(CompTypes.GAUSS_POWER.value),
                       Expected_Component_Nums(CompTypes.SHAPE_POWER.value,
                                       num_shape_basis = SHAPE0_NUM_COEFFS)]

        ##iterate over component details and check they are correct
        for comp_ind, expec_comp in enumerate(expec_comps):
            self.check_single_component(comp_counter, expec_comp, comp_ind)

    # ##CURVED POWER LAW TESTS----------------------------------------------------------
    # ##--------------------------------------------------------------------------------
        
    def test_SinglePointCurved(self):
        """srclist that contains a single point source with power law"""
        
        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_singlepoint_curve.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1, total_comps = 1,
                                         total_point_comps = 1,
                                         num_point_flux_curves = 1)

        expec_source_inds = np.arange(1)
        expec_file_nums = np.arange(1,2)
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)

        ##check the individual components are correct
        expec_comp = Expected_Component_Nums(CompTypes.POINT_CURVE.value)
        self.check_single_component(comp_counter, expec_comp, 0)
        
    def test_SingleGaussianCurved(self):
        """srclist that contains a single gaussian with a power law"""
        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_singlegauss_curve.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1, total_comps = 1,
                                         total_gauss_comps = 1,
                                         num_gauss_flux_curves = 1)
        expec_source_inds = np.arange(1)
        expec_file_nums = np.arange(1,2)
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        ##check the individual components are correct
        expec_comp = Expected_Component_Nums(CompTypes.GAUSS_CURVE.value)
        self.check_single_component(comp_counter, expec_comp, 0)
        
    def test_SingleShapeletCurved(self):
        """srclist that contains a single shapelet with a power law"""
        
        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_singleshape_curve.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1, total_comps = 1,
                                         total_shape_comps = 1,
                                         total_shape_basis = SHAPE0_NUM_COEFFS,
                                         num_shape_flux_curves = 1)
        expec_source_inds = np.arange(1)
        expec_file_nums = np.arange(1,2)
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        ##check the individual components are correct
        expec_comp = Expected_Component_Nums(CompTypes.SHAPE_CURVE.value,
                                             num_shape_basis = SHAPE0_NUM_COEFFS)
        
        self.check_single_component(comp_counter, expec_comp, 0)
        
        
    def test_ThreeSourcesCurved(self):
        """srclist contains three SOURCES, first a point, second a gaussian,
        third a shapelet, all power laws"""

        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_threesources_curve.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 3,
                                         total_comps = 3,
                                         total_point_comps = 1,
                                         num_point_flux_curves = 1,
                                         total_gauss_comps = 1,
                                         num_gauss_flux_curves = 1,
                                         total_shape_comps = 1,
                                         num_shape_flux_curves = 1,
                                         total_shape_basis = SHAPE0_NUM_COEFFS)
        expec_source_inds = np.arange(3)
        expec_file_nums = np.array([1, 15, 33])
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)

        ##All the individual expected components in a list
        expec_comps = [Expected_Component_Nums(CompTypes.POINT_CURVE.value),
                       Expected_Component_Nums(CompTypes.GAUSS_CURVE.value),
                       Expected_Component_Nums(CompTypes.SHAPE_CURVE.value,
                                        num_shape_basis = SHAPE0_NUM_COEFFS)]

        ##iterate over component details and check they are correct
        for comp_ind, expec_comp in enumerate(expec_comps):
            self.check_single_component(comp_counter, expec_comp, comp_ind)

    def test_ThreeComponentsCurved(self):
        """srclist contains one SOURCE with three COMPONENTS, first a point,
        second a gaussian, third a shapelet, all power laws"""

        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_threecomponents_curve.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1,
                                         total_comps = 3,
                                         total_point_comps = 1,
                                         num_point_flux_curves = 1,
                                         total_gauss_comps = 1,
                                         num_gauss_flux_curves = 1,
                                         total_shape_comps = 1,
                                         num_shape_flux_curves = 1,
                                         total_shape_basis = SHAPE0_NUM_COEFFS)
        expec_source_inds = np.zeros(3)
        expec_file_nums = np.array([1, 14, 31])
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        ##All the individual expected components in a list
        expec_comps = [Expected_Component_Nums(CompTypes.POINT_CURVE.value),
                       Expected_Component_Nums(CompTypes.GAUSS_CURVE.value),
                       Expected_Component_Nums(CompTypes.SHAPE_CURVE.value,
                                        num_shape_basis = SHAPE0_NUM_COEFFS)]

        ##iterate over component details and check they are correct
        for comp_ind, expec_comp in enumerate(expec_comps):
            self.check_single_component(comp_counter, expec_comp, comp_ind)


    # ##LIST FLUX TESTS-----------------------------------------------------------------
    # ##--------------------------------------------------------------------------------
        
    def test_SinglePointList(self):
        """srclist that contains a single point source with flux list"""
        
        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_singlepoint_list.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1, total_comps = 1,
                                         total_point_comps = 1,
                                         num_point_flux_lists = 1,
                                         total_point_list_fluxes = POINT0_NUM_LIST_ENTRIES)
        
        expec_source_inds = np.arange(1)
        expec_file_nums = np.arange(1,2)
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)

        ##check the individual components are correct
        expec_comp = Expected_Component_Nums(CompTypes.POINT_LIST.value,
                                             num_list_flux = POINT0_NUM_LIST_ENTRIES)
        self.check_single_component(comp_counter, expec_comp, 0)
        
    def test_SingleGaussianList(self):
        """srclist that contains a single gaussian with a flux list"""
        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_singlegauss_list.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1, total_comps = 1,
                                         total_gauss_comps = 1,
                                         num_gauss_flux_lists = 1,
                                         total_gauss_list_fluxes = GAUSS0_NUM_LIST_ENTRIES)
        expec_source_inds = np.arange(1)
        expec_file_nums = np.arange(1,2)
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        ##check the individual components are correct
        expec_comp = Expected_Component_Nums(CompTypes.GAUSS_LIST.value,
                                             num_list_flux = GAUSS0_NUM_LIST_ENTRIES)
        self.check_single_component(comp_counter, expec_comp, 0)
        
    def test_SingleShapeletList(self):
        """srclist that contains a single shapelet with a flux list"""
        
        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_singleshape_list.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1, total_comps = 1,
                                         total_shape_comps = 1,
                                         total_shape_basis = SHAPE0_NUM_COEFFS,
                                         num_shape_flux_lists = 1,
                                         total_shape_list_fluxes = SHAPE0_NUM_LIST_ENTRIES)
        expec_source_inds = np.arange(1)
        expec_file_nums = np.arange(1,2)
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        ##check the individual components are correct
        expec_comp = Expected_Component_Nums(CompTypes.SHAPE_LIST.value,
                                             num_shape_basis = SHAPE0_NUM_COEFFS,
                                             num_list_flux = SHAPE0_NUM_LIST_ENTRIES)
        
        self.check_single_component(comp_counter, expec_comp, 0)
        
        
    def test_ThreeSourcesList(self):
        """srclist contains three SOURCES, first a point, second a gaussian,
        third a shapelet, all flux lists"""

        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_threesources_list.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 3,
                                         total_comps = 3,
                                         total_point_comps = 1,
                                         num_point_flux_lists = 1,
                                         total_point_list_fluxes = POINT0_NUM_LIST_ENTRIES,
                                         total_gauss_comps = 1,
                                         num_gauss_flux_lists = 1,
                                         total_gauss_list_fluxes = GAUSS0_NUM_LIST_ENTRIES,
                                         total_shape_comps = 1,
                                         num_shape_flux_lists = 1,
                                         total_shape_list_fluxes = SHAPE0_NUM_LIST_ENTRIES,
                                         total_shape_basis = SHAPE0_NUM_COEFFS)
        expec_source_inds = np.arange(3)
        expec_file_nums = np.array([1, 21, 42])
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        
        ##All the individual expected components in a list
        expec_comps = [Expected_Component_Nums(CompTypes.POINT_LIST.value,
                                               num_list_flux = POINT0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.GAUSS_LIST.value,
                                               num_list_flux = GAUSS0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.SHAPE_LIST.value,
                                               num_shape_basis = SHAPE0_NUM_COEFFS,
                                               num_list_flux = SHAPE0_NUM_LIST_ENTRIES)]

        ##iterate over component details and check they are correct
        for comp_ind, expec_comp in enumerate(expec_comps):
            self.check_single_component(comp_counter, expec_comp, comp_ind)

    def test_ThreeComponentsList(self):
        """srclist contains one SOURCE with three COMPONENTS, first a point,
        second a gaussian, third a shapelet, all flux lists"""

        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_threecomponents_list.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 1,
                                         total_comps = 3,
                                         total_point_comps = 1,
                                         num_point_flux_lists = 1,
                                         total_point_list_fluxes = POINT0_NUM_LIST_ENTRIES,
                                         total_gauss_comps = 1,
                                         num_gauss_flux_lists = 1,
                                         total_gauss_list_fluxes = GAUSS0_NUM_LIST_ENTRIES,
                                         total_shape_comps = 1,
                                         num_shape_flux_lists = 1,
                                         total_shape_list_fluxes = SHAPE0_NUM_LIST_ENTRIES,
                                         total_shape_basis = SHAPE0_NUM_COEFFS)
        expec_source_inds = np.zeros(3)
        expec_file_nums = np.array([1, 20, 40])
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        ##All the individual expected components in a list
        expec_comps = [Expected_Component_Nums(CompTypes.POINT_LIST.value,
                                               num_list_flux = POINT0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.GAUSS_LIST.value,
                                               num_list_flux = GAUSS0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.SHAPE_LIST.value,
                                               num_shape_basis = SHAPE0_NUM_COEFFS,
                                               num_list_flux = SHAPE0_NUM_LIST_ENTRIES)]
        
        #iterate over component details and check they are correct
        for comp_ind, expec_comp in enumerate(expec_comps):
            self.check_single_component(comp_counter, expec_comp, comp_ind)


    def test_MultiSourceComponents(self):
        """srclist that contains many combinations of things, in this order:
        SOURCE 0 - singlepoint_power
        SOURCE 2 - singlepoint_list
        SOURCE 3 - singlepoint_curve
        SOURCE 4 - singlegauss_power
        SOURCE 5 - singlegauss_curve
        SOURCE 6 - singlegauss_list
        SOURCE 7 - singleshapelet_curve
        SOURCE 8 - singleshapelet_list
        SOURCE 9 - singleshapelet_power
        SOURCE 10 - all of those single components above appear twice as
                    components in a somewhat random order
        """

        ##sky model used for testing
        skymodel = '{:s}/fits_models/srclist_mulitple_source-components.fits'.format(code_dir)
        
        ##Check the total numbers read in are correct
        expec_nums = Expected_Total_Nums(num_sources = 10,
                                         total_comps = 27,
                                         total_point_comps = 9,
                                         num_point_flux_powers = 3,
                                         num_point_flux_curves = 3,
                                         num_point_flux_lists = 3,
                                         total_point_list_fluxes = 3*POINT0_NUM_LIST_ENTRIES,
                                         total_gauss_comps = 9,
                                         num_gauss_flux_powers = 3,
                                         num_gauss_flux_curves = 3,
                                         num_gauss_flux_lists = 3,
                                         total_gauss_list_fluxes = 3*GAUSS0_NUM_LIST_ENTRIES,
                                         total_shape_comps = 9,
                                         num_shape_flux_powers = 3,
                                         num_shape_flux_curves = 3,
                                         num_shape_flux_lists = 3,
                                         total_shape_list_fluxes = 3*SHAPE0_NUM_LIST_ENTRIES,
                                         total_shape_basis = 9*SHAPE0_NUM_COEFFS)
        
        expec_source_inds = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8,
                                      9, 9, 9, 9, 9, 9, 9, 9, 9,
                                      9, 9, 9, 9, 9, 9, 9, 9, 9])
    
        expec_file_nums = np.array([1, 14, 34, 48, 65,83, 104, 135, 174,
                                    204, 216, 235, 248, 264, 284, 314, 352,
                                    381, 393, 412, 425, 441, 461, 491, 529,
                                    558, 575])
        comp_counter = self.run_and_check_tots_read_fits_radec_count_components(skymodel,
                                    expec_nums, expec_source_inds, expec_file_nums)
        
        ##All the individual expected components in a list
        expec_comps = [Expected_Component_Nums(CompTypes.POINT_POWER.value),
                       Expected_Component_Nums(CompTypes.POINT_LIST.value,
                                               num_list_flux = POINT0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.POINT_CURVE.value),
                       Expected_Component_Nums(CompTypes.GAUSS_POWER.value),
                       Expected_Component_Nums(CompTypes.GAUSS_CURVE.value),
                       Expected_Component_Nums(CompTypes.GAUSS_LIST.value,
                                               num_list_flux = GAUSS0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.SHAPE_CURVE.value,
                                               num_shape_basis = SHAPE0_NUM_COEFFS),
                       Expected_Component_Nums(CompTypes.SHAPE_LIST.value,
                                               num_shape_basis = SHAPE0_NUM_COEFFS,
                                               num_list_flux = SHAPE0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.SHAPE_POWER.value,
                                               num_shape_basis = SHAPE0_NUM_COEFFS),
                       Expected_Component_Nums(CompTypes.POINT_POWER.value),
                       Expected_Component_Nums(CompTypes.POINT_LIST.value,
                                               num_list_flux = POINT0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.POINT_CURVE.value),
                       Expected_Component_Nums(CompTypes.GAUSS_POWER.value),
                       Expected_Component_Nums(CompTypes.GAUSS_LIST.value,
                                               num_list_flux = GAUSS0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.SHAPE_CURVE.value,
                                               num_shape_basis = SHAPE0_NUM_COEFFS),
                       Expected_Component_Nums(CompTypes.SHAPE_LIST.value,
                                               num_shape_basis = SHAPE0_NUM_COEFFS,
                                               num_list_flux = SHAPE0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.SHAPE_POWER.value,
                                               num_shape_basis = SHAPE0_NUM_COEFFS),
                       Expected_Component_Nums(CompTypes.POINT_POWER.value),
                       Expected_Component_Nums(CompTypes.POINT_LIST.value,
                                               num_list_flux = POINT0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.POINT_CURVE.value),
                       Expected_Component_Nums(CompTypes.GAUSS_POWER.value),
                       Expected_Component_Nums(CompTypes.GAUSS_LIST.value,
                                               num_list_flux = GAUSS0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.SHAPE_CURVE.value,
                                               num_shape_basis = SHAPE0_NUM_COEFFS),
                       Expected_Component_Nums(CompTypes.SHAPE_LIST.value,
                                               num_shape_basis = SHAPE0_NUM_COEFFS,
                                               num_list_flux = SHAPE0_NUM_LIST_ENTRIES),
                       Expected_Component_Nums(CompTypes.SHAPE_POWER.value,
                                               num_shape_basis = SHAPE0_NUM_COEFFS),
                       Expected_Component_Nums(CompTypes.GAUSS_CURVE.value),
                       Expected_Component_Nums(CompTypes.GAUSS_CURVE.value)]
                       
    
        #iterate over component details and check they are correct
        for comp_ind, expec_comp in enumerate(expec_comps):
            self.check_single_component(comp_counter, expec_comp, comp_ind)
        
        
    def test_MissingFile(self):
        """If the file doesn't exist, we should has a sys exit so check things
        fail"""
        
        skymodel = "is-not-a-thing-fool.fits"
        
        with self.assertRaises(SystemExit) as cm:
            ##code we are testing
            comp_counter = read_fits_skymodel.read_fits_radec_count_components(skymodel)
            
    def test_StokesV_fracpol(self):
        """
        """
        
        settings = Skymodel_Settings(deg_between_comps = 240,
                                     num_coeff_per_shape = 4,
                                     num_list_values = 4,
                                     comps_per_source = 10,
                                     stokesV_frac_cadence = 5)
        
        fits_skymodel_common.write_full_test_skymodel_fits(settings)
        
        comp_counter = read_fits_skymodel.read_fits_radec_count_components('test_full_skymodel.fits')
        
        comp_counter.print_info()
        
        ##Test we got what we expected
        check_comp_counter(comp_counter, settings)
        
    def test_StokesV_power(self):
        """
        """
        
        settings = Skymodel_Settings(deg_between_comps = 240,
                                     num_coeff_per_shape = 4,
                                     num_list_values = 4,
                                     comps_per_source = 10,
                                     stokesV_pl_cadence = 4)
        
        fits_skymodel_common.write_full_test_skymodel_fits(settings)
        
        comp_counter = read_fits_skymodel.read_fits_radec_count_components('test_full_skymodel.fits')
        
        comp_counter.print_info()
        
        ##Test we got what we expected
        check_comp_counter(comp_counter, settings)
        
    def test_StokesV_curve(self):
        """
        """
        
        settings = Skymodel_Settings(deg_between_comps = 240,
                                     num_coeff_per_shape = 4,
                                     num_list_values = 4,
                                     comps_per_source = 10,
                                     stokesV_cpl_cadence = 3)
        
        fits_skymodel_common.write_full_test_skymodel_fits(settings)
        
        comp_counter = read_fits_skymodel.read_fits_radec_count_components('test_full_skymodel.fits')
        
        comp_counter.print_info()
        
        ##Test we got what we expected
        check_comp_counter(comp_counter, settings)
        
    def test_StokesV_list(self):
        """
        """
        
        settings = Skymodel_Settings(deg_between_comps = 240,
                                     num_coeff_per_shape = 4,
                                     num_list_values = 4,
                                     comps_per_source = 10,
                                     stokesV_list_cadence = 5,
                                     stokesV_num_list = 3)
        
        fits_skymodel_common.write_full_test_skymodel_fits(settings)
        
        comp_counter = read_fits_skymodel.read_fits_radec_count_components('test_full_skymodel.fits')
        
        comp_counter.print_info()
        
        ##Test we got what we expected
        check_comp_counter(comp_counter, settings)
        
        
    def test_linpol_fracpol(self):
        """
        """
        
        settings = Skymodel_Settings(deg_between_comps = 240,
                                     num_coeff_per_shape = 4,
                                     num_list_values = 4,
                                     comps_per_source = 10,
                                     linpol_frac_cadence = 5)
        
        fits_skymodel_common.write_full_test_skymodel_fits(settings)
        
        comp_counter = read_fits_skymodel.read_fits_radec_count_components('test_full_skymodel.fits')
        
        comp_counter.print_info()
        
        ##Test we got what we expected
        check_comp_counter(comp_counter, settings)
        
    def test_linpol_power(self):
        """
        """
        
        settings = Skymodel_Settings(deg_between_comps = 240,
                                     num_coeff_per_shape = 4,
                                     num_list_values = 4,
                                     comps_per_source = 10,
                                     linpol_pl_cadence = 4)
        
        fits_skymodel_common.write_full_test_skymodel_fits(settings)
        
        comp_counter = read_fits_skymodel.read_fits_radec_count_components('test_full_skymodel.fits')
        
        comp_counter.print_info()
        
        ##Test we got what we expected
        check_comp_counter(comp_counter, settings)
        
    def test_linpol_curve(self):
        """
        """
        
        settings = Skymodel_Settings(deg_between_comps = 240,
                                     num_coeff_per_shape = 4,
                                     num_list_values = 4,
                                     comps_per_source = 10,
                                     linpol_cpl_cadence = 3)
        
        fits_skymodel_common.write_full_test_skymodel_fits(settings)
        
        comp_counter = read_fits_skymodel.read_fits_radec_count_components('test_full_skymodel.fits')
        
        comp_counter.print_info()
        
        ##Test we got what we expected
        check_comp_counter(comp_counter, settings)
        
        
    def test_linpol_list(self):
        """
        """
        
        settings = Skymodel_Settings(deg_between_comps = 240,
                                     num_coeff_per_shape = 4,
                                     num_list_values = 4,
                                     comps_per_source = 10,
                                     linpol_list_cadence = 3,
                                     linpol_p_list_cadence = 4,
                                     linpol_num_list = 5,
                                     linpol_num_p_list = 4)
        
        fits_skymodel_common.write_full_test_skymodel_fits(settings)
        
        comp_counter = read_fits_skymodel.read_fits_radec_count_components('test_full_skymodel.fits')
        
        comp_counter.print_info()
        
        ##Test we got what we expected
        check_comp_counter(comp_counter, settings)
    
    
        
        
    def test_all_polarisations(self):
        """
        """
        
        settings = Skymodel_Settings(deg_between_comps = 50,
                                     num_coeff_per_shape = 4,
                                     num_list_values = 4,
                                     comps_per_source = 10,
                                     stokesV_frac_cadence = 5,
                                     stokesV_pl_cadence = 7,
                                     stokesV_cpl_cadence = 11,
                                     stokesV_list_cadence = 13,
                                     stokesV_num_list = 3,
                                     linpol_frac_cadence = 17,
                                     linpol_pl_cadence = 19,
                                     linpol_cpl_cadence = 23,
                                     linpol_list_cadence = 29,
                                     linpol_p_list_cadence = 31,
                                     linpol_num_list = 5,
                                     linpol_num_p_list = 4)
        
        fits_skymodel_common.write_full_test_skymodel_fits(settings)
        
        comp_counter = read_fits_skymodel.read_fits_radec_count_components('test_full_skymodel.fits')
        
        comp_counter.print_info()
        
        ##Test we got what we expected
        check_comp_counter(comp_counter, settings)
        
##Run the test
if __name__ == '__main__':
    unittest.main()
    