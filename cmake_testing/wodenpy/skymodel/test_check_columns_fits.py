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
from astropy.table import Column, Table
from astropy.io import fits

def write_borked_skymodel(unq_source_id=True, name=True, ra=True, dec=True, comp_type=True, major_dc=True,
                          minor_dc=True, pa_dc=True, mod_type=True, norm_comp_pl=True, alpha_pl=True,
                          norm_comp_cpl=True, alpha_cpl=True, curve_cpl=True, v_mod_type=True,
                          use_mod_types='pl', use_comp_types='P', use_v_mod_types='pl', use_lin_mod_types='pl', 
                          v_pol_frac=True, v_norm_comp_pl=True, v_alpha_pl=True, v_norm_comp_cpl=True,
                          v_alpha_cpl=True, v_curve_cpl=True, lin_mod_type=True, rm=True,
                          intr_pol_angle=True, lin_pol_frac=True, lin_norm_comp_pl=True,
                          lin_alpha_pl=True, lin_norm_comp_cpl=True, lin_alpha_cpl=True,
                          lin_curve_cpl=True, int_flux=True, shape_table=True,
                          s_names=True, s_n1s=True, s_n2s=True, s_coeffs=True, shape_max_n=10,
                          p_flux_table=True, p_flux_names=True, p_flux_fluxes=True,
                          q_flux_table=True, q_flux_names=True, q_flux_fluxes=True,
                          u_flux_table=True, u_flux_names=True, u_flux_fluxes=True,
                          v_flux_table=True, v_flux_names=True, v_flux_fluxes=True,
                          use_table_names=True):
    
    num_comps = 3
    main_cols = []
    
    if unq_source_id:
        this_col = Column(data=['bleh']*3, name="UNQ_SOURCE_ID")
        main_cols.append(this_col)

    if name:
        this_col = Column(data=[f'bleh_C{i:2d}' for i in range(num_comps)], name="NAME")
        main_cols.append(this_col)

    if ra:
        this_col = Column(data=np.zeros(num_comps), name="RA")
        main_cols.append(this_col)

    if dec:
        this_col = Column(data=np.zeros(num_comps), name="DEC")
        main_cols.append(this_col)

    if comp_type:
        this_col = Column(data=[use_comp_types]*3, name="COMP_TYPE")
        main_cols.append(this_col)

    if major_dc:
        this_col = Column(data=np.zeros(num_comps), name="MAJOR_DC")
        main_cols.append(this_col)

    if minor_dc:
        this_col = Column(data=np.zeros(num_comps), name="MINOR_DC")
        main_cols.append(this_col)

    if pa_dc:
        this_col = Column(data=np.zeros(num_comps), name="PA_DC")
        main_cols.append(this_col)

    if mod_type:
        this_col = Column(data=[use_mod_types]*num_comps, name="MOD_TYPE")
        main_cols.append(this_col)

    if norm_comp_pl:
        this_col = Column(data=np.zeros(num_comps), name="NORM_COMP_PL")
        main_cols.append(this_col)

    if alpha_pl:
        this_col = Column(data=np.zeros(num_comps), name="ALPHA_PL")
        main_cols.append(this_col)

    if norm_comp_cpl:
        this_col = Column(data=np.zeros(num_comps), name="NORM_COMP_CPL")
        main_cols.append(this_col)

    if alpha_cpl:
        this_col = Column(data=np.zeros(num_comps), name="ALPHA_CPL")
        main_cols.append(this_col)

    if curve_cpl:
        this_col = Column(data=np.zeros(num_comps), name="CURVE_CPL")
        main_cols.append(this_col)

    if v_mod_type:
        this_col = Column(data=[use_v_mod_types]*num_comps, name="V_MOD_TYPE")
        main_cols.append(this_col)

    if v_pol_frac:
        this_col = Column(data=np.zeros(num_comps), name="V_POL_FRAC")
        main_cols.append(this_col)

    if v_norm_comp_pl:
        this_col = Column(data=np.zeros(num_comps), name="V_NORM_COMP_PL")
        main_cols.append(this_col)

    if v_alpha_pl:
        this_col = Column(data=np.zeros(num_comps), name="V_ALPHA_PL")
        main_cols.append(this_col)

    if v_norm_comp_cpl:
        this_col = Column(data=np.zeros(num_comps), name="V_NORM_COMP_CPL")
        main_cols.append(this_col)

    if v_alpha_cpl:
        this_col = Column(data=np.zeros(num_comps), name="V_ALPHA_CPL")
        main_cols.append(this_col)

    if v_curve_cpl:
        this_col = Column(data=np.zeros(num_comps), name="V_CURVE_CPL")
        main_cols.append(this_col)

    if lin_mod_type:
        this_col = Column(data=[use_lin_mod_types]*num_comps, name="LIN_MOD_TYPE")
        main_cols.append(this_col)

    if rm:
        this_col = Column(data=np.zeros(num_comps), name="RM")
        main_cols.append(this_col)

    if intr_pol_angle:
        this_col = Column(data=np.zeros(num_comps), name="INTR_POL_ANGLE")
        main_cols.append(this_col)

    if lin_pol_frac:
        this_col = Column(data=np.zeros(num_comps), name="LIN_POL_FRAC")
        main_cols.append(this_col)

    if lin_norm_comp_pl:
        this_col = Column(data=np.zeros(num_comps), name="LIN_NORM_COMP_PL")
        main_cols.append(this_col)

    if lin_alpha_pl:
        this_col = Column(data=np.zeros(num_comps), name="LIN_ALPHA_PL")
        main_cols.append(this_col)

    if lin_norm_comp_cpl:
        this_col = Column(data=np.zeros(num_comps), name="LIN_NORM_COMP_CPL")
        main_cols.append(this_col)

    if lin_alpha_cpl:
        this_col = Column(data=np.zeros(num_comps), name="LIN_ALPHA_CPL")
        main_cols.append(this_col)

    if lin_curve_cpl:
        this_col = Column(data=np.zeros(num_comps), name="LIN_CURVE_CPL")
        main_cols.append(this_col)
        
    if int_flux:
        this_col = Column(data=np.zeros(num_comps), name="INT_FLUX200.0")
        main_cols.append(this_col)
        
        
    main_table = Table()
    main_table.add_columns(main_cols)
    
    hdulist = [fits.PrimaryHDU(), fits.table_to_hdu(main_table)]
    table_names = ['MAIN']
    
    if shape_table:
        shape_cols = []
        if s_names:
            this_col = Column(data=np.zeros(num_comps), name="NAME")
            shape_cols.append(this_col)
        if s_n1s:
            this_col = Column(data=np.full(num_comps, shape_max_n), name="N1")
            shape_cols.append(this_col)
        if s_n2s:
            this_col = Column(data=np.full(num_comps, shape_max_n), name="N2")
            shape_cols.append(this_col)
        if s_coeffs:
            this_col = Column(data=np.zeros(num_comps), name="COEFF")
            shape_cols.append(this_col)
        
        shape_table = Table()
        shape_table.add_columns(shape_cols)
        table_names.append('SHAPELET')

        hdulist.append(fits.table_to_hdu(shape_table))
        
    if q_flux_table:
        q_flux_cols = []
        if q_flux_names:
            this_col = Column(data=[f'bleh_C{i:2d}' for i in range(num_comps)], name="NAME")
            q_flux_cols.append(this_col)
        if q_flux_fluxes:
            this_col = Column(data=np.zeros(num_comps), name="Q_INT_FLX200.0")
            q_flux_cols.append(this_col)
        
        q_flux_table = Table()
        q_flux_table.add_columns(q_flux_cols)
        table_names.append('Q_LIST_FLUXES')

        hdulist.append(fits.table_to_hdu(q_flux_table))
        
    if u_flux_table:
        u_flux_cols = []
        if u_flux_names:
            this_col = Column(data=[f'bleh_C{i:2d}' for i in range(num_comps)], name="NAME")
            u_flux_cols.append(this_col)
        if u_flux_fluxes:
            this_col = Column(data=np.zeros(num_comps), name="U_INT_FLX200.0")
            u_flux_cols.append(this_col)
        
        u_flux_table = Table()
        u_flux_table.add_columns(u_flux_cols)
        table_names.append('U_LIST_FLUXES')

        hdulist.append(fits.table_to_hdu(u_flux_table))
        
    if p_flux_table:
        p_flux_cols = []
        if p_flux_names:
            this_col = Column(data=[f'bleh_C{i:2d}' for i in range(num_comps)], name="NAME")
            p_flux_cols.append(this_col)
        if p_flux_fluxes:
            this_col = Column(data=np.zeros(num_comps), name="P_INT_FLX200.0")
            p_flux_cols.append(this_col)
        
        p_flux_table = Table()
        p_flux_table.add_columns(p_flux_cols)
        table_names.append('P_LIST_FLUXES')

        hdulist.append(fits.table_to_hdu(p_flux_table))
        
    if v_flux_table:
        v_flux_cols = []
        if v_flux_names:
            this_col = Column(data=[f'bleh_C{i:2d}' for i in range(num_comps)], name="NAME")
            v_flux_cols.append(this_col)
        if v_flux_fluxes:
            this_col = Column(data=np.zeros(num_comps), name="V_INT_FLX200.0")
            v_flux_cols.append(this_col)
        
        v_flux_table = Table()
        v_flux_table.add_columns(v_flux_cols)
        table_names.append('V_LIST_FLUXES')

        hdulist.append(fits.table_to_hdu(v_flux_table))
        
    hdulist = fits.HDUList(hdulist)
    
    if use_table_names:
        for index, name in enumerate(table_names):
            hdulist[index+1].name = name
    
    hdulist.writeto("test_full_skymodel.fits", overwrite=True)
    
    
##Vehicle for running tests
class Test(unittest.TestCase):
    """"""
    
    def assert_check_errors(self):
        """Assert that the check_columns_fits function raises an error"""
        with self.assertRaises(SystemExit) as cm:
            read_fits_skymodel.check_columns_fits('test_full_skymodel.fits')
        
        ##Have a read of the exception code
        print(cm.exception.code)
            
        # read_fits_skymodel.check_columns_fits('test_full_skymodel.fits')
            
    def make_sky_test_fails(self, **kwargs):
        """Given a set of kwargs, write a FITS file with the specified columns,
        and test that it fails the check_columns_fits test"""
        write_borked_skymodel(**kwargs)
        self.assert_check_errors()

    def test_missing_essential(self):
        """Ensure we fail if any of the essential columns are missing"""
        
        self.make_sky_test_fails(unq_source_id=False)
        self.make_sky_test_fails(name=False)
        self.make_sky_test_fails(ra=False)
        self.make_sky_test_fails(dec=False)
        self.make_sky_test_fails(comp_type=False)
        self.make_sky_test_fails(mod_type=False)
        
    def test_missing_gaussian(self):
        """Ensure we fail if we have Gaussians and are missing relevant columns"""
        
        self.make_sky_test_fails(use_comp_types='G', major_dc=False)
        self.make_sky_test_fails(use_comp_types='G', minor_dc=False)
        self.make_sky_test_fails(use_comp_types='G', pa_dc=False)
    
    def test_missing_shapelet(self):
        """Ensure we fail if we have Shapelet and are missing relevant columns
        in main table, missing Shapelet table, or missing Shapelet table columns"""
        
        ##Shapelet columns needed
        self.make_sky_test_fails(use_comp_types='S', major_dc=False)
        self.make_sky_test_fails(use_comp_types='S', minor_dc=False)
        self.make_sky_test_fails(use_comp_types='S', pa_dc=False)
        
        ##Checking the shapelet hdu exists, and contains what it needs
        self.make_sky_test_fails(use_comp_types='S', shape_table=False)
        self.make_sky_test_fails(use_comp_types='S', s_names=False)
        self.make_sky_test_fails(use_comp_types='S', s_n1s=False)
        self.make_sky_test_fails(use_comp_types='S', s_n2s=False)
        self.make_sky_test_fails(use_comp_types='S', s_coeffs=False)
        
        ##Error out if maximum basis order is too high
        self.make_sky_test_fails(use_comp_types='S', shape_max_n=101)
        
        ##Check things still fail the same way if we don't use table names
        self.make_sky_test_fails(use_comp_types='S', use_table_names=False, shape_max_n=101)
        
    def test_missing_stokesI(self):
        """Ensure we fail Stokes I flux models if missing columns"""
        
        ##Shapelet columns needed
        self.make_sky_test_fails(use_mod_types='pl', norm_comp_pl=False)
        self.make_sky_test_fails(use_mod_types='pl', alpha_pl=False)
        self.make_sky_test_fails(use_mod_types='cpl', norm_comp_cpl=False)
        self.make_sky_test_fails(use_mod_types='cpl', alpha_cpl=False)
        self.make_sky_test_fails(use_mod_types='cpl', curve_cpl=False)
        self.make_sky_test_fails(use_mod_types='nan', int_flux=False)
        
    def test_missing_stokesV(self):
        """Ensure we fail Stokes V flux models if missing columns"""
        
        ##Shapelet columns needed
        self.make_sky_test_fails(use_v_mod_types='pl', v_norm_comp_pl=False)
        self.make_sky_test_fails(use_v_mod_types='pl', v_alpha_pl=False)
        self.make_sky_test_fails(use_v_mod_types='cpl', v_norm_comp_cpl=False)
        self.make_sky_test_fails(use_v_mod_types='cpl', v_alpha_cpl=False)
        self.make_sky_test_fails(use_v_mod_types='cpl', v_curve_cpl=False)
        self.make_sky_test_fails(use_v_mod_types='pf', v_pol_frac=False)
        self.make_sky_test_fails(use_v_mod_types='nan', v_flux_table=False)
        self.make_sky_test_fails(use_v_mod_types='nan', v_flux_names=False)
        self.make_sky_test_fails(use_v_mod_types='nan', v_flux_fluxes=False)
        
    def test_missing_lin_pol(self):
        """Ensure we fail linear polarisation models if missing columns"""
        
        ##Shapelet columns needed
        self.make_sky_test_fails(use_lin_mod_types='pl', lin_norm_comp_pl=False)
        self.make_sky_test_fails(use_lin_mod_types='pl', lin_alpha_pl=False)
        self.make_sky_test_fails(use_lin_mod_types='pl', rm=False)
        self.make_sky_test_fails(use_lin_mod_types='cpl', lin_norm_comp_cpl=False)
        self.make_sky_test_fails(use_lin_mod_types='cpl', lin_alpha_cpl=False)
        self.make_sky_test_fails(use_lin_mod_types='cpl', lin_curve_cpl=False)
        self.make_sky_test_fails(use_lin_mod_types='cpl', rm=False)
        self.make_sky_test_fails(use_lin_mod_types='pf', lin_pol_frac=False)
        self.make_sky_test_fails(use_lin_mod_types='pf', rm=False)
        
        self.make_sky_test_fails(use_lin_mod_types='nan', q_flux_table=False)
        self.make_sky_test_fails(use_lin_mod_types='nan', q_flux_names=False)
        self.make_sky_test_fails(use_lin_mod_types='nan', q_flux_fluxes=False)
        self.make_sky_test_fails(use_lin_mod_types='nan', u_flux_table=False)
        self.make_sky_test_fails(use_lin_mod_types='nan', u_flux_names=False)
        self.make_sky_test_fails(use_lin_mod_types='nan', u_flux_fluxes=False)
        
        self.make_sky_test_fails(use_lin_mod_types='p_nan', rm=False)
        self.make_sky_test_fails(use_lin_mod_types='p_nan', p_flux_table=False)
        self.make_sky_test_fails(use_lin_mod_types='p_nan', p_flux_names=False)
        self.make_sky_test_fails(use_lin_mod_types='p_nan', p_flux_fluxes=False)
        
    def test_how_does_this_work(self):
        """"""
        write_borked_skymodel(use_table_names=False)
        
        with fits.open('test_full_skymodel.fits') as hdulist:
            names = [hdu.name for hdu in hdulist]
            # print(dir(hdulist))
            # print(hdulist.info())
            # print(hdulist[1].name)
            # print(hdulist[2].name)
            print(names)
        
        
    
        
##Run the test
if __name__ == '__main__':
    unittest.main()
    