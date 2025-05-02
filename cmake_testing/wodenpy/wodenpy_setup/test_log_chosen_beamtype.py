from sys import path
import os
import unittest
import numpy as np
import numpy.testing as npt
import logging
from wodenpy.wodenpy_setup import run_setup, woden_logger
from wodenpy.use_libwoden.woden_settings import fill_woden_settings_python
from pathlib import Path

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

class Test(unittest.TestCase):
    """Test when we invoke WODEN with different beam settings, the correct
    information is dumped to the log file"""
    
    @classmethod
    def tearDownClass(cls):
        os.remove('test_log.log')

    def make_required_args(self, extra_args=[]):
        """These are arguments that if missing, the parser itself should fail"""
        self.inputs = ['--ra0=0.0', '--dec0=-26.7', '--cat_filename=srclist.txt',
                       '--num_time_steps=16', '--time_res=8.0', '--freq_res=40e+3',
                       '--array_layout={:s}/example_array_layout.txt'.format(code_dir),
                       '--lowest_channel_freq=160e+6', '--date=2019-06-12T13:04:12',
                       "--MWA_FEE_delays=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"]
        
        self.inputs += extra_args
        
        parser = run_setup.get_parser()
        args =  parser.parse_args(self.inputs)
        # args = run_setup.check_args(args)
        
        ##Don't actually run check_args as we want to test things that would
        ##fail it. But set things here that would be set by check_args
        args.num_antennas = 1
        args.dipamps = np.array([1.0])
        args.metafits_filename = "path/to/metafits"
        args.beam_ms_path = "path/to/beam_ms"
        args.coarse_band_width = 1.28e+6
        
        
        if args.primary_beam == 'everybeam_LOFAR' or args.primary_beam == 'everybeam_OSKAR':
            args.eb_ra_point = args.ra0
            args.eb_dec_point = args.dec0
            args.pointed_ms_file_name = Path("path/to/pointed_ms")
            
        if args.primary_beam == 'uvbeam_HERA':
            args.cst_file_list = 'uwotm8.txt'
            args.cst_freqs = [100e+6, 150e+6, 200e+6]
        
        self.args = args
        
    def make_test_logger(self):
        self.logger = logging.getLogger()
        logging.basicConfig(level=logging.DEBUG,
                            format='%(message)s',
                            # datefmt='%Y-%m-%d %H:%M:%S',
                            filename='test_log.log',
                            filemode='w')
        
    def make_woden_settings(self):
        jd_date, lst_deg = 0.0, 0.0
        self.args.band_nums = [1]
        self.woden_settings = fill_woden_settings_python(self.args, jd_date, lst_deg)
        
        ##Fill in some stuff that would be sorted by check_args; we don't want
        ##to run too many other functions here
        self.woden_settings.hdf5_beam_path="path/to/hdf5"
        self.woden_settings.gauss_dec_point = np.radians(self.args.dec0)
        
    def run_log_beamtype(self, extra_args=[]):
        self.make_required_args(extra_args)
        self.make_test_logger()
        self.make_woden_settings()
        woden_logger.log_chosen_beamtype(self.logger, self.woden_settings,
                                         self.args)
        
    def read_lines_check_against_expected(self, expected_lines):
        with open("test_log.log", 'r') as f:
            lines = f.readlines()
        # print(lines[-len(expected_lines):])
        # print(expected_lines)
        self.assertEqual(lines[-len(expected_lines):], expected_lines)
        
    def test_with_no_beam(self):
        self.run_log_beamtype()
        expected_lines = ['No primary beam selected, no beam attenuation will be applied.\n']
        self.read_lines_check_against_expected(expected_lines)
        
    def test_with_gaussian_beam(self):
        self.run_log_beamtype(['--primary_beam=Gaussian'])
        expected_lines = ['Using Gaussian primary beam with the following parameters:\n',
                         '\tLocked to pointing: HA 0.0 deg, Dec -26.7 deg\n',
                         '\tFWHM: 20.0 at reference frequency: 150.0 MHz\n']
        self.read_lines_check_against_expected(expected_lines)
        
        
    def test_with_mwafee_beam(self):
        for primary_beam in ['MWA_FEE', 'MWA_FEE_interp']:
            args = [f'--primary_beam={primary_beam}']
            self.run_log_beamtype(args)
            expected_lines = ['Using MWA primary beam via hyperbeam with the following parameters:\n',
                            '\thdf5 file: path/to/hdf5\n',
                            '\tdelays: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]\n',
                            '\tsetting all dipole amplitudes to 1.0\n',
                            '\twill not flag any dipoles\n']
            self.read_lines_check_against_expected(expected_lines)
            
            args.append('--use_MWA_dipamps')
            self.run_log_beamtype(args)
            expected_lines = ['Using MWA primary beam via hyperbeam with the following parameters:\n',
                            '\thdf5 file: path/to/hdf5\n',
                            '\tdelays: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]\n',
                            '\twill use dipole amplitudes from given metafits file\n',
                            '\twill not flag any dipoles\n',
                            '\tmetafits file: path/to/metafits\n']
            self.read_lines_check_against_expected(expected_lines)
            
            args.append('--use_MWA_dipflags')
            self.run_log_beamtype(args)
            expected_lines = ['Using MWA primary beam via hyperbeam with the following parameters:\n',
                            '\thdf5 file: path/to/hdf5\n',
                            '\tdelays: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]\n',
                            '\twill use dipole amplitudes from given metafits file\n',
                            '\twill use dipole flags from given metafits file\n',
                            '\tmetafits file: path/to/metafits\n']
            self.read_lines_check_against_expected(expected_lines)
        
    def test_with_analy_dipole_beam(self):
        self.run_log_beamtype(['--primary_beam=EDA2'])
        expected_lines = ['Using an analytical dipole primary beam (a.k.a each element is an MWA dipole e.g. EDA2 array).\n']
        self.read_lines_check_against_expected(expected_lines)
        
    def test_with_MWA_analy_beam(self):
        self.run_log_beamtype(['--primary_beam=MWA_analy'])
        expected_lines = ['Using MWA analytic primary beam with:\n',
                          '\tdelays: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]\n']
        self.read_lines_check_against_expected(expected_lines)
        
    def test_with_everybeam_OSKAR(self):
        self.run_log_beamtype(['--primary_beam=everybeam_OSKAR'])
        expected_lines = ['Will run with EveryBeam OSKAR primary beam, based on this measurement set:\n',
                          '\tpath/to/beam_ms\n',
                          f"Created the following minimal MS to point the beam:\n",
                          f"\tpath/to/pointed_ms\n",
                          f"Primary beam is pointed at RA,Dec = {self.args.eb_ra_point}, {self.args.eb_dec_point} deg\n"]
        self.read_lines_check_against_expected(expected_lines)
        
    def test_with_everybeam_LOFAR(self):
        self.run_log_beamtype(['--primary_beam=everybeam_LOFAR'])
        expected_lines = ['Will run with EveryBeam LOFAR primary beam, based on this measurement set:\n',
                          '\tpath/to/beam_ms\n',
                          f"Created the following minimal MS to point the beam:\n",
                          f"\tpath/to/pointed_ms\n",
                          f"Primary beam is pointed at RA,Dec = {self.args.eb_ra_point}, {self.args.eb_dec_point} deg\n"]
        self.read_lines_check_against_expected(expected_lines)
        
    def test_with_everybeam_MWA(self):
        self.run_log_beamtype(['--primary_beam=everybeam_MWA'])
        expected_lines = ['Will run with EveryBeam MWA primary beam, based on this measurement set:\n',
                          '\tpath/to/beam_ms\n',
                          f"Using the following hdf5 file:\n",
                          f"\tpath/to/hdf5\n"]
        self.read_lines_check_against_expected(expected_lines)
        
    def test_with_a_boo_boo(self):
        self.make_required_args()
        self.make_test_logger()
        self.make_woden_settings()
        ##THis should never happen in check_args has been run
        self.woden_settings.beamtype = -1
        
        woden_logger.log_chosen_beamtype(self.logger, self.woden_settings,
                                         self.args)
        
        expected_lines = ["Primary beam type not recognised. This shouldn't be possible if you've used wodenpy.woden_setup.run_setup.check_args().\n"]
        self.read_lines_check_against_expected(expected_lines)
        
    def test_with_uvbeam_HERA(self):
        self.run_log_beamtype(['--primary_beam=uvbeam_HERA'])
        expected_lines = [f"Will create a HERA beam from CST files listed in this file:\n",
                          f"\t{self.args.cst_file_list}\n",
                          f"Have read the following frequencies from this file:\n",
                          f"\t{self.args.cst_freqs}\n"]
        self.read_lines_check_against_expected(expected_lines)
    
if __name__ == '__main__':
    unittest.main()