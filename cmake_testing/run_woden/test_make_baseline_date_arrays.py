from sys import path
import os
import unittest
import numpy as np

##Do some disgusting path finding exercise, there must be a better way
##to do this
fileloc = os.path.realpath(__file__)
path.append('{:s}/../../src'.format(('/').join(fileloc.split('/')[:-1])))

##Code we are testing
import run_woden as rw

##Vehicle for running tests
class Test(unittest.TestCase):
    def test_make_baseline_date_arrays(self):
        """Tests the `rw.make_baseline_date_arrays` function, which should make
        the DATE and BASELINE arrays that are needed to populate a uvfits file"""

        ##Run the code with an example date
        date = "2019-06-12T13:04:12"

        num_antennas = 4
        num_time_steps = 3
        time_res = 2.0

        baselines_array, date_array = rw.make_baseline_date_arrays(num_antennas,
                                                 date, num_time_steps, time_res)

        expected_baselines = np.array([258, 259, 260, 515, 516, 772, 258,
                                        259, 260, 515, 516, 772, 258, 259,
                                        260, 515, 516, 772])

        expected_dates = np.array([0.04458333, 0.04458333, 0.04458333, 0.04458333,
                                   0.04458333, 0.04458333, 0.04460648, 0.04460648,
                                   0.04460648, 0.04460648, 0.04460648, 0.04460648,
                                   0.04462963, 0.04462963, 0.04462963, 0.04462963,
                                   0.04462963, 0.04462963])

        ##Check the outputs are within 1e-10 of expected values
        self.assertTrue(np.allclose(expected_baselines, baselines_array, atol=1e-10))
        self.assertTrue(np.allclose(expected_dates, date_array, atol=1e-10))

##Run the test
if __name__ == '__main__':
   unittest.main()