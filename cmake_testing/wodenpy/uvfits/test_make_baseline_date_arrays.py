from sys import path
import os
import unittest
import numpy as np

from wodenpy.uvfits import wodenpy_uvfits

##Vehicle for running tests
class Test(unittest.TestCase):
    def test_make_baseline_date_arrays_noautos(self):
        """Tests the `wodenpy_uvfits.make_baseline_date_arrays` function, which should make
        the DATE and BASELINE arrays that are needed to populate a uvfits file"""

        ##Run the code with an example date
        date = "2019-06-12T13:04:12"

        num_antennas = 4
        num_time_steps = 3
        time_res = 2.0

        baselines_array, date_array = wodenpy_uvfits.make_baseline_date_arrays(num_antennas,
                                                 date, num_time_steps, time_res)

        expected_baselines = np.array([258, 259, 260, 515, 516, 772, 258,
                                        259, 260, 515, 516, 772, 258, 259,
                                        260, 515, 516, 772])

        expected_dates = np.array([0.04459490726, 0.04459490726, 0.04459490726,
                                   0.04459490726, 0.04459490726, 0.04459490726,
                                   0.04461805541, 0.04461805541, 0.04461805541,
                                   0.04461805541, 0.04461805541, 0.04461805541,
                                   0.04464120356, 0.04464120356, 0.04464120356,
                                   0.04464120356, 0.04464120356, 0.04464120356])

        ##Check the outputs are within 1e-10 of expected values
        self.assertTrue(np.allclose(expected_baselines, baselines_array, atol=1e-10))
        self.assertTrue(np.allclose(expected_dates, date_array, atol=1e-10))

    def test_make_baseline_date_arrays_withautos(self):
        """Tests the `wodenpy_uvfits.make_baseline_date_arrays` function, which should make
        the DATE and BASELINE arrays that are needed to populate a uvfits file"""

        ##Run the code with an example date
        date = "2019-06-12T13:04:12"

        num_antennas = 4
        num_time_steps = 3
        time_res = 2.0

        baselines_array, date_array = wodenpy_uvfits.make_baseline_date_arrays(num_antennas,
                                                 date, num_time_steps, time_res,
                                                 do_autos=True)

        expected_baselines = np.array([257, 258, 259, 260, 514,
                                       515, 516, 771, 772, 1028,
                                       257, 258, 259, 260, 514,
                                       515, 516, 771, 772, 1028,
                                       257, 258, 259, 260, 514,
                                       515, 516, 771, 772, 1028])

        expected_dates = np.array([0.04459490726, 0.04459490726, 0.04459490726,
                                   0.04459490726, 0.04459490726, 0.04459490726,
                                   0.04459490726, 0.04459490726, 0.04459490726,
                                   0.04459490726,
                                   0.04461805541, 0.04461805541, 0.04461805541,
                                   0.04461805541, 0.04461805541, 0.04461805541,
                                   0.04461805541, 0.04461805541, 0.04461805541,
                                   0.04461805541,
                                   0.04464120356, 0.04464120356, 0.04464120356,
                                   0.04464120356, 0.04464120356, 0.04464120356,
                                   0.04464120356, 0.04464120356, 0.04464120356,
                                   0.04464120356])

        ##Check the outputs are within 1e-10 of expected values
        self.assertTrue(np.allclose(expected_baselines, baselines_array, atol=1e-10))
        self.assertTrue(np.allclose(expected_dates, date_array, atol=1e-10))

##Run the test
if __name__ == "__main__":
   unittest.main()
