from sys import path
import os
import unittest

##Do some disgusting path finding exercise, there must be a better way
##to do this
fileloc = os.path.realpath(__file__)
path.append('{:s}/../../../wodenpy/array_layout'.format(('/').join(fileloc.split('/')[:-1])))

import create_array_layout

##Vehicle for running tests
class Test(unittest.TestCase):

    def make_inputs(self):
        """Makes a number of input and output pairs for the tests below"""
        self.b1_1 = 227
        self.b2_1 = 158

        self.b1_2 = 75
        self.b2_2 = 240

        self.b1_3 = 123
        self.b2_3 = 157

        self.b1_4 = 1
        self.b2_4 = 2

        self.b1_5 = 1
        self.b2_5 = 1

        self.code1 = 58270
        self.code2 = 19440
        self.code3 = 31645
        self.code4 = 258
        self.code5 = 257

    def test_RTS_encode_baseline(self):
        """Tests the `create_array_layout.RTS_encode_baseline` function, which should produce
        a single int to encode a pair of antennas based on the ancient ways
        of the AIPS. Just test with a number of pairs"""

        ##Make the inputs
        self.make_inputs()

        ##Test they match expectations
        self.assertEqual(self.code1, create_array_layout.RTS_encode_baseline(self.b1_1, self.b2_1))
        self.assertEqual(self.code2, create_array_layout.RTS_encode_baseline(self.b1_2, self.b2_2))
        self.assertEqual(self.code3, create_array_layout.RTS_encode_baseline(self.b1_3, self.b2_3))
        self.assertEqual(self.code4, create_array_layout.RTS_encode_baseline(self.b1_4, self.b2_4))
        self.assertEqual(self.code5, create_array_layout.RTS_encode_baseline(self.b1_5, self.b2_5))

    def test_RTS_decode_baseline(self):
        """Tests the `create_array_layout.RTS_decode_baseline` function, which should split a
        single baseline code into to antenna numbers, based on the ancient ways
        of the AIPS. Just test with a number of codes"""

        ##Make the inputs
        self.make_inputs()

        ##call the function under test
        b1_1, b2_1 = create_array_layout.RTS_decode_baseline(self.code1)
        b1_2, b2_2 = create_array_layout.RTS_decode_baseline(self.code2)
        b1_3, b2_3 = create_array_layout.RTS_decode_baseline(self.code3)
        b1_4, b2_4 = create_array_layout.RTS_decode_baseline(self.code4)
        b1_5, b2_5 = create_array_layout.RTS_decode_baseline(self.code5)

        ##Test they match expectations
        self.assertEqual(b1_1, self.b1_1)
        self.assertEqual(b1_2, self.b1_2)
        self.assertEqual(b1_3, self.b1_3)
        self.assertEqual(b1_4, self.b1_4)
        self.assertEqual(b2_1, self.b2_1)
        self.assertEqual(b2_2, self.b2_2)
        self.assertEqual(b2_3, self.b2_3)
        self.assertEqual(b2_4, self.b2_4)

        self.assertEqual(b1_5, self.b1_5)
        self.assertEqual(b2_5, self.b2_5)

##Run the test
if __name__ == '__main__':
   unittest.main()
