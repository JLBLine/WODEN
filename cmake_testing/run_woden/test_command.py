from sys import path
import os
import unittest

##Do some disgusting path finding exercise, there must be a better way
##to do this
fileloc = os.path.realpath(__file__)
path.append('{:s}/../../src'.format(('/').join(fileloc.split('/')[:-1])))

##Code we are testing
import run_woden as rw

##Vehicle for running tests
class Test(unittest.TestCase):
    def test_command(self):
        """Tests the `rw.command` function, which should call
        things on the command line. Test by trying to echo the word cheese
        into a text file, and reading the output"""

        ##Run code to be tested, should run this on command line
        cmd = rw.command("echo cheese > example.txt")

        ##Check it's made a file called "example.txt"
        self.assertTrue(os.path.isfile("example.txt"))

        ##Check the contents of example.txt make sense
        with open("example.txt",'r') as infile:
            readwords = infile.read().strip('\n')
            self.assertEqual("cheese", readwords)

##Run test
if __name__ == '__main__':
   unittest.main()
