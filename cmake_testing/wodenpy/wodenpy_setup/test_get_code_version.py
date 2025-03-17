from sys import path
import os
import unittest
import numpy as np

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

# ##Code we are testing
from wodenpy.wodenpy_setup import run_setup
from wodenpy.wodenpy_setup import git_helper



##Vehicle for running tests
class Test(unittest.TestCase):
    """Test whether the args collected by the argument parser are read in
    correctly, and that sanity checks on certain combinations work such that
    we don't feed WODEN arguments that won't work"""

    def test_get_code_version(self):
        """Call `run_setup.code_version` and just check it runs without error.
        Hard to automate a test for this one"""
        version = run_setup.get_code_version()
        print(version)
        
    def test_make_gitdict(self):
        """Call `git_helper.make_gitdict` and just check it runs without error.
        Then try and change current working dir so it returns something False"""
        git_dict = git_helper.make_gitdict()
        
        ##check we can access stuff from the dict
        print(git_dict['describe'])
        print(git_dict['date'])
        print(git_dict['branch'])
        
        cwd = os.getcwd()
        ##now change to dir to something that isn't a git repo
        os.chdir('/')
        git_dict = git_helper.make_gitdict()
        self.assertFalse(git_dict)
        
        ##Get back to cwd
        os.chdir(cwd)
        
    def test_retrieve_gitdict(self):
        """Call `git_helper.retrieve_gitdict`. Needs to have a proper
        installation of WODEN to cover the most lines"""
        git_dict = git_helper.retrieve_gitdict()
        print(git_dict)

##Run the test
if __name__ == '__main__':
    unittest.main()
    