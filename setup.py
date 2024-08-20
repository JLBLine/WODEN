import setuptools
from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
import os
from subprocess import check_output

def get_commandline_output(command_list):
    """
    Takes a command line entry separated into list entries, and returns the
    output from the command line as a string

    Parameters
    ----------
    command_list : list of strings
        list of strings that when combined form a coherent command to input into
        the command line

    Returns
    -------
    output : string
        the output result of running the command

    """
    output = check_output(command_list,universal_newlines=True).strip()
    return output

def make_gitdict():
    """
    Makes a dictionary containing key git information about the repo by running
    specific commands on the command line

    Returns
    -------
    git_dict : dictionary
        A dictionary containing git information with keywords: describe, date,
        branch

    """

    ##Try and get a git version. If this is a release version, it might not
    ##be in a git repo so try, otherwise return False
    try:
        git_dict = {
            'describe': get_commandline_output(["git", "describe", "--always"]),
            'date': get_commandline_output(["git", "log", "-1", "--format=%cd"]),
            'branch': get_commandline_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        }
        
    except:
        git_dict = False

    return git_dict

class GitInfo(setuptools.Command):
  '''A custom command to create a json file containing wodenpy git information.'''

  description = 'Create the file "wodenpy/wodenpy_gitinfo.txt" containing git information '
  user_options = []

  def initialize_options(self):
    '''Set default values for options (this has to be included for
    setuptools.Command to work)'''
    # Each user option must be listed here with their default value.
    self.git_info = True

  def finalize_options(self):
    '''Post-process options (this has to be included for
    setuptools.Command to work)'''
    if self.git_info:
        print('Creating file wodenpy/wodenpy_gitinfo.txt')

  def run(self):
    '''Write the wodenpy git npz file.'''

    ##Find where we are running the pip install from, and add in a sensible
    ##place to save the git dictionary
    save_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),'wodenpy', 'wodenpy_gitinfo.txt')

    git_dict = make_gitdict()
    with open(save_path, "w") as fp:
        for key in git_dict.keys():
            fp.write("{:s},{:s}\n".format(key, git_dict[key]))

class BuildPyCommand(_build_py):
  '''Custom build command to run the gitinfo command during build'''

  def run(self):
    self.run_command('gitinfo')
    _build_py.run(self)

setup(
    name = "wodenpy",
    version = '2.3.0',
    author = "Jack L. B. Line",
    url = "https://github.com/JLBLine/WODEN",
    python_requires=">=3.7",
    description = 'GPU-accelerated code to simulate radio-frequency interferometric observations',
    long_description = open("README.md").read(),
    long_description_content_type = 'text/markdown',
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: Mozilla Public License Version 2.0",
        "Operating System :: Linux",
    ],
    packages = ['wodenpy',
                'wodenpy.wodenpy_setup',
                'wodenpy.use_libwoden',
                'wodenpy.observational',
                'wodenpy.uvfits',
                'wodenpy.array_layout',
                'wodenpy.phase_rotate',
                'wodenpy.skymodel'],
    scripts=["scripts/run_woden.py",
             "scripts/add_woden_uvfits.py",
             "scripts/concat_woden_uvfits.py",
             "scripts/woden_uv2ms.py",
             "scripts/add_instrumental_effects_woden.py"],
    package_data={"wodenpy" : ["libwoden_float.so",
                               "libwoden_double.so",
                               'wodenpy_gitinfo.txt',
                               'bandpass_1kHz.txt']},
    cmdclass={'gitinfo': GitInfo,
              'build_py': BuildPyCommand,
              },
    requirements=['importlib_resources']
)
