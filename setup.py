import setuptools
from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
import os
from wodenpy.wodenpy_setup.git_helper import make_gitdict
import numpy as np

class GitInfo(setuptools.Command):
  '''A custom command to create a json file containing wodenpy git information.'''

  description = 'Create the file "wodenpy/wodenpy_gitinfo.json" containing git information '
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
        print('Creating file wodenpy/wodenpy_gitinfo.npz')

  def run(self):
    '''Write the wodenpy git npz file.'''

    ##Find where we are running the pip install from, and add in a sensible
    ##place to save the git dictionary
    save_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),'wodenpy', 'wodenpy_gitinfo.npz')

    git_dict = make_gitdict()
    np.savez(save_path, **git_dict)


class BuildPyCommand(_build_py):
  '''Custom build command to run the gitinfo command during build'''

  def run(self):
    self.run_command('gitinfo')
    _build_py.run(self)

setup(
    name = "wodenpy",
    version = '2.2.0',
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
                               'wodenpy_gitinfo.npz',
                               'bandpass_1kHz.txt']},
    cmdclass={'gitinfo': GitInfo,
              'build_py': BuildPyCommand,
              },
    install_requires=['numpy<=1.26.0'
                     ]
)
