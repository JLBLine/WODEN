# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import subprocess

sys.path.insert(0, os.path.abspath('../../src/'))
sys.path.insert(0, os.path.abspath('../../wodenpy/'))
sys.path.insert(0, os.path.abspath('../../wodenpy/array_layout'))
sys.path.insert(0, os.path.abspath('../../wodenpy/observational'))
sys.path.insert(0, os.path.abspath('../../wodenpy/phase_rotate'))
sys.path.insert(0, os.path.abspath('../../wodenpy/skymodel'))
sys.path.insert(0, os.path.abspath('../../wodenpy/use_libwoden'))
sys.path.insert(0, os.path.abspath('../../wodenpy/uvfits'))
sys.path.insert(0, os.path.abspath('../../wodenpy/wodenpy_setup'))
sys.path.insert(0, os.path.abspath('../../scripts/'))

read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

if read_the_docs_build:
    sys.path.insert(0, os.path.abspath('../../'))

pwd = os.getcwd()
os.chdir('../doxygen')
print("pwd is {:s}".format(os.getcwd()))

##When updating documentation locally, can be really annoying running
##doxygen when you don't need it, so set
##export NO_RUN_WODEN_DOXYGEN=1
##to switch it off. Can set
##export NO_RUN_WODEN_DOXYGEN=0
##to turn doxygen back on
try:
    norun = os.environ['NO_RUN_WODEN_DOXYGEN']
    if int(norun):
        run_doxygen = False
    else:
        run_doxygen = True
except KeyError:
    run_doxygen = True

if run_doxygen: subprocess.call('doxygen Doxyfile', shell=True)

# -- Project information -----------------------------------------------------

project = 'WODEN'
copyright = '2024, J.L.B. Line'
author = 'J.L.B. Line'

# The full version, including alpha/beta/rc tags
release = '2.3'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinxarg.ext',
    'sphinx.ext.mathjax',
    'sphinx_math_dollar',
    'sphinx.ext.napoleon',
    'sphinx.ext.coverage',
    'sphinx.ext.autosectionlabel',
    'breathe',
    "myst_nb",
]

nb_execution_mode = "off"

todo_include_todos = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'classic'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']


breathe_projects = { "WODEN": "../doxygen/xml/" }
breathe_default_project = "WODEN"

##Add these in so the c++ acknowledges CUDA specific attributes
cpp_id_attributes = ["__global__", "__device__", "_Complex"]
