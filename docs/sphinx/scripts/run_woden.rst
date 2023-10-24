``run_woden.py``
=================

This is the main ``WODEN`` executable. It takes command line arguments, creates and array layout, visibility containers, read sky models, and launches GPU code to calculate the visibilities. Finally, it writes the outputs to ``uvfits`` files.

.. _run_woden command line running options:

*Command line running options*
-------------------------------

.. argparse::
   :filename: ../../scripts/run_woden.py
   :func: get_parser
   :prog: run_woden.py

*Function documentation*
------------------------

.. automodule:: run_woden
   :members:
