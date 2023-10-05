``add_woden_uvfits.py``
========================

Helper script to add visibilities inside of two sets of ``uvfits`` files.
Useful for example if you have a set of 21cm signal visibilities that you want
to add different types of foregrounds to; you only have to run the 21cm
simulation once, then run different foregrounds with the same observational
settings, and just add the foregrounds on top as visibilities are additive.

.. warning:: This script is pretty dumb, so as long as the data have the same shape it will add the contents - it does NOT check if they have the same phase centre, frequencies, time stamps, etc. etc. use at your own risk

.. _add_woden_uvfits command line running options:

*Command line running options*
-------------------------------

.. argparse::
   :filename: ../../scripts/add_woden_uvfits.py
   :func: get_parser
   :prog: add_woden_uvfits.py

*Function documentation*
------------------------

.. automodule:: add_woden_uvfits
   :members:
