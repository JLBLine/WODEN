``concat_woden_uvfits.py``
===========================

Helper script to concatenate a number of uvfits files by frequency into a
single uvfits file. Essential to input ``WODEN`` simulations into ``hyperdrive``.


.. warning:: This script is pretty dumb, so as long as the script can find the uvfits files, it will do the concatenation - it does NOT check if they have the same phase centre, frequencies, time stamps, etc. etc. so use at your own risk

.. _concat_woden_uvfits command line running options:

*Command line running options*
-------------------------------

.. argparse::
   :filename: ../../scripts/concat_woden_uvfits.py
   :func: get_parser
   :prog: concat_woden_uvfits.py

*Function documentation*
------------------------

.. automodule:: concat_woden_uvfits
   :members:
