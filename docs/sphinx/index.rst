.. WODEN documentation master file, created by
   sphinx-quickstart on Mon Mar 22 14:37:03 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. _WSClean: https://sourceforge.net/projects/wsclean/
.. _Mitchell et al. 2008: https://ieeexplore.ieee.org/document/4703504?arnumber=4703504
.. _SHAMFI: https://github.com/JLBLine/SHAMFI
.. _Tingay et al. 2013: https://doi.org/10.1017/pasa.2012.007

The WODEN visibility simulator
=================================

.. note::
  Before ``WODEN`` version 1.4.0, in the output `uvfits` files, the first polarisation (usually called XX) was derived from North-South dipoles, as is the labelling convention according to the IAU. However, most `uvfits` users I've met, as well as the data out of the MWA telescope, define XX as East-West. So although the internal labelling and mathematics within the C/CUDA code is to IAU spec, by default, ``run_woden.py`` now writes out XX as East-West and YY as North-South. From version 1.4.0, a header value of ``IAUORDER = F`` will appear, with ``F`` meaning IAU ordering is False, so the polarisations go EW-EW, NS-NS, EW-NS, NS-EW. If ``IAUORDER = T``, the order is NS-NS, EW-EW, NS-EW, EW-NS. If there is no ``IAUORDER`` at all, assume ``IAUORDER = T``.

``WODEN`` is C / CUDA code designed to be able to simulate low-frequency radio interferometric data. It is written to be simplistic and *fast* to allow all-sky simulations. Although ``WODEN`` was primarily written to simulate Murchinson Widefield Array (MWA, `Tingay et al. 2013`_) visibilities, it is becoming less instrument-specific as time goes on. `WODEN` outputs `uvfits` files.

``WODEN`` has been written with Stokes polarisations in mind. Currently, only Stokes I information is read in, which is then propagated fully through the polarised instrumental response (depending on which primary beam you select), and output into Stokes `XX,XY,YX,YY` polarisations. See :ref:`sky model formats` for more information.

Documentation
-----------------

.. toctree::
   :maxdepth: 2

   installation/installation
   testing/testing
   scripts/scripts
   operating_principles/operating_principles
   examples/example_simulations
   API_reference/API_index
   code_graphs/code_graphs

.. Indices and tables
.. --------------------
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
