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
.. _EveryBeam: https://everybeam.readthedocs.io/en/latest/index.html
.. _Line et al. 2022: https://joss.theoj.org/papers/10.21105/joss.03676
.. _hyperbeam: https://github.com/MWATelescope/mwa_hyperbeam
.. _pyuvdata UVBeam: https://pyuvdata.readthedocs.io/en/latest/uvbeam.html#uvbeam

The WODEN visibility simulator
=================================

``WODEN`` is Python / C / C++ / GPU code designed to be able to simulate low-frequency radio interferometric data. It is written to be simplistic and *fast* to allow all-sky simulations. Although ``WODEN`` was primarily written to simulate Murchinson Widefield Array (MWA, `Tingay et al. 2013`_) visibilities, it is becoming less instrument-specific as time goes on. `WODEN` outputs `uvfits` files. Currently, the ``MWA`` primary beam is supported via `hyperbeam`_; ``LOFAR`` and ``OSKAR``-style primary beams are supported vis the `EveryBeam`_ package; ``HERA``is supported via `pyuvdata UVBeam`_;  ``EDA2``, an analytic ``MWA`` beam, and a simple Gaussian beam are supported natively to ``WODEN``.

``WODEN`` has been written with Stokes polarisations in mind. A fully Stokes IQUV model is propagated through the polarised instrumental response (depending on which primary beam you select), and output into Stokes `XX,XY,YX,YY` polarisations. See :ref:`sky model formats` and :ref:`visibility calculations` for more information.

If you use ``WODEN`` in your research, please cite the JOSS paper `Line et al. 2022`_.

.. Documentation
.. -----------------

.. toctree::
   :maxdepth: 2

   installation/installation
   examples/example_simulations

.. toctree::
   :maxdepth: 3

   operating_principles/operating_principles
   testing/testing
   
.. toctree::
   :maxdepth: 2

   scripts/scripts
   API_reference/API_index
   code_graphs/code_graphs

.. toctree::
   :maxdepth: 4

   developers_guide/developers_guide

.. Indices and tables
.. --------------------
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
