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

``WODEN`` is C / CUDA code designed to be able to simulate low-frequency radio interferometric data. It is written to be simplistic and *fast* to allow all-sky simulations. Although ``WODEN`` was primarily written to simulate Murchinson Widefield Array (MWA, `Tingay et al. 2013`_) visibilities, it is becoming less instrument-specific as time goes on. `WODEN` outputs `uvfits` files.

The unique part of ``WODEN`` is that it can simulate shapelet model sources (along with point and Gaussian) that are compatible with the ``RTS`` (`Mitchell et al. 2008`_). These models are generated with SHApelet Modelling For Interferometers (`SHAMFI`_), specified with the ``--woden_srclist`` option. It also includes a script to convert a multi-scale CLEAN component list out of `WSClean`_ into a ``WODEN``-style srclist (when running ``WSClean`` use the ``-save-source-list`` option). ``WODEN`` can also produce visibilities that can be fed directly into the ``RTS`` to allow testing of calibration and modelling methodologies.

Documentation
-----------------

.. toctree::
   :maxdepth: 2

   installation/installation
   testing/cmake_testing
   testing/script_testing
   operating_principles/operating_principles
   examples/example_simulations
   API_reference/API_index

.. Indices and tables
.. --------------------
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
