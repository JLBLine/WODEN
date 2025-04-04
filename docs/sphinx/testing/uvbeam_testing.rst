.. _everybeam insallation page: https://everybeam.readthedocs.io/en/latest/build-instructions.html

.. _`pyuvdata UVBeam testing`:

EveryBeam Testing and Development
======================================

.. toctree::
   :maxdepth: 2
   :caption: Click these links to go to notebooks of integration tests

   test_uvbeam_MWA.ipynb
   


Work still to be done on UVBeam in WODEN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- tbd

.. Adding new beam models
.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. You need to update the following Python functions to add a new EveryBeam model (in general, look for anywhere ``EB_OSKAR`` is already used, and make sure whatever is done with ``EB_OSKAR`` is also done with your new beam model):

.. - Update the help for ``--primary_beam`` in the parser defined in ``WODEN/wodenpy/wodenpy_setup/run_setup.get_parser`` to show a new recognised beam model. A.k.a add a "my_new_beam" as an accaptable argument to ``--primary_beam``.
.. - Add this "my_new_beam" to the ``eb_args`` list in ``WODEN/wodenpy/wodenpy_setup/run_setup.check_args``.
.. - Add a new beamtype to the ``BeamTypes`` and ``BeamGroups`` enums in ``WODEN/wodenpy/use_libwoden/beam_settings.py``. You must also update the equivalent ``e_beamtype`` in ``WODEN/include/woden_struct_defs.h``. E.g. MWA already has a value of ``EB_MWA``, so try adding in a new ``EB_YOURBEAM``. If the beam model you add in has off-cardinal dipoles, you should add in to the ``BeamGroups.off_cardinal_beam_values`` list.
.. - Make sure your beam model is added to ``woden_settings`` in the ``WODEN/wodenpy/use_libwoden/woden_settings.create_woden_settings`` function.
.. - Add a new function to load in the beam model in ``WODEN/wodenpy/primary_beam/use_everybeam.py`` (similar to ``load_OSKAR_telescope`` or ``load_LOFAR_telescope``).
.. - In ``WODEN/wodenpy/skymodel/read_fits_skymodel.read_fits_skymodel_chunks``, ensure the new beam model initiated correctly.


.. You need to update the following GPU functions:

.. - Ensure your beam model is added to ``beam_settings`` in ``WODEN/src/primary_beam.c``
.. - At the start of ``calculate_visibilities`` in ``WODEN/src/calculate_visibilities.cpp``, ensure ``use_twobeams, num_beams`` are set correctly for new beam.
.. - Wherever ``copy_CPU_beam_gains_to_GPU`` is called in ``WODEN/src/calculate_visibilities.cpp``, ensure that happens for your new beam model.
.. - Add your model to ``WODEN/src/source_components.cpp.free_beam_gains`` function to ensure the leakage terms are freed correctly.
