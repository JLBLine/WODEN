.. _everybeam insallation page: https://everybeam.readthedocs.io/en/latest/build-instructions.html

.. _`EveryBeam testing`:

EveryBeam Testing and Development
======================================

.. toctree::
   :maxdepth: 2
   :caption: Click these links to go to notebooks of integration tests

   test_parallactic_rotation.ipynb
   test_MWA.ipynb
   test_LOFAR_LBA.ipynb
   test_LOFAR_HBA.ipynb
   test_OSKAR_MWA.ipynb
   test_OSKAR_SKA.ipynb
   


Installing EveryBeam
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The build instructions for EveryBeam live on the `everybeam insallation page`_. However, you need to clone my branch as currently it's the only branch where you can easily feed ``az,za`` coords directly into the beam. We want to do that in ``WODEN`` as we precess the array position back to J2000, instead of precessing the ``ra,dec`` coordinates forward to the current date. ``EveryBeam`` will automatically precess any ``ra,dec``. Clone my branch via::

   $ git clone -b mwa_python_wrapper --recursive -j4 https://git.astron.nl/RD/EveryBeam.git

When I ran ``cmake``, I did::

   cmake .. -DCMAKE_INSTALL_PREFIX=/home/jline/software/installed/ 

so I would know where the installation went.


Work still to be done on EveryBeam in WODEN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Only certain frequencies exists in certain models. This means we don't need to call EveryBeam for every single frequency. So create some kind of cache or map for the existing beam models to save compute.
- Getting rid of the reliance on my branch of EveryBeam. This will require a change in the EveryBeam API to allow for feeding in ``az,za`` directly, and working with the authors to get it merged.

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
