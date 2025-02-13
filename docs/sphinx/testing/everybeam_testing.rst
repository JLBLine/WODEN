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
The build instructions for EveryBeam live on the `everybeam insallation page`_. However, you need to clone my branch as currently it's the only branch with an MWA Python wrapper::

   $ git clone -b mwa_python_wrapper --recursive -j4 https://git.astron.nl/RD/EveryBeam.git

When I ran ``cmake``, I did::

   cmake .. -DCMAKE_INSTALL_PREFIX=/home/jline/software/installed/ \
     -DBUILD_WITH_PYTHON=ON


so I would know where the installation went. I also did all of this inside a conda environment called ``woden_dev``, but I think ``cmake`` still found the system Python. Once installed, I had to do::
    
   export PYTHONPATH=$PYTHONPATH:"/home/jline/software/installed/lib/python3.12/site-packages"
   ln -s /home/jline/software/installed/share/everybeam /home/jline/software/anaconda3/envs/woden_dev/share/everybeam

which let me conda environment see everything it needed to. 

When I was running notebooks, which don't load stuff from system, only the conda environment, I had to do::

   conda install -c conda-forge libstdcxx-ng
   conda install hdf5

for certain things to work.


Work still to be done on EveryBeam in WODEN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- I've managed to read the antenna locations from a measurement set, needed to calculate ``u,v,w`` coords. But I haven't worked out where the central latitude/longitude are stored in the measurement set. At the moment you have to put that in via command line arguments. Something to figure out.
- You can only run one ``WODEN`` frequency band at a time. This is because ``WODEN`` is designed to calculate the beam values on GPU, iteratively band by band. It's not the end of the world but a slight design flaw.
- Only certain frequencies exists in certain models. This means we don't need to call EveryBeam for every single frequency. So create some kind of cache or map for the existing beam models to save compute.
- Things are threaded over ``EveryBeam now``, but we still need more optimisation. CPU-only version?
- Unit tests. Mostly everything is an integration test at the moment. Best practice is to break down the testing into units.

Adding new beam models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You need to update the following Python functions to add a new EveryBeam model (in general, look for anywhere ``EB_OSKAR`` is already used, and make sure whatever is done with ``EB_OSKAR`` is also done with your new beam model):

- Update the help for ``--primary_beam`` in the parser defined in ``WODEN/wodenpy/wodenpy_setup/run_setup.get_parser`` to show a new recognised beam model. A.k.a add a "my_new_beam" as an accaptable argument to ``--primary_beam``.
- Add this "my_new_beam" to the ``eb_args`` list in ``WODEN/wodenpy/wodenpy_setup/run_setup.check_args``.
- Add a new beamtype to the ``BeamTypes`` and ``BeamGroups`` enums in ``WODEN/wodenpy/use_libwoden/beam_settings.py``. You must also update the equivalent ``e_beamtype`` in ``WODEN/include/woden_struct_defs.h``. E.g. MWA already has a value of ``EB_MWA``, so try adding in a new ``EB_YOURBEAM``. If the beam model you add in has off-cardinal dipoles, you should add in to the ``BeamGroups.off_cardinal_beam_values`` list.
- Make sure your beam model is added to ``woden_settings`` in the ``WODEN/wodenpy/use_libwoden/woden_settings.create_woden_settings`` function.
- Add a new function to load in the beam model in ``WODEN/wodenpy/primary_beam/use_everybeam.py`` (similar to ``load_OSKAR_telescope`` or ``load_LOFAR_telescope``).
- In ``WODEN/wodenpy/skymodel/read_fits_skymodel.read_fits_skymodel_chunks``, ensure the new beam model initiated correctly.


You need to update the following GPU functions:

- Ensure your beam model is added to ``beam_settings`` in ``WODEN/src/primary_beam.c``
- At the start of ``calculate_visibilities`` in ``WODEN/src/calculate_visibilities.cpp``, ensure ``use_twobeams, num_beams`` are set correctly for new beam.
- Wherever ``copy_CPU_beam_gains_to_GPU`` is called in ``WODEN/src/calculate_visibilities.cpp``, ensure that happens for your new beam model.
- Add your model to ``WODEN/src/source_components.cpp.free_beam_gains`` function to ensure the leakage terms are freed correctly.
