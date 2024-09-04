.. _everybeam insallation page: https://everybeam.readthedocs.io/en/latest/build-instructions.html

EveryBeam Testing and Development
======================================

.. toctree::
   :maxdepth: 2
   :caption: Click these links to go to notebooks of integration tests

   test_OSKAR-MWA.ipynb
   test_LOFAR.ipynb
   test_OSKAR_SKA.ipynb

Installing EveryBeam
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The build instructions for EveryBeam live on the `everybeam insallation page`_.

When I ran ``cmake``, I did::

   cmake .. -DCMAKE_INSTALL_PREFIX=/home/jack-line/software/install \
     -DBUILD_WITH_PYTHON=ON


so I would know where the installation went. I also did all of this inside a conda environment, but I think ``cmake`` still found the system Python. Once installed, I had to do::
    
   export PYTHONPATH=/home/jack-line/software/install/lib/python3.12/site-packages:$PYTHONPATH
   ln -s /home/jack-line/software/install/share/everybeam /home/jack-line/software/anaconda3/envs/everybeam/share/everybeam

which let me conda environment see everything it needed to. When I was running notebooks, which don't load stuff from system, only the conda environment, I had to do::

   conda install -c conda-forge libstdcxx-ng
   conda install hdf5

for certain things to work.


Work still to be done on EveryBeam in WODEN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- I've managed to read the antenna locations from a measurement set, needed to calculate ``u,v,w`` coords. But I haven't worked out where the latitude/longitude are stored in the measurement set. At the moment you have to put that in via command line arguments. Something to figure out.
- You can only run one ``WODEN`` frequency band at a time. This is because ``WODEN`` is designed to calculate the beam values on GPU, iteratively band by band. It's not the end of the world but a slight design flaw.
- Only certain frequencies exists in certain models. This means we don't need to call EveryBeam for every single frequency. So create some kind of cache or map for the existing beam models to save compute.
- Really needs optimisation. We calculate the beam values as we read in chunks of the sky model. This happens in parallel to calculating the visibilities on the GPU. So some investigation is needed into chunking size (controlled by ``--chunking_size`` option) / how many chunks we read in at once (``max_num_chunks=50`` internally hardcoded at the moment), to try and balance the CPU and GPU workloads. ``--chunking_size`` controls the maximum number of visibilities we calculate on the GPU in one go, and is limited by GPU RAM. ``max_num_chunks`` is how many of these chunks we read in at once, before sending them across to the GPU (which will iterate over those chunks). This is limited really in how long it takes to calculate the beam values on the CPU. Likely the most efficient thing to do here is to go fully CPU so we don't have a GPU sitting idle; obviously this will take plenty of development (although could be close to a copy/paste of GPU into CPU code, so only tricky bit would be to make CPU functions multi-threaded).
- Unit tests. Everything is an integration test at the moment. Best practice is to break down the testing into units.

Adding new beam models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You need to update the following Python functions to add a new EveryBeam model (in general, look for anywhere ``EB_OSKAR`` is already used, and make sure whatever is done with ``EB_OSKAR`` is also done with your new beam model):

- Update the help for ``--primary_beam`` in the parser defined in ``WODEN/wodenpy/wodenpy_setup/run_setup.get_parser`` to show the new option
- Add a new beamtype to the ``BeamTypes`` enum in ``WODEN/wodenpy/use_libwoden/beam_settings``. You must also update the equivalent ``e_beamtype`` in ``WODEN/include/woden_struct_defs.h``. E.g. OSKAR already has a value of ``EB_OSKAR``, so if adding an EveryBeam MWA model, you would add ``EB_MWA``.
- Make sure your beam model is added to ``woden_settings`` in the ```WODEN/wodenpy/use_libwoden/woden_settings.create_woden_settings``` function.
- Add your beamtype to the ``eb_beams`` list in ``WODEN/wodenpy/use_libwoden/skymodel_structs.setup_components``
- Add a new function to load in the beam model in ``WODEN/wodenpy/primary_beam/use_everybeam.py`` (similar to ``load_OSKAR_telescope`` or ``load_LOFAR_telescope``).
- In ``WODEN/wodenpy/skymodel/read_fits_skymodel.read_fits_skymodel_chunks``, ensure the new beam model initiated correctly.


You need to update the following GPU functions:

- Ensure your beam model is added to ``beam_settings`` in ``WODEN/src/primary_beam.c``
- At the start of ``calculate_visibilities`` in ``WODEN/src/calculate_visibilities.cpp``, ensure ``use_twobeams, num_beams`` are set correctly for new beam.
- Wherever ``copy_CPU_beam_gains_to_GPU`` is called in ``WODEN/src/calculate_visibilities.cpp``, ensure that happens for your new beam model.
- Add your model to ``WODEN/src/source_components.cpp.free_beam_gains`` function to ensure the leakage terms are freed correctly.