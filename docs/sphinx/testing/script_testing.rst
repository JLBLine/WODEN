Testing installation via scripts
=================================

This is a straight-forward way to check your installation is working; just run some simple small simulations to check various functionality. Easiest way to check is to create images of the results. I've included a second set of scripts to convert the outputs to measurement sets and image them using ``WSClean``.

Running the simulations
------------------------

Once you've built ``WODEN``, and you have variables defined by calling ``WODEN/build/init_WODEN.sh``,
navigate to ``WODEN/test_installation``. To run all the tests immediately, you'll need to have obtained the defined the environment variable ``MWA_FEE_HDF5`` (see :ref:`Post compilation (optional)` for instructions on how to define that). To run all scripted test (including MWA FEE simulations)::

  $ cd WODEN/test_installation
  $ ./run_all_simulations.sh

This should rattle through a number of tests, with various combinations of component types (point, Gaussian, or shapelet), running with and without metafits, and different primary beams. If you don't want to run MWA FEE simulations, you can run::

  $ ./run_all_but_MWAFEE_simulations.sh

If you then want to run MWA FEE tests alone at a later date, you can run::

  $ ./run_only_MWAFEE_simulations.sh

If you want to incrementally run through tests, you can navigate through the ``point_models``, ``gauss_models``, ``shapelet_models``, and ``combo_models`` directories to run each set individually.

Imaging the simulations
------------------------

Dependencies
^^^^^^^^^^^^^

To image the tests, you'll need a way to convert the uvfits files to measurement sets, and then image them. Which means dependencies (again).

+ **CASA** - https://casa.nrao.edu/casa_obtaining.shtml. Download an appropriate tarball (at the time 6.2 is the best version at it's ``python3``) and decompress it::

  $ wget https://casa.nrao.edu/download/distro/casa/release/rhel/casa-6.2.0-124.tar.xz
  $ tar -xvf casa-6.2.0-124.tar.xz

  That's it, it doesn't need installation
+ **WSClean** - https://wsclean.readthedocs.io/en/latest/installation.html. Head to this link to find out how to install ``WSClean``. You can of course use any other CLEANing software you want, but ``WSClean`` is excellent.

I'm assuming if you want to simulate interferometric data, you'll have some kind of FITS file imager already, but if not, ``DS9`` is a good place to start - https://sites.google.com/cfa.harvard.edu/saoimageds9.

Imaging scripts
^^^^^^^^^^^^^^^^

To use the imaging scripts, you **must** set an environment variable ``CASA_DIR`` to point towards the ``casa/bin`` of where you installed your ``casa``. For example, I did this::

  $ export CASA_DIR="/usr/local/casa-6.2.0-124/bin"

This will enable the scripts to convert the uvfits files to measurement sets. Once that variable is set, you can either image all the test outputs with::

  $ ./run_all_imaging.sh

or run the following as required::
  $ ./run_all_but_MWAFEE_imaging.sh
  $ ./run_only_MWAFEE_imaging.sh

Expected outcomes
------------------------


Deleting test outputs
------------------------
If you don't want a bunch of files hanging around on your system for no reason, just run::

  $ ./delete_sim_outputs.sh
  $ ./delete_images.sh

which will nuke the outputs for you.
