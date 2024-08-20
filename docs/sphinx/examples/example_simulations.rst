.. _WODEN demonstrated via examples:

``WODEN`` demonstrated via examples
=====================================

Below are a number of example simulations with various settings. You can find and run these examples in the ``WODEN/examples`` directory. I'll let you know how much storage you'll need for each set of simulation and images (things can get large with radio data).

Two of the sky models are large (a total of about 100 MB), so instead of bundling them into the github, I've added links into the relevant instructions to download them.

.. note:: For all simulation times reported in the below, I used a single NVIDIA GeForce GTX 1080 Ti with 12 GB of RAM.

The examples are:

.. toctree::

   fornaxA_sim
   MWA_EoR1_sim
   eda2_haslam_sim
   dipole_ampflags
   polarisation_examples.ipynb

**Fornax A simulation** - two examples that compare a point/Gaussian model to a shapelet model. This serves as a general introduction on how to simulate an MWA observation using a ``metafits`` file, with a few extra commands. It also serves as a comparison of running with ``woden_float`` and ``woden_double``.

**MWA EoR1 simulation** - demonstrates using a larger (>300,000) source catalogue

**EDA2 Haslam Map simulation** - this demonstrates using ``WODEN`` without a ``metafits`` file, using a text file to describe the array layout, and using the ``EDA2`` beam.

**MWA dipole amplitudes and flags** - this demonstrates using ``WODEN`` when supplying either dipole flags or amplitudes to the MWA FEE beam (via ``hyperbeam``). *Note, currently this only works with the FEE beam, not the analytic MWA beam*

**Polarisation** - demonstrates how to get create polarised sky models, simulate them, and gives examples of how to inspect the results. 

.. warning:: If you have a GPU with small amounts of RAM (say 2GB) some of these simulations won't work to DOUBLE precision, you won't have enough memory. You can add the ``--precision=float`` argument to switch to a lower memory requirement (for the loss of accuracy).
