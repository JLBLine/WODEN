API Reference
==============

Here you can find documentation of all the functions grouped by coding language.
``python`` script documentation also includes a summary of all arguments, which
can be reproduced on the command line using ``script.py --help``.

.. note::  All functions that begin with ``RTS`` are either borrowed or adapted
   from the RTS (`Mitchell et al. 2008`_) calibration package, with permission
   from the original authors. All credit to the original authors. The RTS code
   can currently be found at the `RTS github`_.

.. _Mitchell et al. 2008: https://doi.org/10.1109/JSTSP.2008.2005327
.. _RTS github: https://github.com/ICRAR/mwa-RTS.git

``python`` code
-----------------

.. toctree::
  :maxdepth: 1

  python_code/run_woden
  python_code/convert_WSClean_list_to_WODEN
  python_code/uv2ms

``C`` code
----------------------

.. toctree::
  :maxdepth: 1

  C_code/array_layout
  C_code/chunk_sky_model
  C_code/constants
  C_code/create_sky_model
  C_code/FEE_primary_beam
  C_code/primary_beam
  C_code/print_help
  C_code/shapelet_basis
  C_code/visibility_set
  C_code/woden_settings
  C_code/woden_struct_defs

``CUDA`` code
-------------------------

.. toctree::
  :maxdepth: 1

  CUDA_code/calculate_visibilities
  CUDA_code/cudacheck
  CUDA_code/cudacomplex
  CUDA_code/FEE_primary_beam_cuda
  CUDA_code/fundamental_coords
  CUDA_code/primary_beam_cuda
  CUDA_code/source_components
