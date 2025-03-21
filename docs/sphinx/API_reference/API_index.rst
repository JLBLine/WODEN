API Reference
==============

Here you can find documentation of all the functions grouped by coding language
and/or function. ``python`` script documentation also includes a summary of all
arguments, which can be reproduced on the command line using ``script.py --help``.

.. note::  All functions that begin with ``RTS`` are either borrowed or adapted
   from the RTS (`Mitchell et al. 2008`_) calibration package, with permission
   from the original authors. All credit to the original authors. The RTS code
   can currently be found at the `RTS github`_.

.. _Mitchell et al. 2008: https://doi.org/10.1109/JSTSP.2008.2005327
.. _RTS github: https://github.com/ICRAR/mwa-RTS.git

``python`` code
-----------------

``wodenpy``

.. toctree::
  :maxdepth: 1
  
  python_code/wodenpy/array_layout
  python_code/wodenpy/observational
  python_code/wodenpy/phase_rotate
  python_code/wodenpy/primary_beam
  python_code/wodenpy/skymodel
  python_code/wodenpy/use_libwoden
  python_code/wodenpy/uvfits
  python_code/wodenpy/wodenpy_setup
  

``C`` code
----------------------
This refers to the ``C`` code that has no GPU equivalent; see below for more
details.

The precision of most functions is determined at compilation time, via the
following ``#ifdef`` statement found in ``WODEN/include/woden_precision_defs.h``:

.. code-block:: C

    #ifdef DOUBLE_PRECISION
    /*! If -DDOUBLE_PRECISION flag is added at compilation,
    then user_precision_t is set to double */
    typedef double user_precision_t;
    /*! If -DDOUBLE_PRECISION flag is added at compilation,
    then user_precision_complex_t is set to double _Complex */
    typedef double _Complex user_precision_complex_t;
    #else
    /*! If -DDOUBLE_PRECISION flag is NOT added at compilation,
    then user_precision_t defaults to float */
    typedef float user_precision_t;
    /*! If -DDOUBLE_PRECISION flag is NOT added at compilation,
    then user_precision_complex_t defaults to float _Complex */
    typedef float _Complex user_precision_complex_t;
    #endif

So where you see ``user_precision_t`` in the API, this is either a ``float`` or
``double``, and similarly ``user_precision_complex_t`` is either a ``float _Complex``
or ``double _Complex``.

.. toctree::
  :maxdepth: 1

  C_code/beam_settings
  C_code/call_everybeam_c
  C_code/check_compilation_flags
  C_code/constants
  C_code/logger
  C_code/hyperbeam_error
  C_code/shapelet_basis
  C_code/visibility_set
  C_code/woden
  C_code/woden_precision_defs
  C_code/woden_struct_defs


``GPU`` and ``C`` equivalent functions
-----------------------------------------
These functions call either ``C`` or ``GPU`` code, depending on whether the 
user wants to run on a GPU or not. The ``common`` functions listed here
are designed to minimise the amount of code dupilcation between the ``C`` and
``GPU`` code, either calling a ``cpu`` or ``gpu`` function depending on the
``do_gpu`` boolean.

``common`` code
****************

.. toctree::
  :maxdepth: 1

  common_code/calculate_visibilities_common.rst
  common_code/source_components_common.rst


``C`` code
**************

.. toctree::
  :maxdepth: 1

  CPU_code/calculate_visibilities_cpu
  CPU_code/fundamental_coords_cpu
  CPU_code/primary_beam_cpu
  CPU_code/source_components_cpu

``GPU`` code
**************

All GPU code is either compiled as ``CUDA`` or ``HIP`` code. This is decided at
compilation, depending on whether ``-D__NVCC__`` (for ``CUDA``) or ``-D__HIPCC__``
(for ``HIP``) was passed to the compiler. Various macros found in ``WODEN/include/gpu_macros.h``
are used to drop in the correct GPU calls depending on the language requested.

Similarly to ``C`` code, the precision of most functions is determined at
compilation time, via the following ``#ifdef`` statement found in
``WODEN/include/gpucomplex.h``:

.. code-block:: C

   #ifdef DOUBLE_PRECISION
   /*! If -DDOUBLE_PRECISION flag is added at compilation,
   then cuUserComplex is set to cuDoubleComplex */
   typedef cuDoubleComplex cuUserComplex;
   #else
   /*! If -DDOUBLE_PRECISION flag is NOTE added at compilation,
   then cuUserComplex is set to cuFloatComplex */
   typedef cuFloatComplex cuUserComplex;
   #endif

meaning that ``cuUserComplex`` in the API either means ``cuFloatComplex`` or
``cuDoubleComplex`` depending on compilation flags.

.. toctree::
  :maxdepth: 1

  GPU_code/calculate_visibilities_gpu
  GPU_code/fundamental_coords_gpu
  GPU_code/gpu_macros
  GPU_code/gpucomplex
  GPU_code/primary_beam_gpu
  GPU_code/source_components_gpu
