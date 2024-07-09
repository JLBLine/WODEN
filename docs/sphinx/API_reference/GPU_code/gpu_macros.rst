.. _`astroio`: https://github.com/PaCER-BLINK-Project/astroio
.. _`Sokolowski et al. 2024`: https://ui.adsabs.harvard.edu/abs/2024arXiv240513478S/abstract

``gpu_macros``
==============

.. note::
   These macros were provided by Marcin Sokolowski from the PaCER BLINK `astroio`_ package. See also `Sokolowski et al. 2024`_ for an application of these macros.

API documentation for ``gpu_macros.h``. 

.. doxygenfile:: gpu_macros.h
   :project: WODEN

Macros
---------
See below for a table of the macros employed in ``gpucomplex.h``. Note that many
of the GPU functions are wrapped in the ``GPUErrorCheck`` function (documented 
above as ``docGPUErrorCheck``), which checks for errors and exits with a message
including the line number and file name if an error is detected. If there is error
checking in the macro it is listed below. Again, setting ``-D__NVCC__`` or
``-D__HIPCC__`` at compilations determines whether ``CUDA`` or ``HIP`` functions
are used.

.. list-table:: GPU function macros
   :widths: 30 30 30 10
   :header-rows: 1

   * - Macro
     - ``__NVCC__``
     - ``__HIPCC__``
     - Error wrapped
   * - gpuMalloc
     - cudaMalloc
     - hipMalloc
     - Yes
   * - gpuHostAlloc
     - cudaHostAlloc
     - hipHostMalloc
     - Yes
   * - gpuHostAllocDefault
     - cudaHostAllocDefault
     - 0
     - No
   * - gpuMemcpy
     - cudaMemcpy
     - hipMemcpy
     - Yes
   * - gpuMemcpyAsync
     - cudaMemcpyAsync
     - hipMemcpyAsync
     - Yes
   * - gpuMemset
     - cudaMemset
     - hipMemset
     - Yes
   * - gpuDeviceSynchronize
     - cudaDeviceSynchronize
     - hipDeviceSynchronize
     - No
   * - gpuMemcpyDeviceToHost
     - cudaMemcpyDeviceToHost
     - hipMemcpyDeviceToHost
     - No
   * - gpuMemcpyHostToDevice
     - cudaMemcpyHostToDevice
     - hipMemcpyHostToDevice
     - No
   * - gpuMemcpyDeviceToDevice
     - cudaMemcpyDeviceToDevice
     - hipMemcpyDeviceToDevice
     - No
   * - gpuFree
     - cudaFree
     - hipFree
     - Yes
   * - gpuHostFree
     - cudaFreeHost
     - hipHostFree
     - Yes
   * - gpuStream_t
     - cudaStream_t
     - hipStream_t
     - No
   * - gpuStreamCreate
     - cudaStreamCreate
     - hipStreamCreate
     - Yes
   * - gpuStreamDestroy
     - cudaStreamDestroy
     - hipStreamDestroy
     - Yes
   * - gpuEventCreate
     - cudaEventCreate
     - hipEventCreate
     - Yes
   * - gpuGetDeviceCount
     - cudaGetDeviceCount
     - hipGetDeviceCount
     - Yes
   * - gpuGetLastError
     - cudaGetLastError
     - hipGetLastError
     - No
   * - gpuMemGetInfo
     - cudaMemGetInfo
     - hipMemGetInfo
     - Yes
   * - gpuMallocHost
     - cudaMallocHost
     - hipHostMalloc
     - Yes
   * - gpuFreeHost
     - cudaFreeHost
     - hipFreeHost
     - Yes
   * - gpuGetDeviceProperties
     - cudaGetDeviceProperties
     - hipGetDeviceProperties
     - Yes
   * - gpuDeviceProp
     - cudaDeviceProp
     - hipDeviceProp_t
     - No
   * - gpuPeekAtLastError
     - cudaPeekAtLastError
     - hipPeekAtLastError
     - No


.. list-table:: Complex number operation
   :widths: 33 33 33
   :header-rows: 1

   * - Macro
     - ``__NVCC__``
     - ``__HIPCC__``
   * - gpuCreal
     - cuCreal
     - hipCreal
   * - gpuCrealf
     - cuCrealf
     - hipCrealf
   * - gpuCimag
     - cuCimag
     - hipCimag
   * - gpuCimagf
     - cuCimagf
     - hipCimagf
   * - gpuCadd
     - cuCadd
     - hipCadd
   * - gpuCmul
     - cuCmul
     - hipCmul
   * - gpuCdiv
     - cuCdiv
     - hipCdiv
   * - gpuConj
     - cuConj
     - hipConj
   * - gpuCsub
     - cuCsub
     - hipCsub
   * - gpuCabs
     - cuCabs
     - hipCabs
   * - gpuCaddf
     - cuCaddf
     - hipCaddf
   * - gpuCsubf
     - cuCsubf
     - hipCsubf
   * - gpuCmulf
     - cuCmulf
     - hipCmulf
   * - gpuCdivf
     - cuCdivf
     - hipCdivf
   * - gpuDoubleComplex
     - cuDoubleComplex
     - hipDoubleComplex
   * - gpuFloatComplex
     - cuFloatComplex
     - hipFloatComplex
   * - make_gpuDoubleComplex
     - make_cuDoubleComplex
     - make_hipDoubleComplex
   * - make_gpuFloatComplex
     - make_cuFloatComplex
     - make_hipFloatComplex