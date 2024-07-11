``gpucomplex``
================

This header is in the ``RTS``, and contains useful CUDA operators. All credit to
the original author, R. G. Edgar. I've added in ``double`` definitions, as well
as my `gpuUserComplex`` def, which allows ``float`` or ``double`` to be
selected during compilation using the flag ``--DDOUBLE_PRECISION`` (so even
though in the below it says ``typedef gpuFloatComplex gupuUserComplex``, this
depends on compilation). This code has also been updated with GPU macros, so
it can be used for both CUDA and HIP code in conjunction with the
``gpu_macros.h`` header.

.. doxygenfile:: gpucomplex.h
   :project: WODEN
