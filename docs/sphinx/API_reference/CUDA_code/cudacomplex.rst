``cudacomplex``
================

This header is in the ``RTS``, and contains useful CUDA operators. All credit to
the original author, R. G. Edgar. I've added in ``double`` definitions, as well
as my ``cuUserComplex`` def, which allows ``float`` or ``double`` to be
selected during compilation using the flag ``--DDOUBLE_PRECISION`` (so even
though in the below it says ``typedef cuFloatComplex cuUserComplex``, this
depends on compilation).

.. doxygenfile:: cudacomplex.h
   :project: WODEN
