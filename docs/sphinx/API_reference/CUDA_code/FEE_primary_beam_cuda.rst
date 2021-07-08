``FEE_primary_beam_cuda``
=========================

The MWA Fully Embedded Element primary beam model (`Sokolowski et al. 2017`_) is
at the time of writing the most accurate and advanced MWA primary beam model.
The model is stored in spherical harmonic coefficients, at a frequency
resolution of 1.28MHz, in an ``hdf5`` file called
``mwa_full_embedded_element_pattern.h5``. I have done my best to document the RTS
code I have used. The RTS code makes use of the Legendre Polynomial code
`available here`_ by John Burkardt, and is included in the ``WODEN`` distribution
as ``legendre_polynomial.c``.

.. _Sokolowski et al. 2017: https://doi.org/10.1017/pasa.2017.54
.. _available here: https://people.sc.fsu.edu/~jburkardt/c_src/laguerre_polynomial/laguerre_polynomial.html

.. doxygenfile:: FEE_primary_beam_cuda.h
   :project: WODEN
