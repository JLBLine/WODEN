`FEE_primary_beam`
==================
.. note::  All functions that begin with `RTS` are either borrowed or adapted
   from the RTS (`Mitchell et al. 2008`_) calibration package, with permission
   from the original authors. All credit to the original authors. The RTS code
   can currently be found at the `RTS github`_.

The MWA Fully Embedded Element primary beam model (`Sokolowski et al. 2017`_) is
at the time of writing the most accurate and advanced MWA primary beam model.
The model is stored in spherical harmonic coefficients, at a frequency
resolution of 1.28MHz, in an `hdf5` file called
`mwa_full_embedded_element_pattern.h5`. I have done my best to document the RTS
code I have used. The RTS code makes use of the Legendre Polynomial code
`available here`_ by John Burkardt, and is included in the `WODEN` distribution
as `legendre_polynomial.c`.

.. _Mitchell et al. 2008: https://doi.org/10.1109/JSTSP.2008.2005327
.. _RTS github: https://github.com/ICRAR/mwa-RTS.git
.. _Sokolowski et al. 2017: https://doi.org/10.1017/pasa.2017.54
.. _available here: https://people.sc.fsu.edu/~jburkardt/c_src/laguerre_polynomial/laguerre_polynomial.html

.. doxygenfile:: FEE_primary_beam.h
   :project: WODEN
