``primary_beam_cuda``
=========================
Tests for the functions in ``WODEN/src/primary_beam_cuda.cu``. These functions
calculate the beam responses for the EDA2 and Gaussian beam models.

test_gaussian_beam.c
*********************************
This calls ``primary_beam_cuda::test_kern_gaussian_beam``, which in turn
tests ``primary_beam_cuda::kern_gaussian_beam``, the kernel that calculates
the Gaussian primary beam response. As a Gaussian is an easy function to
calculate, I've setup tests that calculate a north-south and east-west strip
of the beam response, and then compare that to a 1D Gaussian calculation.

As ``kern_gaussian_beam`` just takes in *l,m* coords, these tests just generate
100 *l,m* coords that span from -1 to +1. The tests check whether the kernel
produces the expected coordinates in the *l* and *m* strips, as well as changing
with frequency as expected, by testing 5 input frequencies with a given
reference frequency. For each input frequency :math:`\nu`, the output is
checked against the following calculations:

 - When setting *m* = 0, assert gain = :math:`\exp\left[-\frac{1}{2} \left( \frac{l}{\sigma} \frac{\nu}{\nu_0} \right)^2 \right]`
 - When setting *l* = 0, assert gain = :math:`\exp\left[-\frac{1}{2} \left( \frac{m}{\sigma} \frac{\nu}{\nu_0} \right)^2 \right]`

where :math:`\nu_0` is the reference frequency, and :math:`\sigma_0` the std of
the Gaussian in terms of *lm* coords. The real parts of the beam response are
tested to be within an absolute tolerance of 1e-7, and the imaginary parts
to be equal to zero.

test_analytic_dipole_beam.c
***********************************
This calls ``primary_beam_cuda::test_analytic_dipole_beam``, which in turn
tests ``primary_beam_cuda::calculate_analytic_dipole_beam``, code that copies
az/za angles into GPU memory, calculates an analytic dipole response toward
those directions, and then frees the az/za coords from GPU memory.

Nothing exiting in this test, just call the function for 25 directions on
the sky, for two time steps and two frequencies (a total of 100 beam calculations),
and check that the real beam gains match stored expected values, and the imaginary
values equal zero.
