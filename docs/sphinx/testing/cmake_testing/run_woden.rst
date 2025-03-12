``run_woden``
=========================
Tests for the functions in ``WODEN/src/run_woden.c``. This is basically
``main`` for the dynamic library. 

test_run_woden.c
******************
This is really a token test, as multiple tests are run via ``run_woden.py``
when testing scripts. But this test serves to check everything runs from the 
``C/GPU`` code side. Good for debugging if someone breaks something on the
``Python`` side.

We simply run the same array setup as as when testing ``calculate_visibilities``,
with the interpolated MWA FEE beam model via ``hyperdrive``. One major difference
to the test in ``calculate_visibilities`` is that we run three coarse bands.
We stick a single point source at zenith, and check the output visibilities
by using stored primary beam gains to calculate the expected visibilities.
As we calculate three bands at three different frequencies, we check against
three different sets of beam gains, given the frequency dependence of the
primary beam.