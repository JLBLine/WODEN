``create_sky_model``
=========================
Tests for the functions in ``WODEN/src/create_sky_model.c``. The functions
read in the sky model from a text file, and crop anything below the horizon.

test_read_source_catalogue.c
*********************************
``create_sky_model::read_source_catalogue`` reads in the sky model from a text
file. Each test reads in a sky model, uses ``read_source_catalogue`` to
create a sky model, and tests the correct information has been read in.
``read_source_catalogue`` returns an integer as an error message (0 good, 1 bad),
so some text files below are written to cause failure. See the table below
for each test sky model and the expected result. The way each SOURCE and it's
associated COMPONENTs are stored can affect the way the sky model is cropped,
so all tests check that generated ``source_catalogue_t`` struct is structured
correctly.

.. list-table::
   :widths: 25 50 25
   :header-rows: 1

   * - Sky model
     - Test case
     - Test outcomes
   * - srclist_no-comp_numbers.txt
     - Is missing the line that contains number of COMPONENTs in the SOURCE.
     - Check fails
   * - srclist_badcoeff.txt
     - Has a bad SCOEFF line where one number is sdfasdfasdfasdfasdfasf
     - Check fails
   * - srclist_badspell.txt
     - Contains an incorrect spelling of COMPONENT
     - Check fails
   * - srclist_singlegauss.txt
     - Contains a single GAUSSIAN SOURCE
     - Check sky model values
   * - srclist_singlepoint.txt
     - Contains a single POINT SOURCE
     - Check sky model values
   * - srclist_comment.txt
     - Contains a commented line (and a single POINT SOURCE)
     - Check sky model values
   * - srclist_singleshape.txt
     - Contains a single SHAPELET SOURCE
     - Check sky model values
   * - srclist_empty_line.txt
     - Contains an empty line  (and a single POINT SOURCE)
     - Check sky model values
   * - srclist_threecomponents.txt
     - Contains one SOURCE with three COMPONENTs
     - Check sky model values
   * - srclist_mulitple_source-components.txt
     - Contains multiple SOURCEs each with multiple COMPONENTs
     - Check sky model values
   * - srclist_threesources.txt
     - Contains multiple SOURCEs each with a single COMPONENT
     - Check sky model values


test_horizon_test.c
*********************************
``create_sky_mode::horizon_test`` takes information on a single COMPONENT of a
SOURCE and tests whether it is above or below the horizon. Depending on
whether we are cropping the sky model by SOURCE or by COMPONENT, it updates
various counters that effect the sky model cropping. (cropping by SOURCE
throws away the whole SOURCE if one COMPONENT is below the horizon, cropping
by COMPONENT only throws away the COMPONENTs below the horizon). SHAPELETs are complicating
factors as a single position can match multiple basis function parameters (of
any length) so ``horizon_test`` does some logic to count how many SHAPELET
parameters are being retained.

First set of tests check that cropping on POINT/GAUSSIAN type COMPONENTs work for
both cropping by SOURCE and cropping by COMPONENT. The second set of tests check
that when cropping by COMPONENT, and we have multiple SHAPELET COMPONENTs, that
the correct number of SHAPELET basis functions parameters are retained.

test_crop_sky_model.c
*********************************
``create_sky_mode::crop_sky_model`` calculates the azimuth / zenith angle of
all COMPONENTs in a sky model (for the first LST step of a simulation), and then
crops out either SOURCEs or COMPONENTs that are below the horizon (by using
``horizon_test``). Once cropped, it then calculates the az/za for all time
steps in the simulation for the surviving COMPONENTs.

The tests here are split across 4 sky models - just POINTs, just GAUSSIANs,
just SHAPELETs, and a mix of all three. For each sky model type, the tests are
run for both the cropping by SOURCE and cropping by COMPONENT case. All tests
are run for four different LSTs - a total of 32 tests. The RA/Dec of the
COMPONENTs and the LSTs are chosen in such a way that different combinations
should result in different COMPONENT/SOURCEs being discarded. All tests check
that the correct COMPONENT/SOURCEs are retained (including the SHAPELET basis
function parameters), and that the correct az/za are calculated for each time
step. The sky model setup and expected results are stored in ``test_crop_sky_model.h``.
