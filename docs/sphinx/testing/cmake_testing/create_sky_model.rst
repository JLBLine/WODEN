``create_sky_model``
=========================
Tests for the functions in ``WODEN/src/create_sky_model.c``. The functions
read in the sky model from a text file, and crop anything below the horizon.

test_read_source_catalogue.c
*********************************
``create_sky_model::read_source_catalogue`` reads in the sky model from either a text
file (in the old ``WODEN``-style model format) or a YAML file (in the current
``hyperdrive``-style format), by calling ``read_text_skymodel::read_text_skymodel``
and ``read_yaml_skymodel::read_yaml_skymodel`` respectively. These functions
are tested extensively elsewhere, so only a few tests are run here, by reading
in the following sky models:

.. list-table::
   :widths: 25 25 25
   :header-rows: 1

   * - Sky model
     - Test case
     - Test outcomes
   * - not_a_file.txt
     - Link to a text file that doesn't exist
     - Check fails
   * - not_a_file.yaml
     - Link to a YAML file that doesn't exist
     - Check fails
   * - what_are_you_trying_to_pull.dat
     - Link to a file that doesn't end in either ``.txt`` or ``.yaml``
     - Check fails
   * - srclist_singlepoint.yaml
     - Contains a single POINT SOURCE in YAML format
     - Check ``read_yaml_skymodel`` is called and returns correct values
   * - srclist_singlepoint_power.txt
     - Contains a single POINT SOURCE in text format
     - Check ``read_text_skymodel`` is called and returns correct values

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

Furthermore, each SOURCE contains two each of the POWER_LAW, CURVED_POWER_LAW,
and LIST type flux behaviour for all POINT, GAUSSIAN, and SHAPELET models, for
a total of up 18 COMPONENTS per SOURCE. For each test, all flux entries and
associated values are checked as being copied from the original sources
to the cropped source correctly.

The azimuth and zenith angle outputs are tested to match expectations to within
an absolute tolerance of 1e-6 for the FLOAT compiled code, and 1e-12 for the
DOUBLE.
