``read_yaml_skymodel``
=========================
Tests for functions in ``WODEN/src/read_yaml_skymodel.c``. These functions
read in the sky model from a ``hyperdive``-style text file.

This sky model format allows POINT, GAUSSIAN, and SHAPELET components types as
well as POWER_LAW, CURVED_POWER_LAW, and LIST type flux density behaviours

read_yaml_skymodel.c
*********************************
``read_text_skymodel::read_text_skymodel`` reads in the sky model from a YAML
file. Each test reads in a sky model, uses ``read_yaml_skymodel`` to
create a sky model, and tests the correct information has been read in.
``read_yaml_skymodel`` returns an integer as an error message (0 good, 1 bad),
so some text files below are written to cause failure. See the table below
for each test sky model and the expected result.

.. list-table::
   :widths: 25 25 25
   :header-rows: 1

   * - Sky model
     - Test case
     - Test outcomes
   * - srclist_comment.yaml
     - Model includes commented lines
     - Check passes and values
   * - srclist_empty_line.yaml
     - Model contains empty lines
     - Check passes and values
   * - not_a_file.yaml
     - Path to missing file
     - Check fails
   * - srclist_singlepoint_curve.yaml
     - Single POINT component with CURVED_POWER_LAW
     - Check values
   * - srclist_singlepoint_list.yaml
     - Single POINT component with LIST
     - Check values
   * - srclist_singlepoint_power.yaml
     - Single POINT component with POWER_LAW
     - Check values
   * - srclist_singlegauss_curve.yaml
     - Single GAUSSIAN component with CURVED_POWER_LAW
     - Check values
   * - srclist_singlegauss_list.yaml
     - Single GAUSSIAN component with LIST
     - Check values
   * - srclist_singlegauss_power.yaml
     - Single GAUSSIAN component with POWER_LAW
     - Check values
   * - srclist_singleshape_curve.yaml
     - Single SHAPELET component with CURVED_POWER_LAW
     - Check values
   * - srclist_singleshape_list.yaml
     - Single SHAPELET component with LIST
     - Check values
   * - srclist_singleshape_power.yaml
     - Single SHAPELET component with POWER_LAW
     - Check values
   * - srclist_threecomponents_curve.yaml
     - Single SOURCE with three COMPONENTs with CURVED_POWER_LAW
     - Check values
   * - srclist_threecomponents_list.yaml
     - Single SOURCE with three COMPONENTs with LIST
     - Check values
   * - srclist_threecomponents_power.yaml
     - Single SOURCE with three COMPONENTs with POWER_LAW
     - Check values
   * - srclist_threesources_curve.yaml
     - Three SOURCEs with one COMPONENT with CURVED_POWER_LAW
     - Check values
   * - srclist_threesources_list.yaml
     - Three SOURCEs with one COMPONENT with LIST
     - Check values
   * - srclist_threesources_power.yaml
     - Three SOURCEs with one COMPONENT with POWER_LAW
     - Check values
   * - srclist_mulitple_source-components.yaml
     - Complete mix and match of COMPONENT and FLUX types, in various combinations
     - Check values
