``read_text_skymodel``
=========================
Tests for functions in ``WODEN/src/read_text_skymodel.c``. These functions
read in the sky model from a ``WODEN``-style text file.

.. warning:: From version 1.4 onward, the ``hyperdrive``-style YAML sky model format is preferred over the ``WODEN`` format. The ``WODEN`` format only supports a power-law spectral behaviour and will not be developed any further, whereas the ``hyperdrive`` sky model also supports curved power law and list style flux entries, and will be developed against going forward in ``WODEN``.

The text functionality will remain so any existing catalogues can be read in without extra effort to convert them.

read_text_skymodel.c
*********************************
``read_text_skymodel::read_text_skymodel`` reads in the sky model from a text
file. Each test reads in a sky model, uses ``read_text_skymodel`` to
create a sky model, and tests the correct information has been read in.
``read_text_skymodel`` returns an integer as an error message (0 good, 1 bad),
so some text files below are written to cause failure. See the table below
for each test sky model and the expected result. The way each SOURCE and it's
associated COMPONENTs are stored can affect the way the sky model is cropped,
so all tests check that generated ``source_catalogue_t`` struct is structured
correctly. For all floating point values, when compiling in FLOAT mode, test
asserts that values are within an absolute tolerance of 1e-7, and 1e-15 when
compiling in DOUBLE mode.

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
   * - srclist_linear.txt
     - Contains a single POINT COMPONENT with the LINEAR keyword specifying the SED
     - Check sky model values
