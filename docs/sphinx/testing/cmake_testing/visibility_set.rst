``visibility_set``
=========================
Tests for the functions in ``WODEN/src/visibility_set.c``. These functions handle
a ``visibility_set_t`` struct, which holds the output visibilities. Functions
here include mallocing, filling, and freeing attributes.

test_fill_timefreq_visibility_set.c
*****************************************
This calls ``visibility_set::fill_timefreq_visibility_set``, which uses a
populated ``woden_settings_t`` struct to fill in the following attributes in
a ``visibility_set_t`` struct::

  visibility_set->allsteps_sha0s
  visibility_set->allsteps_cha0s
  visibility_set->allsteps_lsts
  visibility_set->allsteps_wavelengths

Where ``sha0`` and ``cha0`` are sine and cosine of the hour angle of the phase
centre, respectively. This test runs with three time steps with LSTs of
:math:`0`, :math:`\pi/6`, :math:`\pi/4`, so that ``allsteps_sha0s`` and
``allsteps_cha0s`` can be analytically predicted. ``allsteps_wavelengths`` are
calculated using input frequencies, which are set to :math:`c/2`, :math:`3c/4`,
and :math:`c`, meaning the expected output wavelengths should be :math:`2`, :math:`4/3`,
and :math:`1`. The third let's us test the accuracy of the wavelength calculation.

All angles and LSTs are stored at 64bit precision, so both FLOAT and DOUBLE
code versions are tested to within an absolute tolerance of 1e-15. The precision
of ``allsteps_wavelengths`` is set by the user, and is tested to within 1e-7
for FLOAT and 1e-15 for DOUBLE.

test_malloc_and_free.c
*****************************************
Very basic test of ``visibility_set::setup_visibility_set``,
``visibility_set::free_visi_set_inputs``, and ``visibility_set::free_visi_set_outputs``,
which are functions that either ``malloc`` or ``free`` specific attributes in a
``visibility_set_t`` struct. Tests by calling each function and checking that
the following attributes are NOT a NULL if a ``malloc`` was called, and the correct
attributes are NULL if ``free`` was called.

test_write_visi_set_binary.c
*****************************************
Tests ``visibility_set::write_visi_set_binary``, which writes out the contents
of a ``visibility_set_t`` struct to a binary file. Tests by filling a
``visibility_set_t`` struct with some simple non-repeating values, and then
uses ``write_visi_set_binary``. Then reads that binary file in and checks the
contents matched what was in the ``visibility_set_t`` struct.

test_write_visi_set_text.c
*****************************************
Tests ``visibility_set::write_visi_set_text``, which writes out a subset of
the contents of a ``visibility_set_t`` struct to a text file. This test
works similarly to :ref:`test_write_visi_set_binary.c`, by calling
``write_visi_set_text`` with a known set of inputs, and checking the text file
it writes out contains the known inputs.
