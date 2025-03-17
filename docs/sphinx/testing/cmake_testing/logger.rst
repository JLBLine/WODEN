``logger``
=========================
Tests for the functions in ``WODEN/src/logger.c``. ``logger.c`` is has
simple logging functionality that is designed to feed into the controlling
Python logger in ``run_woden.py``. It can work as a standalone logger however.
The WODEN logger is designed to centralise certain errors, so we check here
that we can grab ``hyperbeam`` errors and ``GPU`` errors correctly.

test_hyperbeam_error.c
***************************
This test is fun. Test that we can error hyperbeam and log it correctly. Do this
by calling ``new_fee_beam`` with a bad path to the ``hdf5``. This should cause
and exit, so we redirect the ``exit`` function using some ``setjmp`` shenanigans.
This allows ``Unity`` to check the error message by reading the log, without
exiting and calling things a FAIL.

test_log_gpu_errors.c
***************************
Again, a fun test. This time we test that we can log GPU errors correctly. To do
this, I wrote some bad GPU code in ``bad_gpu_code.cpp``, which tries to copy a
``nullptr`` to a ``num_elements`` double array on the GPU. This time we redirect
the ``exit`` function by ``#define exit mock_exit`` and setting the ``mock_exit``
as something benign. The test asserts that ``GPU ERROR cudaMemcpy`` is in the log.

test_woden_logger.c
***************************
This is the most basic test, just checks we can call the ``logger`` functions,
and that we can set the ``log_callback`` to a different message.
