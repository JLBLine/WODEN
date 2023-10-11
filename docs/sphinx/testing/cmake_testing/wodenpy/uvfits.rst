``uvfits``
=======================


test_RTS_encoding.py
*******************************************************
Tests the ``RTS_encode_baseline`` function, which should take two antenna
numbers and create the ``BASELINE`` number as per AIPS uvfits file convention.
Tests by running with four antenna pairs, and ensuring the output values match
expectation.

Secondly, tests the function ``RTS_decode_baseline``, which should separate
the encoded BASELINE number back into two antenna numbers. Test by decoding the
same four BASELINE numbers and ensuring the correct antenna numbers are found.

test_make_antenna_table.py
*******************************************************
Tests the ``make_antenna_table`` function, which should create
the antenna table that goes into a uvfits file. Test by giving it a
known set of input parameters and checking those values end up in
the correct location and format of output antenna table. The following parameters
are tested as correct:

.. code-block:: python

   ant_table.data['ANNAME']
   ant_table.data['STABXYZ']
   ant_table.data['ORBPARM']
   ant_table.data['NOSTA']
   ant_table.data['MNTSTA']
   ant_table.data['STAXOF']
   ant_table.data['POLTYA']
   ant_table.data['POLAA']
   ant_table.data['POLCALA']
   ant_table.data['POLTYB']
   ant_table.data['POLAB']
   ant_table.data['POLCALB']
   ant_table.header['ARRAYX']
   ant_table.header['ARRAYY']
   ant_table.header['ARRAYZ']
   ant_table.header['FREQ']
   ant_table.header['GSTIA0']
   ant_table.header['DEGPDY']
   ant_table.header['UT1UTC']
   ant_table.header['XYZHAND']
   ant_table.header['FRAME']
   ant_table.header['RDATE']
   ant_table.header['TIMSYS']
   ant_table.header['ARRNAM']
   ant_table.header['NUMORB']
   ant_table.header['NOPCAL']
   ant_table.header['POLTYPE']
   ant_table.header['CREATOR']

test_create_uvfits.py
*******************************************************
Tests the ``create_uvfits`` function, which should take a whole
heap of inputs and write out a uvfits. Test by running function with a known
set of inputs, reading in the created file, and checking contents match
expectations. Along checking all the same parameters in the antenna table as
checked by :ref:`test_make_antenna_table.py`, the following parameters are
checked against the inputs:

.. code-block:: python

    data_table.data['UU']
    data_table.data['VV']
    data_table.data['WW']
    data_table.data['BASELINE']
    ##Astropy automatically adds the header value to the DATE array,
    ##so need to subtract before comparison
    data_table.data['DATE'] - data_table.header['PZERO4']
    ##Check the actual visisbility values are correct
    data_table.data.data
    data_table.header['CTYPE2']
    data_table.header['CRVAL2']
    data_table.header['CRPIX2']
    data_table.header['CDELT2']
    data_table.header['CTYPE3']
    data_table.header['CRVAL3']
    data_table.header['CRPIX3']
    data_table.header['CDELT3']
    data_table.header['CTYPE4']
    data_table.header['CRVAL4']
    data_table.header['CRPIX4']
    data_table.header['CDELT4']
    data_table.header['CTYPE5']
    data_table.header['CRVAL5']
    data_table.header['CRPIX5']
    data_table.header['CDELT5']
    data_table.header['CTYPE6']
    data_table.header['CRVAL6']
    data_table.header['CRPIX6']
    data_table.header['CDELT6']
    data_table.header['PSCAL1']
    data_table.header['PZERO1']
    data_table.header['PSCAL2']
    data_table.header['PZERO2']
    data_table.header['PSCAL3']
    data_table.header['PZERO3']
    data_table.header['PSCAL4']
    data_table.header['PZERO4']
    data_table.header['PSCAL5']
    data_table.header['PZERO5']
    data_table.header['OBJECT']
    data_table.header['OBSRA']
    data_table.header['OBSDEC']
    data_table.header['GITLABEL']
    data_table.header['TELESCOP']
    data_table.header['LAT']
    data_table.header['LON']
    data_table.header['ALT']
    data_table.header['INSTRUME']

test_make_baseline_date_arrays.py
*******************************************************
Tests the ``rw.make_baseline_date_arrays`` function, which should make
the DATE and BASELINE arrays that are needed to populate a uvfits file. Test
by giving the function a known date string, number of antennas, number of time
steps, and time resolution, and checking the output arrays match expectations.

test_read_uvfits_into_pyuvdata.py
*******************************************************
This tests the absolute minimal compliance with `pyuvdata`_. The test calls
``rw.create_uvfits`` with a set of dummy input variables to make a file
called ``unittest_example.uvfits``. It then simply checks that the following
lines don't throw an error:

.. code-block:: python

  from pyuvdata import UVData
  UV = UVData()
  UV.read('unittest_example.uvfits')

This really only tests that the correct keywords and arrays are present in the
output ``unittest_example.uvfits`` to a level that appeases ``pyuvdata``.
The test is setup to skip if the user has not installed ``pyuvdata``.

.. _`pyuvdata`: https://pyuvdata.readthedocs.io/en/latest/index.html