``run_woden``
=========================
Tests for the functions in ``WODEN/src/run_woden.py``. These functions handle:
parsing user arguments; calculating astronomical constants;
reading in variables from an MWA metafits if requested; writing out a ``.json``
file to input into the ``woden`` executable; calling the ``woden`` executable;
reading the binary file written out by ``woden``; creating a ``uvfits`` file;
tidying up after the simulator.

In the following, ``rw`` is short for ``run_woden.py``, and so ``rw.calc_jdcal``
means the function ``calc_jdcal`` in ``run_woden.py``.

test_command.py
*******************************************************
Tests the ``rw.command`` function, which should call things on the command line
from within ``python``. Test by running the command::

   $ echo cheese > example.txt

and then reading the word "cheese" out of the file "example.txt" that should
have been created.

test_calc_jdcal.py
*******************************************************
Tests the ``rw.calc_jdcal`` function, which should calculate the Julian Date and
split into a integer day and fractional day values. Just test by calling
``rw.calc_jdcal`` with two known date strings, and checking the output values
match expectations.

test_get_uvfits_date_and_position_constants.py
*******************************************************
Tests the ``rw.get_uvfits_date_and_position_constants`` function,
which should calculate the LST, GST0 (greenwich sidereal time at 0 hours
of the given date), DEGPDY (rotational speed of the Earth) and UT1UTC (
difference between UT1 and UTC) for a given Long/Lat/Height array location and
UTC date. Test with two combinations of different UTC/Long/Lat/Height and
check the returned values are as expected.

test_RTS_encoding.py
*******************************************************
Tests the ``rw.RTS_encode_baseline`` function, which should take two antenna
numbers and create the ``BASELINE`` number as per AIPS uvfits file convention.
Tests by running with four antenna pairs, and ensuring the output values match
expectation.

Secondly, tests the function ``rw.RTS_decode_baseline``, which should separate
the encoded BASELINE number back into two antenna numbers. Test by decoding the
same four BASELINE numbers and ensuring the correct antenna numbers are found.

test_make_antenna_table.py
*******************************************************
Tests the ``rw.make_antenna_table`` function, which should create
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
Tests the ``rw.create_uvfits`` function, which should take a whole
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

test_enh2xyz.py
*******************************************************
Tests the ``rw.enh2xyz`` function, which should calculate the local X,Y,Z coords
using the local east, north, height. Test using the cases where latitude is 0
and -30 deg, which have analytically predictable outcomes. This runs the same
test as descibred in :ref:`test_RTS_ENH2XYZ_local.c`.

test_load_data.py
*******************************************************
Tests the ``rw.load_data`` function, which should read in a binary
file as output by ``woden_float`` or ``woden_double``, into various arrays.
Test by writing out a binary file with known input params, reading in that
binary file using ``rw.load_data``, and comparing the inputs to outputs. The
test is run in both 32 and 64 bit precision.

test_write_json.py
*******************************************************
Test the ``rw.write_json`` function, which writes an input file to feed into
either ``woden_float`` or ``woden_double``. A number of tests are run, all of
which call ``rw.write_json`` using a minimum set of example input arguments.
The resulting ``.json`` is then read back in, and the following parameters are
checked as correct:

.. code-block:: python

   json_data['ra0']
   json_data['dec0']
   json_data['num_freqs']
   json_data['num_time_steps']
   json_data['cat_filename']
   json_data['time_res']
   json_data['frequency_resolution']
   json_data['chunking_size']
   json_data['jd_date']
   json_data['LST']
   json_data['array_layout']
   json_data['lowest_channel_freq']
   json_data['latitude']
   json_data['coarse_band_width']
   json_data['band_nums']

The following tests run with the following optional arguments:

 - ``test_write_gaussian_beam``: checks that extra arguments that control the Gaussian primary beam are written correctly
 - ``test_write_MWA_FEE_beam``: checks that extra arguments that control the MWA FEE beam are written correctly
 - ``test_write_MWA_FEE_beam_interp``: checks that extra arguments that control the interpolated MWA FEE beam are written correctly
 - ``test_write_MWA_analy_beam``: checks that extra arguments that control analytic MWA beam are written correctly
 - ``test_write_EDA2_beam``: checks that extra arguments that control the EDA2 beam are written correctly
 - ``test_write_no_precession``: checks that the option to turn off precession is added when asked for

test_make_baseline_date_arrays.py
*******************************************************
Tests the ``rw.make_baseline_date_arrays`` function, which should make
the DATE and BASELINE arrays that are needed to populate a uvfits file. Test
by giving the function a known date string, number of antennas, number of time
steps, and time resolution, and checking the output arrays match expectations.

test_remove_phase_tracking.py
*******************************************************
Tests the ``rw.remove_phase_tracking`` function, which should remove
the phase tracking applied to visibilities. The original MWA correlator
did not phase track, so the ``RTS`` expects no phase tracking on the data, so
to input ``WODEN`` simulations into the ``RTS``, have to undo the phase-tracking.
The ``RTS`` calculates it's own ``u,v,w``, so I only fiddle the visibilities
here so be warned.

This test starts by creating a random array layout via:

.. code-block:: python

  num_antennas = 50
  ##Make a random array layout
  east = np.random.uniform(-1000, 1000, num_antennas)
  north = np.random.uniform(-1000, 1000, num_antennas)
  height = np.random.uniform(0, 10, num_antennas)

These coordinates can then be used the calculate *u,v,w* coodinates for a given
array location (I'm using the MWA site) and phase-centre.

First of all, for 10 frequency channels (100MHz to 190MHz at 10MHz resolution),
and for 10 time steps (at a 2s resolution), calculate the phase-tracked
measurement equation:

.. math::

    V_{\textrm{phased}} = \exp\left[2\pi i \left(ul + vm + w(n-1) \right) \right]

where the :math:`u,v,w` and :math:`l,m,n` are calculated with a phase centre of RA, Dec =
:math:`40^\circ, -50^\circ`, and I calculate a single :math:`l,m,n` for a source at
RA, Dec = :math:`10^\circ, -15^\circ` (so in this setup, :math:`u,v,w` change with
time, and :math:`l,m,n` are constant).

I also calculate the  measurement equation without phase tracking, where I calculate
:math:`u_{\mathrm{zen}},v_{\mathrm{zen}},w_{\mathrm{zen}}` and
:math:`l_{\mathrm{zen}},m_{\mathrm{zen}},n_{\mathrm{zen}}`, using the zenith of
the instrument as a coordinate system centre, and use the following
equation:

.. math::

    V_{\textrm{unphased}} = \exp\left[2\pi i \left(u_{\mathrm{zen}}l_{\mathrm{zen}} + v_{\mathrm{zen}}m_{\mathrm{zen}} + w_{\mathrm{zen}}n_{\mathrm{zen}} \right) \right]

(in this setup, :math:`u_{\mathrm{zen}},v_{\mathrm{zen}},w_{\mathrm{zen}}`
are constant with time, and :math:`l_{\mathrm{zen}},m_{\mathrm{zen}},n_{\mathrm{zen}}`
change with time).

I then use :math:`V_{\textrm{phased}}` as an input to ``rw.remove_phase_tracking``
along with :math:`w`, and use that to unwrap the phase tracking. I then assert
that the output of ``rw.remove_phase_tracking`` matches :math:`V_{\textrm{unphased}}`.

test_argument_inputs.py
*******************************************************
These tests run ``rw.get_parser``, which runs the command line parser, and
``rw.check_args``, which checks the ``args`` collected by ``rw.get_parser``
are parsed correctly. It also does sanity checks on certain combinations of args
such that we don't feed WODEN arguments that won't work. The following tests are run
with the expected outcomes:

 - ``test_parser_fails``: There are three required arguments, ``--ra0``, ``--dec0``, and ``--cat_filename``. Check the parser errors if missing.
 - ``test_missing_args_without_metafits_fails``: If the user doesn't supply the ``--metafits`` arg, there are a minimum set of arguments that must be input. Check ``rw.check_args`` errors if they are missing
 - ``test_metafits_read_fails``: Check ``rw.check_args`` errors if there is a bad path to a metafits file
 - ``test_read_metafits_succeeds``: Check the correct values are read in from a known metafits file
 - ``test_EDA2_args_work``: Check the correct arguments are selected for an EDA2 beam simulation
 - ``test_GaussBeam_args_work``: Check that Gaussian beam related arguments work as expected. Iteratively check that if arguments with defaults are not given (e.g. ``--gauss_ra_point``) that they are set to their defaults, and if they *are* supplied, that they match the given value.
 - ``test_MWAFEEBeam_args_work``: Checks that the MWA FEE primary beam is handled by ``ra.check_args`` correctly. The function should error out if certain paths to the hdf5 file that holds the spherical harmonic information is missing, and if the delays have been specified incorrectly. Check that things work when the correct arguments are given.
 - ``test_MWAFEEBeamInterp_args_work``: Same as ``test_MWAFEEBeam_args_work``, but for the interpolated FEE beam.
 - ``test_MWAAnalyBeam_args_work``: Checks for the analytic MWA beam. Same as ``test_MWAFEEBeam_args_work``, but only checking the delays are set correctly, as no need for an hdf5 file for this model


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
