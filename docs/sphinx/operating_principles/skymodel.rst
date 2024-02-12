.. _defined for hyperdrive: https://github.com/MWATelescope/mwa_hyperdrive/wiki/Source-lists
.. _Line et al. 2020: https://doi.org/10.1017/pasa.2020.18
.. _SHAMFI readthedocs: https://shamfi.readthedocs.io/en/latest/
.. _Callingham et al. 2017: https://iopscience.iop.org/article/10.3847/1538-4357/836/2/174/pdf
.. _Lynch et al. 2021: https://doi.org/10.1017/pasa.2021.50

.. _sky model formats:

Sky model
===========================
.. warning:: In version 2.0.0 (the current version) ``WODEN`` only reads Stokes I from the sky model, and ignores all Stokes QUV. This is because a) this makes things way faster as we can ignore QUV in call calculations b) the power-law and curved power-law extrapolations of QUV were unphysical and c) no one was using it anyway. I'm hoping in a future release to implement using a rotation measure to get linear polarisation going, and to get reading


There are currently three sky model formats that ``WODEN`` accepts; the
preferred ``FITS`` format, the ``hyperdrive yaml`` format, and the now deprecated native ``WODEN`` format.
The ``FITS`` and ``yaml`` formats have greater functionality, as curved power law and
list-style flux behaviours are included, beyond a simple power (but as said above, Stokes QUV are all ignored and only I is used). The native
``WODEN`` format only includes power-law behaviour, and will not be developed
any further.

The ``FITS`` format is by far the most efficient in terms of read speed and manipulation however, and so I _strongly_ suggest you suggest you use that format. Below, we'll define how the spectral information is used by ``WODEN``, and then specify the format for each sky model type.

Spectral models
^^^^^^^^^^^^^^^^^^^^

There are three spectral model types: ``power_law``; ``curved_power_law``; ``list``.
The ``power_law`` model is a simple power law with:

.. math::
  S_i = S_0 \left( \frac{\nu_i}{\nu_0} \right)^\alpha

where :math:`S_i` is the flux density extrapolated to frequency :math:`\nu_i`, with a reference flux density :math:`S_0`, reference frequency :math:`\nu_0`, and spectral index  :math:`\alpha`.
An example of these models (these plots are stolen from unit tests of ``WODEN``):

.. image:: test_extrap_power_laws.png
   :width: 400pt

The ``curved_power_law`` model is defined in Equation 2 of `Callingham et al. 2017`_, and is
implemented in ``WODEN`` as:

.. math::
  S_i = S_0 \left( \frac{\nu_i}{\nu_0} \right)^\alpha e^{q\ln(\frac{\nu_i}{\nu_0})^2}

where :math:`q` is a curvature term. This allows for peaked-type SEDs:

.. image:: test_extrap_curve_power_laws.png
   :width: 400pt

Finally, the ``list`` model just takes a list of flux densities, and does a linear
interpolation between them. This is useful for simulating things like the 21cm
signal which should bounce around with frequency.

.. note::

	Most SEDs are assumed to be close to a power-law, and so are linear in log space. The linear interpolation is therefore done in log space. However, for some complicated modelling, negative flux values are required. For these values, the interpolation is done in *linear* space. Hence in the plot below, where the fluxes dip into negative territory, the extrapolated fluxes not longer lie on the straight lines, as these are log-log plots. A 21cm list-style SED has been tested to return the expected power spectra, even with this slightly odd log/linear extrapolation combo when bouncing around zero.

.. image:: test_extrap_list_laws.png
   :width: 400pt

Sky model formats
^^^^^^^^^^^^^^^^^^^^

``LoBES`` FITS sky model format
----------------------------------
This sky model follows (and expands upon) the format of the LoBES catalogue `Lynch et al. 2021`_ and is the preferred format as it's the fastest and easiest to lazy load. There are three COMPONENT types: point source; Gaussian; shapelets. These are all the model types as defined in `Line et al. 2020`_ (including the mathematics of how each model is simulated). You can create any number of SOURCEs, each with any number of COMPONENTs, by using the `UNQ_SOURCE_ID` and `NAME` columns as detailed below. The sky model is a FITS file with at least one HDU with the following columns:



.. list-table:: FITS HDU 0 columns
   :header-rows: 1
   :widths: 10 10 80
   :stub-columns: 1

   *  - Column Name
      - Unit
      - Description
   *  - UNQ_SOURCE_ID
      -
      - Unique source ID. This is used to group COMPONENTs into SOURCEs. If you want to have multiple components in a single source, they must have the same ``UNQ_SOURCE_ID``.
   *  - NAME
      -
      - This is a COMPONENT name, and should read as UNQ_SOURCE_ID_C`number` where \`number\` is a COMPONENT number. For example, if you have a SOURCE with UNQ_SOURCE_ID = FornaxA, and you have two components, you should have two rows with NAME = FornaxA_C000 and FornaxA_C001.
   *  - RA
      - deg
      - Right Ascension (J2000)
   *  - DEC
      - deg
      - Declination (J2000)
   *  - COMP_TYPE
      -
      - Specifies if the component is a point source, Gaussian, or shapelet. Entered as either ``P``, ``G``, or ``S``.
   *  - MAJOR_DC
      - deg
      - Major axis of Gaussian or shapelet model
   *  - MINOR_DC
      - deg
      - Minor axis of Gaussian or shapelet model
   *  - PA_DC
      - deg
      - Position angle of Gaussian or shapelet model
   *  - MOD_TYPE
      -
      - The flux model of this component. Can be either ``pl`` (power-law), ``cpl`` (curved power-law), or ``nan`` (list of flux densities).
   *  - NORM_COMP_PL
      - Jy
      - The reference flux for a power-law ``pl`` component model, *must be at the reference frequency 200MHz.*
   *  - ALPHA_PL
      -
      - The spectral index for a power-law ``pl`` component model
   *  - NORM_COMP_CPL
      - Jy
      - The reference flux for a curved power-law ``cpl`` component model, *must be at the reference frequency 200MHz.*
   *  - ALPHA_CPL
      -
      - The spectral index for a curved power-law ``cpl`` component model, *must be at the reference frequency 200MHz.*
   *  - CURVE_CPL
      -
      - The curvature `q` for a curved power-law ``cpl`` component model
   *  - INT_FLX*frequency*
      - Jy
      - A reference Stokes I flux density, where *frequency* is the frequency in MHz. For a list type flux model, you can include as many INT_FLX*frequency* columns as necessary. For example, if you have three reference fluxes at 100, 150, and 200 MHz, you should have three columns INT_FLX100, INT_FLX150, and INT_FLX200.

If you want to include SHAPELETS, you must include a second HDU that details the shapelet basis functions for each component, using the following columns:

.. list-table:: FITS HDU 1 columns
   :header-rows: 1
   :widths: 10 80
   :stub-columns: 1

   *  - Column Name
      - Description
   *  - NAME
      - The COMPONENT name exactly as appears in the first HDU ``NAME`` column. You can have multiple rows for each COMPONENT, each with a unique combination of ``N1``, ``N2``, and ``COEFF``, to include as many shapelet basis functions as necessary. ``WODEN`` will cross-reference the two HDUs to use these basis functions in conjunction with the position and flux model in the first HDU.
   *  - N1
      - The first shapelet order
   *  - N2
      - The second shapelet order
   *  - COEFF
      - The coefficient to multiply this basis function by

``hyperdrive`` sky model format
----------------------------------
This is the sky model format as `defined for hyperdrive`_. I'll reproduce
some of the documentation to save clicking on the link, but all credit to
Chris Jordan.

There are three COMPONENT types: point source; Gaussian; shapelet. These are
all the model types as defined in `Line et al. 2020`_ (including the mathematics
of how each model is simulated). You can create any number of SOURCEs, each
with any number of COMPONENTs.

.. note:: ``WODEN`` crops out anything below the horizon for a given observation, with ``run_woden.py`` giving the option to either cropby SOURCE (default) or by COMPONENTs (``--sky_crop_components``). The difference is if any COMPONENT in a single SOURCE is below the horizon, when cropping by SOURCE, the whole SOURCE is thrown away, but when cropping by COMPONENT, only the COMPONENT is thrown away.



Read on for how to detail each model in the ``hyperdrive`` format.

Point sources and flux models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example of a single SOURCE with a single point source COMPONENT is::

  source_name:
  - ra: 60.0
    dec: -27.0
    comp_type: point
    flux_type:
      power_law:
        si: -0.8
        fd:
          freq: 170000000.0
          i: 1.0
          q: 2.0
          u: 3.0
          v: 4.0

An explanation of each line and value follows.

::

  source_name:

Initialises the SOURCE, giving it the name ``source_name``.

::

  - ra: 60.0
    dec: -27.0
    comp_type: point

Initialises a new component, specifying the type (either point, gaussian, shapelet) and the RA and DEC (deg, deg). So this line means a point source at RA,DEC = 4h, -27deg.

::

  power_law:
    si: -0.8

This specifies that this is a ``power_law`` type flux behaviour, with a spectral
index of -0.8.

::

  fd:
    freq: 170000000.0
    i: 1.0
    q: 2.0
    u: 3.0
    v: 4.0

This contains the reference flux density information, with a
reference frequency (Hz) of 170MHz, and reference flux densities
(Jy) of the Stokes *I,Q,U,V* of 1,2,3,4 Jy respectively.

To change to a ``curved_power_law`` flux behaviour, use:

::

  point_curve:
  - ra: 15.0
    dec: -30.0
    comp_type: point
    flux_type:
      curved_power_law:
        si: -0.8
        fd:
          freq: 150000000.0
          i: 1.0
          q: 0.0
          u: 0.0
          v: 0.0
        q: 0.2

Where the extra final line ``q: 0.2`` specifies the curvature term. The indentation
becomes important here, otherwise your Stokes Q value and curvature terms
can get mixed up.

To change to a ``list`` flux behaviour, use:

::

  point_list:
  - ra: 15.0
    dec: -30.0
    comp_type: point
    flux_type:
      list:
        - freq: 180000000.0
          i: 10.0
        - freq: 170000000.0
          i: 5.0
          q: 1.0
          u: 2.0
          v: 3.0
        - freq: 190000000.0
          i: 4.0
          u: 3.0
        - freq: 120000000.0
          i: 1.0
          q: -2.0

Which will collect all the listed Stokes parameters inside each new ``-freq``
entry. This example shows you can have missing parameters; these will be filled
in a zero for you. You can also add the frequencies in any order you want; ``WODEN``
will order them as it reads them in. To be explicit, the following information
is read in from this sky model:

.. list-table::
   :widths: 30 30 30 30 30
   :header-rows: 1

   * - Reference freq (MHz)
     - Stokes I (Jy)
     - Stokes Q (Jy)
     - Stokes U (Jy)
     - Stokes V (Jy)
   * - 120
     - 1
     - -2
     - 0
     - 0
   * - 170
     - 5
     - 1
     - 2
     - 3
   * - 180
     - 10
     - 0
     - 0
     - 0
   * - 190
     - 4
     - 0
     - 3
     - 0

``WODEN`` will then perform 4 separate linear interpolations, one for each
Stokes parameter.

Multiple SOURCEs and COMPONENTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To add multiple SOURCEs, simply repeat the process, e.g.:

::

  source1:
  - ra: 60.0
    dec: -27.0
    comp_type: point
    flux_type:
      power_law:
        si: -0.8
        fd:
          freq: 170000000.0
          i: 1.0
          q: 2.0
          u: 3.0
          v: 4.0
  source2:
  - ra: 12.0
    dec: -35.0
    comp_type: point
    flux_type:
      power_law:
        si: -0.1
        fd:
          freq: 120000000.0
          i: 10.0
          q: 0.0
          u: 0.0
          v: 0.0

To put two COMPONENTs into the same SOURCE, just omit the second name. You can
also add comments without breaking the sky model:

::

  one_source_two_components:
  - ra: 60.0
    dec: -27.0
    comp_type: point
    flux_type:
      power_law:
        si: -0.8
        fd:
          freq: 170000000.0
          i: 1.0
          q: 2.0
          u: 3.0
          v: 4.0
  ##Here is a comment
  - ra: 12.0
    dec: -35.0
    comp_type: point
    flux_type:
      power_law:
        si: -0.1
        fd:
          freq: 120000000.0
          i: 10.0
          q: 0.0
          u: 0.0
          v: 0.0

Gaussian sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example srclist containing a single gaussian::

  singlegauss_power:
  - ra: 30.0
    dec: -30.0
    comp_type:
      gaussian:
        maj: 180.
        min: 360.
        pa: -10.
    flux_type:
      power_law:
        si: -0.8
        fd:
          freq: 150000000.0
          i: 2.0
          q: 0.0
          u: 0.0
          v: 0.0

where all lines have the same meaning as the point source, aside from the lines::

  comp_type:
    gaussian:
      maj: 180.
      min: 360.
      pa: -10.

denote Gaussian specific parameters. The FWHM major ``maj`` and minor ``min`` axes
are given in arcseconds, with the position angle (East from North) given in degrees.

Shapelet sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To generate shapelet models compatible with ``WODEN``, use ``SHAMFI`` to fit an
image with the ``--woden_srclist`` option (again see `SHAMFI readthedocs`_.
for more detail). This will ensure all normalisations are correct.

.. warning:: At the time of writing, ``SHAMFI`` spits out either ``RTS`` or ``WODEN`` style sky models. You'll need to use ``hyperdrive`` to convert the outputs into a ``hyperdrive`` style sky model. I have just created an issue on the ``SHAMFI`` github so hopefully it'll get done soonish.

An example sky model (made by hand so the normalisations *won't* be correct) is::

  singleshapelet_power:
  - ra: 45.0
    dec: 20.0
    comp_type:
      shapelet:
        maj: 420.
        min: 300.
        pa: 56.
        coeffs:
          - n1: 0
            n2: 0
            value: 0.48255952
          - n1: 14
            n2: 2
            value: -0.18494293
          - n1: 41
            n2: -15
            value: -0.08973978
          - n1: 37
            n2: 7
            value: -0.22137849
    flux_type:
      power_law:
        si: -0.8
        fd:
          freq: 150000000.0
          i: 10.0
          q: 0.0
          u: 0.0
          v: 0.0

This sky model will generate a single SHAPELET component, with a ``power_law``
type flux behaviour. This SHAPELET model requires 4 basis functions, each detailed
as::

  coeffs:
    - n1: 0
      n2: 0
      value: 0.48255952
    - n1: 14
      n2: 2
      value: -0.18494293
    - n1: 41
      n2: -15
      value: -0.08973978
    - n1: 37
      n2: 7
      value: -0.22137849

where ``n1, n2`` detail the order of the basis function, and ``value`` gives the coefficient
to multiply the basis function by (see `Line et al. 2020`_ for details). You can
write as many ``n1,n2`` paris as necessary, with a maximum order < 100. If you use ``SHAMFI``,
the coefficients will be scaled such that the integrated Stokes I flux density over the
source will be 10 Jy at 150 MHz for this example.

Putting it all together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example skymodel with four sources, the first with all component types, the next three with a single component of each type,  would look something like this::

  point_list_gauss_power_shape_curve:
  - ra: 15.0
    dec: -30.0
    comp_type: point
    flux_type:
      list:
        - freq: 180000000.0
          i: 10.0
        - freq: 170000000.0
          i: 5.0
          q: 1.0
          u: 2.0
          v: 3.0
        - freq: 190000000.0
          i: 4.0
          u: 3.0
        - freq: 120000000.0
          i: 1.0
          q: -2.0
  - ra: 30.0
    dec: -30.0
    comp_type:
      gaussian:
        maj: 180.
        min: 360.
        pa: -10.
    flux_type:
      power_law:
        si: -0.8
        fd:
          freq: 150000000.0
          i: 2.0
          q: 0.0
          u: 0.0
          v: 0.0
  - ra: 15.0
    dec: -30.0
    comp_type:
      shapelet:
        maj: 420.
        min: 300.
        pa: 56.
        coeffs:
          - n1: 0
            n2: 0
            value: 0.48255952
          - n1: 14
            n2: 2
            value: -0.18494293
          - n1: 41
            n2: -15
            value: -0.08973978
          - n1: 37
            n2: 7
            value: -0.22137849
    flux_type:
      curved_power_law:
        si: -0.8
        fd:
          freq: 150000000.0
          i: 1.0
          q: 0.0
          u: 0.0
          v: 0.0
        q: 0.2
  point_list_alone:
  - ra: 15.0
    dec: -30.0
    comp_type: point
    flux_type:
      list:
        - freq: 180000000.0
          i: 10.0
        - freq: 170000000.0
          i: 5.0
          q: 1.0
          u: 2.0
          v: 3.0
        - freq: 190000000.0
          i: 4.0
          u: 3.0
        - freq: 120000000.0
          i: 1.0
          q: -2.0
  gauss_power_alone:
  - ra: 30.0
    dec: -30.0
    comp_type:
      gaussian:
        maj: 180.
        min: 360.
        pa: -10.
    flux_type:
      power_law:
        si: -0.8
        fd:
          freq: 150000000.0
          i: 2.0
          q: 0.0
          u: 0.0
          v: 0.0
  shapelet_curve_alone:
  - ra: 15.0
    dec: -30.0
    comp_type:
      shapelet:
        maj: 420.
        min: 300.
        pa: 56.
        coeffs:
          - n1: 0
            n2: 0
            value: 0.48255952
          - n1: 14
            n2: 2
            value: -0.18494293
          - n1: 41
            n2: -15
            value: -0.08973978
          - n1: 37
            n2: 7
            value: -0.22137849
    flux_type:
      curved_power_law:
        si: -0.8
        fd:
          freq: 150000000.0
          i: 1.0
          q: 0.0
          u: 0.0
          v: 0.0
        q: 0.2

``WODEN`` sky model format
-------------------------------

.. note Really, just use the ``hyperdrive`` format from now on.

The ``WODEN`` source catalogue is a modified version of the ``RTS`` srclist. In the current version of ``WODEN``, you create one single SOURCE which can include as many COMPONENTS as desired, each of type ``POINT``, ``GAUSSIAN`` or ``SHAPELET``. A ``POINT`` is a dirac delta point source model, a GAUSSIAN is a 2D Gaussian model (with a major, minor, and position angle), and a ``SHAPELET`` model uses multiple 'shapelet' basis functions to build a model. For details on the model types, see `Line et al. 2020`_. If you want to build a shapelet model, you can use the software ``SHAMFI``, which you can read about on the `SHAMFI readthedocs`_.

Currently, every source is given a simple power-law frequency behaviour as:

.. math::
  S = S_0 \left( \frac{\nu_0}{\nu} \right)^\alpha

where :math:`S` is the flux density at frequency :math:`\nu`, with a reference flux density :math:`S_0`, reference frequency :math:`\nu_0`, and spectral index  :math:`\alpha`.

Point sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example of a single SOURCE with a single point source COMPONENT is::

  SOURCE source_name P 1 G 0 S 0 0
  COMPONENT POINT 4.0 -27.0
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  ENDCOMPONENT
  ENDSOURCE

An explanation of each line and value follows.

::

  SOURCE source_name P 1 G 0 S 0 0

Initialises the SOURCE, giving it the name ``source_name``, and specifying the number and type of components (P = point, G = gaussian, S = shapelet). For shapelet, the two numbers are total number of coefficients and total number of components. Read on further for more explanation of shapelets.

::

  COMPONENT POINT 4.0 -27.0

Initialises a component, specifying the type (either POINT, GAUSSIAN, SHAPELET) and the RA and DEC (hours, deg). So this line means a point source at RA,DEC = 4h, -27deg.

::

  LINEAR 1.8e+08 10.0 0 0 0 -0.8

Specifies a reference Stokes flux density as *LINEAR Freq I Q U V SI*, where the Freq is in Hz, Stokes params *I,Q,U,V* are all in units of Jy, and SI is the spectral index. It's labelled ``LINEAR`` as a power-law is linear in log-log space. This example line specifies we have a source that has a flux density of purely Stokes I of 10 Jy at 180 MHz, with a spectral index if -0.8.

::

  ENDCOMPONENT

This line ends the component.

::

  ENDSOURCE

This line ends the source.


To add multiple point sources, simply repeat the ``COMPONENT`` / ``ENDCOMPONENT`` sections with new details, i.e.

::

  SOURCE multi_point P 3 G 0 S 0 0
  COMPONENT POINT 4.0 -27.0
  LINEAR 1.8e+08 10.0 0 0 0 -0.4
  ENDCOMPONENT
  COMPONENT POINT 3.0 -37.0
  LINEAR 1.3e+08 1.0.0 0 0 0 -0.786
  ENDCOMPONENT
  COMPONENT POINT 5.0 -47.0
  LINEAR 3.9e+08 0.04 0 0 0 .02
  ENDCOMPONENT
  ENDSOURCE

noting that at the very top line, I have updated ``P 3`` to reflect there are now three point sources. These numbers are used to quickly allocate memory, that's why they re included.

.. note:: ``WODEN`` crops everything below the horizon out of the sky model. It can do this one of two ways - either by ``COMPONENT`` or by ``SOURCE``. In the example above, we have three COMPONENT in one SOURCE. If you ask ``WODEN`` to crop by ``SOURCE``, if just one of the ``COMPONENTS`` is below the horizon, it'll crop the *entire* source.

Gaussian sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example srclist containing a single gaussian::

  SOURCE gaussian_source P 0 G 1 S 0 0
  COMPONENT GAUSSIAN 3.378 -37.2
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  GPARAMS 45.0000000000 6.0 3.0
  ENDCOMPONENT
  ENDSOURCE

where all lines have the same meaning as the point source, and the meaning of the extra line::

  GPARAMS 45.0000000000 6.0 3.0

which specifies the Gaussian parameters as ``GPARAMS pa(deg) major_axis(arcmin) minor_axis(arcmin)``. The major and minor axes are specified as FWHM. Note this line needs to sit in between the lines starting with ```COMPONENT GAUSSIAN`` and ```ENDCOMPONENT``.

Shapelet sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To generate shapelet models compatible with WODEN, simply use ``SHAMFI`` to fit an image with the ``--woden_srclist`` option (again see `SHAMFI readthedocs`_. for more detail). This will ensure all normalisations are correct. An example sky model (made by hand so the normalisations *won't* be correct) is::

  SOURCE shapelet_source P 0 G 0 S 1 3
  COMPONENT SHAPELET 3.378 -37.2
  FREQ 1.8e+08 10.0 0 0 0
  SPARAMS 45.0000000000 6.0 3.0
  SCOEFF 0 0 0.92342
  SCOEFF 1 10 0.0002354
  SCOEFF 4 5 0.004567
  ENDCOMPONENT
  ENDSOURCE

which generates a single shapelet component, including 3 shapelet basis functions, hence ``S 1 3`` in the first line. The ``SPARAMS`` line is similar to the ``GAUSSIAN`` line with ``SPARAMS pa(deg) major_axis(arcmin) minor_axis(arcmin)``. The extra lines like::

  SCOEFF 0 0 0.92342

encode the order of the shapelet basis function (see `Line et al. 2020`_ for details) and fitted coefficient as ``SCOEFF p1 p2 coeff_value``. You can add as many ``SCOEFF`` lines as necessary, with a maximum order < 100. If you use ``SHAMFI``, the coefficients will be scaled such that the Stokes I flux density of the full source will be 10 Jy at 180 MHz for this example. You may have noticed the SED information here is different::

  FREQ 1.8e+08 10.0 0 0 0

This line will still assume a power-law frequency behaviour, with a reference flux of 10 Jy at 180 MHz, but use a default SI = -0.8.

Putting it all together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example skymodel with four sources, the first with all component types, the next three with a single component of each type,  would look something like this::

  SOURCE multi_sources P 3 G 1 S 2 7
  COMPONENT SHAPELET 3.378 -37.2
  FREQ 1.8e+08 10.0 0 0 0
  SPARAMS 45.0000000000 6.0 3.0
  SCOEFF 0 0 0.92342
  SCOEFF 1 10 0.0002354
  SCOEFF 4 5 0.004567
  ENDCOMPONENT
  COMPONENT SHAPELET 3.12 -32.2
  FREQ 1.8e+08 3.1 0 0 0
  SPARAMS 56.0000000000 9.0 3.0
  SCOEFF 0 0 0.02345
  SCOEFF 3 0 -0.234234
  SCOEFF 21 34 0.82342
  SCOEFF 31 5 -0.00876234
  ENDCOMPONENT
  COMPONENT GAUSSIAN 3.378 -37.2
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  GPARAMS 45.0000000000 6.0 3.0
  ENDCOMPONENT
  COMPONENT POINT 4.0 -27.0
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  ENDCOMPONENT
  COMPONENT POINT 3.0 -37.0
  LINEAR 1.8e+08 0.6 0 0.2 0 -0.8
  ENDCOMPONENT
  COMPONENT POINT 5.0 -47.0
  LINEAR 70E+6 87.0 0 0 0 -0.8
  ENDCOMPONENT
  ENDSOURCE
  SOURCE source_name P 1 G 0 S 0 0
  COMPONENT POINT 4.0 -27.0
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  ENDCOMPONENT
  ENDSOURCE
  SOURCE gaussian_source P 0 G 1 S 0 0
  COMPONENT GAUSSIAN 3.378 -37.2
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  GPARAMS 45.0000000000 6.0 3.0
  ENDCOMPONENT
  ENDSOURCE
  SOURCE shapelet_source P 0 G 0 S 1 3
  COMPONENT SHAPELET 3.378 -37.2
  LINEAR 1.1e+08 10.0 2.0 0 0.8 -0.7
  SPARAMS 45.0000000000 6.0 3.0
  SCOEFF 0 0 0.92342
  SCOEFF 1 10 0.0002354
  SCOEFF 4 5 0.004567
  ENDCOMPONENT
  ENDSOURCE
