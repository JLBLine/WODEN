``WODEN`` sky model format
============================

.. _Line et al. 2020: https://doi.org/10.1017/pasa.2020.18
.. _SHAMFI readthedocs: https://shamfi.readthedocs.io/en/latest/


The ``WODEN`` source catalogue is a modified version of the ``RTS`` srclist. In the current version of ``WODEN``, you create one single SOURCE which can include as many COMPONENTS as desired, each of type ``POINT``, ``GAUSSIAN`` or ``SHAPELET``. A ``POINT`` is a dirac delta point source model, a GAUSSIAN is a 2D Gaussian model (with a major, minor, and position angle), and a ``SHAPELET`` model uses multiple 'shapelet' basis functions to build a model. For details on the model types, see `Line et al. 2020`_. If you want to build a shapelet model, you can use the software ``SHAMFI``, which you can read about on the `SHAMFI readthedocs`_.

Currently, every source is given a simple power-law frequency behaviour as:

.. math::
  S = S_0 \left( \frac{\nu_0}{\nu} \right)^\alpha

where :math:`S` is the flux density at frequency :math:`\nu`, with a reference flux density :math:`S_0`, reference frequency :math:`\nu_0`, and spectral index  :math:`\alpha`.

Point sources
^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^

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
