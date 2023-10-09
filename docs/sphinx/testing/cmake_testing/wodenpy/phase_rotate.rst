``phase_rotate``
=========================
Tests for the functions in ``phase_rotate.phase_tracking``, which simply attempt to undo phase tracking on visibilities.

test_remove_phase_tracking.py
*******************************************************
Tests the ``phase_tracking.remove_phase_tracking`` function, which should remove
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

I then use :math:`V_{\textrm{phased}}` as an input to ``phase_rotate.remove_phase_tracking``
along with :math:`w`, and use that to unwrap the phase tracking. I then assert
that the output of ``phase_rotate.remove_phase_tracking`` matches :math:`V_{\textrm{unphased}}`.

