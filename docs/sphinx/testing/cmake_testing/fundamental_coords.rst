``fundamental_coords``
=========================
Tests for the functions in ``WODEN/src/fundamental_coords.cu``.

test_lmn_coords.c
*********************************
This runs ``fundamental_coords::test_kern_calc_lmn``, which tests
``fundamental_coords::kern_calc_lmn``, which calculates *l,m,n* coords.

This runs two control tests, both that generate analytically predictable
outcomes. Both set the phase centre to *RA*:math:`_{\textrm{phase}}`, *Dec*:math:`_{\textrm{phase}}` = :math:`0^\circ, 0^\circ`. One
test holds *Dec* = :math:`0^\circ`, and varies *RA*, the other holds
*RA* = :math:`0^\circ`, and varies *Dec*.  Under these settings the following
should be true:

.. list-table:: Outcomes when *Dec* = :math:`0^\circ`
   :widths: 25 25 25 25
   :header-rows: 1

   * - *RA*
     - *l*
     - *m*
     - *n*
   * - :math:`\frac{3\pi}{2}`
     - :math:`-1`
     - :math:`0`
     - :math:`0`
   * - :math:`\frac{5\pi}{3}`
     - :math:`-\frac{\sqrt{3}}{2}`
     - :math:`0`
     - :math:`0.5`
   * - :math:`\frac{7\pi}{4}`
     - :math:`-\frac{\sqrt{2}}{2}`
     - :math:`0`
     - :math:`\frac{\sqrt{2}}{2}`
   * - :math:`\frac{11\pi}{6}`
     - :math:`-0.5`
     - :math:`0`
     - :math:`\frac{\sqrt{3}}{2}`
   * - :math:`0`
     - :math:`0`
     - :math:`0`
     - :math:`1`
   * - :math:`\frac{\pi}{6}`
     - :math:`0.5`
     - :math:`0`
     - :math:`\frac{\sqrt{3}}{2}`
   * - :math:`\frac{\pi}{4}`
     - :math:`\frac{\sqrt{2}}{2}`
     - :math:`0`
     - :math:`\frac{\sqrt{2}}{2}`
   * - :math:`\frac{\pi}{3}`
     - :math:`\frac{\sqrt{3}}{2}`
     - :math:`0`
     - :math:`0.5`
   * - :math:`\frac{\pi}{2}`
     - :math:`1.0`
     - :math:`0`
     - :math:`0`

.. list-table:: Outcomes when *RA* = :math:`0^\circ`
   :widths: 25 25 25 25
   :header-rows: 1

   * - *Dec*
     - *l*
     - *m*
     - *n*
   * - :math:`-\frac{\pi}{2}`
     - :math:`0`
     - :math:`-1`
     - :math:`0`
   * - :math:`-\frac{\pi}{3}`
     - :math:`0`
     - :math:`-\frac{\sqrt{3}}{4}`
     - :math:`0.5`
   * - :math:`-\frac{\pi}{4}`
     - :math:`0`
     - :math:`-\frac{\sqrt{2}}{2}`
     - :math:`\frac{\sqrt{2}}{2}`
   * - :math:`-\frac{\pi}{6}`
     - :math:`0`
     - :math:`-0.5`
     - :math:`\frac{\sqrt{3}}{2}`
   * - :math:`0`
     - :math:`0`
     - :math:`0`
     - :math:`1`
   * - :math:`\frac{\pi}{6}`
     - :math:`0`
     - :math:`0.5`
     - :math:`\frac{\sqrt{3}}{2}`
   * - :math:`\frac{\pi}{4}`
     - :math:`0`
     - :math:`\frac{\sqrt{2}}{2}`
     - :math:`\frac{\sqrt{2}}{2}`
   * - :math:`\frac{\pi}{3}`
     - :math:`0`
     - :math:`\frac{\sqrt{3}}{2}`
     - :math:`0.5`
   * - :math:`\frac{\pi}{2}`
     - :math:`0`
     - :math:`1.0`
     - :math:`0`

The tests ensure that the inputs yield the outputs as expected, and in
the process check that the execution of the kernel yields the correct number
of outputs. The tests require that the *n* coords be within an absolute tolerance
of 2e-7 to pass. This is a result of the accuracy afforded by the ``sinf``,
``cosf`` float cosine functions.

test_uvw_coords.c
*********************************
This runs both ``fundamental_coords::test_kern_calc_uvw`` as well as
``fundamental_coords::test_kern_calc_uvw_shapelet``, which in turn test
``fundamental_coords::kern_calc_uvw`,
``fundamental_coords::kern_calc_uvw_shapelet`` respectively.

Both kernels calculate *u,v,w* coords (in wavelenths) in slightly different circumstances.
``kern_calc_uvw`` calculates *u,v,w* coords towards a specified *RA,Dec* phase centre,
for a given set of baseline lengths :math:`X_{\mathrm{diff}}, Y_{\mathrm{diff}}, Z_{\mathrm{diff}}`, for a number of LSTs and frequencies (meaning *u,v,w*
change with time and frequency).
``kern_calc_uvw_shapelet`` does the above for a number of *u,v,w* coordinate systems,
each centred on a different SHAPELET component.

Both kernels are tested for scaling by wavelength, and changing by time. To generate
analytically predictable outcomes, the phase centre is again set to
*RA*:math:`_{\textrm{phase}}`, *Dec*:math:`_{\textrm{phase}}` = :math:`0^\circ, 0^\circ`.
Under these conditions, the following is true:

.. math::

   \begin{eqnarray}
   u & = & \left[\sin(H_{\textrm{phase}}) X_{\mathrm{diff}} + \cos(H_{\textrm{phase}}) Y_{\mathrm{diff}} \right] / \lambda \\
   v & = & Z_{\mathrm{diff}} / \lambda \\
   w & = & \left[\cos(H_{\textrm{phase}}) X_{\mathrm{diff}} - \sin(H_{\textrm{phase}}) Y_{\mathrm{diff}} \right] / \lambda
   \end{eqnarray}

where :math:`H_{\textrm{phase}}` is the hour angle of the phase centre. These tests
check that this holds true over multiple time and frequency steps. In the case
of ``kern_calc_uvw_shapelet``, this is checked for each SHAPELET component,
making sure that the outputs are ordered as expected.
