``fundamental_coords``
=========================
Tests for the functions in ``WODEN/src/fundamental_coords_cpu.c`` and
``WODEN/src/fundamental_coords_gpu.cpp``. These functions calculate ``lmn``
and ``uvw`` coordinates, on the cpu and gpu respectively. The same tests are
run on both the cpu and gpu code, to ensure consistency.

``test_lmn_coords*.c``
*********************************
This either runs ``fundamental_coords_cpu.c::calc_lmn_cpu`` or
``fundamental_coords_gpu.cpp::calc_lmn_for_components_gpu``, which calculate *l,m,n* coords.

We run two control tests, both that generate analytically predictable
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
of outputs. Note that this function is entirely 64-bit whether in FLOAT or
DOUBLE compile mode. The absolute tolerance for outputs vs expectation tabled
above is 1e-15 (this function is a good 'un).

test_uvw_coords*.c
*********************************
This runs either ``fundamental_coords_gpu.cpp::calc_uvw_gpu`` and 
``fundamental_coords_gpu.cpp::calc_uv_shapelet_gpu``, or 
``fundamental_coords_cpu.c::calc_uvw_cpu`` and
``fundamental_coords_cpu.c::calc_uv_shapelet_cpu``.

Both functions calculate *u,v,w* coords in slightly different circumstances.
``calc_uvw`` calculates *u,v,w* coords towards a specified *RA,Dec* phase centre,
for a given set of baseline lengths :math:`X_{\mathrm{diff}}, Y_{\mathrm{diff}}, Z_{\mathrm{diff}}`,
for a number of LSTs and frequencies (meaning *u,v,w* change with time and frequency).

``calc_uv_shapelet`` calculates a number of *u,v* coordinate systems,
each centred on a different SHAPELET component, and does *not* scale by wavelength;
everything remains in metres. This is to save on memory, at the cost of
having to divide by wavelength inside kernels at later times.

Both functions are tested for scaling as changing by time, and
``calc_uvw`` is tested to change with wavelength. To generate
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
of ``calc_uvw_shapelet``, this is checked for each SHAPELET component,
making sure that the outputs are ordered as expected. 
