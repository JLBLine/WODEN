``add_instrumental_effects_woden.py``
======================================

Use this script to add instrumental effects to a WODEN ``uvfits`` output. The script is designed to be run on a single ``uvfits``, so I suggest running ``add_woden_uvfits.py`` if you have 24 separate ``uvfits`` files. Read on for details on how different types of instrumental effects are added.

Adding noise
-------------
.. warning:: Everything is set for MWA highband simulations currently. This script needs developing. The noise is only correctly added for cross-correlations. I don't know how to add noise for auto-correlations yet, so they are left alone.

The noise on cross-correlations is drawn for a normal distribution with a mean of zero and a standard deviation of $\sigma$, which is calculated via the following (see Equation 6.50 in TMS 3rd edition for details):

$$
\sigma = \frac{\sqrt{2}k_b(T_{sky} + T_{rec})}{A_{eff}\sqrt{\Delta\nu\Delta t}},
$$

where variables and there defaults are listed below

.. list-table:: 
   :widths: 25 25 50
   :header-rows: 1

   * - Variable
     - Default
     - Description
   * - $T_{sky}$
     - $228(\frac{\nu}{150\mathrm{MHz}})^{-2.53}$
     - The sky temperature in Kelvin; estimated by extrapolating values from 150 MHz
   * - $T_{rec}$
     - 50
     - The receiver temperature in Kelvin (estimate for MWA highband)
   * - $A_{eff}$
     - 20.35
     - The effective area of the antenna in square metres (estimate for MWA highband)
   * - $\Delta\nu$
     - Value in ``uvfits`` header
     - The frequency resolution in Hz
   * - $\Delta t$
     - Inferred from ``uvfits``
     - The time integration in seconds

The noise is drawn randomly and added separately for the real and imaginary parts of each visibility, and only added to the XX and YY polarisations. If noise is to be added, it is done before adding any antenna gains or leakages, which should naturally include any noise from XX,YY into XY,YX.

Adding antenna (tile) gains and leakages
-----------------------------------------
Each receiving element (often called a tile for MWA, called an antenna hereon) in a dual polarisation interferometer can have a gain ($g_x, g_y$) and leakage ($D_x,D_y$) term for each polarisation. Each visibility is a correlation of two antennas. ``WODEN`` implements any antenna gains and leakages by multiplying the visibility between tiles 1 and 2 through the following operation:

.. math::
   \begin{bmatrix}
   V^`_{12\,XX} & V^`_{12\,XY} \\
   V^`_{12\,YX} & V^`_{12\,YY}
   \end{bmatrix} =
   \begin{bmatrix}
   g_{x1} & D_{x1} \\
   D_{y1} & g_{y1}
   \end{bmatrix}
   \begin{bmatrix}
   V_{12\,XX} & V_{12\,XY} \\
   V_{12\,YX} & V_{12\,YY}
   \end{bmatrix}
   \begin{bmatrix}
   g_{x2}^{\ast} & D_{x2}^{\ast} \\
   D_{y2}^{\ast} & g_{y2}^{\ast}
   \end{bmatrix},

where $V$ is a visibility in the ``uvfits`` file, $\ast$ means the complex conjugate, and $V^`$ is the visibility after the antenna gains and leakages have been applied. In ``WODEN`` 2.0.0, you can add a random frequency independent amplitude and a phase dependent slope to each antenna gain (see ``--ant_gain_amp_error`` and ``--ant_gain_phase_error``below for more detail).  The leakage terms are calculated via Equation A4.5 from TMS third edition where

$
\begin{align}
D_x = \Psi - j \chi \\
D_y = -\Psi + j \chi
\end{align}
$

where $\Psi, \chi$ are alignment errors of the dipoles. This equation is really designed for single antennas, but in the MWA case, you could imagine that all the dipoles in a tile are aligned aligned perfectly to the mesh, and the mesh is slightly offset, so the alignment errors for all dipoles are the same. The parameters. $\Psi, \chi$ are set by the user via ``--ant_leak_errs`` so you can tune however much leakage you want.

*Command line running options*
-------------------------------

.. argparse::
   :filename: ../../scripts/add_instrumental_effects_woden.py
   :func: get_parser
   :prog: add_instrumental_effects_woden.py

*Function documentation*
------------------------

.. automodule:: add_instrumental_effects_woden
   :members:
