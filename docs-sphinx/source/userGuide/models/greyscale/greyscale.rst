======================================================================================
Single-Phase Multi-Scale Porous Media (Greyscale Model - lbpm_greyscale_simulator)
======================================================================================

The LBPM model to simulate single-phase multi-scale porous media, also called greyscale model, 
is constructed to approximate the solution of the Darcy-Brinkman equations in grey regions, 
coupled to a Navier-Stokes solution in open regions. The generalized model for isothermal fluid 
flow in porous media employs volume-averaged Navier-Stokes equations that account for both linear 
and nonlinear drag effects:

.. math::
   :nowrap:

   $$
   \nabla \cdot \boldsymbol{u} = 0 \;,
   $$

.. math::
   :nowrap:

   $$
   \frac{\partial \boldsymbol{u}}{\partial t} + (\boldsymbol{u} \cdot \nabla) \left(\frac{\boldsymbol{u}}{\epsilon}\right) = -\frac{1}{\rho} \nabla(\epsilon p) + \nu_eff \nabla^2 \boldsymbol{u}  -\frac{\epsilon \nu}{K} \boldsymbol{u} - \frac{\epsilon F_e}{\sqrt{K}} |\boldsymbol{u}| \boldsymbol{u} + \epsilon \boldsymbol{g} \;.
   $$

In this formulation, :math:`p` and :math:`\boldsymbol{u}` denote the volume-averaged pressure and velocity, 
respectively, and :math:`\rho` represents the fluid density. The model incorporates :math:`\nu` as the 
kinematic viscosity of the fluid and :math:`\nu_eff` as an effective viscosity parameter. External force 
effects are captured through the acceleration term :math:`\boldsymbol{g}`. The geometric function :math:`F_e`, 
derived from Ergun's experimental investigations, characterizes the porous medium structure through the relation:

.. math::
   :label: Forchheimer-eq-greysingle

   F_e = \frac{1.75}{\sqrt{150\epsilon^3}}.

****************************
LBM Formulation
****************************

.. raw:: html

   <details>
   <summary><strong>Click here to open: LBM Formulation</strong></summary>

A multiple-relaxation-time (MRT) lattice Boltzmann equation with D3Q19 discretization describes momentum transport for the phase-averaged equations: 

.. math::
   :nowrap:

   $$
   f_q(\boldsymbol{x}_i + \boldsymbol{\xi}_q \delta_t,t + \delta_t) - f_q(\boldsymbol{x}_i,t) = \sum^{Q-1}_{k=0} M^{-1}_{qk} S_{kk} (m_k^{eq}-m_k)  + \sum^{Q-1}_{k=0} M^{-1}_{qk} \left(1-\frac{S_{kk}}{2}\right) \hat{F}_q \;.
   $$

The total force :math:`\boldsymbol{F}` is imposed based on the discrete force construction developed by Guo et al. (2002):

.. math::
   :nowrap:

   $$
   \boldsymbol{F} = - \frac{\epsilon \nu \rho}{K} \boldsymbol{u} - \frac{\epsilon F_{e}}{\sqrt{K}} \rho |\boldsymbol{u}|\boldsymbol{u} +  \epsilon \rho \boldsymbol{g} \;.
   $$

The fluid velocity is computed through the relation for the forcing scheme:

.. math::
   :nowrap:

   $$
   \rho \boldsymbol{u} = \sum_i \boldsymbol{\xi}_i f_i + \frac{\delta_t}{2}  \boldsymbol{F} \;.
   $$

However, to handle the nonlinearity arising from the velocity-dependent force term, an explicit solution for the velocity is employed:

.. math::
   :nowrap:

   $$
   \boldsymbol{u} = \frac{\boldsymbol{v}}{c_0 + \sqrt{c_0^2 + c_1 |\boldsymbol{v}|}} \;,
   $$

that utilizes a temporal velocity field :math:`\boldsymbol{v}`:

.. math::
   :nowrap:

   $$
   \rho \boldsymbol{v} = \sum_i \boldsymbol{\xi}_i f_i + \frac{\delta_t}{2} \rho \epsilon \boldsymbol{g} \;,
   $$

and the coefficients governing this nonlinear transformation are:

.. math::
   :nowrap:

   $$
   c_0 = \frac{1}{2} \left(  1 + \frac{\delta_t}{2}\frac{\epsilon \nu}{K} \right), \qquad c_1 = \frac{\delta_t}{2}\frac{\epsilon F_e}{\sqrt{K}} \;.
   $$

where :math:`c_s^2 = 1/3` is the speed of sound and :math:`\boldsymbol{\xi}_q` are the unit lattice vectors for the D3Q19 lattice discretization:

.. math::
   :nowrap:

   $$
   \boldsymbol{\xi}_q=\xi_{q,\alpha} =\left\{
   \begin{array}{llll}
      \alpha=0, & \{0,0,0\}, & \text{for } q = 0, & w_{0}=1/3 \\
      \alpha=x, & \{\pm 1,0,0\}, & \text{for } q = 1,2, & w_{1,2}=1/18 \\
      \alpha=y, & \{0,\pm 1,0\}, & \text{for } q = 3,4, & w_{3,4}=1/18 \\
      \alpha=z, & \{0,0,\pm 1\}, & \text{for } q = 5,6, & w_{5,6}=1/18 \\
      \alpha=x,y,& \{\pm 1,\pm 1,0\}, & \text{for } q = 7,8,9,10, & w_{7,8,9,10}=1/36 \\
      \alpha=x,z,& \{\pm 1,0,\pm 1\}, & \text{for } q = 11,12,13,14, & w_{11,12,13,14}=1/36 \\
      \alpha=y,z,& \{0,\pm 1,\pm 1\}, & \text{for } q = 15,16,17,18, & w_{15,16,17,18}=1/36 .
   \end{array}\right.
   $$

The moments defined in the Gram-Schmidt basis are linearly independent functions of the distribution functions:

.. math::
   :nowrap:
   
   $$
   m_k = \sum_{q=0}^{18} M_{qk} f_q\; \quad \rightarrow \quad M_{qk} = \left(\tiny
   \begin{array}{rrrrrrrrrrrrrrrrrrr}
   1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\
   -30 & -11 & -11 & -11 & -11 & -11 & -11 & 8 & 8 & 8 & 8 & 8 & 8 & 8 & 8 & 8 & 8 & 8 & 8 \\
   12 & -4 & -4 & -4 & -4 & -4 & -4 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\
   0 & 1 & -1 & 0 & 0 & 0 & 0 & 1 & -1 & 1 & -1 & 0 & 0 & 1 & -1 & 1 & -1 & 0 & 0 \\
   0 & -4 & 4 & 0 & 0 & 0 & 0 & 1 & -1 & 1 & -1 & 0 & 0 & 1 & -1 & 1 & -1 & 0 & 0 \\
   0 & 0 & 0 & 1 & -1 & 0 & 0 & 1 & -1 & 0 & 0 & 1 & -1 & -1 & 1 & 0 & 0 & 1 & -1 \\
   0 & 0 & 0 & -4 & 4 & 0 & 0 & 1 & -1 & 0 & 0 & 1 & -1 & -1 & 1 & 0 & 0 & 1 & -1 \\
   0 & 0 & 0 & 0 & 0 & 1 & -1 & 0 & 0 & 1 & -1 & 1 & -1 & 0 & 0 & -1 & 1 & -1 & 1 \\
   0 & 0 & 0 & 0 & 0 & -4 & 4 & 0 & 0 & 1 & -1 & 1 & -1 & 0 & 0 & -1 & 1 & -1 & 1 \\
   0 & 2 & 2 & -1 & -1 & -1 & -1 & 1 & 1 & 1 & 1 & -2 & -2 & 1 & 1 & 1 & 1 & -2 & -2 \\
   0 & -4 & -4 & 2 & 2 & 2 & 2 & 1 & 1 & 1 & 1 & -2 & -2 & 1 & 1 & 1 & 1 & -2 & -2 \\
   0 & 0 & 0 & 1 & 1 & -1 & -1 & 1 & 1 & -1 & -1 & 0 & 0 & 1 & 1 & -1 & -1 & 0 & 0 \\
   0 & 0 & 0 & -2 & -2 & 2 & 2 & 1 & 1 & -1 & -1 & 0 & 0 & 1 & 1 & -1 & -1 & 0 & 0 \\
   0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & -1 & -1 & 0 & 0 & 0 & 0 \\
   0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & -1 & -1 \\
   0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & -1 & -1 & 0 & 0 \\
   0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & -1 & -1 & 1 & 0 & 0 & 1 & -1 & -1 & 1 & 0 & 0 \\
   0 & 0 & 0 & 0 & 0 & 0 & 0 & -1 & 1 & 0 & 0 & 1 & -1 & 1 & -1 & 0 & 0 & 1 & -1 \\
   0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & -1 & -1 & 1 & 0 & 0 & -1 & 1 & 1 & -1
   \end{array}
   \right)
   $$

The equilibrium formulation is adapted from the generalized volume-averaged macroscopic equations. As 
established in the appendix of Guo and Zhao (2002), the definition of the equilibrium distribution 
function yields the following moments, where the presence of the porous matrix modifies the second-order 
velocity terms by a factor of :math:`1/\epsilon`:

.. math::
   :nowrap:

   $$
   f_{q}^{(eq)} = w_{q}\rho \left( 1 + \frac{\xi_{q,\alpha} u_{\alpha}}{c_{s}^{2}} + \frac{u_{\alpha}u_{\beta}\left(\xi_{q,\alpha}\xi_{q,\beta} - c_{s}^{2}\delta_{\alpha\beta}\right)}{2\epsilon c_{s}^{2}} \right) ,
   $$

Applying this to the D3Q19 MRT basis yields the corresponding equilibrium moments, given by:

.. math::

   m_0^{eq} = \rho \;,

.. math::

   m_1^{eq} = 19 \frac{\rho\left(u_{x}^2+u_{y}^2+u_{z}^2\right)}{\epsilon} - 11\rho \;, \qquad m_2^{eq} = 3\rho - \frac{11}{2} \frac{\rho\left(u_{x}^2+u_{y}^2+u_{z}^2\right)}{\epsilon} \;,

.. math::

   m_3^{eq} = \rho u_x \;, \qquad m_5^{eq} = \rho u_y \;, \qquad m_7^{eq} = \rho u_z \;,

.. math::

   m_4^{eq} = -\frac{2}{3} \rho u_x \;, \qquad m_6^{eq} = -\frac{2}{3} \rho u_y \;, \qquad m_8^{eq} = -\frac{2}{3} \rho u_z \;,

.. math::

   m_9^{eq} = \frac{\rho\left(2u_x^2-u_y^2-u_z^2\right)}{\epsilon} \;, \qquad m_{10}^{eq} = -\frac{\rho\left(2u_x^2-u_y^2-u_z^2\right)}{2\epsilon} \;,

.. math::

   m_{11}^{eq} = \frac{\rho\left(u_y^2-u_z^2\right)}{\epsilon} \;, \qquad m_{12}^{eq} = -\frac{\rho\left(u_y^2-u_z^2\right)}{2 \epsilon} \;,

.. math::

   m_{13}^{eq} = \frac{\rho(u_x u_y)}{\epsilon} \;, \qquad m_{14}^{eq} = \frac{\rho(u_y u_z)}{\epsilon} \;, \qquad m_{15}^{eq} = \frac{\rho(u_x u_z)}{\epsilon} \;,

.. math::

   m_{16}^{eq} = 0 \;, \qquad m_{17}^{eq} = 0 \;, \qquad m_{18}^{eq} = 0 \;.


Although the theoretical derivation of the macroscopic equations yields the non-zero 
expressions shown above for :math:`m_{10}^{eq}` and :math:`m_{12}^{eq}`, the actual 
LBPM MRT implementation enforces these moments to be strictly zero (i.e., :math:`m_{10}^{eq} = 0` 
and :math:`m_{12}^{eq} = 0`) following moments-based approach proposed by :cite:t:`d2002multiple`. 

The relaxation parameters for the momentum diffusion (collision operator) are determined 
based on the effective relaxation time :math:`\tau_{eff}`:

.. math::

   \lambda_1 =  \lambda_2 =  \lambda_9 = \lambda_{10} = \lambda_{11} = \lambda_{12} = \lambda_{13} = \lambda_{14} = \lambda_{15} = s_\nu = \frac{1}{\tau_{eff}} \;,

.. math::

   \lambda_4 = \lambda_6 = \lambda_8 = \lambda_{16} = \lambda_{17} = \lambda_{18} = \frac{8(2-s_\nu)}{8-s_\nu} \;.

where the effective relaxation time :math:`\tau_{eff}` is related to the effective kinematic 
viscosity of the fluid by :math:`\nu_{eff} = (\tau_{eff} - 1/2)/3`. Note that in the generalized 
volume-averaged formulation, this effective viscosity :math:`\nu_{eff}` dictates the viscous 
momentum diffusion (Brinkman term), while the linear drag term (Darcy) depends strictly on the 
thermodynamic fluid kinematic viscosity :math:`\nu`, defined by the base relaxation time :math:`\tau`. 
The relation imposed on the even relaxation times enforces the no-slip velocity condition exactly 
at a position one-half lattice spacing from the solid boundary. Consequently, it leads to a 
viscosity-independent numerical error in the permeability estimate :cite:p:`d2009viscosity`.

.. raw:: html

   </details>

.. raw:: html

   <br>


************************************
Run lbpm_greyscale_simulator
************************************

Launch the simulator from the directory containing the input database and ``.raw`` image:

.. code-block:: bash

   mpirun -np $NUMPROCS lbpm_greyscale_simulator input.db

where ``$NUMPROCS`` must equal the product of the three values in ``nproc``.
For example, ``nproc = 2, 2, 2`` requires eight MPI processes. The executable must
be on the shell search path. The input database must contain four sections: ``Greyscale``, ``Domain``,
``Analysis``, and ``Visualization``. The last two may be empty to use their
defaults.

------------------------------
Model parameters
------------------------------

Set the flow and collision parameters in the ``Greyscale`` class section:

.. list-table:: Greyscale class parameters
   :header-rows: 1
   :widths: 22 20 58
   :align: center

   * - Key
     - Default
     - Meaning
   * - ``collision``
     - ``"MRT"``
     - This routine allow choose differente collision models, on the purpose of 
       tests for LBM developers. The default collision model is set ``"MRT"`` 
       and it recommended model. However, by set ``collision = "MRT"``,
       ``collision = "IMRT"``, and ``collision = "BGK"``.
   * - ``tau``
     - ``1.0``
     - Fluid relaxation time, determining the viscosity used in porous drag
       and the reported permeability.
   * - ``ViscosityRatio``
     - ``1.0``
     - Ratio of effective viscosity to fluid viscosity (:math:`\mathrm{ViscosityRatio} = \displaystyle\frac{\nu_{eff}}{\nu}`), used to derive
       ``tau_eff`` when it is not supplied directly.
   * - ``tau_eff``
     - Derived
     - Effective relaxation time controlling viscous momentum diffusion in
       the MRT collision operator. An explicit value overrides the ratio.
   * - ``F``
     - ``0.0, 0.0, 0.0``
     - External force vector per cubic meter in lattice units, corresponding to
       :math:`F=\rho\boldsymbol{g}=\Delta p / L` in the formulation. Set a nonzero value for
       force-driven periodic flow (``BC = 0``).
   * - ``Forchheimer``
     - ``false``
     - Enable the nonlinear drag term (:eq:`Forchheimer-eq-greysingle`) when ``true``. Linear Darcy drag is
       included in grey voxels in either case.
   * - ``timestepMax``
     - ``100000``
     - Maximum number of lattice timesteps.
   * - ``tolerance``
     - ``0.01``
     - Relative change in the volume-averaged flow between analysis samples
       at which the run stops.
   * - ``din``, ``dout``
     - ``1.0``, ``1.0``
     - Inlet and outlet density values for the pressure boundaries used
       with ``BC = 3``.

The fluid and effective kinematic viscosities in lattice units are

.. math::
   :nowrap:

   $$
   \nu = \frac{\tau - 1/2}{3}, \qquad
   \nu_{\mathrm{eff}} = \frac{\tau_{\mathrm{eff}} - 1/2}{3},
   $$

with the default relation

.. math::
   :nowrap:

   $$
   \tau_{\mathrm{eff}} = \frac{\nu_{eff}}{\nu} (\tau - 1/2) + 1/2.
   $$

Both relaxation times must exceed ``0.5`` for positive viscosity. Values close
to ``0.5`` can be unstable numerically; ``tau = 1.0`` and ``tau_eff = 1.0``
provide a starting point for checking a new geometry.

----------------------------------------
Image labels and grey-region properties
----------------------------------------

The microtomography image must be segmented into impermeable solid regions, 
resolved open pore space, and one or more grey regions. In LBPM, the geometry 
labels are assigned by the ``Domain`` class through the mapping from 
``ReadValues`` to ``WriteValues``, as follows:

- Labels less than or equal to zero identify solid voxels.
- Positive labels identify active voxels, including both open and grey regions.
- Assign porosity ``1.0`` to resolved open pore space. The collision kernels
  disable porous drag at these voxels.
- Assign porosity strictly between zero and one to grey regions, together with
  a positive local permeability.

There are two ways to supply the porous-medium properties.

~~~~~~~~~~~~~~~~~~~~~~~~~
Properties by image label
~~~~~~~~~~~~~~~~~~~~~~~~~

Supply three paired lists in ``Greyscale``:

.. code-block:: text

   ComponentLabels = 1, 2
   PorosityList = 1.0, 0.3
   PermeabilityList = 1.0, 0.01

These values describe open pore space with label ``1`` and a grey region with
label ``2``. All three lists must have equal lengths, and every active label
must appear, including labels introduced by inlet or outlet layers. Each
porosity must satisfy :math:`0 < \epsilon \leq 1`, and each permeability must
be positive. The positive permeability supplied for an open-pore label is a
placeholder because drag is disabled when its porosity is exactly one.

The permeability values specified in ``PermeabilityList`` are expressed in 
square micrometres (:math:`\mu m^{2}`) and internally converted to lattice 
units using ``voxel_length``, which represents the image resolution in 
micrometres per voxel:

.. math::

   K_{\mathrm{lattice}} = \frac{K_{\mathrm{input}}}{h^2},

where :math:`h` is ``voxel_length``. For example, a permeability of
``0.01`` at a resolution of ``2.0`` micrometres per voxel becomes ``0.0025``
in lattice units. Local grey-region permeability is an input property;
the bulk sample permeability is computed from the simulated flow.

~~~~~~~~~~~~~~~~~~~~~~~~~~~
Properties from voxel maps
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, replacing ``ComponentLabels``, ``PorosityList`` and 
``PermeabilityList`` with ``.raw`` input file maps:

.. code-block:: text

   FileVoxelPorosityMap = "porosity.raw", "double"
   FileVoxelPermeabilityMap = "permeability.raw", "double"

each file is a global, headerless array of double-precision values with the
same dimensions, ordering, and alignment as the geometry image. The sequence
of data in input map follow ``index = x + Nx*y + Nx*Ny*z``. Rank zero reads 
the arrays and distributes them using
``Domain.N``, ``Domain.n``, ``Domain.nproc``, and ``Domain.offset``. Only the
datatype ``"double"`` is implemented for these property maps. Permeability
values use the same physical units as ``PermeabilityList``, and the validity
checks apply at every active voxel. When using porosity and permeability maps, 
``ComponentLabels`` is not required.
However, ``ReadValues`` and ``WriteValues`` must still be specified in the
``Domain`` section, even when the geometry contains many labels. To preserve
the existing labels, use identical values in both lists. Porosity and
permeability are assigned to each active voxel using the values at the
corresponding spatial position in the maps.

If both voxel maps (``FileVoxelPorosityMap`` 
and ``FileVoxelPermeabilityMap``) and parameter lists (``PorosityList`` and 
``PermeabilityList``) are provided, which is not recommended, the porosity and 
permeability values will be read exclusively from the voxel maps. The corresponding 
parameter lists will be ignored.

Use property maps with ``BC = 0`` and without generated inlet or outlet mixing
layers. The map reader does not generate property values for added reservoir
or mixing layers; use the label-based method for those configurations.

When ``Domain.Filename`` is supplied, the simulator decomposes the global
geometry image internally. Set ``N`` to the original image dimensions,
``n`` to each rank's interior dimensions, and ``nproc`` to the process grid.
For an uncropped image, the componentwise product of ``n`` and ``nproc``
matches ``N``. ``offset`` selects the starting voxel when reading a subregion.
See :doc:`../domain/domain` for the decomposition settings. Without
``Filename``, the simulator instead expects pre-decomposed ``ID.00000``,
``ID.00001``, and subsequent rank files, including their halo layers.

------------------------------
Boundary Conditions
------------------------------

Set ``BC`` in the ``Domain`` section. The greyscale driver supports the
following flow boundary conditions:

- ``BC = 0`` -- periodic boundaries. Drive flow with a nonzero ``Greyscale.F``.
- ``BC = 3`` -- fixed pressure at the lower and upper z faces, with periodic
  boundaries in x and y. Set ``Greyscale.din`` and ``Greyscale.dout``.

For MRT, the pressure relation is :math:`p = \rho/3`, so ``din`` and ``dout``
are density values rather than pressures in physical units. For example,
``din = 1.001`` and ``dout = 1.0`` impose a lattice pressure difference of
:math:`0.001/3`.

For periodic images, mismatched inlet and outlet pores can restrict flow
(as illustrated in :doc:`../../../userGuide/models/mrt/mrt` and
:doc:`../../../examples/SinglePhasePoreScale/3D-DigitalRocks/3D-DigitalRocks-Cases`).
``InletLayers`` and ``OutletLayers`` can modify the geometry near the domain
faces to improve connectivity. Specify three counts, in x, y, z order.
In the greyscale model, zero or negative labels represent solid voxels,
while positive labels represent active pore or grey regions.

When ``Domain.Filename`` is supplied, the following settings select the
first five and last five existing z planes for modification:

.. code-block:: text

   InletLayers = 0, 0, 5
   OutletLayers = 0, 0, 5

The ``InletLayers`` and ``OutletLayers`` command only works for z planes.
These settings modify existing planes; they do not append additional
slices. The two zeros disable layer generation along x and y. LBPM 
changes only solid voxels whose corresponding
position on the opposite z face has a positive label. It copies that
positive label from the same x and y coordinates. Existing positive labels
remain unchanged, as illustrated below:

.. list-table:: Label changes
   :header-rows: 1
   :align: center

   * - Current label
     - Opposite-face label
     - Result
   * - ``0``
     - ``3``
     - ``3``
   * - ``2``
     - ``3``
     - ``2`` (unchanged)
   * - ``0``
     - ``0``
     - ``0`` (unchanged)

With label-based properties, copied grey labels receive the porosity and
permeability associated with their new labels.

For ``BC = 3``, provide open-pore inlet and outlet reservoirs in the geometry
and use the label-based property method.

.. note::

   The built-in ``AbsPerm`` calculation uses the magnitude of ``F`` as the
   driving acceleration. It does not infer a pressure gradient from ``din``
   and ``dout``. Use the periodic, force-driven setup below for the reported
   permeability; pressure-driven runs require a separate permeability
   calculation from the imposed pressure drop and measured bulk flux.

------------------------------------
Input File Example for ``BC = 0``
------------------------------------

The following ``input.db`` assumes a ``64 x 64 x 64`` raw 8-bit image named
``geometry.raw``. Label ``0`` is solid, label ``1`` is open pore space, and
label ``2`` is a grey region. The property values are illustrative and should
be replaced with values appropriate to the sample.

.. code-block:: text

   Greyscale {
      collision = "MRT"
      tau = 1.0
      ViscosityRatio = 1.0
      Den = 1.0
      F = 0.0, 0.0, 1.0e-6
      Forchheimer = false
      timestepMax = 100000
      tolerance = 0.001
      ComponentLabels = 1, 2
      PorosityList = 1.0, 0.3
      PermeabilityList = 1.0, 0.01 // square micrometres
   }
   Domain {
      Filename = "geometry.raw"
      ReadType = "8bit"
      N = 64, 64, 64        // original image dimensions
      nproc = 1, 1, 1      // MPI process grid
      n = 64, 64, 64       // interior voxels per MPI rank
      offset = 0, 0, 0
      voxel_length = 1.0   // micrometres per voxel
      ReadValues = 0, 1, 2
      WriteValues = 0, 1, 2
      BC = 0
   }
   Analysis {
      analysis_interval = 1000
      visualization_interval = 10000
      restart_interval = 1000000
   }
   Visualization {
      format = "vtk"
   }

Run this example with one MPI process:

.. code-block:: bash

   mpirun -np 1 lbpm_greyscale_simulator input.db

To use eight processes for the same image, change ``nproc`` to ``2, 2, 2``,
``n`` to ``32, 32, 32``, and launch with ``-np 8``. Keep ``N`` unchanged.
The example places checkpoint output beyond the maximum timestep; visualization
and analysis remain active.

--------------------------------------
Assessing Steady State Permeability
--------------------------------------

The ``Analysis`` section controls the sampling and output intervals:

- ``analysis_interval = 1000`` -- compute averaged flow, convergence, and
  permeability every 1000 timesteps by default.
- ``visualization_interval = 1000`` -- write field data every 1000 timesteps
  by default.
- ``restart_interval = 10000`` -- write checkpoint files every 10000 timesteps
  by default.

All three intervals belong in ``Analysis``, including
``visualization_interval``. Use positive, even intervals and an even
``timestepMax`` because the run loop advances an odd/even pair before checking
analysis and output conditions.

Defining :math:`\overline{u}_{\parallel}` as the bulk-volume
average velocity projected along ``F``:

.. math::

   \overline{u}_{\parallel} = \overline{\mathbf{u}}\cdot\frac{\mathbf{F}}{|\mathbf{F}|} =
   \frac{\overline{u}_{x}F_x+\overline{u}_{y}F_y+\overline{u}_{z}F_{z}}{\sqrt{F_x^2+F_y^2+F_z^2}}.

For example, if the force acts only in the **positive z direction**, 
then :math:`\overline{u}_{\parallel}=\overline{u}_z`. 

The convergence error employs :math:`\overline{u}_{\parallel}` is

.. math::

   E_n = \frac{\left|\overline{u}_{\parallel,n}
         - \overline{u}_{\parallel,n-1}\right|}
         {\left|\overline{u}_{\parallel,n}\right|}.

The simulation stops when this sampled relative change is no greater than
``tolerance`` or when ``timestepMax`` is reached. This criterion measures the
change in averaged flow, so inspect the permeability history before treating
the final value as a steady-state measurement. A connected flow path and a
nonzero driving force are needed for a meaningful periodic-flow measurement.
A zero or nonfinite flow, or a reported NaN error, indicates a failed
measurement and should not be interpreted as convergence.

For the periodic, force-driven case, the reported permeability is

.. math::

   K_{\mathrm{sample}} = h^2 \nu
       \frac{\overline{u}_{\parallel}}{|\boldsymbol{F}|}.

The output is in square micrometres when ``voxel_length`` is in micrometres.
The averaged velocity is normalized by the full sample volume, including
solid voxels. With nonlinear Forchheimer drag enabled, this ratio is an
apparent permeability that depends on the driving force.

Rank zero appends the results to ``Permeability.csv``. Despite the extension,
the file is whitespace-delimited and has the columns

.. list-table::
   :widths: auto
   :align: center

   * - ``timestep``
     - ``Fx``
     - ``Fy``
     - ``Fz``
     - ``mu``
     - ``Vs``
     - ``As``
     - ``Hs``
     - ``Xs``
     - ``vax``
     - ``vay``
     - ``vaz``
     - ``AbsPerm``

``AbsPerm`` is the bulk permeability in square micrometres. ``vax``, ``vay``,
and ``vaz`` are the bulk-volume average velocity components in lattice units;
``mu`` contains the kinematic viscosity :math:`(\tau - 1/2)/3`.
``Vs``, ``As``, ``Hs``, and ``Xs`` contain the geometric Minkowski measures.
Use a separate working directory for each run to keep appended histories
from different simulations separate.

For example, plot the permeability history with Python, pandas, and Matplotlib:

.. code-block:: python

   import pandas as pd
   from matplotlib import pyplot as plt

   data = pd.read_csv("Permeability.csv", sep=r"\s+")
   permeability_md = data["AbsPerm"] / 0.0009869233

   plt.plot(data["timestep"], permeability_md)
   plt.xlabel("Timestep")
   plt.ylabel("Permeability (mD)")
   plt.tight_layout()
   plt.show()

If the permeability is still drifting, rerun with a larger ``timestepMax``
and a smaller ``tolerance``. To assess intrinsic permeability in the Darcy
regime, also check that reducing the driving force gives a consistent result.

--------------------------------------
Visualization and Checkpoint Output
--------------------------------------

The ``Visualization.format`` setting selects the field-output format. The
default is ``"silo"``; ``"vtk"`` provides VTK image data, and ``"hdf5"``
requires HDF5 support in the build. Silo output similarly requires a build
with Silo support.

The simulator exports ``SignDist``, ``Velocity_x``, ``Velocity_y``,
``Velocity_z``, and ``Pressure`` on the ``domain`` mesh. Velocity and pressure
are in lattice units; for MRT the pressure is :math:`p = \sum_q f_q/3`.
An output directory such as ``vis10000`` contains the fields for that timestep.
With ``format = "vtk"``, open its ``summary.pvti`` in ParaView. VTK coordinates
use unit voxel spacing; apply ``voxel_length`` when converting lengths to
micrometres. The executable also writes the final fields after the run ends,
even if termination falls between visualization intervals.

****************************
Benchmark Cases
****************************

.. list-table::
   :header-rows: 1
   :widths: 40 40 40
   :align: center

   * - .. centered:: :doc:`../../../examples/SinglePhaseMultiScale/Benchmarks/parallel-series`
     - .. centered:: :doc:`../../../examples/SinglePhaseMultiScale/Benchmarks/poiseuille`
     - .. centered:: :doc:`../../../examples/SinglePhaseMultiScale/Benchmarks/Gray-Multiscale-Mircromodels`

   * - .. image:: ../../../_static/images/permeability_parallel_clean.svg
          :width: 300px
          :align: center

     - .. image:: ../../../_static/images/grayscale_poiseuille_flow_bw.svg
          :width: 300px
          :align: center

     - .. image:: ../../../_static/images/M2-Gray.png
          :width: 200px
          :align: center
