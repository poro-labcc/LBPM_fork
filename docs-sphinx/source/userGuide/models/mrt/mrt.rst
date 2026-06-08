###############################################################################
Single-Phase Pore-Scale Model
###############################################################################

The LBPM single fluid model is implemented by a multi-relaxation time (MRT) collision model in a D3Q19
lattice representation of lattice Boltzmann equation (LBE) to solve the fluid momentum transport, recovering the Navier-Stokes
equations. Beside the MRT scheme a set two-relaxation time (TRT) is used to avoid half-way bounce-back slip effects :cite:p:`pan2006`. The implemented 
modeling is used to assess the permeability of digital rock images in either the Darcy or non-Darcy flow regimes. 

****************************
LBM Formulation
****************************

The LBE governing momentum transport is defined based on a MRT scheme based on the D3Q19 discrete
velocity set:

.. math::
   :nowrap:

   $$
      f_q(\boldsymbol{x}_i + \boldsymbol{\xi}_q \delta t,t + \delta t) - f_q(\boldsymbol{x}_i,t) = \sum^{Q-1}_{k=0} M^{-1}_{qk} \lambda_{k} (m_k^{eq}-m_k) + w_q \boldsymbol{\xi}_q \cdot \frac{\boldsymbol{F}}{c_s^2} \;,
   $$
where :math:`\boldsymbol{F}` is an external body force, :math:`c_s^2 = 1/3` is the speed of sound and :math:`\boldsymbol{\xi}_q` are the unit lattice vectors for the Lattice discretization:

.. math::
   :nowrap:

   $$
      \boldsymbol{\xi}_q=\xi_{q,\alpha} =\left\{
      \begin{array}{llll}
         \alpha=0, & \{0,0,0\}, & \text{for } q = 0, & w_{0}=1/3 \\[4pt]
         \alpha=x, & \{\pm 1,0,0\}, & \text{for } q = 1,2, & w_{1,2}=1/18 \\[4pt]
         \alpha=y, & \{0,\pm 1,0\}, & \text{for } q = 3,4, & w_{3,4}=1/18 \\[4pt]
         \alpha=z, & \{0,0,\pm 1\}, & \text{for } q = 5,6, & w_{5,6}=1/18 \\[4pt]
         \alpha=x,y,& \{\pm 1,\pm 1,0\}, & \text{for } q = 7,8,9,10, & w_{7,8,9,10}=1/36 \\[4pt]
         \alpha=x,z,& \{\pm 1,0,\pm 1\}, & \text{for } q = 11,12,13,14, & w_{11,12,13,14}=1/36 \\[4pt]
         \alpha=y,z,& \{0,\pm 1,\pm 1\}, & \text{for } q = 15,16,17,18, & w_{15,16,17,18}=1/36 .
      \end{array}\right.
   $$


The moments are linearly indepdendent functions of the distributions (Present Grand-Schimidt base):

.. math::
   :nowrap:

   $$
      m_k = \sum_{q=0}^{18} M_{qk} f_q\; \qquad \rightarrow \qquad M_{qk} = \left(
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
      \right),
   $$

and the non-conserved equilibrium moments are then obtained following the moment-based approach proposed by :cite:t:`d2002multiple`, in which selected equilibrium moments are explicitly enforced:

.. math::
   :nowrap:

   $$
     m_1^{eq} = 19\frac{j_x^2+j_y^2+j_z^2}{\rho} - 11\rho \;, \qquad m_2^{eq} = 3\rho - \frac{11}{2} \frac{j_x^2+j_y^2+j_z^2}{\rho} \;,
   $$

.. math::
   :nowrap:

   $$
     m_4^{eq} = -\frac 2 3 j_x \;, \qquad m_6^{eq} = -\frac 2 3 j_y \;, \qquad m_8^{eq} = -\frac 2 3 j_z \;,
   $$

.. math::
   :nowrap:

   $$
     m_9^{eq} = \frac{2j_x^2-j_y^2-j_z^2}{\rho}\;, \qquad m_{10}^{eq} = -\frac{2j_x^2-j_y^2-j_z^2}{2\rho} \;,
   $$

.. math::
   :nowrap:

   $$
     m_{11}^{eq} = \frac{j_y^2-j_z^2}{\rho} \;, \qquad m_{12}^{eq} = -\frac{j_y^2-j_z^2}{2\rho} \;,
   $$

.. math::
   :nowrap:

   $$
     m_{13}^{eq} = \frac{j_x j_y}{\rho} \;, \qquad m_{14}^{eq} = \frac{j_y j_z}{\rho} \;, \qquad m_{15}^{eq} = \frac{j_x j_z}{\rho} \;.
   $$

The relaxation parameters are determined based on the relaxation time :math:`\tau`:

.. math::
   :nowrap:

   $$
     \lambda_1 =  \lambda_2=  \lambda_9 = \lambda_{10}= \lambda_{11}= \lambda_{12}= \lambda_{13}= \lambda_{14}= \lambda_{15} = s_\nu = \frac{1}{\tau} \;,
   $$
.. math::
   :nowrap:
      
   $$
     \lambda_{4}= \lambda_{6}= \lambda_{8} = \lambda_{16} = \lambda_{17} = \lambda_{18}= \frac{8(2-s_\nu)}{8-s_\nu} \;.
   $$

where :math:`\tau` is related to the kinematic viscosity of the fluid by
:math:`\nu = (\tau - 1/2)/3`. The relation imposed on the even relaxation
times enforces the no-slip velocity condition exactly at a position
one-half lattice spacing from the solid boundary. Consequently, it
leads to a viscosity-independent numerical error in the permeability
estimate :cite:p:`d2009viscosity`.

****************************
Run lbpm_permeability_simulator
****************************

A typical command to launch the LBPM single-phase simulator is as follows

```
mpirun -np $NUMPROCS lbpm_permeability_simulator input.db
```

where ``$NUMPROCS`` is the number of MPI processors to be used and ``input.db`` is
the name of the input database that provides the simulation parameters.
Note that the specific syntax to launch MPI tasks may vary depending on your system.
For additional details please refer to your local system documentation.

***************************
Model parameters
***************************

The essential model parameters for the single-phase MRT model are

- ``tau`` -- control the fluid viscosity -- :math:`0.5 < \tau < \infty`

The kinematic viscosity in LBM scale is given by

.. math::
   :nowrap:

   $$
      \nu = \frac{1}{3} \left( \tau - \frac 12 \right).
   $$

Numerical simulations may become unstable for values of :math:`\tau` close to 0.5, corresponding to vanishing viscosity. The stability 
range of :math:`\tau` is investigated in the benchmark tests :doc:`../../../examples/SinglePhasePoreScale/bcc/bcc` and XX. Additionally, 
the parameters governing fluid flow through the medium depend on the selected boundary condition, as discussed in the Boundary Conditions section below.

****************************
Boundary Conditions
****************************

The following external boundary conditions are supported by ``lbpm_permeability_simulator``
and can be set by setting the ``BC`` key values in the ``Domain`` section of the
input file database

- ``BC = 0`` -- fully periodic boundary conditions
- ``BC = 3`` -- constant pressure boundary condition
- ``BC = 4`` -- constant volumetric flux boundary condition

For ``BC = 0`` any mass that exits on one side of the domain will re-enter at the other
side. If the pore-structure for the image is tight, the mismatch between the inlet and
outlet can artificially reduce the permeability of the sample due to the blockage of
flow pathways at the boundary. LBPM includes an internal utility that will reduce the impact
of the boundary mismatch by eroding the solid labels within the inlet and outlet layers
(https://doi.org/10.1007/s10596-020-10028-9) to create a mixing layer.
The number mixing layers to use can be set using the key values in the ``Domain`` section
of the input database

- ``InletLayers  = 5`` -- set the number of mixing layers to ``5`` voxels at the inlet
- ``OUtletLayers  = 5`` -- set the number of mixing layers to ``5`` voxels at the outlet

For the other boundary conditions a thin reservoir of fluid  (default ``3`` voxels)
is established at either side of the domain. The inlet is defined as the boundary face
where ``z = 0`` and the outlet is the boundary face where ``z = nprocz*nz``. By default a
reservoir of fluid A is established at the inlet and a reservoir of fluid B is established at
the outlet, each with a default thickness of three voxels. To over-ride the default label at
the inlet or outlet, the ``Domain`` section of the database may specify the following key values

- ``InletLayerPhase = 2`` -- establish a reservoir of component B at the inlet
- ``OutletLayerPhase = 1`` -- establish a reservoir of component A at the outlet

****************************
Example Input File
****************************

.. code-block:: c

   MRT {
      tau = 1.0
      F = 0.0, 0.0, 1.0e-5
      timestepMax = 2000
      tolerance = 0.01
   }
   Domain {
      Filename = "Bentheimer_LB_sim_intermediate_oil_wet_Sw_0p37.raw"  
      ReadType = "8bit"      // data type
      N = 900, 900, 1600     // size of original image
      nproc = 2, 2, 2        // process grid
      n = 200, 200, 200      // sub-domain size
      offset = 300, 300, 300 // offset to read sub-domain
      voxel_length = 1.66    // voxel length (in microns)
      ReadValues = 0, 1, 2   // labels within the original image
      WriteValues = 0, 1, 2  // associated labels to be used by LBPM
      InletLayers = 0, 0, 10 // specify 10 layers along the z-inlet
      BC = 0                 // boundary condition type (0 for periodic)
   }
   Visualization {
   }

****************************
Benchmark Cases
****************************

.. list-table:: Benchmarks
   :header-rows: 1
   :widths: 30 30 30

   * - :doc:`../../../examples/SinglePhasePoreScale/bcc/bcc`
     - 3D digital Rocks
     - Multiscale Micromodels

   * - .. image:: ../../../_static/images/bcc-bench.png
          :width: 150px
          :align: center

     - .. image:: ../../../_static/images/bentheimer-3d.png
          :width: 150px
          :align: center

     - .. image:: ../../../_static/images/M2-P.png
          :width: 150px
          :align: center
