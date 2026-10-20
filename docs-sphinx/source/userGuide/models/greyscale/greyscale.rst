=============================================
Single-Phase Multi-Scale Porous Media (Greyscale Model - lbpm_greyscale_simulator)
=============================================

The LBPM model to simulate single-phase multi-scale porous media, also called greyscale model, 
is constructed to approximate the solution of the Darcy-Brinkman equations in grey regions, 
coupled to a Navier-Stokes solution in open regions. The generalized model for isothermal fluid 
flow in porous media employs volume-averaged Navier-Stokes equations that account for both linear 
and nonlinear drag effects:

.. math::

   \nabla \cdot \boldsymbol{u} = 0 \;,

.. math::

   \frac{\partial \boldsymbol{u}}{\partial t} + (\boldsymbol{u} \cdot \nabla) \left(\frac{\boldsymbol{u}}{\epsilon}\right) = -\frac{1}{\rho} \nabla(\epsilon p) + \nu_eff \nabla^2 \boldsymbol{u}  -\frac{\epsilon \nu}{K} \boldsymbol{u} - \frac{\epsilon F_e}{\sqrt{K}} |\boldsymbol{u}| \boldsymbol{u} + \epsilon \boldsymbol{g} \;.

In this formulation, :math:`p` and :math:`\boldsymbol{u}` denote the volume-averaged pressure and velocity, 
respectively, and :math:`\rho` represents the fluid density. The model incorporates :math:`\nu` as the 
kinematic viscosity of the fluid and :math:`\nu_eff` as an effective viscosity parameter. External force 
effects are captured through the acceleration term :math:`\boldsymbol{g}`. The geometric function :math:`F_e`, 
derived from Ergun's experimental investigations, characterizes the porous medium structure through the relation:

.. math::

   F_e = \frac{1.75}{\sqrt{150\epsilon^3}} \;.

****************************
LBM Formulation
****************************

.. raw:: html

   <details>
   <summary><strong>Click here to open: LBM Formulation</strong></summary>

A multiple-relaxation-time (MRT) lattice Boltzmann equation with D3Q19 discretization describes momentum transport for the phase-averaged equations: 

.. math::

   f_q(\boldsymbol{x}_i + \boldsymbol{\xi}_q \delta_t,t + \delta_t) - f_q(\boldsymbol{x}_i,t) = \sum^{Q-1}_{k=0} M^{-1}_{qk} S_{kk} (m_k^{eq}-m_k)  + \sum^{Q-1}_{k=0} M^{-1}_{qk} \left(1-\frac{S_{kk}}{2}\right) \hat{F}_q \;.

The total force :math:`\boldsymbol{F}` is imposed based on the discrete force construction developed by Guo et al. (2002):

.. math::

   \boldsymbol{F} = - \frac{\epsilon \nu \rho}{K} \boldsymbol{u} - \frac{\epsilon F_{e}}{\sqrt{K}} \rho |\boldsymbol{u}|\boldsymbol{u} +  \epsilon \rho \boldsymbol{g} \;.

The fluid velocity is computed through the relation for the forcing scheme:

.. math::

   \rho \boldsymbol{u} = \sum_i \boldsymbol{\xi}_i f_i + \frac{\delta_t}{2}  \boldsymbol{F} \;.

However, to handle the nonlinearity arising from the velocity-dependent force term, an explicit solution for the velocity is employed:

.. math::

   \boldsymbol{u} = \frac{\boldsymbol{v}}{c_0 + \sqrt{c_0^2 + c_1 |\boldsymbol{v}|}} \;,

that utilizes a temporal velocity field :math:`\boldsymbol{v}`:

.. math::

   \rho \boldsymbol{v} = \sum_i \boldsymbol{\xi}_i f_i + \frac{\delta_t}{2} \rho \epsilon \boldsymbol{g} \;,

and the coefficients governing this nonlinear transformation are:

.. math::

   c_0 = \frac{1}{2} \left(  1 + \frac{\delta_t}{2}\frac{\epsilon \nu}{K} \right), \qquad c_1 = \frac{\delta_t}{2}\frac{\epsilon F_e}{\sqrt{K}} \;.

where :math:`c_s^2 = 1/3` is the speed of sound and :math:`\boldsymbol{\xi}_q` are the unit lattice vectors for the D3Q19 lattice discretization:

.. math::

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

The moments defined in the Gram-Schmidt basis are linearly independent functions of the distribution functions:

.. math::

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
   \right)

The equilibrium formulation is adapted from the generalized volume-averaged macroscopic equations. As established in the appendix of Guo and Zhao (2002), the definition of the equilibrium distribution function yields the following moments, where the presence of the porous matrix modifies the second-order velocity terms by a factor of :math:`1/\epsilon`:

.. math::

   \sum_i f_i^{(eq)} = \rho \;,

.. math::

   \sum_i \xi_{i\alpha} f_i^{(eq)} = \rho u_\alpha \;,

.. math::

   \sum_i \xi_{i\alpha} \xi_{i\beta} f_i^{(eq)} = c_s^2 \rho \delta_{\alpha\beta} + \frac{1}{\epsilon} \rho u_\alpha u_\beta \;,

.. math::

   \sum_i \xi_{i\alpha} \xi_{i\beta} \xi_{i\gamma} f_i^{(eq)} = c_s^2 \rho (u_\alpha \delta_{\beta\gamma} + u_\beta \delta_{\alpha\gamma} + u_\gamma \delta_{\alpha\beta}) \;.

Applying this to the D3Q19 MRT basis yields the corresponding equilibrium moments, given by:

.. math::

   m_0^{eq} = \rho \;,

.. math::

   m_1^{eq} = \frac{19}{\rho \epsilon} \left(j_x^2+j_y^2+j_z^2\right) - 11\rho \;, \qquad m_2^{eq} = 3\rho - \frac{11}{2 \rho \epsilon}\left(j_x^2+j_y^2+j_z^2\right) \;,

.. math::

   m_3^{eq} = j_x = \rho u_x \;, \qquad m_5^{eq} = j_y = \rho u_y \;, \qquad m_7^{eq} = j_z = \rho u_z \;,

.. math::

   m_4^{eq} = -\frac{2}{3} j_x \;, \qquad m_6^{eq} = -\frac{2}{3} j_y \;, \qquad m_8^{eq} = -\frac{2}{3} j_z \;,

.. math::

   m_9^{eq} = \frac{1}{\rho \epsilon}\left(2j_x^2-j_y^2-j_z^2\right) \;, \qquad m_{10}^{eq} = -\frac{1}{2 \rho \epsilon}\left(2j_x^2-j_y^2-j_z^2\right) \;,

.. math::

   m_{11}^{eq} = \frac{1}{\rho \epsilon} \left(j_y^2-j_z^2\right) \;, \qquad m_{12}^{eq} = -\frac{1}{2 \rho \epsilon}\left(j_y^2-j_z^2\right) \;,

.. math::

   m_{13}^{eq} = \frac{1}{\rho \epsilon}(j_x j_y) \;, \qquad m_{14}^{eq} = \frac{1}{\rho \epsilon}(j_y j_z) \;, \qquad m_{15}^{eq} = \frac{1}{\rho \epsilon} (j_x j_z) \;,

.. math::

   m_{16}^{eq} = 0 \;, \qquad m_{17}^{eq} = 0 \;, \qquad m_{18}^{eq} = 0 \;.

.. note::
   Although the theoretical derivation of the macroscopic equations yields the non-zero expressions shown above for :math:`m_{10}^{eq}` and :math:`m_{12}^{eq}`, the actual LBPM MRT implementation enforces these moments to be strictly zero (i.e., :math:`m_{10}^{eq} = 0` and :math:`m_{12}^{eq} = 0`). The theoretical values remain commented out in the source code.

The relaxation parameters for the momentum diffusion (collision operator) are determined based on the effective relaxation time :math:`\tau_{eff}`:

.. math::

   \lambda_1 =  \lambda_2 =  \lambda_9 = \lambda_{10} = \lambda_{11} = \lambda_{12} = \lambda_{13} = \lambda_{14} = \lambda_{15} = s_\nu = \frac{1}{\tau_{eff}} \;,

.. math::

   \lambda_4 = \lambda_6 = \lambda_8 = \lambda_{16} = \lambda_{17} = \lambda_{18} = \frac{8(2-s_\nu)}{8-s_\nu} \;.

where the effective relaxation time :math:`\tau_{eff}` is related to the effective kinematic viscosity of the fluid by :math:`\nu_{eff} = (\tau_{eff} - 1/2)/3`. Note that in the generalized volume-averaged formulation, this effective viscosity :math:`\nu_{eff}` dictates the viscous momentum diffusion (Brinkman term), while the linear drag term (Darcy) depends strictly on the thermodynamic fluid kinematic viscosity :math:`\nu`, defined by the base relaxation time :math:`\tau`. The relation imposed on the even relaxation times enforces the no-slip velocity condition exactly at a position one-half lattice spacing from the solid boundary. Consequently, it leads to a viscosity-independent numerical error in the permeability estimate.

.. raw:: html

   </details>

.. raw:: html

   <br>


************************************
Run lbpm_greyscale_simulator
************************************

To use the lbpm_greyscale_simulator, the input image
should be segmented into "grey" and open regions. Each particular "grey" label should be assigned 
both a porosity and permeability value. 

A typical command to launch the LBPM color simulator is as follows

```
mpirun -np $NUMPROCS lbpm_greyscale_simulator input.db
```

where ``$NUMPROCS`` is the number of MPI processors to be used and ``input.db`` is
the name of the input database that provides the simulation parameters.
Note that the specific syntax to launch MPI tasks may vary depending on your system.
For additional details please refer to your local system documentation.

***************************
Model parameters
***************************

The essential model parameters for the single phase greyscale model are

- ``tau`` -- control the fluid viscosity -- :math:`0.7 < \tau < 1.5`

The kinematic viscosity is given by