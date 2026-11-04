*********************************************
Parameters Scaling and ``Input.db`` Generator
*********************************************

Convert physical fluid properties to lattice units before configuring a Color
model simulation. The calculator below treats fractional flow, core flooding,
and centrifuge separately because their driving conditions and capillary-number
conventions differ. Physical quantities carry a superscript :math:`*`.

Shared Color model scaling
==========================

Enter the voxel resolution in micrometres; the calculator converts it to
:math:`\Delta x` in metres. Choose a reference density
:math:`\rho_{\mathrm{ref}}` in :math:`\mathrm{kg/m^3}` and a lattice interfacial
tension :math:`\sigma`. These choices define the mass unit and time step:

.. math::

   m_{\mathrm{ref}} = \rho_{\mathrm{ref}}\Delta x^3,
   \qquad
   \Delta t = \sqrt{\frac{\sigma\rho_{\mathrm{ref}}\Delta x^3}{\sigma^*}},
   \qquad
   \alpha = \frac{\sigma}{6}.

For each fluid :math:`i \in \{A,B\}`,

.. math::

   \rho_i = \frac{\rho_i^*}{\rho_{\mathrm{ref}}},
   \qquad
   \nu_i = \frac{\nu_i^*\Delta t}{\Delta x^2},
   \qquad
   \mu_i = \rho_i\nu_i,
   \qquad
   \tau_i = \frac{1}{2} + 3\nu_i.

The remaining conversions are

.. math::

   L = \frac{L^*}{\Delta x},
   \qquad
   u = \frac{u^*\Delta t}{\Delta x},
   \qquad
   Q = \frac{Q^*\Delta t}{\Delta x^3},

.. math::

   p = \frac{p^*\Delta t^2}{\rho_{\mathrm{ref}}\Delta x^2},
   \qquad
   F_z = \frac{G_z^*\Delta t^2}{\rho_{\mathrm{ref}}\Delta x},
   \qquad
   G_z^* = -\frac{\partial p^*}{\partial z^*}.

Here :math:`Q` is a volumetric flow rate, :math:`\mu` is dynamic viscosity,
and :math:`F_z` is a body-force density, equivalent to a driving pressure
gradient. Lattice pressure units correspond to
:math:`p_0 = \rho_{\mathrm{ref}}\Delta x^2/\Delta t^2
= \sigma^*/(\sigma\Delta x)` pascals. Lattice interfacial tension is expressed
in lattice mass and time units, not :math:`\mathrm{N/m}`.

Increasing :math:`\sigma` increases :math:`\Delta t`, but it also changes the
relaxation times. There is no universal stability guarantee from the scaling
alone. The :doc:`Color model guide <../../userGuide/models/color/index>`
gives approximate working ranges :math:`\alpha \lesssim 0.01` and
:math:`0.7 \lesssim \tau_i \lesssim 1.5`; check the selected parameters and
flow speeds for your geometry.

Protocol-specific scaling
=========================

For flow along :math:`z`, let :math:`A=N_xN_y` be the full cross-sectional
area in lattice units and :math:`\phi` the fluid porosity. Enter the dimensions
of the simulated region, including any padding that contributes to that area.
For true volumetric flow rates, the single-phase Reynolds number is
:math:`Re=L|Q|/(A\nu)`.

Interactive calculator
======================

Each Color model protocol has independent inputs and outputs. Choose **Compute**
to update its results or **Clear** to reset only its output fields. Warnings
identify parameter choices outside the approximate working ranges above. The
``input.db`` block below each panel is populated from its calculated parameters.

.. raw:: html

   <style>
   #color-scaling-calculator .scaling-columns {display:flex; gap:2rem; flex-wrap:wrap; align-items:flex-start;}
   #color-scaling-calculator .scaling-column {flex:1 1 28rem; min-width:0;}
   #color-scaling-calculator .scaling-field {display:grid; grid-template-columns:minmax(12rem,1fr) 11rem; gap:0.75rem; align-items:center; margin-bottom:0.6rem;}
   #color-scaling-calculator input {min-width:0; width:100%; padding:0.35rem;}
   #color-scaling-calculator input[type="radio"] {width:auto;}
   #color-scaling-calculator input[readonly] {background:#f3f5f7;}
   #color-scaling-calculator .scaling-options {display:flex; gap:1rem; flex-wrap:wrap; margin:1rem 0;}
   #color-scaling-calculator .scaling-options label {display:flex; gap:0.4rem; align-items:center;}
   #color-scaling-calculator .scaling-actions {display:flex; gap:0.7rem; margin:1rem 0;}
   #color-scaling-calculator .scaling-actions button {padding:0.4rem 0.8rem;}
   #color-scaling-calculator .scaling-status {white-space:pre-line; margin-top:1rem;}
   #color-scaling-calculator .scaling-error {color:#a11;}
   @media (max-width:600px) {#color-scaling-calculator .scaling-field {grid-template-columns:1fr;}}
   </style>
   <div id="color-scaling-calculator">
   <h3>Select Routine</h3>
   <div class="scaling-options">
   <label><input type="radio" name="routine" value="routine4" onchange="changeRoutine()">Routine: Permeability Simulator</label>
   <label><input type="radio" name="routine" value="routine5" onchange="changeRoutine()">Routine: Gray-Scale Permeability Simulator</label>
   <label><input type="radio" name="routine" value="routine6" onchange="changeRoutine()">Routine: Morphological Drainage</label>
   </div>
   <h3>Routine: Color Model</h3>
   <div class="scaling-options">
   <label><input type="radio" name="routine" value="routine1" checked onchange="changeRoutine()">Protocol: Fractional Flow</label>
   <label><input type="radio" name="routine" value="routine2" onchange="changeRoutine()">Protocol: Core Flooding</label>
   <label><input type="radio" name="routine" value="routine3" onchange="changeRoutine()">Protocol: Centrifuge</label>
   </div>
   <div id="routine1Box" class="routine-box">
   <h4>Fractional Flow Protocol</h4>

.. container:: protocol-explanation

   The :doc:`fractional-flow protocol
   <../../userGuide/models/color/protocols/fractionalFlow>` uses periodic
   boundaries (``BC = 0``). Define the conventional, viscosity-weighted capillary
   number using superficial phase velocities :math:`q_i=Q_i/A`:

   .. math::

      Ca = \frac{|\mu_A Q_A + \mu_B Q_B|}{A\sigma}.

   For an assumed co-current volumetric split :math:`f_B=Q_B/(Q_A+Q_B)`, the
   corresponding flow-rate magnitudes are

   .. math::

      |Q| = \frac{A\sigma Ca}{\mu_A(1-f_B)+\mu_B f_B},
      \qquad
      |Q_A|=(1-f_B)|Q|,
      \qquad
      |Q_B|=f_B|Q|.

   The single-phase endpoints are :math:`|Q_A|_{\max}=A\sigma Ca/\mu_A` and
   :math:`|Q_B|_{\max}=A\sigma Ca/\mu_B`. Porosity does not appear because these
   velocities are superficial. The split is a planning assumption, not a constraint
   imposed by the fractional-flow algorithm; the actual split depends on the
   saturation and pore geometry.

   In this checkout, the fractional-flow controller in ``models/ColorModel.cpp``
   normalizes by :math:`5.796\alpha`, while the interfacial-tension conversion uses
   :math:`\sigma=6\alpha`. To achieve the conventional input :math:`Ca`, set

   .. math::

      \mathtt{Color.capillary\_number} = Ca\frac{6}{5.796}.

   The calculator returns this adjusted database value. An initial driving pressure
   gradient is a separate input: a capillary number alone cannot determine the
   body force without a permeability or a measured flow response. LBPM subsequently
   adjusts the force toward its target capillary number.

.. raw:: html

   <div class="scaling-columns">
   <div class="scaling-column">
   <h4>Input Parameters</h4>
   <h5>Physical properties and lattice scale</h5>
   <div class="scaling-field"><label for="viscA_r1">Fluid A kinematic viscosity <span class="math">\(\nu_A^*\)</span> (m²/s)</label><input id="viscA_r1" aria-label="Fluid A kinematic viscosity" type="number" value="0.00000125" step="any"></div>
   <div class="scaling-field"><label for="viscB_r1">Fluid B kinematic viscosity <span class="math">\(\nu_B^*\)</span> (m²/s)</label><input id="viscB_r1" aria-label="Fluid B kinematic viscosity" type="number" value="0.00000125" step="any"></div>
   <div class="scaling-field"><label for="rhoA_real_r1">Fluid A density <span class="math">\(\rho_A^*\)</span> (kg/m³)</label><input id="rhoA_real_r1" aria-label="Fluid A density" type="number" value="1000" step="any"></div>
   <div class="scaling-field"><label for="rhoB_real_r1">Fluid B density <span class="math">\(\rho_B^*\)</span> (kg/m³)</label><input id="rhoB_real_r1" aria-label="Fluid B density" type="number" value="1000" step="any"></div>
   <div class="scaling-field"><label for="rho_ref_r1">Reference density <span class="math">\(\rho_{\mathrm{ref}}\)</span> (kg/m³)</label><input id="rho_ref_r1" aria-label="Reference density" type="number" value="1000" step="any"></div>
   <div class="scaling-field"><label for="int_ten_r1">Physical interfacial tension <span class="math">\(\sigma^*\)</span> (N/m)</label><input id="int_ten_r1" aria-label="Physical interfacial tension" type="number" value="0.040" step="any"></div>
   <div class="scaling-field"><label for="dx_r1">Voxel resolution <span class="math">\(\Delta x\)</span> (µm)</label><input id="dx_r1" aria-label="Voxel resolution" type="number" value="1.66" step="any"></div>
   <div class="scaling-field"><label for="int_ten_lbm_r1">Lattice interfacial tension <span class="math">\(\sigma\)</span> (lattice units)</label><input id="int_ten_lbm_r1" aria-label="Lattice interfacial tension" type="number" value="0.06" step="any"></div>
   <h5>Protocol and domain parameters</h5>
   <div class="scaling-field"><label for="Ca_r1">Weighted superficial capillary number <span class="math">\(Ca\)</span> (dimensionless)</label><input id="Ca_r1" aria-label="Weighted superficial capillary number" type="number" value="0.001" step="any"></div>
   <div class="scaling-field"><label for="fB_r1">Assumed co-current flow fraction of B <span class="math">\(f_B\)</span> (0 to 1)</label><input id="fB_r1" aria-label="Assumed co-current flow fraction of B" type="number" value="0.5" step="any"></div>
   <div class="scaling-field"><label for="Lx_r1">Simulated x dimension <span class="math">\(N_x\)</span> (voxels)</label><input id="Lx_r1" aria-label="Simulated x dimension" type="number" value="900" step="any"></div>
   <div class="scaling-field"><label for="Ly_r1">Simulated y dimension <span class="math">\(N_y\)</span> (voxels)</label><input id="Ly_r1" aria-label="Simulated y dimension" type="number" value="900" step="any"></div>
   <div class="scaling-field"><label for="Lz_r1">Simulated z dimension <span class="math">\(N_z\)</span> (voxels)</label><input id="Lz_r1" aria-label="Simulated z dimension" type="number" value="1600" step="any"></div>
   <div class="scaling-field"><label for="gradient_r1">Initial driving pressure gradient <span class="math">\(G_z^*\)</span> (Pa/m)</label><input id="gradient_r1" aria-label="Initial driving pressure gradient" type="number" value="1000" step="any"></div>
   <div class="scaling-actions"><button type="button" onclick="compute()">Compute</button><button type="button" onclick="clearFields()">Clear</button></div>
   </div>
   <div class="scaling-column">
   <h4>Output Parameters</h4>
   <h5>Shared Color model parameters</h5>
   <div class="scaling-field"><label for="dt_r1">Physical time step <span class="math">\(\Delta t\)</span> (s)</label><input id="dt_r1" aria-label="Physical time step" type="text" readonly></div>
   <div class="scaling-field"><label for="tauA_r1">Fluid A relaxation time <span class="math">\(\tau_A\)</span> (time steps)</label><input id="tauA_r1" aria-label="Fluid A relaxation time" type="text" readonly></div>
   <div class="scaling-field"><label for="tauB_r1">Fluid B relaxation time <span class="math">\(\tau_B\)</span> (time steps)</label><input id="tauB_r1" aria-label="Fluid B relaxation time" type="text" readonly></div>
   <div class="scaling-field"><label for="alpha_r1">Color model alpha <span class="math">\(\alpha\)</span> (lattice units)</label><input id="alpha_r1" aria-label="Color model alpha" type="text" readonly></div>
   <div class="scaling-field"><label for="rhoAlb_r1">Fluid A density <span class="math">\(\rho_A\)</span> (lattice units)</label><input id="rhoAlb_r1" aria-label="Fluid A density" type="text" readonly></div>
   <div class="scaling-field"><label for="rhoBlb_r1">Fluid B density <span class="math">\(\rho_B\)</span> (lattice units)</label><input id="rhoBlb_r1" aria-label="Fluid B density" type="text" readonly></div>
   <div class="scaling-field"><label for="viscAlb_r1">Fluid A kinematic viscosity <span class="math">\(\nu_A\)</span> (lu²/ts)</label><input id="viscAlb_r1" aria-label="Fluid A kinematic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="viscBlb_r1">Fluid B kinematic viscosity <span class="math">\(\nu_B\)</span> (lu²/ts)</label><input id="viscBlb_r1" aria-label="Fluid B kinematic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="muAlb_r1">Fluid A dynamic viscosity <span class="math">\(\mu_A\)</span> (lattice units)</label><input id="muAlb_r1" aria-label="Fluid A dynamic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="muBlb_r1">Fluid B dynamic viscosity <span class="math">\(\mu_B\)</span> (lattice units)</label><input id="muBlb_r1" aria-label="Fluid B dynamic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="pressure_scale_r1">Physical pressure per lattice unit <span class="math">\(p_0\)</span> (Pa)</label><input id="pressure_scale_r1" aria-label="Physical pressure per lattice unit" type="text" readonly></div>
   <h5>Protocol parameters</h5>
   <div class="scaling-field"><label for="capillary_number_r1">Adjusted Color.capillary_number (dimensionless)</label><input id="capillary_number_r1" aria-label="Adjusted Color.capillary_number" type="text" readonly></div>
   <div class="scaling-field"><label for="Qtotal_r1">Total volumetric flow <span class="math">\(Q\)</span> (lu³/ts)</label><input id="Qtotal_r1" aria-label="Total volumetric flow" type="text" readonly></div>
   <div class="scaling-field"><label for="QA_r1">Fluid A volumetric flow <span class="math">\(Q_A\)</span> (lu³/ts)</label><input id="QA_r1" aria-label="Fluid A volumetric flow" type="text" readonly></div>
   <div class="scaling-field"><label for="QB_r1">Fluid B volumetric flow <span class="math">\(Q_B\)</span> (lu³/ts)</label><input id="QB_r1" aria-label="Fluid B volumetric flow" type="text" readonly></div>
   <div class="scaling-field"><label for="QAmax_r1">Single-phase A endpoint magnitude <span class="math">\(|Q_A|_{\max}\)</span> (lu³/ts)</label><input id="QAmax_r1" aria-label="Single-phase A endpoint magnitude" type="text" readonly></div>
   <div class="scaling-field"><label for="QBmax_r1">Single-phase B endpoint magnitude <span class="math">\(|Q_B|_{\max}\)</span> (lu³/ts)</label><input id="QBmax_r1" aria-label="Single-phase B endpoint magnitude" type="text" readonly></div>
   <div class="scaling-field"><label for="dpdl_r1">Initial body-force density <span class="math">\(F_z\)</span> (lattice units)</label><input id="dpdl_r1" aria-label="Initial body-force density" type="text" readonly></div>
   </div>
   </div>
   <p id="status_r1" class="scaling-status" role="status" aria-live="polite"></p>

.. rubric:: Generated input.db

Click **Compute** to populate this block, then copy the text to ``input.db``.
The image filename, label mapping, wetting settings, run length, and
single-process decomposition are example settings to customize.

.. code-block:: text
   :name: input-db-r1
   :class: generated-input-db

   // Click Compute to generate input.db for this protocol.

.. raw:: html

   </div>
   <div id="routine2Box" class="routine-box" style="display:none;">
   <h4>Core Flooding Protocol</h4>

.. container:: protocol-explanation

   The :doc:`core-flooding protocol
   <../../userGuide/models/color/protocols/coreFlooding>` uses ``BC = 4``.
   The local implementation defines its target using fluid **B's dynamic
   viscosity** and an intrinsic pore velocity:

   .. math::

      Ca = \frac{\mu_B u_{\mathrm{pore}}}{\sigma},
      \qquad
      Q_z = \phi A u_{\mathrm{pore}}
          = \frac{\phi A\sigma Ca}{\mu_B}.

   This convention differs from fractional flow. The calculator follows LBPM's
   ``capillary_number`` initialization, even when the two fluids have different
   densities or viscosities. Fluid A is the default inlet phase; its viscosity
   must not replace :math:`\mu_B` in this initialization formula.

   The number of steps to inject one pore volume is rounded up:

   .. math::

      n_{\mathrm{PV}} = \left\lceil\frac{\phi N_xN_yN_z}{Q_z}\right\rceil.

   The runtime estimate assumes the entered MLUPS counts **active fluid-node
   updates across all MPI processes**, not all image voxels. Sum per-process
   throughputs before entering the value. With :math:`M` in million updates per
   second, the estimate is

   .. math::

      t_{\mathrm{wall}} \approx
      \frac{\phi N_xN_yN_z\,n_{\mathrm{PV}}}{10^6 M}.

   This excludes initialization, analysis, and output overhead. The physical
   injection time is :math:`n_{\mathrm{PV}}\Delta t`, which is distinct from wall
   clock time.

.. raw:: html

   <figure><img src="../../_static/images/core-flood-illus.svg" alt="Core flooding with imposed volumetric flux along the z-axis." style="max-width:90%; height:auto;"><figcaption>Core flooding with imposed volumetric flux.</figcaption></figure>
   <div class="scaling-columns">
   <div class="scaling-column">
   <h4>Input Parameters</h4>
   <h5>Physical properties and lattice scale</h5>
   <div class="scaling-field"><label for="viscA_r2">Fluid A kinematic viscosity <span class="math">\(\nu_A^*\)</span> (m²/s)</label><input id="viscA_r2" aria-label="Fluid A kinematic viscosity" type="number" value="0.00000125" step="any"></div>
   <div class="scaling-field"><label for="viscB_r2">Fluid B kinematic viscosity <span class="math">\(\nu_B^*\)</span> (m²/s)</label><input id="viscB_r2" aria-label="Fluid B kinematic viscosity" type="number" value="0.00000125" step="any"></div>
   <div class="scaling-field"><label for="rhoA_real_r2">Fluid A density <span class="math">\(\rho_A^*\)</span> (kg/m³)</label><input id="rhoA_real_r2" aria-label="Fluid A density" type="number" value="1000" step="any"></div>
   <div class="scaling-field"><label for="rhoB_real_r2">Fluid B density <span class="math">\(\rho_B^*\)</span> (kg/m³)</label><input id="rhoB_real_r2" aria-label="Fluid B density" type="number" value="1000" step="any"></div>
   <div class="scaling-field"><label for="rho_ref_r2">Reference density <span class="math">\(\rho_{\mathrm{ref}}\)</span> (kg/m³)</label><input id="rho_ref_r2" aria-label="Reference density" type="number" value="1000" step="any"></div>
   <div class="scaling-field"><label for="int_ten_r2">Physical interfacial tension <span class="math">\(\sigma^*\)</span> (N/m)</label><input id="int_ten_r2" aria-label="Physical interfacial tension" type="number" value="0.040" step="any"></div>
   <div class="scaling-field"><label for="dx_r2">Voxel resolution <span class="math">\(\Delta x\)</span> (µm)</label><input id="dx_r2" aria-label="Voxel resolution" type="number" value="1.66" step="any"></div>
   <div class="scaling-field"><label for="int_ten_lbm_r2">Lattice interfacial tension <span class="math">\(\sigma\)</span> (lattice units)</label><input id="int_ten_lbm_r2" aria-label="Lattice interfacial tension" type="number" value="0.06" step="any"></div>
   <h5>Protocol and domain parameters</h5>
   <div class="scaling-field"><label for="Ca_r2">Pore-velocity capillary number (fluid B) <span class="math">\(Ca\)</span> (dimensionless)</label><input id="Ca_r2" aria-label="Pore-velocity capillary number (fluid B)" type="number" value="0.001" step="any"></div>
   <div class="scaling-field"><label for="Lx_r2">Simulated x dimension <span class="math">\(N_x\)</span> (voxels)</label><input id="Lx_r2" aria-label="Simulated x dimension" type="number" value="900" step="any"></div>
   <div class="scaling-field"><label for="Ly_r2">Simulated y dimension <span class="math">\(N_y\)</span> (voxels)</label><input id="Ly_r2" aria-label="Simulated y dimension" type="number" value="900" step="any"></div>
   <div class="scaling-field"><label for="Lz_r2">Simulated z dimension <span class="math">\(N_z\)</span> (voxels)</label><input id="Lz_r2" aria-label="Simulated z dimension" type="number" value="1600" step="any"></div>
   <div class="scaling-field"><label for="phi_r2">Fluid porosity <span class="math">\(\phi\)</span> (0 to 1)</label><input id="phi_r2" aria-label="Fluid porosity" type="number" value="0.23" step="any"></div>
   <div class="scaling-field"><label for="mlups_lbm_r2">Aggregate active-fluid throughput <span class="math">\(M\)</span> (MLUPS)</label><input id="mlups_lbm_r2" aria-label="Aggregate active-fluid throughput" type="number" value="1000" step="any"></div>
   <div class="scaling-actions"><button type="button" onclick="compute()">Compute</button><button type="button" onclick="clearFields()">Clear</button></div>
   </div>
   <div class="scaling-column">
   <h4>Output Parameters</h4>
   <h5>Shared Color model parameters</h5>
   <div class="scaling-field"><label for="dt_r2">Physical time step <span class="math">\(\Delta t\)</span> (s)</label><input id="dt_r2" aria-label="Physical time step" type="text" readonly></div>
   <div class="scaling-field"><label for="tauA_r2">Fluid A relaxation time <span class="math">\(\tau_A\)</span> (time steps)</label><input id="tauA_r2" aria-label="Fluid A relaxation time" type="text" readonly></div>
   <div class="scaling-field"><label for="tauB_r2">Fluid B relaxation time <span class="math">\(\tau_B\)</span> (time steps)</label><input id="tauB_r2" aria-label="Fluid B relaxation time" type="text" readonly></div>
   <div class="scaling-field"><label for="alpha_r2">Color model alpha <span class="math">\(\alpha\)</span> (lattice units)</label><input id="alpha_r2" aria-label="Color model alpha" type="text" readonly></div>
   <div class="scaling-field"><label for="rhoAlb_r2">Fluid A density <span class="math">\(\rho_A\)</span> (lattice units)</label><input id="rhoAlb_r2" aria-label="Fluid A density" type="text" readonly></div>
   <div class="scaling-field"><label for="rhoBlb_r2">Fluid B density <span class="math">\(\rho_B\)</span> (lattice units)</label><input id="rhoBlb_r2" aria-label="Fluid B density" type="text" readonly></div>
   <div class="scaling-field"><label for="viscAlb_r2">Fluid A kinematic viscosity <span class="math">\(\nu_A\)</span> (lu²/ts)</label><input id="viscAlb_r2" aria-label="Fluid A kinematic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="viscBlb_r2">Fluid B kinematic viscosity <span class="math">\(\nu_B\)</span> (lu²/ts)</label><input id="viscBlb_r2" aria-label="Fluid B kinematic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="muAlb_r2">Fluid A dynamic viscosity <span class="math">\(\mu_A\)</span> (lattice units)</label><input id="muAlb_r2" aria-label="Fluid A dynamic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="muBlb_r2">Fluid B dynamic viscosity <span class="math">\(\mu_B\)</span> (lattice units)</label><input id="muBlb_r2" aria-label="Fluid B dynamic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="pressure_scale_r2">Physical pressure per lattice unit <span class="math">\(p_0\)</span> (Pa)</label><input id="pressure_scale_r2" aria-label="Physical pressure per lattice unit" type="text" readonly></div>
   <h5>Protocol parameters</h5>
   <div class="scaling-field"><label for="Qz_lb_r2">LBPM volumetric flux <span class="math">\(Q_z\)</span> (lu³/ts)</label><input id="Qz_lb_r2" aria-label="LBPM volumetric flux" type="text" readonly></div>
   <div class="scaling-field"><label for="Qz_real_r2">Physical volumetric flow <span class="math">\(Q_z^*\)</span> (m³/s)</label><input id="Qz_real_r2" aria-label="Physical volumetric flow" type="text" readonly></div>
   <div class="scaling-field"><label for="pv_ft_r2">Steps for one pore volume (rounded up) <span class="math">\(n_{\mathrm{PV}}\)</span> (time steps)</label><input id="pv_ft_r2" aria-label="Steps for one pore volume (rounded up)" type="text" readonly></div>
   <div class="scaling-field"><label for="pv_time_r2">Physical injection time for one PV (s)</label><input id="pv_time_r2" aria-label="Physical injection time for one PV" type="text" readonly></div>
   <div class="scaling-field"><label for="time_sim_seconds_r2">Estimated wall time for one PV (s)</label><input id="time_sim_seconds_r2" aria-label="Estimated wall time for one PV" type="text" readonly></div>
   <div class="scaling-field"><label for="time_sim_r2">Estimated wall time for one PV (h / min / s)</label><input id="time_sim_r2" aria-label="Estimated wall time for one PV" type="text" readonly></div>
   </div>
   </div>
   <p id="status_r2" class="scaling-status" role="status" aria-live="polite"></p>

.. rubric:: Generated input.db

Click **Compute** to populate this block, then copy the text to ``input.db``.
The image filename, label mapping, wetting settings, run length, and
single-process decomposition are example settings to customize.

The template sets ``flux`` explicitly. Its ``capillary_number`` line is
commented out so LBPM preserves the calculated flux.

.. code-block:: text
   :name: input-db-r2
   :class: generated-input-db

   // Click Compute to generate input.db for this protocol.

.. raw:: html

   </div>
   <div id="routine3Box" class="routine-box" style="display:none;">
   <h4>Centrifuge Protocol</h4>

.. container:: protocol-explanation

   The :doc:`centrifuge protocol
   <../../userGuide/models/color/protocols/centrifuge>` uses ``BC = 3`` with equal
   reservoir settings ``din = dout = 1.0``. A body force drives displacement;
   the reservoir pressure difference remains zero.

   Enter a signed, physical **equivalent driving pressure span** :math:`P_c^*`
   across the length :math:`N_z\Delta x`. Its lattice value and corresponding
   initial body force are

   .. math::

      P_c = \frac{P_c^*}{p_0}
          = \frac{P_c^*\sigma\Delta x}{\sigma^*},
      \qquad
      G_z^* = \frac{P_c^*}{N_z\Delta x},
      \qquad
      F_z = \frac{P_c}{N_z}.

   Positive values drive fluid A into the sample; negative values reverse the
   direction. Zero gives zero body force. Porosity does not enter this conversion.
   The pressure span specifies the imposed forcing scale, not a guaranteed final
   phase-average capillary pressure. Infer the latter from the simulated phase
   pressures and convert using :math:`p_c^*=p_0(p_A-p_B)`.

.. raw:: html

   <div class="scaling-columns">
   <div class="scaling-column">
   <h4>Input Parameters</h4>
   <h5>Physical properties and lattice scale</h5>
   <div class="scaling-field"><label for="viscA_r3">Fluid A kinematic viscosity <span class="math">\(\nu_A^*\)</span> (m²/s)</label><input id="viscA_r3" aria-label="Fluid A kinematic viscosity" type="number" value="0.00000125" step="any"></div>
   <div class="scaling-field"><label for="viscB_r3">Fluid B kinematic viscosity <span class="math">\(\nu_B^*\)</span> (m²/s)</label><input id="viscB_r3" aria-label="Fluid B kinematic viscosity" type="number" value="0.00000125" step="any"></div>
   <div class="scaling-field"><label for="rhoA_real_r3">Fluid A density <span class="math">\(\rho_A^*\)</span> (kg/m³)</label><input id="rhoA_real_r3" aria-label="Fluid A density" type="number" value="1000" step="any"></div>
   <div class="scaling-field"><label for="rhoB_real_r3">Fluid B density <span class="math">\(\rho_B^*\)</span> (kg/m³)</label><input id="rhoB_real_r3" aria-label="Fluid B density" type="number" value="1000" step="any"></div>
   <div class="scaling-field"><label for="rho_ref_r3">Reference density <span class="math">\(\rho_{\mathrm{ref}}\)</span> (kg/m³)</label><input id="rho_ref_r3" aria-label="Reference density" type="number" value="1000" step="any"></div>
   <div class="scaling-field"><label for="int_ten_r3">Physical interfacial tension <span class="math">\(\sigma^*\)</span> (N/m)</label><input id="int_ten_r3" aria-label="Physical interfacial tension" type="number" value="0.040" step="any"></div>
   <div class="scaling-field"><label for="dx_r3">Voxel resolution <span class="math">\(\Delta x\)</span> (µm)</label><input id="dx_r3" aria-label="Voxel resolution" type="number" value="1.66" step="any"></div>
   <div class="scaling-field"><label for="int_ten_lbm_r3">Lattice interfacial tension <span class="math">\(\sigma\)</span> (lattice units)</label><input id="int_ten_lbm_r3" aria-label="Lattice interfacial tension" type="number" value="0.06" step="any"></div>
   <h5>Protocol and domain parameters</h5>
   <div class="scaling-field"><label for="P_real_r3">Signed equivalent driving pressure span <span class="math">\(P_c^*\)</span> (Pa)</label><input id="P_real_r3" aria-label="Signed equivalent driving pressure span" type="number" value="10000" step="any"></div>
   <div class="scaling-field"><label for="Lx_r3">Simulated x dimension <span class="math">\(N_x\)</span> (voxels)</label><input id="Lx_r3" aria-label="Simulated x dimension" type="number" value="900" step="any"></div>
   <div class="scaling-field"><label for="Ly_r3">Simulated y dimension <span class="math">\(N_y\)</span> (voxels)</label><input id="Ly_r3" aria-label="Simulated y dimension" type="number" value="900" step="any"></div>
   <div class="scaling-field"><label for="Lz_r3">Simulated z length <span class="math">\(N_z\)</span> (voxels)</label><input id="Lz_r3" aria-label="Simulated z length" type="number" value="1600" step="any"></div>
   <div class="scaling-actions"><button type="button" onclick="compute()">Compute</button><button type="button" onclick="clearFields()">Clear</button></div>
   </div>
   <div class="scaling-column">
   <h4>Output Parameters</h4>
   <h5>Shared Color model parameters</h5>
   <div class="scaling-field"><label for="dt_r3">Physical time step <span class="math">\(\Delta t\)</span> (s)</label><input id="dt_r3" aria-label="Physical time step" type="text" readonly></div>
   <div class="scaling-field"><label for="tauA_r3">Fluid A relaxation time <span class="math">\(\tau_A\)</span> (time steps)</label><input id="tauA_r3" aria-label="Fluid A relaxation time" type="text" readonly></div>
   <div class="scaling-field"><label for="tauB_r3">Fluid B relaxation time <span class="math">\(\tau_B\)</span> (time steps)</label><input id="tauB_r3" aria-label="Fluid B relaxation time" type="text" readonly></div>
   <div class="scaling-field"><label for="alpha_r3">Color model alpha <span class="math">\(\alpha\)</span> (lattice units)</label><input id="alpha_r3" aria-label="Color model alpha" type="text" readonly></div>
   <div class="scaling-field"><label for="rhoAlb_r3">Fluid A density <span class="math">\(\rho_A\)</span> (lattice units)</label><input id="rhoAlb_r3" aria-label="Fluid A density" type="text" readonly></div>
   <div class="scaling-field"><label for="rhoBlb_r3">Fluid B density <span class="math">\(\rho_B\)</span> (lattice units)</label><input id="rhoBlb_r3" aria-label="Fluid B density" type="text" readonly></div>
   <div class="scaling-field"><label for="viscAlb_r3">Fluid A kinematic viscosity <span class="math">\(\nu_A\)</span> (lu²/ts)</label><input id="viscAlb_r3" aria-label="Fluid A kinematic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="viscBlb_r3">Fluid B kinematic viscosity <span class="math">\(\nu_B\)</span> (lu²/ts)</label><input id="viscBlb_r3" aria-label="Fluid B kinematic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="muAlb_r3">Fluid A dynamic viscosity <span class="math">\(\mu_A\)</span> (lattice units)</label><input id="muAlb_r3" aria-label="Fluid A dynamic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="muBlb_r3">Fluid B dynamic viscosity <span class="math">\(\mu_B\)</span> (lattice units)</label><input id="muBlb_r3" aria-label="Fluid B dynamic viscosity" type="text" readonly></div>
   <div class="scaling-field"><label for="pressure_scale_r3">Physical pressure per lattice unit <span class="math">\(p_0\)</span> (Pa)</label><input id="pressure_scale_r3" aria-label="Physical pressure per lattice unit" type="text" readonly></div>
   <h5>Protocol parameters</h5>
   <div class="scaling-field"><label for="Pc_lb_r3">Equivalent lattice pressure span <span class="math">\(P_c\)</span> (lattice units)</label><input id="Pc_lb_r3" aria-label="Equivalent lattice pressure span" type="text" readonly></div>
   <div class="scaling-field"><label for="gradient_real_r3">Physical driving pressure gradient <span class="math">\(G_z^*\)</span> (Pa/m)</label><input id="gradient_real_r3" aria-label="Physical driving pressure gradient" type="text" readonly></div>
   <div class="scaling-field"><label for="dpdl_r3">Initial body-force density <span class="math">\(F_z\)</span> (lattice units)</label><input id="dpdl_r3" aria-label="Initial body-force density" type="text" readonly></div>
   </div>
   </div>
   <p id="status_r3" class="scaling-status" role="status" aria-live="polite"></p>

.. rubric:: Generated input.db

Click **Compute** to populate this block, then copy the text to ``input.db``.
The image filename, label mapping, wetting settings, run length, and
single-process decomposition are example settings to customize.

.. code-block:: text
   :name: input-db-r3
   :class: generated-input-db

   // Click Compute to generate input.db for this protocol.

.. raw:: html

   </div>
   </div>
   <script>
   function getSelectedRoutine() {
     return document.querySelector('input[name="routine"]:checked').value;
   }

   function changeRoutine() {
     const selected = getSelectedRoutine();
     for (let i = 1; i <= 6; i++) {
       const box = document.getElementById("routine" + i + "Box");
       if (box) box.style.display = selected === "routine" + i ? "block" : "none";
     }
     const box = document.getElementById(selected + "Box");
     if (box && window.MathJax && window.MathJax.typesetPromise) {
       window.MathJax.typesetPromise([box]).catch(err => console.error(err));
     }
   }

   function readScalingInput(stem, suffix, options = {}) {
     const field = document.getElementById(stem + suffix);
     const label = field.getAttribute("aria-label");
     const text = field.value.trim();
     const value = Number(text);
     if (text === "" || !Number.isFinite(value)) {
       throw new Error(label + " must be a finite number.");
     }
     if (!options.signed && !options.zero && value <= 0) {
       throw new Error(label + " must be greater than zero.");
     }
     if (options.zero && value < 0) {
       throw new Error(label + " must not be negative.");
     }
     if (options.integer && !Number.isSafeInteger(value)) {
       throw new Error(label + " must be a positive whole number.");
     }
     if (options.max !== undefined && value > options.max) {
       throw new Error(label + " must not exceed " + options.max + ".");
     }
     return value;
   }

   function sharedColorScaling(suffix) {
     const dx = readScalingInput("dx", suffix) * 1e-6;
     const sigmaSI = readScalingInput("int_ten", suffix);
     const sigma = readScalingInput("int_ten_lbm", suffix);
     const rhoRef = readScalingInput("rho_ref", suffix);
     const dt = Math.sqrt(sigma * rhoRef * dx * dx * dx / sigmaSI);
     const rhoA = readScalingInput("rhoA_real", suffix) / rhoRef;
     const rhoB = readScalingInput("rhoB_real", suffix) / rhoRef;
     const nuA = readScalingInput("viscA", suffix) * dt / (dx * dx);
     const nuB = readScalingInput("viscB", suffix) * dt / (dx * dx);
     const muA = rhoA * nuA;
     const muB = rhoB * nuB;
     const tauA = 0.5 + 3 * nuA;
     const tauB = 0.5 + 3 * nuB;
     const alpha = sigma / 6;
     const pressureScale = rhoRef * dx * dx / (dt * dt);
     if (![dx, dt, rhoA, rhoB, nuA, nuB, muA, muB, pressureScale].every(v => Number.isFinite(v) && v > 0)
         || !Number.isFinite(tauA) || !Number.isFinite(tauB) || tauA <= 0.5 || tauB <= 0.5) {
       throw new Error("The selected scales exceed numerical precision; check the physical inputs and lattice tension.");
     }
     return {dx, sigmaSI, sigma, rhoRef, dt, rhoA, rhoB, nuA, nuB, muA, muB, tauA, tauB, alpha, pressureScale};
   }

   function computeRoutine1(scales, suffix) {
     const Ca = readScalingInput("Ca", suffix);
     const fB = readScalingInput("fB", suffix, {zero: true, max: 1});
     const area = readScalingInput("Lx", suffix, {integer: true}) * readScalingInput("Ly", suffix, {integer: true});
     const gradient = readScalingInput("gradient", suffix, {signed: true});
     if (gradient === 0) throw new Error("The initial driving pressure gradient must be nonzero for fractional flow.");
     const weightedMu = scales.muA * (1 - fB) + scales.muB * fB;
     const flow = Math.sign(gradient) * area * scales.sigma * Ca / weightedMu;
     return {
       // The local fractional-flow controller divides by 5.796 * alpha.
       capillary_number: Ca * 6 / 5.796,
       Qtotal: flow,
       QA: (1 - fB) * flow,
       QB: fB * flow,
       QAmax: area * scales.sigma * Ca / scales.muA,
       QBmax: area * scales.sigma * Ca / scales.muB,
       dpdl: gradient * scales.dt * scales.dt / (scales.rhoRef * scales.dx)
     };
   }

   function computeRoutine2(scales, suffix) {
     const Ca = readScalingInput("Ca", suffix);
     const nx = readScalingInput("Lx", suffix, {integer: true});
     const ny = readScalingInput("Ly", suffix, {integer: true});
     const nz = readScalingInput("Lz", suffix, {integer: true});
     const phi = readScalingInput("phi", suffix, {max: 1});
     const mlups = readScalingInput("mlups_lbm", suffix);
     // Match ScaLBL_ColorModel::Initialize: dynamic viscosity of fluid B.
     const flow = phi * nx * ny * scales.sigma * Ca / scales.muB;
     const poreVolume = phi * nx * ny * nz;
     const steps = Math.ceil(poreVolume / flow);
     if (!Number.isSafeInteger(steps) || steps < 1) {
       throw new Error("The pore-volume step count exceeds numerical precision; check the capillary number and domain size.");
     }
     return {
       Qz_lb: flow,
       Qz_real: flow * scales.dx * scales.dx * scales.dx / scales.dt,
       pv_ft: steps,
       pv_time: steps * scales.dt,
       time_sim_seconds: poreVolume * steps / (mlups * 1e6)
     };
   }

   function computeRoutine3(scales, suffix) {
     const pressure = readScalingInput("P_real", suffix, {signed: true});
     const nz = readScalingInput("Lz", suffix, {integer: true});
     const latticePressure = pressure / scales.pressureScale;
     return {
       Pc_lb: latticePressure,
       gradient_real: pressure / (nz * scales.dx),
       dpdl: latticePressure / nz
     };
   }

   function setInputDbText(routine, text) {
     const block = document.getElementById("input-db-r" + routine.slice(-1));
     block.querySelector("pre").textContent = text;
   }

   function generateInputDb(routine, scales, outputs) {
     const suffix = "_r" + routine.slice(-1);
     const dimensions = ["Lx", "Ly", "Lz"].map(stem =>
       readScalingInput(stem, suffix, {integer: true, max: 2147483645}));
     const voxelLength = readScalingInput("dx", suffix);
     const protocols = {routine1: "fractional flow", routine2: "core flooding", routine3: "centrifuge"};
     const boundary = {routine1: 0, routine2: 4, routine3: 3}[routine];
     const force = routine === "routine2" ? 0 : outputs.dpdl;
     let driving;
     if (routine === "routine1") {
       driving = `    capillary_number = ${outputs.capillary_number}\n`;
     } else if (routine === "routine2") {
       driving = `    flux = ${outputs.Qz_lb}\n` +
         `    // capillary_number = ${readScalingInput("Ca", suffix)}\n` +
         "    // Leave capillary_number commented to preserve the flux above.\n";
     } else {
       driving = "    din = 1.0\n    dout = 1.0\n";
     }
     const flowAdaptor = routine === "routine1"
       ? "    min_steady_timesteps = 100000\n" +
         "    max_steady_timesteps = 1000000\n" +
         "    mass_fraction_factor = 0.006\n" +
         "    fractional_flow_increment = 0.05\n" +
         "    endpoint_threshold = 0.1\n"
       : "";
     return `// Generated from the ${protocols[routine]} scaling calculator.
   // Physical time step: ${scales.dt} s; physical interfacial tension: ${scales.sigmaSI} N/m.
   // Replace the example image filename and customize image labels and wetting.
   Domain {
       Filename = "segmented.raw"
       ReadType = "8bit"
       N = ${dimensions.join(", ")}
       nproc = 1, 1, 1                 // Example: one MPI process
       n = ${dimensions.join(", ")}  // For MPI, keep n[i] * nproc[i] equal to N[i]
       offset = 0, 0, 0
       voxel_length = ${voxelLength}  // Micrometres
       ReadValues = 0, 1, 2            // Solid, fluid A, fluid B
       WriteValues = 0, 1, 2
       InletLayers = 0, 0, 0
       OutletLayers = 0, 0, 0
       BC = ${boundary}
   }
   Color {
       protocol = "${protocols[routine]}"
       timestepMax = 1000000           // Example run length; adjust for your experiment
       alpha = ${scales.alpha}
       rhoA = ${scales.rhoA}
       rhoB = ${scales.rhoB}
       tauA = ${scales.tauA}
       tauB = ${scales.tauB}
       F = 0, 0, ${force}
   ${driving}    WettingConvention = "SCAL"
       ComponentLabels = 0
       ComponentAffinity = 0.0        // Example: neutral solid wetting
       Restart = false
   }
   Analysis {
       analysis_interval = 1000
       subphase_analysis_interval = 5000
       visualization_interval = 100000
       restart_interval = 100000
       restart_file = "Restart"
   }
   Visualization {
       format = "vtk"
       write_silo = false
       save_8bit_raw = true
       save_phase_field = true
       save_pressure = true
       save_velocity = true
   }
   FlowAdaptor {
   ${flowAdaptor}}\n`;
   }

   function formatTimeNice(seconds) {
     if (seconds === 0) return "0 s";
     if (seconds < 1) return seconds.toExponential(3) + " s";
     if (seconds > Number.MAX_SAFE_INTEGER) return seconds.toExponential(3) + " s";
     seconds = Math.round(seconds * 10) / 10;
     const h = Math.floor(seconds / 3600);
     const m = Math.floor((seconds % 3600) / 60);
     const s = (seconds % 60).toFixed(1);
     if (h > 0) return `${h}h ${m}m ${s}s`;
     if (m > 0) return `${m}m ${s}s`;
     return `${s}s`;
   }

   function compute() {
     const routine = getSelectedRoutine();
     const suffix = "_r" + routine.slice(-1);
     const status = document.getElementById("status" + suffix);
     if (!status) return;
     clearFields();
     try {
       const scales = sharedColorScaling(suffix);
       const calculate = {routine1: computeRoutine1, routine2: computeRoutine2, routine3: computeRoutine3}[routine];
       const outputs = {
         dt: scales.dt, tauA: scales.tauA, tauB: scales.tauB, alpha: scales.alpha,
         rhoAlb: scales.rhoA, rhoBlb: scales.rhoB, viscAlb: scales.nuA, viscBlb: scales.nuB,
         muAlb: scales.muA, muBlb: scales.muB, pressure_scale: scales.pressureScale,
         ...calculate(scales, suffix)
       };
       if (!Object.values(outputs).every(Number.isFinite)) {
         throw new Error("The result exceeds numerical precision; check the input magnitudes.");
       }
       const inputDb = generateInputDb(routine, scales, outputs);
       for (const [stem, value] of Object.entries(outputs)) {
         document.getElementById(stem + suffix).value = stem === "pv_ft" ? String(value) : value.toExponential(8);
       }
       if (routine === "routine2") {
         document.getElementById("time_sim" + suffix).value = formatTimeNice(outputs.time_sim_seconds);
       }
       setInputDbText(routine, inputDb);
       const warnings = [];
       if (scales.alpha > 0.01) warnings.push("Lattice tension implies alpha > 0.01; consider a smaller lattice tension.");
       if (scales.tauA < 0.7 || scales.tauA > 1.5 || scales.tauB < 0.7 || scales.tauB > 1.5) {
         warnings.push("A relaxation time is outside the approximate 0.7–1.5 working range; adjust the lattice scale or grid resolution and verify stability.");
       }
       status.className = "scaling-status";
       status.textContent = warnings.length ? "Scaling computed.\n" + warnings.join("\n") : "Scaling computed.";
     } catch (error) {
       status.className = "scaling-status scaling-error";
       status.textContent = error.message;
     }
   }

   function clearFields() {
     const routine = getSelectedRoutine();
     const box = document.getElementById(routine + "Box");
     if (!box) return;
     box.querySelectorAll("input[readonly]").forEach(field => {field.value = "";});
     setInputDbText(routine, "// Click Compute to generate input.db for this protocol.");
     const status = document.getElementById("status_r" + routine.slice(-1));
     status.textContent = "";
     status.className = "scaling-status";
   }
   </script>
