###############################################################################
Representative Elementary Volume (REV) Analysis (lbpm_rev_analysis)
###############################################################################

The degree of heterogeneity observed in digital rocks is an important factor in characterizing porous media. Advances 
in the resolution and sample size achievable with X-ray microtomography have made it increasingly feasible to digitally 
investigate larger volumes that capture the intrinsic heterogeneities of porous structures. These developments enable the 
identification of representative elementary volumes (REVs), beyond which the influence of local heterogeneities becomes 
negligible and specific effective properties of the porous medium can be reliably characterized.

The ``lbpm_rev_analysis`` routine provides an automated REV analysis for a range of 
geometrical and hydrodynamic properties, as summarized in :numref:`Table %s <rev-properties>`. By combining deterministic and statistical approaches, 
the routine enables a comprehensive characterization of these properties, providing valuable insights for assessing the 
representativeness of microtomography images and estimating effective sample volumes for subsequent numerical analyses.

.. _rev-properties:

.. list-table:: Properties available for REV analysis.
   :widths: 25 25 25 25
   :header-rows: 1
   :align: center

   * - 
     - Property 1
     - Property 2
     - Property 3
   * - Geometrical
     - Porosity (:math:`\phi`)
     - Tortuosity (:math:`\tau`)
     - Superficial Area (:math:`A_{s}`)
   * - Hydrodynamic 
     - Absolute Permeability (:math:`K_{z}`)
     - Hydraulic Tortuosity (:math:`\tau_{h}`)
     -

****************************
Methodology
****************************

The implemented REV routine combines deterministic and statistical approaches into a **hybrid methodology**, as illustrated 
in the workflow diagram presented in :numref:`rev-workflow`. The deterministic stage provides an initial estimate of the REV 
size at a relatively low computational cost. This estimate is subsequently used to initialize the statistical stage, which 
evaluates the spatial variability of the property over randomly distributed subvolumes. Both deterministic and statistical 
analyses employ a default convergence threshold of :math:`\gamma = 0.1`  for their respective criteria. The final REV is accepted only when the 
statistical convergence conditions are satisfied. For a detailed explanation of the steps illustrated in :numref:`rev-workflow`, 
expand the toggle sections below for the deterministic and statistical analyses.

.. list-table::
   :widths: 50 50
   :class: borderless
   :align: center

   * - .. _rev-workflow:

       .. figure:: ../../../_static/images/rev-workflow.svg
          :width: 70%
          :align: center

          Workflow of the hybrid REV analysis, combining deterministic
          (steps 1–5) and statistical (steps 6–9) convergence assessments.

     - .. figure:: ../../../_static/images/REVHybridAnimation.gif
          :width: 100%
          :align: center

          GIF Animation illustrating the hybrid REV analysis for porosity.

----------------------
Deterministic Analysis
----------------------

.. raw:: html

   <details>
   <summary><strong>Click here to open: Deterministic Analysis</strong></summary>

The **deterministic approach** evaluates a property of interest within a
subvolume centered in the computational domain. The analysis is repeated
for progressively larger subvolumes until the measured property becomes
approximately independent of the subvolume edge length. The smallest
subvolume satisfying the prescribed convergence criteria is identified
as the deterministic REV.

The main advantage of this approach is its relatively low computational
cost, as only one subvolume is evaluated for each edge length. However,
the results may be sensitive to local fluctuations, including maxima,
minima, and inflection points, potentially leading to false convergence.

The deterministic procedure, illustrated in the upper part of the
workflow diagram, consists of the following steps:

#. **Initialize the subvolume size:** The initial edge length is defined
   as twice the average pore diameter,

   .. math::

      \ell_0 = 2\overline{d}_p.

#. **Evaluate the domain-size constraint:** Before each evaluation,
   the subvolume dimensions are checked against the available computational
   domain. If the prescribed maximum size is exceeded before convergence,
   no deterministic REV is identified.

#. **Measure the property:** The property of interest, :math:`\eta_i`,
   is evaluated within the centered subvolume.

#. **Compute the convergence criteria:** Two complementary metrics are
   calculated to quantify the variation of the measured property with
   increasing subvolume size.

   The first metric is the relative error (RE), defined following
   :cite:t:`li2009permeability` as

   .. math::
      :label: rev-relative-error

      RE_i =
      \left|
      \frac{\eta_{i+1}-\eta_{i-1}}
      {\eta_{i+1}+\eta_{i-1}}
      \frac{V_{i-1}}{V_{i+1}}
      \right| < \gamma,

   where :math:`RE_i` is the relative error at iteration :math:`i`,
   :math:`\eta_i` represents the measured property in the :math:`i`-th
   subvolume, and :math:`V_i` denotes its corresponding volume.
   This metric quantifies the local variation of the property as the
   evaluated subvolume increases in size.

   The second metric is the convergence criterion (CC), introduced by
   :cite:t:`Marafini2020`:

   .. math::
      :label: rev-convergence-criterion

      CC_i =
      \frac{
      \sigma(\eta_{i-m+1},\ldots,\eta_i)
      }{
      \overline{\eta}(\eta_{i-m+1},\ldots,\eta_i)
      } < \gamma,

   where :math:`\sigma` denotes the standard deviation,
   :math:`\overline{\eta}` represents the arithmetic mean, and
   :math:`m` is the number of consecutive measurements considered.
   By evaluating the dispersion of the most recent measurements,
   this metric reduces the likelihood of false convergence associated
   with local fluctuations. The default value is :math:`m=10`.

#. **Assess convergence and update the subvolume size:** The deterministic
   REV is identified when both convergence criteria satisfy the prescribed
   tolerance. Otherwise, the subvolume edge length is increased by 6%,
   according to

   .. math::
      :label: rev-edge-growth

      \ell_{i+1}=1.06\ell_i,

   and the procedure is repeated until convergence is achieved or
   the maximum admissible subvolume size is reached.

.. raw:: html

   </details>

.. raw:: html

   <br>

--------------------
Statistical Analysis
--------------------

.. raw:: html

   <details>
   <summary><strong>Click here to open: Statistical Analysis</strong></summary>

The **statistical approach** evaluates the spatial variability of a
property by sampling multiple subvolumes distributed throughout the
computational domain. For a given edge length, the property is measured
in a set of randomly selected subvolumes, and the dispersion of the
resulting measurements is assessed.

If the dispersion exceeds the prescribed tolerance, the subvolume size
is increased and the procedure is repeated until the measurements become
sufficiently consistent.

Compared with the deterministic approach, statistical analysis is
computationally more expensive because multiple subvolumes must be
evaluated at each edge length. Nevertheless, it provides a more robust
assessment of spatial heterogeneity and reduces the risk of identifying
false convergence caused by localized features.

An important limitation arises as the subvolume size approaches the
dimensions of the computational domain. In this situation, increasing
overlap between sampled subvolumes can artificially reduce the measured
dispersion, potentially leading to an overestimation of convergence.

The statistical procedure, represented in the lower part of the workflow
diagram, follows these steps:

6. **Initialize the statistical analysis:** The edge length corresponding
   to the deterministic REV is adopted as the initial subvolume size,
   avoiding the repeated evaluation of smaller volumes that have already
   been assessed.

7. **Select and evaluate random subvolumes:** For the prescribed edge
   length, :math:`N` subvolumes are randomly distributed within the
   computational domain, and the property of interest is evaluated
   independently in each subvolume.

8. **Evaluate statistical convergence:** The dispersion of the measurements
   is quantified using the coefficient of variation (CV), following
   :cite:t:`Zhang2000`:

   .. math::
      :label: rev-coefficient-variation

      CV = \frac{\sigma(\eta)}{\overline{\eta}} < \gamma,

   where :math:`\sigma(\eta)` is the standard deviation of the property
   measurements and :math:`\overline{\eta}` is their arithmetic mean.
   Normalization by the mean provides a dimensionless measure of relative
   dispersion.

   The number of sampled subvolumes is determined adaptively by monitoring
   the variation of the CV as additional subvolumes are included.
   Sampling continues until the estimated dispersion stabilizes, ensuring
   that the statistical measure is sufficiently representative of the
   selected subvolume size.

9. **Determine the statistical REV:** Once the CV has stabilized, its
   value is compared with the prescribed convergence tolerance.
   If the criterion is satisfied, the current edge length is identified
   as the statistical REV. Otherwise, the subvolume edge length is
   increased by 6%, and the statistical evaluation is repeated.
   If the admissible domain size is exceeded before convergence,
   no statistical REV is identified.

.. raw:: html

   </details>

.. raw:: html

   <br>

************************************
Run REV Analysis
************************************

A typical command to launch the LBPM single-phase simulator is as follows. The REV analysis can be run 
in two different routines: on the ``lbpm_rev_analysis`` and inside of ``lbpm_permeability_simulator``.

-------------------------------
lbpm_permeability_simulator
-------------------------------

Inside of ``lbpm_permeability_simulator``, the proceddure follow the typical command to launch the LBPM 
single-phase simulator is as follows. However, inside ``input.db`` you have activate REV inside of MRT 
class by set ``REV = True``. With REV analysis activated the REV class of parameters need to be set, where 
we have:

.. code-block:: c

   REV {
      Geomtrical = True
      Hydraulic = True
      gamma = 0.1
      SizeStep = 0.06
      CCInterval = 10
      StatSamples = 50
   }

where ``Geomtrical`` and ``Hidraulic`` equal Ture means that REV analysis for each list of geometrical and hydraulic
properties will be done and ``False`` they will note be performed. The gamma 


By default, when ``REV = True`` inside of ``lbpm_permeability_simulator`` the 
``Geomtrical`` and ``Hidraulic`` will be equal True beeing not necessary activate

------------------------------
Model parameters
------------------------------

.. code-block:: text

   REV {
      Geometrical = True
      Hydraulic = True
      gamma = 0.1
      SizeStep = 0.06
      CCInterval = 10
      StatSamples = 50
   }

``Geometrical`` and ``Hydraulic`` are boolean settings read into the REV class.
Both default to ``True`` when omitted and accept ``True`` or ``False``.

--------------------------------------
Assessing Steady State Permeability
--------------------------------------



****************************
Benchmark Cases
****************************


