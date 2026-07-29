*****************
3D Digital Rocks
*****************

The use of X-ray microtomography to generate digital rock images has significantly advanced the characterization of 
porous materials. Three-dimensional images provide valuable insights into complex pore networks, grain arrangements, 
and fluid-flow behavior, which can be investigated through numerical simulations.

In this benchmark, three-dimensional digital rock images available in the literature are used with the 
``lbpm_permeability_simulator`` to estimate absolute permeability. :numref:`Figs. %s <3d-berea>` - :numref:`%s <3d-bentheimer>` illustrates the pore 
structures of some digital rocks considered in this section. The datasets are available from the following sources:

- **Berea sandstone** (a): https://figshare.com/articles/dataset/Berea_Sandstone/1153794/2

- **LV60A sandpack** (b): https://figshare.com/articles/dataset/LV60A_sandpack/1153795/2

- **Bentheimer sandstone** (c): https://www.digitalrocksportal.org/projects/218


.. list-table::
   :widths: 50 50 50
   :align: center

   * - .. figure:: ../../../_static/images/3d-berea.png
          :width: 60%
          :align: center
          :name: 3d-berea

          Berea sandstone (a) :cite:p:`dong2008berea`.

     - .. figure:: ../../../_static/images/3d-lv60a.png
          :width: 60%
          :align: center
          :name: 3d-lv60a

          LV60A sandpack (b) :cite:p:`dong2008lv60a`.

     - .. figure:: ../../../_static/images/3d-bentheimer.png
          :width: 65%
          :align: center
          :name: 3d-bentheimer

          Bentheimer sandstone (c) :cite:p:`dalton2019bentheimer`.

:numref:`Tab. %s <3d-parameters-samples>` summarizes the parameters of the digital porous media analyzed, including 
the image sizes, the resolution of the 3D images (voxel size), as well as the corresponding porosity values.


.. list-table:: Properties of the 3D digital-rock samples.
   :header-rows: 1
   :align: center
   :name: 3d-parameters-samples
   :widths: 20 35 25 20

   * - **Rock type**
     - **Size** (:math:`N_x \times N_y \times N_z`)
     - **Resolution** (:math:`\mu\text{m}`)
     - **Porosity**
   * - Berea (a)
     - :math:`400 \times 400 \times 400`
     - 5.345
     - 0.198
   * - LV60A (b)
     - :math:`450 \times 450 \times 450`
     - 10.002
     - 36.8
   * - Bentheimer (c)
     - :math:`900 \times 900 \times 1600`
     - 1.66
     - 23.55

-------------------------------------------
Absolute Permeability - Numerical Setup
-------------------------------------------

The digital rocks considered in the present analysis are simulated using both the symmetry and mixing-layer schemes to compare 
their effects on the results. As discussed in section :doc:`Boundary Conditions <../../../userGuide/models/mrt/mrt>`, symmetric images 
provide more accurate results in periodic domains (``BC=0``) than those obtained using the mixing-layer scheme. Therefore, the 
present analysis aims to assess the influence of the mixing-layer scheme on the computed absolute permeability.


The LBPM setup for the ``.db`` file follows the structure demonstrated in the toogles below, the first one for mixing layers scheme and
the second one for the symmetric image. The symmetric porous domain needs be created outside of LBPM framework (mirrored in z-axis) and 
inputed as ``.raw`` file in a format of ``uint8`` or ``int8`` or ``uint16`` or ``int16``.

.. raw:: html

   <details>
   <summary><strong>Click here to open: .db example with mixing layers scheme</strong></summary>

.. code-block:: c

   MRT {
      tau = 1.0
      F = 0.0, 0.0, 1.0e-6
      timestepMax = 60000
      tolerance = 0.001
   }
   Domain {
      Filename = "Berea.raw"
      ReadType = "8bit"      // data type
      N = 400, 400, 400     // size of original image
      nproc = 2, 2, 2        // process grid
      n = 200, 200, 200      // sub-domain size
      voxel_length = 5.345    // voxel length (in microns)
      ReadValues = 0, 1, 2   // labels within the original image
      WriteValues = 1, 0, 2  // associated labels to be used by LBPM
      InletLayers = 0, 0, 5 // specify 10 layers along the z-inlet
      OutletLayers = 0, 0, 5 // specify 10 layers along the z-inlet
      BC = 0                 // boundary condition type (0 for periodic)
   }
   Visualization {
      write_silo = true     // write SILO databases with assigned variables
      save_pressure = true    // save pressure field within SILO database
      save_velocity = true    // save velocity field within SILO database
   }

.. raw:: html

   </details>

.. raw:: html

   <br>


To illustrate the aplication of  mixing layers in 3D Digitla Rokcs,  :numref:`Figures %s <berea-mix>` - 
:numref:`%s <bentheimer-mix>`  show the scheme where the superposition 
of the fluid nodes present in the original inlet and outlet planes are applied. Specifically, if a position :math:`(y, z)` 
corresponds to a fluid node at either the inlet or the outlet, the same position in the added layers 
is defined as a fluid node. Conversely, if the position :math:`(y, z)` corresponds to a solid node in both the inlet 
and the outlet planes, it is assigned as a solid node in the added layers. The resulting mixed layers are 
illustrated in the :numref:`Figures %s <berea-mix-res>` - 
:numref:`%s <bentheimer-mix-res>`. This procedure allows the inlet and outlet to be 
connected without the need to reflect the domain sample.

.. figure:: ../../../_static/images/label-color.png
   :width: 84%
   :align: center

.. list-table::
   :widths: 50 50 50
   :align: center

   * - .. figure:: ../../../_static/images/berea-mix.png
          :width: 65%
          :align: center
          :name: berea-mix

          Berea mixing layers.

     - .. figure:: ../../../_static/images/lv60a-mix.png
          :width: 65%
          :align: center
          :name: lv60a-mix

          LV60A mixing layers.

     - .. figure:: ../../../_static/images/bentheimer-mix.png
          :width: 65%
          :align: center
          :name: bentheimer-mix

          Bentheimer mixing layers.

.. figure:: ../../../_static/images/label-color-2.png
   :width: 40%
   :align: center

.. list-table::
   :widths: 50 50 50
   :align: center

   * - .. figure:: ../../../_static/images/berea-mix-res.png
          :width: 65%
          :align: center
          :name: berea-mix-res

          Berea mixed-layer.

     - .. figure:: ../../../_static/images/lv60a-mix-res.png
          :width: 65%
          :align: center
          :name: lv60a-mix-res

          LV60A mixed-layer.

     - .. figure:: ../../../_static/images/bentheimer-mix-res.png
          :width: 65%
          :align: center
          :name: bentheimer-mix-res

          Bentheimer mixed-layer.

:numref:`Figure %s <berea-symmetry>` illustrates the symmetry applied along the :math:`z`-axis, in comparison 
with the original image shown in :numref:`Figure %s <3d-berea-2>`. As indicated in the toggle below, the image size in the input ``.db`` file 
is doubled, which consequently doubles the computational cost.


.. raw:: html

   <details>
   <summary><strong>Click here to open: .db example for symmetric image scheme</strong></summary>

.. code-block:: c

   MRT {
      tau = 1.0
      F = 0.0, 0.0, 1.0e-6
      timestepMax = 60000
      tolerance = 0.001
   }
   Domain {
      Filename = "Berea.raw"
      ReadType = "8bit"      // data type
      N = 400, 400, 400     // size of original image
      nproc = 2, 2, 2        // process grid
      n = 200, 200, 200      // sub-domain size
      voxel_length = 5.345    // voxel length (in microns)
      ReadValues = 0, 1, 2   // labels within the original image
      WriteValues = 1, 0, 2  // associated labels to be used by LBPM
      InletLayers = 0, 0, 5 // specify 10 layers along the z-inlet
      OutletLayers = 0, 0, 5 // specify 10 layers along the z-inlet
      BC = 0                 // boundary condition type (0 for periodic)
   }
   Visualization {
      write_silo = true     // write SILO databases with assigned variables
      save_pressure = true    // save pressure field within SILO database
      save_velocity = true    // save velocity field within SILO database
   }

.. raw:: html

   </details>

.. raw:: html

   <br>


.. list-table::
   :widths: 50 50
   :align: center

   * - .. figure:: ../../../_static/images/3d-berea.png
          :width: 45%
          :align: center
          :name: 3d-berea-2

          Berea sandstone.

     - .. figure:: ../../../_static/images/berea-symmetry.png
          :width: 60%
          :align: center
          :name: berea-symmetry

          Berea sandstone symmetric in z-axis.

~~~~~~~~~~~~~~~~~~~~~~~~~
Results for ``BC = 0``
~~~~~~~~~~~~~~~~~~~~~~~~~

Table XX presents the results obtained for the samples listed in Table YY. The average percentage deviation between the symmetry 
and mixing-layer schemes is approximately ZZ%, indicating good agreement between the two approaches. Additionally, isotropy 
analyses were performed for samples XX, YY, and ZZ, the corresponding mean values and standard deviations are reported based on
values absolute permeability observed for each direction.


.. raw:: html

   <table id="abs-perm-results"
          class="docutils align-center"
          style="border-collapse: collapse; width: 60%;">

     <caption>
       Absolute-permeability results obtained using the symmetric and mixing-layer schemes.
     </caption>

     <thead>
       <tr>
         <th rowspan="2"
             style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: left; vertical-align: middle;">
           Rock type
         </th>

         <th colspan="3"
             style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           \(K_{zz}\,[\mathrm{Da}]\)
         </th>

         <th colspan="3"
             style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           \(K_{yy}\,[\mathrm{Da}]\)
         </th>

         <th colspan="3"
             style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           \(K_{xx}\,[\mathrm{Da}]\)
         </th>
       </tr>

       <tr>
         <th style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           Symmetric
         </th>
         <th style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           Mixing layer
         </th>
         <th style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           \(\Delta_{\%}\)
         </th>

         <th style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           Symmetric
         </th>
         <th style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           Mixing layer
         </th>
         <th style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           \(\Delta_{\%}\)
         </th>

         <th style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           Symmetric
         </th>
         <th style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           Mixing layer
         </th>
         <th style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;">
           \(\Delta_{\%}\)
         </th>
       </tr>
     </thead>

     <tbody>
       <tr>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: left;">
           Berea (a)
         </td>

         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>

         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>

         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
       </tr>

       <tr style="background-color: #f3f4f4;">
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: left;">
           LV60A (b)
         </td>

         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>

         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>

         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
       </tr>

       <tr>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: left;">
           Bentheimer (c)
         </td>

         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>

         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>

         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
         <td style="border: 1px solid #d8d8d8; padding: 8px 12px; text-align: center;"></td>
       </tr>
     </tbody>
   </table>