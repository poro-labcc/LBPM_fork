*****************
3D Digital Rocks
*****************

The use of x-ray microtomography imaging in generating digital rocks has significantly advanced the characterization of 
porous properties. For instance, the 3D images obtained offer valuable insights into the complex pore network, grain 
arrangement, fluid flow properties through numerical simulations, and the overall rock morphology.

In this section, we simulate the 3D images of Berea sandstone, sandpack LV60A and Bentheimer to validate the lbpm_permeability_simulator 
routine for determining absolute permeability in 3D digital rocks. The results obtained are then compared with those reported by 
{cite:t}`michels2021` and {cite:t}`mcclure2021lbpm`. Figure \ref{3DR-PoreStruc} illustrates the pore structure of these 3D digital 
rocks. The raw images of the digital rocks used in this section are available at the following source:

- {numref}`BE-LV-BE-2` (a)  Berea sandstone (https://figshare.com/articles/dataset/Berea_Sandstone/1153794/2);

- {numref}`BE-LV-BE-2` (b)  Sandpack LV60A (https://figshare.com/articles/dataset/LV60A_sandpack/1153795/2);

- {numref}`BE-LV-BE-2` (c)  Bentheimer (https://www.digitalrocksportal.org/projects/218).


```{figure} BE-LV-BE-2.png
---
scale: 90%
align: center
name: BE-LV-BE-2
---
Digital Rocks: a)  Berea sandstone, b)  Sandpack LV60A and c) Bentheimer.

------------------------------
Absolute Permeability - Numerical Setup
------------------------------

{numref}`3DPPM` provides a summary of the parameters for the digital porous media that were employed in the simulations. This table includes detailed information such as the physical and image sizes, the resolution of the 3D images (voxel size) and the porosity. The LBPM setup follows the same structure demonstrated in the Code below. However, it has been adapted to the dimensions of each raw data file.

```{table} - Digital porous media parameters. The superscript $*$ and $\star$ indicate the results reported by Michels et al. (2021) and McClure et al. (2021),  respectively.
:name: 3DPPM
:align: center
| Parameter |Berea Sandstone | LV60A Sandpack | Bentheimer |
|:-:|:-:|:-:|:-:|
| Physical size $(\mu m)$ |$1603.5$ | $3000.6$ | $494^{2}\times 2656$ |
| Image size (pixels) | $300^{3}$ | $300^{3}$  | $900^{2}\times 1600$ |
| Voxel size $(\mu m)$ |$5.345$ | $10.002$ | $1.66$ |
| Porosity $(\%)$ | $19.8$ | $36.8$  | $23.55$  |
| Absolute Permeability (Darcy) | $1.55$  | $36.58$ | $2.92$  |

The LBPM setup code, using as example the Berea sandstone, is given by


```c
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
   offset = 0, 0, 0 // offset to read sub-domain
   voxel_length = 5.345    // voxel length (in microns)
   ReadValues = 0, 1, 2   // labels within the original image
   WriteValues = 1, 0, 2  // associated labels to be used by LBPM
   InletLayers = 0, 0, 5 // specify 10 layers along the z-inlet
   OutletLayers = 0, 0, 5 // specify 10 layers along the z-inlet
   BC = 0                 // boundary condition type (0 for periodic)
}
Visualization {
   write_silo = true     // write SILO databases with assigned variables
   save_8bit_raw = true  // write labeled 8-bit binary files with phase assignments
   save_phase_field = true  // save phase field within SILO database
   save_pressure = true    // save pressure field within SILO database
   save_velocity = true    // save velocity field within SILO database
}
```

In all cases, the domain set of " InletLayers = 0, 0, 5" and "OutletLayers = 0, 0, 5" are used to conect inlet and outlet fo the periodic domain.

~~~~~~~~~~~~~~~~~~~~~
Results for ``BC = 0``
~~~~~~~~~~~~~~~~~~~~~

The results derived from this approach are detailed in {numref}`3DR-R` and {numref}`3DPPM`. {numref}`3DR-R` visually illustrates the streamlines 
observed in each case study, effectively demonstrating the fluid flow through connected pores in z-axis. {numref}`3DPPM` provides a comparative 
analysis of the permeability values obtained from our LBPM setup against the numerical results previously published by {cite:t}`michels2021` 
and {cite:t}`mcclure2021lbpm`.

```{figure} 3DR-R.png
---
scale: 80%
align: center
name: 3DR-R
---
Stream lines: a)  Berea sandstone, b)  Sandpack LV60A and c) Bentheimer.
```

The permeability values obtained numerically for Berea sandstone and sandpack LV60A align closely with the reported 
by {cite:t}`michels2021`. Conversely, for Bentheimer sandstone, the permeability value obtained LBPM was significantly 
higher than the value reported by {cite:t}`mcclure2021lbpm`.