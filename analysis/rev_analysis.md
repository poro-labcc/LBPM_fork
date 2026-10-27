# REV analysis

`ScaLBL_REV` supports standalone analysis and analysis of existing model data.

## Standalone executable

Build the `lbpm_rev_analysis` target, then run it with one MPI process:

```sh
mpirun -np 1 /path/to/build/tests/lbpm_rev_analysis params.db
```

The executable reads the geometry configured in `Domain`, computes its signed
distance field, and runs deterministic and statistical REV analysis. Optional
permeability subdomain selection also runs flow simulations. The API sequence is:

```cpp
ScaLBL_REV rev(rank, nprocs, comm);
rev.ReadParams(filename, false);
rev.SetDomain();
rev.DetRevAnalysis();
rev.StatRevAnalysis();
rev.sub_domain_select(); // runs only when REV.SubdomainSelect = true
```

Keep the usual `Domain` settings, including the geometry filename, dimensions,
labels, voxel length, and process grid. Set `Domain.nproc = 1, 1, 1`.

The standalone executable defaults to `Hydraulic = false`, so geometry-only
analysis requires no saved velocity output. Explicit settings in `REV` override
this default. For geometric analysis, use this REV section:

```text
REV {
    Geometrical = true
    Hydraulic = false
    SizeStep = 0.1
    CCInterval = 10
    gamma = 0.1
    StatSamples = 50
    Seed_Anchor = true
}
```

Use `Geometrical` for the geometrical setting. The presence of the `REV` section
enables analysis, independently of any `MRT.REV` setting.

`Seed_Anchor` controls random subvolume placement in the statistical analysis
and defaults to `true` in both REV executables. With `Seed_Anchor = true`, the
random generator uses the fixed seed `0`, making sampling reproducible for the
same input and configuration. With `Seed_Anchor = false`, it uses a seed from
`std::random_device`, so sampled locations can vary between runs. This setting
does not affect the deterministic analysis.

## Saved velocity fields

For hydraulic analysis, set `Hydraulic = true` and supply output from an earlier
simulation. The standalone executable requires explicit, nonempty values for
`SimulationParams`, `VelocityPath`, and `VelocityTimestep`. If any are missing or
blank, it lists the missing settings and exits with a nonzero status before
reading saved inputs or geometry. Geometry-only analysis does not require them.

```text
REV {
    Geometrical = true
    Hydraulic = true
    Seed_Anchor = true
    SimulationParams = "/path/to/simulation/params.db"
    VelocityPath = "/path/to/simulation/results"
    VelocityTimestep = "vis3000"
}
```

`SimulationParams` points to the database used to run the simulation. The
referenced database must contain an `MRT` section, from which REV reads `tau`
and `F`. A `BGK` section is not accepted as a substitute. You do not need to
repeat these values in the analysis database. Relative `SimulationParams`
paths are resolved against the directory containing the analysis database;
absolute paths are also accepted. For example,
`SimulationParams = "../simulation/params.db"` selects a sibling simulation
directory.

The referenced file supplies the flow parameters only: keep `Domain` and the
REV settings in the analysis database. Its `MRT.REV` value does not enable or
disable the analysis. The referenced flow parameters take precedence over the
analysis database's `MRT` values. Missing `tau` or `F` entries in that referenced
`MRT` section use the usual solver defaults; the standalone executable requires
the `SimulationParams` path itself whenever hydraulic analysis is enabled.

`VelocityPath` is the output directory containing timestep directories.
Relative paths are resolved from the current working directory; use `"."` to
select it explicitly. `VelocityTimestep` is a directory name, quoted as a
string. Both settings must be supplied for standalone hydraulic analysis; the
executable does not select a default directory or timestep.

The selected timestep must contain `LBM.summary` and a `domain` mesh with scalar
cell fields named `Velocity_x`, `Velocity_y`, and `Velocity_z`. LBPM's existing
IO reader supports Silo and HDF5 when the corresponding libraries are enabled
in the build. The saved output must have one domain with dimensions matching
the input geometry. Use the same geometry, voxel length, relaxation time, and
forcing as the simulation that generated the velocity fields.

Setting `Hydraulic = false` requires no saved simulation output. Disabled
geometric or hydraulic output is not opened or written; existing CSV files
are retained.

## Representative subdomain selection

Set `SubdomainSelect = true` in `REV` to select representative subdomains after
the statistical analysis. It defaults to `false` in both REV executables. For
example, this selects a `64 x 64 x 64` geometric subdomain without hydraulic
input, using the statistical REV properties as the selection targets:

```text
REV {
    Geometrical = true
    Hydraulic = false
    SubdomainSelect = true
    subdomain_size = 64
    save_subdomain_raw = false
    gamma = 0.1
    StatSamples = 50
    Seed_Anchor = true
}
```

`save_subdomain_raw` controls whether selection writes the selected `.raw` crops
and defaults to `false`. Set it to `true` in `input.db` to enable raw exports for
both geometric and hydraulic selection. With raw export disabled, the selection
CSV and domain databases are still written, and any existing raw crops or
temporary raw files are left untouched. The CSV's `raw_file` field is blank.

Enable `Hydraulic` and supply the saved-flow settings described above to also
select a permeability subdomain. The statistical analysis retains the first
size at which a property's coefficient of variation satisfies
`0 <= CV < gamma`, together with its mean and dimensions. Porosity and
permeability retain independent results. Selection uses these results directly;
it does not repeat the REV analysis. The mean specific surface area is also
retained at the porosity REV size for geometric scoring, and mean hydraulic
tortuosity is retained at the permeability REV size for hydraulic scoring.
Each target pair therefore comes from the same statistical window dimensions,
even if the surface-area or tortuosity REV would have a different size.

`subdomain_size` optionally overrides the selection window dimensions in
voxels. A single positive integer, such as `subdomain_size = 64`, selects a
`64 x 64 x 64` cube. Three positive integers, such as
`subdomain_size = 64, 48, 32`, specify the exact x, y, and z dimensions. The
configured dimensions apply to both geometric and hydraulic selection and
must fit within the analyzed `Domain.n`. If omitted, each property uses its
own statistical REV dimensions: cubic for a cubic image, or proportionally
scaled and rounded as in the statistical analysis for a rectangular image.

The override changes only the selection window. Statistical REV sizes and
all retained target means, including specific surface area and tortuosity,
remain unchanged. Selection compares candidates of the requested size against
those targets; it does not repeat REV analysis at the custom size. A converged
statistical REV and valid target properties are still required. A custom size
does not provide missing reference values or enable an unconverged property.

Every valid window is evaluated with a one-voxel increment in each direction,
including the last position at each boundary. Porosity uses the analysis pore
convention, `Distance > 0`, and summed voxel counts. Specific surface area
`S_v` uses the existing geometric measure: count each internal solid-pore face
once along x, y, and z, using `S_v = N_internal_faces / (N_voxels * h)`, where
`N_voxels` is the number of voxels in the window and `h` is the voxel length.
Its units are inverse micrometres. Window exterior faces and
periodic wraparound faces are excluded. Permeability and hydraulic tortuosity
use the same fresh MRT simulation for each candidate, with periodic boundaries,
the same `tau` and force vector, and the existing Darcy conversion to
millidarcies for permeability. Statistical hydraulic means still come from
the previously computed full-domain velocity fields; the candidate values
come from simulations of the cropped geometry.
These boundary conditions can produce differences from the statistical mean.

The geometric selection minimizes a fixed weighted score, with lower values
indicating a closer match to the statistical targets:

```text
E_phi = abs(phi_local - mean_porosity) / abs(mean_porosity)
E_surface = abs(S_v_local - mean_surface) / abs(mean_surface)
selection_score = 0.6 * E_phi + 0.4 * E_surface
```

The hydraulic selection similarly minimizes:

```text
E_k = abs(k_local - mean_permeability) / abs(mean_permeability)
E_tortuosity = abs(T_local - mean_tortuosity) / abs(mean_tortuosity)
selection_score = 0.6 * E_k + 0.4 * E_tortuosity
```

The weights apply automatically: 60% porosity and 40% specific surface for the
geometric region, and 60% permeability and 40% tortuosity for the hydraulic
region. These fixed weights apply to candidates at the actual selection
window dimensions, including a configured `subdomain_size`.
For each component, a zero statistical mean gives zero deviation for an exact
zero and infinite deviation otherwise. Equal
scores select the first candidate in traversal order: x changes fastest,
followed by y and z. The traversal remains exhaustive after an exact match.

Hydraulic tortuosity follows the existing displayed analysis convention,
`T = sum_pores(|u|) / sum_pores(u_z)`, using the converged numerical velocities.
The statistical routine internally stores `T - 1`; the retained selection
target is `T`, matching `stat_tort_mean`. Its denominator uses z velocity,
rather than the force-projected velocity used for permeability, so use flow
along positive z for the usual positive tortuosity interpretation. A candidate
with no flow or undefined net z transport has undefined tortuosity and cannot
win the hydraulic score. Every window is still evaluated and counted. If no
candidate has usable tortuosity, the hydraulic row has status
`no_valid_tortuosity`. If the statistical tortuosity target is unavailable at
the permeability REV size, its status is `unavailable_tortuosity`. Neither case
produces a hydraulic raw crop or domain database; tortuosity is never replaced
with a geometric estimate or an arbitrary zero.

Nonzero `Domain.InletLayers` and `Domain.OutletLayers` are supported. In this
case selection rereads the original source geometry to build the candidate
pore mask before synthetic layer modifications. Raw labels use the existing
`ReadValues`/`WriteValues` convention, including the signed-byte narrowing
performed by the current `16bit` geometry reader. This classification does not
change the original bytes exported below. `ID.00000` labels are already mapped
and are used directly.

The layer settings, full-domain analysis, saved velocities, and retained
statistical REV sizes and means remain unchanged. Those statistical targets
can include synthetic layer effects, so deviations from candidates measured
on the original rock can be larger.

Selection requires one MPI process. Permeability selection additionally
requires periodic flow (`BC = 0`) and a nonzero body force. It inherits
`timestepMax`, `analysis_interval`, and
`tolerance` from the referenced `SimulationParams` MRT block, or the local MRT
block when no reference is supplied. Their defaults are `100000`, `1000`, and
`1e-8`, respectively. A candidate that fails to converge stops selection with
an error. There are `(Nx-nx+1)*(Ny-ny+1)*(Nz-nz+1)` candidates, where
`Nx,Ny,Nz` are the analyzed image dimensions from `Domain.n` and `nx,ny,nz`
are the actual selection window dimensions. Exhaustive permeability selection
runs a simulation for each candidate and can be substantially more expensive
than statistical REV analysis.

Results are written to the current working directory:

- When `save_subdomain_raw = true`, `subdomain_porosity_<nx>x<ny>x<nz>.raw` and
  `subdomain_permeability_<nx>x<ny>x<nz>.raw` contain the selected windows with
  x varying fastest and no halos. Here `nx`, `ny`, and `nz` are the selected
  window dimensions in voxels along x, y, and z, respectively. For example,
  a `64 x 64 x 64` geometric window is saved as
  `subdomain_porosity_64x64x64.raw`. Each property's filename uses its own
  selected dimensions. For `Domain.Filename` with
  `ReadType = "8bit"` or `"16bit"`, the original bytes are copied without
  relabeling or data-type conversion. `ID.00000` input is exported as its
  original one-byte labels, with halos removed, and requires `BC = 0`.
  Selection does not support `GridFile` or `swc` input.
- `subdomain_selection.csv` is comma-separated and records each property's
  status, REV size in voxels and micrometres, window dimensions and physical
  lengths, statistical mean, selected value, relative deviation, candidate
  count, raw filename in `raw_file` (blank when raw export is disabled), and
  source layout. The raw filename includes the selected window dimensions.
  `rev_size_vx` and `rev_size_um` retain the statistical reference size, even
  when `subdomain_size` is set. `nx,ny,nz`, physical lengths, coordinates,
  candidate count, and the raw filename describe the actual selection window.
  Coordinates are zero-based with exclusive upper bounds. `x0,y0,z0` and
  `x1,y1,z1` refer to the halo-free
  analysis image; the `source_*` coordinates additionally account for
  `Domain.offset` (and input halos for `ID.00000`), so the crop can be reproduced
  from the source file. Source dimensions and bytes per voxel are also recorded.
  The geometric row keeps `statistical_mean`, `selected_value`, and
  `relative_deviation` for porosity. It also records `statistical_mean_surface`,
  `selected_surface`, `relative_deviation_surface`, `porosity_weight`, and
  `surface_weight`; these fields are blank for the hydraulic row. The hydraulic
  row keeps `statistical_mean`, `selected_value`, and `relative_deviation` for
  permeability, and records `statistical_mean_tortuosity`, `selected_tortuosity`,
  `relative_deviation_tortuosity`, `permeability_weight`, and `tortuosity_weight`;
  these fields are blank for the geometric row. `selection_score` records the
  applicable weighted score for either selected region.
- `geometrical-rev-domain.db` locates the weighted geometric selection, and
  `hydraulic-rev-domain.db` locates the weighted hydraulic selection. Each file
  contains only a `Domain` block referencing the original full sample: `Filename` is the
  absolute source path, `N` is its full size, `n` is the actual selected window
  size (including any `subdomain_size` override), and `offset` is its starting
  position in the source image, including any
  original analysis offset. The CSV's `domain_file` column identifies this file.

The geometric row retains the property name `porosity`, and its optional raw
crop uses `subdomain_porosity_<nx>x<ny>x<nz>.raw`; both identify the same region
selected using the combined porosity and specific-surface score.
The hydraulic row retains the property name `permeability`, and
the optional `subdomain_permeability_<nx>x<ny>x<nz>.raw` crop identifies the
region selected using the combined permeability and tortuosity score.

The generated domain blocks preserve the raw input type, label mapping, and
voxel length. They use one MPI process, `BC = 0`, and zero inlet/outlet layers
so loading them reproduces the selected original rock region. Physical `L`
lengths, when present, are updated for the selected dimensions. For `ID.00000`
input, the file is referenced as an `8bit` raw image with its full haloed size,
an offset that includes the one-voxel halo, and an identity label mapping.
Load an `ID.00000`-based domain from a separate working directory: the normal
decomposition writes its own `ID.00000` in the working directory.

For example, a selected region within a `400 x 400 x 400` sample could produce:

```text
Domain {
    Filename = "/path/to/Berea.raw"
    ReadType = "8bit"
    ReadValues = 0, 1
    WriteValues = 1, 0
    N = 400, 400, 400
    n = 64, 64, 64
    offset = 100, 120, 80
    nproc = 1, 1, 1
    voxel_length = 5.345
    BC = 0
    InletLayers = 0, 0, 0
    OutletLayers = 0, 0, 0
}
```

Reuse this `Domain` block in a subsequent input database to load the selected
region from the full sample. The generated file has no `MRT`, `REV`, or
`Visualization` blocks; add the other settings needed for the intended solver.

A disabled property or a property with no converged statistical REV is recorded
with a status instead of a selected window; no raw crop or domain database is
produced for it.

## Existing model data

The methods also accept any compatible model through templates:

```cpp
rev.PoreSize(model);
rev.DetRevAnalysis(model);
rev.StatRevAnalysis(model);
rev.sub_domain_select(model);
```

Call `rev.ReadParams(filename)` first; this keeps the hydraulic default enabled
for the combined simulation workflow. The model must expose `Nx`, `Ny`, `Nz`,
`Dm->voxel_length`, `Distance`, `Velocity_x`, `Velocity_y`, `Velocity_z`, `tau`,
`Fx`, `Fy`, and `Fz`. The arrays use LBPM's `DoubleArray` type and include the
usual one-cell halos. MRT and BGK models meet this interface.

These calls borrow the model's arrays; they do not copy them or load saved
velocities. They also use the model's `tau` and `F`, while `SimulationParams`
supplies the parameters for standalone analysis. The calculations remain
implemented in `rev_analysis.cpp`, while
the template wrappers in the header provide the common data interface.
