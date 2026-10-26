# REV analysis

`ScaLBL_REV` supports standalone analysis and analysis of existing model data.

## Standalone executable

Build the `lbpm_rev_analysis` target, then run it with one MPI process:

```sh
mpirun -np 1 /path/to/build/tests/lbpm_rev_analysis params.db
```

The executable reads the geometry configured in `Domain`, computes its signed
distance field, and runs deterministic and statistical REV analysis. It does not
run a flow solver. The API sequence is:

```cpp
ScaLBL_REV rev(rank, nprocs, comm);
rev.ReadParams(filename, false);
rev.SetDomain();
rev.DetRevAnalysis();
rev.StatRevAnalysis();
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

## Existing model data

The methods also accept any compatible model through templates:

```cpp
rev.PoreSize(model);
rev.DetRevAnalysis(model);
rev.StatRevAnalysis(model);
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
