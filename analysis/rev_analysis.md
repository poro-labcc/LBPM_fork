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
rev.ReadParams(filename);
rev.SetDomain();
rev.DetRevAnalysis();
rev.StatRevAnalysis();
```

Keep the usual `Domain` settings, including the geometry filename, dimensions,
labels, voxel length, and process grid. Set `Domain.nproc = 1, 1, 1`.

For geometric analysis, use this REV section:

```text
REV {
    Geometrical = true
    Hydraulic = false
    SizeStep = 0.1
    CCInterval = 10
    gamma = 0.1
    StatSamples = 50
}
```

The existing `Geomtrical` spelling is also accepted. If both spellings are
present, `Geomtrical` takes precedence for compatibility. If an `MRT` section is
present, `MRT.REV` controls whether REV analysis runs. Without an `MRT` section,
the presence of the `REV` section enables analysis.

## Saved velocity fields

For hydraulic analysis, set `Hydraulic = true` and supply output from an earlier
simulation:

```text
REV {
    Geometrical = true
    Hydraulic = true
    SimulationParams = "/path/to/simulation/params.db"
    VelocityPath = "/path/to/simulation/results"
    VelocityTimestep = "vis3000"
}
```

`SimulationParams` points to the database used to run the simulation. REV reads
`tau` and `F` from its `MRT` section, or its `BGK` section if `MRT` is absent.
You do not need to repeat these values in the analysis database. Relative
`SimulationParams` paths are resolved against the directory containing the
analysis database; absolute paths are also accepted. For example,
`SimulationParams = "../simulation/params.db"` selects a sibling simulation
directory.

The referenced file supplies the flow parameters only: keep `Domain` and the
REV settings in the analysis database. Its `MRT.REV` value does not enable or
disable the analysis. Without `SimulationParams`, the analysis continues to
read `tau` and `F` from its own `MRT` section, with the usual solver defaults
for omitted values. When `SimulationParams` is present, the referenced flow
parameters take precedence over the analysis database's `MRT` values.

`VelocityPath` is the output directory containing `LBM.visit` and timestep
directories. It defaults to the current working directory. `VelocityTimestep`
is a directory name, quoted as a string. If omitted, the last entry in the
output's timestep list is used.

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

Call `rev.ReadParams(filename)` first. The model must expose `Nx`, `Ny`, `Nz`,
`Dm->voxel_length`, `Distance`, `Velocity_x`, `Velocity_y`, `Velocity_z`, `tau`,
`Fx`, `Fy`, and `Fz`. The arrays use LBPM's `DoubleArray` type and include the
usual one-cell halos. MRT and BGK models meet this interface.

These calls borrow the model's arrays; they do not copy them or load saved
velocities. They also use the model's `tau` and `F`, while `SimulationParams`
supplies the parameters for standalone analysis. The calculations remain
implemented in `rev_analysis.cpp`, while
the template wrappers in the header provide the common data interface.
