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
MRT {
    REV = true
    tau = 1.0
    F = 0.0, 0.0, 1.0e-5
}
REV {
    Geometrical = true
    Hydraulic = true
    VelocityPath = "/path/to/simulation/results"
    VelocityTimestep = "vis3000"
}
```

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
velocities. The calculations remain implemented in `rev_analysis.cpp`, while
the template wrappers in the header provide the common data interface.
