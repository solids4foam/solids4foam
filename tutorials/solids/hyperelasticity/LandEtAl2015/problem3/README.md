---
sort: 10
---

# Inflation and active contraction of an idealised ventricle: `LandEtAl2015/problem3`

## Tutorial Aims

- Demonstrate an anisotropic hyperelastic material with a spatially varying
  fibre direction.
- Demonstrate the `electroMechanicalLaw`, which adds an active fibre tension
  to a passive mechanical law.
- Provide the set-up for problem 3 of the cardiac mechanics benchmark of Land
  et al. [1].

```note
This case requires solids4foam to be built with PETSc, and it only runs with
OpenFOAM.com versions. Otherwise, the `Allrun` script exits without running the
case.
```

## Case Overview

The geometry is the idealised left ventricle of Land et al. [1]: a truncated
ellipsoid, with an endocardial surface of semi-axes `7 mm` and `17 mm`, an
epicardial surface of semi-axes `10 mm` and `20 mm`, and a base plane at
`z = 5 mm`. The same geometry is used, without fibres or active tension, in the
[`idealisedVentricle`](../../idealisedVentricle/README.md) tutorial (problem 2
of the benchmark).

The base (`fixed` patch) is held fixed. A pressure is applied to the
endocardium (`inside` patch), increasing linearly from zero to `15 kPa` over
the `1 s` of the simulation. The epicardium (`outside` patch) is
traction-free.

The material is defined in `constant/mechanicalProperties`:

- The passive behaviour is the transversely isotropic `GuccioneElastic` law,
  with `k = 2 kPa`, `cf = 8`, `ct = 2` and `cfs = 4`.
- The `electroMechanicalLaw` adds an active tension along the fibre
  direction, ramped linearly to `60 kPa` over `rampTime = 1 s`.
- The fibre direction varies through the wall, from `90` degrees at the
  endocardium to `-90` degrees at the epicardium. The fibre field is created
  by the `setFibreField` utility before the solver is run.

The mesh is created in three steps: `blockMesh` creates a slice of the
ventricle wall, `extrudeMesh` rotationally extrudes it through `360` degrees,
and `createPatch` removes the leftover empty patches. Six mesh densities are
provided in `verification/mesh`, from the coarsest (`1`, the default) to the
finest (`6`).

Two solid formulations are provided, both using the total Lagrangian
finite-strain approach and solved with PETSc SNES:

- `displacement` (default): the displacement-based formulation;
- `pressure`: the mixed displacement-pressure formulation.

## Expected Results

The ventricle deforms under the combined action of the applied pressure and
the active tension along the fibres. The displacement of a point in the wall
is written to the `postProcessing` directory by the `solidPointDisplacement`
function object. On the coarsest mesh with the displacement formulation, the
final displacement magnitude at this point is approximately `0.85 mm`.

```note
This README is a short summary. A comparison with the benchmark results of
Land et al. [1] has not yet been added.
```

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/hyperelasticity/LandEtAl2015/problem3`. Run the
included script from this directory:

```bash
./Allrun
```

The script accepts the following options, which can be combined:

```bash
# Use the mixed displacement-pressure formulation
./Allrun pressure

# Run in parallel (see numberOfSubdomains in system/decomposeParDict)
./Allrun parallel

# Use a finer mesh: 1 is the coarsest and 6 is the finest
MESH_NUMBER=2 ./Allrun
```

The case can be cleaned with `./Allclean`.

## References

[1]
[S. Land, V. Gurev, S. Arens, et al., "Verification of cardiac mechanics
software: benchmark problems and solutions for testing active and passive
material behaviour", Proceedings of the Royal Society A, 471(2184), 20150641,
2015.](https://doi.org/10.1098/rspa.2015.0641)
