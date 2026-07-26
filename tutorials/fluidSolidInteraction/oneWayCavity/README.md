---
sort: 3
---

# One-way fluid-solid interaction in a lid-driven cavity: `oneWayCavity`

## Tutorial Aims

- Demonstrate a one-way fluid-solid interaction (FSI) workflow.
- Show how pre-computed fluid pressures and tractions can be transferred to a
  solid without feeding the solid deformation back to the fluid.

## Case Overview

The fluid domain is a two-dimensional square cavity with side length
`0.1 m`. The top wall moves at `1 m/s`, while the other walls are fixed. The
fluid is incompressible, with density `1000 kg/m3` and kinematic viscosity
`1e-3 m2/s`.

The cavity is surrounded on its left, bottom, and right sides by a solid. The
solid is linear-elastic steel with density `8000 kg/m3`, Young's modulus
`200 GPa`, and Poisson's ratio `0.3`. Its outer top edges are fixed. The fluid
loads on the stationary cavity walls are transferred directly to the solid
`interface` patch.

This is a one-way FSI case: the fluid solution is calculated on the undeformed
fluid mesh, and the resulting fluid loading drives the solid deformation. The
solid motion therefore does not alter the fluid solution.

The case uses the `oneWayFsiFluid` fluid model and the `oneWayCoupling`
fluid-solid interface with `directMap` transfer.

## Running the Case

Run the included script from this directory:

```bash
./Allrun
```

The script performs the following steps:

1. Creates the fluid mesh and runs `icoFoam` to completion.
2. Creates the solid mesh.
3. Links the completed fluid fields into the solid case as a `fluid` region.
4. Runs `solids4Foam`, which reads the fluid loading and solves the solid
   response.

The `Allclean` script removes generated meshes, time directories, logs, and
the links created for the fluid results:

```bash
./Allclean
```

## Results

After the run, fluid velocity and pressure fields are available in the main
case time directories. The solid displacement field is written under the
`solid` region. These can be viewed together in ParaView by opening the case
and selecting the required regions.

For a broader description of the one-way coupling implementation, see the
[fluid-solid interfaces README](../../../src/solids4FoamModels/fluidSolidInterfaces/README.md).
