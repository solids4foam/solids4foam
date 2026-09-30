---
sort: 3
---

# Flow over a heated plate: `flowOverHeatedPlate`

## Tutorial Aims

- Demonstrate a conjugate heat transfer (thermal fluid-solid interaction)
  workflow, where the fluid and solid temperature fields are coupled at a
  shared interface.
- Show the use of the `thermal` fluid-solid interface with `thermalRobin`
  boundary conditions on both sides of the interface.

## Case Overview

A fluid flows through a two-dimensional channel over a plate that is heated
from below. The set-up resembles the forced-convection conjugate heat transfer
problem of Vynnycky et al. [1].

The fluid region is a channel of length `3.5 m` ($$-0.5 \le x \le 3$$ m) and
height `0.5 m`. The fluid enters at the `inlet` with a uniform velocity of
`1 m/s` and a temperature of `300 K`. The fluid has a density of `1 kg/m3`, a
kinematic viscosity of `2e-4 m2/s`, a specific heat capacity of `250 J/(kg K)`
and a thermal conductivity of `5 W/(m K)`. The flow is laminar and buoyancy is
neglected (`beta` is zero), so the temperature is transported as a passive
scalar.

The solid region is a plate of length `1 m` and thickness `0.25 m`, located
below the channel at $$0 \le x \le 1$$ m. Its `bottom` surface is held at
`350 K`, its `left` and `right` sides are insulated, and its `top` surface is
the interface with the fluid. The plate has a thermal conductivity of
`100 W/(m K)`, a specific heat capacity of `100 J/(kg K)` and a density of
`1 kg/m3`.

The channel floor upstream of the plate (`slip-bottom`) is a slip surface,
while the plate surface (`interface`) and the floor downstream of the plate
(`bottom`) are no-slip walls. All fluid boundaries apart from the inlet and
the interface are adiabatic.

The case uses:

- the `pimpleFluid` fluid model with `solveEnergyEq true`;
- the `thermalSolid` solid model, which solves only for the temperature in the
  plate, i.e. the plate does not deform;
- the `thermal` fluid-solid interface in `constant/fsiProperties`, which
  couples the solid `top` patch to the fluid `interface` patch.

The simulation is transient, with a time step of `0.01 s` and an end time of
`2 s`. Within each time step, the fluid and solid are solved repeatedly until
the interface residual falls below `outerCorrTolerance`.

## Expected Results

A thermal boundary layer develops in the fluid above the plate, and the plate
cools from its top surface, where heat is conducted into the fluid. The
interface residuals are written to the `residuals` directory. The temperature
fields in both regions can be viewed in ParaView.

```note
This README is a short summary. Reference results and a comparison with the
literature have not yet been added.
```

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/thermoFluidSolidInteraction/flowOverHeatedPlate`. Run
the included script from this directory:

```bash
./Allrun
```

The script creates the solid and fluid meshes with `blockMesh` and then runs
the `solids4Foam` solver. The case can be cleaned with `./Allclean`.

## References

[1]
[M. Vynnycky, S. Kimura, K. Kanev, and I. Pop, "Forced convection heat
transfer from a flat plate: the conjugate problem", International Journal of
Heat and Mass Transfer, 41(1), 45-59,
1998.](https://doi.org/10.1016/S0017-9310(97)00113-0)
