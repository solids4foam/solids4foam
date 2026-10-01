---
sort: 2
---

# Bi-material cylinder with a prescribed temperature field: `hotCylinderPredefinedTFieldMultipleMaterials`

## Tutorial Aims

- Demonstrate a thermo-elastic analysis where the temperature field is not
  solved for, but is instead read from a separate case.
- Demonstrate the use of multiple materials, where each material is assigned
  to a cell zone.

## Case Overview

This case is a variant of the
[`hotCylinder`](../hotCylinder/README.md) tutorial. The geometry is the same
plane strain quarter model of a thick-walled cylinder with an inner radius of
`0.5 m` and an outer radius of `0.7 m`, discretised with `10` cells in the
radial direction and `60` cells in the circumferential direction. The inner and
outer surfaces are traction-free.

There are two differences from `hotCylinder`.

**Two materials.** The cylinder is split at a radius of `0.6 m` into an inner
steel layer and an outer aluminium layer. The `steel` and `aluminium` cell
zones are created by `setSet` (using `batch.setSet`) and `setsToZones`, and
each entry in `constant/mechanicalProperties` is applied to the cell zone of
the same name:

| Material  | Young's modulus | Poisson's ratio | Thermal expansion |
| --------- | --------------- | --------------- | ----------------- |
| steel     | `200 GPa`       | `0.3`           | `1e-5 1/K`        |
| aluminium | `70 GPa`        | `0.2`           | `2.3e-5 1/K`      |

**Prescribed temperature.** The `linearGeometryTotalDisplacement` solid model
is used, which solves for the displacement only. The `thermoMechanicalLaw`
reads the temperature field at each time step from the time directories of the
case given by its `TcaseDirectory` entry:

```c++
steel
{
    type            thermoMechanicalLaw;
    alpha           alpha [0 0 0 -1 0 0 0] 1e-05;
    T0              T0 [0 0 0 1 0 0 0] 0;
    TcaseDirectory  "hotCylinderTemperatureField";
    mechanicalLaw
    {
        type            linearElastic;
        ...
    }
}
```

The `hotCylinderTemperatureField` directory holds a uniform temperature field
for each time: `50 K` at time `0`, rising by `10 K` per time step to `90 K` at
time `4`. The stress-free reference temperature is `T0 = 0 K`.

## Expected Results

The temperature is uniform at each time. The stresses reflect both the
plane-strain constraint and the mismatch in thermal expansion between the
two materials: the aluminium expands more than the steel. Even a
single-material cylinder can develop thermal stress under the plane-strain
constraint. At the final time (`90 K`), the maximum equivalent (von Mises)
stress reported by the solver is approximately `182 MPa` and the maximum
equivalent strain is approximately `1.5e-3`.

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials`.
Run the included script from this directory:

```bash
./Allrun
```

The script creates the mesh with `blockMesh`, creates the material cell zones
with `setSet` and `setsToZones`, and runs the `solids4Foam` solver. The case
can be cleaned with `./Allclean`.

To use a different temperature history, replace the `T` files in the
`hotCylinderTemperatureField` time directories, e.g. with the results of a
separate heat-transfer simulation on the same mesh.
