---
sort: 15
---

# Flexible flag behind a block in a wind flow: `foilInWind`

Prepared by Philip Cardiff

---

## Tutorial Aims

- Demonstrate a three-dimensional fluid-solid interaction case with a thin,
  compliant structure in an air flow.
- Demonstrate the partitioned (IQNILS and Aitken) and quasi-monolithic
  (Newton) coupling approaches on the same case.
- Demonstrate the `decomposeParMonolithic` utility, which decomposes the fluid
  and solid regions together for a parallel quasi-monolithic run.

---

## Case Overview

A thin flexible foil (flag) is clamped to the downstream face of a rigid block
and exposed to a uniform wind flow. The foil is thin and compliant
(rubber-like), so the wind loading can cause significant bending deflection.
The case exercises large-displacement nonlinear structural response coupled to
incompressible laminar fluid flow.

The air domain is $$20 \times 11 \times 11$$ cm ($$x \times y \times z$$). The
rigid block occupies $$5 \le x \le 6$$ cm, $$5 \le y \le 6$$ cm and
$$4 \le z \le 7$$ cm, and the foil occupies $$6 \le x \le 10$$ cm,
$$5.47 \le y \le 5.53$$ cm and $$4 \le z \le 7$$ cm, i.e. it is 4 cm long, 3
cm wide and 0.6 mm thick. The foil is fixed on the face attached to the block
(`fixed` patch) and loaded by the air on its other faces. The fluid and solid
meshes share a common interface (`flag`/`interface` patch pair) across which
displacement and traction are exchanged.

Air enters at $$x = 0$$ with a uniform velocity of 1 m/s, the outlet has a zero
pressure, the four lateral walls are slip walls, and the block is a no-slip
wall. The Reynolds number based on the block side is about 650. The flow is
laminar, and the fluid velocity is initialised with a potential flow solution
(`potentialFoam`).

### Physical parameters

| Property | Value |
|----------|-------|
| Fluid kinematic viscosity | 1.542e-5 m²/s |
| Fluid density | 1.18 kg/m³ |
| Solid density | 2000 kg/m³ |
| Young's modulus | 2 MPa |
| Poisson's ratio | 0.35 |
| Material model | St. Venant-Kirchhoff nonlinear elastic |

### Simulation parameters

| Parameter | Value |
|-----------|-------|
| End time | 12 s |
| Time step | 0.01 s |
| Total steps | 1200 |

The fluid mesh is created with `blockMesh` and refined once around the foil
with `setSet` and `refineMesh`, so that the fluid `flag` patch conforms to the
solid `interface` patch, giving 176 248 fluid cells. The solid mesh has
$$32 \times 6 \times 24$$ cells. The interface quantities are transferred with
the `directMap` method.

The partitioned approach uses the `pimpleFluid` fluid model, a
`nonLinearGeometryTotalLagrangianTotalDisplacement` solid model solved with
PETSc SNES, and IQNILS (default) or Aitken coupling. The quasi-monolithic
approach (`fluidSolidInterface NewtonQuasiMonolithic`) solves the fluid
velocity and pressure (`newtonIcoFluid`) and the solid velocity
(`nonLinearGeometryTotalLagrangianVelocity`, with the same St. Venant-Kirchhoff
law) together with PETSc SNES in every time step, so there are no FSI
sub-iterations; the fluid mesh motion is a separate linear elastic solid model
in the fluid region.

---

## Running the Case

The tutorial case can be run using the included `Allrun` script:

```bash
./Allrun                      # partitioned IQNILS, serial
./Allrun parallel             # partitioned IQNILS, parallel
./Allrun aitken               # partitioned Aitken, serial
./Allrun aitken parallel      # partitioned Aitken, parallel
./Allrun monolithic           # quasi-monolithic Newton, serial
./Allrun monolithic parallel  # quasi-monolithic Newton, parallel
```

The default files are the partitioned IQNILS setup. The Aitken and
quasi-monolithic files are stored beside them as `*.aitken` and `*.monolithic`
files, which `Allrun` links in and `Allclean` removes again. The case needs
OpenFOAM.com and a solids4foam build with PETSc; with other OpenFOAM variants,
or without PETSc, `Allrun` exits without running. `Allrun` also writes
`deflection.pdf` and `force.pdf` with gnuplot, showing the displacement of the
mid-span point of the foil trailing edge and the force on the foil.

In the partitioned parallel run, the fluid and solid regions are decomposed
independently with the dictionaries in `system/fluid` and `system/solid`. In
the quasi-monolithic parallel run, the
[`decomposeParMonolithic`](https://github.com/solids4foam/solids4foam/tree/master/applications/utilities/decomposeParMonolithic)
utility decomposes both regions together as one graph, connected across the
`flag`/`interface` patch pair, using `system/decomposeParDict.monolithic`:

```c++
numberOfSubdomains 4;

method          scotch;

monolithicCoeffs
{
    regions (fluid solid);

    regionWeights
    {
        fluid   4;
        solid   3;
    }

    interfaces
    (
        {
            regionA     fluid;
            patchA      flag;
            regionB     solid;
            patchB      interface;
        }
    );
}
```

The PETSc options of the quasi-monolithic system are in
`system/fvSolution.monolithic`: a Newton method with a backtracking line
search and matrix-free Jacobian-vector products, and an LGMRES linear solver
with a multiplicative fluid/solid field-split preconditioner. The fluid block
uses a Schur complement split of velocity and pressure, with ILU(0) for the
velocity and hypre BoomerAMG for the pressure; the solid block is gathered
onto one rank and factorised with LU (`redundant`), which keeps the solid
preconditioner exact for any number of ranks.

---

## Expected Results

The geometry and the inflow are symmetric about the mid-plane of the foil, so
the foil is only loaded by the drag of the air flowing along it until the
wake of the block becomes unsteady. Over the first $$0.2$$ s, the foil
trailing edge moves less than $$10^{-7}$$ m, and the drag force on the foil
decays from about $$2 \times 10^{-5}$$ N, during the start-up transient, to
about $$2 \times 10^{-6}$$ N. The partitioned and quasi-monolithic approaches
give a maximum streamwise trailing-edge displacement within $$1\%$$ of each
other, and drag forces within about $$10\%$$; at these tiny loads, the
cross-stream displacement (below $$2 \times 10^{-8}$$ m) is at the level of
the solver tolerances and differs between the approaches. The flapping of
the foil, if it develops later in the run, has not been verified.

The full run is expensive: on an Apple M1 Ultra, the first $$0.2$$ s took
about 13 min for the partitioned IQNILS setup on 2 processes (5 to 15
coupling iterations per time step) and about 29 min for the quasi-monolithic
setup on 2 processes (3 to 14 Newton iterations per time step).
