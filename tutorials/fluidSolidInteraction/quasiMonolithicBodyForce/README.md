# quasiMonolithicBodyForce

A small check that the quasi-monolithic (QM) coupling applies body forces
given as finite-volume options (`fvOptions`) in the fluid and in the solid.

## Case

- **Fluid:** an enclosed unit box of fluid (`newtonIcoFluid`, 8 x 8 cells),
  with no-slip walls on three sides.
- **Solid:** a clamped neo-Hookean layer above the box
  (`nonLinearGeometryTotalLagrangianVelocity`, 8 x 2 cells).
- **Coupling:** the fluid and the solid meet at a flat interface at y = 1.

`./Allrun` runs one of three variants:

- `hydrostatic` (default): a uniform fluid body force `g = (0 -10 0)`
  (`vectorSemiImplicitSource` in `constant/fluid/fvOptions`). The fluid stays
  at rest with `p = 10 (1 - y)`, which is zero at the interface, and the solid
  is unloaded and does not deform.
- `gravity`: solid gravity from `constant/solid/g`, the reference solid
  deformation.
- `bodyForce`: the same solid load as an fvOptions body force, `rho g` per unit
  reference volume. It must give the same deformation as `gravity`.

The `hydrostatic` result is an exact discrete solution: the pressure is linear,
and the fluid pressure reference is set to its exact value. It holds only if
the body force enters both:

- the fluid momentum residual;
- the wall-normal momentum balance that sets the `fixedGradient` wall and
  interface pressure gradient.

The case needs PETSc and OpenFOAM.com.

## Regression test

`./regressionTest.sh` runs the three variants for three time steps and checks:

- `hydrostatic`: the maximum fluid velocity is below 1e-9, the pressure
  differs from `10 (1 - y)` by less than 1e-8, and the solid displacement is
  below 1e-12;
- `gravity` and `bodyForce`: the solid displacements agree to 1e-8 relative.
