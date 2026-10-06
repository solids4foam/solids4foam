# 3dTube: cross-platform `u_z,min` investigation

Between two platforms the 3dTube verification gave the same radial
displacement, arrival time and pulse speed (to 0.1-0.2%) but an axial trough
`u_z,min(A)` that differed by 1.5% on level 1, growing to 3.0% on level 3:

| Platform | OpenFOAM build | solids4foam compile flags | u_z,min(A), level 1 (mm) |
|---|---|---|---:|
| xenosim (AMD EPYC 9684X) | v2412 and v2512 Ubuntu packages, GCC 11.4 | `-O3` (wmake default) | -0.08954 |
| MeluXina (AMD EPYC 7H12) | v2412 EasyBuild foss-2024a, GCC 13.3 | `-O2 -fno-tree-vectorize -march=znver2` (EasyBuild rules) | -0.08819 |
| Apple M1 Ultra (earlier README) | v2412, Clang | | -0.08823 |

Each platform was deterministic and independent of MPI ranks, preconditioner,
tolerances, coupling (Robin-Neumann or IQN-ILS), mesh source and
`FOAM_SETNAN`.

## Root causes

Two independent defects, found in this order. Both are fixed on this branch.

### 1. Aliasing miscompilation of the fluid wall viscous force (level 1 onwards)

**An aliasing miscompilation of OpenFOAM's field inner product, triggered by
the way solids4foam evaluated the fluid wall viscous force.**
(Classification: compiler-dependent undefined behaviour in a field operator;
a force-evaluation, not a discretisation, difference.)

`pimpleFluid::patchViscousForce` (and the other fluid models) computed

```cpp
rho*(mesh().boundary()[patchID].nf() & (-devReff().boundaryField()[patchID]))
```

`nf()` returns a `tmp<vectorField>`. For `tmp<Field<vector>> &
tmp<Field<symmTensor>>`, OpenFOAM reuses the storage of the first tmp for the
vector result (`reuseTmpTmp`), and its field loops (`ListLoopM.H`) access the
result and the operands through `__restrict__` pointers. The result therefore
aliases an input that the compiler is told is not aliased. With `-O3`, GCC 11
and GCC 13 compute the x and y components, store them, and then form the z
component from the overwritten normal:

`z = (n&(-D))_x (-D_xz) + (n&(-D))_y (-D_yz)` instead of
`n_x (-D_xz) + n_y (-D_yz) + n_z (-D_zz)`.

For a wall whose normal has no axial component, the x and y components are
right and the axial (z) traction is wrong. On the first time step of level 1
the integrated axial viscous force on the wall was -2.597e-6 N instead of
-6.519e-6 N; individual faces differed by up to 1.05 Pa. The radial force was
unaffected (bit-identical across platforms).

## Evidence (smallest reproducer first)

1. **Inputs.** Bit-identical case directories (copied between machines,
   meshes identical to 2e-17 m) still gave the split, so case generation is
   excluded.
2. **Plain OpenFOAM agrees.** `pimpleFoam` on the level-1 fluid mesh with a
   rigid wall, one time step: axial viscous force 5.4626698912391e-06 N on
   both platforms (v2412 and v2512 on xenosim bit-identical; MeluXina agrees
   to 1e-15).
3. **Smallest solids4foam reproducer.** Level 1, first time step, one FSI
   iteration, one PIMPLE corrector (`allowUnconvergedCoupling yes`):
   axial interface force -2.597e-6 (xenosim) vs -6.519e-6 N (MeluXina),
   radial components equal to 16 digits.
4. **Fields are identical.** After that solve, U, p, phi, Uf, meshPhi and the
   mesh-motion fields agree to 1e-15 (relative) in every cell and patch face
   (`compareFoamFields.py`).
5. **Wall-force decomposition** (`traction-diagnostics.patch`): the pressure
   contribution to the axial force is 7e-19 N on both; the difference is
   entirely viscous. `snGrad(U)` on the wall and the cell gradient are
   identical on both platforms (the wall cells are orthogonal: `|k| < 1e-19`,
   so the registered `grad(U)` used by `elasticWallVelocity::snGrad` plays no
   role). The wall `devReff` is identical too, and on xenosim
   `rho*(nf & -devReff_b)` evaluated on a *named* copy of the operands, or by
   hand, gives MeluXina's value (-6.519e-6 N); only the one-line expression
   on tmp operands gives -2.597e-6 N.
6. **Standalone OpenFOAM reproducer** (`Test-tmpDotAlias.C`,
   `buildAliasTest.sh`; no solids4foam code): `tmp<vectorField> &
   (-tmp<symmTensorField>)` against a hand loop, 1600 entries:

   | Flags | GCC 11.4 (xenosim) | GCC 13.3 (MeluXina) |
   |---|---:|---:|
   | `-O3` | wrong (max error 7.0e-4) | wrong (7.0e-4) |
   | `-O3 -ffp-contract=off` | wrong | wrong |
   | `-O3 -D__restrict__=` | exact | exact |
   | `-O3 -march=native` | exact | exact |
   | `-O2`, `-O1`, `-O0` | exact | exact |

   The named-operand form is exact in every case. Floating-point contraction
   is not involved; removing `__restrict__` removes the error, so this is
   aliasing undefined behaviour, exposed by the `-O3` code generation for the
   default x86-64 target. The OpenFOAM release does not matter (v2412 and
   v2512 behave the same); the compile flags of the solids4foam library do.
7. **Fix verified.** With the normal as a named field (`const vectorField
   nf(...)`) the xenosim build (`-O3`) gives -6.5190215264804e-6 N on the
   reproducer (MeluXina -6.5190215264811e-6 N) and the full level-1 run gives
   MeluXina's QoIs (u_r,max 0.15978 mm, u_z,min -0.08819 mm, t_arr 5.873 ms,
   c_p 4.654 m/s).

### 2. Degenerate barycentric weights in the FSI point transfer (level 2 onwards)

With the first fix, level 1 agreed between the platforms in every digit, but
level 2 did not (`u_r,max(A)` 0.16000 against 0.15943 mm), and on xenosim it
changed with the number of MPI ranks (0.16000 on 8, 0.16029 on 32), while
MeluXina was rank-independent (0.15944, 0.15943).

The two-step level-2 reproducer (fully converged FSI) located it: the first
solid solve was identical on both platforms (SNES norms to 1e-11), the second
started from a residual 0.6% different, and between them the solid interface
displacement increment, identical on input to 1e-10, came out of
`amiInterfaceToInterfaceMapping::transferPointsZoneToZone` different between
the platforms (0.15-0.2% in the sums) and asymmetric in x and y on both,
although the quarter tube is symmetric about x = y.

The point transfer interpolates with the barycentric weights of the
fan triangle (two face vertices and the face centre) that a target point
projects onto, from `triangle::pointToBarycentric`. That function returns the
"degenerate" weights (1/3, 1/3, 1/3) when `d00*d11 - d01^2 < SMALL`; the
quantity is four times the squared triangle area, in m^4, and `SMALL` is
1e-15. The fan triangles of the interface faces have `4 A^2` of about 1.5e-14
on level 1 (faces 0.39 x 0.63 mm), 9.4e-16 on level 2 (0.20 x 0.31 mm) and
6e-17 on level 3, so **every point of levels 2 and 3 received the weights
(1/3, 1/3, 1/3)** (all 6601 points of level 2 in the dump), while level 1 was
interpolated correctly. The transferred displacement was then a smeared
three-point average, broke the x-y symmetry (solid wall force x and y 0.3%
apart), and, because the average depends on which triangle the addressing
search picks, it changed with decomposition and build.

This is a refinement-dependent error of the coupling itself: it entered the
mesh study between levels 1 and 2.

`triangleWeights.H` computes the same weights with a degeneracy test relative
to the triangle (squared sine of the angle below `SMALL`), used by
`amiInterfaceToInterfaceMapping` and `amiZoneInterpolation` (commit
`0ff7f462b`). With it, level 2 is rank-independent on xenosim (1.5e-6 in the
wall displacement after two steps, against 1% before), the solid wall force
is symmetric to 1e-8, and the platforms agree: `u_r,max(A)` 0.15946 /
0.15946 mm, `u_z,min(A)` -0.08659 / -0.08659 mm, `t_arr(A)` 5.843 / 5.843 ms,
`c_p` 4.635 / 4.634 m/s (xenosim 8 ranks / MeluXina 32 ranks). Level 1 is
unchanged in every digit.

## Fixes

Commit `635ff464f`: the normal is a named field in `patchViscousForce` of
`pimpleFluid` (all three forks), `newtonIcoFluid`, `pimpleOversetFluid` and
`interFluid`, and the same pattern is removed from `newtonIcoFluid`'s
deformed-mesh variant and the `solidTractions` function object. The
underlying hazard (a reused tmp result aliasing an input through
`__restrict__` in compound inner products such as `vector & tensor`,
`tensor & vector`, `tensor & tensor`) is in OpenFOAM itself and can affect any
code, including OpenFOAM's own, compiled with `-O3`; it should be reported
upstream. Other solids4foam expressions of the form `tmp<vectorField> &
tensorField` were not audited.

## Consequences

- Before fix 1, the xenosim results carry a wrong axial wall shear:
  `u_z,min(A)` by 1.5-3% (levels 1-3), `u_r,max(A)` by about 0.1%, the
  arrival time and the late trough by up to 0.3%. The MeluXina and Apple M1
  results were not affected by it.
- Before fix 2, levels 2 and 3 on every platform used the degenerate
  point transfer, and level 1 did not: the earlier mesh studies compared
  differently coupled problems between level 1 and the finer levels. The
  non-monotone `u_r,max(A)` sequence reported before is not evidence about
  the discretisation.
- The corrected mesh study (both fixes) is in `../README.md`.

## Files

- `Test-tmpDotAlias.C`, `buildAliasTest.sh`: standalone OpenFOAM
  reproducer and the flag matrix.
- `compareFoamFields.py`: field-by-field comparison of two ASCII OpenFOAM
  fields (max absolute, max relative, relative L2, location).
- `traction-diagnostics.patch`: the temporary instrumentation (wall-force
  decomposition, `snGrad` lookup, per-face `devReff` and traction), applied to
  the commit before the fix and enabled with `S4F_DEBUG_TRACTION=1`.
