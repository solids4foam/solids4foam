# beamInCrossFlow Richter verification

This opt-in study verifies the solids4foam discretisation against the
stationary three-dimensional benchmark introduced by Richter (2012). It is
separate from the tutorial regression test and is not run by normal CI.

## Commands

After sourcing a supported OpenFOAM environment and building solids4foam with
PETSc, run from this directory:

```bash
# Richter-consistent graded spatial family
./Allverify --study mesh --family graded

# Individual factors for scheduled runs
./Allverify --study mesh --family graded --levels 4,8 --cores 128

# L3 time-step and coupling-tolerance controls
./Allverify --study temporal --family graded --cores 128
./Allverify --study coupling --family graded --cores 128

# Static structural diagnostic
python3 scripts/solid_discretisation.py --levels 1,2,4,8 --cores 8
```

The driver creates isolated copies under `work/` and compact CSV output under
`postProcessing/`; both directories are ignored. Versioned reference CSVs in
`reference/` contain the completed production evidence. No raw transient
fields are versioned.

## Exact Richter definition

The primary source is Thomas Richter, *Goal-oriented error estimation for
fluid-structure interaction problems*, Computer Methods in Applied Mechanics
and Engineering 223--224 (2012), 28--42,
doi:10.1016/j.cma.2012.02.014, section 7.2 and Tables 5--6.

<!-- markdownlint-disable MD013 -->

| Item | Richter definition | solids4foam production definition |
| --- | --- | --- |
| Fluid domain | `(0,1.5) x (0,0.4) x (-0.4,0.4) m` | Same; negative-`z` symmetry half |
| Solid | `(0.4,0.5) x (0,0.2) x (-0.2,0.2) m` | Same; negative-`z` half |
| Symmetry | `x-y` plane | `z=0` fluid and solid symmetry |
| Inlet | bi-parabolic profile with explicit coefficient and peak `0.3 m/s` | Same, after a numerical ramp |
| Reynolds-number statement | average speed `0.2 m/s` and `Re=40`, inconsistent with the printed profile | Retained as a documented source inconsistency, not used to alter the inlet equation |
| Fluid | `rho=1000 kg/m3`, `nu=0.001 m2/s` | Same |
| Outlet | do-nothing / zero traction | zero gauge kinematic pressure and zero-gradient velocity |
| Other fluid walls | no slip | Same |
| Solid law | compressible St Venant--Kirchhoff | Same |
| Solid constants | shear modulus `0.5 MPa`, Poisson ratio `0.4` | `E=1.4 MPa`, `nu=0.4`, hence the same `mu` and `lambda=2 MPa` |
| Clamp | solid base at `y=0` | Same |
| Point A | `(0.45,0.15,0.15) m` | symmetry-equivalent `(0.45,0.15,-0.15) m` |
| Drag | fluid traction integral on the half-solid interface | pressure plus viscous traction on the same half-interface |
| Problem type | stationary monolithic ALE FSI | partitioned transient-to-steady continuation |

<!-- markdownlint-enable MD013 -->

The inlet field is

```text
u_x = 0.3 y(0.4-y)(0.4^2-z^2)/(0.2^2 0.4^2),  u_y=u_z=0.
```

The coefficient and the immediately following statement `vmax=0.3` are the
most explicit definition of the boundary condition in the primary source and
are therefore used without tuning to Table 6. In the same paragraph Richter
states an average inlet speed of `0.2 m/s` and `Re=40`; however, direct
integration of the printed bi-parabolic profile gives an area average of
`(4/9) vmax = 0.1333... m/s`. This inconsistency means that the tabulated
reference cannot yet be assumed to represent the printed boundary condition.
The custom inlet condition implements the printed equation exactly. A
one-second cosine ramp is used only to start the partitioned calculation
robustly. The reported state must satisfy the steady criterion after the
boundary value has been constant; the ramp is not part of the benchmark
physics. Richter does not state a solid density for this stationary example.
Solids4foam uses
`1000 kg/m3` during continuation, but density drops out of the converged
stationary balance. A temporal control checks that the continuation does not
pollute the reported quantities.

The old solids4foam definition placed the solid at `x=0.45...0.55 m` and
sampled its upstream face rather than its mid-thickness. Those material changes
have been replaced. The old definition is not retained as another production
variant; git history preserves the earlier study.

## Richter reference provenance

Richter solves a stationary, monolithic ALE formulation with equal-order
piecewise-linear finite elements and local-projection stabilisation. Table 5
reports five uniformly refined levels. The finest has 7,600,775 algebraic
unknowns, `F_x=1.3380 N`, and `u_x(A)=5.9202e-5 m`. Table 6 extrapolates the
finest three levels:

| Quantity | Reference | Stated accuracy |
| --- | ---: | ---: |
| `u_x(A)` | `5.924e-5 m` | `+/-1e-7 m` |
| `F_x` | `1.327 N` | `+/-0.01 N` |

The values are tabulated, not digitised. Richter says their relative accuracy
is at most about 1%. Re-entrant solid corners prevent the expected improvement
from piecewise-quadratic elements, so this is a strong independent source with
finite and explicitly retained numerical uncertainty, not an exact solution.
Richter defines interface drag and evaluates an equivalent solid-base/residual
functional. The solids4foam surface-traction integral is physically equivalent
at equilibrium but not algebraically identical, which is retained as a small
comparison-method qualification.

## Graded family

The conformal 11-block fluid topology and one-block solid topology are fixed.
Every cell count is multiplied by the level factor. Fixed total block
expansion ratios `xUp=0.125`, `xDown=8`, `yOuter=6`, and `zOuter=1/6` cluster
cells on the beam faces, free corner, and near wake while coarsening toward the
outer boundaries. Solid and fluid interface counts match at every level.

<!-- markdownlint-disable MD013 -->

| Level | Factor | Fluid | Solid | Near `h` (m) | Far `h` (m) | Through thickness | Interface `y x z` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| L0 | 1 | 7,584 | 256 | 0.00960629 | 0.0975879 | 4 | `8 x 8` |
| L1 | 2 | 60,672 | 2,048 | 0.00487940 | 0.0491578 | 8 | `16 x 16` |
| L2 | 4 | 485,376 | 16,384 | 0.00245789 | 0.0246677 | 16 | `32 x 32` |
| L3 | 8 | 3,883,008 | 131,072 | 0.00123339 | 0.0123558 | 32 | `64 x 64` |
| L4 | 16 | 31,064,064 | 1,048,576 | 0.000617792 | 0.00618338 | 64 | `128 x 128` |

<!-- markdownlint-enable MD013 -->

The table values are calculated from the geometric-series block definitions;
final `checkMesh` measurements and quality results are recorded with the
production evidence.

## Numerical controls

All definitive coupled levels use one MeluXina build and software stack. The
spatial family holds a documented timestep selected by the fine-mesh startup
and temporal controls. Production IQN-ILS uses direct mapping,
`outerCorrTolerance=1e-6`, requires both normalized residual measures, permits
at most 200 iterations, uses prediction and no reused modes, and applies
relative QR filtering `0.01`. A fine level is repeated at tighter tolerance to
quantify iterative error.

The continuation runs to `t=8 s`. A primary quantity is treated as steady only
when its relative change over `t=7...8 s` is below `0.1%`; otherwise the run is
extended. Temporal and coupling changes must be materially smaller than the
last retained spatial change.

## Static solid diagnostic

The exact three-dimensional Richter solid is loaded by a uniform `100 Pa`
horizontal traction on its upstream face, with all other wetted faces
traction-free. Four systematically refined solid meshes use the coupled
family's solid/interface resolution. The diagnostic records displacement at
the exact material point A and compares it with the slender Euler--Bernoulli
value for the corresponding full-width uniform load. The analytical value is
diagnostic only because it omits three-dimensional Poisson and clamp/end
effects.

<!-- SOLID_RESULTS -->

## Mesh quality

<!-- MESH_QUALITY -->

## MeluXina environment

<!-- ENVIRONMENT -->

## Spatial convergence

<!-- SPATIAL_RESULTS -->

## Temporal and coupling controls

<!-- CONTROL_RESULTS -->

## Numerical uncertainty and Richter comparison

<!-- UNCERTAINTY_RESULTS -->

## Limitations

<!-- LIMITATIONS -->
