---
sort: 23
---

# GuccioneElastic

The Guccione transversely isotropic exponential hyperelastic law, intended
for passive myocardium. The runtime type is:

```text
GuccioneElastic
```

The response is anisotropic about a fibre direction `f0`, given either as a
field or as a single vector. The underlying constitutive model is due to
J. M. Guccione, A. D. McCulloch and L. K. Waldman (1991); the formulation
follows E. Garcia-Blanco, R. Ortigosa, A. J. Gil, C. H. Lee and J. Bonet
(2019).

---

## User Guide

### What it computes

This is a **finite-strain** law: it implements the finite-strain update only,
so it is used with a nonlinear geometry solid model. A linear geometry solid
model asks for the small-strain update, and the run stops with an error.

The strain energy is split into an isochoric and a volumetric part:

```text
W = 0.5*k*(exp(Q) - 1) + U(J)

Fbar = J^(-1/3)*F
E    = 0.5*(Fbar^T & Fbar - I)        (isochoric Green-Lagrange strain)

I1 = tr(E)
I2 = 0.5*(sqr(tr(E)) - tr(E & E))
I4 = E && (f0*f0)
I5 = (E & E) && (f0*f0)

Q  = ct*sqr(I1) - 2*ct*I2
   + (cf - 2*cfs + ct)*sqr(I4)
   + 2*(cfs - ct)*I5
```

The Cauchy stress is

```text
S     = 0.5*k*exp(Q)*dQ/dE
sigma = dev(symm(Fbar & S & Fbar^T))/J + dU/dJ*I

U(J)  = 0.25*K*(J^2 - 1 - 2*ln(J))
dU/dJ = 0.5*K*(J^2 - 1)/J
```

where `K` is `bulkModulus`, a penalty enforcing near incompressibility.
Setting `cf = ct = cfs` makes `Q` isotropic, and the fibre direction then has
no effect.

The sheet and sheet-normal directions are not used: the invariant form above
is the only form of `Q` implemented, so the response is transversely
isotropic about `f0`.

A deformation gradient with `J <= sqrt(SMALL)` stops the run.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `k` | yes | - | Stress scale of the exponential, pressure units |
| `cf` | yes | - | Fibre coefficient of `Q`, dimensionless scalar |
| `ct` | yes | - | Transverse coefficient of `Q`, dimensionless scalar |
| `cfs` | yes | - | Fibre-sheet coefficient of `Q`, dimensionless scalar |
| `bulkModulus` | yes | - | Volumetric penalty `K`, pressure units, `> 0` |
| `uniformFibreField` | no | `no` | Take the default `f0` from this dictionary |
| `f0` | if uniform | - | Fibre vector; required if `uniformFibreField yes` |

`bulkModulus` must be positive; zero or a negative value is a fatal error.
`f0` need not be of unit length: it is normalised before use.

### Fibre direction

`f0` is a prescribed state variable (see below). When the manager sets up the
law's state it looks for a `volVectorField` named `f0`, in this order:

1. a file `f0` in the current (start) time directory;
2. a file `f0` in `0`;
3. a field `f0` already registered by another model.

If one is found it is used, and the cell values are mapped to the law's
integration points (interpolated to faces for face-based topologies).
Only if none is found is the dictionary value used: `f0` when
`uniformFibreField` is `yes`, otherwise a zero vector.

A zero-length fibre stops the run at the first evaluation, naming the two
ways of supplying one. So a case must either provide an `f0` field (the
`setFibreField` utility writes one for the `LandEtAl2015` benchmarks) or set
`uniformFibreField yes` with `f0` in the dictionary. Note that an `f0` field
in the time directory takes precedence over the dictionary vector even when
`uniformFibreField` is `yes`.

### State variables

| Name | Type | Role | Description |
| --- | --- | --- | --- |
| `f0` | vector | prescribed | Reference fibre direction |

A prescribed variable is read once and never written back, so `f0` is not
written with the results and is not part of the restart files; the `f0` file
in `0` is found again on restart. The law has no persistent history.

### Tangents

- `scalar`: `(4/3)*mu + K`, with `mu = 0.5*k*(cf + cfs + ct)/3`, the
  small-strain shear modulus implied by `k` and the coefficients. It is a
  constant stiffness estimate for the segregated solver, not a tangent of the
  exponential energy.
- `scalarDeviatoric`: `(4/3)*mu`, for mixed displacement-pressure
  formulations.
- `fourthOrderFiniteDifference`: supported, through the base class.
- `fourthOrder` (analytical): not implemented; requesting it is a fatal
  error.

### Mixed displacement-pressure formulation

The law provides the volumetric split: because `Q` is written on the
isochoric strain, it can return the isochoric stress and `dU/dJ` separately.
It can therefore be used with a mixed formulation, such as
[`coupledPressureDisplacementSolid`](../../../solidModels/coupledPressureDisplacementSolid/README.md)
or
[`nonLinearGeometryTotalLagrangianTotalDisplacement`](../../../solidModels/nonLinGeomTotalLagTotalDispSolid/README.md)
with `solvePressure true`, where the solved pressure replaces `dU/dJ`.
`coupledPressureDisplacementSolid` uses `1/bulkModulus` as the
compressibility of its pressure equation, so a very large `bulkModulus` (the
tutorials use `1e15`) gives a fully incompressible material.

The mixed formulation is chosen by the solid model, not by this dictionary.

### Differences from the legacy law

- `pressureDisplacement`, `impKcoeff`,
  `calculateStressInLocalCoordinateSystem`, `writeS0N0R` and `tangentEps` are
  no longer read. They are not reported as errors, so remove them from
  migrated cases. The local coordinate system variant of `Q` is gone.
- `bulkModulus` is always required, including for
  `coupledPressureDisplacementSolid`, where the legacy law set it to `GREAT`.
- The law never reads the solid model's `p` field: the pressure enters
  through the volumetric split.
- `Q` is taken on the isochoric strain rather than the full Green-Lagrange
  strain. The two agree as `J` approaches one; on the
  `LandEtAl2015/problem3` tutorial they differ by about `4e-4`.
- `f0` is normalised, and a separate face field `f0f` is no longer needed.
- The pressure-displacement implicit stiffness, a per-step perturbation
  estimate of an effective shear modulus, is replaced by the constant
  scalar tangents above.

### Example

From the `heartTissueBeam` tutorial, used with
`coupledPressureDisplacementSolid`:

```text
planeStress no;

mechanical
(
    rubber
    {
        type        GuccioneElastic;
        k           k [ 1 -1 -2 0 0 0 0 ] 2000;
        cf          8;
        ct          2;
        cfs         4;
        rho         rho [ 1 -3 0 0 0 0 0 ] 1000;
        bulkModulus bulkModulus [ 1 -1 -2 0 0 0 0 ] 1e15;
    }
);
```

This case gives no `f0` in the dictionary; the fibre field is read from
`0/f0`. With a uniform fibre direction instead:

```text
        uniformFibreField yes;
        f0                (0 0 1);
```

---

## Tutorials

- [`solids/hyperelasticity/heartTissueBeam`](../../../../../tutorials/solids/hyperelasticity/heartTissueBeam)
- [`solids/hyperelasticity/idealisedVentricle`](../../../../../tutorials/solids/hyperelasticity/idealisedVentricle),
  in both `caseOptions/petsc` and `caseOptions/pressureDisplacement`
- [`solids/hyperelasticity/LandEtAl2015/problem3`](../../../../../tutorials/solids/hyperelasticity/LandEtAl2015/problem3),
  as the passive law inside
  [`electroMechanicalLaw`](../electroMechanicalLawMechanicalConstitutiveLaw/README.md)

---

## References

- J. M. Guccione, A. D. McCulloch and L. K. Waldman, _Passive material
  properties of intact ventricular myocardium determined from a cylindrical
  model_, Journal of Biomechanical Engineering, 113 (1991) 42-55.
- E. Garcia-Blanco, R. Ortigosa, A. J. Gil, C. H. Lee and J. Bonet, _A new
  computational framework for electro-activation in cardiac mechanics_,
  Computer Methods in Applied Mechanics and Engineering, 348 (2019) 796-845.
