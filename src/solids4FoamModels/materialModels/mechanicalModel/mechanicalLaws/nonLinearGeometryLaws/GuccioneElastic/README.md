---
sort: 23
---

# GuccioneElastic

This page documents the Guccione transversely isotropic exponential
hyperelastic law, intended for passive myocardium. The runtime type is:

```text
GuccioneElastic
```

The law is anisotropic about a local fibre direction `f0`, which is either a
single uniform vector read from `mechanicalProperties` or a full
`volVectorField` read from a time directory. The implementation follows the
formulation used in E. Garcia-Blanco, R. Ortigosa, A. J. Gil, C. H. Lee and
J. Bonet, _A new computational framework for electro-activation in cardiac
mechanics_, Computer Methods in Applied Mechanics and Engineering, 348 (2019)
796-845. The underlying constitutive model is due to J. M. Guccione,
A. D. McCulloch and L. K. Waldman (1991).

---

## User Guide

### Constitutive relation

The strain energy density is exponential in a quadratic form `Q` of the
Green-Lagrange strain `E`:

```text
W = 0.5*k*(exp(Q) - 1)

S = dW/dE = 0.5*k*exp(Q)*dQ/dE
```

where `S` is the deviatoric part of the 2nd Piola-Kirchhoff stress and `k` has
units of pressure. Two equivalent expressions for `Q` are implemented, selected
by `calculateStressInLocalCoordinateSystem`.

With `calculateStressInLocalCoordinateSystem no` (the default), `Q` is written
in terms of invariants of `E` and the structural tensor `f0*f0`:

```text
I1 = tr(E)
I2 = 0.5*(sqr(tr(E)) - tr(E & E))
I4 = E && (f0*f0)
I5 = (E & E) && (f0*f0)

Q  = ct*sqr(I1) - 2*ct*I2
   + (cf - 2*cfs + ct)*sqr(I4)
   + 2*(cfs - ct)*I5
```

With `calculateStressInLocalCoordinateSystem yes`, `E` is first rotated into
the orthonormal fibre frame `(f0, s0, n0)` to give `EStar`, and `Q` takes the
classical component form:

```text
Q = cf*sqr(E11)
  + ct*(sqr(E22) + sqr(E33) + 2*sqr(E23))
  + cfs*(2*sqr(E12) + 2*sqr(E13))
```

The Cauchy stress is then assembled differently in the two solver modes. In
the standard (displacement-only) mode a penalty bulk term enforces near
incompressibility:

```text
sigma = dev(symm(F & S & F.T))/J + sigmaHyd*I

sigmaHyd = 0.5*bulkModulus*(sqr(J) - 1)/J
```

In pressure-displacement mode the hydrostatic part comes from the solved
pressure field `p` instead:

```text
sigma = dev(symm(F & S & F.T)) - p*I
```

### Fibre, sheet and sheet-normal directions

Only `f0` is user input. The constructor normalises it, then constructs the
sheet direction `s0` by removing the `f0` component from the global `x`
direction (falling back to `y` where that degenerates), normalises it, and sets
`n0 = f0 ^ s0`. The rotation tensor `R` has `f0`, `s0` and `n0` as its columns
and is used for the local-coordinate-system branch. Because `s0` and `n0` are
generated automatically rather than supplied, the transverse response is
isotropic in the plane normal to `f0` unless `ct` and `cfs` are chosen to
differ.

### Model options

Required entries in the material sub-dictionary:

| Entry | Type | Description |
| --- | --- | --- |
| `k` | `dimensionedScalar` | Exponential stress scale, `[1 -1 -2 0 0 0 0]` |
| `cf` | `scalar` | Fibre-direction stiffness coefficient |
| `ct` | `scalar` | Transverse stiffness coefficient |
| `cfs` | `scalar` | Fibre-transverse shear coefficient |
| `rho` | `dimensionedScalar` | Density, read by the base `mechanicalLaw` |

`bulkModulus` is also read with `lookup()` and is therefore required, unless
`pressureDisplacement` is `yes`, in which case it is silently set to `GREAT`
and never read.

Optional entries:

| Entry | Default | Description |
| --- | --- | --- |
| `pressureDisplacement` | `false` | Pressure-displacement formulation |
| `uniformFibreField` | `false` | Read `f0` as a single vector |
| `f0` | none | Fibre vector; required if `uniformFibreField` |
| `calculateStressInLocalCoordinateSystem` | `false` | Use `EStar` form of `Q` |
| `impKcoeff` | `1.0` | Scales `impK`, pressure-displacement only |
| `writeS0N0R` | `false` | Write `s0`, `n0` and `R` at construction |

`tangentEps` is read with `lookup()` inside `materialTangentField()`. It is
only needed by solid models that request a material tangent, but it is a hard
requirement for those.

When `uniformFibreField` is `false`, a `volVectorField` named `f0` must exist
in the start time directory or in `0`; the code tries the current time
directory first and falls back to `0`, and aborts if neither has it.

### Compatibility

`correct(volSymmTensorField&)` is implemented for both modes, so the law works
with the cell-centred nonlinear-geometry solid models, both total and updated
Lagrangian (`updateF()` handles both).

`correct(surfaceSymmTensorField&)` is implemented **only** when
`pressureDisplacement` is enabled; otherwise it hits `notImplemented`. The
face-based path is what `coupledPressureDisplacementSolid` uses, and that model
also supplies the `p` and `pf` fields that the pressure-displacement branch
looks up from the registry. Selecting `pressureDisplacement yes` with a solid
model that does not create `p` will fail at the first `correct()` call.

### Recommended dictionary setup

Standard penalty-bulk-modulus form, for example with
`nonLinearGeometryTotalLagrangianTotalDisplacement`:

```text
mechanical
(
    myocardium
    {
        type            GuccioneElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;

        k               k [1 -1 -2 0 0 0 0] 2000;
        cf              8;
        ct              2;
        cfs             4;

        // Penalty bulk modulus; aim for ~100 x the shear modulus
        bulkModulus     bulkModulus [1 -1 -2 0 0 0 0] 1e6;

        uniformFibreField yes;
        f0              (0 0 1);

        calculateStressInLocalCoordinateSystem no;
    }
);
```

Pressure-displacement form, for `coupledPressureDisplacementSolid`:

```text
mechanical
(
    myocardium
    {
        type            GuccioneElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        nu              nu [0 0 0 0 0 0 0] 0.5;

        pressureDisplacement yes;
        impKcoeff       2;

        k               k [1 -1 -2 0 0 0 0] 2000;
        cf              8;
        ct              2;
        cfs             4;

        calculateStressInLocalCoordinateSystem true;
    }
);
```

Note that `nu` is not read by this law; it appears in the tutorials because
other parts of the case setup expect it.

### Field glossary

- `f0`: reference fibre direction, unit length after construction.
- `s0`, `n0`: generated sheet and sheet-normal directions; written only if
  `writeS0N0R` is enabled.
- `R`, `Rf`: rotation from the local fibre frame to global Cartesian axes.
- `f0f0`, `f0f0f`: structural tensor `sqr(f0)` at cells and faces.
- `S2PK`, `S2PKf`: the deviatoric 2nd Piola-Kirchhoff stress `S`.
- `muEff`: effective shear modulus used for `impK` in pressure-displacement
  mode; recomputed at the end of every time step.
- `expQf`: face field holding `exp(Q)`, under-relaxed between outer iterations
  in the pressure-displacement face path.

---

## Developer Notes

### Class role

`GuccioneElastic` derives from `mechanicalLaw`. It carries two parallel sets of
fields, cell-centred and face-centred, so the same law can serve both the
standard cell-centred solid models and `coupledPressureDisplacementSolid`.

### Construction

The constructor:

1. Reads `bulkModulus` via the file-local helper `readBulkModulus()`, which
   short-circuits to `GREAT` when `pressureDisplacement` is set.
2. Reads `k`, `cf`, `ct`, `cfs` and forms the linearised shear modulus
   `mu = 0.5*k*(cf + cfs + ct)/3`.
3. Builds `f0` and `f0f` through `makeF0()`/`makeF0f()`, either from a uniform
   vector or from a field on disk.
4. Stores old-time copies of `F` and `Ff`.
5. Normalises `f0`, constructs `s0` and `n0`, and assembles `R` and `Rf`.
6. If `pressureDisplacement`, calls `calcInitialShearModulus()` and
   `calcEffectiveShearModulus()`.

### `correct()` walkthrough

`correct(volSymmTensorField&)` first calls `updateF()`. If that returns `true`
the solver has enforced linearised elasticity for stability and the function
returns immediately. Otherwise:

- `C = symm(F.T() & F)` and `E = 0.5*(C - I)`;
- `Q` and `dQ/dE` are formed by one of the two branches above;
- `S_ = dQdE*0.5*k*exp(Q)` (rotated back to global axes in the local branch);
- the deviatoric Cauchy stress is `dev(symm(F & S & F.T))/J`;
- `updateSigmaHyd()` supplies the hydrostatic term, with implicit stiffness
  `(4/3)*mu + bulkModulus`.

In the pressure-displacement branch the `/J` division is omitted and the
hydrostatic term is `-p*I`, taken from the registry field `p`.

`correct(surfaceSymmTensorField&)` mirrors this for faces, but only in
pressure-displacement mode. It additionally stores and under-relaxes `expQf_`
between outer iterations, which damps the strong exponential nonlinearity.

### Implicit stiffness `impK()`

```text
pressureDisplacement:  impK = impKcoeff*muEff
otherwise:             impK = (4/3)*mu + bulkModulus
```

`muEff` is not a material constant. `calcEffectiveShearModulus()` perturbs the
displacement gradient by `0.001` in each of the three shear planes of the
principal stress frame, evaluates the deviatoric Cauchy stress through
`calcDevCauchy()`, and takes the largest resulting secant shear stiffness. It
does this for both the old-time and current deformation gradients and keeps the
maximum. `updateTotalFields()` refreshes it once per time step, so `impK`
stiffens automatically as the exponential response steepens.

`calcInitialShearModulus()` runs the same perturbation procedure once at
construction to replace the analytic small-strain `mu`.

### Material tangent

`materialTangentField()` builds a `mat66` per face by finite-differencing the
face stress with respect to each component of `grad(D)f`, using the step size
`tangentEps`. It calls the private `calculateStress()` helper, which is
currently `notImplemented`, so the tangent path is not usable as it stands.

### Extension points

- Supply `s0` explicitly rather than generating it, if genuine orthotropy
  (`cf`, `ct`, `cfs` all distinct with a meaningful sheet plane) is wanted.
- Implement `calculateStress(surfaceSymmTensorField&, const surfaceTensorField&)`
  to enable the material tangent and the non-pressure-displacement face path.
- Replace the penalty `bulkModulus` with a different volumetric energy in the
  `updateSigmaHyd()` call.

---

## Tutorials

Cases that select `GuccioneElastic`:

- `solids/hyperelasticity/heartTissueBeam` (pressure-displacement)
- `solids/hyperelasticity/idealisedVentricle/caseOptions/pressureDisplacement`
- `solids/hyperelasticity/idealisedVentricle/caseOptions/petsc`, where it is
  the passive law nested inside `electroMechanicalLaw`
