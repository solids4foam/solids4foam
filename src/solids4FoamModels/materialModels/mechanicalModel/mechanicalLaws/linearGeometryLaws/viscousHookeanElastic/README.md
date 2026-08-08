---
sort: 9
---

# viscousHookeanElastic

This page documents the small-strain generalised Maxwell law, with an elastic
volumetric response and a relaxing deviatoric response. The runtime type is:

```text
viscousHookeanElastic
```

---

## User Guide

### What it computes

The law places a relaxed spring in parallel with a list of Maxwell arms. It
first forms the instantaneous modulus and relative moduli:

```text
E0        = EInfinity + sum(E[i])
gammaInf  = EInfinity/E0
gamma[i]  = E[i]/E0
mu        = E0/(2*(1 + nu))

lambdaInf = nu*EInfinity/((1 + nu)*(1 - 2*nu))  // plane strain / 3-D
lambdaInf = nu*EInfinity/((1 + nu)*(1 - nu))    // planeStress yes
k         = lambdaInf + (2/3)*gammaInf*mu
```

For a non-incremental solid model, the trial deviatoric stress and elastic
volumetric stress are:

```text
e       = dev(symm(grad(D)))
s       = 2*mu*e
sigmaVol = k*tr(grad(D))*I
```

For an incremental solid model, the implemented update is:

```text
De       = dev(symm(grad(DD)))
s        = s.oldTime() + 2*mu*De
sigmaVol = tr(sigma.oldTime())*I + k*tr(grad(DD))*I
```

Each Maxwell-arm stress is then advanced using the midpoint exponential
update. With `aT = 1` unless temperature shifting is enabled, the update and
total stress are:

```text
aT = 10^(-C1*(T - Tref)/(C2 + (T - Tref)))

h[i] = exp(-deltaT/(aT*tau[i]))*h[i].oldTime()
       + exp(-deltaT/(2*aT*tau[i]))*(s - s.oldTime())

sigma = sigmaVol + gammaInf*s + sum(gamma[i]*h[i])
```

Only the deviatoric response relaxes; `k` is constant. The law implements
`correct(volSymmTensorField&)` and `correct(surfaceSymmTensorField&)`. The
point-symmetric-tensor overload inherits `notImplemented`. On OpenFOAM.com and
OpenFOAM.org, the `CompactListList` quadrature-point overload also inherits
`notImplemented`; that overload is not compiled for foam-extend.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `type` | yes | Runtime type; `viscousHookeanElastic` |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `EInfinity` | yes | Relaxed modulus, `[1 -1 -2 0 0 0 0]` |
| `E` | yes | Maxwell-arm modulus list, in pressure units |
| `relaxationTimes` | yes | Maxwell-arm relaxation-time list |
| `nu` | yes | Relaxed Poisson's ratio, dimensionless |
| `WilliamsLandelFerry` | no | Enable WLF shifting; default `no` |
| `WilliamsLandelFerryCoeffs` | with WLF | Dictionary containing WLF data |
| `C1` | with WLF | Dimensionless first WLF coefficient |
| `C2` | with WLF | Second WLF coefficient, `[0 0 0 1 0 0 0]` |
| `Tref` | with WLF | Reference temperature, `[0 0 0 1 0 0 0]` |
| `regionName` | no | Base mesh region selected by name |
| `solvePressureEqn` | no | Base switch; default `no`, unused here |
| `pressureSmoothingScaleFactor` | no | Default `100`; unused here |

`E` and `relaxationTimes` are read as scalar lists and must have equal lengths.
The code uses each relaxation time with `deltaTValue()`, so the values have the
case time unit. Enabling `WilliamsLandelFerry` also requires a cell-centred `T`
field for the cell correction and a face-centred `T` field for the face
correction.

The base class selects `solid` as `regionName` when that mesh exists, then
`region0`; otherwise `regionName` must be supplied. Although the base class
reads the two pressure-equation options, this law does not call
`updateSigmaHyd()`, so neither option changes its stress update.

### Recommended dictionary setup

The following polymer values are used by the `viscoTube` tutorial:

```text
planeStress     no;

mechanical
(
    polymer
    {
        type            viscousHookeanElastic;
        rho             rho [1 -3 0 0 0 0 0] 7850;
        EInfinity       EInfinity [1 -1 -2 0 0 0 0] 39.58e9;
        nu              nu [0 0 0 0 0 0 0] 0.33;
        E               (2.9318518519e9 5.8637037037e9
                         6.5966666667e9 18.3240740741e9);
        relaxationTimes (30 300 3000 12000);

        // Optional temperature shift
        // WilliamsLandelFerry yes;
        // WilliamsLandelFerryCoeffs
        // {
        //     C1   C1 [0 0 0 0 0 0 0] 17.44;
        //     C2   C2 [0 0 0 1 0 0 0] 51.6;
        //     Tref Tref [0 0 0 1 0 0 0] 293.15;
        // }
    }
);
```

### Field glossary

- `s`, `sf`: instantaneous deviatoric trial stress at cell centres and faces.
  Their old-time values are stored, but the fields are neither read nor
  written.
- `h0`, `h1`, and so on: cell-centred Maxwell-arm internal stresses. Each field
  stores old-time and previous-iteration values and is neither read nor
  written.
- `hf0`, `hf1`, and so on: face-centred counterparts of the Maxwell-arm
  stresses, with the same time and iteration storage.
- `T`: temperature looked up when WLF shifting is enabled. Its location must
  match the selected cell- or face-centred correction.
- `K`: temporary cell field containing the constant relaxed bulk modulus `k`.
- `impK`: temporary cell field containing the implicit stiffness used by the
  solid model's Laplacian term.

---

## Developer Notes

### Class role

`viscousHookeanElastic` derives directly from `mechanicalLaw` and is registered
in the `linGeomMechLaw` runtime-selection table. It stores `EInfinity`, the
Maxwell moduli and relaxation times, their relative moduli, `nu`, `lambda`,
`mu`, `k`, the cell and face internal-stress lists, the cell and face trial
deviatoric stresses, and the WLF switch and coefficients.

The class uses `#ifdef OPENFOAM_NOT_EXTEND` for OpenFOAM-specific field headers
and for `primitiveField()` in its residual calculation; foam-extend uses
`internalField()`. It appears in both `Make/files.openfoam` and
`Make/files.foamextend`.

### Construction

The base constructor first reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, adding their defaults to the stored dictionary,
and selects the base region from `regionName`, `solid`, or `region0`. Failure to
find a region is fatal. Density is read later by `rhoScalar()`, rather than by
the constructor.

The derived constructor then reads `EInfinity`, `E`, `relaxationTimes`, and
`nu`. It reads `WilliamsLandelFerry` with a default of `false`; when enabled,
it requires `C1`, `C2`, and `Tref` from `WilliamsLandelFerryCoeffs`. It checks
that the two lists have equal lengths, calculates `E0` and the relative moduli,
and rejects a relaxation time or Maxwell modulus below `SMALL` with a fatal
error.

It creates zero-valued `s`, `sf`, `h[i]`, and `hf[i]` fields and stores their
old-time values. Finally, it calculates `mu`, rejects `nu` outside
`[-1, 0.5]`, and calculates `lambda` and `k`. For `nu == 0.5`, it prints an
incompressibility message and sets `lambda` and `k` to `GREAT`. Construction
also prints the relative moduli and the state of WLF shifting. It emits no
warnings.

### Key methods

- `impK()` forms a time-step-dependent algorithmic stiffness. Its scale factor
  is `gammaInf + sum(gamma[i]*exp(-deltaT/(2*tau[i])))`; the returned field is
  that factor times `2*mu`, plus `lambda` unless `nu == 0.5` or a cell field
  named `p` exists. This diffusivity affects the outer-iteration convergence
  rate rather than the converged answer. WLF shifting is not included in this
  scale factor.
- `correct(volSymmTensorField&)` and its face-field counterpart update the
  trial stress, Maxwell-arm stresses, and total stress. They select `grad(D)`
  or `grad(DD)` according to the solid model's incremental setting and apply
  WLF shifting when enabled.
- `residual()` returns the largest relative previous-iteration change over all
  Maxwell-arm stresses. It uses face fields when the base mesh has an `Ff`
  surface field and cell fields otherwise.
- `bulkModulus()` delegates to `K()`, which returns a uniform field containing
  `k`.
- `materialTangent()` is not overridden. The base implementation aborts with
  `notImplemented`, so solid models requiring a material tangent cannot use
  this law.

### Extension points

A related linear viscoelastic law can copy this class and replace the
Maxwell-arm update in both `correct` overloads. New internal variables need
matching cell and face fields, old-time storage, previous-iteration storage if
they contribute to `residual()`, and corresponding terms in `impK()`. A law
needed by point- or quadrature-based solid models must also implement those
`correct` overloads, and a law needed by Newton-based solid models must provide
`materialTangent()`.

The source is at
[viscousHookeanElastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/viscousHookeanElastic/viscousHookeanElastic.C).

---

## Tutorials

- `solids/viscoelasticity/viscoTube`
