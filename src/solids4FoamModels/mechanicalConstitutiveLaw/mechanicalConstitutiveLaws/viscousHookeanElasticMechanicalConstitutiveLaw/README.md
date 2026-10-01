---
sort: 9
---

# viscousHookeanElastic

Small-strain linear viscoelasticity: a generalised Maxwell (Prony series)
model with an elastic volumetric response and a relaxing deviatoric response.
The runtime type is:

```text
viscousHookeanElastic
```

The dictionary layout common to all laws is described in the
[material models page](../../../materialModels/README.md).

---

## User Guide

### What it computes

A relaxed spring (`EInfinity`) acts in parallel with a list of Maxwell arms
(`E[i]`, `relaxationTimes[i]`). The law first forms the instantaneous modulus
and relative moduli:

```text
E0        = EInfinity + sum(E[i])
gammaInf  = EInfinity/E0
gamma[i]  = E[i]/E0
mu        = E0/(2*(1 + nu))

lambdaInf = nu*EInfinity/((1 + nu)*(1 - 2*nu))  // plane strain / 3-D
lambdaInf = nu*EInfinity/((1 + nu)*(1 - nu))    // planeStress yes
k         = lambdaInf + (2/3)*gammaInf*mu
```

At each integration point, with `epsilon = symm(grad(D))` from the solid
model and `deltaT` the time step:

```text
s     = 2*mu*dev(epsilon)
h[i]  = exp(-deltaT/tau[i])*h[i]_0 + exp(-deltaT/(2*tau[i]))*(s - s_0)
sigma = k*tr(epsilon)*I + gammaInf*s + sum(gamma[i]*h[i])
```

following Simo and Hughes (1998), eq. 10.3.12. Only the deviatoric response
relaxes; the bulk response `k` is the long-term one and constant.

With the Williams-Landel-Ferry (WLF) shift enabled, each relaxation time is
multiplied by

```text
aT = 10^(-C1*(T - Tref)/(C2 + (T - Tref)))
```

evaluated with the temperature `T` at the integration point.

The law is small-strain only; a nonlinear geometry solid model stops with an
error.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `EInfinity` | yes | - | Relaxed modulus, `[1 -1 -2 0 0 0 0]` |
| `E` | yes | - | List of Maxwell-arm moduli, in Pa |
| `relaxationTimes` | yes | - | List of relaxation times, in s |
| `nu` | yes | - | Poisson's ratio, dimensionless |
| `WilliamsLandelFerry` | no | `no` | Shift relaxation times with `T` |
| `WilliamsLandelFerryCoeffs` | with WLF | - | Sub-dictionary of WLF data |

`WilliamsLandelFerryCoeffs` contains:

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `C1` | yes | - | First WLF coefficient, dimensionless |
| `C2` | yes | - | Second WLF coefficient, `[0 0 0 1 0 0 0]` |
| `Tref` | yes | - | Reference temperature, `[0 0 0 1 0 0 0]` |

The following are fatal errors at construction:

- `EInfinity` without pressure dimensions, or `nu` not dimensionless;
- `E` and `relaxationTimes` of different lengths, or empty;
- any `E[i]` or `relaxationTimes[i]` below `SMALL`;
- `nu` at or above `0.5` (to within `SMALL`).

`E` and `relaxationTimes` are plain scalar lists, without dimensions; the
relaxation times are in the case's time unit.

With `WilliamsLandelFerry yes`, the law asks the framework for a live
temperature input `T`. It is taken from a registered `T` field (for example
one solved by a thermal solid model), from a case directory given under
`inputCaseDirectories` in the law's dictionary, or from a `T` file in the
current or initial time directory; if none is found the run stops. See the
[framework README](../../README.md) for the lookup order.

### Tangents and volumetric split

- Scalar tangent: `factor*2*mu + lambdaInf`, and `factor*(4/3)*mu` for the
  deviatoric scalar tangent, with
  `factor = gammaInf + sum(gamma[i]*exp(-deltaT/(2*tau[i])))`. It depends on
  the time step. The WLF shift is not included in it.
- Fourth-order tangent: `fourthOrderFiniteDifference` only. Asking for the
  analytical `fourthOrder` tangent is a fatal error.
- Volumetric split: no. The law cannot be used with `solvePressureEqn` or with
  the mixed displacement-pressure formulations.

### State variables

| Name | Type | Role | Description |
| --- | --- | --- | --- |
| `s` | symmTensor | persistent | Instantaneous deviatoric stress |
| `h0`, `h1`, ... | symmTensor | persistent | Maxwell-arm internal stress |

There is one `h<i>` per entry of `E`. At each write time, for cell-centred
solid models, each persistent variable is written as a volField named after
it (`s`, `h0`, `h1`, ...) for viewing. The restart data are written alongside
as `<material>_<topology>_<variable>`, for example
`polymer_cellCentredIntegrationPointTopology_h0`, together with a
`<material>_<topology>_integrationPoints` file. To restart a run from a later
time, set `restart yes;` in the solid model's coefficients dictionary from the
start of the run.

### Example

From the `viscoTube` tutorial:

```text
planeStress     no;

mechanical
(
    polymer
    {
        type            viscousHookeanElastic;
        rho             rho [1 -3 0 0 0 0 0] 7850;
        EInfinity       EInfinity [1 -1 -2 0 0 0 0] 39.58e+9;
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

### Migrating from the legacy law

- `regionName` is no longer read. The legacy law silently ignored
  `solvePressureEqn`; the solid model now refuses it for this law.
- `nu = 0.5` is now a fatal error; the legacy law accepted it and set
  `lambda` and `k` to `GREAT`.
- The update is always written on the total displacement gradient and its
  old-time value. The legacy law had a separate incremental form, based on
  `grad(DD)` and the old-time total stress.
- The WLF temperature is gathered by the framework from any of the sources
  above; the legacy law required a `T` field at cells or at faces, matching
  the stress location.
- The scalar tangent no longer drops `lambda` when a field named `p` exists.
- The Maxwell-arm stresses are now written (`h0`, ...) and restartable; the
  legacy `h`, `hf`, `s` and `sf` fields were neither read nor written.
- The law now serves every integration-point location, and a fourth-order
  tangent is available by finite differences.

---

## Tutorials

- [viscoTube](../../../../../tutorials/solids/viscoelasticity/viscoTube/README.md),
  whose regression tests also cover restart

---

## References

- J.C. Simo and T.J.R. Hughes, _Computational Inelasticity_, Springer, 1998.
