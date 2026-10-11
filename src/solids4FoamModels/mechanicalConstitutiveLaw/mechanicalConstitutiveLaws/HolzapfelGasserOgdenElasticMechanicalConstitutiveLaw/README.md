---
sort: 24
---

# HolzapfelGasserOgdenElastic

The Holzapfel-Gasser-Ogden (HGO) anisotropic hyperelastic law for
fibre-reinforced soft tissue, in particular arterial wall. The runtime type
is:

```text
HolzapfelGasserOgdenElastic
```

The constitutive model is due to G. A. Holzapfel, T. C. Gasser and
R. W. Ogden (2000).

---

## User Guide

### What it computes

This is a **finite-strain** law: it implements the finite-strain update only,
so it is used with a nonlinear geometry solid model. A linear geometry solid
model asks for the small-strain update, and the run stops with an error.

The response is a neo-Hookean ground substance plus two exponential fibre
families, wound symmetrically about the circumferential direction `Ec`:

```text
W = 0.5*mu*(I1bar - 3)
  + (k1/(2*k2))*sum_i (exp(k2*sqr(I_i - 1)) - 1)
  + U(J)

N4 =  cos(fibreAngle)*Ec + sin(fibreAngle)*Ea
N6 = -cos(fibreAngle)*Ec + sin(fibreAngle)*Ea
```

All invariants are taken on the isochoric deformation `Fbar = J^(-1/3)*F`:

```text
bbar = symm(Fbar & Fbar^T)
Cbar = symm(Fbar^T & Fbar)
I4   = Cbar && symm(N4*N4),   n4 = Fbar & N4
I6   = Cbar && symm(N6*N6),   n6 = Fbar & N6

T     = mu*bbar
      + 2*k1*(I4 - 1)*exp(k2*sqr(I4 - 1))*symm(n4*n4)
      + 2*k1*(I6 - 1)*exp(k2*sqr(I6 - 1))*symm(n6*n6)

sigma = dev(T)/J + dU/dJ*I

U(J)  = 0.25*K*(J^2 - 1 - 2*ln(J))
dU/dJ = 0.5*K*(J^2 - 1)/J
```

where `K` is `bulkModulus`, a penalty enforcing near incompressibility.

```warning
The fibre terms are active in compression as well as tension. The standard
HGO model switches them off when `I4 < 1`; this implementation does not, and
has no fibre dispersion parameter.
```

A deformation gradient with `J <= sqrt(SMALL)` stops the run.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `mu` | yes | - | Ground-substance shear modulus, pressure units, `> 0` |
| `k1` | yes | - | Fibre stress scale, pressure units |
| `k2` | yes | - | Fibre exponent; dimensioned, dimensionless |
| `fibreAngle` | yes | - | Angle from `Ec`, in **degrees**; dimensioned |
| `bulkModulus` | yes | - | Volumetric penalty `K`, pressure units, `> 0` |
| `uniformLocalBasis` | no | `no` | Default `Ec`, `Ea` from this dictionary |
| `Ec` | if uniform | - | Circumferential direction vector |
| `Ea` | if uniform | - | Axial direction vector |

`k2` and `fibreAngle` are read as dimensioned scalars, e.g.
`k2 k2 [0 0 0 0 0 0 0] 1.465;`, and only their values are used. `mu` and
`bulkModulus` must be positive; otherwise construction stops with an error.

### Local basis

`Ec` and `Ea` are prescribed state variables. When the manager sets up the
law's state it looks, for each, for a `volVectorField` of that name:

1. a file in the current (start) time directory;
2. a file in `0`;
3. a field already registered by another model.

If one is found it is used, mapped from the cells to the law's integration
points. Otherwise the dictionary vector is used when `uniformLocalBasis` is
`yes`, and a zero vector when it is not. A zero-length `Ec` or `Ea` stops the
run at the first evaluation. Both are normalised before use, so they need not
be unit vectors. A field in the time directory takes precedence over the
dictionary vector even when `uniformLocalBasis` is `yes`.

In the `ratCarotid` tutorial the fields are written by the case-local
`calcLocCoordinates` utility, which `Allrun` compiles and runs before the
solver.

### State variables

| Name | Type | Role | Description |
| --- | --- | --- | --- |
| `Ec` | vector | prescribed | Reference circumferential direction |
| `Ea` | vector | prescribed | Reference axial direction |

Prescribed variables are read once and never written back, so they are not
written with the results or to the restart files. The law has no persistent
history.

### Tangents

- `scalar`: `(4/3)*mu + K`. The fibre stiffness is not included; it is a
  constant stiffness estimate for the segregated solver, not a tangent of the
  energy.
- `scalarDeviatoric`: `(4/3)*mu`, for mixed displacement-pressure
  formulations.
- `fourthOrderFiniteDifference`: supported, through the base class.
- `fourthOrder` (analytical): not implemented; requesting it is a fatal
  error.

### Mixed displacement-pressure formulation

The law provides the volumetric split, so it can be used with a mixed
formulation such as
[`coupledPressureDisplacementSolid`](../../../solidModels/coupledPressureDisplacementSolid/README.md)
or
[`nonLinearGeometryTotalLagrangianTotalDisplacement`](../../../solidModels/nonLinGeomTotalLagTotalDispSolid/README.md)
with `solvePressure true`. The solid model asks for the isochoric stress and
`dU/dJ` separately and replaces `dU/dJ` with its solved pressure. A large
`bulkModulus` approaches the incompressible limit; `ratCarotid` uses one
thousand times `mu`. The law can also be run in the displacement-only
formulation, where the penalty alone controls the volume change.

### Differences from the legacy law

- The shear modulus is given as `mu` only. `E` and `nu` are no longer read,
  and the `nu = 0.5` requirement is gone. The given `mu` is used as it is:
  the legacy law overwrote it at start-up with a perturbation estimate.
- `bulkModulus` is new and required. The legacy law was incompressible only
  and took the spherical stress from the solid model's `p` and `pf` fields;
  this law never reads `p`, and the pressure enters through the volumetric
  split instead, so it is no longer tied to `coupledPressureDisplacementSolid`.
- The invariants are taken on `Fbar`; the legacy law assumed `J = 1`.
- Only `Ec` and `Ea` are needed, as cell fields. `Er`, `Ecf`, `Eaf` and `Erf`
  are no longer read. A uniform basis can be given in the dictionary with
  `uniformLocalBasis`.
- `impKcoeff` is no longer read, and the per-step effective shear modulus
  used for the implicit stiffness is replaced by the constant scalar
  tangents above. Unused legacy entries are ignored rather than reported.

### Example

From the `ratCarotid` tutorial:

```text
planeStress no;

mechanical
(
    pipe
    {
        type        HolzapfelGasserOgdenElastic;
        rho         rho [1 -3 0 0 0 0 0] 1200;
        mu          mu [1 -1 -2 0 0 0 0] 44.23e3;
        fibreAngle  fibreAngle [0 0 0 0 0 0 0] 39.76;
        k1          k1 [1 -1 -2 0 0 0 0] 0.206e3;
        k2          k2 [0 0 0 0 0 0 0] 1.465;
        bulkModulus bulkModulus [1 -1 -2 0 0 0 0] 44.23e6;
    }
);
```

`Ec` and `Ea` are read from the fields written by `calcLocCoordinates`. For
a uniform basis instead:

```text
        uniformLocalBasis yes;
        Ec                (1 0 0);
        Ea                (0 1 0);
```

---

## Tutorials

- [`solids/hyperelasticity/ratCarotid`](../../../../../tutorials/solids/hyperelasticity/ratCarotid):
  by default with `nonLinearGeometryTotalLagrangianTotalDisplacement` and
  `solvePressure true`; `./Allrun pressureDisplacement` selects
  `coupledPressureDisplacementSolid` (foam-extend only).

---

## References

- G. A. Holzapfel, T. C. Gasser and R. W. Ogden, _A new constitutive
  framework for arterial wall mechanics and a comparative study of material
  models_, Journal of Elasticity, 61 (2000) 1-48.
