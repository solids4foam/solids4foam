# circulationDdtSchemes

Third- and fourth-order backward differentiation formula (BDF) time schemes
for OpenFOAM v2412, built as `libcirculationDdtSchemes`.

| `ddtSchemes` entry | Formula (uniform step) |
|---|---|
| `BDF3` | `(11/6 f[n+1] - 3 f[n] + 3/2 f[n-1] - 1/3 f[n-2]) / dt` |
| `BDF4` | `(25/12 f[n+1] - 4 f[n] + 3 f[n-1] - 4/3 f[n-2] + 1/4 f[n-3]) / dt` |

Both are hard-coded in one shared implementation,
`BDFDdtScheme/BDFDdtScheme.{H,C}`, which differs between the two only in its
maximum order. `BDF3DdtScheme` and `BDF4DdtScheme` register the two names.
Every `ddtScheme` operation is implemented: `fvm::ddt` and `fvc::ddt` in all
their constant/variable-density forms, and the `ddtCorr` flux corrections that
`pimpleFoam`'s pressure equation uses, as the direct generalisation of the
second-order `backward` scheme's.

## Use

```
// system/controlDict
libs ( "libcirculationDdtSchemes.so" );

// system/fvSchemes
ddtSchemes
{
    default         BDF4;
}
```

## Behaviour to know about

**Start-up.** An order-q formula needs q old time levels. Each level's time
index records which step filled it, and a level created on demand shares the
index of the one above. The scheme counts the distinct levels, so a fresh
start ramps Euler -> BDF2 -> BDF3 -> BDF4 over its first steps. This is the
same fallback OpenFOAM's `backward` uses on its first step. Those start-up
steps have lower-order local errors, and they can limit the global order. To
avoid that, write the exact old-time fields into the start time (`U_0`,
`U_0_0`, `U_0_0_0`); OpenFOAM reads them back with distinct time indices, so
the run is at full order from step one. A restart finds the old-time fields
the previous run wrote and resumes at full order.

**Constant time step.** The coefficients are for a uniform step. A change of
`deltaT` during a run is detected on the step it happens and is fatal, rather
than silently losing the order. A restart with a different `deltaT` cannot be
detected; don't do it.

**Static meshes only.** `fvMesh` stores `V0` and `V00`, fewer old volumes than
these formulas need; construction on a moving mesh is fatal.

**Stability.** BDF3 and BDF4 are not A-stable (their stability regions leave
out a sector around the imaginary axis, wider for BDF4). They are standard for
diffusion-dominated, stiff problems. In strongly convection-dominated flow at
large Courant numbers they can be less robust than `backward`.

## Verification

`tutorials/testCases/windkesselHigherOrder/ddtVerification`: a one-cell ODE
`dT/dt = 2 pi cos(2 pi t)` whose only error is the ddt scheme's, over twelve
timesteps from 0.1 to 0.002:

| Scheme | Fitted order | Pairwise order at finest pair |
|---|---:|---:|
| Euler    | 1.000 | 1.000 |
| backward | 1.984 | 2.000 |
| BDF3     | 2.982 | 3.000 |
| BDF4     | 4.049 | 4.000 |

These use exact old-time fields. From a fresh start both measure 2.96: BDF3
barely changes, but BDF4 loses a full order because its Euler and BDF2
start-up steps set the global error. This is the start-up effect described
above, shown directly. (It is 3 rather than 2 here only because this solution
has T''(0) = 0.)
