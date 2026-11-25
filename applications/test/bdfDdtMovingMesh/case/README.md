# BDF ddt scheme on a moving mesh

`Test-bdfDdtMovingMesh` checks the `BDF` ddt scheme of orders 1 to 6
(`src/solids4FoamModels/numerics/bdfDdtScheme`) on a moving mesh. The points
of a box of 20 cells along `x` are moved every step by

```text
x = X + eps L sin(pi X/L) g(t),   g = sin(2 pi t)
```

so that every cell volume `V(t)` is known in closed form, and the field
`psi = f(t) (1 + X_c)` is set from an exact function each step. On a moving
mesh the finite volume time derivative is that of the cell integral,
`(1/V) d(V psi)/dt`, which `BDF q` must reproduce to order `q`.

`regressionTest.sh` runs the application in a copy of the case and checks,
for each order:

1. the observed order of `fvc::ddt(psi)` against the exact derivative, from
   the time steps `0.1`, `0.05`, `0.025` and `0.0125`: the order from the two
   smallest steps must be at least `q - 0.5`;
2. `fvm::ddt(psi)` evaluated at `psi` against `fvc::ddt(psi)`, and the
   geometric conservation law, `fvc::ddt(1)` against `div(meshPhi)`: both to
   `1e-9`;
3. exactness, to `1e-9`, when the motion and `f` are polynomials of degree
   `q` and `q - 1` (`-poly`);
4. with the mesh moved 2% past its position before the ddt terms of each
   step and to its position after them (`-moveAfterDdt`), as the FSI
   interfaces' final mesh correction (`additionalMeshCorrection`) moves it
   after the last ddt call of a step: the implicit form and the conservation
   law, to `1e-9`;
5. in every run, that `fvc::ddt(psi)` rebuilt from the weights returned by
   `fieldWeights(psi)` and the scheme's volume levels agrees with
   `fvc::ddt(psi)` to `1e-9`, and that the order of the final step is `q`.

The errors are relative to the largest value of the quantity compared, and
the derivative errors are measured over the second half of the run, after
the order has ramped up from Euler. The logs, and the table of errors and
observed orders in `log.checks`, are in `regressionTests/main`.

Run it with:

```bash
./regressionTest.sh
```

`tutorials/Alltest-regression` runs it with the tutorial regression tests.
The `BDF` scheme is built for OpenFOAM.com only, so the check is skipped on
the other forks.
