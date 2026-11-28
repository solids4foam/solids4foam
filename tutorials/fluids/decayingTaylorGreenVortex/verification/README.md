# decayingTaylorGreenVortex reconstruction check

This opt-in study checks the data that a high-order (k-exact) finite volume
scheme for `newtonIcoFluid` needs, before any solver change. It does not run
the solver.

A k-exact scheme treats the cell values as cell averages and reconstructs a
polynomial of degree k in each cell. The cell average differs from the
cell-centre value by a term of second order in the cell size, so second-order
initial fields or error norms would hide any gain in accuracy. The boundary
data must also be exact at the face quadrature points, not only at the face
centres.

It is separate from `regressionTest.sh`, and nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## What is checked

`checkTaylorGreenVortexReconstruction` sets U and p to the exact cell averages
of the Taylor-Green solution, computed with a cell quadrature that is exact
for polynomials of degree 6. It reconstructs both fields with the
`kExactLeastSquares` scheme in `system/reconstructionCheckDict`, for
`polynomialOrder` 1, 2 and 3, and prints the L2 and maximum errors of

- the values and gradients at the cell centres,
- the values (from both sides) and gradients at the internal-face quadrature
  points,
- the values and gradients at the boundary-face quadrature points.

The velocity walls are in the U stencils, and the
`decayingTaylorGreenVortexVelocity` condition gives the exact velocity at their
quadrature points. The pressure walls have a fixed gradient, so the pressure
reconstruction is one-sided there.

The driver runs the utility on the uniform and graded meshes of the tutorial
(levels 2 to 6, from 100 to 25,600 cells) and computes the observed order of
each L2 error, with the cell size taken as $$h \propto N^{-1/2}$$. On the two
finest meshes, the orders of the cell and internal-face errors must be at
least $$k + 1 - 0.2$$ for the values and $$k - 0.2$$ for the gradients. The
boundary-face errors are reported but not checked.

## Running the study

```bash
./Allverify
./Allverify --families uniform --levels 2 3 4
./Allclean
```

`Allverify` builds the utility and the boundary condition library, runs the
meshes in `work/`, and writes all errors to
`postProcessing/reconstruction_errors.csv`. It exits with 1 if a check fails.
The full study takes under a minute.

## Results

Observed orders of the L2 errors between the two finest meshes (6,400 and
25,600 cells), written as uniform / graded.

| k | Field | Cell value | Cell gradient | Face value | Face gradient |
| - | ----- | ---------- | ------------- | ---------- | ------------- |
| 1 | U | 2.00 / 2.00 | 1.51 / 1.59 | 2.31 / 2.10 | 1.52 / 1.55 |
| 1 | p | 2.00 / 2.00 | 1.50 / 1.62 | 2.08 / 2.03 | 1.51 / 1.60 |
| 2 | U | 3.53 / 3.67 | 1.99 / 1.98 | 2.99 / 2.99 | 2.00 / 1.98 |
| 2 | p | 4.35 / 3.92 | 2.01 / 1.98 | 3.02 / 2.98 | 2.01 / 1.98 |
| 3 | U | 3.98 / 3.98 | 3.50 / 3.67 | 4.04 / 4.00 | 3.07 / 2.99 |
| 3 | p | 4.22 / 3.99 | 3.52 / 3.58 | 4.27 / 4.03 | 3.28 / 3.02 |

Every measure reaches the k-exact order, and several exceed it. For k = 1 and
3, the cell-centre gradients converge at about k + 1/2, and for k = 1 so do
the face gradients. For k = 2, the cell-centre values converge at 3.5 to 4.4.
This is most likely because the stencils on these smooth meshes are nearly
symmetric, so part of the leading error cancels. For k = 2 and 3, the face
values and gradients are close to the plain k-exact orders.

## Stencil size

With `faceStencilExtraCells 15`, some stencils in the stretched cells near
the walls of the graded meshes are rank-deficient for `polynomialOrder 3`
(condition numbers near $$10^{16}$$ on `graded.4`), and the gradients in those
cells are wrong by many orders of magnitude while the cell values look
correct. The check uses 20 extra cells, which removes them on all levels. The
scheme does not warn about such stencils.
