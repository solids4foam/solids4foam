# `leastSquares`

Point-value least-squares reconstruction with one polynomial per cell,
selected with `type leastSquares;` in the `displacement` (or `pressure`)
sub-dictionary of `highOrderCoeffs`. It is unrelated to the OpenFOAM gradient
scheme of the same name, which lives in a different run-time selection table.

## What is stored and what is reconstructed

The stored cell unknown $u_P$ is the field value at the cell centre. For every
cell one polynomial of order $p$ (1, 2 or 3), centred at the cell centre and
with its constant fixed to $u_P$, is fitted by weighted least squares to the
point values of the stencil cells:

$$
u_P^R(\boldsymbol{x})=u_P
+\sum_{1\le|\alpha|\le p}\frac{\partial^{\alpha}u_P}{\alpha!}
\,\boldsymbol{r}^{\alpha},
\qquad
\boldsymbol{r}=\boldsymbol{x}-\boldsymbol{x}_P,
$$

so that $u_P^R(\boldsymbol{x}_P)=u_P$. Each stencil cell $N$ supplies one
equation $u_N-u_P=\sum_\alpha \boldsymbol{d}_{PN}^{\alpha}/\alpha!\;
\partial^\alpha u_P$. The cell polynomial is then reused for values and
derivatives at the cell centre, the face centres, the face quadrature points
and arbitrary points.

This differs from

- `movingLeastSquares`, which fits a new polynomial at every evaluation point
  and does not fix the constant to the cell value, and
- `kExactLeastSquares`, whose stored values are cell averages and whose basis
  monomials are volume-averaged and corrected with the cell moments.

## Relation to `kExactLeastSquares`

`leastSquares` is the k-exact formulation with every positive-order central
cell moment set to zero: the averaged monomial of a neighbour reduces to the
point monomial $\boldsymbol{d}_{PN}^{\alpha}$, the owner-cell averages of the
non-constant terms vanish, the mirrored observations on symmetry patches
become monomials of the mirrored offset, the prescribed-value boundary
observation $u_b(\boldsymbol{x}_q)-u_P$ loses its moment term, and the value
coefficients at a point lose the derivative-moment contractions. No cell
moments are computed or exchanged. Everything else, from the stencils and
weights to the face-gradient coefficients, the boundary treatment and the
processor communication, is the shared implementation in `cellLeastSquares`,
documented function by function in
[`../kExactLeastSquares/README.md`](../kExactLeastSquares/README.md).

## Dictionary entries

The entries are those of `kExactLeastSquares`:

```text
displacement
{
    type                  leastSquares;
    polynomialOrder       3;
    cellStencilExtraCells 15;   // faceStencilExtraCells is accepted as well
    weightFunctionCoeffs
    {
        type radiallySymmetricExponential;
        k    6;
    }
    calcConditionNumber   false;
}
```

`haloDepthScale` (default `2.5*polynomialOrder`) controls the parallel halo
depth of the stencil, as for the other schemes.

## Representation query

Code outside the reconstruction classes asks
`leastSquaresScheme::cellAverageUnknowns()` whether a scheme stores cell
averages; `leastSquares` and `movingLeastSquares` return `false`,
`kExactLeastSquares` returns `true`. The tutorial analytical-solution function
objects and `Test-highOrderGrad` use this to compare against point values or
cell averages of the exact solution.
