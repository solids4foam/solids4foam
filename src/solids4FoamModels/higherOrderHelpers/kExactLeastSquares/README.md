# `kExactLeastSquares`

## Class principal functions

This document provides explanation of the principal functions documented below
and their roles in the reconstruction procedure.

- `valueAtPoint()`
  - Evaluates a scalar or vector reconstruction associated with one cell at a
    requested point $\boldsymbol{x}$.
- `cellValueCoeffsAtPoint()`
  - Constructs scalar coefficients that multiply the owner and stencil-cell
    averages to obtain the reconstructed value at $\boldsymbol{x}$.
  - It converts the derivative-coefficient rows stored for cell $P$ into scalar
    value coefficients at the requested point $$\boldsymbol{x}$$.
- `cellGradCoeffsAtPoint()`
  - Constructs vector coefficients that multiply stencil-minus-owner cell
    average differences to obtain the reconstructed gradient
    at $$\boldsymbol{x}$$.
  - It converts the derivative-coefficient rows stored for cell $P$ into vector
    coefficients at the requested point $$\boldsymbol{x}$$.
- `faceGradStencil()` and `makeFaceGradStencil()`
  - Lazily provide and construct the merged global-cell addressing required
    for face-gradient reconstruction.
- `calcFaceGradCoeffs()`
  - Calculates gradient coefficients at face quadrature points for
    internal and boundary faces.
- `calcCellCoeffs()`
  - Solves the weighted k-exact least-squares systems and stores the
    cell-centred first-, second- and third-derivative coefficient rows.

### General

The derivation is presented for one scalar field component $u$. For the
displacement field, the same reconstruction and the same geometrical
coefficients are applied independently to $u_x$, $u_y$, and $u_z$.

### Stored cell value

The value stored in cell $P$ is the cell volume average

$$
\overline{u}_P
=
\frac{1}{|\Omega_P|}
\int_{\Omega_P}u(\boldsymbol{x})\,\mathrm{d}\Omega.
$$

It is generally not equal to the value at the cell centroid:

$$
\overline{u}_P \neq u(\boldsymbol{x}_P).
$$

The polynomial $u_P^R(\boldsymbol{x})$ is the reconstruction associated
with cell $P$. Introduce the coordinates relative to the centroid of $P$:

$$
\boldsymbol{r}
=
\boldsymbol{x}-\boldsymbol{x}_P,
\qquad
r_x=x-x_P,
\qquad
r_y=y-y_P,
\qquad
r_z=z-z_P.
$$

The reconstruction is constructed so that its average over $P$ is exactly
the stored value:

$$
\frac{1}{|\Omega_P|}
\int_{\Omega_P}u_P^R(\boldsymbol{x})\,\mathrm{d}\Omega
=
\overline{u}_P.
$$

### Cell moments

The volume-normalised central moments of cell $P$ are denoted by $M^P$.
For example,
$$
M_{xx}^P
=
\frac{1}{|\Omega_P|}
\int_{\Omega_P}r_x^2\,\mathrm{d}\Omega,
\qquad
M_{xy}^P
=
\frac{1}{|\Omega_P|}
\int_{\Omega_P}r_xr_y\,\mathrm{d}\Omega,
$$

and

$$
M_{xxy}^P
=
\frac{1}{|\Omega_P|}
\int_{\Omega_P}r_x^2r_y\,\mathrm{d}\Omega.
$$

The first-order moments are zero because $\boldsymbol{x}_P$ is the volume
centroid:

$$
M_x^P=M_y^P=M_z^P=0.
$$

This is an analytical property of a volume centroid. `calcCellCoeffs()` uses
this property directly and sets first-order central moments to zero. The
`cellMoments` test utility can be used to measure the numerical quadrature
error in this identity on a particular mesh.

### Explicit mean-free reconstruction

Using derivative components instead of compact second- and third-order
derivative tensor notation, the reconstruction through third order is
$$
\begin{aligned}
u_P^R(\boldsymbol{x})={}&\overline{u}_P
+u_{x,P}r_x+u_{y,P}r_y+u_{z,P}r_z
\\
&+\frac{1}{2}u_{xx,P}\left(r_x^2-M_{xx}^P\right)
+u_{xy,P}\left(r_xr_y-M_{xy}^P\right)
\\
&+u_{xz,P}\left(r_xr_z-M_{xz}^P\right)
+\frac{1}{2}u_{yy,P}\left(r_y^2-M_{yy}^P\right)
\\
&+u_{yz,P}\left(r_yr_z-M_{yz}^P\right)
+\frac{1}{2}u_{zz,P}\left(r_z^2-M_{zz}^P\right)
\\
&+\frac{1}{6}u_{xxx,P}\left(r_x^3-M_{xxx}^P\right)
+\frac{1}{2}u_{xxy,P}\left(r_x^2r_y-M_{xxy}^P\right)
\\
&+\frac{1}{2}u_{xxz,P}\left(r_x^2r_z-M_{xxz}^P\right)
+\frac{1}{2}u_{xyy,P}\left(r_xr_y^2-M_{xyy}^P\right)
\\
&+u_{xyz,P}\left(r_xr_yr_z-M_{xyz}^P\right)
+\frac{1}{2}u_{xzz,P}\left(r_xr_z^2-M_{xzz}^P\right)
\\
&+\frac{1}{6}u_{yyy,P}\left(r_y^3-M_{yyy}^P\right)
+\frac{1}{2}u_{yyz,P}\left(r_y^2r_z-M_{yyz}^P\right)
\\
&+\frac{1}{2}u_{yzz,P}\left(r_yr_z^2-M_{yzz}^P\right)
+\frac{1}{6}u_{zzz,P}\left(r_z^3-M_{zzz}^P\right).
\end{aligned}
$$

### Compact tensor form

The same reconstruction can be written more compactly as

$$
\begin{aligned}
u_P^R(\boldsymbol{x})={}&\overline{u}_P
+\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}\boldsymbol{\nabla}u_P
+\frac{1}{2}
\left[
\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}
\boldsymbol{H}_P\boldsymbol{r}
-\boldsymbol{H}_P:\boldsymbol{M}^{(2)}_P
\right]
+\frac{1}{6}
\left[
\boldsymbol{T}_P:\boldsymbol{r}^{3}
-\boldsymbol{T}_P:\boldsymbol{M}^{(3)}_P
\right].
\end{aligned}
$$

or like this:
$$
\begin{aligned}
  u_P^R(\boldsymbol{x})={}&\overline{u}_P
  +\boldsymbol r\boldsymbol{\cdot}(\nabla u)_P
 +\frac{1}{2}
  (\nabla\nabla u)_P:
  \left(
  \boldsymbol r\otimes\boldsymbol r-\boldsymbol M_P^{(2)}
  \right)
 +\frac{1}{6}
  (\nabla\nabla\nabla u)_P:
  \left(
  \boldsymbol r\otimes\boldsymbol r\otimes\boldsymbol r
  -\boldsymbol M_P^{(3)}
  \right).
  \end{aligned}
$$

The second-order derivative tensor contains all second derivatives:
$$
\boldsymbol{H}_P
=
\begin{bmatrix}
u_{xx,P} & u_{xy,P} & u_{xz,P}\\
u_{xy,P} & u_{yy,P} & u_{yz,P}\\
u_{xz,P} & u_{yz,P} & u_{zz,P}
\end{bmatrix}.
$$

The third-order derivative tensor contains

$$
\left(\boldsymbol{T}_P\right)_{ijk}=u_{ijk,P},
\qquad i,j,k\in\{x,y,z\}.
$$

The moment tensors and the compact notation $\boldsymbol{r}^{3}$ mean

$$
\boldsymbol{M}^{(2)}_P
=
\frac{1}{|\Omega_P|}
\int_{\Omega_P}
\boldsymbol{r}\otimes\boldsymbol{r}\,\mathrm{d}\Omega,
$$

$$
\boldsymbol{M}^{(3)}_P
=
\frac{1}{|\Omega_P|}
\int_{\Omega_P}
\boldsymbol{r}\otimes\boldsymbol{r}\otimes\boldsymbol{r}
\,\mathrm{d}\Omega,
\qquad
\boldsymbol{r}^{3}
=
\boldsymbol{r}\otimes\boldsymbol{r}\otimes\boldsymbol{r}.
$$

The contractions in the compact equation are

$$
\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}
\boldsymbol{H}_P\boldsymbol{r}
=
\sum_{i}\sum_{j}r_i u_{ij,P}r_j,
$$

$$
\boldsymbol{H}_P:\boldsymbol{M}^{(2)}_P
=
\sum_i\sum_j u_{ij,P}M_{ij}^P,
$$

and

$$
\boldsymbol{T}_P:\boldsymbol{r}^{3}
=
\sum_i\sum_j\sum_k u_{ijk,P}r_ir_jr_k,
\qquad
\boldsymbol{T}_P:\boldsymbol{M}^{(3)}_P
=
\sum_i\sum_j\sum_k u_{ijk,P}M_{ijk}^P.
$$

For the third-order tensors, the colon therefore denotes contraction over all
three indices. Because the derivatives are symmetric, the sums contain three
equal contributions for terms such as $u_{xxy,P}$ and six equal contributions
for $u_{xyz,P}$. After multiplication by $1/6$, their factors in the explicit
equation are consequently $1/2$ and $1$, respectively.

For $p=1$, only the first line is retained. For $p=2$, the linear and
quadratic terms are retained. For $p=3$, all terms are retained. In a 2D
case, every term containing $z$ is omitted.

Every bracket containing a second- or third-order monomial has zero average
over cell $P$. Consequently, the average of the complete reconstruction is
exactly $\overline{u}_P$.

### Stencil-coefficient notation used in this README

Uppercase $C$ denotes the geometry-dependent coefficients stored by the
code. For stencil cell $N\in\mathcal{S}_P$, the derivative coefficients are
written as
$$
C_{x,N}^P,
\quad
C_{y,N}^P,
\quad
C_{xx,N}^P,
\quad
C_{xy,N}^P,
\quad
C_{xxx,N}^P,
\quad\ldots
$$

They multiply the difference between the neighbour and owner cell averages.
For example,

$$
u_{x,P}
=
\sum_{N\in\mathcal{S}_P}
C_{x,N}^P
\left(\overline{u}_N-\overline{u}_P\right),
$$

$$
u_{xy,P}
=
\sum_{N\in\mathcal{S}_P}
C_{xy,N}^P
\left(\overline{u}_N-\overline{u}_P\right),
$$

and

$$
u_{xxx,P}
=
\sum_{N\in\mathcal{S}_P}
C_{xxx,N}^P
\left(\overline{u}_N-\overline{u}_P\right).
$$

The same rule applies to every derivative component. In the code:

- `cellGradCoeffs()` stores $C_{x,N}^P,C_{y,N}^P,C_{z,N}^P$, $$\to \mathbf{C}_x^P=(C_{x,N}^P,C_{y,N}^P,C_{z,N}^P)$$;
- `cellSecondGradCoeffs()` stores $C_{xx,N}^P,C_{xy,N}^P,\ldots,C_{zz,N}^P$,
   $$\to \mathbf{C}_{xx}^P=(C_{xx,N}^P,C_{xy,N}^P,\ldots,C_{zz,N}^P)$$;
- `cellThirdGradCoeffs()` stores $C_{xxx,N}^P,C_{xxy,N}^P,\ldots,C_{zzz,N}^P$ ,
   $$\to \mathbf{C}_{xxx}^P=(C_{xxx,N}^P,C_{xxy,N}^P,\ldots,C_{zzz,N}^P)$$.

## `valueAtPoint`

The common reconstruction interface evaluates a scalar or vector field from a
specified cell at an arbitrary point:

```cpp
reconstruction.valueAtPoint(field, cellID, x);
```

The stored value is the cell average. With the position offset

$$
\boldsymbol{r}=\boldsymbol{x}-\boldsymbol{x}_P,
$$

the mean-free basis gives

$$
\begin{aligned}
u_P^R(\boldsymbol{x})={}&\overline{u}_P
+\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}(\nabla u)_P
+\frac{1}{2}(\nabla\nabla u)_P:
\left(
\boldsymbol{r}\otimes\boldsymbol{r}
-\boldsymbol{M}^{(2)}_P
\right)
+\frac{1}{6}(\nabla\nabla\nabla u)_P:
\left(
\boldsymbol{r}\otimes\boldsymbol{r}\otimes\boldsymbol{r}
-\boldsymbol{M}^{(3)}_P
\right).
\end{aligned}
$$

The moment corrections are essential: without them, the expansion would treat
the cell average as though it were a point value at the centroid. For a vector
field, the same equation is applied independently to each component.

To evaluate this equation, the scalar and vector `valueAtPoint()` overloads
call the templated `evaluateAtPoint()` function. `evaluateAtPoint()` then calls
`cellValueCoeffsAtPoint(cellID, x, valueCoeffs)`, which combines the first-,
second- and third-derivative coefficient rows and the cell moments into scalar
value coefficients at the requested point. `evaluateAtPoint()` multiplies
these coefficients by the owner value and the stencil-cell values to obtain
$u_P^R(\boldsymbol{x})$.

`valueAtPoint()` does not construct `grad()`, `secondGrad()` or `thirdGrad()`
field objects over the complete mesh. However, on first use, it may construct
and cache the geometry-dependent coefficient and moment tables for all local
cells. Once these mesh-wide tables are available, evaluation uses only the
coefficient rows associated with the selected cell. The same calculation works
for scalar and vector fields.

In a parallel run, the operation must be called by every processor because the
selected cell stencil may contain remote cells. The processor that owns the
requested cell passes its local `cellID`; all other processors pass `-1` and
receive zero. All processors participate in the remote-field exchange, but
only the owning processor evaluates a cell reconstruction. The resulting value
can then be combined with a global sum. `solidPointDisplacement` follows this
pattern.

## `cellValueCoeffsAtPoint`

```cpp
void kExactLeastSquares::cellValueCoeffsAtPoint
(
    const label cellI,
    const point& x,
    UList<scalar>& coeffs
) const
```

This function does not use field values. It converts the derivative-coefficient
rows stored for cell $P$ into scalar value coefficients at the requested point $$\boldsymbol{x}$$.

Let

$$
\boldsymbol{r}=\boldsymbol{x}-\boldsymbol{x}_P.
$$

For stencil cell $N\in\mathcal{S}_P$, the stored derivative coefficients are:

- $\boldsymbol{C}_{x,N}^P$, the first-derivative coefficient vector;
- $\boldsymbol{C}_{xx,N}^P$, the second-derivative coefficient tensor;
- $\boldsymbol{C}_{xxx,N}^P$, the third-derivative coefficient tensor.

For example,

$$
\boldsymbol{C}_{x,N}^P
=
\left(
C_{x,N}^P,
C_{y,N}^P,
C_{z,N}^P
\right).
$$

### Value coefficient for stencil cell \(N\)

The scalar coefficient multiplying the value stored in stencil cell $N$ is

$$
\begin{aligned}
C_N^P(\boldsymbol{x})={}&
\boldsymbol{r}\mathbin{\boldsymbol{\cdot}}\boldsymbol{C}_{x,N}^P
+\frac{1}{2}\boldsymbol{C}_{xx,N}^P:
\left(
\boldsymbol{r}\otimes\boldsymbol{r}
-\boldsymbol{M}^{(2)}_P
\right)
+\frac{1}{6}\boldsymbol{C}_{xxx,N}^P:
\left(
\boldsymbol{r}\otimes\boldsymbol{r}\otimes\boldsymbol{r}
-\boldsymbol{M}^{(3)}_P
\right).
\end{aligned}
$$

In the implementation:

- `gradCoeffs[stencilI]` is $\boldsymbol{C}_{x,N}^P$;
- `secondDerivativeCoeff` is $\boldsymbol{C}_{xx,N}^P$;
- `thirdDerivativeCoeff` is $\boldsymbol{C}_{xxx,N}^P$;
- `coeffs[stencilI]` receives the resulting $C_N^P(\boldsymbol{x})$.

The moment-contraction helper functions calculate only the moment parts of the
second- and third-order contributions. The corresponding point-dependent parts
are calculated by `r & (secondDerivativeCoeff & r)` and
`cubicForm(thirdDerivativeCoeff, r)`.

### Difference form and owner coefficient

The derivative coefficient rows multiply differences between stencil and owner
cell averages. Substitution into the mean-free polynomial therefore gives

$$
u_P^R(\boldsymbol{x})
=
\overline{u}_P
+
\sum_{N\in\mathcal{S}_P}
C_N^P(\boldsymbol{x})
\left(\overline{u}_N-\overline{u}_P\right).
$$

Expanding the differences gives

$$
u_P^R(\boldsymbol{x})
=
\left(
1-\sum_{N\in\mathcal{S}_P}C_N^P(\boldsymbol{x})
\right)
\overline{u}_P
+
\sum_{N\in\mathcal{S}_P}
C_N^P(\boldsymbol{x})\overline{u}_N.
$$

The owner coefficient is consequently

$$
C_P^P(\boldsymbol{x})
=
1-\sum_{N\in\mathcal{S}_P}C_N^P(\boldsymbol{x}).
$$

The absolute-value reconstruction is

$$
u_P^R(\boldsymbol{x})
=
C_P^P(\boldsymbol{x})\overline{u}_P
+
\sum_{N\in\mathcal{S}_P}
C_N^P(\boldsymbol{x})\overline{u}_N.
$$

### Output ordering

The caller must provide `coeffs` with

$$
\texttt{stencil().cellsStencil()[cellI].size()}+1
$$

entries. If the stencil ordering is
$\mathcal{S}_P=[N_1,N_2,\ldots,N_n]$, the output ordering is

$$
\left[
C_{N_1}^P(\boldsymbol{x}),
C_{N_2}^P(\boldsymbol{x}),
\ldots,
C_{N_n}^P(\boldsymbol{x}),
C_P^P(\boldsymbol{x})
\right].
$$

Thus, the first $n$ entries retain the stencil ordering and the last entry is
the owner coefficient. All coefficients are scalar, so the same list
reconstructs a scalar field or each component of a vector field.

The absolute coefficients satisfy

$$
C_P^P(\boldsymbol{x})
+
\sum_{N\in\mathcal{S}_P}C_N^P(\boldsymbol{x})
=1.
$$

Consequently, a constant field is reconstructed exactly up to round-off. The
individual stencil coefficients $C_N^P(\boldsymbol{x})$ are not generally zero;
constant preservation follows from the sum of all absolute coefficients being
one, or equivalently from
$\overline{u}_N-\overline{u}_P=0$ in the difference form.

### How the coefficients are used

Once the required local and remote field values are available, evaluation can
be written schematically as

```cpp
Type value = coeffs[cellStencil.size()]*ownerFieldValue;

forAll(cellStencil, stencilI)
{
    value += coeffs[stencilI]*stencilFieldValue[stencilI];
}
```

---

## `cellGradCoeffsAtPoint`

```cpp
void kExactLeastSquares::cellGradCoeffsAtPoint
(
    const label cellI,
    const point& x,
    UList<vector>& coeffs
) const;
```

This function does not use field values. It converts the derivative-coefficient
rows stored for cell $P$ into vector gradient coefficients at the requested
point $\boldsymbol{x}$.

Let

$$
\boldsymbol{r}=\boldsymbol{x}-\boldsymbol{x}_P.
$$

Differentiating the cell reconstruction with respect to
$\boldsymbol{x}$ gives

$$
\begin{aligned}
\nabla u_P^R(\boldsymbol{x})={}&(\nabla u)_P
+(\nabla\nabla u)_P\mathbin{\boldsymbol{\cdot}}\boldsymbol{r}
+\frac{1}{2}
(\nabla\nabla\nabla u)_P:
\left(\boldsymbol{r}\otimes\boldsymbol{r}\right).
\end{aligned}
$$

The moment corrections do not appear in this equation because
$\boldsymbol{M}_P^{(2)}$ and $\boldsymbol{M}_P^{(3)}$ are constant for cell
$P$. Their derivatives with respect to $\boldsymbol{x}$ are zero.

### Gradient coefficient for stencil cell \(N\)

For each stencil cell $N\in\mathcal{S}_P$, the function calculates

$$
\begin{aligned}
\boldsymbol{C}_N^P(\boldsymbol{x})={}&
\boldsymbol{C}_{x,N}^P
+\boldsymbol{C}_{xx,N}^P
\mathbin{\boldsymbol{\cdot}}\boldsymbol{r}
+\frac{1}{2}
\boldsymbol{C}_{xxx,N}^P:
\left(\boldsymbol{r}\otimes\boldsymbol{r}\right).
\end{aligned}
$$

Here:

- $\boldsymbol{C}_{x,N}^P$ is the first-derivative coefficient vector;
- $\boldsymbol{C}_{xx,N}^P$ is the second-derivative coefficient tensor;
- $\boldsymbol{C}_{xxx,N}^P$ is the third-derivative coefficient tensor;
- $\boldsymbol{C}_N^P(\boldsymbol{x})$ is the resulting vector coefficient
  for the gradient at $\boldsymbol{x}$.

In the implementation:

- `gradCoeffs[cI]` is $\boldsymbol{C}_{x,N}^P$;
- `(*secondGradCoeffsPtr)[cellI][cI] & r` evaluates $\boldsymbol{C}_{xx,N}^P\boldsymbol{\cdot}\boldsymbol{r}$;
- `0.5*(((*thirdGradCoeffsPtr)[cellI][cI] & r) & r)` evaluates $\frac{1}{2}\boldsymbol{C}_{xxx,N}^P:
  (\boldsymbol{r}\otimes\boldsymbol{r})$;
- `coeffs[cI]` receives $\boldsymbol{C}_N^P(\boldsymbol{x})$.

For $p=1$, only the first-derivative coefficient is used, so the reconstructed
gradient is constant inside cell $P$. For $p=2$, the second-derivative term
makes the gradient linear in $\boldsymbol{r}$. For $p=3$, the third-derivative
term adds quadratic variation.

### Difference form

The returned coefficients retain the neighbour-minus-owner form:

$$
\nabla u_P^R(\boldsymbol{x})
=
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_N^P(\boldsymbol{x})
\left(
\overline{u}_N-\overline{u}_P
\right).
$$

Unlike `cellValueCoeffsAtPoint()`, this function does not append an owner
coefficient. It returns one vector coefficient for every cell in
`stencil().cellsStencil()[cellI]`. If an absolute-value form is required, the
owner gradient coefficient is

$$
\boldsymbol{C}_P^P(\boldsymbol{x})
=
-\sum_{N\in\mathcal{S}_P}\boldsymbol{C}_N^P(\boldsymbol{x}).
$$

The absolute gradient coefficients consequently sum to the zero vector:

$$
\boldsymbol{C}_P^P(\boldsymbol{x})
+
\sum_{N\in\mathcal{S}_P}\boldsymbol{C}_N^P(\boldsymbol{x})
=\boldsymbol{0}.
$$

This ensures that the reconstructed gradient of a constant field is zero.

#### Output ordering for newly calculated coefficients

The caller must provide `coeffs` with exactly

$$
\texttt{stencil().cellsStencil()[cellI].size()}
$$

entries. If the stencil ordering is
$\mathcal{S}_P=[N_1,N_2,\ldots,N_n]$, the output ordering is

$$
\left[
\boldsymbol{C}_{N_1}^P(\boldsymbol{x}),
\boldsymbol{C}_{N_2}^P(\boldsymbol{x}),
\ldots,
\boldsymbol{C}_{N_n}^P(\boldsymbol{x})
\right].
$$

There is no additional owner entry.

#### Scalar and vector fields

For a scalar field, each returned vector coefficient is multiplied by a
cell-average scalar difference:

$$
\nabla u_P^R(\boldsymbol{x})
=
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_N^P(\boldsymbol{x})
\left(\overline{u}_N-\overline{u}_P\right).
$$

For displacement, the same operation is applied component-wise and produces
the displacement-gradient tensor:
$$
\nabla\boldsymbol{D}_P^R(\boldsymbol{x})
=
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_N^P(\boldsymbol{x})
\otimes
\left(
\overline{\boldsymbol{D}}_N-\overline{\boldsymbol{D}}_P
\right).
$$

Although the helper accepts any point, it is intended for points inside or
close to owner cell $P$, particularly face quadrature points.

## `faceGradStencil`

The entries are global cell IDs. This is necessary because a cell
reconstruction stencil can contain cells owned by another processor.

### Grad stencil for internal faces

For an internal face $f$ with owner cell $P$ and neighbour cell $N$,
`makeFaceGradStencil()` constructs

$$
\mathcal{S}_f
=
\{P\}
\cup\mathcal{S}_P
\cup\{N\}
\cup\mathcal{S}_N.
$$

The global IDs of $P$ and $N$ are inserted explicitly, followed by all entries
from the two cell stencils. A `labelHashSet` removes duplicate global IDs.
The resulting IDs are sorted before storage, so the face-stencil ordering is
deterministic but does not separately identify owner-side and neighbour-side
entries.

### Fixed-value boundary faces

For a boundary patch selected by `includePatchInStencils_`, the face stencil
contains the owner cell and its reconstruction stencil:

$$
\mathcal{S}_f
=
\{P\}\cup\mathcal{S}_P.
$$

Prescribed values at boundary quadrature points are known data rather than cell
unknowns. They therefore use the separate `faceBoundaryDataStencil()`
addressing and are not inserted into `faceGradStencil()`.

The face stencil contains addressing only. The corresponding vector
coefficients at every face quadrature point are calculated later by
`calcFaceGradCoeffs()` and stored in `faceGradCoeffsPtr_`.

## `calcCellCoeffs`

`calcCellCoeffs()` is the workhorse that constructs the geometry-dependent
coefficients used to recover the first, second and third derivatives in every
cell. It does not use values of the reconstructed field. Once calculated, the
same coefficients can be applied to any scalar field or independently to every
component of a vector field.

```cpp
void kExactLeastSquares::calcCellCoeffs() const;
```

The function fills:

- `cellGradCoeffsPtr_` for every polynomial order;
- `cellSecondGradCoeffsPtr_` when $p\ge2$;
- `cellThirdGradCoeffsPtr_` when $p\ge3$.

Each stored row follows the ordering of
`stencil().cellsStencil()[cellI]`. The owner cell is not appended because the
derivative reconstruction uses stencil-minus-owner cell-average differences.

### Polynomial basis

A polynomial basis is simply the list of polynomial terms that the
reconstruction is allowed to contain. Instead of writing a separate case for
every term, the code describes each term using three integer powers.

For a multi-index

$$
\boldsymbol{\alpha}=(i,j,k),
\qquad
|\boldsymbol{\alpha}|=i+j+k,
\qquad
\boldsymbol{\alpha}!=i!j!k!,
$$

define

$$
\boldsymbol{r}^{\boldsymbol{\alpha}}
=r_x^i r_y^j r_z^k.
$$

#### How to read the multi-index notation

The multi-index is only a compact bookkeeping notation. The three numbers
$(i,j,k)$ are the powers of $r_x$, $r_y$ and $r_z$. The same three numbers
also identify the derivative associated with that polynomial term:

$$
D^{(i,j,k)}u
=
\frac{
\partial^{i+j+k}u
}{
\partial x^i\partial y^j\partial z^k
}.
$$

Examples are:

| $(i,j,k)$ | Polynomial term | Derivative | $i!j!k!$ |
|---:|---:|---:|---:|
| $(1,0,0)$ | $r_x$ | $u_x$ | 1 |
| $(2,0,0)$ | $r_x^2$ | $u_{xx}$ | 2 |
| $(1,1,0)$ | $r_xr_y$ | $u_{xy}$ | 1 |
| $(2,1,0)$ | $r_x^2r_y$ | $u_{xxy}$ | 2 |
| $(1,1,1)$ | $r_xr_yr_z$ | $u_{xyz}$ | 1 |

For example, $\boldsymbol{\alpha}=(2,1,0)$ simply means:

$$
\boldsymbol{r}^{\boldsymbol{\alpha}}=r_x^2r_y,
\qquad
D^{\boldsymbol{\alpha}}u=u_{xxy},
\qquad
\boldsymbol{\alpha}!=2!1!0!=2.
$$

In practical terms, the polynomial basis is just the following list of terms.
In 2D:

- $p=1$: $r_x,\ r_y$;
- $p=2$: add $r_x^2,\ r_xr_y,\ r_y^2$;
- $p=3$: add $r_x^3,\ r_x^2r_y,\ r_xr_y^2,\ r_y^3$.

In 3D:

- $p=1$: $r_x,\ r_y,\ r_z$;
- $p=2$: add $r_x^2,\ r_xr_y,\ r_xr_z,\ r_y^2,\ r_yr_z,\ r_z^2$;
- $p=3$: add $r_x^3,\ r_x^2r_y,\ r_x^2r_z,\ r_xr_y^2,\$
  $r_xr_yr_z,r_xr_z^2,\ r_y^3,\ r_y^2r_z,\ r_yr_z^2,\ r_z^3$.

These lists show only the monomial part of each basis function. The
second- and third-order cell moments are still subtracted to make the basis
mean-free over owner cell $P$.

`generateExponents()` generates all non-constant monomials satisfying

$$
1\le|\boldsymbol{\alpha}|\le p.
$$

Each generated $(i,j,k)$ identifies one row of the least-squares matrix and
one derivative that will be reconstructed.

Terms are grouped by total polynomial degree. In 2D, $k=0$ and all terms
containing $z$ are omitted. The number of non-constant terms, returned by
`minNn()`, is

$$
N_p^{2D}
=
\frac{(p+1)(p+2)}{2}-1,
$$

and

$$
N_p^{3D}
=
\frac{(p+1)(p+2)(p+3)}{6}-1.
$$

| Order | Number of 2D terms | Number of 3D terms |
|---:|---:|---:|
| $p=1$ | 2 | 3 |
| $p=2$ | 5 | 9 |
| $p=3$ | 9 | 19 |

The constant monomial is excluded. Its coefficient is fixed by the requirement
that the reconstruction average over owner cell $P$ equals
$\overline{u}_P$.

`rowOf()` returns the matrix row associated with a requested exponent
$(i,j,k)$. `calcDerivativeRows()` records the rows required to construct the
second- and third-derivative tensors. Consequently, the code does not depend on
manually assumed row positions.

#### Cell-average equations

The mean-free reconstruction associated with owner cell $P$ can be written as

$$
u_P^R(\boldsymbol{x})
=
\overline{u}_P
+
\sum_{\boldsymbol{\alpha}}
\frac{
\left(D^{\boldsymbol{\alpha}}u\right)_P
}{
\boldsymbol{\alpha}!
}
\left[
\boldsymbol{r}^{\boldsymbol{\alpha}}
-M_P^{\boldsymbol{\alpha}}
\right],
$$

where

$$
\boldsymbol{r}=\boldsymbol{x}-\boldsymbol{x}_P
$$

and

$$
M_P^{\boldsymbol{\alpha}}
=
\frac{1}{|\Omega_P|}
\int_{\Omega_P}
\boldsymbol{r}^{\boldsymbol{\alpha}}\,\mathrm{d}\Omega.
$$

For stencil cell $N$, take the volume average of this reconstruction over
$\Omega_N$ and subtract the owner-cell average:

$$
\overline{u}_N-\overline{u}_P
=
\sum_{\boldsymbol{\alpha}}
\frac{
\left(D^{\boldsymbol{\alpha}}u\right)_P
}{
\boldsymbol{\alpha}!
}
\left[
\left\langle
\boldsymbol{r}^{\boldsymbol{\alpha}}
\right\rangle_N
-M_P^{\boldsymbol{\alpha}}
\right].
$$

The notation

$$
\left\langle
\boldsymbol{r}^{\boldsymbol{\alpha}}
\right\rangle_N
=
\frac{1}{|\Omega_N|}
\int_{\Omega_N}
\boldsymbol{r}^{\boldsymbol{\alpha}}\,\mathrm{d}\Omega
$$

denotes an average over stencil cell $N$, but $\boldsymbol{r}$ is still
measured from the owner centroid $\boldsymbol{x}_P$. This distinction is why
both the neighbour-cell centre offset and the neighbour central moments are
needed.

#### The averageMonomial calculation

Let

$$
\boldsymbol{d}_{PN}
=
\boldsymbol{x}_N-\boldsymbol{x}_P,
\qquad
\boldsymbol{\xi}_N
=
\boldsymbol{x}-\boldsymbol{x}_N.
$$

Inside cell $N$,

$$
\boldsymbol{r}
=
\boldsymbol{d}_{PN}+\boldsymbol{\xi}_N.
$$

This change of origin is the exact place where the binomial expansion is
needed. The reconstruction polynomial uses coordinates measured from
$\boldsymbol{x}_P$, but the moments of cell $N$ are stored using coordinates
measured from $\boldsymbol{x}_N$. The binomial expansion converts between
these two coordinate systems.

For a linear term,

$$
\begin{aligned}
\left\langle r_x\right\rangle_N
&=
\left\langle d_{PN,x}+\xi_{N,x}\right\rangle_N
=
d_{PN,x}+\left\langle\xi_{N,x}\right\rangle_N
=d_{PN,x},
\end{aligned}
$$

because the first central moment $\left\langle\xi_{N,x}\right\rangle_N$ is zero.

For a pure quadratic term,

$$
\begin{aligned}
\left\langle r_x^2\right\rangle_N
&=
\left\langle
\left(d_{PN,x}+\xi_{N,x}\right)^2
\right\rangle_N
=
d_{PN,x}^2
+2d_{PN,x}\left\langle\xi_{N,x}\right\rangle_N
+\left\langle\xi_{N,x}^2\right\rangle_N
=
d_{PN,x}^2+M_{xx}^N.
\end{aligned}
$$

For a mixed quadratic term,

$$
\begin{aligned}
\left\langle r_xr_y\right\rangle_N
=
\left\langle
\left(d_{PN,x}+\xi_{N,x}\right)
\left(d_{PN,y}+\xi_{N,y}\right)
\right\rangle_N
=
d_{PN,x}d_{PN,y}+M_{xy}^N.
\end{aligned}
$$

A cubic example shows why more combinations appear:

$$
\begin{aligned}
\left\langle r_x^2r_y\right\rangle_N
={}
d_{PN,x}^2d_{PN,y}
+2d_{PN,x}M_{xy}^N
+d_{PN,y}M_{xx}^N
+M_{xxy}^N.
\end{aligned}
$$

The code must perform the same expansion for every combination of powers
$(i,j,k)$. This is why `averageMonomial()` contains one loop for each
coordinate direction.

Applying the binomial expansion independently in the $x$, $y$ and $z$
directions gives

$$
\begin{aligned}
\left\langle r_x^i r_y^j r_z^k\right\rangle_N
={}&
\sum_{a=0}^{i}
\sum_{b=0}^{j}
\sum_{c=0}^{k}
\binom{i}{a}
\binom{j}{b}
\binom{k}{c}
\\
&\quad
d_{PN,x}^{\,i-a}
d_{PN,y}^{\,j-b}
d_{PN,z}^{\,k-c}
M_N^{(a,b,c)}.
\end{aligned}
$$

This equation explains the three nested loops in `averageMonomial()`:

- index $a$ expands the power in the $x$ direction;
- index $b$ expands the power in the $y$ direction;
- index $c$ expands the power in the $z$ direction.

For each stencil cell and each polynomial term, `calcCellCoeffs()`
calls `averageMonomial()` twice:

1. with cell $N$ and $\boldsymbol{d}_{PN}$ to calculate the neighbour-cell
   average;
2. with cell $P$ and $\boldsymbol{d}_{PP}=\boldsymbol{0}$ to calculate the
   owner-cell average.

Their difference becomes one entry of the matrix $\boldsymbol{Q}_P$ after
division by the distance scaling and factorial.

#### Complete 2D quadratic example

For $p=2$ in 2D, the cell-average equation associated with one stencil cell
$N$ is

$$
\begin{aligned}
\overline{u}_N-\overline{u}_P
={}&
u_{x,P}d_{PN,x}
+u_{y,P}d_{PN,y}
\\
&+\frac{1}{2}u_{xx,P}
\left(
d_{PN,x}^2+M_{xx}^N-M_{xx}^P
\right)
\\
&+u_{xy,P}
\left(
d_{PN,x}d_{PN,y}+M_{xy}^N-M_{xy}^P
\right)
\\
&+\frac{1}{2}u_{yy,P}
\left(
d_{PN,y}^2+M_{yy}^N-M_{yy}^P
\right).
\end{aligned}
$$

This is the explicit, non-multi-index form of one least-squares equation.
The corresponding column of the matrix $\boldsymbol{Q}_P$ is

$$
\boldsymbol{Q}_{P,\cdot N}
=
\begin{bmatrix}
d_{PN,x}/h_P
\\
d_{PN,y}/h_P
\\
\left(d_{PN,x}^2+M_{xx}^N-M_{xx}^P\right)/(2h_P^2)
\\
\left(d_{PN,x}d_{PN,y}+M_{xy}^N-M_{xy}^P\right)/h_P^2
\\
\left(d_{PN,y}^2+M_{yy}^N-M_{yy}^P\right)/(2h_P^2)
\end{bmatrix}.
$$

It multiplies the scaled derivative vector

$$
\boldsymbol{s}_P
=
\begin{bmatrix}
h_Pu_{x,P}
&
h_Pu_{y,P}
&
h_P^2u_{xx,P}
&
h_P^2u_{xy,P}
&
h_P^2u_{yy,P}
\end{bmatrix}^{T}.
$$

The product

$$
\boldsymbol{Q}_{P,\cdot N}^{T}\boldsymbol{s}_P
$$

reproduces the complete equation for
$\overline{u}_N-\overline{u}_P$ above. Every additional stencil cell supplies
another column of $\boldsymbol{Q}_P$, or equivalently another equation in the
least-squares system. For $p=3$, the procedure is identical, but the cubic
terms additionally use the third-order moments.

For the owner cell, $\boldsymbol{d}_{PP}=\boldsymbol{0}$, so the same function
returns $M_P^{(i,j,k)}$. First-order central moments are set to zero because
cell centres are volume centroids. Second- and third-order moments are read
from `fvMeshQuadrature` and exchanged when a stencil contains remote cells.

#### Distance scaling and matrix \(Q\)

For each owner cell, the code defines

$$
h_P
=
2\max_{N\in\mathcal{S}_P}
\left|
\boldsymbol{x}_N-\boldsymbol{x}_P
\right|.
$$

The derivatives are scaled as

$$
s_{\boldsymbol{\alpha}}^P
=
h_P^{|\boldsymbol{\alpha}|}
\left(D^{\boldsymbol{\alpha}}u\right)_P.
$$

For exponent row $\boldsymbol{\alpha}=(i,j,k)$ and stencil column $N$, the
matrix assembled in the code is

$$
Q_{\boldsymbol{\alpha}N}
=
\frac{
\left\langle
\boldsymbol{r}^{\boldsymbol{\alpha}}
\right\rangle_N
-M_P^{\boldsymbol{\alpha}}
}{
h_P^{|\boldsymbol{\alpha}|}\boldsymbol{\alpha}!
}.
$$

Collecting all stencil-cell equations gives

$$
\Delta\overline{\boldsymbol{u}}_P
=
\boldsymbol{Q}_P^T\boldsymbol{s}_P,
$$

where

$$
\left(
\Delta\overline{\boldsymbol{u}}_P
\right)_N
=
\overline{u}_N-\overline{u}_P.
$$

The code stores $\boldsymbol{Q}_P$ with $N_p$ rows and
$N_n=|\mathcal{S}_P|$ columns. It requires

$$
N_n\ge N_p.
$$

The factor $h_P^{|\boldsymbol{\alpha}|}$ makes the matrix rows associated with
different derivative orders dimensionally comparable. The derivatives are
unscaled again before their coefficients are stored.

#### Weights and QR solution

The diagonal weight matrix contains one spatial weight for every stencil cell:

$$
\left(\boldsymbol{W}_P\right)_{NN}
=
w
\left(
\left|\boldsymbol{x}_N-\boldsymbol{x}_P\right|,
\max_{K\in\mathcal{S}_P}
\left|\boldsymbol{x}_K-\boldsymbol{x}_P\right|
\right).
$$

The scaled derivatives are obtained from the weighted least-squares problem

$$
\min_{\boldsymbol{s}_P}
\left\|
\boldsymbol{W}_P^{1/2}
\left(
\boldsymbol{Q}_P^T\boldsymbol{s}_P
-\Delta\overline{\boldsymbol{u}}_P
\right)
\right\|_2.
$$

`QRSolve()` returns the matrix $\boldsymbol{A}_P$ such that

$$
\boldsymbol{s}_P
=
\boldsymbol{A}_P
\Delta\overline{\boldsymbol{u}}_P.
$$

It performs a Householder QR decomposition of

$$
\boldsymbol{W}_P^{1/2}\boldsymbol{Q}_P^T
$$

and solves the resulting reduced square system. The normal-equation matrix is
not formed explicitly. If `calcConditionNumber` is enabled, the condition
number of the square upper factor is calculated and written to the
`cellConditionNumber` field.

#### Unscaling and storing derivative coefficients

For derivative $\boldsymbol{\alpha}$, row $\boldsymbol{\alpha}$ of
$\boldsymbol{A}_P$ gives

$$
\left(D^{\boldsymbol{\alpha}}u\right)_P
=
\sum_{N\in\mathcal{S}_P}
\frac{
\left(\boldsymbol{A}_P\right)_{\boldsymbol{\alpha}N}
}{
h_P^{|\boldsymbol{\alpha}|}
}
\left(
\overline{u}_N-\overline{u}_P
\right).
$$

Therefore, the stored coefficient is

$$
C_{\boldsymbol{\alpha},N}^P
=
\frac{
\left(\boldsymbol{A}_P\right)_{\boldsymbol{\alpha}N}
}{
h_P^{|\boldsymbol{\alpha}|}
}.
$$

The first-derivative rows are divided by $h_P$ and stored as

$$
\boldsymbol{C}_{x,N}^P
=
\left(
C_{x,N}^P,
C_{y,N}^P,
C_{z,N}^P
\right)
$$

in `cellGradCoeffsPtr_`.

For $p\ge2$, the second-derivative rows are divided by $h_P^2$ and stored as
the `symmTensor` $\boldsymbol{C}_{xx,N}^P$ in `cellSecondGradCoeffsPtr_`.

For $p\ge3$, the third-derivative rows are divided by $h_P^3$ and stored as
the `symmTensor3rdOrder` $\boldsymbol{C}_{xxx,N}^P$ in `cellThirdGradCoeffsPtr_`.

In 2D, every stored component containing a $z$ derivative is explicitly set
to zero. Each coefficient table has one row per local owner cell, and each row
has exactly the size and ordering of that cell's stencil. No additional owner
entry is stored.

Prescribed boundary quadrature-point values are not used by
`calcCellCoeffs()`. This function constructs the unconstrained cell-average
reconstruction coefficients. Fixed-value boundary data are introduced later
when `calcFaceGradCoeffs()` constructs the one-sided boundary-face gradient
operators.

## `calcFaceGradCoeffs`

`calcFaceGradCoeffs()` constructs the geometry-dependent operators used by
`fGrad()` to evaluate gradients at face quadrature points.

```cpp
void kExactLeastSquares::calcFaceGradCoeffs() const;
```

The function does not use values of the reconstructed field. It constructs and
caches:

- `faceGradCoeffsPtr_`, containing coefficients for unknown cell averages;
- `faceBoundaryDataStencilPtr_`, containing addresses of prescribed boundary
  quadrature-point values;
- `faceBoundaryDataCoeffsPtr_`, containing coefficients multiplying those
  prescribed values.

Internal face and fixed-value boundary faces require different constructions.

### Internal faces

Consider an internal face $f$ with owner cell $P$, neighbour cell $R$, and
quadrature point $\boldsymbol{x}_{f,q}$. A generic cell in either
reconstruction stencil is denoted by $N$.

The owner-side reconstruction gives

$$
\nabla u_P^R(\boldsymbol{x}_{f,q})
=
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_N^P(\boldsymbol{x}_{f,q})
\left(
\overline{u}_N-\overline{u}_P
\right).
$$

The neighbour-side reconstruction gives

$$
\nabla u_R^R(\boldsymbol{x}_{f,q})
=
\sum_{N\in\mathcal{S}_R}
\boldsymbol{C}_N^R(\boldsymbol{x}_{f,q})
\left(
\overline{u}_N-\overline{u}_R
\right).
$$

The face gradient is the arithmetic average of these two reconstructions:

$$
\nabla u_f(\boldsymbol{x}_{f,q})
=
\frac{1}{2}
\nabla u_P^R(\boldsymbol{x}_{f,q})
+
\frac{1}{2}
\nabla u_R^R(\boldsymbol{x}_{f,q}).
$$

`cellGradCoeffsAtPoint()` calculates $$\boldsymbol{C}_N^P(\boldsymbol{x}_{f,q})$$
and $$\boldsymbol{C}_N^R(\boldsymbol{x}_{f,q})$$.

#### Converting the internal-face differences to absolute values

For the owner side,

$$
\begin{aligned}
\frac{1}{2}\nabla u_P^R(\boldsymbol{x}_{f,q})
={}
\frac{1}{2}
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_N^P(\boldsymbol{x}_{f,q})
\overline{u}_N
-
\frac{1}{2}
\left[
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_N^P(\boldsymbol{x}_{f,q})
\right]
\overline{u}_P.
\end{aligned}
$$

The code first inserts the positive stencil-cell contributions:

```cpp
coeffs[faceStencilI] += 0.5*coeff;
ownSum += coeff;
```

It then inserts the coefficient multiplying the owner value:

```cpp
coeffs[ownIndex] -= 0.5*ownSum;
```

The same operation is performed for the neighbour-side reconstruction:

```cpp
coeffs[faceStencilI] += 0.5*coeff;
neiSum += coeff;

coeffs[neiIndex] -= 0.5*neiSum;
```

The factor `0.5` comes directly from the arithmetic average of the two
reconstructed gradients.

All contributions are accumulated into the merged face stencil

$$
\mathcal{S}_f
=
\{P\}
\cup\mathcal{S}_P
\cup\{R\}
\cup\mathcal{S}_R.
$$

If a global cell appears in both cell stencils, its contributions are added to
the same face-stencil entry.

After assembly, the face gradient has the absolute-value form

$$
\nabla u_f(\boldsymbol{x}_{f,q})
=
\sum_{N\in\mathcal{S}_f}
\boldsymbol{C}_N^{f,q}\overline{u}_N.
$$

The stored coefficients satisfy

$$
\sum_{N\in\mathcal{S}_f}
\boldsymbol{C}_N^{f,q}
=\boldsymbol{0},
$$

so the face gradient of a constant field is zero.

#### Internal-face storage

For each internal face:

- `faceGradStencil()[faceI][stencilI]` stores the global cell ID;
- `faceGradCoeffs()[faceI][qpI][stencilI]` stores the corresponding $\boldsymbol{C}_N^{f,q}$;
- the shared `stencilI` index associates each global cell with its coefficient.

The owner-side and neighbour-side labels are needed only while constructing
the operator. After the contributions have been merged, `fGrad()` only needs
the global cell ID and its final coefficient.

#### Why fixed-value boundary faces are different

A fixed-value boundary face has no neighbour cell whose reconstruction can be
averaged with the owner reconstruction. Instead, prescribed values at boundary
quadrature points are added as weighted observations to the owner-cell
least-squares reconstruction.

The resulting boundary-aware reconstruction is the only reconstruction used
at that face:

$$
\nabla u_f(\boldsymbol{x}_{f,q})
=
\nabla u_P^R(\boldsymbol{x}_{f,q}).
$$

There is no factor $1/2$ because the prescribed boundary value is an
observation used to determine the owner polynomial, not a second gradient to
be averaged with the owner gradient.

For every owner cell $P$ touching a selected fixed-value patch, the function
collects all prescribed quadrature points on all selected faces of that cell.
One address is stored as

```text
{meshFaceID, quadraturePointID}
```

Using all points belonging to $P$ means that a corner cell touching two
fixed-value faces uses the prescribed data from both faces in the same
reconstruction.

At this stage, only the quadrature-point locations and their future addresses
are known. `calcFaceGradCoeffs()` does not read the actual prescribed field
values. `fGrad()` obtains those values later through
`patchFaceQuadValues()`.

---

### Boundary equations

Let $c$ denote one prescribed boundary quadrature point with location
$\boldsymbol{x}_c$ and prescribed value $u_c^D$. Define
$$
b_c
=
u_c^D-\overline{u}_P.
$$

The boundary observation equation is

$$
b_c
=
\sum_{\boldsymbol{\alpha}}
D_{c\boldsymbol{\alpha}}
s_{\boldsymbol{\alpha}}^P,
$$

where $\boldsymbol{s}_P$ is the scaled-derivative vector introduced in the
`calcCellCoeffs()` section. The boundary matrix entry is

$$
D_{c\boldsymbol{\alpha}}
=
\frac{
\boldsymbol{r}_c^{\boldsymbol{\alpha}}
-M_P^{\boldsymbol{\alpha}}
}{
h_P^{|\boldsymbol{\alpha}|}\boldsymbol{\alpha}!
},
\qquad
\boldsymbol{r}_c
=
\boldsymbol{x}_c-\boldsymbol{x}_P.
$$

This is the mean-free polynomial evaluated at a point. It differs from a
neighbour-cell equation because no volume average over the boundary point is
required.

For a quadratic reconstruction in 2D, define

$$
r_{c,x}=x_c-x_P,
\qquad
r_{c,y}=y_c-y_P.
$$

The boundary observation equation written with ordinary, unscaled derivatives
is

$$
\begin{aligned}
u_c^D-\overline{u}_P
={}&
r_{c,x}u_{x,P}
+r_{c,y}u_{y,P}
\\
&+\frac{1}{2}
\left(r_{c,x}^2-M_{xx}^P\right)u_{xx,P}
\\
&+\left(r_{c,x}r_{c,y}-M_{xy}^P\right)u_{xy,P}
\\
&+\frac{1}{2}
\left(r_{c,y}^2-M_{yy}^P\right)u_{yy,P}.
\end{aligned}
$$

The left-hand side is known: it is the prescribed value at boundary point
$c$ minus the stored average of owner cell $P$. The quantities multiplying the
unknown derivatives form one boundary observation.

The code solves for the scaled derivative vector

$$
\boldsymbol{s}_P
=
\begin{bmatrix}
h_Pu_{x,P}
&
h_Pu_{y,P}
&
h_P^2u_{xx,P}
&
h_P^2u_{xy,P}
&
h_P^2u_{yy,P}
\end{bmatrix}^{T}.
$$

Therefore, the corresponding row of the boundary matrix is

$$
\boldsymbol{D}_{c,P}
=
\begin{bmatrix}
r_{c,x}/h_P
&
r_{c,y}/h_P
&
\left(r_{c,x}^2-M_{xx}^P\right)/(2h_P^2)
&
\left(r_{c,x}r_{c,y}-M_{xy}^P\right)/h_P^2
&
\left(r_{c,y}^2-M_{yy}^P\right)/(2h_P^2)
\end{bmatrix}.
$$

The same observation equation can consequently be written as

$$
u_c^D-\overline{u}_P
=
\boldsymbol{D}_{c,P}\boldsymbol{s}_P.
$$

The factors $1/2$ appear only for the pure second derivatives $u_{xx}$ and
$u_{yy}$ because $2!=2$. The mixed derivative $u_{xy}$ has no factor $1/2$
because its factorial denominator is $1!1!=1$.

In 3D, the same row additionally contains the $z$, $xz$, $yz$, and $zz$
terms.

Every collected prescribed quadrature point supplies one row of
$\boldsymbol{D}_P$ and one entry of the known-data vector
$\boldsymbol{b}_P$.

#### Combined cell and boundary least-squares system

The cell equations used in `calcCellCoeffs()` are rebuilt for the boundary
cell. The prescribed boundary-point equations are then appended to them. The
complete known-data vector is
$$
\boldsymbol{y}_P
=
\begin{bmatrix}
\Delta\overline{\boldsymbol{u}}_P
\\
\boldsymbol{b}_P
\end{bmatrix},
$$

where $\Delta\overline{\boldsymbol{u}}_P$ contains the stencil-cell average
differences and $\boldsymbol{b}_P$ contains the prescribed boundary-value
differences. The combined polynomial matrix and weight matrix are
$$
\boldsymbol{Q}_P^{\mathrm{aug}}
=
\begin{bmatrix}
\boldsymbol{Q}_P
&
\boldsymbol{D}_P^T
\end{bmatrix},
$$

and

$$
\boldsymbol{W}_P^{\mathrm{aug}}
=
\begin{bmatrix}
\boldsymbol{W}_P & \boldsymbol{0}
\\
\boldsymbol{0} & \boldsymbol{W}_P^D
\end{bmatrix}.
$$

The first $N_n$ columns of $\boldsymbol{Q}_P^{\mathrm{aug}}$ represent
stencil-cell average equations. The following $N_c$ columns represent
prescribed boundary-point equations. This column-oriented storage is the
orientation expected by `QRSolve()`.

The boundary-aware reconstruction solves

$$
\begin{aligned}
\min_{\boldsymbol{s}_P}\quad
&\left\|
\left(\boldsymbol{W}_P^{\mathrm{aug}}\right)^{1/2}
\left(
\left(\boldsymbol{Q}_P^{\mathrm{aug}}\right)^T
\boldsymbol{s}_P
-\boldsymbol{y}_P
\right)
\right\|_2^2
\end{aligned}
$$

Thus, prescribed boundary values are weighted least-squares observations.
They are not enforced as exact equality constraints, and no arbitrarily large
penalty is used. The same distance-based `weightFunction` is applied to both
cell and boundary observations.

`QRSolve()` directly calculates the combined coefficient matrix
$\boldsymbol{A}_P^{\mathrm{aug}}$. The scaled derivatives are therefore

$$
\boldsymbol{s}_P
=
\boldsymbol{A}_P^{\mathrm{aug}}
\boldsymbol{y}_P.
$$

The implementation uses a weighted QR decomposition. Splitting
the columns according to their input data gives
$$
\boldsymbol{A}_P^{\mathrm{aug}}
=
\begin{bmatrix}
\boldsymbol{A}_P^{\mathrm{cell}}
&
\boldsymbol{A}_P^D
\end{bmatrix}.
$$

The first part maps the stencil-cell differences, while the second part maps
the prescribed boundary-value differences.

#### Gradient evaluation at a boundary quadrature point

For every target quadrature point $\boldsymbol{x}_{f,q}$, the matrix
$\boldsymbol{L}_{f,q}$ differentiates the scaled mean-free polynomial:

$$
\nabla u_P^R(\boldsymbol{x}_{f,q})
=
\boldsymbol{L}_{f,q}\boldsymbol{s}_P.
$$

The complete gradient map is

$$
\boldsymbol{C}^{f,q}
=
\boldsymbol{L}_{f,q}\boldsymbol{A}_P^{\mathrm{aug}}.
$$

This matrix is named `gradientMap` in the code. Its first $N_n$ columns
multiply stencil-cell differences, and its following $N_c$ columns multiply
prescribed boundary-value differences.

Before conversion to absolute values, the gradient is

$$
\begin{aligned}
\nabla u_P^R(\boldsymbol{x}_{f,q})
={}&
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_{N}^{f,q}
\left(
\overline{u}_N-\overline{u}_P
\right)
\\
&+
\sum_{c=1}^{N_c}
\boldsymbol{C}_{c,D}^{f,q}
\left(
u_c^D-\overline{u}_P
\right).
\end{aligned}
$$

Expanding the differences gives

$$
\begin{aligned}
\nabla u_P^R(\boldsymbol{x}_{f,q})
={}&
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_{N}^{f,q}\overline{u}_N
+
\sum_{c=1}^{N_c}
\boldsymbol{C}_{c,D}^{f,q}u_c^D
\\
&-
\left[
\sum_{N\in\mathcal{S}_P}
\boldsymbol{C}_{N}^{f,q}
+
\sum_{c=1}^{N_c}
\boldsymbol{C}_{c,D}^{f,q}
\right]
\overline{u}_P.
\end{aligned}
$$

This is why the implementation accumulates both the cell and prescribed-data
coefficients in `coefficientSum` and then executes

```cpp
coefficients[ownIndex] -= coefficientSum;
```

The resulting absolute coefficients again give zero gradient for a constant
field when the cell averages and prescribed values contain the same constant.

#### Boundary-face storage and use by fGrad

For a selected fixed-value boundary face:

- `faceGradStencil()[faceI]` stores the owner cell and its cell stencil;
- `faceGradCoeffs()[faceI][qpI]` stores coefficients multiplying unknown cell
  averages;
- `faceBoundaryDataStencil()[faceI]` stores the
  `{meshFaceID, quadraturePointID}` address of every prescribed value used by
  the owner-cell reconstruction;
- `faceBoundaryDataCoeffs()[faceI][qpI]` stores the vector coefficient
  multiplying every prescribed value.

During `fGrad()` evaluation, the unknown-cell contribution is

$$
\sum_{N\in\mathcal{S}_f}
\boldsymbol{C}_N^{f,q}\overline{u}_N,
$$

and the known-data contribution is

$$
\sum_{c=1}^{N_c}
\boldsymbol{C}_{c,D}^{f,q}u_c^D.
$$

For a vector field, multiplication of each vector coefficient by the
corresponding vector value produces the reconstructed gradient tensor.

Boundary patches not selected by `includePatchInStencils_` are not processed
by this boundary augmentation. Processor-face gradients are handled
separately in `fGrad()` by exchanging the two one-sided reconstructed
gradients and averaging them.
