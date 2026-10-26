---
sort: 3
---

# Block with a flexible tail: `blockWithFlexibleTail`

---

Prepared by Philip Cardiff and the solids4foam contributors

---

## Tutorial Aims

- Demonstrates the Wall--Ramm fluid-solid interaction benchmark in which
  vortex shedding from a rigid square block excites a flexible cantilever.
- Provides partitioned IQN-ILS and quasi-monolithic Newton coupling variants.
- Records the benchmark definition and separates the two parameter sets that
  are commonly conflated in later publications.

---

## Case Overview

A 1 cm square rigid block is centred in a two-dimensional channel. Its
upstream face is 4.5 cm from the inlet. A 4 cm long and 0.06 cm thick elastic
tail is clamped to the centre of its downstream face. The channel is 19.5 cm
long and 12 cm high, leaving 10 cm between the tail tip and the outlet. The
top and bottom are slip boundaries, the block and tail are no-slip, and the
outlet pressure is zero. These dimensions are taken from Figure 6.6 of
Cardiff's accessible thesis [2]. Wall's original thesis is available [1], but
the scan did not provide searchable benchmark values and its German figure
could not be checked independently in this work.

The fluid starts from rest and has a uniform 31.5 cm/s inlet speed. A transverse
inlet velocity of only 0.01 cm/s is applied for the first 0.1 s to break the
otherwise exact symmetry; the benchmark inlet is used thereafter. The solid
uses a St. Venant--Kirchhoff law and the plane-stress assumption. Table 1 gives
the values checked against Sections 6.2--6.3 and Figure 6.6 of [2]. The same
values and Reynolds number are reported in the associated paper [3].

**Table 1: Physical and geometric parameters used by this tutorial.**

| Parameter | Value in source units | SI value | Source |
| :-- | --: | --: | :-- |
| Channel length | 19.5 cm | 0.195 m | Figure 6.6 of [2] |
| Channel height | 12 cm | 0.12 m | Figure 6.6 of [2] |
| Inlet to block upstream face | 4.5 cm | 0.045 m | Figure 6.6 of [2] |
| Square block side | 1 cm | 0.01 m | Figure 6.6 of [2] |
| Tail length | 4 cm | 0.04 m | Figure 6.6 of [2] |
| Tail thickness | 0.06 cm | 0.0006 m | [2, 3] |
| Fluid density, $$\rho_F$$ | 1.18e-3 g/cm$$^3$$ | 1.18 kg/m$$^3$$ | [2, 3] |
| Dynamic viscosity, $$\mu_F$$ | 1.82e-4 g/(cm s) | 1.82e-5 Pa s | [2, 3] |
| Inlet speed, $$U$$ | 31.5 cm/s | 0.315 m/s | [2, 3] |
| Solid density, $$\rho_S$$ | 2.0 g/cm$$^3$$ | 2000 kg/m$$^3$$ | [2, 3] |
| Young's modulus, $$E_S$$ | 2.0e6 g/(cm s$$^2$$) | 0.2 MPa | [2, 3] |
| Poisson's ratio, $$\nu_S$$ | 0.35 | 0.35 | [2, 3] |
| Reynolds number, $$\rho_F U H/\mu_F$$ | 204 | 204 | [2, 3] |

There are at least two widely cited Wall--Ramm tail cases. This tutorial uses
the dense, stiffer tail and 31.5 cm/s inflow above, which excites the second
bending mode. The more frequently tabulated case instead uses 51.3 cm/s,
$$\rho_S = 0.1$$ g/cm$$^3$$ and $$E_S = 2.5\times10^6$$
g/(cm s$$^2$$); its roughly 1.0--1.34 cm amplitudes must not be used as
references for this tutorial. The distinction is explicit in [2] and is also
discussed by Sánchez et al. [5].

Jaiman and Joshi's 2022 textbook [6] was bibliographically verified, but the
relevant pages were not openly accessible, so no numerical value here is
attributed to it. The suggested Wood et al., Habchi et al. and Olivier et al.
citations were not used: accessible summaries concern the lighter 51.3 cm/s
case or did not establish that they used this exact parameter set.

The fluid mesh has 13,020 cells and the tail mesh has 320 cells. Both have 80
faces along each long tail surface and four faces at the tip, so the interface
is conformal. This is a shake-down mesh, not a mesh-converged discretisation.

---

## Running the Case

Run the partitioned IQN-ILS case, with an implicit-segregated solid solve,
using:

```bash
./Allrun
```

Use the cubic high-order finite-volume solid discretisation with:

```bash
./Allrun highOrder
```

This variant requires PETSc and is unavailable on foam-extend. The
implementation supports polynomial orders 1, 2 and 3; this tutorial selects
the cubic option. The high-order residual and Jacobian are not implemented for
the velocity-based solid model used by the quasi-monolithic coupling, so
`highOrder` cannot be combined with `monolithic`.

A controlled solid-only gravity-step study used PETSc SNES, `alpha`
stabilisation with scale 0.1, and the high-order residual and Jacobian for
every polynomial order. The Euler--Bernoulli references are 0.607 Hz and
3.804 Hz for the first two bending modes.

| Tail mesh | Order | First mode (Hz) | Error | Second mode (Hz) | Error |
| :-- | --: | --: | --: | --: | --: |
| 80 x 2 | 1 | 0.4740 | -21.9% | 3.5035 | -7.9% |
| 80 x 2 | 2 | 1.0388 | +71.1% | 6.2595 | +64.6% |
| 80 x 2 | 3 | 0.7918 | +30.4% | 4.9145 | +29.2% |
| 80 x 4 | 1 | 0.5907 | -2.7% | 3.7180 | -2.3% |
| 80 x 4 | 2 | 0.7900 | +30.1% | 4.8875 | +28.5% |
| 80 x 4 | 3 | 0.6068 | -0.0% | 3.7795 | -0.6% |

The 80 x 4 cubic result motivates the mesh and order selected here. The
unexpectedly large errors for order 2 are not yet explained. Earlier
segregated order-1 calculations gave errors of +39.2%/+37.5% on 80 x 2 and
+19.4%/+18.2% on 80 x 4. Their disagreement with the controlled order-1 rows
shows that stabilisation or the solution path, rather than the mesh alone,
strongly affects the bending response.

Run the quasi-monolithic Newton variant with:

```bash
./Allrun monolithic
```

The quasi-monolithic solid uses the same 80 x 4 linear mesh and nominal
`alpha` stabilisation scale 0.1 as the controlled order-1 path. It remains a
distinct velocity-based formulation without the partitioned model's
high-order residual/Jacobian machinery. Its solid-only frequencies are 0.7245
and 4.4990 Hz, 19.4% and 18.3% above the beam references. The compact
velocity residual therefore remains too stiff on this mesh; a stabilisation
and mesh-refinement study is required before its FSI result can be compared
with the partitioned calculation.

Append `parallel` to either command for a four-rank run. A parallel
quasi-monolithic run uses `decomposeParMonolithic`, with the coupled regions
and conformal tail patches declared in `system/decomposeParDict`.

The default time step is 0.001 s, matching the second-mode study in [2, 3],
and the case runs to 10 s. Restart state is written every 1 s and the four most
recent checkpoints are retained. The `solidPointDisplacement` function object writes
the displacement of the initial tail-tip point `(0.095 0.06 0.0005)`. If
Gnuplot is installed, `Allrun` writes `deflection.pdf`.

---

## Expected Results

Cardiff [2] and Cardiff et al. [3] report the mesh and time-step study in
Table 2. Their amplitude is half the peak-to-peak vertical tail-tip motion.
Hübner et al. and Xia et al. values in the final two rows are reproduced from
that accessible comparison rather than independently extracted from the
original papers.

**Table 2: Published second-mode tail-tip response for the parameter set in
Table 1.**

| Fluid mesh | Solid mesh | Time step (s) | Amplitude (cm) | Frequency (Hz) |
| --: | :--: | --: | --: | --: |
| 6,000 | 40 x 1 | 0.001 | 0.78 | 2.98 |
| 6,000 | 80 x 2 | 0.001 | 0.78 | 2.98 |
| 25,000 | 120 x 1 | 0.001 | 0.76 | 2.97 |
| 50,000 | 160 x 1 | 0.001 | 0.76 | 2.97 |
| 6,000 | 40 x 1 | 0.0005 | 0.78 | 2.97 |
| Hübner et al. [4], as reported by [2, 3] | - | - | 0.80 | 3.10 |
| Xia et al., as reported by [2, 3] | - | - | 0.81 | 3.30 |

Table 3 summarises the current four-rank shake-down results. The controlled
order-1 run uses the same PETSc SNES, high-order residual/Jacobian machinery
and `alpha` stabilisation as the modal study, with `polynomialOrder 1`.

**Table 3: Current computed tail-tip response. Amplitudes are half the
peak-to-peak vertical displacement.**

| Solid setup | Time window (s) | Amplitude (cm) | Frequency (Hz) | Status |
| :-- | :--: | --: | --: | :-- |
| 80 x 4, controlled order 1 | 8--10 | 0.650 | 3.24 | Stopped at 11.135 s |
| 80 x 4, order 3 | 8--10 | 0.651 | 3.25 | Stopped at 10.316 s |
| 80 x 4, revised order 1 | 8--10 | 0.650 | 3.18 | Stopped at 10.655 s |
| 80 x 4, revised order 3 | 8--9.353 | 0.644 | 3.30 | Stopped at 9.353 s |
| Literature [2--4] | Limit cycle | 0.76--0.81 | 2.97--3.3 | Reference |

Both orders now reproduce the published frequency and give essentially the
same 8--10 s amplitude. Growth is slowing, but 10 s is not a demonstrated
limit cycle and neither attempted 25 s run completed. The order-1 run reached
0.689 cm before IQN-ILS stalled at 11.135 s. The cubic run reached a similar
amplitude before a large IQN-ILS correction and solid linear-solver failure at
10.316 s.

Mesh checks identify severe deformation of the fluid grid near the tail tip.
At 10.3 s in the cubic run the maximum non-orthogonality was 91.5 degrees, the
maximum skewness was 4.95, the minimum cell volume had fallen to 3.7% of its
initial value, and one face was incorrectly oriented. Changing to unsquared
inverse-distance diffusion, filtering IQN-ILS modes, guarding mode reuse and
rebuilding the solid preconditioner did not delay the failure in fresh runs;
those changes have therefore been removed.

The `solids4foamSolidModel` motion solver with a total-Lagrangian
`diffusionHyperElastic` law is a promising replacement and is already used by
the quasi-monolithic variant. It cannot currently be selected by the
partitioned interface, however: partitioned coupling prescribes
`pointMotionU`, whereas this solver owns `pointD`, and the interface aborts at
its first mesh update. The tutorial retains `velocityLaplacian` until that
solver-interface gap is fixed. Restart checkpoints remain enabled so a later
large-amplitude run can continue without repeating its first 10 s.

Diagnostics point to excessive numerical bending stiffness in the original
shake-down solid discretisation. A solid-only gravity-step test gave a first
natural frequency of 1.031 Hz on the original 80 x 2 segregated mesh, versus
0.607 Hz from beam theory and the finite-element reference. In contrast,
probes in the computed wake give approximately 4.0 Hz, close to the published
uncoupled-fluid value of 3.71 Hz. The controlled modal study above supports
the current 80 x 4 mesh with the cubic discretisation. Its repeated FSI run
and a mesh-convergence study are still required.

The controlled order-1 10 s run took 3,574 s on four ranks. The earlier 10 s
partitioned and quasi-monolithic runs took 2,140 s and 8,336 s respectively;
their approximately 1.2 Hz first-mode response used the original, excessively
stiff solid setup and is superseded by Table 3. That quasi-monolithic run used
an 80 x 2 tail and `diffStencilLaplacian` stabilisation at scale 0.5; its
80 x 4, `alpha`-stabilised replacement is under test, but its solid-only modal
error is approximately 19%. The old run averaged 1.77 solid and 2.00 fluid
SNES Newton iterations per time step. These timings are indicative only.

The reusable analysis script reports half peak-to-peak amplitude, peak
magnitude, detrended upward-zero-crossing frequency and Hann-windowed FFT
frequency over a selected interval. For example:

```bash
python3 verification/analyseTipDisplacement.py \
    --case order1 postProcessing 8 10
```

---

## References

[1] W. A. Wall, [*Fluid-Struktur-Interaktion mit stabilisierten Finiten
Elementen*](https://edocs.tib.eu/files/e01dd01/367181061l.pdf), PhD thesis,
Universität Stuttgart, 1999.

[2] P. N. Cardiff, [*Development of parallel strongly coupled hybrid
fluid-structure interaction technology involving thin geometrically
non-linear structures*](https://repository.up.ac.za/bitstreams/250ee98a-041b-4185-8b54-6d9a320acb6f/download),
PhD thesis, University of Cape Town, 2012, Sections 6.2--6.3.

[3] P. Cardiff, A. Karac and A. Ivanković, ["A matrix free, partitioned
solution of fluid-structure interaction problems using finite volume and
finite element methods"](https://doi.org/10.1016/j.euromechflu.2014.10.002),
*European Journal of Mechanics B/Fluids*, 49 (2015), 272--286.

[4] B. Hübner, E. Walhorn and D. Dinkler, ["A monolithic approach to
fluid-structure interaction using space-time finite
elements"](https://doi.org/10.1016/j.cma.2004.01.024), *Computer Methods in
Applied Mechanics and Engineering*, 193 (2004), 2087--2104.

[5] R. Sánchez, J. Thomas, R. Palacios and S. Terrapon,
[*Fluid-Structure Interaction Problems using SU2 and External Finite-Element
Solvers*](https://economon.github.io/docs/AIAA-2016-0205.pdf), AIAA 2016-0205,
2016.

[6] R. K. Jaiman and V. Joshi, [*Computational Mechanics of Fluid-Structure
Interaction*](https://doi.org/10.1007/978-981-16-5355-1), Springer, 2022.

---
