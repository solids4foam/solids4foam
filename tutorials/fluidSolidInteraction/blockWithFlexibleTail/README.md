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

The fluid mesh has 12,920 cells and the tail mesh has 160 cells. Both have 80
faces along each long tail surface and two faces at the tip, so the interface
is conformal. This is a shake-down mesh, not a mesh-converged discretisation.

---

## Running the Case

Run the partitioned IQN-ILS case, with an implicit-segregated solid solve,
using:

```bash
./Allrun
```

Run the quasi-monolithic Newton variant with:

```bash
./Allrun monolithic
```

Append `parallel` to either command for a four-rank run. A parallel
quasi-monolithic run uses `decomposeParMonolithic`, with the coupled regions
and conformal tail patches declared in `system/decomposeParDict`.

The default time step is 0.001 s, matching the second-mode study in [2, 3],
and the case runs to 10 s. The `solidPointDisplacement` function object writes
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

The partitioned case completed a three-step smoke test with the stored mesh and
settings. The three steps needed 6, 6 and 5 IQN-ILS iterations and 4.5 s in
serial, but this is far too short to measure an oscillation. A long shake-down
could not be submitted because the test sandbox prevented Slurm socket access;
PETSc's OpenMPI initialisation was blocked for the same reason, so the
quasi-monolithic time stepping was not exercised. Its case selection, 162
conformal interface-face pairs and four-rank combined decomposition were
checked with `decomposeParMonolithic`. No solids4foam amplitude, frequency or
quasi-monolithic Newton-iteration count is therefore claimed here.

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
