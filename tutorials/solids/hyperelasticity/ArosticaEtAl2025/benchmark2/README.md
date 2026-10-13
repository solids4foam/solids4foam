---
sort: 12
---

# Cardiac elastodynamics benchmark, idealised biventricle: `ArosticaEtAl2025/benchmark2`

---

Prepared by Chanikya Valeti (UCD), with Claude Opus 5.5

---

## Tutorial Aims

- Reproduce benchmark 2 (biventricular mechanics) of Aróstica et al. [1].
- Use the same `viscoActiveHolzapfelOgden` law and `solidSpringDashpot`
  condition as `benchmark1`, with two cavity pressures.

## Case Overview

**Geometry and mesh** (Section 4.2, Table 6). The idealised biventricle is the
tetrahedral mesh distributed with the benchmark (cardiac-benchmark-toolkit
[2], `data/biv_mesh_and_fibers_{coarse,fine}_p1`). `Allrun` downloads it, at
a fixed commit of the toolkit and checked against SHA-256 sums, into `data/`,
and `../scripts/convertBiVMesh.py` converts it to an OpenFOAM `polyMesh`. This
needs `python3` with `numpy` and `h5py`, and internet access (or the two
`.h5` files placed in `data/<mesh>` beforehand). Patches: `base` (marker 10),
`endo_rv` (20), `endo_lv` (30), `epi` (40).

| Mesh | Tetrahedra | Nodes |
| --- | --- | --- |
| `coarse` (default) | 45 304 | 11 444 |
| `fine` | 121 133 | 27 807 |

**Fibres** (Section 4.3). The BT-LDRB fibres (∓60°) of the benchmark, computed
by the participants with life^x and distributed as nodal (P1) vectors. The
converter averages them to cell and face centres and re-orthonormalises
(f0, then s0 by Gram–Schmidt, n0 = f0 × s0). The toolkit's `sheet` vector is
used as s0, as in the reference implementation.

**Material.** As in benchmark 1 (Tables 1–2), with the active stress of
Table 8 (σ0 = 1.5·10⁵ Pa, t_sys = 0.163 s, t_dias = 0.5 s; τ_max ≈ 120.1 kPa).

**Pressures** (Table 9): p_lv(t) on `endo_lv` (max 16 490 Pa) and p_rv(t) on
`endo_rv` (max 4 163 Pa), both follower pressures, tabulated by `Allrun` in
`constant/pressureLV` and `constant/pressureRV` with
`../scripts/makeTimeSeries.py --pressure biv_lv|biv_rv`.

**Boundary conditions** (Table 7): `base` spring-dashpot α = 10⁶ Pa/m,
β = 5·10³ Pa s/m; `epi` normal spring-dashpot α = 10⁸ Pa/m, β = 5·10³ Pa s/m.

**Output**: displacement at p0 = (0.025, 0.03, 0), p1 = (0, 0.03, 0) and
p2 = (0.025, 0, 0.072) m.

**Results so far** (coarse mesh, pressure formulation, serial, full beat
t ∈ [0, 1] s; every step converges in 3 Newton iterations, about 4 h; run
with this PR's first spring-dashpot condition, before the switch to the
`solidSpringDashpot` of #528, so they are to be re-run). The
timing of contraction and relaxation follows the participants' at all three
points, and the deformed shape at t = 0.3 and 0.5 s matches Fig. 11 of [1].
Of the nine displacement components, 6 are inside the participants' range at
t = 0.17 s (end of filling), 6 at 0.3 s, 5 at 0.5 s, 7 at 0.7 s and 6 at
1 s; all of p2 is inside during systole. Outside it:

- the long-axis shortening in systole is too large: p0 x 9.6 vs a median of
  8.2 mm (max 8.9) at t = 0.5 s, and p1 x just above the range, as in the
  step 1 of `benchmark1`;
- at low pressure the response is too stiff: during filling some components
  are about a third of the participants' (p0 x 0.08 vs 0.45 mm, p0 z 0.41 vs
  1.03 mm at t = 0.17 s), and in late relaxation the displacements return to
  zero by t ≈ 0.75 s, while the participants keep 0.2-0.5 mm at t = 1 s;
- p1 y overshoots to +0.5 mm at the onset of contraction and reaches
  −1.45 mm at t = 0.5 s (range −1.22 to −0.71 mm).

A polyhedral mesh made from the coarse tetrahedra by `polyDualMesh` (12 573
cells, fibres mapped from the nearest tetrahedron) was not better: it moved
p0 and p1 x into the range, but over-contracted p2, and its filling phase was
as stiff. The stiffness at low load is therefore attributed to the momentum
stabilisation (`diffStencilLaplacian`, `scaleFactor 100`) rather than to the
tetrahedral cells; a study of the stabilisation scale and the `fine` mesh is
still to be run.

The displacement formulation, with the benchmark's κ = 10⁶ Pa in the momentum
equation, follows the pressure formulation during filling (within 0.14 mm),
but lags it at the onset of contraction, by up to 1.2 mm at t = 0.22 s: it
is stiffer, as expected from volumetric locking and from its stabilisation
being scaled by κ. Lowering κ would change the benchmark's material, so the
pressure formulation is the default.

## Running the Case

```bash
./Allrun                  # coarse mesh, mixed (pressure) formulation
./Allrun coarse displacement
./Allrun fine             # pressure formulation: serial only
./Allrun coarse displacement parallel
python3 ../scripts/compareWithReference.py coarse   # also run by Allrun
./Allclean
```

The reference curves in `verification/reference/{coarse,fine}` are the
participants' submitted solutions from the toolkit [2] (MIT licence).
`compareWithReference.py` plots all of them, and takes the median over one
solution per group.

## References

[1] R. Aróstica et al., "A software benchmark for cardiac elastodynamics",
CMAME 435 (2025) 117485. <https://doi.org/10.1016/j.cma.2024.117485>

[2] <https://github.com/Reidmen/cardiac_benchmark_toolkit>
