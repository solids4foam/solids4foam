---
sort: 12
---

# Cardiac elastodynamics benchmark, idealised biventricle: `ArosticaEtAl2025/benchmark2`

## Tutorial Aims

- Reproduce benchmark 2 (biventricular mechanics) of Aróstica et al. [1].
- Use the same `viscoActiveHolzapfelOgden` law and `solidSpringDashpot`
  condition as `benchmark1`, with two cavity pressures.

## Case Overview

**Geometry and mesh** (Section 4.2, Table 6). The idealised biventricle is the
tetrahedral mesh distributed with the benchmark (cardiac-benchmark-toolkit
[2], `data/biv_mesh_and_fibers_{coarse,fine}_p1`), converted to an OpenFOAM
`polyMesh` by `scripts/convertBiVMesh.py`. Patches: `base` (marker 10),
`endo_rv` (20), `endo_lv` (30), `epi` (40).

| Mesh | Tetrahedra | Nodes |
| --- | --- | --- |
| `coarse` (default, shipped) | 45 304 | 11 444 |
| `fine` (converted on request, see `Allrun`) | 121 133 | 27 807 |

**Fibres** (Section 4.3). The BT-LDRB fibres (∓60°) of the benchmark, computed
by the participants with life^x and distributed as nodal (P1) vectors. The
converter averages them to cell and face centres and re-orthonormalises
(f0, then s0 by Gram–Schmidt, n0 = f0 × s0). The toolkit's `sheet` vector is
used as s0, as in the reference implementation.

**Material.** As in benchmark 1 (Tables 1–2), with the active stress of
Table 8 (σ0 = 1.5·10⁵ Pa, t_sys = 0.163 s, t_dias = 0.5 s; τ_max ≈ 120.1 kPa).

**Pressures** (Table 9): p_lv(t) on `endo_lv` (max 16 490 Pa) and p_rv(t) on
`endo_rv` (max 4 163 Pa), both follower pressures, tabulated in
`constant/timeSeries` by `scripts/makeTimeSeries.py --pressure biv_lv|biv_rv`.

**Boundary conditions** (Table 7): `base` spring-dashpot α = 10⁶ Pa/m,
β = 5·10³ Pa s/m; `epi` normal spring-dashpot α = 10⁸ Pa/m, β = 5·10³ Pa s/m.

**Output**: displacement at p0 = (0.025, 0.03, 0), p1 = (0, 0.03, 0) and
p2 = (0.025, 0, 0.072) m.

**Results so far** (coarse mesh, pressure formulation, serial, run to t = 0.51
s; every step converges in 3 Newton iterations): the deformed shape at t = 0.3
and 0.5 s matches the behaviour in Fig. 11 of [1], and 11 of the 18 values of
the three points at those times are inside the participants' range; all of p2
is. Outside it are the long-axis x of p0 and p1 (by 0.1 to 0.7 mm; p0 x is 9.6
vs a median of 8.2 mm at t = 0.5 s), p0 y (by 0.04 to 0.06 mm) and p1 y at t =
0.5 s (−1.45 vs −1.22 to −0.71 mm). The filling phase (t < 0.17 s) is too stiff:
some components are about a third of the participants' (p0 uz 0.4 vs about 1
mm), which the hexahedral benchmark 1 does not show, so the tetrahedral cells
are the main suspect.

## Running the Case

```bash
./Allrun                  # coarse mesh, mixed (pressure) formulation
./Allrun coarse displacement
./Allrun fine                    # pressure formulation: serial only
./Allrun coarse displacement parallel
python3 scripts/compareWithReference.py coarse   # also run by Allrun
python3 scripts/checkJacobian.py                 # J = det(F) over time
```

The reference curves in `verification/reference/{coarse,fine}` are the
participants' submitted solutions from the toolkit [2] (MIT licence).

```note
Tetrahedra are the least favourable cells for cell-centred finite volumes.
If the coarse tet mesh proves too stiff or too noisy, `polyDualMesh` turns it
into a polyhedral mesh with the same boundary (the fibres must then be
re-sampled at the new cell centres).
```

## References

[1] R. Aróstica et al., "A software benchmark for cardiac elastodynamics",
CMAME 435 (2025) 117485. <https://doi.org/10.1016/j.cma.2024.117485>

[2] <https://github.com/Reidmen/cardiac_benchmark_toolkit>
