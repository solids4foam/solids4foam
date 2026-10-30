---
sort: 12
---

# Cardiac elastodynamics benchmark, idealised biventricle: `ArosticaEtAl2025/benchmark2`

---

Prepared by Chanikya Valeti and Philip Cardiff (UCD), with Claude Opus 5.5

---

## Tutorial Aims

- Reproduce benchmark 2 (biventricular mechanics) of Aróstica et al. [1].
- Use the same `viscoActiveHolzapfelOgden` law and `solidSpringDashpot`
  condition as `benchmark1`, with two cavity pressures.

## Case Overview

**Geometry and mesh** (Section 4.2, Table 6). The idealised biventricle is the
tetrahedral mesh distributed with the benchmark (cardiac-benchmark-toolkit
[2], `data/biv_mesh_and_fibers_{coarse,fine}_p1`), converted once to OpenFOAM
format. `Allrun` downloads the converted meshes and fibres from the
solids4foam website (`ArosticaEtAl2025.zip`, checked against its SHA-256 sum)
into `data/`, with `curl` or `wget`, and `unzip`. Without internet access,
place `ArosticaEtAl2025.zip` in `data/` first. The
archive also holds the converter, `convertBiVMesh.py`, and the provenance of
the toolkit's files. Patches: `base`, `endo_rv`, `endo_lv` and `epi`.

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
`constant/pressureLV` and `constant/pressureRV` by the
`makeArosticaTimeSeries` utility shared with `benchmark1`
(`-pressure bivLV|bivRV`).

**Boundary conditions** (Table 7): `base` spring-dashpot α = 10⁶ Pa/m,
β = 5·10³ Pa s/m; `epi` normal spring-dashpot α = 10⁸ Pa/m, β = 5·10³ Pa s/m.

**Output**: displacement at p0 = (0.025, 0.03, 0), p1 = (0, 0.03, 0) and
p2 = (0.025, 0, 0.072) m.

**Formulation and solver.** As in `benchmark1`: the mixed
displacement-pressure formulation (default, serial only) or the displacement
formulation, with PETSc SNES, BDF2 and Δt = 1 ms over t ∈ [0, 1] s.

## Expected Results

The participants' solutions in [1] spread over 1-6 mm, depending on the
component, over the beat. At t = 0.5 s (end of contraction), their median
displacements are about 8 mm (p0 x), 13 mm (p1 x) and 15 mm (p2 x).
If `gnuplot` is available, `Allrun` plots the displacements of p0, p1 and p2
against all the participants' solutions in `verification/<mesh>.pdf`.

```note
The results depend on the momentum stabilisation (`diffStencilLaplacian`,
`scaleFactor` in `constant/solidProperties.pressure`): with too large a
scale the wall is too stiff at low load, which shows in the filling phase
(t < 0.17 s) and in late relaxation. The displacement formulation, with the
benchmark's κ = 10⁶ Pa in the momentum equation, is stiffer still at the onset
of contraction (volumetric locking), so the pressure formulation is the
default.
```

## Running the Case

```bash
./Allrun                  # coarse mesh, mixed (pressure) formulation
./Allrun coarse displacement
./Allrun fine             # pressure formulation: serial only
./Allrun coarse displacement parallel
./Allclean
```

The reference curves in `verification/reference/{coarse,fine}` are the
participants' submitted solutions from the toolkit [2] (MIT licence).

## References

[1] R. Aróstica et al., "A software benchmark for cardiac elastodynamics",
CMAME 435 (2025) 117485. <https://doi.org/10.1016/j.cma.2024.117485>

[2] <https://github.com/Reidmen/cardiac_benchmark_toolkit>
