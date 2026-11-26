# multiMaterialPlate verification study

This opt-in study has two parts:

- a **thickness study**, which runs the tutorial mesh for the four mid-layer
  thicknesses, `h/W = 0`, `0.002`, `0.005` and `0.01`, and plots the normal
  stress along the interfaces against the boundary element results of Ioka et
  al. (2007);
- a **mesh study**, which coarsens and refines the tutorial mesh for
  `h = 0.01W` and checks that the normal stress along the interfaces
  converges.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study
compares it with published results and checks mesh convergence. Nothing here
is run by `tutorials/Alltest` or `tutorials/Alltest-regression`.

## Running

Source an OpenFOAM environment with solids4foam built (with PETSc for the
`petscSnes` variant), and run:

```bash
cd tutorials/solids/multiMaterial/multiMaterialPlate/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick                 # two thicknesses and two levels only
./Allverify --studies thickness     # one of the two studies
./Allverify --variants segregated   # implicit segregated solution algorithm
./Allverify --levels 2,3,4          # a subset of the mesh family
./Allverify --thicknesses 0,0.01    # a subset of the thicknesses
./Allverify --reuse                 # resume a sweep without re-running cases
```

Each run is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. Results are written to the ignored
`verification/postProcessing/` directory as `thickness_study.csv`,
`mesh_convergence.csv`, `verification_summary.md`, the sampled profiles under
`profiles/`, and, when `gnuplot` is available, the plots
`multiMaterialPlate_thickness_interface13.pdf` and `..._interface23.pdf` of
the thickness study, and `multiMaterialPlate_mesh_interface13.pdf`,
`..._interface23.pdf` and `..._midLine.pdf` of the mesh study.

## Reference curves

The files `reference/ioka_sigmaxx_interface13_hOverW_<h/W>.dat` and
`reference/ioka_sigmaxx_interface23_hOverW_<h/W>.dat` hold the normal stress
along the two interfaces from Ioka et al. (2007), as `r/W` and
`sigma_xx/sigma_0`. They were digitised from the published figure with
WebPlotDigitizer, so they are drawn in the plots for a visual comparison only
and take no part in the acceptance criteria below. Without a mid-layer the
two are the same interface, the one between materials 1 and 2, so the curve
for `h/W = 0` is stored once, under interface 13, and drawn for both.
The `plot.gnuplot` script of the tutorial draws the curve for the thickness
that was run in `sigmaxx.png` as well.

The plot scripts choose the terminal from the extension of the output file.
The figures of the four thicknesses in the tutorial README are the plots of
the thickness study written as `.png`.

## Mesh levels

The tutorial ships level 3 of the family. Every level scales the block
divisions in the plane by a power of two and retains the tutorial grading,
with a ratio of 10000 between the largest and the smallest cell:

| Level | Cells along each outer material | Cells across the mid-layer | Cells |
|---:|---:|---:|---:|
| 1 | 30 x 60 | 125 | 11 100 |
| 2 | 60 x 120 | 250 | 44 400 |
| 3 | 120 x 240 | 500 | 177 600 |
| 4 | 240 x 480 | 1000 | 710 400 |

## Acceptance criteria

- **Equilibrium**: the mean of `sigma_xx/sigma_0` over each interface is one,
  whatever the materials. Its error must be below 1e-3 for every thickness on
  the tutorial mesh, and on the finest mesh of the mesh study.
- **Order of the singularity**: for the bi-material plate, `h = 0`, the order
  fitted to the stress in `1e-4 < r/W < 3e-4` must be within 5% of the
  analytical order of Bogy (1970). With a mid-layer, the
  fitted and analytical orders are reported but not judged, as the range in
  which the singular term dominates narrows with the mid-layer thickness.
- **Convergence**: `profile_rms_relative_change`, the RMS relative change of
  the stress along an interface between successive mesh levels in
  `1e-4 <= r/W <= 1`, must decrease monotonically on both interfaces.

A `--quick` run only exercises two thicknesses and the two coarsest meshes,
so it checks equilibrium and the order of the bi-material plate, but not
convergence.

`verification/Allclean` removes the generated `work/` and `postProcessing/`
directories, as does `./Allclean all` in the case directory; the plain case
`./Allclean` keeps them.

## Reference results

Recorded with OpenFOAM v2412 and the `petscSnes` variant on an Apple M4 Pro.
The whole study takes about one minute.

Thickness study, on the tutorial mesh:

| h/W | Interface | Cells | Mean stress error | Fitted order | Analytical order |
|---:|---:|---:|---:|---:|---:|
| 0 | 12 | 57 600 | 5.3e-06 | 0.0916 | 0.0934 |
| 0.002 | 13 | 117 600 | 2.9e-07 | 0.2661 | 0.2823 |
| 0.002 | 23 | 117 600 | 5.8e-06 | 0.0473 | 0.2392 |
| 0.005 | 13 | 117 600 | 4.6e-07 | 0.2894 | 0.2823 |
| 0.005 | 23 | 117 600 | 4.6e-06 | 0.1754 | 0.2392 |
| 0.01 | 13 | 177 600 | 2.3e-06 | 0.2876 | 0.2823 |
| 0.01 | 23 | 177 600 | 2.5e-06 | 0.2191 | 0.2392 |

The computed stress along both interfaces follows the digitised curves of Ioka
et al. for the four thicknesses.

Mesh study, for `h = 0.01W`:

| Level | Cells | Interface | Fitted order | Profile change | Time (s) |
|---:|---:|---:|---:|---:|---:|
| 1 | 11 100 | 13 | 0.3835 | – | <1 |
| 1 | 11 100 | 23 | 0.2985 | – | <1 |
| 2 | 44 400 | 13 | 0.3171 | 1.30e-02 | 1 |
| 2 | 44 400 | 23 | 0.2383 | 1.04e-02 | 1 |
| 3 | 177 600 | 13 | 0.2876 | 3.24e-03 | 5 |
| 3 | 177 600 | 23 | 0.2191 | 2.12e-03 | 5 |
| 4 | 710 400 | 13 | 0.2898 | 6.06e-04 | 21 |
| 4 | 710 400 | 23 | 0.2222 | 6.85e-04 | 21 |

The change in the stress along the interfaces falls by a factor of 4 to 5 for
each halving of the cell size. Along the mid-line, the stress across the
mid-layer changes by at most 1.5e-5 of the applied stress between the two
finest meshes.

## References

S. Ioka, K. Masuda, S. Kubo, Singular stress field near the edge of interface
of bonded dissimilar materials with an interlayer, *International Journal of
Solids and Structures*, 44, 6232–6238, 2007.

D. B. Bogy, On the problem of edge-bonded elastic quarter-planes loaded at the
boundary, *International Journal of Solids and Structures*, 6, 1287–1313,
1970.

Ž. Tuković, A. Ivanković, A. Karač, Finite-volume stress analysis in
multi-material linear elastic body, *International Journal for Numerical
Methods in Engineering*, 93(4), 400–419, 2013.
