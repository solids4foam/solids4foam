# Robin coupling study: choosing the added-mass coefficient `hs`

This directory holds the scripts used to investigate how the Robin coefficient
`rho_s*hs` of the `elasticWallPressure` boundary condition can be chosen
automatically, robustly and (close to) optimally. It is not a runnable case:
the scripts copy existing FSI tutorials, modify them and run them in
`tutorialsTest-robinHs/` at the repository root (ignored by git).

## Background

The Robin-Neumann coupling solves the fluid with

    p + (rho_s hs/rho_f) dp/dn = p_prev - rho_s hs a_s,n

and the solid with the resulting traction. For a Fourier mode with fluid
added-mass impedance `S_f` and time-discrete solid impedance `S_s` (both in
kg/m^2, acceleration form), each FSI iteration multiplies the error by

    rate(alpha) = (S_f/S_s) |S_s - alpha| / (S_f + alpha),   alpha = rho_s hs

(Badia, Nobile and Vergara 2008; Gerardo-Giorda, Nobile and Vergara 2010).
Therefore:

- the optimum is `alpha = S_s`, the solid's own interface impedance;
- for strong added mass (`S_f >> S_s`) the rate is about `|1 - alpha/S_s|`:
  `alpha > 2 S_s` diverges while `alpha < S_s` only slows the iterations;
- `S_s` depends on the mode: soft (bending) modes have a small impedance,
  local (grid-scale) and stiff modes a large one, so a single `alpha` must
  balance them.

For an implicit time step the solid through-thickness response decays over the
length `l = c_p dt_eff`. A slab of thickness `H` with a free back face has the
impedance `rho_s l tanh(H/l)`, which tends to the half-space value `rho_s l`
(the original default `hs = c_p dt_eff`) for thick walls and to the wall
inertia `rho_s H` for thin walls. Stiffness adds `dt_eff^2 K`, and for lateral
wavenumbers `k` the impedance grows to about `rho_s l sqrt(1 + (k l)^2)`.

## Models available in `elasticWallPressure`

`hsModel` selects the coefficient:

- `pWaveSpeed`: `hs = c_p dt_eff`, with `c_p = sqrt(impK/rho_s)`;
- `constant`: `hs = constantHs`, or a nonuniform `hs` field;
- `thicknessLimited`: `hs = l tanh(H/l)`, where `H` is the local wall
  thickness found by ray casting through the solid, halved for walls wetted on
  both sides;
- `secant` (default): starts from `seedModel` (default `thicknessLimited`)
  and rescales it using the solid and fluid impedances measured from the
  changes in interface pressure and acceleration between FSI iterations.

Secant settings: `secantUpdate` (`iteration` or `timeStep`), `secantFit`
(`minimax` over the recorded impedance samples, or the `median` of the
per-iteration fits), `secantMemory`, `secantMinFactor`, `secantMaxFactor`,
`secantMaxIncrease`, `secantMaxChange`, `secantSpatial`. Other settings:
`hsScale`, `thicknessBlend`, `twoSidedHalving`, `waveSpeed`,
`divergenceSafeguard` (a-priori models only), `writeDiagnostics`.

With `writeDiagnostics yes` the boundary condition writes
`postProcessing/robinCoefficient_<patch>.dat` with, for each FSI iteration, the
coefficient, the measured solid and fluid impedances, and the predicted and
observed contraction factors.

## Scripts

- `scripts/robin_hs_study.py`: study driver.
  - `run --case <case> --variant "tag:key=value,..." [--common ...] [--jobs N]`
    copies and runs variants of a tutorial;
  - `summarize --case <case>` tabulates the FSI iteration statistics;
  - `compare --case <case> --ref <tag>` compares the output time series
    against a reference run;
  - `report` writes Markdown tables for all cases.
- `scripts/robin_toy_model.py`: Fourier model of the iteration for the
  tutorial parameters (predicted optimum and divergence limit).

Cases: `3dTube`, `beamInCrossFlow` (modified form), `beamInCrossFlowOriginal`,
`fillingElasticContainer`, `cerebralAneurysm`, and Robin variants of
`HronTurekFsi3` (`HronTurekFsi3Robin`) and `flexibleDamBreak`
(`flexibleDamBreakRobin`).

Example:

    cd scripts
    secant="bc.hsModel=secant,bc.secantFit=minimax,bc.secantUpdate=iteration"
    ./robin_hs_study.py run --case 3dTube --jobs 4 \
        --common "fsi.robinFluxTolerance=1,bc.writeDiagnostics=yes" \
        --variant "default" \
        --variant "tl:bc.hsModel=thicknessLimited" \
        --variant "secant:$secant"
    ./robin_hs_study.py summarize --case 3dTube

Part 1 of the study (the Robin coefficient) set `robinFluxTolerance 1`: with
the original interface flux formulation the leakage levels off at a
case-dependent floor, so with its default tolerance every time step runs to
`nOuterCorr` whatever `hs` is. Part 2 fixes the leakage and revises the
convergence criteria.

## Results

See `RESULTS.md`.
