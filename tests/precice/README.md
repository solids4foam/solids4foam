# preCICE coupling tests

solids4foam can be coupled to OpenFOAM fluid solvers, and to other solvers,
through [preCICE](https://precice.org). The coupling is provided by the generic
[OpenFOAM preCICE adapter](https://github.com/precice/openfoam-adapter), loaded
as a function object; there is no preCICE code in solids4foam itself.

The cases that exercise this live upstream in
[precice/tutorials](https://github.com/precice/tutorials), not in this
repository. The tests here clone that repository at a pinned commit, run the
registered participant combinations, and compare the results against stored
reference values.

The preCICE team already run these cases in their own system tests, but those
pin a *released* solids4foam Docker image, so they cannot catch a regression
introduced on a branch. These tests are deliberate redundancy on top of that,
giving per-branch coverage ahead of a release.

## Running locally

You need preCICE and the OpenFOAM preCICE adapter installed, with
`libpreciceAdapterFunctionObject.so` in `$FOAM_USER_LIBBIN`. See the
[preCICE quickstart](https://precice.org/quickstart.html) and the
[adapter documentation](https://precice.org/adapter-openfoam-overview.html).

```bash
./tests/precice/Alltest-precice              # all registered cases
./tests/precice/Alltest-precice perpendicular-flap   # one case
./tests/precice/Alltest-precice --keep       # reuse an existing tutorials clone
```

If preCICE or the adapter is missing, the script prints a message and exits
successfully, so it is safe to run anywhere.

The tutorials clone and all logs are written to `../tutorialsTest-precice`,
alongside the `../tutorialsTest` and `../tutorialsTest-regression` directories
used by the other test scripts.

## Running in CI

`.github/workflows/preciceTest.yml` builds preCICE, the adapter and solids4foam
in a container, then runs `Alltest-precice`. Building solids4foam dominates the
run time, so the job does not run on every pull request. It runs:

- always, for pull requests targeting `master`, which precede a release;
- on any pull request carrying the `test-precice` label — remove and re-add the
  label to trigger it again, the same convention the preCICE adapter
  repositories use with their `trigger-system-tests` label;
- manually, via `gh workflow run preciceTest.yml`.

## Layout

```text
tests/precice/
  Alltest-precice          driver
  generateReferences.sh    prints the values the checks compare against
  versions.txt             pinned preCICE, adapter and tutorials versions
  <case>/caseSetup         which tutorial and participants this case runs
  <case>/regressionTest.sh the checks, with reference values and tolerances
```

Every immediate subdirectory containing a `caseSetup` is a registered case, so
the driver needs no central list.

## Adding a case

When a solids4foam participant lands in another upstream preCICE tutorial:

1. bump `TUTORIALS_REF` in `versions.txt` to a commit containing it;
2. add `tests/precice/<tutorial>/caseSetup`, declaring `TUTORIAL`,
   `PARTICIPANTS`, `MAX_TIME` and `TIMEOUT`;
3. add `tests/precice/<tutorial>/regressionTest.sh`, following
   `perpendicular-flap/regressionTest.sh`;
4. generate and commit its reference values, as below.

Choose `MAX_TIME` so that the truncated run still spans behaviour worth
checking. For `perpendicular-flap`, 1.0 s covers the first oscillation peak and
the reversal, which makes the checks sensitive to a period or phase error; a
shorter window would only see a monotonic rise.

No change to the driver or the workflow is needed.

## Updating the reference values

The reference values are tied to the versions in `versions.txt` and to the
upstream case configuration, so they must be regenerated whenever either
changes deliberately.

```bash
./tests/precice/generateReferences.sh
```

This runs the cases and prints the quantities the checks compare against. It
does not edit the `regressionTest.sh` scripts: read the values, run it twice to
see the run-to-run spread, and set the tolerances by hand with headroom above
that spread.

Generate the committed values in the same container CI uses, since results can
differ between OpenFOAM versions:

```bash
docker run --rm --user 1001 -e HOME=/work -v "$PWD:/work/ws" -w /work \
  philippic/openfoam-v2512:ubuntu-22.04-petsc-no-openmp-system-blas-gmsh \
  bash -lc '...'   # see .github/workflows/preciceTest.yml for the full steps
```

## What is checked

For each case, the preCICE watch-point output and the implicit-coupling
iteration log. For `perpendicular-flap` that is the flap-tip x-displacement —
the quantity used to compare this tutorial's solid participants against each
other, see [precice/tutorials#909](https://github.com/precice/tutorials/pull/909)
— plus the total number of coupling iterations, which catches a convergence
regression that leaves the displacement history close enough to pass.

This is deliberately a small set of scalars rather than a field-by-field
comparison. It matches the idiom of `tutorials/Alltest-regression` and the
per-case `regressionTest.sh` scripts used elsewhere in this repository. The
upstream system tests already do the full `fieldcompare` comparison against
Git-LFS reference results.

A case that does not reach its requested end time is reported as a skip rather
than a failure, matching the convention of the other regression tests; the
driver reports the participant exit status separately.
