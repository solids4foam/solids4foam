#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Usage: ./regressionTest.sh [--check-only]
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_DIR="$SCRIPT_DIR/regressionTests/main"
CHECK_ONLY=false
for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run) CHECK_ONLY=true ;;
        *) echo "Unknown option $arg"; exit 1 ;;
    esac
done

source solids4FoamScripts.sh

if [[ "$CHECK_ONLY" == false ]]; then
    rm -rf "$CASE_DIR"
    mkdir -p "$CASE_DIR"
    for item in "$SCRIPT_DIR"/*; do
        name=$(basename "$item")
        case "$name" in
            regressionTests|postProcessing|log.*) continue ;;
        esac
        cp -a "$item" "$CASE_DIR/"
    done
    (cd "$CASE_DIR" && ./Allrun > log.Allrun 2>&1) || true
fi

if solids4Foam::regressionCaseSkipped "$CASE_DIR/log.Allrun"; then
    echo "SKIP: pressureStabilisationFourier is not supported in this OpenFOAM flavour"
    exit 0
fi

failures=0
SOLVER_LOG="$CASE_DIR/log.stabilisationFourierCheck.N16"
# The utility compares every measured eigenvalue with its analytical symbol,
# in both the raw and the spectrally normalised campaigns, and the equivalent
# operator pairs with each other, and counts the failures
if [[ -f "$SOLVER_LOG" ]] && grep -q '^Failures: 0$' "$SOLVER_LOG" \
&& grep -q '^End' "$SOLVER_LOG"; then
    echo "PASS: $(grep -c ' PASS$' "$SOLVER_LOG") eigenvalue and pair checks, 0 failures"
else
    echo "FAIL: stabilisationFourierCheck did not complete with 0 failures"
    failures=$((failures + 1))
fi
for campaign in RAW NORMALISED; do
    if ! grep -q "^$campaign Fourier check" "$SOLVER_LOG" 2>/dev/null; then
        echo "FAIL: the $campaign campaign did not run"
        failures=$((failures + 1))
    fi
done
# Field output: every pressure mode and model response is written, and every
# mode is an eigenmode of every model
LOG_FIELDS="$CASE_DIR/log.stabilisationFourierCheck.fields"
missing=""
for field in fourierP_lowX fourierP_nyquist \
    fourierS_laplacian_lowX fourierS_laplacian_nyquist \
    fourierS_gen1Normalised_lowX fourierS_gen1Normalised_nyquist; do
    if ! grep -q '^internalField' "$CASE_DIR/0/$field" 2>/dev/null; then
        missing="$missing $field"
    fi
done
if [[ -f "$LOG_FIELDS" ]] && grep -q '^Failures: 0$' "$LOG_FIELDS" \
&& grep -q '^End' "$LOG_FIELDS" && [[ -z "$missing" ]]; then
    echo "PASS: field output: 2 modes and 2 models written, all eigenmodes"
else
    echo "FAIL: field output (missing fields:${missing:- none})"
    failures=$((failures + 1))
fi
# The 3-D level: three active directions, referenceNyquistDirections 3
LOG_3D="$CASE_DIR/log.stabilisationFourierCheck.3D"
if [[ -f "$LOG_3D" ]] && grep -q '^Failures: 0$' "$LOG_3D" \
&& grep -q '^End' "$LOG_3D" && grep -q '^NORMALISED Fourier check' "$LOG_3D"; then
    echo "PASS: 3-D level: $(grep -c ' PASS$' "$LOG_3D") eigenvalue and pair checks, 0 failures"
else
    echo "FAIL: the 3-D level did not complete with 0 failures"
    failures=$((failures + 1))
fi
nRejected=$(grep -c '^PASS rejection' "$CASE_DIR/log.Allrun" || true)
if [[ "$nRejected" -eq 9 ]] && ! grep -q '^FAIL rejection' "$CASE_DIR/log.Allrun"; then
    echo "PASS: 9 invalid configurations rejected with their messages"
else
    echo "FAIL: parameter validation ($nRejected of 9 rejected as expected)"
    failures=$((failures + 1))
fi
# A renamed, reversed pair differing in the last bit fails as it must
if grep -q '^PASS bitwise pair' "$CASE_DIR/log.Allrun"; then
    echo "PASS: a renamed, reversed pair must agree bitwise"
else
    echo "FAIL: the bitwise requirement of a renamed, reversed pair"
    failures=$((failures + 1))
fi
# The keyword spectralNormalisation is accepted
if grep -q '^PASS accepted: acceptSpectralNormalisation' "$CASE_DIR/log.Allrun"; then
    echo "PASS: spectralNormalisation accepted"
else
    echo "FAIL: spectralNormalisation was not accepted"
    failures=$((failures + 1))
fi

if ((failures)); then
    echo "Regression test FAILED ($failures checks)"
    exit 1
fi
echo "Regression test PASSED"
