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
nRejected=$(grep -c '^PASS rejection' "$CASE_DIR/log.Allrun" || true)
if [[ "$nRejected" -eq 10 ]] && ! grep -q '^FAIL rejection' "$CASE_DIR/log.Allrun"; then
    echo "PASS: 10 invalid configurations rejected with their messages"
else
    echo "FAIL: parameter validation ($nRejected of 10 rejected as expected)"
    failures=$((failures + 1))
fi
# spectralNormalisation accepted silently, its alias normalise with a warning
nAccepted=$(grep -c '^PASS accepted' "$CASE_DIR/log.Allrun" || true)
if [[ "$nAccepted" -eq 2 ]] && ! grep -q '^FAIL accepted' "$CASE_DIR/log.Allrun"; then
    echo "PASS: spectralNormalisation and its alias normalise accepted as documented"
else
    echo "FAIL: keyword spelling ($nAccepted of 2 accepted as expected)"
    failures=$((failures + 1))
fi

if ((failures)); then
    echo "Regression test FAILED ($failures checks)"
    exit 1
fi
echo "Regression test PASSED"
