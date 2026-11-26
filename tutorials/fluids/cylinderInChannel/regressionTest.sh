#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# GNU sed, for the in-place edits below
solids4Foam::requireGnuSed

# ============================================================
# cylinderInChannel regression test
# main: the force history as a cheap fluid benchmark check.
# shortOutlet, shortOutletNewton: the channel cut 4 D behind the
# cylinder, so that the starting wake crosses the outlet, with the
# backflowStabilisedOutletVelocity outlet (beta = 0.2), run with
# pimpleFluid and, when PETSc is available, newtonIcoFluid. The
# condition damps the kinetic energy carried in by backflow, which
# the pimpleFluid arm checks; with newtonIcoFluid, which builds the
# outlet flux from the boundary velocity, it also reduces the peak
# outlet backflow (0.81% with zero gradient)
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

REG_END_TIME=5
DRAG_MIN=0.622
DRAG_MAX=0.624
LIFT_MIN=0.005
LIFT_MAX=0.012

SHORT_END_TIME=8
SHORT_DRAG_MIN=0.623
SHORT_DRAG_MAX=0.626
SHORT_LIFT_MIN=0.020
SHORT_LIFT_MAX=0.030

# The newtonIcoFluid on development is not the latest version, which is on
# feature-petsc-snes-quasi-monolithic. This version gives a different
# solution from pimpleFluid for this case whatever the outlet condition
# (final drag 0.6198 and lift 0.0182 at t = 5 in the full channel), so its
# short-outlet forces have their own bands. With the latest newtonIcoFluid,
# use the pimpleFluid force bands and a peak backflow band of
# [0.0080, 0.0089] (0.847%; zero gradient 0.92%). The lower lift limit
# covers foam-extend (0.0476; 0.0529 with OpenFOAM), where the wake does not
# reverse the outlet flow (issue #733)
NEWTON_DRAG_MIN=0.621
NEWTON_DRAG_MAX=0.624
NEWTON_LIFT_MIN=0.045
NEWTON_LIFT_MAX=0.058

# pimpleFluid: peak kinetic energy influx through the outlet, relative to
# the outflux, -sum(0.5*|U|^2*min(phi, 0))/sum(0.5*|U|^2*max(phi, 0)), over
# the written times: 1.7e-05 (OpenFOAM-v2512), against 3.3e-05 with zero
# gradient; beta = 0.5 and 1 give 1.0e-05 and 0.5e-05. The traction does
# not reduce the backflow flux itself with pimpleFluid (0.99%, against
# 0.96% with zero gradient), so that is not checked.
# newtonIcoFluid: peak backflow flux, -sum(phi < 0)/sum(phi > 0), 0.747%,
# against 0.812% with zero gradient (OpenFOAM-v2512).
# With foam-extend the wake does not reverse the outlet flow, for either
# outlet or fluid model, so the backflow checks are skipped there (issue
# #733)
SHORT_ENERGY_BACKFLOW_MIN=1.2e-05
SHORT_ENERGY_BACKFLOW_MAX=2.5e-05
NEWTON_BACKFLOW_MIN=0.0070
NEWTON_BACKFLOW_MAX=0.0078
NEWTON_MAX_ITERATIONS=6

ALLRUN_LOGFILE="log.Allrun"
FOAM_FLAVOUR=$(solids4Foam::foamFlavour)

echo "============================================================"
echo "cylinderInChannel regression test"
echo "main: end time ${REG_END_TIME}," \
    "final drag in [${DRAG_MIN}, ${DRAG_MAX}]," \
    "final lift in [${LIFT_MIN}, ${LIFT_MAX}]"
echo "shortOutlet: end time ${SHORT_END_TIME}," \
    "final drag in [${SHORT_DRAG_MIN}, ${SHORT_DRAG_MAX}]," \
    "final lift in [${SHORT_LIFT_MIN}, ${SHORT_LIFT_MAX}]," \
    "peak energy backflow in" \
    "[${SHORT_ENERGY_BACKFLOW_MIN}, ${SHORT_ENERGY_BACKFLOW_MAX}]"
echo "shortOutletNewton: end time ${SHORT_END_TIME}," \
    "final drag in [${NEWTON_DRAG_MIN}, ${NEWTON_DRAG_MAX}]," \
    "final lift in [${NEWTON_LIFT_MIN}, ${NEWTON_LIFT_MAX}]," \
    "peak backflow in" \
    "[${NEWTON_BACKFLOW_MIN}, ${NEWTON_BACKFLOW_MAX}]," \
    "at most ${NEWTON_MAX_ITERATIONS} Newton iterations per step"
echo "============================================================"
echo

prepare_case() {
    local case_dir=$1
    local end_time=$2
    local item base_item

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done

    "${SOLIDS4FOAM_SED}" -i \
        "s/^endTime[[:space:]]\+50;/endTime         ${end_time};/" \
        "${case_dir}/system/controlDict"
}

find_force_file() {
    local case_dir=$1
    local candidate
    for candidate in \
        "${case_dir}/forces/0/forces.dat" \
        "${case_dir}/forces/0/force.dat" \
        "${case_dir}/postProcessing/forces/0/force.dat" \
        "${case_dir}/postProcessing/forces/0/forces.dat" \
        "${case_dir}/postProcessing/fluid/forces/0/force.dat" \
        "${case_dir}/postProcessing/fluid/forces/0/forces.dat"
    do
        if [[ -f "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done
    return 1
}

# Final total force component: 1 for drag, 2 for lift
final_force() {
    local force_file=$1
    local component=$2

    awk -v c="${component}" '
        END {
            gsub(/[()]/, "", $0)

            if (NF >= 13)
            {
                print $(1 + c) + $(4 + c)
            }
            else if (NF >= 9)
            {
                print $(1 + c)
            }
        }
    ' "${force_file}"
}

# Backflow fraction, -sum(phi < 0)/sum(phi > 0), of the outlet patch in a
# phi file; reads the uniform and nonuniform forms
outlet_backflow() {
    local phi_file=$1

    awk -v patch=outlet '
        function add(v)
        {
            if (v + 0 > 0) { pos += v } else { neg -= v }
        }
        BEGIN { inBf = 0; inPatch = 0; state = 0; found = 0 }
        /^boundaryField/ { inBf = 1; next }
        inBf && !inPatch && $1 == patch { inPatch = 1; next }
        inPatch && state == 0 && $1 == "value" {
            if ($2 == "uniform")
            {
                v = $3; sub(/;$/, "", v); add(v)
                found = 1; exit
            }
            if (index($0, "("))
            {
                line = substr($0, index($0, "(") + 1)
                sub(/\).*$/, "", line)
                n = split(line, a, " ")
                for (i = 1; i <= n; i++) { add(a[i]) }
                found = 1; exit
            }
            state = 1; next
        }
        state == 1 && /^ *\( *$/ { state = 2; next }
        state == 2 && /^ *\)/ { found = 1; exit }
        state == 2 { add($1) }
        END {
            if (!found || pos <= 0) { exit 1 }
            printf "%.6g\n", neg/pos
        }
    ' "${phi_file}"
}

# Kinetic energy influx through the outlet patch relative to the outflux,
# -sum(0.5*|U|^2*min(phi, 0))/sum(0.5*|U|^2*max(phi, 0)), from the phi and
# U files of a time directory; reads the uniform and nonuniform forms
outlet_energy_backflow() {
    local time_dir=$1

    awk -v patch=outlet '
        # Append the numbers of the outlet value entry of the current file
        # to vals, setting uni if the entry is uniform
        FNR == 1 {
            f++; inBf = 0; inPatch = 0; state = 0; text = ""
        }
        state == 3 { next }
        /^boundaryField/ { inBf = 1; next }
        inBf && !inPatch && $1 == patch { inPatch = 1; next }
        inPatch && state == 0 && $1 == "value" {
            if ($2 == "uniform")
            {
                uni[f] = 1
                text = substr($0, index($0, "uniform") + 7)
                state = 3
                store(f, text)
                next
            }
            text = $0
            if (index($0, ";")) { state = 3; store(f, list(text)); next }
            state = 1; next
        }
        state == 1 {
            text = text " " $0
            if ($0 ~ /^ *; *$/ || $0 ~ /\) *; *$/)
            {
                state = 3; store(f, list(text))
            }
            next
        }
        function list(t)
        {
            sub(/^[^(]*\(/, "", t)
            return t
        }
        function store(i, t,    n, a, k)
        {
            gsub(/[();]/, " ", t)
            n = split(t, a, " ")
            nv[i] = 0
            for (k = 1; k <= n; k++) { vals[i, ++nv[i]] = a[k] }
        }
        END {
            # File 1: phi (scalars), file 2: U (vectors)
            if (nv[1] == 0 || nv[2] == 0) { exit 1 }
            nFaces = uni[1] ? nv[2]/3 : nv[1]
            if (uni[1] && uni[2]) { nFaces = 1 }
            for (i = 1; i <= nFaces; i++)
            {
                phi = uni[1] ? vals[1, 1] : vals[1, i]
                j = uni[2] ? 0 : 3*(i - 1)
                ke = 0.5*(vals[2, j + 1]^2 + vals[2, j + 2]^2 \
                    + vals[2, j + 3]^2)
                if (phi + 0 < 0) { eIn -= ke*phi } else { eOut += ke*phi }
            }
            if (eOut <= 0) { exit 1 }
            printf "%.6g\n", eIn/eOut
        }
    ' "${time_dir}/phi" "${time_dir}/U"
}

# Peak over the written times of a case of a per-time measure: the
# arguments are the case directory and outlet_backflow or
# outlet_energy_backflow
peak_over_times() {
    local case_dir=$1
    local measure=$2
    local peak=""
    local dir name value

    for dir in "${case_dir}"/*; do
        [[ -f "${dir}/phi" && -f "${dir}/U" ]] || continue
        name="${dir##*/}"
        [[ "${name}" =~ ^[0-9.eE+-]+$ ]] || continue
        if [[ "${measure}" == "outlet_backflow" ]]; then
            value=$(outlet_backflow "${dir}/phi") || return 1
        else
            value=$(outlet_energy_backflow "${dir}") || return 1
        fi
        if [[ -z "${peak}" ]] \
            || awk "BEGIN {exit !(${value} > ${peak})}"; then
            peak=${value}
        fi
    done

    [[ -n "${peak}" ]] || return 1
    echo "${peak}"
}

# Largest number of Newton iterations in a time step, from the PETSc SNES
# monitor
max_newton_iterations() {
    local log_file=$1

    awk '
        /^Time = / { if (n > max) max = n; n = 0 }
        /SNES Function norm/ { n++ }
        END { if (n > max) max = n; print max - 1 }
    ' "${log_file}"
}

failures=0

check_range() {
    local label=$1
    local value=$2
    local min=$3
    local max=$4

    if [[ -z "${value}" ]]; then
        echo "FAIL: ${label}: no value"
        failures=$((failures + 1))
    elif awk "BEGIN {exit !(${value} >= ${min} && ${value} <= ${max})}"; then
        printf "PASS: %s = %.6g\n" "${label}" "${value}"
    else
        printf "FAIL: %s = %.6g, outside [%s, %s]\n" \
            "${label}" "${value}" "${min}" "${max}"
        failures=$((failures + 1))
    fi
}

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------

CHECK_ONLY=false

for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run)
            CHECK_ONLY=true
            ;;
        *)
            ;;
    esac
done

if [ "$CHECK_ONLY" = true ]; then
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# Run a case: the arguments are the case name, the end time and the Allrun
# arguments. Returns 1 when the tutorial skipped in this environment
run_case() {
    local name=$1
    local end_time=$2
    shift 2
    local case_dir="${REGRESSION_ROOT}/${name}"

    if [ "$CHECK_ONLY" = false ]; then
        prepare_case "${case_dir}" "${end_time}"
        ( cd "${case_dir}" && ./Allclean > /dev/null 2>&1 ) || true
        if ! ( cd "${case_dir}" && ./Allrun "$@" > "${ALLRUN_LOGFILE}" 2>&1 )
        then
            echo "FAIL: ${name}: Allrun failed;" \
                "see ${case_dir}/${ALLRUN_LOGFILE}"
            failures=$((failures + 1))
        fi
    fi

    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"
    then
        echo "SKIP: ${name}: the tutorial skipped in this environment"
        return 1
    fi

    return 0
}

check_forces() {
    local name=$1
    local drag_min=$2
    local drag_max=$3
    local lift_min=$4
    local lift_max=$5
    local case_dir="${REGRESSION_ROOT}/${name}"
    local force_file=""

    if ! force_file=$(find_force_file "${case_dir}"); then
        echo "FAIL: ${name}: could not find force history output"
        failures=$((failures + 1))
        return 0
    fi

    check_range "${name}: final drag" \
        "$(final_force "${force_file}" 1)" "${drag_min}" "${drag_max}"
    check_range "${name}: final lift" \
        "$(final_force "${force_file}" 2)" "${lift_min}" "${lift_max}"
}

check_end_time() {
    local name=$1
    local end_time=$2
    local case_dir="${REGRESSION_ROOT}/${name}"
    local latest

    latest=$(solids4Foam::latestTime "${case_dir}")
    if [[ -n "${latest}" ]] \
        && awk "BEGIN {exit !(${latest} == ${end_time})}"; then
        echo "PASS: ${name}: reached t = ${latest}"
    else
        echo "FAIL: ${name}: latest time '${latest}', expected ${end_time}"
        failures=$((failures + 1))
    fi
}

clean_case() {
    local name=$1

    if [ "$CHECK_ONLY" = false ]; then
        ( cd "${REGRESSION_ROOT}/${name}" && ./Allclean > /dev/null 2>&1 ) \
            || true
    fi
}

# ------------------------------------------------------------
# main: full channel
# ------------------------------------------------------------

if run_case main "${REG_END_TIME}"; then
    check_forces main "${DRAG_MIN}" "${DRAG_MAX}" "${LIFT_MIN}" "${LIFT_MAX}"
    clean_case main
fi

# ------------------------------------------------------------
# shortOutlet: backflow through the outlet, pimpleFluid
# ------------------------------------------------------------

if run_case shortOutlet "${SHORT_END_TIME}" shortOutlet pimpleFluid; then
    check_end_time shortOutlet "${SHORT_END_TIME}"
    check_forces shortOutlet "${SHORT_DRAG_MIN}" "${SHORT_DRAG_MAX}" \
        "${SHORT_LIFT_MIN}" "${SHORT_LIFT_MAX}"
    if [[ "${FOAM_FLAVOUR}" == "foamextend" ]]; then
        echo "SKIP: shortOutlet: peak outlet energy backflow: no backflow" \
            "with foam-extend (issue #733)"
    else
        check_range "shortOutlet: peak outlet energy backflow" \
            "$(peak_over_times "${REGRESSION_ROOT}/shortOutlet" \
            outlet_energy_backflow || true)" \
            "${SHORT_ENERGY_BACKFLOW_MIN}" "${SHORT_ENERGY_BACKFLOW_MAX}"
    fi
    clean_case shortOutlet
fi

# ------------------------------------------------------------
# shortOutletNewton: as shortOutlet, with newtonIcoFluid (PETSc)
# ------------------------------------------------------------

NEWTON_CASE_DIR="${REGRESSION_ROOT}/shortOutletNewton"

if run_case shortOutletNewton "${SHORT_END_TIME}" \
    shortOutlet newtonIcoFluid; then
    check_end_time shortOutletNewton "${SHORT_END_TIME}"
    check_forces shortOutletNewton "${NEWTON_DRAG_MIN}" \
        "${NEWTON_DRAG_MAX}" "${NEWTON_LIFT_MIN}" "${NEWTON_LIFT_MAX}"
    if [[ "${FOAM_FLAVOUR}" == "foamextend" ]]; then
        echo "SKIP: shortOutletNewton: peak outlet backflow: no backflow" \
            "with foam-extend (issue #733)"
    else
        check_range "shortOutletNewton: peak outlet backflow" \
            "$(peak_over_times "${NEWTON_CASE_DIR}" outlet_backflow \
            || true)" \
            "${NEWTON_BACKFLOW_MIN}" "${NEWTON_BACKFLOW_MAX}"
    fi
    check_range "shortOutletNewton: max Newton iterations per step" \
        "$(max_newton_iterations "${NEWTON_CASE_DIR}/log.solids4Foam")" \
        1 "${NEWTON_MAX_ITERATIONS}"
    clean_case shortOutletNewton
fi

echo
if (( failures == 0 )); then
    echo "============================================================"
    echo "Regression test PASSED"
    echo "============================================================"
    exit 0
else
    echo "============================================================"
    echo "Regression test FAILED (${failures} checks)"
    echo "============================================================"
    exit 1
fi
