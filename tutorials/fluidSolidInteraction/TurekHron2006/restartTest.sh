#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

# ============================================================
# Restart consistency of the quasi-monolithic coupling
# ============================================================
# A run restarted from a written time must continue as the uninterrupted run
# does (solids4foam #740). The quasi-monolithic FSI3 case is run in parallel to
# END_TIME, and again from RESTART_TIME:
#   - on the same decomposition: the tip displacement must agree with the
#     continuous run to SAME_REL_TOL (the solver tolerance), as the fields
#     are written with full precision;
#   - decomposed again region by region (decomposePar -region), which does not
#     copy <time>/uniform: solids4Foam must stop with an error rather than
#     continue as a new run;
#   - decomposed again for all regions in one call with -copyUniform: the
#     result must agree to REDECOMP_REL_TOL (the decomposition changes the
#     Krylov iterates, not the converged solution)
# Needs OpenFOAM.com and PETSc; it is skipped otherwise.

RESTART_TIME=0.02
END_TIME=0.03
N_PROCS=2
N_PROCS_REDECOMP=3
SAME_REL_TOL=1e-6
REDECOMP_REL_TOL=1e-5

TEST_ROOT="${SCRIPT_DIR}/regressionTests/restart"

if [[ "${WM_PROJECT_VERSION:-}" != v* || -z "${PETSC_DIR:-}" ]]; then
    echo "Skipping the restart test: it needs OpenFOAM.com and PETSc"
    exit 0
fi

# GNU sed, for the in-place edits below
solids4Foam::requireGnuSed

echo "============================================================"
echo "Quasi-monolithic restart consistency test"
echo "Restart time                   = ${RESTART_TIME}"
echo "End time                       = ${END_TIME}"
echo "Same decomposition tolerance   < ${SAME_REL_TOL}"
echo "Re-decomposition tolerance     < ${REDECOMP_REL_TOL}"
echo "============================================================"
echo

set_subdomains() {
    local nProcs="$1"
    local dict
    for dict in system/decomposeParDict system/fluid/decomposeParDict \
        system/solid/decomposeParDict
    do
        "${SOLIDS4FOAM_SED}" -i \
            -e "s/^\(numberOfSubdomains[[:space:]]*\).*/\1${nProcs};/" \
            -e "s/^\(method[[:space:]]*\).*/\1scotch;/" \
            "${dict}"
    done
}

# Remove the processor time directories later than the restart time
remove_later_times() {
    local dir time
    for dir in processor*/*
    do
        time=$(basename "${dir}")
        if [[ "${time}" =~ ^[0-9.eE+-]+$ ]] \
        && awk "BEGIN {exit !(${time} > ${RESTART_TIME} + 1e-12)}"
        then
            rm -rf "${dir}"
        fi
    done
}

# Final tip displacement magnitude, from the history written by a run that
# started at START (0 or RESTART_TIME); empty unless that run reached END_TIME
final_tip_disp() {
    local start="$1"
    local file="postProcessing/${start}/solidPointDisplacement_pointDisp.dat"
    [[ -f "${file}" ]] || return 0
    awk -v e="${END_TIME}" '($1 + 0) == $1 {t = $1; d = $5}
        END {if (d != "" && t >= e - 1e-9) print d}' "${file}"
}

prepare_case() {
    local caseDir="$1"
    rm -rf "${caseDir}"
    mkdir -p "${caseDir}"
    local item base_item
    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" \
           || "${base_item}" == "verification" ]]; then
            continue
        fi
        cp -a "${item}" "${caseDir}/"
    done
    (
        cd "${caseDir}"
        # A copy of a tutorial that has been run holds its output
        ./Allclean > /dev/null 2>&1
        solids4Foam::restoreCaseVariant .
        ln -nsf fsiProperties.monolithic constant/fsiProperties
        solids4Foam::selectCaseVariant . monolithic
        solids4Foam::convertCaseFormat .
        "${SOLIDS4FOAM_SED}" -i \
            -e "s/^\(endTime[[:space:]]*\).*/\1${END_TIME};/" \
            -e "s/^\(writeControl[[:space:]]*\).*/\1timeStep;/" \
            -e "s/^\(writeInterval[[:space:]]*\).*/\110;/" \
            -e "s/^\(writePrecision[[:space:]]*\).*/\116;/" \
            system/controlDict
        # Couple from the start, so that the plate and the fluid mesh move
        # before the restart time (the benchmark holds the plate fixed until
        # t = 2), with a gently ramped inflow so that the coupled start-up
        # stays well conditioned
        "${SOLIDS4FOAM_SED}" -i \
            -e "s/^\([[:space:]]*maxValue[[:space:]]*\).*/\10.3;/" \
            -e "s/^\([[:space:]]*transitionPeriod[[:space:]]*\).*/\10.05;/" \
            0/fluid/U
        "${SOLIDS4FOAM_SED}" -i \
            -e "s/^\([[:space:]]*coupled[[:space:]]*\).*/\1yes;/" \
            -e "/^[[:space:]]*couplingStartTime[[:space:]]/d" \
            constant/fsiProperties
        # Both solids write the fields a restart needs
        local props
        for props in constant/solid/solidProperties \
            constant/fluid/solidProperties
        do
            "${SOLIDS4FOAM_SED}" -i \
                -e '/Coeffs$/{n;s/{/{\n    restart yes;/}' "${props}"
        done
        set_subdomains "${N_PROCS}"
        blockMesh -region solid > log.blockMesh.solid 2>&1
        blockMesh -region fluid > log.blockMesh.fluid 2>&1
    ) > /dev/null
}

failures=0

check_rel() {
    local label="$1" value="$2" reference="$3" tolerance="$4"
    if awk -v v="${value}" -v r="${reference}" -v t="${tolerance}" \
        'BEGIN {d = v - r; if (d < 0) d = -d; a = (r < 0 ? -r : r);
         exit !(d <= t*a + 1e-30)}'
    then
        printf "PASS: %s tip displacement = %.12g (continuous %.12g)\n" \
            "${label}" "${value}" "${reference}"
    else
        printf "FAIL: %s tip displacement = %.12g (continuous %.12g)\n" \
            "${label}" "${value}" "${reference}"
        failures=$((failures + 1))
    fi
}

# Continuous run
CONT="${TEST_ROOT}/continuous"
prepare_case "${CONT}"
(
    cd "${CONT}"
    decomposePar -regions '(fluid solid)' > log.decomposePar 2>&1
    solids4Foam::runParallel -o solids4Foam > /dev/null
)
ref_disp=$(cd "${CONT}" && final_tip_disp 0)
if [[ -z "${ref_disp}" ]]; then
    echo "FAIL: the continuous run did not complete (see ${CONT}/log.solids4Foam)"
    exit 1
fi

# Restart on the same decomposition
SAME="${TEST_ROOT}/sameDecomposition"
rm -rf "${SAME}"
cp -a "${CONT}" "${SAME}"
(
    cd "${SAME}"
    rm -rf postProcessing
    remove_later_times
    "${SOLIDS4FOAM_SED}" -i \
        "s/^\(startFrom[[:space:]]*\).*/\1latestTime;/" system/controlDict
    solids4Foam::runParallel -o solids4Foam > /dev/null
)
check_rel "same decomposition" "$(cd "${SAME}" && final_tip_disp "${RESTART_TIME}")" \
    "${ref_disp}" "${SAME_REL_TOL}"

# Reconstruct the restart time, for the two re-decompositions
reconstruct_restart_time() {
    local caseDir="$1"
    rm -rf "${caseDir}"
    cp -a "${CONT}" "${caseDir}"
    (
        cd "${caseDir}"
        rm -rf postProcessing
        reconstructPar -region fluid -time "${RESTART_TIME}" \
            > log.reconstructPar.fluid 2>&1
        reconstructPar -region solid -time "${RESTART_TIME}" \
            > log.reconstructPar.solid 2>&1
        rm -rf processor*
        set_subdomains "${N_PROCS_REDECOMP}"
        "${SOLIDS4FOAM_SED}" -i \
            "s/^\(startFrom[[:space:]]*\).*/\1latestTime;/" system/controlDict
    )
}

# Re-decomposed region by region: <time>/uniform is not copied, and the run
# must stop rather than continue as a new run
BYREGION="${TEST_ROOT}/redecomposedByRegion"
reconstruct_restart_time "${BYREGION}"
(
    cd "${BYREGION}"
    decomposePar -region fluid -time "${RESTART_TIME}" \
        > log.decomposePar.fluid 2>&1
    decomposePar -region solid -time "${RESTART_TIME}" \
        > log.decomposePar.solid 2>&1
    solids4Foam::runParallel -o solids4Foam > /dev/null || true
)
if grep -q "uniform/time is missing" "${BYREGION}/log.solids4Foam"; then
    echo "PASS: region-by-region re-decomposition stops with the missing"
    echo "      uniform/time error"
else
    echo "FAIL: region-by-region re-decomposition did not stop with the"
    echo "      missing uniform/time error (see ${BYREGION}/log.solids4Foam)"
    failures=$((failures + 1))
fi

# Re-decomposed for all regions in one call, copying <time>/uniform
ALLREGIONS="${TEST_ROOT}/redecomposedAllRegions"
reconstruct_restart_time "${ALLREGIONS}"
(
    cd "${ALLREGIONS}"
    decomposePar -regions '(fluid solid)' -copyUniform \
        -time "${RESTART_TIME}" > log.decomposePar 2>&1
    solids4Foam::runParallel -o solids4Foam > /dev/null
)
check_rel "re-decomposed with uniform" \
    "$(cd "${ALLREGIONS}" && final_tip_disp "${RESTART_TIME}")" "${ref_disp}" "${REDECOMP_REL_TOL}"

echo
if (( failures == 0 )); then
    echo "Restart test PASSED"
    exit 0
else
    echo "Restart test FAILED (${failures} checks)"
    exit 1
fi
