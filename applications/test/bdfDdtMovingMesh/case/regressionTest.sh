#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# BDF ddt scheme on a moving mesh
# Runs Test-bdfDdtMovingMesh, which moves the points of a box every step so
# that the cell volumes are known in closed form: the observed order of BDF 1
# to 6 from four time steps, exactness for a polynomial motion and field,
# and, with the mesh moved again after the last ddt call of each step, the
# implicit form, the geometric conservation law and the weights returned by
# fieldWeights; the conservation law also with meshPhi first called late in
# the run. See README.md
# ============================================================

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
APP="Test-bdfDdtMovingMesh"

# Source solids4Foam scripts. From PATH, as a tutorial's Allrun does, since
# this case is run both in place and from a copy at a different depth
source solids4FoamScripts.sh

echo "============================================================"
echo "BDF ddt scheme on a moving mesh"
echo "============================================================"

# The BDF ddt scheme is only built for OpenFOAM.com
solids4Foam::caseDoesNotRunWithFoamExtend
solids4Foam::caseDoesNotRunWithOpenFOAMOrg

# A skip where the application is not built, and a failure in CI, where it
# always is
solids4Foam::requireTestApp "${APP}" || exit $(( $? - 1 ))

CASE_DIR="${SCRIPT_DIR}/regressionTests/main"
rm -rf "${CASE_DIR}"
mkdir -p "${CASE_DIR}/constant"
cp -a "${SCRIPT_DIR}/system" "${CASE_DIR}/"
cd "${CASE_DIR}"

function failed()
{
    echo "FAIL: $1"
    echo "============================================================"
    echo "Regression test FAILED"
    echo "============================================================"
    exit 1
}

solids4Foam::runApplication blockMesh || failed "blockMesh"

# Convergence: four time steps, each half the previous
DELTA_TS=(0.1 0.05 0.025 0.0125)
for dt in "${DELTA_TS[@]}"
do
    solids4Foam::runApplication -s "dt${dt}" "${APP}" -deltaT "${dt}" \
        || failed "${APP} -deltaT ${dt}"
done

# Exactness for a polynomial motion and field
solids4Foam::runApplication -s poly "${APP}" -deltaT 0.05 -poly \
    || failed "${APP} -poly"

# The mesh moved again after the last ddt call of each step
MOVE_DTS=(0.05 0.025)
for dt in "${MOVE_DTS[@]}"
do
    solids4Foam::runApplication -s "moveAfterDdt.dt${dt}" "${APP}" \
        -deltaT "${dt}" -moveAfterDdt \
        || failed "${APP} -deltaT ${dt} -moveAfterDdt"
done

# meshPhi first called in the second half of the run
LATE_DTS=(0.05)
for dt in "${LATE_DTS[@]}"
do
    solids4Foam::runApplication -s "lateMeshPhi.dt${dt}" "${APP}" \
        -deltaT "${dt}" -lateMeshPhi \
        || failed "${APP} -deltaT ${dt} -lateMeshPhi"
done

LOGS=()
for dt in "${DELTA_TS[@]}"
do
    LOGS+=("log.${APP}.dt${dt}")
done
LOGS+=("log.${APP}.poly")
for dt in "${MOVE_DTS[@]}"
do
    LOGS+=("log.${APP}.moveAfterDdt.dt${dt}")
done
for dt in "${LATE_DTS[@]}"
do
    LOGS+=("log.${APP}.lateMeshPhi.dt${dt}")
done

# Each log has one line per order,
# "BDF q: [ddt error e] fvm-fvc f GCL g weights w order o".
# Checks, for each order q:
#   - the observed order from the two smallest time steps is at least
#     q - 0.5, unless the error has reached the round-off floor (1e-9);
#   - fvm-fvc and GCL are below 1e-9 in the convergence runs;
#   - the polynomial run has a ddt error, fvm-fvc and GCL below 1e-9;
#   - the runs with the mesh moved after the ddt terms have fvm-fvc and GCL
#     below 1e-9;
#   - the runs with meshPhi first called late have fvm-fvc and GCL below
#     1e-9;
#   - in every run, fvc::ddt rebuilt from the weights of fieldWeights agrees
#     with fvc::ddt to 1e-9 and the order of the final step is q
# The checks fail if a log does not have all six orders.
# The time steps as one string, for the table heading
DT_LIST=$(printf '%s ' "${DELTA_TS[@]}")

awk \
    -v nConv="${#DELTA_TS[@]}" \
    -v nMove="${#MOVE_DTS[@]}" \
    -v nLate="${#LATE_DTS[@]}" \
    -v dtList="${DT_LIST% }" \
'
function value(name,    i)
{
    for (i = 3; i < NF; i++)
    {
        if ($i == name)
        {
            return $(i + 1) + 0
        }
    }
    return 0
}

function check(ok, message)
{
    if (ok)
    {
        print "PASS: " message
        nPass++
    }
    else
    {
        print "FAIL: " message
        nFail++
    }
}

function maxOver(values, first, last, q,    f, m)
{
    m = 0
    for (f = first; f <= last; f++)
    {
        if (values[f, q] > m)
        {
            m = values[f, q]
        }
    }
    return m
}

FNR == 1 { file++ }

/^BDF [1-6]:/ {
    q = $2 + 0
    seen[file, q] = 1
    ddt[file, q] = value("error")
    fvm[file, q] = value("fvm-fvc")
    gcl[file, q] = value("GCL")
    wts[file, q] = value("weights")
    ord[file, q] = value("order")
}

END {
    tol = 1e-9
    roundOff = 1e-9
    poly = nConv + 1
    nRuns = nConv + 1 + nMove + nLate
    gsub(" ", ", ", dtList)

    printf "order | ddt error at dt = %s | observed orders\n", dtList

    for (q = 1; q <= 6; q++)
    {
        complete = 1
        for (f = 1; f <= nRuns; f++)
        {
            if (!seen[f, q])
            {
                complete = 0
            }
        }
        if (!complete)
        {
            check(0, "BDF " q ": missing from a log")
            continue
        }

        line = sprintf("  %d   |", q)
        for (f = 1; f <= nConv; f++)
        {
            line = line sprintf(" %.2e", ddt[f, q])
        }
        line = line " |"
        for (f = 1; f < nConv; f++)
        {
            p = 0
            if (ddt[f + 1, q] > 0)
            {
                p = log(ddt[f, q]/ddt[f + 1, q])/log(2)
            }
            line = line sprintf(" %.2f", p)
        }
        print line

        # p is now the order from the two smallest time steps
        check(p >= q - 0.5 || ddt[nConv, q] < roundOff, \
            sprintf("BDF %d: observed order %.2f >= %.1f", q, p, q - 0.5))

        mf = maxOver(fvm, 1, nConv, q)
        mg = maxOver(gcl, 1, nConv, q)
        check(mf < tol && mg < tol, \
            sprintf("BDF %d: fvm-fvc %.1e, GCL %.1e < %g", q, mf, mg, tol))

        check(ddt[poly, q] < tol && fvm[poly, q] < tol && gcl[poly, q] < tol, \
            sprintf("BDF %d: polynomial ddt error %.1e < %g", \
                q, ddt[poly, q], tol))

        mf = maxOver(fvm, poly + 1, poly + nMove, q)
        mg = maxOver(gcl, poly + 1, poly + nMove, q)
        check(mf < tol && mg < tol, \
            sprintf("BDF %d: moved after ddt, fvm-fvc %.1e, GCL %.1e < %g", \
                q, mf, mg, tol))

        mf = maxOver(fvm, poly + nMove + 1, nRuns, q)
        mg = maxOver(gcl, poly + nMove + 1, nRuns, q)
        check(mf < tol && mg < tol, \
            sprintf("BDF %d: meshPhi first called late, fvm-fvc %.1e, " \
                "GCL %.1e < %g", q, mf, mg, tol))

        mw = maxOver(wts, 1, nRuns, q)
        orderReached = 1
        for (f = 1; f <= nRuns; f++)
        {
            if (ord[f, q] != q)
            {
                orderReached = 0
            }
        }
        check(mw < tol && orderReached, \
            sprintf("BDF %d: fieldWeights, ddt from weights %.1e < %g, " \
                "final order %d", q, mw, tol, ord[poly, q]))
    }

    print "BDF moving-mesh checks: " (nPass + 0) " passed, " \
        (nFail + 0) " failed"
    exit (nFail > 0 || nPass != 36)
}
' "${LOGS[@]}" > log.checks || status=$?

cat log.checks

if [[ ${status:-0} -ne 0 ]]
then
    failed "$(grep -E "^BDF moving-mesh checks" log.checks || echo "checks")"
fi

echo "============================================================"
echo "Regression test PASSED"
echo "============================================================"
