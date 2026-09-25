#!/bin/bash
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Script
#     testGnuSedResolution.sh
#
# Description
#     Focused shell tests for the centralised GNU sed resolution in
#     applications/scripts/solids4FoamScripts.sh.
#
#     The tests run entirely in a temporary directory and do not require an
#     OpenFOAM installation, a build of solids4foam, or a macOS machine: the
#     macOS situation, where the system "sed" is BSD sed and GNU sed is
#     installed by Homebrew as "gsed", is emulated with stub commands placed
#     first in the PATH.
#
#     Tests:
#       1. GNU sed exposed as "sed"
#       2. GNU sed exposed only as "gsed", with a BSD-like "sed" that fails on
#          in-place editing (the macOS case)
#       3. A clear error when neither provides GNU sed
#       4. SOLIDS4FOAM_SED is always set, so "set -u" is safe
#       5. Both case-format conversion directions, under "set -u", with GNU sed
#          as "sed" and as "gsed"
#       6. The D -> DD and DD -> D rename paths in runSolidModel(), under
#          "set -u", with GNU sed as "sed" and as "gsed"
#
# Usage
#     ./testGnuSedResolution.sh
#
# License
#     GNU Lesser General Public License, version 3.
#     https://www.gnu.org/licenses/lgpl-3.0.en.html
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

set -u

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="${TEST_DIR}/../solids4FoamScripts.sh"

if [[ ! -f ${SCRIPTS} ]]
then
    echo "Cannot find ${SCRIPTS}"
    exit 1
fi

if ! sed --version 2>/dev/null | grep -q "GNU sed"
then
    echo "These tests require GNU sed to be available as 'sed' on the test"
    echo "machine, as it is used to build the 'gsed' stub"
    exit 1
fi

WORK_DIR="$(mktemp -d)"
trap 'rm -rf "${WORK_DIR}"' EXIT

REAL_SED="$(command -v sed)"

N_PASS=0
N_FAIL=0

function pass()
{
    echo "  PASS: $1"
    N_PASS=$((N_PASS + 1))
}

function fail()
{
    echo "  FAIL: $1"
    N_FAIL=$((N_FAIL + 1))
}

function check()
{
    # 1: description, 2: expected, 3: actual
    if [[ "$2" == "$3" ]]
    then
        pass "$1"
    else
        fail "$1: expected '$2' but got '$3'"
    fi
}

function checkFileContains()
{
    # 1: description, 2: file, 3: pattern
    if [[ -f "$2" ]] && grep -qF -- "$3" "$2"
    then
        pass "$1"
    else
        fail "$1: '$3' not found in '$2'"
    fi
}

function checkFileLacks()
{
    # 1: description, 2: file, 3: pattern
    if [[ -f "$2" ]] && ! grep -qF -- "$3" "$2"
    then
        pass "$1"
    else
        fail "$1: '$3' unexpectedly found in '$2'"
    fi
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Stub PATHs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# macOS-like: "sed" is BSD sed, which does not understand "--version" and
# rejects "-i" without a backup suffix; GNU sed is available as "gsed"
MACOS_STUB="${WORK_DIR}/stubs-macos"
mkdir -p "${MACOS_STUB}"

cat > "${MACOS_STUB}/sed" <<EOF
#!/bin/bash
# BSD-like sed stub: no --version, and "-i" requires a backup suffix
for arg in "\$@"
do
    if [[ \${arg} == "--version" ]]
    then
        echo "sed: illegal option -- -" 1>&2
        exit 1
    fi
    if [[ \${arg} == "-i" ]]
    then
        echo "sed: -i may not be used without a backup suffix" 1>&2
        exit 1
    fi
done
exec "${REAL_SED}" "\$@"
EOF
chmod +x "${MACOS_STUB}/sed"

cat > "${MACOS_STUB}/gsed" <<EOF
#!/bin/bash
exec "${REAL_SED}" "\$@"
EOF
chmod +x "${MACOS_STUB}/gsed"

# No GNU sed at all: both "sed" and "gsed" are BSD-like
NOGNU_STUB="${WORK_DIR}/stubs-nognu"
mkdir -p "${NOGNU_STUB}"
\cp "${MACOS_STUB}/sed" "${NOGNU_STUB}/sed"
\cp "${MACOS_STUB}/sed" "${NOGNU_STUB}/gsed"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fixture case
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function makeCase()
{
    # 1: case directory to create
    local DIR="$1"

    mkdir -p "${DIR}"/0 "${DIR}"/system "${DIR}"/constant/polyMesh

    cat > "${DIR}"/0/D <<'EOF'
FoamFile { version 2.0; format ascii; class volVectorField; object D; }
boundaryField
{
    front { type symmetry; }
}
EOF

    cat > "${DIR}"/constant/polyMesh/boundary <<'EOF'
(
    front { type symmetry; }
)
EOF

    cat > "${DIR}"/system/blockMeshDict <<'EOF'
boundary
(
    symmetry front ((0 1 2 3));
);
EOF

    echo "solidModel linearGeometryTotalDisplacement;" \
        > "${DIR}"/constant/solidProperties

    cat > "${DIR}"/system/fvSchemes <<'EOF'
gradSchemes
{
    default extendedLeastSquares 0;
    grad(D) pointCellsLeastSquares;
}
EOF

    cat > "${DIR}"/system/mirrorMeshDict <<'EOF'
planeType pointAndNormal;
pointAndNormalDict
{
    point (0 0 0);
    normal (0 1 0);
}
EOF

    cat > "${DIR}"/force.gnuplot <<'EOF'
plot "postProcessing/sample/0/force.dat" u 1:($2) w l, "" u 1:($3) w l
EOF

    cat > "${DIR}"/plot.gnuplot <<'EOF'
plot "postProcessing/sample/0/data.xy", "postProcessing/sample.surfaces/x.raw"
EOF
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# solidModelDicts fixture, for runSolidModel
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function makeDicts()
{
    # 1: parent directory, 2: model name, 3: displacement name
    local PARENT="$1"
    local MODEL="$2"
    local DISP="$3"
    local DICTS="${PARENT}/solidModelDicts/${MODEL}"

    mkdir -p "${DICTS}"
    echo "solidModel ${MODEL};" > "${DICTS}"/solidProperties
    echo "solvers {}" > "${DICTS}"/fvSolution
    echo "gradSchemes {}" > "${DICTS}"/fvSchemes
    echo "${DISP}" > "${DICTS}"/displacementName
}

# A copy of solids4FoamScripts.sh whose DICTS_PARENT_DIR points at the fixture
function scriptsWithDicts()
{
    # 1: dicts parent dir, 2: output script path
    "${REAL_SED}" \
        "s|^    DICTS_PARENT_DIR=.*|    DICTS_PARENT_DIR=$1|" \
        "${SCRIPTS}" > "$2"
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Tests 1 to 4: resolution of the GNU sed command
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

echo "Test 1: GNU sed exposed as 'sed'"
RESULT=$(
    set -u
    unset SOLIDS4FOAM_SED
    # shellcheck disable=SC1090
    source "${SCRIPTS}"
    solids4Foam::requireGnuSed > /dev/null
    echo "${SOLIDS4FOAM_SED}"
)
check "SOLIDS4FOAM_SED resolves to sed" "sed" "${RESULT}"

echo "Test 2: GNU sed exposed only as 'gsed' (macOS-like)"
RESULT=$(
    set -u
    unset SOLIDS4FOAM_SED
    export PATH="${MACOS_STUB}:${PATH}"
    # shellcheck disable=SC1090
    source "${SCRIPTS}"
    solids4Foam::requireGnuSed > /dev/null
    echo "${SOLIDS4FOAM_SED}"
)
check "SOLIDS4FOAM_SED resolves to gsed" "gsed" "${RESULT}"

echo "Test 3: neither 'sed' nor 'gsed' provides GNU sed"
OUTPUT=$(
    set -u
    unset SOLIDS4FOAM_SED
    export PATH="${NOGNU_STUB}:${PATH}"
    # shellcheck disable=SC1090
    source "${SCRIPTS}"
    solids4Foam::requireGnuSed
    echo "SHOULD-NOT-REACH-HERE"
) && STATUS=0 || STATUS=$?
check "requireGnuSed exits with status 1" "1" "${STATUS}"
if echo "${OUTPUT}" | grep -q "requires GNU sed" &&
    echo "${OUTPUT}" | grep -q "brew install gnu-sed" &&
    ! echo "${OUTPUT}" | grep -q "SHOULD-NOT-REACH-HERE"
then
    pass "requireGnuSed prints a clear installation message"
else
    fail "requireGnuSed error message: ${OUTPUT}"
fi

echo "Test 4: SOLIDS4FOAM_SED is set on sourcing, so 'set -u' is safe"
RESULT=$(
    set -u
    unset SOLIDS4FOAM_SED
    # shellcheck disable=SC1090
    source "${SCRIPTS}"
    echo "unset-safe:[${SOLIDS4FOAM_SED}]"
) && STATUS=0 || STATUS=$?
check "sourcing under set -u succeeds" "0" "${STATUS}"
check "SOLIDS4FOAM_SED is empty but set" "unset-safe:[]" "${RESULT}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Tests 5 and 6: both conversion directions, GNU sed as 'sed' and as 'gsed'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function testConversion()
{
    # 1: label, 2: extra PATH prefix (may be empty), 3: expected sed command
    local LABEL="$1"
    local PATH_PREFIX="$2"
    local EXPECTED_SED="$3"

    local CASE_DIR="${WORK_DIR}/case-${LABEL}"
    local LOG="${WORK_DIR}/log-${LABEL}"
    local PRISTINE="${WORK_DIR}/pristine-${LABEL}"
    local STATUS

    makeCase "${CASE_DIR}"
    \cp -r "${CASE_DIR}" "${PRISTINE}"

    # Convert to the foam extend format and back, under 'set -u'
    (
        set -u
        unset SOLIDS4FOAM_SED
        if [[ -n ${PATH_PREFIX} ]]
        then
            export PATH="${PATH_PREFIX}:${PATH}"
        fi
        export WM_PROJECT="foam"
        export WM_PROJECT_VERSION="4.1"
        # shellcheck disable=SC1090
        source "${SCRIPTS}"

        solids4Foam::convertCaseFormat "${CASE_DIR}"

        echo "RESOLVED_SED=${SOLIDS4FOAM_SED}"

        # Check the forward conversion before restoring
        grep -q "symmetryPlane;" "${CASE_DIR}/0/D" ||
            { echo "FORWARD-FAIL: 0/D"; exit 1; }
        grep -q "symmetryPlane front" \
            "${CASE_DIR}/constant/polyMesh/blockMeshDict" ||
            { echo "FORWARD-FAIL: blockMeshDict"; exit 1; }
        grep -q " leastSquares;" "${CASE_DIR}/system/fvSchemes" ||
            { echo "FORWARD-FAIL: fvSchemes"; exit 1; }
        grep -q "basePoint" "${CASE_DIR}/system/mirrorMeshDict" ||
            { echo "FORWARD-FAIL: mirrorMeshDict basePoint"; exit 1; }
        grep -q "normalVector" "${CASE_DIR}/system/mirrorMeshDict" ||
            { echo "FORWARD-FAIL: mirrorMeshDict normalVector"; exit 1; }
        grep -q "forces.dat" "${CASE_DIR}/force.gnuplot" ||
            { echo "FORWARD-FAIL: force.gnuplot"; exit 1; }
        grep -q 'postProcessing/sets/' "${CASE_DIR}/plot.gnuplot" ||
            { echo "FORWARD-FAIL: plot.gnuplot"; exit 1; }
        echo "FORWARD-OK"

        solids4Foam::restoreCaseFormat "${CASE_DIR}"
        echo "RESTORE-DONE"
    ) > "${LOG}" 2>&1 && STATUS=0 || STATUS=$?

    check "${LABEL}: conversion round trip exits cleanly" "0" "${STATUS}"
    check "${LABEL}: uses '${EXPECTED_SED}'" \
        "RESOLVED_SED=${EXPECTED_SED}" \
        "$(grep '^RESOLVED_SED=' "${LOG}" || true)"
    check "${LABEL}: forward conversion (com -> foam extend)" \
        "FORWARD-OK" "$(grep '^FORWARD-OK' "${LOG}" || true)"

    # Reverse conversion: the case should be back in the stored format
    checkFileContains "${LABEL}: reverse conversion restores symmetry in 0/D" \
        "${CASE_DIR}/0/D" "symmetry;"
    checkFileLacks "${LABEL}: no symmetryPlane remains in 0/D" \
        "${CASE_DIR}/0/D" "symmetryPlane"
    checkFileContains "${LABEL}: reverse conversion restores fvSchemes" \
        "${CASE_DIR}/system/fvSchemes" "pointCellsLeastSquares;"
    checkFileContains "${LABEL}: reverse conversion restores mirrorMeshDict" \
        "${CASE_DIR}/system/mirrorMeshDict" "    point ("
    checkFileLacks "${LABEL}: no basePoint remains in mirrorMeshDict" \
        "${CASE_DIR}/system/mirrorMeshDict" "basePoint"
    checkFileContains "${LABEL}: reverse conversion restores force.gnuplot" \
        "${CASE_DIR}/force.gnuplot" "force.dat"
    checkFileContains "${LABEL}: reverse conversion restores plot.gnuplot" \
        "${CASE_DIR}/plot.gnuplot" "postProcessing/sample/"

    if diff -r "${PRISTINE}/0" "${CASE_DIR}/0" > /dev/null 2>&1 &&
        diff "${PRISTINE}/system/fvSchemes" "${CASE_DIR}/system/fvSchemes" \
            > /dev/null 2>&1 &&
        diff "${PRISTINE}/system/mirrorMeshDict" \
            "${CASE_DIR}/system/mirrorMeshDict" > /dev/null 2>&1 &&
        diff "${PRISTINE}/force.gnuplot" "${CASE_DIR}/force.gnuplot" \
            > /dev/null 2>&1 &&
        diff "${PRISTINE}/plot.gnuplot" "${CASE_DIR}/plot.gnuplot" \
            > /dev/null 2>&1
    then
        pass "${LABEL}: round trip is lossless for the edited files"
    else
        fail "${LABEL}: round trip changed the edited files"
    fi
}

echo "Test 5: case-format conversion, GNU sed as 'sed'"
testConversion "sed" "" "sed"

echo "Test 6: case-format conversion, GNU sed as 'gsed' (macOS-like)"
testConversion "gsed" "${MACOS_STUB}" "gsed"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Tests 7 to 10: the D/DD rename paths in runSolidModel
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function testRunSolidModel()
{
    # 1: label, 2: PATH prefix, 3: expected sed, 4: field present, 5: wanted
    local LABEL="$1"
    local PATH_PREFIX="$2"
    local EXPECTED_SED="$3"
    local HAVE="$4"
    local WANT="$5"

    local BASE="${WORK_DIR}/rsm-${LABEL}"
    local CASE_DIR="${BASE}/case"
    local LOG="${WORK_DIR}/log-rsm-${LABEL}"
    local STATUS

    mkdir -p "${CASE_DIR}/0" "${CASE_DIR}/constant" "${CASE_DIR}/system"
    makeDicts "${BASE}" "testModel" "${WANT}"
    scriptsWithDicts "${BASE}" "${BASE}/solids4FoamScripts.sh"

    cat > "${CASE_DIR}/0/${HAVE}" <<EOF
FoamFile
{
    class volVectorField;
    object ${HAVE};
}
internalField uniform (0 0 0);
EOF

    # Stub Allrun, so that no OpenFOAM solver is run
    echo '#!/bin/bash' > "${CASE_DIR}/Allrun"
    echo 'echo "stub Allrun"' >> "${CASE_DIR}/Allrun"
    chmod +x "${CASE_DIR}/Allrun"

    (
        set -u
        unset SOLIDS4FOAM_SED
        if [[ -n ${PATH_PREFIX} ]]
        then
            export PATH="${PATH_PREFIX}:${PATH}"
        fi
        cd "${CASE_DIR}" || exit 1
        # shellcheck disable=SC1090
        source "${BASE}/solids4FoamScripts.sh"
        solids4Foam::runSolidModel "${CASE_DIR}" testModel
        echo "RESOLVED_SED=${SOLIDS4FOAM_SED}"
    ) > "${LOG}" 2>&1 && STATUS=0 || STATUS=$?

    check "${LABEL}: runSolidModel exits cleanly under set -u" "0" "${STATUS}"
    check "${LABEL}: uses '${EXPECTED_SED}'" \
        "RESOLVED_SED=${EXPECTED_SED}" \
        "$(grep '^RESOLVED_SED=' "${LOG}" || true)"

    if [[ -f "${CASE_DIR}/0/${WANT}" && ! -f "${CASE_DIR}/0/${HAVE}" ]]
    then
        pass "${LABEL}: 0/${HAVE} renamed to 0/${WANT}"
    else
        fail "${LABEL}: 0/${HAVE} was not renamed to 0/${WANT}"
    fi

    checkFileContains "${LABEL}: object entry updated to ${WANT}" \
        "${CASE_DIR}/0/${WANT}" "object ${WANT};"
}

echo "Test 7: runSolidModel DD -> D rename, GNU sed as 'sed'"
testRunSolidModel "DD-to-D-sed" "" "sed" "DD" "D"

echo "Test 8: runSolidModel D -> DD rename, GNU sed as 'sed'"
testRunSolidModel "D-to-DD-sed" "" "sed" "D" "DD"

echo "Test 9: runSolidModel DD -> D rename, GNU sed as 'gsed' (macOS-like)"
testRunSolidModel "DD-to-D-gsed" "${MACOS_STUB}" "gsed" "DD" "D"

echo "Test 10: runSolidModel D -> DD rename, GNU sed as 'gsed' (macOS-like)"
testRunSolidModel "D-to-DD-gsed" "${MACOS_STUB}" "gsed" "D" "DD"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "| ${N_PASS} passed, ${N_FAIL} failed"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

if [[ ${N_FAIL} -ne 0 ]]
then
    exit 1
fi
