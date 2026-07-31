#!/bin/bash
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Library
#     solids4Foam bash functions for converting a case to the appropriate
#     OpenFOAM format
#
#     The script broadly follows the style given at
#     https://google.github.io/styleguide/shellguide.html as well as the OpenFOAM
#     coding style at https://openfoam.org/dev/coding-style-guide
#     The script is checked with https://www.shellcheck.net
#
# Authors
#     Philip Cardiff, UCD
#
# License
#     GNU Lesser General Public License, version 3.
#     https://www.gnu.org/licenses/lgpl-3.0.en.html
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# On macOS, ensure OpenFOAM libraries remain discoverable for child processes.
case "$(uname -s)" in
Darwin)
    export DYLD_LIBRARY_PATH="${FOAM_LIBBIN}:${FOAM_USER_LIBBIN}${DYLD_LIBRARY_PATH:+:${DYLD_LIBRARY_PATH}}"
    ;;
esac

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Case format
#
#     Tutorial cases are stored in the OpenFOAM.com (ESI) format. At run time,
#     solids4Foam::convertCaseFormat converts a case to the format required by
#     the loaded OpenFOAM flavour; solids4Foam::restoreCaseFormat converts it
#     back to the stored format.
#
#     Files that cannot be converted with a text substitution are stored twice,
#     with the unsuffixed name always holding the OpenFOAM.com version, e.g.
#     system/functions and system/functions.foamextend. During conversion the
#     unsuffixed file is backed up with a .openfoam suffix and the flavour
#     specific file is copied into its place; restoring moves the backup back.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# foamFlavour
#     Echoes the loaded OpenFOAM flavour: foamextend, com or org
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::foamFlavour()
{
    if [[ ${WM_PROJECT} == "foam" ]]
    then
        echo "foamextend"
    elif [[ ${WM_PROJECT_VERSION} == *"v"* ]]
    then
        echo "com"
    else
        echo "org"
    fi
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# blockMeshDictDir
#     Echoes the directory blockMesh reads blockMeshDict from, for the loaded
#     OpenFOAM flavour. Intended for cases that generate their blockMeshDict,
#     e.g. with m4.
# Arguments:
#     1: optional region name
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::blockMeshDictDir()
{
    local REGION="${1:-}"

    if [[ $(solids4Foam::foamFlavour) == "foamextend" ]]
    then
        if [[ -n ${REGION} ]]
        then
            echo "constant/${REGION}/polyMesh"
        else
            echo "constant/polyMesh"
        fi
    else
        if [[ -n ${REGION} ]]
        then
            echo "system/${REGION}"
        else
            echo "system"
        fi
    fi
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# requireGnuSed
#     Exits if GNU sed is not available, as the conversion functions rely on it
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::requireGnuSed()
{
    if sed --version 2>/dev/null | grep -q "GNU sed"
    then
        echo "GNU sed detected"
    else
        echo "Error: This script requires GNU sed. Please install it (e.g."
        echo "via Homebrew: 'brew install gnu-sed') and use 'gsed' instead."
        exit 1
    fi
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# fieldDirs
#     Echoes the initial field directories present in a case, e.g. 0 and 0.orig
# Arguments:
#     1: CASE_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::fieldDirs()
{
    local CASE_DIR="$1"
    local DIR

    for DIR in "${CASE_DIR}"/0 "${CASE_DIR}"/0.orig "${CASE_DIR}"/0.org
    do
        if [[ -d ${DIR} ]]
        then
            echo "${DIR}"
        fi
    done
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# convertCaseFormat
#     Converts a case from the stored OpenFOAM.com format to the format required
#     by the loaded OpenFOAM flavour. No changes are applied when OpenFOAM.com
#     is loaded.
# Arguments:
#     1: CASE_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::convertCaseFormat()
{
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4Foam::convertCaseFormat start                                |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo

    # Check number of input parameters is correct
    if [[ $# -ne 1 ]]
    then
        solids4Foam::err "convertCaseFormat: incorrect number of parameters"
    fi

    # Give sensible names to the argument
    local CASE_DIR="$1"

    case "$(solids4Foam::foamFlavour)" in
    com)
        echo "OpenFOAM.com loaded: cases are stored in this format"
        echo "No changes made"; echo
        ;;
    org)
        solids4Foam::requireGnuSed
        solids4Foam::applyOpenFOAMOrgTweaks "${CASE_DIR}"
        ;;
    foamextend)
        solids4Foam::requireGnuSed
        solids4Foam::applyFoamExtendTweaks "${CASE_DIR}"
        ;;
    esac

    echo
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4Foam::convertCaseFormat end                                  |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# restoreCaseFormat
#     Restores a case to the stored OpenFOAM.com format, regardless of which
#     OpenFOAM flavour is loaded. The operation is idempotent.
# Arguments:
#     1: CASE_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::restoreCaseFormat()
{
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4Foam::restoreCaseFormat start                                |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo

    # Check number of input parameters is correct
    if [[ $# -ne 1 ]]
    then
        solids4Foam::err "restoreCaseFormat: incorrect number of parameters"
    fi

    # Give sensible names to the argument
    local CASE_DIR="$1"

    solids4Foam::requireGnuSed

    solids4Foam::undoFoamExtendTweaks "${CASE_DIR}"
    solids4Foam::undoOpenFOAMOrgTweaks "${CASE_DIR}"

    echo
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4Foam::restoreCaseFormat end                                  |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# convertCaseFormatFoamExtend
#     Deprecated: retained so that existing user case scripts keep working.
#     Use solids4Foam::restoreCaseFormat instead.
# Arguments:
#     1: CASE_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::convertCaseFormatFoamExtend()
{
    echo "Note: solids4Foam::convertCaseFormatFoamExtend is deprecated as cases"
    echo "      are now stored in the OpenFOAM.com format; please use"
    echo "      solids4Foam::restoreCaseFormat instead"
    echo

    solids4Foam::restoreCaseFormat "$@"
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# applyFoamExtendTweaks
#     Converts a case from the stored OpenFOAM.com format to the foam extend
#     format
# Arguments:
#     1: CASE_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::applyFoamExtendTweaks()
{
    local CASE_DIR="$1"
    local FILE
    local DIR

    # 1. symmetry in OpenFOAM becomes symmetryPlane in foam extend
    #    Patches are given as "symmetry;" in field files, boundary files and
    #    createPatchDict, and as "symmetry <patchName>" in blockMeshDict

    while IFS= read -r -d '' FILE
    do
        if grep -qE "symmetry;|^[[:space:]]*symmetry[[:space:]]" "${FILE}"
        then
            echo "Changing symmetry to symmetryPlane in ${FILE}"
            sed -i 's|symmetry;|symmetryPlane;|g' "${FILE}"
            sed -i -E 's|^([[:space:]]*)symmetry([[:space:]]+)|\1symmetryPlane\2|' \
                "${FILE}"
        fi
    done < <(find "${CASE_DIR}" \( -name 'blockMeshDict*' -o -name boundary \
        -o -name createPatchDict \) -print0)

    while IFS= read -r DIR
    do
        while IFS= read -r -d '' FILE
        do
            if grep -q "symmetry;" "${FILE}"
            then
                echo "Changing symmetry to symmetryPlane in ${FILE}"
                sed -i 's|symmetry;|symmetryPlane;|g' "${FILE}"
            fi
        done < <(find "${DIR}" -type f -print0)
    done < <(solids4Foam::fieldDirs "${CASE_DIR}")

    # 2. foam extend reads the blockMeshDict from constant/polyMesh
    for DIR in "" solid fluid
    do
        if [[ -z ${DIR} ]]
        then
            local SRC="${CASE_DIR}/system/blockMeshDict"
            local DST="${CASE_DIR}/constant/polyMesh"
        else
            local SRC="${CASE_DIR}/system/${DIR}/blockMeshDict"
            local DST="${CASE_DIR}/constant/${DIR}/polyMesh"
        fi

        # Note: -f is false for a dangling symlink, so -L is also checked, as
        # a case may link its blockMeshDict to a variant file
        if [[ -f ${SRC} || -L ${SRC} ]]
        then
            echo "Moving ${SRC} to ${DST}"
            mkdir -p "${DST}"
            \mv "${SRC}" "${DST}"
        fi
    done

    # 2b. Use the foam extend functions file, if present
    if [[ -f "${CASE_DIR}"/system/functions.foamextend ]]
    then
        echo "Replacing system/functions with system/functions.foamextend"
        \cp "${CASE_DIR}"/system/functions \
            "${CASE_DIR}"/system/functions.openfoam
        \cp -f "${CASE_DIR}"/system/functions.foamextend \
            "${CASE_DIR}"/system/functions
    fi

    # 3. Rename turbulence model
    if [[ -n $(find "${CASE_DIR}" -name turbulenceProperties) ]]
    then
        echo "Changing RAS to RASModel in turbulenceProperties"
        find "${CASE_DIR}" -name turbulenceProperties \
            -exec sed -i "s/RAS;/RASModel;/g" {} +
    fi

    # 4. Use the foam extend boundaryData, if present
    if [[ -d "${CASE_DIR}"/constant/boundaryData.foamextend ]]
    then
        echo "Moving constant/boundaryData to constant/boundaryData.openfoam"
        \mv "${CASE_DIR}"/constant/boundaryData \
            "${CASE_DIR}"/constant/boundaryData.openfoam

        echo "Moving constant/boundaryData.foamextend to constant/boundaryData"
        \mv "${CASE_DIR}"/constant/boundaryData.foamextend \
            "${CASE_DIR}"/constant/boundaryData
    fi

    # 6. foam extend uses timeVaryingUniformFixedValue instead of
    #    uniformFixedValue
    if [[ -n $(find "${CASE_DIR}" -name p) ]]
    then
        echo "Changing uniformFixedValue to timeVaryingUniformFixedValue in p"
        find "${CASE_DIR}" -name p \
            -exec sed -i "s|^\([[:space:]]*\)type\(.*\)uniformFixedValue;|\1//type\2uniformFixedValue;|g" {} +
        find "${CASE_DIR}" -name p \
            -exec sed -i "s|^\([[:space:]]*\)//type\(.*\)timeVaryingUniformFixedValue;|\1type\2timeVaryingUniformFixedValue;|g" {} +
    fi

    # 7. Use the foam extend changeDictionaryDict, if present
    if [[ -f "${CASE_DIR}/system/changeDictionaryDict.foamextend" ]]
    then
        echo "Moving system/changeDictionaryDict to system/changeDictionaryDict.openfoam"
        \mv "${CASE_DIR}/system/changeDictionaryDict" \
            "${CASE_DIR}/system/changeDictionaryDict.openfoam"
        echo "Moving system/changeDictionaryDict.foamextend to system/changeDictionaryDict"
        \mv "${CASE_DIR}/system/changeDictionaryDict.foamextend" \
            "${CASE_DIR}/system/changeDictionaryDict"
    fi

    # 8. Use the foam extend createPatchDict, if present
    if [[ -f "${CASE_DIR}/system/createPatchDict.foamextend" ]]
    then
        echo "Moving system/createPatchDict to system/createPatchDict.openfoam"
        \mv "${CASE_DIR}/system/createPatchDict" \
            "${CASE_DIR}/system/createPatchDict.openfoam"
        echo "Moving system/createPatchDict.foamextend to system/createPatchDict"
        \mv "${CASE_DIR}/system/createPatchDict.foamextend" \
            "${CASE_DIR}/system/createPatchDict"
    fi

    # 9. pointCellsLeastSquares is used as the gradScheme for the solid in
    #    OpenFOAM, as it is consistent with boundary non-orthogonal correction;
    #    foam extend uses leastSquares
    if [[ -f "${CASE_DIR}"/constant/solidProperties ]]
    then
        echo "foam extend specific: replacing 'pointCellsLeastSquares' with"
        echo "'leastSquares' in system/fvSchemes"
        sed -i "s/ pointCellsLeastSquares;/ leastSquares;/g" \
            "${CASE_DIR}"/system/fvSchemes
    elif [[ -f "${CASE_DIR}"/constant/solid/solidProperties ]]
    then
        echo "foam extend specific: replacing 'pointCellsLeastSquares' with"
        echo "'leastSquares' in system/solid/fvSchemes"
        sed -i "s/ pointCellsLeastSquares;/ leastSquares;/g" \
            "${CASE_DIR}"/system/solid/fvSchemes
    fi

    # 10. foam extend writes forces.dat rather than force.dat
    solids4Foam::useForcesDat "${CASE_DIR}"

    # 11. foam extend samples to postProcessing/sets
    # 12. foam extend writes surfaces to postProcessing/surfaces
    for FILE in $(find "${CASE_DIR}" -name plot.gnuplot)
    do
        echo "Updating ${FILE}"
        sed -i "s@postProcessing/sample/@postProcessing/sets/@g" "${FILE}"
        sed -i "s@postProcessing/sample.surfaces/@postProcessing/surfaces/@g" "${FILE}"
    done

    # 13. foam extend uses basePoint and normalVector in mirrorMeshDict
    for FILE in "${CASE_DIR}"/system/mirrorMeshDict \
        "${CASE_DIR}"/system/solid/mirrorMeshDict
    do
        if [[ -f ${FILE} ]]
        then
            echo "foam extend specific: replacing 'point' and 'normal' with"
            echo "'basePoint' and 'normalVector' in ${FILE}"
            sed -i -E 's/^([[:space:]]*)point([[:space:]])/\1basePoint\2/g' "${FILE}"
            sed -i -E 's/^([[:space:]]*)normal([[:space:]])/\1normalVector\2/g' "${FILE}"
        fi
    done
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# undoFoamExtendTweaks
#     Reverses applyFoamExtendTweaks, i.e. converts a case to the stored
#     OpenFOAM.com format. The operation is idempotent and is applied regardless
#     of which OpenFOAM flavour is loaded.
# Arguments:
#     1: CASE_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::undoFoamExtendTweaks()
{
    local CASE_DIR="$1"
    local FILE
    local DIR

    # 1. symmetryPlane in foam extend becomes symmetry in OpenFOAM
    while IFS= read -r -d '' FILE
    do
        if grep -q "symmetryPlane" "${FILE}"
        then
            echo "Changing symmetryPlane to symmetry in ${FILE}"
            sed -i 's|symmetryPlane|symmetry|g' "${FILE}"
        fi
    done < <(find "${CASE_DIR}" \( -name 'blockMeshDict*' -o -name boundary \
        -o -name createPatchDict \) -print0)

    while IFS= read -r DIR
    do
        while IFS= read -r -d '' FILE
        do
            if grep -q "symmetryPlane;" "${FILE}"
            then
                echo "Changing symmetryPlane to symmetry in ${FILE}"
                sed -i 's|symmetryPlane;|symmetry;|g' "${FILE}"
            fi
        done < <(find "${DIR}" -type f -print0)
    done < <(solids4Foam::fieldDirs "${CASE_DIR}")

    # 2. OpenFOAM reads the blockMeshDict from system
    for DIR in "" solid fluid
    do
        if [[ -z ${DIR} ]]
        then
            local SRC="${CASE_DIR}/constant/polyMesh/blockMeshDict"
            local DST="${CASE_DIR}/system"
        else
            local SRC="${CASE_DIR}/constant/${DIR}/polyMesh/blockMeshDict"
            local DST="${CASE_DIR}/system/${DIR}"
        fi

        # Note: -f is false for a dangling symlink, so -L is also checked, as
        # a case may link its blockMeshDict to a variant file
        if [[ -f ${SRC} || -L ${SRC} ]]
        then
            mkdir -p "${DST}"

            if [[ -f "${DST}/blockMeshDict" || -L "${DST}/blockMeshDict" ]]
            then
                echo "Removing generated ${SRC}; ${DST}/blockMeshDict exists"
                \rm "${SRC}"
            else
                echo "Moving ${SRC} to ${DST}"
                \mv "${SRC}" "${DST}"
            fi

            # Remove the polyMesh directory if the conversion created it
            rmdir "$(dirname "${SRC}")" 2>/dev/null || true
        fi
    done

    # 2b. Restore the OpenFOAM functions file, if it was replaced
    if [[ -f "${CASE_DIR}"/system/functions.openfoam ]]
    then
        echo "Restoring system/functions from system/functions.openfoam"
        \mv -f "${CASE_DIR}"/system/functions.openfoam \
            "${CASE_DIR}"/system/functions
    fi

    # 3. Rename turbulence model
    if [[ -n $(find "${CASE_DIR}" -name turbulenceProperties) ]]
    then
        echo "Changing RASModel to RAS in turbulenceProperties"
        find "${CASE_DIR}" -name turbulenceProperties \
            -exec sed -i "s/RASModel;/RAS;/g" {} +
    fi

    # 4. Restore the OpenFOAM boundaryData, if it was replaced
    if [[ -d "${CASE_DIR}"/constant/boundaryData.openfoam ]]
    then
        echo "Moving constant/boundaryData to constant/boundaryData.foamextend"
        \rm -rf "${CASE_DIR}"/constant/boundaryData.foamextend
        \mv "${CASE_DIR}"/constant/boundaryData \
            "${CASE_DIR}"/constant/boundaryData.foamextend

        echo "Moving constant/boundaryData.openfoam to constant/boundaryData"
        \mv "${CASE_DIR}"/constant/boundaryData.openfoam \
            "${CASE_DIR}"/constant/boundaryData
    fi

    # 6. OpenFOAM uses uniformFixedValue instead of timeVaryingUniformFixedValue
    if [[ -n $(find "${CASE_DIR}" -name p) ]]
    then
        echo "Changing timeVaryingUniformFixedValue to uniformFixedValue in p"
        find "${CASE_DIR}" -name p \
            -exec sed -i "s|^\([[:space:]]*\)//type\(.*\)uniformFixedValue;|\1type\2uniformFixedValue;|g" {} +
        find "${CASE_DIR}" -name p \
            -exec sed -i "s|^\([[:space:]]*\)type\(.*\)timeVaryingUniformFixedValue;|\1//type\2timeVaryingUniformFixedValue;|g" {} +
    fi

    # 7. Restore the OpenFOAM changeDictionaryDict, if it was replaced
    if [[ -f "${CASE_DIR}/system/changeDictionaryDict.openfoam" ]]
    then
        echo "Moving system/changeDictionaryDict to system/changeDictionaryDict.foamextend"
        \mv "${CASE_DIR}/system/changeDictionaryDict" \
            "${CASE_DIR}/system/changeDictionaryDict.foamextend"
        echo "Moving system/changeDictionaryDict.openfoam to system/changeDictionaryDict"
        \mv "${CASE_DIR}/system/changeDictionaryDict.openfoam" \
            "${CASE_DIR}/system/changeDictionaryDict"
    fi

    # 8. Restore the OpenFOAM createPatchDict, if it was replaced
    if [[ -f "${CASE_DIR}/system/createPatchDict.openfoam" ]]
    then
        echo "Moving system/createPatchDict to system/createPatchDict.foamextend"
        \mv "${CASE_DIR}/system/createPatchDict" \
            "${CASE_DIR}/system/createPatchDict.foamextend"
        echo "Moving system/createPatchDict.openfoam to system/createPatchDict"
        \mv "${CASE_DIR}/system/createPatchDict.openfoam" \
            "${CASE_DIR}/system/createPatchDict"
    fi

    # 9. Either pointCellsLeastSquares or edgeCellsLeastSquares should be used
    #    as the gradScheme for the solid in OpenFOAM, as these are the only
    #    schemes consistent with boundary non-orthogonal correction
    if [[ -f "${CASE_DIR}"/constant/solidProperties ]]
    then
        echo "OpenFOAM specific: replacing 'leastSquares' with"
        echo "'pointCellsLeastSquares' in system/fvSchemes"
        sed -i "s/ leastSquares;/ pointCellsLeastSquares;/g" \
            "${CASE_DIR}"/system/fvSchemes
    elif [[ -f "${CASE_DIR}"/constant/solid/solidProperties ]]
    then
        echo "OpenFOAM specific: replacing 'leastSquares' with"
        echo "'pointCellsLeastSquares' in system/solid/fvSchemes"
        sed -i "s/ leastSquares;/ pointCellsLeastSquares;/g" \
            "${CASE_DIR}"/system/solid/fvSchemes
    fi

    # 10. OpenFOAM.com writes force.dat rather than forces.dat
    solids4Foam::useForceDat "${CASE_DIR}"

    # 11. OpenFOAM samples to postProcessing/sample
    # 12. OpenFOAM writes surfaces to postProcessing/sample.surfaces
    for FILE in $(find "${CASE_DIR}" -name plot.gnuplot)
    do
        echo "Updating ${FILE}"
        sed -i "s@postProcessing/sets/@postProcessing/sample/@g" "${FILE}"
        sed -i "s@postProcessing/surfaces/@postProcessing/sample.surfaces/@g" "${FILE}"
    done

    # 13. OpenFOAM uses point and normal in mirrorMeshDict
    for FILE in "${CASE_DIR}"/system/mirrorMeshDict \
        "${CASE_DIR}"/system/solid/mirrorMeshDict
    do
        if [[ -f ${FILE} ]]
        then
            echo "OpenFOAM specific: replacing 'basePoint' and 'normalVector'"
            echo "with 'point' and 'normal' in ${FILE}"
            sed -i -E 's/\bbasePoint\b/point/g' "${FILE}"
            sed -i -E 's/\bnormalVector\b/normal/g' "${FILE}"
        fi
    done
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# applyOpenFOAMOrgTweaks
#     Converts a case from the stored OpenFOAM.com format to the OpenFOAM.org
#     format
# Arguments:
#     1: CASE_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::applyOpenFOAMOrgTweaks()
{
    local CASE_DIR="$1"

    # 5. OpenFOAM.org renamed the uniform and face sampledSet types
    if [[ -f "${CASE_DIR}"/system/sample ]]
    then
        echo "OpenFOAM.org specific: replacing 'uniform' with 'lineUniform' in"
        echo "system/sample"
        sed -i "s/type.*uniform;/type lineUniform;/g" "${CASE_DIR}"/system/sample

        echo "OpenFOAM.org specific: replacing 'face' with 'lineFace' in"
        echo "system/sample"
        sed -i "s/type.*face;/type lineFace;/g" "${CASE_DIR}"/system/sample
    fi

    # 10. OpenFOAM.org writes forces.dat rather than force.dat
    solids4Foam::useForcesDat "${CASE_DIR}"
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# undoOpenFOAMOrgTweaks
#     Reverses applyOpenFOAMOrgTweaks. The operation is idempotent and is
#     applied regardless of which OpenFOAM flavour is loaded.
# Arguments:
#     1: CASE_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::undoOpenFOAMOrgTweaks()
{
    local CASE_DIR="$1"

    # 5. OpenFOAM.org renamed the uniform and face sampledSet types
    if [[ -f "${CASE_DIR}"/system/sample ]]
    then
        echo "Replacing 'lineUniform' with 'uniform' in system/sample"
        sed -i "s/type.*lineUniform;/type uniform;/g" "${CASE_DIR}"/system/sample

        echo "Replacing 'lineFace' with 'face' in system/sample"
        sed -i "s/type.*lineFace;/type face;/g" "${CASE_DIR}"/system/sample
    fi

    # 10. OpenFOAM.com writes force.dat rather than forces.dat
    solids4Foam::useForceDat "${CASE_DIR}"
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# useForcesDat, useForceDat
#     Set the forces functionObject output expected by the gnuplot scripts.
#
#     OpenFOAM.com writes force.dat, with the total force in columns 2 to 4,
#     followed by the pressure and viscous contributions, so the total force is
#     plotted directly as ($2) and ($3).
#
#     foam extend and OpenFOAM.org write forces.dat, with only the pressure and
#     viscous contributions, followed by the moments, so the total force is
#     plotted as their sum, ($2+$5) and ($3+$6).
# Arguments:
#     1: CASE_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::useForcesDat()
{
    local CASE_DIR="$1"
    local FILE

    for FILE in $(find "${CASE_DIR}" -name force.gnuplot)
    do
        echo "Changing force.dat to forces.dat in ${FILE}"
        sed -i "s|force\.dat|forces.dat|g" "${FILE}"

        echo "Summing the pressure and viscous forces in ${FILE}"
        sed -i 's|(\$2)|($2+$5)|g' "${FILE}"
        sed -i 's|(\$3)|($3+$6)|g' "${FILE}"
    done
}

function solids4Foam::useForceDat()
{
    local CASE_DIR="$1"
    local FILE

    for FILE in $(find "${CASE_DIR}" -name force.gnuplot)
    do
        echo "Changing forces.dat to force.dat in ${FILE}"
        sed -i "s|forces\.dat|force.dat|g" "${FILE}"

        echo "Using the total force in ${FILE}"
        sed -i 's|(\$2+\$5)|($2)|g' "${FILE}"
        sed -i 's|(\$3+\$6)|($3)|g' "${FILE}"
    done
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Print error message to stderr
# Arguments:
#     1. error message
#     2. optional: log file that will be copied to errorCommandLog.txt
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::err()
{
    echo; echo "ERROR: see error.txt"

    # Error message
    errMsg="[$(date +'%Y-%m-%dT%H:%M:%S%z')]: $*"

    # Print to stderr
    echo "${errMsg}" >&2

    # Print error to error.txt file
    echo "${errMsg}" > error.txt

    # Copy log file to errorCommandLog.txt file
    if [[ $# -gt 1 ]]
    then
        \cp -f "${2}" errorCommandLog.txt
        echo "       see errorCommandLog.txt"
    fi

    echo

    # Stop with error
    exit 1
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# caseOnlyRunsWithFoamExtend
#     Exit if OpenFOAM version is not foam-extend
# Arguments:
#     None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::caseOnlyRunsWithFoamExtend()
{
    if [[ $WM_PROJECT != "foam" ]]
    then
        echo; echo "This case currently only runs in foam-extend"; echo
        exit 0
    fi
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# caseOnlyRunsWithOpenFOAM
#     Exit if OpenFOAM version is foam-extend
# Arguments:
#     None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::caseDoesNotRunWithFoamExtend()
{
    if [[ $WM_PROJECT == "foam" ]]
    then
        echo; echo "This case currently does not run with foam-extend"; echo
        exit 0
    fi
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# caseDoesNotRunWithOpenFOAMOrg
#     Exit if OpenFOAM version is OpenFOAM.org
# Arguments:
#     None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::caseDoesNotRunWithOpenFOAMOrg()
{
    if [[ $WM_PROJECT == "OpenFOAM" ]]
    then
        if [[ $WM_PROJECT_VERSION != v* ]]
        then
            echo; echo "This case currently does not run with OpenFOAM.org"; echo
            exit 0
        fi
    fi
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# requirePetscOrExitSilently
#     Exit cleanly (exit 0) if PETSc is not available, i.e. PETSC_DIR is unset.
#     Used by tutorials whose (selected) solution approach requires PETSc. This
#     function is intentionally unaware of case-specific approach names: the
#     caller decides whether the current approach needs PETSc before calling it.
# Arguments:
#     None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::requirePetscOrExitSilently()
{
    if [[ -z "${PETSC_DIR:-}" ]]
    then
        echo
        echo "Skipping this case as PETSc is not installed"
        echo "If you would like to run this case, solids4foam should be rebuilt with PETSc"
        echo "(set the PETSC_DIR variable and rebuild solids4foam)"
        echo
        exit 0
    fi
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# regressionCaseSkipped
#     Return success when an Allrun log contains a known case-skip message.
#     Regression scripts use this to exit cleanly when a tutorial is not
#     intended to run in the current OpenFOAM flavour.
# Arguments:
#     1: LOG_FILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::regressionCaseSkipped()
{
    local LOG_FILE=$1

    if [[ ! -f "${LOG_FILE}" ]]
    then
        return 1
    fi

    if grep -Eq "This case currently only runs in foam-extend|This case currently does not run with foam-extend|This case currently does not run with OpenFOAM.org|Skipping this case as it does not currently (working|run) with OpenFOAM.org|OpenFOAM-v[0-9]+ or a newer version is required|Skipping this case as PETSc is not installed" "${LOG_FILE}"
    then
        return 0
    fi

    return 1
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# removeEmptyDirs
#     Ported from preCICE toolbox
#     Remove empty time directories that are generated when running FSI cases
#     with preCICE
# Arguments:
#     None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::removeEmptyDirs()
{
    (
        set -e -u
        echo "Removing time directories without results"

        for f in [0-9]* [0-9]*.[0-9]*; do
            if ! [ -f "${f}/U" ] && ! [ -f "${f}/T" ] && ! [ -f "${f}/U.gz" ] && ! [ -f "${f}/T.gz" ] && ! [ -f "${f}/D" ] && ! [ -f "${f}/pointD" ] && ! [ -f "${f}/DD" ] && ! [ -f "${f}/pointDD" ] && ! [ -f "${f}/D.gz" ] && ! [ -f "${f}/pointD.gz" ] && ! [ -f "${f}/DD.gz" ] && ! [ -f "${f}/pointDD.gz" ]; then
                rm -rfv "${f}"
            fi
        done
        if [ -d processor0 ]; then
            for d in processor*; do
                cd "${d}"
                for f in [0-9]* [0-9]*.[0-9]*; do
                    if ! [ -f "${f}/U" ] && ! [ -f "${f}/T" ] && ! [ -f "${f}/U.gz" ] && ! [ -f "${f}/T.gz" ] && ! [ -f "${f}/D" ] && ! [ -f "${f}/pointD" ] && ! [ -f "${f}/DD" ] && ! [ -f "${f}/pointDD" ] && ! [ -f "${f}/D.gz" ] && ! [ -f "${f}/pointD.gz" ] && ! [ -f "${f}/DD.gz" ] && ! [ -f "${f}/pointDD.gz" ]; then
                        rm -rfv "${f}"
                    fi
                done
                cd ..
            done
        fi
    )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# runApplication
#     runApplication that works that same irrespective of the OpenFOAM version
#     Copied from OpenFOAM-v2012
# Arguments:
#     None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::runApplication()
{
    local appName="" appRun="" logFile="" logMode="" appArgs=""

    # Parse options until executable is encountered
    while [ "$#" -gt 0 ] && [ -z "$appRun" ]
    do
        case "$1" in
            -a | -append)
                logMode=append
                ;;
            -o | -overwrite)
                logMode=overwrite
                ;;
            -s | -suffix)
                logFile=".$2"
                shift
                ;;
            -decomposeParDict)
                appArgs="$appArgs $1 $2"
                shift
                ;;
            '')
                ;;
            *)
                appRun="$1"
                ;;
        esac
        shift
    done

    appName="${appRun##*/}"
    logFile="log.$appName$logFile"

    if [ -f "$logFile" ] && [ -z "$logMode" ]
    then
        echo "$appName already run on $PWD:" \
             "remove log file '$logFile' to re-run"
    else
        echo "Running $appRun on $PWD"
        if [ "$logMode" = append ]
        then
            $appRun $appArgs "$@" >> $logFile 2>&1 || echo "ERROR" >> $logFile
        else
            $appRun $appArgs "$@" > $logFile 2>&1 || echo "ERROR" >> $logFile
        fi
    fi
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# getNumberOfProcessors
#      Extract 'numberOfSubdomains' from system/decomposeParDict
#      (or alternative location). Copied from OpenFOAM-v2012 and adapted to
#      work with all versions.
# Arguments:
#     1:decomposeParDict location
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::getNumberOfProcessors()
{
    local dict="${1:-system/decomposeParDict}"

    numberOfSubdomains=$(grep -i numberOfSubdomains "$dict" | awk '{print $2}' | sed 's/;//g')

    if [ "$#" -eq 1 ]
    then
        echo "$numberOfSubdomains"
    else
        echo "Error getting 'numberOfSubdomains' from '$dict'" 1>&2
        echo 1      # Fallback is 1 proc (serial)
        return 1
    fi
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# runParallel
#     runParallel that works that same irrespective of the OpenFOAM version
#     Copied from OpenFOAM-v2012
# Arguments:
#     None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::runParallel()
{
    local appName appRun logFile logMode nProcs

    # Any additional parsed arguments (eg, decomposeParDict)
    local appArgs="-parallel"

    local mpirun="mpirun"
    if [ "$FOAM_MPI" = msmpi ]
    then
        mpirun="mpiexec"
    fi

    # Parse options until executable is encountered
    while [ "$#" -gt 0 ] && [ -z "$appRun" ]
    do
        case "$1" in
            -a | -append)
                logMode=append
                ;;
            -o | -overwrite)
                logMode=overwrite
                ;;
            -s | -suffix)
                logFile=".$2"
                shift
                ;;
            -n | -np)
                nProcs="$2"
                shift
                ;;
            -decomposeParDict)
                appArgs="$appArgs $1 $2"
                nProcs=$(solids4Foam::getNumberOfProcessor "$2")
                shift
                ;;
            '')
                ;;
            *)
                appRun="$1"
                ;;
        esac
        shift
    done

    [ -n "$nProcs" ] || nProcs=$(solids4Foam::getNumberOfProcessors "system/decomposeParDict")

    appName="${appRun##*/}"
    logFile="log.$appName$logFile"

    if [ -f "$logFile" ] && [ -z "$logMode" ]
    then
        echo "$appName already run on $PWD:" \
             "remove log file '$logFile' to re-run"
    else
        echo "Running $appRun ($nProcs processes) on $PWD "
        # Options '-n' and '-np' are synonymous, but msmpi only supports '-n'
        if [ "$logMode" = append ]
        then
        (
            $mpirun -n $nProcs $appRun $appArgs "$@" </dev/null >> $logFile 2>&1
        )
        else
        (
            $mpirun -n $nProcs $appRun $appArgs "$@" </dev/null > $logFile 2>&1
        )
        fi
    fi
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# runSolidModel
#     Run the case using the specified solid model. This script modifies the
#     case to work with the specified solid model and then runs the case.
#     Modification are made to the solidProperties, fvSchemes, fvSolution,
#     boundary conditions, etc., as needed.
# Arguments:
#     1: Address of the case directory
#     2: Name of the solid model to use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::runSolidModel()
{
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4Foam::runSolidModel start                                   |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo

    # Check number of input parameters is correct
    if [[ $# -ne 2 ]]
    then
        solids4Foam::err "runSolidModel: incorrect number of parameters"
    fi

    # Give sensible names to the argument
    CASE_DIR=$1
    SOLID_MODEL=$2

    # Directory containing the template solid model dicts
    # Note: the Allwmake script replaces /Users/philipc/OpenFOAM/philipc-v2012/solids4foam/applications/scripts with $(pwd)
    DICTS_PARENT_DIR=/Users/philipc/OpenFOAM/philipc-v2012/solids4foam/applications/scripts
    DICTS_DIR=${DICTS_PARENT_DIR}/solidModelDicts/${SOLID_MODEL}

    # Check default dictionaries exist for the specified solid model
    if [ ! -d "${DICTS_DIR}" ]; then
        solids4Foam::err "Cannot find ${DICTS_DIR}. Do default dicts exist for this solid model?"
    fi

    # Copy default dictionaries for the specified solid model

    echo "Replacing solidProperties"
    \cp ${DICTS_DIR}/solidProperties ${CASE_DIR}/constant/

    echo "Replacing fvSolution"
    \cp ${DICTS_DIR}/fvSolution ${CASE_DIR}/system/

    echo "Replacing fvSchemes"
    \cp ${DICTS_DIR}/fvSchemes ${CASE_DIR}/system/

    # Add files to case, if the addFilesToCase directory is present
    if [ -d "${DICTS_DIR}/addFilesToCase" ]; then
        for f in "${DICTS_DIR}/addFilesToCase/"*; do
            if [ -f "$f" ]; then  # Ensure it's a file
                echo "Copying ${f} to ${CASE_DIR}"
                \cp "${f}" "${CASE_DIR}"
            fi
        done
    fi

    # Check the displacement file is correct
    if [ ! -f "${DICTS_DIR}/displacementName" ]; then
        solids4Foam::err "Cannot find ${DICTS_DIR}/displacementName. This file should contain the name of displacement field primary unknown"
    fi
    DISP_NAME_FILE="${DICTS_DIR}"/displacementName
    DISP=$(cat "${DISP_NAME_FILE}")

    if [ ! -f "${DICTS_DIR}/0/${DISP}" ]; then
        # The DISP field is not present. Check if we can copy a field that is
        # present

        if [ "${DISP}" = "D" ] && [ -f "${CASE_DIR}/0/DD" ]; then
            echo "Renaming ${CASE_DIR}/0/DD to ${CASE_DIR}/0/D"
            \mv "${CASE_DIR}/0/DD" "${CASE_DIR}/0/D"
            sed -i "s/object.*DD;/object D;/g" "${CASE_DIR}/0/D"
        fi

        if [ "${DISP}" = "DD" ] && [ -f "${CASE_DIR}/0/D" ]; then
            echo "Renaming ${CASE_DIR}/0/D to ${CASE_DIR}/0/DD"
            \mv "${CASE_DIR}/0/D" "${CASE_DIR}/0/DD"
            sed -i "s/object.*D;/object DD;/g" "${CASE_DIR}/0/DD"
        fi

        if [ "${DISP}" != "D" ] && [ "${DISP}" != "DD" ]; then
            solids4Foam::err "Cannot find ${DICTS_DIR}/0/${DISP} and it cannot be copied from D or DD"
        fi
    fi

    # Other checks - to-do
    # 1. solidModel may use D, DD, pointD or other
    # 2. Specify forms of boundary condition, e.g. solidTraction vs
    #    blockSolidTraction

    # Run the case
    if [ ! -f "${CASE_DIR}/Allrun" ]; then
        solids4Foam::err "Cannot find ${CASE_DIR}/Allrun"
    fi
    echo "Running ./Allrun on ${CASE_DIR}"
    ./Allrun &> log.Allrun

    echo
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4Foam::runSolidModel end                                     |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo
}
