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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# convertCaseFormat
#     Converts a case from foam extend format to OpenFOAM format. No changes are
#     applied if foam extend is loaded.
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
    CASE_DIR=$1

    # Exit if foam extend is loaded
    if [[ $WM_PROJECT = "foam" ]]
    then
        echo "foam-extend loaded: no changes made"; echo
        return 0
    fi

    # 1. symmetryPlane in foam extend becomes symmetry in OpenFOAM

    if [[ -n $(find "${CASE_DIR}" -name blockMeshDict*) ]]
    then
        echo "Changing symmetryPlane to symmetry in blockMeshDict"; echo
        find "${CASE_DIR}" -name blockMeshDict* | xargs sed -i 's\symmetryPlane\symmetry\g'
    fi

    if [[ -n $(find "${CASE_DIR}" -name boundary) ]]
    then
        echo "Changing symmetryPlane to symmetry in boundary"; echo
        find "${CASE_DIR}" -name boundary | xargs sed -i 's\symmetryPlane\symmetry\g'
    fi

    # Check all files in the 0 directory
    for FILE in $(find ./[0-9]* -type f)
    do
        if [[ -f "${FILE}" ]]
        then
            if grep -q "symmetryPlane;" "${FILE}"
            then
                echo "Changing symmetryPlane to symmetry in ${FILE}"; echo
                sed -i 's\symmetryPlane;\symmetry;\g' "${FILE}"
            fi
        fi
    done

    # 2. If found, move the blockMeshDict to the system directory
    if [[ -f "${CASE_DIR}"/constant/polyMesh/blockMeshDict ]]
    then
        echo "Moving constant/polyMesh/blockMeshDict to system"
        \mv "${CASE_DIR}"/constant/polyMesh/blockMeshDict "${CASE_DIR}"/system/
    fi
    if [[ -f "${CASE_DIR}"/constant/solid/polyMesh/blockMeshDict ]]
    then
        echo "Moving constant/solid/polyMesh/blockMeshDict to system/solid"
        \mv "${CASE_DIR}"/constant/solid/polyMesh/blockMeshDict "${CASE_DIR}"/system/solid
    fi
    if [[ -f "${CASE_DIR}"/constant/fluid/polyMesh/blockMeshDict ]]
    then
        echo "Moving constant/fluid/polyMesh/blockMeshDict to system/fluid"
        \mv "${CASE_DIR}"/constant/fluid/polyMesh/blockMeshDict "${CASE_DIR}"/system/fluid
    fi

    # Replace the functions file
    if [[ -f "${CASE_DIR}"/system/functions ]]
    then
        echo "Replacing system/functions with system/functions.openfoam"
        \cp "${CASE_DIR}"/system/functions \
            "${CASE_DIR}"/system/functions.foamextend
        \cp -f "${CASE_DIR}"/system/functions.openfoam \
            "${CASE_DIR}"/system/functions
    fi

    # 3. Rename turbulence model
    if [[ -n $(find "${CASE_DIR}" -name turbulenceProperties) ]]
    then
        echo "Changing RASModel to RAS in turbulenceProperties"
        find "${CASE_DIR}" -name turbulenceProperties | xargs sed -i "s/RASModel;/RAS;/g"
    fi

    # 4. Check for boundaryData
    if [[ -d "${CASE_DIR}"/constant/boundaryData && -d "${CASE_DIR}"/constant/boundaryData.openfoam ]]
    then
        echo "Moving constant/boundaryData to constant/boundaryData.foam-extend"
        \mv "${CASE_DIR}"/constant/boundaryData "${CASE_DIR}"/constant/boundaryData.foam-extend

        echo "Moving constant/boundaryData.openfoam to constant/boundaryData"
        \mv "${CASE_DIR}"/constant/boundaryData.openfoam "${CASE_DIR}"/constant/boundaryData
    fi

    echo
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4Foam::convertCaseFormat end                                  |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# convertCaseFormatFoamExtend
#     Converts a case to the FOAM EXTEND format, regardless of what OpenFOAM
#     version is loaded
# Arguments:
#     1: CASE_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4Foam::convertCaseFormatFoamExtend()
{
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4Foam::convertCaseFormatFoamExtend start                      |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo

    # Check number of input parameters is correct
    if [[ $# -ne 1 ]]
    then
        solids4Foam::err "convertCaseFormatFoamExtend: incorrect number of parameters"
    fi

    # Give sensible names to the argument
    CASE_DIR=$1

    # Un-do changes made in convertCaseFormat, if any

    # 1. symmetryPlane in foam extend becomes symmetry in OpenFOAM

    if [[ -n $(find "${CASE_DIR}" -name blockMeshDict) ]]
    then
        echo "Changing symmetry to symmetryPlane in blockMeshDict"; echo
        find "${CASE_DIR}" -name blockMeshDict | xargs sed -i 's\symmetry \symmetryPlane \g'
    fi

    if [[ -n $(find "${CASE_DIR}" -name boundary) ]]
    then
    echo "Changing symmetry to symmetryPlane in boundary"; echo
        find "${CASE_DIR}" -name boundary | xargs sed -i 's\symmetry;\symmetryPlane;\g'
    fi

    for FILE in $(find ./[0-9]* -type f)
    do
        if [[ -f "${FILE}" ]]
        then
            if grep -q "symmetry;" "${FILE}"
            then
                echo "Changing symmetry to symmetryPlane in ${FILE}"; echo
                sed -i 's\symmetry;\symmetryPlane;\g' "${FILE}"
            fi
        fi
    done

    # 2. If found, move the blockMeshDict to the system directory
    if [[ -f "${CASE_DIR}"/system/blockMeshDict ]]
    then
        echo "Moving system/blockMeshDict to constant/polyMesh"
        mkdir -p "${CASE_DIR}"/constant/polyMesh
        \mv "${CASE_DIR}"/system/blockMeshDict "${CASE_DIR}"/constant/polyMesh
    fi
    if [[ -f "${CASE_DIR}"/system/solid/blockMeshDict ]]
    then
        echo "Moving system/solid/blockMeshDict to constant/solid/polyMesh"
        mkdir -p "${CASE_DIR}"/constant/solid/polyMesh
        \mv "${CASE_DIR}"/system/solid/blockMeshDict "${CASE_DIR}"/constant/solid/polyMesh
    fi
    if [[ -f "${CASE_DIR}"/system/fluid/blockMeshDict ]]
    then
        echo "Moving system/fluid/blockMeshDict to constant/fluid/polyMesh"
        mkdir -p "${CASE_DIR}"/constant/fluid/polyMesh
        \mv "${CASE_DIR}"/system/fluid/blockMeshDict "${CASE_DIR}"/constant/fluid/polyMesh
    fi

    if [[ -f "${CASE_DIR}"/system/functions.foamextend ]]
    then
        echo "Replacing system/functions with system/functions.openfoam"
        \mv -f "${CASE_DIR}"/system/functions.foamextend \
            "${CASE_DIR}"/system/functions
    fi

    # 3. Rename turbulence model
    if [[ -n $(find "${CASE_DIR}" -name turbulenceProperties) ]]
    then
        echo "Changing RAS to RASModel in turbulenceProperties"
        find "${CASE_DIR}" -name turbulenceProperties | xargs sed -i "s/RAS;/RASModel;/g"
    fi

    # 4. Check for boundaryData
    if [[ -d "${CASE_DIR}"/constant/boundaryData && -d "${CASE_DIR}"/constant/boundaryData.foam-extend ]]
    then
        echo "Moving constant/boundaryData to constant/boundaryData.openfoam"
        \mv "${CASE_DIR}"/constant/boundaryData "${CASE_DIR}"/constant/boundaryData.openfoam

        echo "Moving constant/boundaryData.foam-extend to constant/boundaryData"
        \mv "${CASE_DIR}"/constant/boundaryData.foam-extend "${CASE_DIR}"/constant/boundaryData
    fi

    echo
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4Foam::convertCaseFormatFoamExtend end                        |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo
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
#     Give error if OpenFOAM version is not foam-extend
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
