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
        echo "FOAM EXTEND loaded: no changes made"; echo
        exit 0
    fi

    # 1. symmetryPlane in foam extend becomes symmetry in OpenFOAM

    echo "Changing symmetryPlane to symmetry in D*"; echo
    if [[ -n $(find "${CASE_DIR}" -name "D*") ]]
    then
        find "${CASE_DIR}" -name "D*" | xargs sed -i 's\symmetryPlane\symmetry\g'
    fi

    echo "Changing symmetryPlane to symmetry in pointD*"; echo
    if [[ -n $(find "${CASE_DIR}" -name "pointD*") ]]
    then
        find "${CASE_DIR}" -name "pointD*" | xargs sed -i 's\symmetryPlane\symmetry\g'
    fi

    echo "Changing symmetryPlane to symmetry in U"; echo
    if [[ -n $(find "${CASE_DIR}" -name "U") ]]
    then
        find "${CASE_DIR}" -name "U" | xargs sed -i 's\symmetryPlane\symmetry\g'
    fi

    echo "Changing symmetryPlane to symmetry in p"; echo
    if [[ -n $(find "${CASE_DIR}" -name "p") ]]
    then
        find "${CASE_DIR}" -name "p" | xargs sed -i 's\symmetryPlane\symmetry\g'
    fi

    echo "Changing symmetryPlane to symmetry in pointMotionU"; echo
    if [[ -n $(find "${CASE_DIR}" -name "pointMotionU") ]]
    then
        find "${CASE_DIR}" -name "pointMotionU" | xargs sed -i 's\symmetryPlane\symmetry\g'
    fi

    echo "Changing symmetryPlane to symmetry in blockMeshDict"; echo
    if [[ -n $(find "${CASE_DIR}" -name blockMeshDict) ]]
    then
        find "${CASE_DIR}" -name blockMeshDict | xargs sed -i 's\symmetryPlane\symmetry\g'
    fi

    echo "Changing symmetryPlane to symmetry in boundary"; echo
    if [[ -n $(find "${CASE_DIR}" -name boundary) ]]
    then
        find "${CASE_DIR}" -name boundary | xargs sed -i 's\symmetryPlane\symmetry\g'
    fi

    # 2. If found, move the blockMeshDict to the system directory 
    echo "2. Checking if the blockMeshDict exists"
    if [[ -f "${CASE_DIR}"/constant/polyMesh/blockMeshDict ]]
    then
        echo "    Moving constant/polyMesh/blockMeshDict to system"
        \mv "${CASE_DIR}"/constant/polyMesh/blockMeshDict "${CASE_DIR}"/system/
    fi
    if [[ -f "${CASE_DIR}"/constant/solid/polyMesh/blockMeshDict ]]
    then
        echo "    Moving constant/solid/polyMesh/blockMeshDict to system/solid"
        \mv "${CASE_DIR}"/constant/solid/polyMesh/blockMeshDict "${CASE_DIR}"/system/solid
    fi
    if [[ -f "${CASE_DIR}"/constant/fluid/polyMesh/blockMeshDict ]]
    then
        echo "    Moving constant/fluid/polyMesh/blockMeshDict to system/fluid"
        \mv "${CASE_DIR}"/constant/fluid/polyMesh/blockMeshDict "${CASE_DIR}"/system/fluid
    fi
    echo

    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4Foam::convertCaseFormat end                                  |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# convertCaseFormatFoamExtend
#     Converts a case from foam extend format to FOAM EXTEND format, regardless
#     of what OpenFOAM version is loaded
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

    echo "Changing symmetry to symmetryPlane in D*"; echo
    if [[ -n $(find "${CASE_DIR}" -name "D*") ]]
    then
        find "${CASE_DIR}" -name "D*" | xargs sed -i 's\symmetry;\symmetryPlane;\g'
    fi
 
    echo "Changing symmetry to symmetryPlane in pointD*"; echo
    if [[ -n $(find "${CASE_DIR}" -name "pointD*") ]]
    then
        find "${CASE_DIR}" -name "pointD*" | xargs sed -i 's\symmetry;\symmetryPlane;\g'
    fi

    echo "Changing symmetry to symmetryPlane in U"; echo
    if [[ -n $(find "${CASE_DIR}" -name "U") ]]
    then
        find "${CASE_DIR}" -name "U" | xargs sed -i 's\symmetry;\symmetryPlane;\g'
    fi

    echo "Changing symmetry to symmetryPlane in p"; echo
    if [[ -n $(find "${CASE_DIR}" -name "p") ]]
    then
        find "${CASE_DIR}" -name "p" | xargs sed -i 's\symmetry;\symmetryPlane;\g'
    fi

    echo "Changing symmetry to symmetryPlane in pointMotionU"; echo
    if [[ -n $(find "${CASE_DIR}" -name "pointMotionU") ]]
    then
        find "${CASE_DIR}" -name "pointMotionU" | xargs sed -i 's\symmetry;\symmetryPlane;\g'
    fi

    echo "Changing symmetry to symmetryPlane in blockMeshDict"; echo
    if [[ -n $(find "${CASE_DIR}" -name blockMeshDict) ]]
    then
        find "${CASE_DIR}" -name blockMeshDict | xargs sed -i 's\symmetry \symmetryPlane \g'
    fi

    echo "Changing symmetry to symmetryPlane in boundary"; echo
    if [[ -n $(find "${CASE_DIR}" -name boundary) ]]
    then
        find "${CASE_DIR}" -name boundary | xargs sed -i 's\symmetry;\symmetryPlane;\g'
    fi

    # 2. If found, move the blockMeshDict to the system directory 
    echo "2. Checking if the blockMeshDict exists"
    if [[ -f "${CASE_DIR}"/system/blockMeshDict ]]
    then
        echo "    Moving system/blockMeshDict to constant/polyMesh"
        mkdir "${CASE_DIR}"/constant/polyMesh
        \mv "${CASE_DIR}"/system/blockMeshDict "${CASE_DIR}"/constant/polyMesh
    fi
    if [[ -f "${CASE_DIR}"/system/solid/blockMeshDict ]]
    then
        echo "    Moving system/solid/blockMeshDict to constant/solid/polyMesh"
        mkdir "${CASE_DIR}"/constant/solid/polyMesh
        \mv "${CASE_DIR}"/system/solid/blockMeshDict "${CASE_DIR}"/constant/solid/polyMesh
    fi
    if [[ -f "${CASE_DIR}"/system/fluid/blockMeshDict ]]
    then
        echo "    Moving system/fluid/blockMeshDict to constant/fluid/polyMesh"
        mkdir "${CASE_DIR}"/constant/fluid/polyMesh
        \mv "${CASE_DIR}"/system/fluid/blockMeshDict "${CASE_DIR}"/constant/fluid/polyMesh
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
