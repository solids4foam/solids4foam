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
    for FILE in $(find "${CASE_DIR}" -type f)
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

    # 5. Check for sample for foundation version
    if [[ "${WM_PROJECT_VERSION}" != *"v"* ]]
    then
        if [[ -f "${CASE_DIR}"/system/sample ]]
        then
           echo "OpenFOAM.org specific: replacing 'uniform' with 'lineUniform' in system/sample"
           sed -i "s/type.*uniform;/type lineUniform;/g" "${CASE_DIR}"/system/sample

           echo "OpenFOAM.org specific: replacing 'face' with 'lineFace' in system/sample"
           sed -i "s/type.*face;/type lineFace;/g" "${CASE_DIR}"/system/sample
        fi
    fi

    # 6. Check for timeVaryingUniformFixedValue
    if [[ -n $(find "${CASE_DIR}" -name p) ]]
    then
        echo "Changing timeVaryingUniformFixedValue to uniformValue in p"
        find "${CASE_DIR}" -name p | \
            xargs sed -i "s|//type.*uniformFixedValue;|type          uniformFixedValue;|g"
        find "${CASE_DIR}" -name p | \
            xargs sed -i "s|type.*timeVaryingUniformFixedValue;|//type        timeVaryingUniformFixedValue;|g"
    fi

    # 7. Check for changeDictionaryDict.openfoam
    if [[ -f "${CASE_DIR}/system/changeDictionaryDict.openfoam" ]]
    then
        echo "Moving ${CASE_DIR}/system/changeDictionaryDict to system/changeDictionaryDict.foamextend"
        mv "${CASE_DIR}/system/changeDictionaryDict" "system/changeDictionaryDict.foamextend"
        echo "Moving ${CASE_DIR}/system/changeDictionaryDict.openfoam to system/changeDictionaryDict"
        mv "${CASE_DIR}/system/changeDictionaryDict.openfoam" "system/changeDictionaryDict"
    fi

    # 8. Check for createPatchDict.openfoam
    if [[ -f "${CASE_DIR}/system/createPatchDict.openfoam" ]]
    then
        echo "Moving ${CASE_DIR}/system/createPatchDict to system/createPatchDict.foamextend"
        mv "${CASE_DIR}/system/createPatchDict" "system/createPatchDict.foamextend"
        echo "Moving ${CASE_DIR}/system/createPatchDict.openfoam to system/createPatchDict"
        mv "${CASE_DIR}/system/createPatchDict.openfoam" "system/createPatchDict"
    fi

    # 9. Either pointCellsLeastSquares or edgeCellsLeastSquares should be used
    #    the gradScheme for the solid in OpenFOAM.com, as these are the only
    #    schemes consistent with boundary non-orthogonal correction
    if [[ "${WM_PROJECT_VERSION}" == *"v"* ]]
    then
        if [[ -f "${CASE_DIR}"/constant/solidProperties ]]
        then
            echo "OpenFOAM.com specific: replacing 'leastSquares' with 'pointCellsLeastSquares' in system/fvSchemes"
            sed -i "s/ leastSquares;/ pointCellsLeastSquares;/g" "${CASE_DIR}"/system/fvSchemes
        elif [[ -f "${CASE_DIR}"/constant/solid/solidProperties ]]
        then
            echo "OpenFOAM.com specific: replacing 'leastSquares' with 'pointCellsLeastSquares' in system/solid/fvSchemes"
            sed -i "s/ leastSquares;/ pointCellsLeastSquares;/g" "${CASE_DIR}"/system/solid/fvSchemes
        fi
    fi

    # 10. Resolve force post-processing path from foam-extend
    if  [[ -n $(find "${CASE_DIR}" -name force.gnuplot) ]]
    then
        if [[ $WM_PROJECT_VERSION == *"v"* ]]
        then
            echo "Modifying force.gnuplot in consistent with ESI version"
            sed -i "s|forces.dat|force.dat|g" force.gnuplot
        fi
    fi

    # 11. Resolve sample post-processing path from foam-extend
    if  [[ -n $(find "${CASE_DIR}" -name plot.gnuplot) ]]
    then
        echo "Updating plot.gnuplot"
        sed -i "s|postProcessing/sets/|postProcessing/sample/|g" plot.gnuplot
    fi

    # 12. Resolve sampleDict post-processing path from foam-extend
    if  [[ -n $(find "${CASE_DIR}" -name plot.gnuplot) ]]
    then
        echo "Updating plot.gnuplot"
        sed -i  "s|postProcessing/surfaces/|postProcessing/sample.surfaces/|g" plot.gnuplot
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

    for FILE in $(find "${CASE_DIR}" -type f)
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

    # 5. Check for sample
    if [[ "${WM_PROJECT_VERSION}" != *"v"* ]]
    then
        if [[ -f "${CASE_DIR}"/system/sample ]]
        then
           echo "OpenFOAM.org specific: replacing 'lineUniform' with 'uniform' in system/sample"
           sed -i "s/type.*lineUniform;/type uniform;/g" "${CASE_DIR}"/system/sample

           echo "OpenFOAM.org specific: replacing 'lineFace' with 'face' in system/sample"
           sed -i "s/type.*lineFace;/type face;/g" "${CASE_DIR}"/system/sample
        fi
    fi

    # 6. Check for timeVaryingUniformFixedValue
    if [[ -n $(find "${CASE_DIR}" -name p) ]]
    then
        echo "Changing uniformValue to timeVaryingUniformFixedValue in p"
        find "${CASE_DIR}" -name p | \
            xargs sed -i "s|type.*uniformFixedValue;|//type          uniformFixedValue;|g"
        find "${CASE_DIR}" -name p | \
            xargs sed -i "s|//type.*timeVaryingUniformFixedValue;|type        timeVaryingUniformFixedValue;|g"

        # Remove any //// that were introdued
        find "${CASE_DIR}" -name p | xargs sed -i "s|////|//|g"
    fi

    # 7. Check for changeDictionaryDict.openfoam
    if [[ -f "${CASE_DIR}/system/changeDictionaryDict.foamextend" ]]
    then
        echo "Moving ${CASE_DIR}/system/changeDictionaryDict to system/changeDictionaryDict.openfoam"
        mv "${CASE_DIR}/system/changeDictionaryDict" "system/changeDictionaryDict.openfoam"
        echo "Moving ${CASE_DIR}/system/changeDictionaryDict.foamextend to system/changeDictionaryDict"
        mv "${CASE_DIR}/system/changeDictionaryDict.foamextend" "system/changeDictionaryDict"
    fi

    # 8. Check for createPatchDict.openfoam
    if [[ -f "${CASE_DIR}/system/createPatchDict.foamextend" ]]
    then
        echo "Moving ${CASE_DIR}/system/createPatchDict to system/createPatchDict.openfoam"
        mv "${CASE_DIR}/system/createPatchDict" "system/createPatchDict.openfoam"
        echo "Moving ${CASE_DIR}/system/createPatchDict.foamextend to system/createPatchDict"
        mv "${CASE_DIR}/system/createPatchDict.foamextend" "system/createPatchDict"
    fi

    # 9. Either pointCellsLeastSquares or edgeCellsLeastSquares should be used
    #    the gradScheme for the solid in OpenFOAM.com, as these are the only
    #    schemes consistent with boundary non-orthogonal correction
    if [[ "${WM_PROJECT_VERSION}" == *"v"* ]]
    then
        if [[ -f "${CASE_DIR}"/constant/solidProperties ]]
        then
            echo "OpenFOAM.com specific: replacing 'pointCellsLeastSquares' with 'leastSquares' in system/fvSchemes"
            sed -i "s/ pointCellsLeastSquares;/ leastSquares;/g" "${CASE_DIR}"/system/fvSchemes
        elif [[ -f "${CASE_DIR}"/constant/solid/solidProperties ]]
        then
            echo "OpenFOAM.com specific: replacing 'pointCellsLeastSquares' with 'leastSquares' in system/solid/fvSchemes"
            sed -i "s/ pointCellsLeastSquares;/ leastSquares;/g" "${CASE_DIR}"/system/solid/fvSchemes
        fi
    fi

    # 10. Resolve force post-processing path for foam-extend
    if  [[ -n $(find "${CASE_DIR}" -name force.gnuplot) ]]
    then
        if [[ $WM_PROJECT_VERSION == *"v"* ]]
        then
            echo "Reverting force.gnuplot from ESI version to foam-extend or .org "
            sed -i "s|force.dat|forces.dat|g" force.gnuplot
        fi
    fi

    # 11. Resolve sample post-processing path for foam-extend
    if  [[ -n $(find "${CASE_DIR}" -name plot.gnuplot) ]]
    then
        echo "Updating plot.gnuplot"
        sed -i "s|postProcessing/sample/|postProcessing/sets/|g" plot.gnuplot
    fi

    # 12. Resolve sampleDict post-processing path for foam-extend
    if  [[ -n $(find "${CASE_DIR}" -name plot.gnuplot) ]]
    then
        echo "Updating plot.gnuplot"
        sed -i "s|postProcessing/sample.surfaces/|postProcessing/surfaces/|g" plot.gnuplot
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
# caseOnlyRunsWithOpenFOAM
#     Give error if OpenFOAM version is foam-extend
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
    local appName appRun logFile logMode

    # Any additional parsed arguments (eg, decomposeParDict)
    local appArgs

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
                nProcs=$(getNumberOfProcessors "$2")
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

    [ -n "$nProcs" ] || nProcs=$(getNumberOfProcessors system/decomposeParDict)

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
