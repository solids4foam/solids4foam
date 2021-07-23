#!/bin/bash
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Library
#     Bash functions for building and testing the solids4foam toolbox using
#     docker containers, modified to be used by Jenkins
#     The script broadly follows the style given at
#     https://google.github.io/styleguide/shellguide.html as well as the OpenFOAM
#     coding style at https://openfoam.org/dev/coding-style-guide
#     The script is checked with https://www.shellcheck.net
#
# Authors
#     Philip Cardiff, UCD, February 2021
#
# License
#     GNU Lesser General Public License, version 3.
#     https://www.gnu.org/licenses/lgpl-3.0.en.html
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# createContainer
#     Create and start a docker container
# Arguments:
#     1: VERSION
#     2: IMAGE
#     3: NAME_OF_CONTAINER_TO_BE_CREATED
#     4: CONTAINER_HOME_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4foam::createContainer()
{
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4foam::createContainer start                                 |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo

    # Check number of input parameters is correct
    if [[ $# -ne 4 ]]
    then
        echo; echo "solids4foam::createContainer incorrect number of parameters"
        echo
        exit 1
    fi

    # Give sensible names to the arguments and local variables
    VERSION="${1}"
    IMAGE="${2}"
    CONTAINER="${3}"
    HOME_DIR="${4}"

    echo "Version: ${VERSION}"
    echo "Image: ${IMAGE}"
    echo "Container: ${CONTAINER}"
    echo "Container home directory: ${HOME_DIR}"
    echo

    # Pull latest image
    echo; echo "Pulling latest docker image: ${IMAGE}"; echo
    docker pull "${IMAGE}"

    # Create temporary container
    echo; echo "Creating temporary docker container"; echo
    CWD=$(pwd)
    docker create --network=host --entrypoint /bin/bash --mount \
        src="${CWD}",target="${HOME_DIR}"/solids4foam-release,type=bind \
        --name "${CONTAINER}" -it "${IMAGE}"

    # Start the temporary container
    echo "Starting the temporary docker container"; echo
    docker start "${CONTAINER}"

    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4foam::createContainer end                                   |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# build
#     Compile solids4foam using the given docker container
# Arguments:
#     1: CONTAINER
#     2: VERSION
#     3: FOAM_BASH_ADDRESS_IN_CONTAINER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4foam::build()
{
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4foam::build start                                           |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo

    # Check number of input parameters is correct
    if [[ $# -ne 3 ]]
    then
        echo; echo "solids4foam::build incorrect number of parameters"; echo
        exit 1
    fi

    # Give sensible names to the arguments and local variables
    CONTAINER="${1}"
    VERSION="${2}"
    BASH_ADDRESS="${3}"

    echo "Container: ${CONTAINER}"
    echo "Version: ${VERSION}"
    echo "Bash address: ${BASH_ADDRESS}"
    echo

    echo "Compiling solids4foam with ${VERSION} via docker"; echo
    buildPassed=false
    BUILD_LOG_FILE="Allwmake.${VERSION}.log"
    if docker exec "${CONTAINER}" bash -c \
       "source ${BASH_ADDRESS} &> /dev/null \
       && cd solids4foam-release \
       && ./Allwmake &> ${BUILD_LOG_FILE}"
    then
        # Also check the log file for the string "Error " if the build fails
        if [[ $(grep -c "Error " "${BUILD_LOG_FILE}") -gt 0 ]]
        then
            echo; echo "solids4foam ${VERSION}: Allwmake failed *************"; echo
            buildPassed=false
        else
            echo; echo "solids4foam ${VERSION}: Allwmake passed"; echo
            buildPassed=true
        fi
    else
        echo; echo "solids4foam ${VERSION}: Alltest failed **************"; echo
        buildPassed=false
    fi

    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4foam::build end                                              |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo

    if [ "${buildPassed}" = true ]
    then
        return 0
    else
        return 1
    fi
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# test
#     Run solids4foam Alltest script
# Arguments:
#     1: CONTAINER
#     2: VERSION
#     3: FOAM_BASH_ADDRESS_IN_CONTAINER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4foam::test()
{
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4foam::test start                                            |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo

    # Check number of input parameters is correct
    if [[ $# -ne 4 ]]
    then
        echo; echo "solids4foam::test incorrect number of parameters"; echo
        exit 1
    fi

    # Give sensible names to the arguments and local variables
    CONTAINER="${1}"
    VERSION="${2}"
    BASH_ADDRESS="${3}"
    ALLTEST_SCRIPT="${4}"

    echo "Container: ${CONTAINER}"
    echo "Version: ${VERSION}"
    echo "Bash address: ${BASH_ADDRESS}"
    echo "Alltest script: ${ALLTEST_SCRIPT}"
    echo

    # Print the number of tutorials and the number only supported by foam extend
    NUM_TUTS=$(find . -name Allrun -type f | wc -l)
    NUM_TUTS_FE_ONLY=$(find . -name Allrun -type f | xargs grep "solids4Foam::caseOnlyRunsWithFoamExtend" | wc -l)
    echo "Number of tutorials: ${NUM_TUT}"
    echo "Number of tutorials limit to foam extend: ${NUM_TUTS_FE_ONLY}"
    echo

    # Run Alltest and check if all tutorials passed
    echo "Running ${ALLTEST_SCRIPT}"
    docker exec "${CONTAINER}" bash -c \
        "source ${BASH_ADDRESS} &> /dev/null \
        && cd solids4foam-release/tutorials \
        && ./${ALLTEST_SCRIPT} &> ../Alltest.${VERSION}.log \
        && cd .. && chmod -R 777 . &> chmod.${VERSION}.log"

    # The generated testLoopReport will contain the word "FATAL" for any
    # failed commands
    testPassed=false
    TEST_REPORT_FILE="tutorialsTest/testLoopReport"
    if [[ -f  "${TEST_REPORT_FILE}" ]]
    then
        mkdir -p logs
        cp "${TEST_REPORT_FILE}" "logs/testLoopReport.${VERSION}.log"
        cp tutorialsTest/logs "logs/logs.${VERSION}.log"

        if [[ $(grep -c "FATAL" "${TEST_REPORT_FILE}") -gt 0 ]]
        then
            echo; echo "solids4foam ${VERSION}: Alltest failed **************"; echo
        else
            echo; echo "solids4foam ${VERSION}: Alltest passed"; echo
            testPassed=true
        fi
    else
        echo; echo "solids4foam ${VERSION}: Alltest failed (no testLoopReport)";
        echo
    fi

    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4foam::test end                                              |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo

    if [ "${testPassed}" = true ]
    then
        return 0
    else
        return 1
    fi
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# removeContainer
#     Remove docker container and remove solids4foam directory
# Arguments:
#     1: CONTAINER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function solids4foam::removeContainer()
{
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4foam::removeContainer start                                 |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo

    # Check number of input parameters is correct
    if [[ $# -ne 1 ]]
    then
        echo; echo "solids4foam::removeContainer incorrect number of parameters"; echo
        exit 1
    fi

    # Give sensible names to the arguments and local variables
    CONTAINER="${1}"

    echo "Container: ${CONTAINER}"
    echo

    # Clean
    echo "Running Allwclean"
    docker exec "${CONTAINER}" bash -c \
        "source ${BASH_ADDRESS} &> /dev/null \
        && cd solids4foam-release \
        && ./Allwclean &> Allwclean.${VERSION}.log \
        && git clean -fd --exclude='*.log' --exclude='logs' \
        && rm -rf tutorialsTest &> rmTutorialsTest.${VERSION}.log"
    rm -rf tutorialsTest &> "rmTutorialsTest.sys.${VERSION}.log"
    git clean -fd --exclude='*.log' --exclude='logs'

    # Remove container
    echo "Removing ${CONTAINER} container"; echo
    docker stop "${CONTAINER}"
    docker rm --force "${CONTAINER}"
    
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "| solids4foam::removeContainer end                                   |"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo
}
