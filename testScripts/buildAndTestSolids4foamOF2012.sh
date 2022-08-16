#!/bin/bash +x

echo "****************************************"
echo "* User build script: start"
echo "****************************************"
echo


# Load solids4foam bash functions
testScriptsDir="${0%/*}"
source ${testScriptsDir}/solids4FoamFunctions.sh

# Change permissions to allow docker to write in the current directory
chmod -R o+w .

# Create container
CONTAINER="autobuild-solids4foam-openfoam-v2012"
solids4foam::createContainer "openfoam-v2012" \
    "philippic/openfoam-v2012-centos73.gfortran:release" \
    "${CONTAINER}" \
    "/home/dockeruser"

# Build
if solids4foam::build "${CONTAINER}" "openfoam-v2012" \
       "/opt/OpenFOAM/OpenFOAM-v2012/etc/bashrc"
then
    BUILD_STATUS=true
else
    BUILD_STATUS=false
fi

# Test
if solids4foam::test "${CONTAINER}" "openfoam-v2012" \
       "/opt/OpenFOAM/OpenFOAM-v2012/etc/bashrc" "Alltest.openfoam"
then
    TEST_STATUS=true
else
    TEST_STATUS=false
fi

# Remove container
solids4foam::removeContainer "${CONTAINER}"


# Write status to files so it can be read by the automated emails
echo "Create BUILD/TEST STATUS FILES"
echo "${BUILD_STATUS}" > ${WORKSPACE}/BUILD_STATUS
echo "${TEST_STATUS}" > ${WORKSPACE}/TEST_STATUS


echo "****************************************"
echo "* User build script: end"
echo "****************************************"
echo


if [[ "${BUILD_STATUS}" == "false" ]] || [[ "${TEST_STATUS}" == "false" ]]
then
    echo "Exit: build and/or test failed"; echo
    exit 1
fi
