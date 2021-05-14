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
CONTAINER="autobuild-solids4foam-foam-extend-4.0"
solids4foam::createContainer "foam-extend-4.0" \
    "philippic/foam-extend-4.0-ubuntu18.04:latest" \
    "${CONTAINER}" \
    "/home/app"


# Build
if solids4foam::build "${CONTAINER}" "foam-extend-4.0" \
       "/home/app/foam/foam-extend-4.0/etc/bashrc"
then
    BUILD_STATUS=true
else
    BUILD_STATUS=false
fi


# Test
if solids4foam::test "${CONTAINER}" "foam-extend-4.0" \
       "/home/app/foam/foam-extend-4.0/etc/bashrc"
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
