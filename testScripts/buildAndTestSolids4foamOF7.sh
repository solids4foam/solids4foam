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
CONTAINER="autobuild-solids4foam--openfoam-v1812"
solids4foam::createContainer "openfoam-7" \
    "openfoam/openfoam7-graphical-apps" \
    "${CONTAINER}" \
    "/home/openfoam"

# Build
if solids4foam::build "${CONTAINER}" "openfoam-7" \
       "/opt/openfoam7/etc/bashrc"
then
    BUILD_STATUS=true
else
    BUILD_STATUS=false
fi

# Test
# if solids4foam::test "${CONTAINER}" "openfoam-7" \
#        "/opt/openfoam7/etc/bashrc"
# then
#     TEST_STATUS=true
# else
#     TEST_STATUS=false
# fi

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
