#!/bin/bash
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Library
#     Bash function to perform automated parallel scaling study on a given case.
#     The case to run and the numbers of cores to run the case on are given as
#     arguments (see example below). The script then performs the following
#     steps:
#     - Makes a copy of the case to a directory parallelScalingCase/nCoresX
#       where X is the number of cores
#     - Updates the number of cores in system/decomposeParDict, and also
#       system/solid/decomposeParDict and system/fluid/decomposeParDict if they
#       are present
#     - Execute the Allrun script passing the number of cores as an argument
#       e.g. "Allrun 4"
#     - Writes the execution time from log.solids4foam to the std out and
#       parallelScalingCase/scaling.txt

#     If a run.pbs file is found then this is submitted (using sbatch run.pbs)
#     and the Allrun script is ignored. In that case, the execution times are
#     not written.
#
#     The script assumes the following:
#     - An Allrun script is present in the case
#     - The Allrun script takes the number of cores as an argument
#     - The Allrun script will create a log called log.solids4foam
#
#     The script has been checked with https://www.shellcheck.net.
#
#     Example usage:
#
#         ./parallelScalingStudy plateHole 1 2 4 8 16 32 64
#
# Authors
#     Philip Cardiff, UCD, October 2021
#
# License
#     GNU Lesser General Public License, version 3.
#     https://www.gnu.org/licenses/lgpl-3.0.en.html
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "| parallelScalingStudy start                                         |"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo

# Check that at least two arguments has been provided
if [[ $# -lt 2 ]]
then
    echo "Error: incorrect number of parameters";
    echo "Exmaple Usage: ./parallelScalingStudy plateHole 1 2 4 8 16 32 64"
    echo
    exit 1
fi

# Give a sensible name to the first argument
CASE_DIR="${1}"
echo "Case: ${CASE_DIR}"; echo

# Create the directory for the parallel scaling cases
SCALING_DIR="parallelScalingCases"
if [[ -d "${SCALING_DIR}" ]]
then
    echo "Error: the ${SCALING_DIR} directory already exists! Please remove it"
    echo
    exit 1
fi
mkdir "${SCALING_DIR}"

# Create scaling file for results
SCALING_FILE="${SCALING_DIR}/scaling.txt"
echo "# nCores timeInSeconds" >> "${SCALING_FILE}"

# Loop through the list of core counts
for N_CORES in "${@:2}"
do
    echo "Number of cores: ${N_CORES}"

    # Set numberOfSubdomains in decomposeParDict
    if [[ -f "${CASE_DIR}/system/decomposeParDict" ]]
    then
        echo "Setting numberOfSubdomains to ${N_CORES} in system/decomposeParDict"
        sed -i "s/numberOfSubdomains .*$/numberOfSubdomains ${N_CORES};/g" \
            "${CASE_DIR}/system/decomposeParDict"
    else
        echo; echo "Error: cannot find system/decomposeParDict"
        echo
        exit 1
    fi

    # Set numberOfSubdomains in solid/decomposeParDict
    if [[ -f "${CASE_DIR}/system/solid/decomposeParDict" ]]
    then
        echo "Setting numberOfSubdomains to ${N_CORES} in system/solid/decomposeParDict"
        sed -i "s/numberOfSubdomains .*$/numberOfSubdomains ${N_CORES};/g" \
            "${CASE_DIR}/system/solid/decomposeParDict"
    fi

    # Set numberOfSubdomains in fluid/decomposeParDict
    if [[ -f "${CASE_DIR}/system/fluid/decomposeParDict" ]]
    then
        echo "Setting numberOfSubdomains to ${N_CORES} in system/fluid/decomposeParDict"
        sed -i "s/numberOfSubdomains .*$/numberOfSubdomains ${N_CORES};/g" \
            "${CASE_DIR}/system/fluid/decomposeParDict"
    fi

    # Check if an Allrun script exists
    if [[ -f "${CASE_DIR}/run.pbs" ]]
    then
        PBS_FILE="{CASE_DIR}/run.pbs"
    elif [[ ! -f "${CASE_DIR}/Allrun" ]]
    then
        echo; echo "Error: cannot find the Allrun script"
        echo
        exit 1
    fi

    # Make a copy of the case
    SCALING_CASE="${SCALING_DIR}/nCores${N_CORES}"
    mkdir "${SCALING_CASE}"
    echo "Copying the case to ${SCALING_CASE}"
    rsync -a --exclude "parallelScalingCases" * "${SCALING_CASE}/"

    # Execute the Allrun script
    if [[ -z "${PBS_FILE}" ]]
    then
        echo "Executing ${SCALING_CASE}/Allrun ${N_CORES}"
        if ! (cd "${SCALING_CASE}" && ./Allrun "${N_CORES}" &> log.Allrun)
        then
            echo; echo "Error: Allrun exited with an error"; echo
            exit 1
        else
            echo "Allrun exited cleanly"
        fi
    else
        echo "Submitting ${SCALING_CASE}/run.pbs"
        if ! (cd "${SCALING_CASE}" && sbatch run.pbs)
        then
            echo; echo "Error: 'sbatch run.pbs' returned an error"; echo
            exit 1
        else
            echo "run.pbs was submitted cleanly"

            # Skip to the next case
            continue
        fi
    fi

    # Print the execution time from log.solids4foam
    if [[ -f "${SCALING_CASE}/log.solids4foam" ]]
    then
        if [[ $OSTYPE == 'darwin'* ]]
        then
            EXEC_TIME=$(tail -r "${SCALING_CASE}/log.solids4Foam" | grep -m 1 "ExecutionTime" | awk '{print $3}')
        else
            EXEC_TIME=$(tac "${SCALING_CASE}/log.solids4Foam" | grep -m 1 "ExecutionTime" | awk '{print $3}')
        fi

        # Print to std out
        echo "Execution time = ${EXEC_TIME} s"

        # Write to file
        echo "${N_CORES} ${EXEC_TIME}" >> "${SCALING_FILE}"
    else
        echo; echo "${SCALING_CASE}/log.solids4foam does not exist: skipping"
    fi
    echo
done


echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "| parallelScalingStudy end                                           |"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo
