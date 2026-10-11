#!/bin/bash
#
# Check that every tutorial directory containing an Allrun script also has a
# README.md file beside it.
#
# A small number of directories are documented elsewhere. They are listed
# explicitly below, each with the README that documents it and the reason it
# is exempt. Directories are not exempted because some ancestor has a
# README, so a new tutorial without a README is still reported.

set -eu

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

# Exceptions, one per line: <directory>|<documenting README>|<reason>
exceptions="\
tutorials/solids/hyperelasticity/rigidRotation/rotatingBlock|\
tutorials/solids/hyperelasticity/rigidRotation/README.md|\
sub-case documented by the shared rigidRotation README
tutorials/solids/hyperelasticity/rigidRotation/rotatingCylinder|\
tutorials/solids/hyperelasticity/rigidRotation/README.md|\
sub-case documented by the shared rigidRotation README
tutorials/solids/hyperelasticity/rigidRotation/rotatingSphere|\
tutorials/solids/hyperelasticity/rigidRotation/README.md|\
sub-case documented by the shared rigidRotation README
tutorials/solids/hyperelasticity/cylindricalPressureVessel/caseOptions/\
pressureDisplacement/cylinderLin/analyticalSolution/calculation|\
tutorials/solids/hyperelasticity/cylindricalPressureVessel/README.md|\
helper that generates the analytical solution, not a separate tutorial"

failed=0

# Returns 0 if the given directory is listed as an exception
is_exception()
{
    local dir="$1"
    local exc_dir exc_readme exc_reason

    while IFS='|' read -r exc_dir exc_readme exc_reason
    do
        if [ "$dir" = "$exc_dir" ]
        then
            return 0
        fi
    done <<< "$exceptions"

    return 1
}

# Validate the exception list itself, so that it cannot silently go stale
while IFS='|' read -r exc_dir exc_readme exc_reason
do
    if [ ! -f "$exc_dir/Allrun" ]
    then
        echo "Stale README exception (no Allrun): $exc_dir"
        failed=1
    elif [ -f "$exc_dir/README.md" ]
    then
        echo "Unneeded README exception (README.md exists): $exc_dir"
        failed=1
    fi

    if [ ! -f "$exc_readme" ]
    then
        echo "Missing documenting README for exception: $exc_readme"
        echo "    (exempts $exc_dir: $exc_reason)"
        failed=1
    fi
done <<< "$exceptions"

while IFS= read -r allrun_file
do
    case_dir="$(dirname "$allrun_file")"

    if [ -f "$case_dir/README.md" ] || is_exception "$case_dir"
    then
        continue
    fi

    echo "Missing README.md: $case_dir"
    failed=1
done < <(find tutorials -mindepth 2 -name Allrun | sort)

if [ "$failed" -ne 0 ]
then
    echo
    echo "One or more tutorial cases are missing README.md files, or the"
    echo "exception list in .github/scripts/checkTutorialReadmes.sh is out of"
    echo "date."
    echo "Add a README.md to each reported case. If a directory is a helper or"
    echo "a sub-case documented elsewhere, add it to the exception list with"
    echo "the documenting README and a reason."
    exit 1
fi

echo "All tutorial cases have README.md files."
