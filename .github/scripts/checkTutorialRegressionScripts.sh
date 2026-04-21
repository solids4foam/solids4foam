#!/bin/bash

set -eu

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

missing=0

while IFS= read -r allrun_file
do
    case_dir="$(dirname "$allrun_file")"
    regression_file="$case_dir/regressionTest.sh"

    if [ ! -f "$regression_file" ]
    then
        echo "Missing regressionTest.sh: $case_dir"
        missing=1
    fi
done < <(find tutorials -mindepth 2 -name Allrun | sort)

if [ "$missing" -ne 0 ]
then
    echo "One or more tutorial cases are missing regressionTest.sh files."
    exit 1
fi

echo "All tutorial cases have regressionTest.sh files."
