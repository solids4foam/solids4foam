#!/bin/bash

set -eu

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

missing=0

while IFS= read -r allrun_file
do
    case_dir="$(dirname "$allrun_file")"
    readme_file="$case_dir/README.md"

    if [ ! -f "$readme_file" ]
    then
        echo "Missing README.md: $case_dir"
        missing=1
    fi
done < <(find tutorials -mindepth 2 -name Allrun | sort)

if [ "$missing" -ne 0 ]
then
    echo "One or more tutorial cases are missing README.md files."
    exit 1
fi

echo "All tutorial cases have README.md files."
