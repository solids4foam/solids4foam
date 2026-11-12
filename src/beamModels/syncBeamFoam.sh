#!/bin/bash
#------------------------------------------------------------------------------
# License
#     This file is part of solids4foam.
#
#     solids4foam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     solids4foam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     syncBeamFoam.sh
#
# Description
#     Copy the beam library and beam mesh utilities from a beamFoam commit
#     into solids4foam, unchanged.
#
#     The sources are taken with git archive, so only files tracked at that
#     commit are copied. Before copying, the files that beamFoam tracked at the
#     previously synchronised commit, recorded in
#     src/beamModels/BEAMFOAM_VERSION, are removed; any other file in the
#     synchronised directories is kept. The build files are solids4foam's own
#     and are never replaced: Make/, Allwmake, BEAMFOAM_VERSION and this
#     script. Editor configuration (.vscode) is not copied. The new commit is
#     recorded in src/beamModels/BEAMFOAM_VERSION.
#
#     The beamFoam repository must contain both the previous and the new
#     commit.
#
# Usage
#     src/beamModels/syncBeamFoam.sh <beamFoam git repository> [commit]
#
#     The commit defaults to HEAD of the given repository.
#
#------------------------------------------------------------------------------

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]
then
    echo "Usage: $0 <beamFoam git repository> [commit]"
    exit 1
fi

beamFoamRepo="$1"
commit=$(git -C "$beamFoamRepo" rev-parse --verify "${2:-HEAD}^{commit}")
s4fRoot=$(cd "${0%/*}/../.." && pwd)

# The commit that the current copy comes from
previous=$(awk '$1 == "commit" { print $2 }' \
    "$s4fRoot/src/beamModels/BEAMFOAM_VERSION")

if ! git -C "$beamFoamRepo" cat-file -e "${previous}^{commit}" 2> /dev/null
then
    echo "The previously synchronised commit '$previous' (BEAMFOAM_VERSION)"
    echo "is not in $beamFoamRepo; fetch it first"
    exit 1
fi

# beamFoam path -> solids4foam path
syncPaths=(
    "src/wireBunchingModels:src/beamModels"
    "applications/utilities/createBeamMesh:applications/utilities/createBeamMesh"
    "applications/utilities/setInitialPositionBeam:applications/utilities/setInitialPositionBeam"
)

# Files owned by solids4foam inside the synchronised directories
keep=(Make Allwmake BEAMFOAM_VERSION syncBeamFoam.sh)

# beamFoam files that are not copied, at any depth (editor configuration)
exclude=(.vscode)

# Whether a path relative to a synchronised directory is solids4foam's own or
# excluded, so that the sync neither removes nor copies it
skipPath()
{
    local path="$1" name component components

    for name in "${keep[@]}"
    do
        [[ "${path%%/*}" == "$name" ]] && return 0
    done

    IFS=/ read -r -a components <<< "$path"
    for component in "${components[@]}"
    do
        for name in "${exclude[@]}"
        do
            [[ "$component" == "$name" ]] && return 0
        done
    done

    return 1
}

tmpDir=$(mktemp -d)
trap 'rm -rf "$tmpDir"' EXIT

for pair in "${syncPaths[@]}"
do
    from="${pair%%:*}"
    to="$s4fRoot/${pair##*:}"

    echo "Syncing $from -> ${pair##*:}"

    rm -rf "$tmpDir/archive"
    mkdir -p "$tmpDir/archive" "$to"
    git -C "$beamFoamRepo" archive "$commit" "$from" \
        | tar -x -C "$tmpDir/archive"

    for name in "${exclude[@]}"
    do
        find "$tmpDir/archive" -name "$name" -prune -exec rm -rf {} +
    done

    # Remove the files of the previous copy, keeping any other file
    while IFS= read -r -d '' path
    do
        path="${path#"$from"/}"
        if ! skipPath "$path"
        then
            rm -f "$to/$path"
        fi
    done < <(git -C "$beamFoamRepo" ls-tree -r -z --name-only \
        "$previous" -- "$from")

    find "$to" -mindepth 1 -type d -empty -delete

    # Copy the new sources, skipping the beamFoam build files
    (
        cd "$tmpDir/archive/$from"
        find . -mindepth 1 -maxdepth 1 \
            $(printf -- '-not -name %s ' "${keep[@]}") \
            -exec cp -a {} "$to"/ \;
    )

    # Report beamFoam source list changes, which Make/files must follow
    if [[ -f "$tmpDir/archive/$from/Make/files" && -f "$to/Make/files" ]]
    then
        if ! diff -qwB \
            <(grep -v '^ *\(LIB\|EXE\) *=' "$tmpDir/archive/$from/Make/files") \
            <(grep -v '^ *\(LIB\|EXE\) *=' "$to/Make/files") > /dev/null
        then
            echo "    NOTE: ${pair##*:}/Make/files differs from beamFoam:"
            diff -wB \
                <(grep -v '^ *\(LIB\|EXE\) *=' "$tmpDir/archive/$from/Make/files") \
                <(grep -v '^ *\(LIB\|EXE\) *=' "$to/Make/files") \
                | sed 's/^/        /' || true
        fi
    fi
done

cat > "$s4fRoot/src/beamModels/BEAMFOAM_VERSION" << EOF
# beamFoam commit synchronised by syncBeamFoam.sh
# Repository: https://github.com/solids4foam/beamFoam
commit $commit
date $(git -C "$beamFoamRepo" log -1 --format=%cs "$commit")
EOF

echo "Synchronised beamFoam $commit"

#------------------------------------------------------------------------------
