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
#     commit are copied. The build files are solids4foam's own and are kept:
#     Make/, Allwmake, BEAMFOAM_VERSION and this script. The commit is
#     recorded in src/beamModels/BEAMFOAM_VERSION.
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

# beamFoam path -> solids4foam path
syncPaths=(
    "src/wireBunchingModels:src/beamModels"
    "applications/utilities/createBeamMesh:applications/utilities/createBeamMesh"
    "applications/utilities/setInitialPositionBeam:applications/utilities/setInitialPositionBeam"
)

# Files owned by solids4foam inside the synchronised directories
keep=(Make Allwmake BEAMFOAM_VERSION syncBeamFoam.sh)

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

    # Remove the previous copy, keeping the solids4foam files
    find "$to" -mindepth 1 -maxdepth 1 \
        $(printf -- '-not -name %s ' "${keep[@]}") \
        -exec rm -rf {} +

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
