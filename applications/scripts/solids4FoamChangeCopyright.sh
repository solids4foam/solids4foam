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
#     solids4FoamChangeCopyright
#     adapted from foamChangeCopyright in foam-extend-4.1
#
# Description
#     Updates the header of .H and .C files.
#
#     To run the script on multiple files, use the following command:
#
#     $ find . \( -name "*.C" -o -name "*.H" \) | xargs -n 1 "./scriptPath.sh"
#
#------------------------------------------------------------------------------

# Check if the correct number of arguments are provided
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <file>"
    echo "Wrong number of arguments, expected 1 found $#"
    exit 1
fi

# Check if the file exists
if [[ ! -f "$1" ]]; then
    echo "File not found: $1"
    exit 1
fi

# Read original file
orig_file=$(< "$1")

perl -0777 -p -i -e '
s!(([#% ]*) *)(  =========.*?\n)\n.*?\n[#%]*\n.*?\n[#%]*\n.*?\n[#%]*\n!$1License
$1    This file is part of solids4foam.
$2
$1    solids4foam is free software: you can redistribute it and/or modify it
$1    under the terms of the GNU General Public License as published by the
$1    Free Software Foundation, either version 3 of the License, or (at your
$1    option) any later version.
$2
$1    solids4foam is distributed in the hope that it will be useful, but
$1    WITHOUT ANY WARRANTY; without even the implied warranty of
$1    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
$1    General Public License for more details.
$2
$1    You should have received a copy of the GNU General Public License
$1    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.
$2
!s;
s/[ \t]+$//mg;
' $1

# Read modified file
mod_file=$(< "$1")

# Compare original and modified file
if [ "$orig_file" != "$mod_file" ]; then
    echo "The header has been updated!"
    exit 1;
else
    echo "No changes have been made!"
    exit 0;
fi
