#!/bin/bash
# This script was adapted from foamChangeCopyright in foam-extend-4.1

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

