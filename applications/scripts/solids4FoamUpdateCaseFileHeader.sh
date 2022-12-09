#!/bin/sh
#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     | Version:     4.1
#   \\  /    A nd           | Web:         http://www.foam-extend.org
#    \\/     M anipulation  | For copyright notice see file Copyright
#------------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     solids4foamUpdateCaseFileHeader
#     adapted from foamUpdateCaseFileHeader in foam-extend-4.1
#
# Description
#     Updates the header of application files.
#     By default, writes current version in the header.
#     Alternatively version can be specified with -v option.
#     Also removes consecutive blank lines from file.
#
#------------------------------------------------------------------------------
solids4foamVersion=v2.0

usage() {
    cat<<USAGE

Usage: ${0##*/} [OPTION] <file1> ... <fileN>

options:
  -v VER  specifies the version to be written in the header
  -h      help

  Updates the header of application files and removes consecutive blank lines.
  By default, writes version ${solids4foamVersion} in the header.
  An alternative version can be specified with the -v option.

USAGE
    exit 1
}

printHeader() {
    cat<<HEADER
/*--------------------------------*- C++ -*----------------------------------*\\
| solids4foam: solid mechanics and fluid-solid interaction simulations        |
| Version:     ${solids4foamVersion}                                                           |
| Web:         https://solids4foam.github.io                                  |
| Disclaimer:  This offering is not approved or endorsed by OpenCFD Limited,  |
|              producer and distributor of the OpenFOAM software via          |
|              www.openfoam.com, and owner of the OPENFOAM® and OpenCFD®      |
|              trade marks.                                                   |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ${FORMAT};
    class       ${CLASS};
HEADER

    if [ -n "${NOTE}" ];
    then
        cat<<HEADER
    note        ${NOTE};
HEADER
    fi

    if [ -n "${LOCATION}" ];
    then
        cat<<HEADER
    location    ${LOCATION};
HEADER
    fi

        cat<<HEADER
    object      ${OBJECT};
}
HEADER
}

#
# extract attribute '$1' from file '$2'
#
FoamFileAttribute() {
    sed -n -e 's/[ ;]*$//' -e "s/^ *$1 *//p" $2
}

#
# OPTIONS
#
opts=$(getopt hv: $*)
if [ $? -ne 0 ]
then
    echo "Aborting due to invalid option"
    usage
fi
eval set -- '$opts'
while [ "$1" != "--" ]
do
    case $1 in
    -v)
        foamVersion=$2
        shift
        ;;
    -h)
        usage
        ;;
    esac
    shift
done
shift

[ $# -ge 1 ] || usage

# constant width for version
foamVersion=$(printf %-33s $foamVersion)

#
# MAIN
#

for caseFile
do
    if [ ! -x "$caseFile" ] && (grep "^ *FoamFile" $caseFile >/dev/null 2>&1)
    then
        echo "Updating case file: $caseFile"
        sed -n '/FoamFile/,/}/p' $caseFile > FoamFile.tmp

        FORMAT=$(FoamFileAttribute format FoamFile.tmp)
        CLASS=$(FoamFileAttribute  class  FoamFile.tmp)
        NOTE=$(FoamFileAttribute note FoamFile.tmp)
        LOCATION=$(FoamFileAttribute location FoamFile.tmp)
        OBJECT=$(FoamFileAttribute object FoamFile.tmp)

        printHeader > FoamFile.tmp
        sed '1,/}/d' $caseFile | sed '/./,/^$/!d' | sed 's/ *$//g' >> FoamFile.tmp
        #sed '1,/}/d' $caseFile >> FoamFile.tmp

        # use cat to avoid removing/replace soft-links
        [ -s FoamFile.tmp ] && cat FoamFile.tmp >| $caseFile
        rm -f FoamFile.tmp 2>/dev/null
    else
        echo " Invalid case file: $caseFile" 1>&2
    fi
done

#------------------------------------------------------------------------------
