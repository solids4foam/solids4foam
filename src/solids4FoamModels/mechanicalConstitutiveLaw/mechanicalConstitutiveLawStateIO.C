/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mechanicalConstitutiveLawStateIO.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::mechanicalConstitutiveLawStateIO::fieldName
(
    const word& lawName,
    const word& topologyName,
    const wordList& childPath,
    const word& variableName
)
{
    // Every part is needed to make the name unique. The law, because two
    // materials may each declare epsilonP and mean different things by it. The
    // topology, because the same law may hold state at cell centres and at
    // faces in the same run, and the two are different lengths. The child
    // path, because a composite law and its sub-law each have their own
    // unqualified namespace and may both use the same name
    word name(lawName + ':' + topologyName);

    forAll(childPath, i)
    {
        name += ':' + childPath[i];
    }

    return name + ':' + variableName;
}


Foam::label Foam::mechanicalConstitutiveLawStateIO::totalSize
(
    const stateParts& parts
)
{
    label n = 0;

    forAll(parts, partI)
    {
        n += parts[partI]->size();
    }

    return n;
}


// ************************************************************************* //
