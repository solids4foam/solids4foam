/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "readFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Mesh, class GeoField>
void Foam::readFields
(
    const Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeoField>& fields
)
{
    // Search list of objects for volScalarFields
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

    // Remove the cellDist field
    IOobjectList::iterator celDistIter = fieldObjects.find("cellDist");
    if (celDistIter != fieldObjects.end())
    {
        fieldObjects.erase(celDistIter);
    }

    // Construct the vol scalar fields
    fields.setSize(fieldObjects.size());

    label fieldI = 0;
    for
    (
        IOobjectList::iterator iter = fieldObjects.begin();
        iter != fieldObjects.end();
        ++iter
    )
    {
        fields.set
        (
            fieldI++,
            new GeoField
            (
                *iter(),
                mesh
            )
        );
    }
}


// ************************************************************************* //
