/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#ifdef OPENFOAMESIORFOUNDATION

#include "amiInterfaceToInterfaceMapping.H"

namespace Foam
{

namespace interfaceToInterfaceMappings
{

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

template<class Type>
void amiInterfaceToInterfaceMapping::transferFacesZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<Type>& fromField,    // from field
    Field<Type>& toField             // to field
) const
{
    Info<< "Interpolating face values using AMI" << endl;

    // Check field sizes are correct
    interfaceToInterfaceMapping::checkFieldSizes
    (
        fromZone.size(), toZone.size(), fromField.size(), toField.size()
    );

    // Check if fromZone is zoneA or zoneB by checking the memory address
    if (&fromZone == &zoneA() && &toZone == &zoneB())
    {
        // fromZone is the master (zoneA); toZone is the slave (zoneB)
        toField = interpolator().interpolateToTarget(fromField);
    }
    else if (&toZone == &zoneA() && &fromZone == &zoneB())
    {
        // toZone is the master (zoneA); fromZone is the slave (zoneB)
        toField = interpolator().interpolateToSource(fromField);
    }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp< Field<Type> >\n"
            "Foam::amiInterfaceToInterfaceMapping::transferFacesZoneToZone\n"
            "(\n"
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
            ") const"
        )   << "Unknown zones!" << abort(FatalError);
    }
}


template<class Type>
void amiInterfaceToInterfaceMapping::transferPointsZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<Type>& fromField,    // from field
    Field<Type>& toField             // to field
) const
{
    Info<< "Interpolating point values using AMI" << endl;

    // Check field sizes are correct
    interfaceToInterfaceMapping::checkFieldSizes
    (
        fromZone.nPoints(), toZone.nPoints(), fromField.size(), toField.size()
    );

    // Check if fromZone is zoneA or zoneB by checking the memory address
    if (&fromZone == &zoneA() && &toZone == &zoneB())
    {
        // fromZone is the master (zoneA); toZone is the slave (zoneB)
        //toField = interpolator().interpolateToTargetPoints(fromField);

        // Escape the interpolation if there are no faces in the target patch
        if (zoneB().nPoints() == 0)
        {
            return;
        }

        Field<Type>& result = toField;

        const List<face>& zoneAFaces = zoneA().localFaces();

        const List<labelPair>& addr = zoneBPointAddr();

        const FieldField<Field, double>& weights = zoneBPointWeights();

        forAll(result, pointI)
        {
            if (addr[pointI].first() > -1)
            {
                const face& hitFace = zoneAFaces[addr[pointI].first()];

                const label pI = addr[pointI].second();

                const Type ctrF = average(Field<Type>(fromField, hitFace));

                result[pointI] =
                    weights[pointI][0]*fromField[hitFace[pI]]
                  + weights[pointI][1]*fromField[hitFace.nextLabel(pI)]
                  + weights[pointI][2]*ctrF;
            }
            else
            {
                FatalErrorIn
                (
                    "void\n"
                    "Foam::amiInterfaceToInterfaceMapping::transferPointsZoneToZone\n"
                    "(\n"
                    "    ...\n"
                    ") const"
                )   << "zoneB point addressing is not correct"
                    << abort(FatalError);
            }
        }
    }
    else if (&toZone == &zoneA() && &fromZone == &zoneB())
    {
        // toZone is the master (zoneA); fromZone is the slave (zoneB)
        //toField = interpolateToSourcePoints(interpolator(), fromField);

        // Escape the interpolation if there are no faces in the target patch
        if (zoneA().nPoints() == 0)
        {
            return;
        }

        Field<Type>& result = toField;

        const List<face>& zoneBFaces = zoneB().localFaces();

        const List<labelPair>& addr = zoneAPointAddr();

        const FieldField<Field, double>& weights = zoneAPointWeights();

        forAll(result, pointI)
        {
            if (addr[pointI].first() > -1)
            {
                const face& hitFace = zoneBFaces[addr[pointI].first()];

                const label pI = addr[pointI].second();

                const Type ctrF = average(Field<Type>(fromField, hitFace));

                result[pointI] =
                    weights[pointI][0]*fromField[hitFace[pI]]
                  + weights[pointI][1]*fromField[hitFace.nextLabel(pI)]
                  + weights[pointI][2]*ctrF;
            }
            else
            {
                FatalErrorIn
                (
                    "void\n"
                    "Foam::amiInterfaceToInterfaceMapping::transferPointsZoneToZone\n"
                    "(\n"
                    "    ...\n"
                    ") const"
                )   << "zoneA point addressing is not correct"
                    << abort(FatalError);
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void\n"
            "Foam::amiInterfaceToInterfaceMapping::transferPointsZoneToZone\n"
            "(\n"
            "    ...\n"
            ") const"
        )   << "Unknown zones!" << abort(FatalError);
    }
}

// ************************************************************************* //

} // End namespace interfaceToInterfaceMappings

} // End namespace Foam

#endif // end of #ifdef OPENFOAMESIORFOUNDATION

// ************************************************************************* //
