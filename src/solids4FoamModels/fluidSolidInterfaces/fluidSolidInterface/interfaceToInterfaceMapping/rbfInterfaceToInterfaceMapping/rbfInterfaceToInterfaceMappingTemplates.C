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

#include "rbfInterfaceToInterfaceMapping.H"

namespace Foam
{

namespace interfaceToInterfaceMappings
{

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

template<class Type>
void rbfInterfaceToInterfaceMapping::transferFacesZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<Type>& fromField,    // from field
    Field<Type>& toField             // to field
) const
{
    Info<< "Interpolating face values using RBF" << endl;

    // Check field sizes are correct
    interfaceToInterfaceMapping::checkFieldSizes
    (
        fromZone.size(), toZone.size(), fromField.size(), toField.size()
    );

    // Transfer fromField into matrix format
    matrix fromRbfField(fromField.size(), int(pTraits<Type>::nComponents));
    matrix toRbfField(toField.size(), int(pTraits<Type>::nComponents));
    for(int cmptI = 0; cmptI < pTraits<Type>::nComponents; cmptI++)
    {
        const scalarField fromFieldCmptI(fromField.component(cmptI));

        forAll(fromField, faceI)
        {
            fromRbfField(faceI, cmptI) = fromFieldCmptI[faceI];
        }
    }

    // Check if fromZone is zoneA or zoneB by checking the memory address
    if (&fromZone == &zoneA() && &toZone == &zoneB())
    {
        // fromZone is zoneA; toZone is zoneB
        zoneAToZoneBInterpolator()->interpolate(fromRbfField, toRbfField);
    }
    else if (&toZone == &zoneA() && &fromZone == &zoneB())
    {
        // toZone is zoneA; fromZone is zoneB
        zoneBToZoneAInterpolator()->interpolate(fromRbfField, toRbfField);
    }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp< Field<Type> >\n"
            "Foam::rbfInterfaceToInterfaceMapping::"
            "transferFacesZoneToZone\n"
            "(\n"
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
            ") const"
        )   << "Unknown zones!" << abort(FatalError);
    }

    // Transfer toField from matrix format into field
    for(int cmptI = 0; cmptI < pTraits<Type>::nComponents; cmptI++)
    {
        scalarField toFieldCmptI(toField.size(), 0.0);

        forAll(toField, faceI)
        {
            toFieldCmptI[faceI] = toRbfField(faceI, cmptI);
        }

        toField.replace(cmptI, toFieldCmptI);
    }
}


template<class Type>
void rbfInterfaceToInterfaceMapping::transferPointsZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<Type>& fromField,    // from field
    Field<Type>& toField             // to field
) const
{
    Info<< "Interpolating point values using RBF" << endl;

    // Check field sizes are correct
    interfaceToInterfaceMapping::checkFieldSizes
    (
        fromZone.nPoints(), toZone.nPoints(), fromField.size(), toField.size()
    );

    // Transfer fromField into matrix format
    matrix fromRbfField(fromField.size(), int(pTraits<Type>::nComponents));
    matrix toRbfField(toField.size(), int(pTraits<Type>::nComponents));
    for(int cmptI = 0; cmptI < pTraits<Type>::nComponents; cmptI++)
    {
        const scalarField fromFieldCmptI(fromField.component(cmptI));

        forAll(fromField, pointI)
        {
            fromRbfField(pointI, cmptI) = fromFieldCmptI[pointI];
        }
    }

    // Check if fromZone is zoneA or zoneB by checking the memory address
    if (&fromZone == &zoneA() && &toZone == &zoneB())
    {
        // fromZone is zoneA; toZone is zoneB
        zoneAToZoneBInterpolator()->interpolate(fromRbfField, toRbfField);
    }
    else if (&toZone == &zoneA() && &fromZone == &zoneB())
    {
        // toZone is zoneA; fromZone is zoneB
        zoneBToZoneAInterpolator()->interpolate(fromRbfField, toRbfField);
    }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp< Field<Type> >\n"
            "Foam::rbfInterfaceToInterfaceMapping::"
            "transferPointsZoneToZone\n"
            "(\n"
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
            ") const"
        )   << "Unknown zones!" << abort(FatalError);
    }

    // Transfer toField from matrix format into field
    for(int cmptI = 0; cmptI < pTraits<Type>::nComponents; cmptI++)
    {
        scalarField toFieldCmptI(toField.size(), 0.0);

        forAll(toField, pointI)
        {
            toFieldCmptI[pointI] = toRbfField(pointI, cmptI);
        }

        toField.replace(cmptI, toFieldCmptI);
    }
}

// ************************************************************************* //

} // End namespace interfaceToInterfaceMappings

} // End namespace Foam


// ************************************************************************* //
