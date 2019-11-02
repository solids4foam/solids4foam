/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fluidSolidInterface::transferFacesZoneToZone
(
    const word& fromRegion,                    // from region name
    const word& toRegion,                      // to region name
    const PtrList<globalPolyPatch>& fromZones, // from zones
    const PtrList<globalPolyPatch>& toZones,   // to zones
    const List<Field<Type> >& fromFields,      // from fields
    List<Field<Type> >& toFields               // to fields
) const
{
    forAll(solid().globalPatches(), i)
    {
        const standAlonePatch& fromZone = fromZones[i].globalPatch();
        const standAlonePatch& toZone = toZones[i].globalPatch();

        const Field<Type>& fromField = fromFields[i];
        Field<Type>& toField = toFields[i];

        // Check field sizes are correct
        if (fromField.size() != fromZone.size())
        {
            FatalErrorIn
            (
                "Foam::tmp< Field<Type> >\n"
                "Foam::fluidSolidInterface::transferFacesZoneToZone\n"
                "(\n"
                "    const word& fromRegion,\n"
                "    const word& toRegion,\n"
                "    const PtrList<globalPolyPatch>& fromZones,\n"
                "    const PtrList<globalPolyPatch>& toZones,\n"
                "    const List<Field<Type> >& fromFields,\n"
                "    List<Field<Type> >& toFields\n"
                ") const"
            )   << "fromField is wrong size!" << nl
                << "fromField size: " << fromField.size()
                << ", fromZone size: " << fromZone.size()
                << abort(FatalError);
        }

        if (toField.size() != toZone.size())
        {
            FatalErrorIn
            (
                "Foam::tmp< Field<Type> >\n"
                "Foam::fluidSolidInterface::transferFacesZoneToZone\n"
                "(\n"
                "    const word& fromRegion,\n"
                "    const word& toRegion,\n"
                "    const PtrList<globalPolyPatch>& fromZones,\n"
                "    const PtrList<globalPolyPatch>& toZones,\n"
                "    const List<Field<Type> >& fromFields,\n"
                "    List<Field<Type> >& toFields\n"
                ") const"
            )   << "toField is wrong size!" << nl
                << "toField size: " << toField.size()
                << ", toZone size: " << toZone.size()
                << abort(FatalError);
        }

        if (transferMethod_ == directMap)
        {
            Info<< "Interpolating from " << fromRegion << " to " << toRegion
                << " using direct mapping" << endl;

            if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
            {
                const labelList& fluidToSolidMap = fluidToSolidFaceMaps()[i];
                forAll(toField, faceI)
                {
                    toField[faceI] = fromField[fluidToSolidMap[faceI]];
                }
            }
            else if
            (
                fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
            )
            {
                const labelList& solidToFluidMap = solidToFluidFaceMaps()[i];
                forAll(toField, faceI)
                {
                    toField[faceI] = fromField[solidToFluidMap[faceI]];
                }
            }
            else
            {
                FatalErrorIn
                (
                    "Foam::tmp< Field<Type> >\n"
                    "Foam::fluidSolidInterface::transferFacesZoneToZone\n"
                    "(\n"
                    "    const word& fromRegion,\n"
                    "    const word& toRegion,\n"
                    "    const PtrList<globalPolyPatch>& fromZones,\n"
                    "    const PtrList<globalPolyPatch>& toZones,\n"
                    "    const List<Field<Type> >& fromFields,\n"
                    "    List<Field<Type> >& toFields\n"
                    ") const"
                )   << "Unknown regions:  " << fromRegion
                    << " and/or " << toRegion << abort(FatalError);
            }
        }
        else if (transferMethod_ == RBF)
        {
            Info<< "Interpolating from " << fromRegion << " to " << toRegion
                << " using RBF interpolation" << endl;

            matrix fromRbfField(fromField.size(), int(pTraits<Type>::nComponents));
            matrix toRbfField(toField.size(), int(pTraits<Type>::nComponents));

            for(int cmptI = 0; cmptI < pTraits<Type>::nComponents; cmptI++)
            {
                const scalarField fromFieldCmptI = fromField.component(cmptI);

                forAll(fromField, faceI)
                {
                    fromRbfField(faceI, cmptI) = fromFieldCmptI[faceI];
                }
            }

            if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
            {
                rbfFluidToSolid()[i]->interpolate(fromRbfField, toRbfField);
            }
            else if
            (
                fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
            )
            {
                rbfSolidToFluid()[i]->interpolate(fromRbfField, toRbfField);
            }
            else
            {
                FatalErrorIn
                (
                    "Foam::tmp< Field<Type> >\n"
                    "Foam::fluidSolidInterface::transferFacesZoneToZone\n"
                    "(\n"
                    "    const word& fromRegion,\n"
                    "    const word& toRegion,\n"
                    "    const PtrList<globalPolyPatch>& fromZones,\n"
                    "    const PtrList<globalPolyPatch>& toZones,\n"
                    "    const List<Field<Type> >& fromFields,\n"
                    "    List<Field<Type> >& toFields\n"
                    ") const"
                )   << "Unknown regions:  " << fromRegion
                    << " and/or " << toRegion << abort(FatalError);
            }

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
        else if (transferMethod_ == GGI)
        {
            Info<< "Interpolating from " << fromRegion << " to " << toRegion
                << " using GGI/AMI interpolation" << endl;

            if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
            {
                // fluid is the master; solid is the slave
                toField = ggiInterpolators()[i].masterToSlave(fromField);
            }
            else if
            (
                fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
            )
            {
                toField = ggiInterpolators()[i].slaveToMaster(fromField);
            }
            else
            {
                FatalErrorIn
                (
                    "Foam::tmp< Field<Type> >\n"
                    "Foam::fluidSolidInterface::transferFacesZoneToZone\n"
                    "(\n"
                    "    const word& fromRegion,\n"
                    "    const word& toRegion,\n"
                    "    const PtrList<globalPolyPatch>& fromZones,\n"
                    "    const PtrList<globalPolyPatch>& toZones,\n"
                    "    const List<Field<Type> >& fromFields,\n"
                    "    List<Field<Type> >& toFields\n"
                    ") const"
                )   << "Unknown regions:  " << fromRegion
                    << " and/or " << toRegion << abort(FatalError);
            }
        }
        else
        {
            FatalErrorIn
            (
                "Foam::tmp< Field<Type> >\n"
                "Foam::fluidSolidInterface::transferFacesZoneToZone\n"
                "(\n"
                "    const word& fromRegion,\n"
                "    const word& toRegion,\n"
                "    const PtrList<globalPolyPatch>& fromZones,\n"
                "    const PtrList<globalPolyPatch>& toZones,\n"
                "    const List<Field<Type> >& fromFields,\n"
                "    List<Field<Type> >& toFields\n"
                ") const"
            )   << "Unknown transferMethod:  "
                << interfaceTransferMethodNames_[transferMethod_] << nl
                << "Available transfer methods are: "
                << interfaceTransferMethodNames_
                << abort(FatalError);
        }
    }
}


template<class Type>
void Foam::fluidSolidInterface::transferPointsZoneToZone
(
    const word& fromRegion,                    // from region name
    const word& toRegion,                      // to region name
    const PtrList<globalPolyPatch>& fromZones, // from zones
    const PtrList<globalPolyPatch>& toZones,   // to zones
    const List<Field<Type> >& fromFields,      // from fields
    List<Field<Type> >& toFields               // to fields
) const
{
    forAll(solid().globalPatches(), i)
    {
        const standAlonePatch& fromZone = fromZones[i].globalPatch();
        const standAlonePatch& toZone = toZones[i].globalPatch();

        const Field<Type>& fromField = fromFields[i];
        Field<Type>& toField = toFields[i];

        // Check field sizes are correct
        if (fromField.size() != fromZone.nPoints())
        {
            FatalErrorIn
            (
                "Foam::tmp< Field<Type> >\n"
                "Foam::fluidSolidInterface::transferPointsZoneToZone\n"
                "(\n"
                "    const word& fromRegion,\n"
                "    const word& toRegion,\n"
                "    const PtrList<globalPolyPatch>& fromZones,\n"
                "    const PtrList<globalPolyPatch>& toZones,\n"
                "    const List<Field<Type> >& fromFields,\n"
                "    List<Field<Type> >& toFields\n"
                ") const"
            )   << "fromField is wrong size!" << nl
                << "fromField size: " << fromField.size()
                << ", fromZone size: " << fromZone.nPoints()
                << abort(FatalError);
        }

        if (toField.size() != toZone.nPoints())
        {
            FatalErrorIn
            (
                "Foam::tmp< Field<Type> >\n"
                "Foam::fluidSolidInterface::transferPointsZoneToZone\n"
                "(\n"
                "    const word& fromRegion,\n"
                "    const word& toRegion,\n"
                "    const PtrList<globalPolyPatch>& fromZones,\n"
                "    const PtrList<globalPolyPatch>& toZones,\n"
                "    const List<Field<Type> >& fromFields,\n"
                "    List<Field<Type> >& toFields\n"
                ") const"
            )   << "toField is wrong size!" << nl
                << "toField size: " << toField.size()
                << ", toZone size: " << toZone.nPoints()
                << abort(FatalError);
        }

        if (transferMethod_ == directMap)
        {
            Info<< "Interpolating from " << fromRegion << " to " << toRegion
                << " using direct mapping" << endl;

            if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
            {
                const labelList& fluidToSolidMap = fluidToSolidPointMaps()[i];
                forAll(toField, pointI)
                {
                    toField[pointI] = fromField[fluidToSolidMap[pointI]];
                }
            }
            else if
            (
                fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
            )
            {
                const labelList& solidToFluidMap = solidToFluidPointMaps()[i];
                forAll(toField, pointI)
                {
                    toField[pointI] = fromField[solidToFluidMap[pointI]];
                }
            }
            else
            {
                FatalErrorIn
                (
                    "Foam::tmp< Field<Type> >\n"
                    "Foam::fluidSolidInterface::transferPointsZoneToZone\n"
                    "(\n"
                    "    const word& fromRegion,\n"
                    "    const word& toRegion,\n"
                    "    const PtrList<globalPolyPatch>& fromZones,\n"
                    "    const PtrList<globalPolyPatch>& toZones,\n"
                    "    const List<Field<Type> >& fromFields,\n"
                    "    List<Field<Type> >& toFields\n"
                    ") const"
                )   << "Unknown regions:  " << fromRegion
                    << " and/or " << toRegion << abort(FatalError);
            }
        }
        else if (transferMethod_ == RBF)
        {
            Info<< "Interpolating from " << fromRegion << " to " << toRegion
                << " using RBF interpolation" << endl;

            matrix fromRbfField(fromField.size(), int(pTraits<Type>::nComponents));
            matrix toRbfField(toField.size(), int(pTraits<Type>::nComponents));

            for(int cmptI = 0; cmptI < pTraits<Type>::nComponents; cmptI++)
            {
                const scalarField fromFieldCmptI = fromField.component(cmptI);

                forAll(fromField, faceI)
                {
                    fromRbfField(faceI, cmptI) = fromFieldCmptI[faceI];
                }
            }

            if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
            {
                rbfFluidToSolid()[i]->interpolate(fromRbfField, toRbfField);
            }
            else if
            (
                fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
            )
            {
                rbfSolidToFluid()[i]->interpolate(fromRbfField, toRbfField);
            }
            else
            {
                FatalErrorIn
                (
                    "Foam::tmp< Field<Type> >\n"
                    "Foam::fluidSolidInterface::transferPointsZoneToZone\n"
                    "(\n"
                    "    const word& fromRegion,\n"
                    "    const word& toRegion,\n"
                    "    const PtrList<globalPolyPatch>& fromZones,\n"
                    "    const PtrList<globalPolyPatch>& toZones,\n"
                    "    const List<Field<Type> >& fromFields,\n"
                    "    List<Field<Type> >& toFields\n"
                    ") const"
                )   << "Unknown regions:  " << fromRegion
                    << " and/or " << toRegion << abort(FatalError);
            }

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
        else if (transferMethod_ == GGI)
        {
            Info<< "Interpolating from " << fromRegion << " to " << toRegion
                << " using GGI/AMI interpolation" << endl;

            if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
            {
                // fluid is the master; solid is the slave
                toField =
                    ggiInterpolators()[i].masterToSlavePointInterpolate(fromField);
            }
            else if
            (
                fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
            )
            {
                toField =
                    ggiInterpolators()[i].slaveToMasterPointInterpolate(fromField);
            }
            else
            {
                FatalErrorIn
                (
                    "Foam::tmp< Field<Type> >\n"
                    "Foam::fluidSolidInterface::transferPointsZoneToZone\n"
                    "(\n"
                    "    const word& fromRegion,\n"
                    "    const word& toRegion,\n"
                    "    const PtrList<globalPolyPatch>& fromZones,\n"
                    "    const PtrList<globalPolyPatch>& toZones,\n"
                    "    const List<Field<Type> >& fromFields,\n"
                    "    List<Field<Type> >& toFields\n"
                    ") const"
                )   << "Unknown regions:  " << fromRegion
                    << " and/or " << toRegion << abort(FatalError);
            }
        }
        else
        {
            FatalErrorIn
            (
                "Foam::tmp< Field<Type> >\n"
                "Foam::fluidSolidInterface::transferPointsZoneToZone\n"
                "(\n"
                "    const word& fromRegion,\n"
                "    const word& toRegion,\n"
                "    const PtrList<globalPolyPatch>& fromZones,\n"
                "    const PtrList<globalPolyPatch>& toZones,\n"
                "    const List<Field<Type> >& fromFields,\n"
                "    List<Field<Type> >& toFields\n"
                ") const"
            )   << "Unknown transferMethod:  "
                << interfaceTransferMethodNames_[transferMethod_] << nl
                << "Available transfer methods are: "
                << interfaceTransferMethodNames_
                << abort(FatalError);
        }
    }
}


// ************************************************************************* //
