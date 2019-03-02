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
    const word& fromRegion,          // from region name
    const word& toRegion,            // to region name
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<Type>& fromField,    // from field
    Field<Type>& toField             // to field
) const
{
    // Check field sizes are correct

    if (fromField.size() != fromZone.size())
    {
        FatalErrorIn
        (
            "Foam::tmp< Field<Type> >\n"
            "Foam::fluidSolidInterface::transferFacesZoneToZone\n"
            "(\n"
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
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
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
            ") const"
        )   << "toField is wrong size!" << nl
            << "toField size: " << toField.size()
            << ", toZone size: " << toZone.size()
            << abort(FatalError);
    }

    if (transferMethod_ == directMap)
    {
        if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
        {
            const labelList& fluidToSolidMap = fluidToSolidFaceMap();
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
            const labelList& solidToFluidMap = solidToFluidFaceMap();
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
                "    const standAlonePatch& fromZone,\n"
                "    const standAlonePatch& toZone,\n"
                "    const Field<Type>& fromField,\n"
                "    Field<Type>& toField\n"
                ") const"
            )   << "Unknown regions:  " << fromRegion
                << " and/or " << toRegion << abort(FatalError);
        }
    }
    else if (transferMethod_ == RBF)
    {
        Info<< "Interpolating from " << fromRegion << " to " << toRegion
            << " using RBF interpolation" << endl;

        matrix fromRbfField(fromField.size(), 3);
        matrix toRbfField(toField.size(), 3);

        forAll(fromField, faceI)
        {
            fromRbfField(faceI, 0) = fromField[faceI].x();
            fromRbfField(faceI, 1) = fromField[faceI].y();
            fromRbfField(faceI, 2) = fromField[faceI].z();
        }

        if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
        {
            rbfFluidToSolid()->interpolate(fromRbfField, toRbfField);
        }
        else if
        (
            fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
        )
        {
            rbfSolidToFluid()->interpolate(fromRbfField, toRbfField);
        }
        else
        {
            FatalErrorIn
            (
                "Foam::tmp< Field<Type> >\n"
                "Foam::fluidSolidInterface::transferFacesZoneToZone\n"
                "(\n"
                "    const standAlonePatch& fromZone,\n"
                "    const standAlonePatch& toZone,\n"
                "    const Field<Type>& fromField,\n"
                "    Field<Type>& toField\n"
                ") const"
            )   << "Unknown regions:  " << fromRegion
                << " and/or " << toRegion << abort(FatalError);
        }

        forAll(toField, faceI)
        {
            toField[faceI].x() = toRbfField(faceI, 0);
            toField[faceI].y() = toRbfField(faceI, 1);
            toField[faceI].z() = toRbfField(faceI, 2);
        }
    }
    else if (transferMethod_ == GGI)
    {
        Info<< "Interpolating from " << fromRegion << " to " << toRegion
            << " using GGI/AMI interpolation" << endl;

        if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
        {
            // fluid is the master; solid is the slave
            toField = ggiInterpolator().masterToSlave(fromField);
        }
        else if
        (
            fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
        )
        {
            toField = ggiInterpolator().slaveToMaster(fromField);
        }
        else
        {
            FatalErrorIn
            (
                "Foam::tmp< Field<Type> >\n"
                "Foam::fluidSolidInterface::transferFacesZoneToZone\n"
                "(\n"
                "    const standAlonePatch& fromZone,\n"
                "    const standAlonePatch& toZone,\n"
                "    const Field<Type>& fromField,\n"
                "    Field<Type>& toField\n"
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
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
            ") const"
        )   << "Unknown transferMethod:  "
            << interfaceTransferMethodNames_[transferMethod_] << nl
            << "Available transfer methods are: "
            << interfaceTransferMethodNames_
            << abort(FatalError);
    }
}


template<class Type>
void Foam::fluidSolidInterface::transferPointsZoneToZone
(
    const word& fromRegion,          // from region name
    const word& toRegion,            // to region name
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<Type>& fromField,    // from field
    Field<Type>& toField             // to field
) const
{
    // Check field sizes are correct

    if (fromField.size() != fromZone.nPoints())
    {
        FatalErrorIn
        (
            "Foam::tmp< Field<Type> >\n"
            "Foam::fluidSolidInterface::transferPointsZoneToZone\n"
            "(\n"
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
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
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
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
            const labelList& fluidToSolidMap = fluidToSolidPointMap();
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
            const labelList& solidToFluidMap = solidToFluidPointMap();
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
                "Foam::fluidSolidInterface::transferPoiubtZoneToZone\n"
                "(\n"
                "    const standAlonePatch& fromZone,\n"
                "    const standAlonePatch& toZone,\n"
                "    const Field<Type>& fromField,\n"
                "    Field<Type>& toField\n"
                ") const"
            )   << "Unknown regions:  " << fromRegion
                << " and/or " << toRegion << abort(FatalError);
        }
    }
    else if (transferMethod_ == RBF)
    {
        Info<< "Interpolating from " << fromRegion << " to " << toRegion
            << " using RBF interpolation" << endl;

        matrix fromRbfField(fromField.size(), 3);
        matrix toRbfField(toField.size(), 3);

        forAll(fromField, faceI)
        {
            fromRbfField(faceI, 0) = fromField[faceI].x();
            fromRbfField(faceI, 1) = fromField[faceI].y();
            fromRbfField(faceI, 2) = fromField[faceI].z();
        }

        if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
        {
            rbfFluidToSolid()->interpolate(fromRbfField, toRbfField);
        }
        else if
        (
            fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
        )
        {
            rbfSolidToFluid()->interpolate(fromRbfField, toRbfField);
        }
        else
        {
            FatalErrorIn
            (
                "Foam::tmp< Field<Type> >\n"
                "Foam::fluidSolidInterface::transferPointsZoneToZone\n"
                "(\n"
                "    const standAlonePatch& fromZone,\n"
                "    const standAlonePatch& toZone,\n"
                "    const Field<Type>& fromField,\n"
                "    Field<Type>& toField\n"
                ") const"
            )   << "Unknown regions:  " << fromRegion
                << " and/or " << toRegion << abort(FatalError);
        }

        forAll(toField, faceI)
        {
            toField[faceI].x() = toRbfField(faceI, 0);
            toField[faceI].y() = toRbfField(faceI, 1);
            toField[faceI].z() = toRbfField(faceI, 2);
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
                ggiInterpolator().masterToSlavePointInterpolate(fromField);
        }
        else if
        (
            fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
        )
        {
            toField =
                ggiInterpolator().slaveToMasterPointInterpolate(fromField);
        }
        else
        {
            FatalErrorIn
            (
                "Foam::tmp< Field<Type> >\n"
                "Foam::fluidSolidInterface::transferPointsZoneToZone\n"
                "(\n"
                "    const standAlonePatch& fromZone,\n"
                "    const standAlonePatch& toZone,\n"
                "    const Field<Type>& fromField,\n"
                "    Field<Type>& toField\n"
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
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
            ") const"
        )   << "Unknown transferMethod:  "
            << interfaceTransferMethodNames_[transferMethod_] << nl
            << "Available transfer methods are: "
            << interfaceTransferMethodNames_
            << abort(FatalError);
    }
}


// ************************************************************************* //
