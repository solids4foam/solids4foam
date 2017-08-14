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

Application
    volPointInterpolate

Description
    Interpolate the given field from cell centres to mesh vertices using least
    squares interpolation procedure.

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volPointInterpolateField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("fieldName");
    argList::validOptions.insert("noMeshUpdate", "");

#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    // Read arguments
    const word fieldName(args.additionalArgs()[0]);
    const bool noMeshUpdate = args.optionFound("noMeshUpdate");

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"

    for (label i = startTime; i < endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        if (!noMeshUpdate)
        {
            Info << "    Reading mesh" << endl;
            mesh.readUpdate();
        }

        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Create pointMesh
        pointMesh pMesh(mesh);

        // Create interpolator
        newLeastSquaresVolPointInterpolation volToPointInterp(mesh);

        // Check sigma exists
        if (fieldHeader.headerOk())
        {
            mesh.readUpdate();

            bool processed = false;

            volPointInterpolateField<scalar>
            (
                fieldHeader, mesh, pMesh, volToPointInterp, processed
            );
            volPointInterpolateField<vector>
            (
                fieldHeader, mesh, pMesh, volToPointInterp, processed
            );
            volPointInterpolateField<sphericalTensor>
            (
                fieldHeader, mesh, pMesh, volToPointInterp, processed
            );
            volPointInterpolateField<symmTensor>
            (
                fieldHeader, mesh, pMesh, volToPointInterp, processed
            );
            volPointInterpolateField<tensor>
            (
                fieldHeader, mesh, pMesh, volToPointInterp, processed
            );

            if (!processed)
            {
                FatalError
                    << "Unable to process " << fieldName << nl
                    << "Not implemented for fields of type "
                    << fieldHeader.headerClassName() << nl << nl
                    << exit(FatalError);
            }
        }
        else
        {
            Info<< "    No field " << fieldName << endl;
        }

        Info<< endl;
    }

    Info<< "End" << endl;

    return(0);
}


// ************************************************************************* //
