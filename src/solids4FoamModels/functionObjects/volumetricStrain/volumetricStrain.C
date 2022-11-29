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

\*----------------------------------------------------------------------------*/

#include "volumetricStrain.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(volumetricStrain, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        volumetricStrain,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::volumetricStrain::writeData()
{
    if (runTime_.outputTime())
    {
        if (mesh_.foundObject<volTensorField>("grad(D)"))
        {
            // Lookup D field
            const volTensorField& gradD =  mesh_.lookupObject<volTensorField>("grad(D)");

            const volScalarField volEpsilon("volEpsilon", tr(gradD));

            // Write fields
            Info<< "    Writing volEpsilon" << endl;

            volEpsilon.write();
        }
        else 
        {
            Info<< name_ << ": grad(D) not found!" << endl;
        }

    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volumetricStrain::volumetricStrain
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    runTime_(t),
    mesh_
    (
        runTime_.lookupObject<fvMesh>
        (
            dict.lookupOrDefault<word>("region", "region0")
        )
    )
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::volumetricStrain::start()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


#if FOAMEXTEND
bool Foam::volumetricStrain::execute(const bool forceWrite)
#else
bool Foam::volumetricStrain::execute()
#endif
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


bool Foam::volumetricStrain::read(const dictionary& dict)
{
    return true;
}

#ifdef OPENFOAMESIORFOUNDATION
bool Foam::volumetricStrain::write()
{
    return writeData();
}
#endif

// ************************************************************************* //
