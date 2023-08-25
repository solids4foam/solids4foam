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
    if (runTime_.outputTime() && mesh_.foundObject<volTensorField>("grad(D)"))
    {
        // Lookup the displacement gradient field
        const volTensorField& gradD =
            mesh_.lookupObject<volTensorField>("grad(D)");

        // Calculate volumetric strain
        volScalarField volEpsilon("volEpsilon", tr(gradD));

        if (compressionPositive_)
        {
            volEpsilon = -volEpsilon;
        }

        // Write the field
        Info<< name_ << ": writing volEpsilon" << endl;
        volEpsilon.write();
    }
    else
    {
        Info<< name_ << ": grad(D) not found!" << endl;
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
    ),
    compressionPositive_
    (
        dict.lookupOrDefault("compressionPositive", false)
    )
{
    Info<< "Creating " << this->name() << " function object" << nl
        << "    compressionPositive: " << compressionPositive_ << endl;
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

#ifdef OPENFOAM_NOT_EXTEND
bool Foam::volumetricStrain::write()
{
    return false;
}
#endif

// ************************************************************************* //
