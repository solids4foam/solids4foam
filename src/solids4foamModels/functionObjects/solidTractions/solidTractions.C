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

#include "solidTractions.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidTractions, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidTractions,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidTractions::writeData()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>("region0");

    if (mesh.foundObject<volSymmTensorField>(stressName_))
    {
        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>(stressName_);

        volVectorField traction
        (
            IOobject
            (
                "traction",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimForce/dimArea, vector::zero)
        );

        const surfaceVectorField n = mesh.Sf()/mesh.magSf();

        if (nonLinearTL_)
        {
            const volTensorField Finv =
                hinv(mesh.lookupObject<volTensorField>("F"));

            forAll(traction.boundaryField(), patchi)
            {
                if (!traction.boundaryField()[patchi].coupled())
                {
                    traction.boundaryField()[patchi] =
                        (
                            Finv.boundaryField()[patchi].T()
                            & n.boundaryField()[patchi]
                        ) & sigma.boundaryField()[patchi];
                }
            }
        }
        else
        {
            forAll(traction.boundaryField(), patchi)
            {
                if (!traction.boundaryField()[patchi].coupled())
                {
                    traction.boundaryField()[patchi] =
                        n.boundaryField()[patchi]
                        & sigma.boundaryField()[patchi];
                }
            }
        }

        traction.write();
    }
    else
    {
        InfoIn(this->name() + " function object constructor")
            << stressName_ << " not found" << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidTractions::solidTractions
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    stressName_
    (
        dict.found("stressName")
        ? word(dict.lookup("stressName"))
        : word("sigma")
    ),
    nonLinearTL_
    (
        dict.found("nonLinearTL")
        ? bool(dict.lookup("nonLinearTL"))
        : bool(false)
    )
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidTractions::start()
{
    return writeData();
}


bool Foam::solidTractions::execute()
{
    return writeData();
}


bool Foam::solidTractions::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
