/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

InClass
    Foam::mechanicalLaw

\*---------------------------------------------------------------------------*/

#include "mechanicalLaw.H"
#include "volFields.H"
#include "fvc.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mechanicalLaw, 0);
    defineRunTimeSelectionTable(mechanicalLaw, dictionary);
}

// * * * * * * * * * * * * * * Protected Members * * * * * * * * * * * * * * //

bool Foam::mechanicalLaw::planeStress() const
{
    if (mesh_.foundObject<IOdictionary>("mechanicalProperties"))
    {
        return
            Switch
            (
                mesh_.lookupObject<IOdictionary>
                (
                    "mechanicalProperties"
                ).lookup("planeStress")
            );
    }
    else
    {
        // It is not straight-forward to lookup the mechanicalProperties from
        // here as we only have access to a subMesh fvMesh objectRegistry
        // We will read it here again; this switch only gets called at the start
        // of a simulation so it is not a problem
        IOdictionary mechProp
        (
            IOobject
            (
                "mechanicalProperties",
                "constant",
                mesh_.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        return Switch(mechProp.lookup("planeStress"));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalLaw::mechanicalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    baseMeshRegionName_()
{
    // Set the base mesh region name
    // For an FSI case, the region will be called solid, else it will be called
    // region0.
    if (mesh.db().foundObject<fvMesh>("solid"))
    {
        baseMeshRegionName_ = "solid";
    }
    else if (mesh.db().foundObject<fvMesh>("region0"))
    {
        baseMeshRegionName_ = "region0";
    }
    else
    {
        FatalErrorIn
        (
            "Foam::mechanicalLaw::mechanicalLaw\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        ) << "solid region name not found" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * //


Foam::tmp<Foam::surfaceScalarField> Foam::mechanicalLaw::impKf() const
{
    return fvc::interpolate(impK());
}


void Foam::mechanicalLaw::correct(surfaceSymmTensorField&)
{
    notImplemented
    (
        type() + "::correct(surfaceSymmTensorField&)\n"
        "The correct(surfaceSymmTensorField&) function is not implemented\n"
        " for the " + type() + " mechanical law"
    );
}


Foam::scalar Foam::mechanicalLaw::residual()
{
    // Default to zero; this can be overwritten by any derived mechanical law
    return 0.0;
}


// ************************************************************************* //
