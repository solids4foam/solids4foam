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

\*---------------------------------------------------------------------------*/

#include "mechanicalModel.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mechanicalModel, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalModel::mechanicalModel(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "mechanicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    mechanicalLawPtr_(mechanicalLaw::New("law", mesh, subDict("mechanical"))),
    planeStress_(lookup("planeStress"))
{
    Info<< "Creating mechanical model" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //


Foam::mechanicalModel::~mechanicalModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::fvMesh& Foam::mechanicalModel::mesh() const
{
    return mesh_;
}


const Foam::Switch& Foam::mechanicalModel::planeStress() const
{
    return planeStress_;
}


const Foam::mechanicalLaw& Foam::mechanicalModel::law() const
{
    return mechanicalLawPtr_();
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::rho() const
{
    return mechanicalLawPtr_->rho();
}


void Foam::mechanicalModel::correct(Foam::volSymmTensorField& sigma)
{
    mechanicalLawPtr_->correct(sigma);
}


void Foam::mechanicalModel::correct(Foam::surfaceSymmTensorField& sigma)
{
    mechanicalLawPtr_->correct(sigma);
}


Foam::scalar Foam::mechanicalModel::residual()
{
    return mechanicalLawPtr_->residual();
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::impK() const
{
    return mechanicalLawPtr_->impK();
}


Foam::tmp<Foam::surfaceScalarField> Foam::mechanicalModel::impKf() const
{
    return fvc::interpolate(mechanicalLawPtr_->impK());
}


void Foam::mechanicalModel::updateTotalFields()
{
    mechanicalLawPtr_->updateTotalFields();
}


// Foam::tmp<Foam::volScalarField>
// Foam::mechanicalModel::plasticDissipationRate() const
// {
//     return mechanicalLawPtr_->plasticDissipationRate();
// }


bool Foam::mechanicalModel::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
