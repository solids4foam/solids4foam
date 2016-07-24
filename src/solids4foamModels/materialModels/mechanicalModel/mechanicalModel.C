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
#include "crackerFvMesh.H"

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
    // cohesiveDictPtr_(NULL),
    // cohesiveLawPtr_(NULL),
    planeStress_(lookup("planeStress"))
{
    Info<< "Creating mechanical model" << endl;

    // Create cohesiveLaw if the mesh is a crackerFvMesh
    // if (isA<crackerFvMesh>(mesh))
    // {
    //   Info<< "Reading cohesiveProperties because the mesh is a crackerFvMesh"
    //         << endl;

    //     cohesiveDictPtr_ =
    //         new IOdictionary
    //         (
    //             IOobject
    //             (
    //                 "cohesiveProperties",
    //                 mesh.time().constant(),
    //                 mesh,
    //                 IOobject::MUST_READ,
    //                 IOobject::NO_WRITE
    //             )
    //        );

    //     cohesiveLawPtr_ =
    //         cohesiveLaw::New
    //         (
    //             "law", mesh, cohesiveDictPtr_->subDict("cohesive")
    //         );
    // }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //


Foam::mechanicalModel::~mechanicalModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::mechanicalModel::viscoActive() const
{
    return mechanicalLawPtr_->viscoActive();
}


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


// const Foam::cohesiveLaw& Foam::mechanicalModel::cohLaw() const
// {
//     return cohesiveLawPtr_();
// }


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::rho() const
{
    return mechanicalLawPtr_->rho();
}

Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::rho(scalar t) const
{
    return mechanicalLawPtr_->rho(t);
}


void Foam::mechanicalModel::correct(Foam::volSymmTensorField& sigma)
{
    mechanicalLawPtr_->correct(sigma);
}


void Foam::mechanicalModel::correct
(
    Foam::volSymmTensorField& sigma, const int flag
)
{
    if (flag == 0)
    {
        mechanicalLawPtr_->correct(sigma);
    }
    else
    {
        mechanicalLawPtr_->correct(sigma, flag);
    }
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


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::mu() const
{
    const volScalarField lawE = mechanicalLawPtr_->E();
    const volScalarField lawNu = mechanicalLawPtr_->nu();

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            lawE/(2.0*(1.0 + lawNu))
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::lambda() const
{
    const volScalarField lawE = mechanicalLawPtr_->E();
    const volScalarField lawNu = mechanicalLawPtr_->nu();

    // Check if incompressible i.e. nu = 0.5
    if (mag(gMax(lawNu) - 0.5) < SMALL)
    {
        FatalErrorIn("tmp<volScalarField> mechanicalModel::lambda() const")
            << "At least one cell is incompressible: lambda is infinity!"
            << abort(FatalError);
    }

    if (planeStress())
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - lawNu))
            )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - 2.0*lawNu))
            )
        );
    }
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::threeK() const
{
    const volScalarField lawE = mechanicalLawPtr_->E();
    const volScalarField lawNu = mechanicalLawPtr_->nu();

    if (planeStress())
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "threeK",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawE/(1 - lawNu)
            )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "threeK",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawE/(1 - 2*lawNu)
            )
        );
    }
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::mu(Foam::scalar t) const
{
    const volScalarField lawE = mechanicalLawPtr_->E(t);
    const volScalarField lawNu = mechanicalLawPtr_->nu(t);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            lawE/(2.0*(1.0 + lawNu))
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::mechanicalModel::mu(const Foam::volScalarField& epsilonEq) const
{
    const volScalarField lawE = mechanicalLawPtr_->E(epsilonEq);
    const volScalarField lawNu = mechanicalLawPtr_->nu();

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            lawE/(2.0*(1.0 + lawNu))
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::mechanicalModel::lambda(Foam::scalar t) const
{
    const volScalarField lawE = mechanicalLawPtr_->E(t);
    const volScalarField lawNu = mechanicalLawPtr_->nu(t);

    if (planeStress())
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - lawNu))
            )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - 2.0*lawNu))
            )
        );
    }
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::lambda
(
    const Foam::volScalarField& epsilonEq
) const
{
    const volScalarField lawE = mechanicalLawPtr_->E(epsilonEq);
    const volScalarField lawNu = mechanicalLawPtr_->nu();

    if (planeStress())
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - lawNu))
            )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - 2.0*lawNu))
            )
        );
    }
}


Foam::tmp<Foam::volDiagTensorField> Foam::mechanicalModel::K() const
{
    return mechanicalLawPtr_->K();
}


Foam::tmp<Foam::volSymmTensor4thOrderField> Foam::mechanicalModel::C() const
{
    return mechanicalLawPtr_->C();
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::bulkModulus() const
{
    return mechanicalLawPtr_->bulkModulus();
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::viscosity() const
{
    return mechanicalLawPtr_->viscosity();
}


void Foam::mechanicalModel::updateYieldStress()
{
    mechanicalLawPtr_->updateYieldStress();
}


Foam::tmp<Foam::volScalarField>
Foam::mechanicalModel::plasticDissipationRate() const
{
    return mechanicalLawPtr_->plasticDissipationRate();
}


bool Foam::mechanicalModel::read()
{
    if (regIOobject::read()) // && cohesiveDictPtr_->regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
