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

#include "poroMechanicalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poroMechanicalLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, poroMechanicalLaw, linGeomMechLaw
    );
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::poroMechanicalLaw::makeP0f() const
{
    if (p0fPtr_)
    {
        FatalErrorIn("void Foam::poroMechanicalLaw::makeP0f() const")
            << "pointer already set" << abort(FatalError);
    }

    p0fPtr_ =
        new surfaceScalarField
        (
            "p0f",
            fvc::interpolate(p0_)
        );
}


const Foam::surfaceScalarField& Foam::poroMechanicalLaw::p0f() const
{
    if (!p0fPtr_)
    {
        makeP0f();
    }

    return *p0fPtr_;
}


const Foam::volScalarField& Foam::poroMechanicalLaw::lookupPressureField() const
{
    if (mesh().thisDb().parent().foundObject<objectRegistry>(pRegion_))
    {
        return mesh().thisDb().parent().subRegistry
        (
            pRegion_
        ).lookupObject<volScalarField>(pName_);
    }
    else if
    (
        mesh().thisDb().parent().foundObject<objectRegistry>("solid")
    )
    {
        return mesh().thisDb().parent().subRegistry
        (
            "solid"
        ).lookupObject<volScalarField>(pName_);
    }
    else
    {
        FatalErrorIn("Foam::poroMechanicalLaw::lookupPressureField()")
            << "Cannot find " << pName_ << " field in " << pRegion_
            << " or in 'solid'" << abort(FatalError);
    }

    // Keep compiler happy
    return mesh().lookupObject<volScalarField>("null");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::poroMechanicalLaw::poroMechanicalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    effectiveStressMechLawPtr_
    (
        mechanicalLaw::NewLinGeomMechLaw
        (
            word(dict.subDict("effectiveStressMechanicalLaw").lookup("type")),
            mesh,
            dict.subDict("effectiveStressMechanicalLaw"),
            nonLinGeom
        )
    ),
    b_
    (
        dict.lookupOrDefault<dimensionedScalar>
        (
            "biotCoeff", dimensionedScalar("0", dimless, 1.0)
        )
    ),
    pName_(dict.lookupOrDefault<word>("pressureFieldName", "p")),
    pRegion_(dict.lookupOrDefault<word>("pressureFieldRegion", "region0")),
    p0_
    (
        IOobject
        (
            "p0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dict.lookupOrDefault<dimensionedScalar>
        (
            "p0",
            dimensionedScalar("zero", dimPressure, 0.0)
        )
    ),
    p0fPtr_(NULL)
{
    if (gMax(mag(p0_)()) > SMALL)
    {
        Info<< "Reading p0 initial/residual pore-pressure field" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::poroMechanicalLaw::~poroMechanicalLaw()
{
    deleteDemandDrivenData(p0fPtr_);
}


Foam::tmp<Foam::volScalarField> Foam::poroMechanicalLaw::impK() const
{
    return effectiveStressMechLawPtr_->impK();
}


void Foam::poroMechanicalLaw::correct(volSymmTensorField& sigma)
{
    // Calculate effective stress
    // Note that we could just pass "sigma" here but we use a separate field
    // called sigmaEff just for post-processing visualisation of the effective
    // stress
    effectiveStressMechLawPtr_->correct(sigmaEff());

    // Lookup the pressure field
    const volScalarField& p = lookupPressureField();

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure
    sigma = sigmaEff() - b_*(p + p0_)*symmTensor(I);
}


void Foam::poroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    // Calculate effective stress
    effectiveStressMechLawPtr_->correct(sigma);

    // Lookup the pressure field
    const volScalarField& p = lookupPressureField();

    // Interpolate pressure to the faces
    const surfaceScalarField pf(fvc::interpolate(p));

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure
    sigma -= b_*(pf + p0f())*symmTensor(I);
}


// ************************************************************************* //
