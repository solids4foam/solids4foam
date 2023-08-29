/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "poroMechanicalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"
#include "demandDrivenData.H"

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

bool Foam::poroMechanicalLaw::checkSigmaEffReady(const volSymmTensorField& sigma, const volScalarField& p)
{
    if (sigmaEff_.valid())
    {
        return true;
    }
    
    sigmaEff_.set(
        new volSymmTensorField{
            IOobject
            (
                "sigmaEff",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
           sigma + b_*(p + p0_)*symmTensor(I)
        }
    );
    return true;
}

bool Foam::poroMechanicalLaw::checkSigmaEffReady(const surfaceSymmTensorField& sigma, const surfaceScalarField& p)
{
    if (sigmaEfff_.valid())
    {
        return true;
    }
    
    sigmaEfff_.set(
        new surfaceSymmTensorField{
            IOobject
            (
                "sigmaEfff",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
           sigma + b_*(p + p0f())*symmTensor(I)
        }
    );
    return true;
}

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
    sigmaEff_(),
    sigmaEfff_(),
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
    // Lookup the pressure field
    const volScalarField& p = lookupPressureField();

    // check if sigmaEff has been initialized (should be done only once per calculation)
    checkSigmaEffReady(sigma,p);

    // Calculate effective stress
    //-- Note that we could just pass "sigma" here but we use a separate field
    //-- called sigmaEff just for post-processing visualisation of the effective
    //-- stress <-- for stress state dependent material laws like Mohr-Coulomb it is important to use sigmaEff since the strength depends on tr(sigmaEff)
    effectiveStressMechLawPtr_->correct(sigmaEff_());

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure
    sigma = sigmaEff_() - b_*(p + p0_)*symmTensor(I);
}


void Foam::poroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    // Lookup the pressure field
    const volScalarField& p = lookupPressureField();

    // Interpolate pressure to the faces
    const surfaceScalarField pf(fvc::interpolate(p));

    // check if sigmaEff has been initialized (should be done only once per calculation)
    checkSigmaEffReady(sigma,pf);

    // Calculate effective stress
    effectiveStressMechLawPtr_->correct(sigmaEfff_());

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure
    sigma = sigmaEfff_() - b_*(pf + p0f())*symmTensor(I);
}


// ************************************************************************* //
