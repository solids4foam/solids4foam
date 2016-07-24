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

#include "neoHookeanElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformGeometricField.H"
#include "fvc.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanElastic, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanElastic::neoHookeanElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mechanicalLaw(name, mesh, dict),
    restarted_
    (
        IOobject
        (
            "sigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        ).headerOk()
    ),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    K_
    (
        mesh.lookupObject<IOdictionary>
        (
            "mechanicalProperties"
        ).lookup("planeStress")
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    ),
    bEbar_
    (
        IOobject
        (
            "bEbar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    )
{
    // On restart, some fields may be defaulted to zero when they should default
    // to I
    if (gMin(mag(bEbar_.internalField())) < SMALL)
    {
        WarningIn("neoHookeanElastic::neoHookeanElastic()")
            << "Resseting zero in bEbar fields to I" << endl;

        symmTensorField& bEbarI = bEbar_.internalField();

        forAll(bEbarI, cellI)
        {
            if (mag(bEbarI[cellI]) < SMALL)
            {
                bEbarI[cellI] = 1.0*I;
            }
        }

        forAll(bEbar_.boundaryField(), patchI)
        {
            forAll(bEbar_.boundaryField()[patchI], faceI)
            {
                if (mag(bEbar_.boundaryField()[patchI][faceI]) < SMALL)
                {
                    bEbar_.boundaryField()[patchI][faceI] = 1.0*I;
                }
            }
        }

        bEbar_.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanElastic::~neoHookeanElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::E() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            E_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::nu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            nu_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::Ep() const
{
    notImplemented
    (
        "Foam::tmp<Foam::volScalarField> "
        "Foam::neoHookeanElastic::Ep() const"
    );

    // Keep compiler happy
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "undefined",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            0.0,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


void Foam::neoHookeanElastic::correct(volSymmTensorField& sigma)
{
    const fvMesh& mesh = this->mesh();

    // Lookup the total deformation gradient from the solver
    const volTensorField& F = mesh.lookupObject<volTensorField>("F");

    // Lookup the Jacobian of the deformation gradient from the solver
    const volScalarField& J = mesh.lookupObject<volScalarField>("J");

    // Update left Cauchy Green strain tensor with volumetric term removed
    bEbar_ = pow(J, -2.0/3.0)*symm(F & F.T());

    // Calculate deviatoric stress
    volSymmTensorField s = mu_*dev(bEbar_);

    // Calculate the Cauchy stress
    sigma = (1.0/J)*0.5*K_*(pow(J, 2) - 1)*I + s;
}


void Foam::neoHookeanElastic::correct(volSymmTensorField& tau, const int flag)
{
    correct(tau);
}


void Foam::neoHookeanElastic::correct(surfaceSymmTensorField& tau)
{
    notImplemented("neoHookeanElastic::correct(surfaceSymmTensorField& tau)");

    // const fvMesh& mesh = this->mesh();

    // // Lookup relative deformation gradient from the solver
    // const surfaceTensorField& relFbar =
    //     mesh.lookupObject<surfaceTensorField>("relFbarf");

    // // Update left Cauchy Green strain tensor with volumetric term removed
    // bEbarf_ = transform(relFbar, bEbarf_.oldTime());

    // // Lookup the Jacobian of the deformation gradient from the solver
    // const surfaceScalarField& J =
    // mesh.lookupObject<surfaceScalarField>("Jf");

    // // Calculate deviatoric stress
    // surfaceSymmTensorField s = mu_*dev(bEbarf_);

    // // Calculate new Kirchhoff stress

    // surfaceSymmTensorField newTau = 0.5*K_*(pow(J, 2) - 1)*I + s;

    // // Add thermal component, if active
    // if (mesh.foundObject<thermalModel>("thermalProperties"))
    // {
    //     const surfaceScalarField& threeKalpha =
    //         mesh.lookupObject<surfaceScalarField>("threeKalphaf");

    //     const surfaceScalarField& DT =
    //         mesh.lookupObject<surfaceScalarField>("DTf");

    //     newTau += threeKalpha*DT*symmTensor(I);
    // }

    // // Assign Kirchhoff stress
    // // For now, to deal with multi-materials, we will multiply by curMaterial
    // // index field so only cells in the current material are calculated:
    // // we should be able to do this in a better way

    // const surfaceScalarField curMat =
    //     pos(fvc::interpolate(curMaterial()) - SMALL);

    // tau = curMat*newTau + (1.0 - curMat)*tau;
}


Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElastic::plasticDissipationRate() const
{
    // Elastic material dissipates zero plastic energy
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "plasticDissipationRate",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimForce/(dimArea*dimTime), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


void Foam::neoHookeanElastic::setMaterialIndex(label curMatIndex)
{
    // Set current material index
    curMaterialIndex() = curMatIndex;

    // Rename fields to avoid conflicts with other mechanical laws
    curMaterial().rename(curMaterial().name() + '_' + Foam::name(curMatIndex));
    bEbar_.rename("bEbar_mat" + Foam::name(curMatIndex));

    // Temporary fix: fields are called sigmaY_mat0 for example but we read in
    // sigmaY in the constructor which does not exist; so we need to read them
    // here now that we know the material ID.
    // A better way would be the materials knows its ID from the constructor.

    const fvMesh& mesh = this->mesh();

    // Check for restarted case (if sigmaY_matID exists)
    if
    (
        IOobject
        (
            "bEbar_mat" + Foam::name(curMatIndex),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        ).headerOk()
    )
    {
        Info<< "Restarted case: correcting fields" << endl;

        bEbar_ =
            volSymmTensorField
            (
                IOobject
                (
                    "bEbar_mat" + Foam::name(curMatIndex),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("I", dimless, I)
            );

        // Reset cells in new roller regions which may have fields incorrectly
        // initialised to zero

        symmTensorField& bEbarI = bEbar_.internalField();

        forAll(bEbarI, cellI)
        {
            if (mag(bEbarI[cellI]) < SMALL)
            {
                bEbarI[cellI] = I;
            }
        }

        // Check all boundary values
        forAll(bEbar_.boundaryField(), patchI)
        {
            symmTensorField& pBEbar = bEbar_.boundaryField()[patchI];

            forAll(pBEbar, faceI)
            {
                if (mag(pBEbar[faceI]) < SMALL)
                {
                    pBEbar[faceI] = I;
                }
            }
        }

        bEbar_.correctBoundaryConditions();
    }
}


// ************************************************************************* //
