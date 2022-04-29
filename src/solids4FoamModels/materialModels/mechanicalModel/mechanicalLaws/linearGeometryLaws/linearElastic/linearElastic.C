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

#include "linearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearElastic, linGeomMechLaw
    );
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::linearElastic::makeSigma0f() const
{
    if (sigma0fPtr_)
    {
        FatalErrorIn("void Foam::linearElastic::makeSigma0f() const")
            << "pointer already set" << abort(FatalError);
    }

    sigma0fPtr_ =
        new surfaceSymmTensorField
        (
            "sigma0f",
            fvc::interpolate(sigma0_)
        );
}


const Foam::surfaceSymmTensorField& Foam::linearElastic::sigma0f() const
{
    if (!sigma0fPtr_)
    {
        makeSigma0f();
    }

    return *sigma0fPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElastic::linearElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu_("mu", dimPressure, 0.0),
    K_("K", dimPressure, 0.0),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    lambda_("lambda", dimPressure, 0.0),
    sigma0_
    (
        IOobject
        (
            "sigma0",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dict.lookupOrDefault<dimensionedSymmTensor>
        (
            "sigma0",
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
        )
    ),
    sigma0fPtr_(NULL)
{
    // Store old times
    epsilon().storeOldTime();
    epsilonf().storeOldTime();

    // Read elastic parameters
    // The user can specify E and nu or mu and K
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        E_ = dimensionedScalar(dict.lookup("E"));

        // Read the Poisson's ratio
        nu_ = dimensionedScalar(dict.lookup("nu"));

        // Set the shear modulus
        mu_ = E_/(2.0*(1.0 + nu_));

        // Set the bulk modulus
        if (nu_.value() < 0.5)
        {
            K_ = (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_.value() = GREAT;
        }
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        // Read shear modulus
        mu_ = dimensionedScalar(dict.lookup("mu"));

        // Read bulk modulus
        K_ = dimensionedScalar(dict.lookup("K"));

        // Calculate Young's modulus
        E_ = 9.0*K_*mu_/(3.0*K_ + mu_);

        // Calculate Poisson's ratio
        nu_ = (3.0*K_ - 2.0*mu_)/(2.0*(3.0*K_ + mu_));
    }
    else
    {
        FatalErrorIn
        (
            "linearElasticMisesPlastic::linearElasticMisesPlastic::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    // Set first Lame parameter
    if (nu_.value() < 0.5)
    {
        lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_));
    }
    else
    {
        lambda_.value() = GREAT;
    }

    // Check for physical Poisson's ratio
    if (nu_.value() < -1.0 || nu_.value() > 0.5)
    {
        FatalErrorIn
        (
            "Foam::linearElastic::linearElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Unphysical Poisson's ratio: nu should be >= -1.0 and <= 0.5"
            << abort(FatalError);
    }

    // Check for incompressibility or quasi-incompressibility
    if (nu_.value() > 0.49 && !solvePressureEqn())
    {
        WarningIn(type() + "::" + type())
            << "Poisson's ratio is greater than 0.49: "
            << "consider setting 'solvePressureEqn' to 'yes'!" << endl;
    }

    if (gMax(mag(sigma0_)()) > SMALL)
    {
        Info<< "Reading sigma0 initial/residual stress field" << endl;
    }

    if (gMax(mag(sigma0_)()) > SMALL)
    {
        Info<< "Reading sigma0 initial/residual stress field" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElastic::~linearElastic()
{
    deleteDemandDrivenData(sigma0fPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField> Foam::linearElastic::bulkModulus() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "bulkModulus",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            K_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

#ifdef OPENFOAMESIORFOUNDATION
    tresult.ref().correctBoundaryConditions();
#else
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::linearElastic::impK() const
{
    if (nu_.value() == 0.5)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "impK",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                2.0*mu_
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
                    "impK",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                2.0*mu_ + lambda_
            )
        );
    }

}


const Foam::dimensionedScalar& Foam::linearElastic::mu() const
{
    return mu_;
}


const Foam::dimensionedScalar& Foam::linearElastic::K() const
{
    return K_;
}


const Foam::dimensionedScalar& Foam::linearElastic::E() const
{
    return E_;
}


const Foam::dimensionedScalar& Foam::linearElastic::nu() const
{
    return nu_;
}


const Foam::dimensionedScalar& Foam::linearElastic::lambda() const
{
    return lambda_;
}


void Foam::linearElastic::correct(volSymmTensorField& sigma)
{
    // Update epsilon
    updateEpsilon();

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn
            (
                "void Foam::linearElasticMisesPlastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << abort(FatalError);
        }

        epsilon().replace
        (
            symmTensor::ZZ,
           -(nu_/E_)
           *(sigma.component(symmTensor::XX) + sigma.component(symmTensor::YY))
        );
    }

    // Hooke's law : standard form
    //sigma = 2.0*mu_*epsilon() + lambda_*tr(epsilon())*I + sigma0_;

    // Hooke's law : partitioned deviatoric and dilation form

    // Calculate hydrostatic stress
    updateSigmaHyd(K_*tr(epsilon()), 2*mu_ + lambda_);

    // Add deviatoric and initial stresses
    sigma = 2.0*mu_*dev(epsilon()) + sigmaHyd()*I + sigma0_;
}


void Foam::linearElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update epsilon
    updateEpsilonf();

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn
            (
                "void Foam::linearElasticMisesPlastic::"
                "correct(surfaceSymmTensorField& sigma)"
            )   << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << abort(FatalError);
        }

        epsilonf().replace
        (
            symmTensor::ZZ,
           -(nu_/E_)
           *(sigma.component(symmTensor::XX) + sigma.component(symmTensor::YY))
        );
    }

    // Hooke's law : standard form
    //sigma = 2.0*mu_*epsilonf() + lambda_*tr(epsilonf())*I + sigma0f();

    // Calculate hydrostatic stress
    surfaceScalarField* sigmaHydfPtr = NULL;
    surfaceScalarField& sigmaHydf = *sigmaHydfPtr;
    if (solvePressureEqn())
    {
        // Solve pressure equation at cells
        updateEpsilon();
        updateSigmaHyd(K_*tr(epsilon()), 2*mu_ + lambda_);

        // Interpolate to faces
        sigmaHydf = fvc::interpolate(sigmaHyd());
    }
    else
    {
        // Calculate hydrostatic stress at the faces
        sigmaHydf = K_*tr(epsilonf());
    }

    // Add deviatoric and initial stresses
    sigma = 2.0*mu_*dev(epsilonf()) + sigmaHydf*I + sigma0f();
}


// ************************************************************************* //
