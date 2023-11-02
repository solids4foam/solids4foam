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
    lambda_("lambda", dimPressure, 0.0)
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
        if (planeStress())
        {
            lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - nu_));

            if (solvePressureEqn())
            {
                FatalErrorIn
                (
                    "linearElasticMisesPlastic::linearElasticMisesPlastic::()"
                )   << "planeStress must be 'off' when solvePressureEqn is "
                    << "enabled" << abort(FatalError);
            }
        }
        else
        {
            lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_));
        }
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
            "    dictionary& dict\n"
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

    // Read the initial stress
    if (dict.found("sigma0"))
    {
        Info<< "Reading sigma0 from the dict" << endl;
        sigma0() = dimensionedSymmTensor(dict.lookup("sigma0"));
    }
    else if (gMax(mag(sigma0())()) > SMALL)
    {
        Info<< "Reading sigma0 stress field" << endl;
    }
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

#ifdef OPENFOAM_NOT_EXTEND
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


#ifdef OPENFOAM_NOT_EXTEND
Foam::scalarSquareMatrix
Foam::linearElastic::materialTangent() const
{
    scalarSquareMatrix matTang(6, 0.0);

    matTang(0,0) = 2*mu_.value() + lambda().value();
    matTang(0,1) = lambda().value();
    matTang(0,2) = lambda().value();

    matTang(1,0) = lambda().value();
    matTang(1,1) = 2*mu_.value() + lambda().value();
    matTang(1,2) = lambda().value();

    matTang(2,0) = lambda().value();
    matTang(2,1) = lambda().value();
    matTang(2,2) = 2*mu_.value() + lambda().value();

    matTang(3,3) = mu_.value();
    matTang(4,4) = mu_.value();
    matTang(5,5) = mu_.value();

    return matTang;
}
#endif


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

    if (solvePressureEqn())
    {
        // Calculate hydrostatic stress
        updateSigmaHyd(K_*tr(epsilon()), 2*mu_ + lambda_);

        // Hooke's law: partitioned deviatoric and dilation form
        sigma = 2.0*mu_*dev(epsilon()) + sigmaHyd()*I + sigma0();
    }
    else
    {
        // Hooke's law: standard form
        sigma = 2.0*mu_*epsilon() + lambda_*tr(epsilon())*I + sigma0();

        // Update sigmaHyd variable
        sigmaHyd() = -K_*tr(epsilon());
    }
}


void Foam::linearElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update epsilon
    updateEpsilonf();

    if (solvePressureEqn())
    {
        // Calculate hydrostatic stress at the cell-centres
        // Solve pressure equation at cells
        updateEpsilon();
        updateSigmaHyd(K_*tr(epsilon()), 2*mu_ + lambda_);

        // Interpolate to faces
        const surfaceScalarField sigmaHydf(fvc::interpolate(sigmaHyd()));

        // Add deviatoric and initial stresses
        sigma = 2.0*mu_*dev(epsilonf()) + sigmaHydf*I + sigma0f();
    }
    else
    {
        // Hooke's law : standard form
        sigma = 2.0*mu_*epsilonf() + lambda_*tr(epsilonf())*I + sigma0f();
    }
}


// ************************************************************************* //
