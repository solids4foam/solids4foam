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

#include "rotationInvariantLinearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "eig3Field.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rotationInvariantLinearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, rotationInvariantLinearElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::rotationInvariantLinearElastic::rotationInvariantLinearElastic
(
    const word& name,
    const fvMesh& mesh,
    dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    biotStrain_
    (
        IOobject
        (
            "biotStrain_" + name,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("0", dimless, symmTensor::zero)
    ),
    mu_("mu", dimPressure, 0.0),
    K_("K", dimPressure, 0.0),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    lambda_("lambda", dimPressure, 0.0)
{
    // Read mechanical properties
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
        FatalErrorIn(type())
            << "Either E and nu or mu and K elastic parameters should be "
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
        notImplemented
        (
            "Not implemented for incompressible behaviour: "
            "try a smaller value of Poisson's ratio, nu"
        );
    }

    // Store old F
    F().storeOldTime();
    Ff().storeOldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rotationInvariantLinearElastic::~rotationInvariantLinearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::rotationInvariantLinearElastic::impK() const
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
            (4.0/3.0)*mu_ + K_ // == 2*mu + lambda
        )
    );
}


void Foam::rotationInvariantLinearElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F()));

    // Calculate the right Cauchy Green tensor
    const volSymmTensorField C(symm(F().T() & F()));

    // Eigen value field of C
    volVectorField lambda
    (
        IOobject
        (
            "eigenVal(" + C.name() + ")",
            C.time().timeName(),
            C.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C.mesh(),
        dimensionedVector("zero", C.dimensions(), vector::zero)
    );

    // Eigen vectors will be store in the rows i.e. the first eigen vector
    // is (eigenVec.xx() eigenVec.xy() eigenVec.xz())
    volTensorField eigVec
    (
        IOobject
        (
            "eigenVec(" + C.name() + ")",
            C.time().timeName(),
            C.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C.mesh(),
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    // Calculate eigen values and eigen vectors of C
    eig3Field(C, eigVec, lambda);

    // Create a volVectorField for each principal direction
    const vector i(1, 0, 0);
    const vector j(0, 1, 0);
    const vector k(0, 0, 1);
    const volVectorField eigVec1
    (
        "eigVec1",
        eigVec.component(tensor::XX)*i
      + eigVec.component(tensor::XY)*j
      + eigVec.component(tensor::XZ)*k
    );
    const volVectorField eigVec2
    (
        "eigVec2",
        eigVec.component(tensor::YX)*i
      + eigVec.component(tensor::YY)*j
      + eigVec.component(tensor::YZ)*k
    );
    const volVectorField eigVec3
    (
        "eigVec3",
        eigVec.component(tensor::ZX)*i
      + eigVec.component(tensor::ZY)*j
      + eigVec.component(tensor::ZZ)*k
    );

    // Calculate the Biot strain
    biotStrain_ =
        sqrt(lambda.component(vector::X))*sqr(eigVec1)
      + sqrt(lambda.component(vector::Y))*sqr(eigVec2)
      + sqrt(lambda.component(vector::Z))*sqr(eigVec3)
      - I;

    // Calculate the 2nd Piola-Kirchhoff stress using Hooke's law
    const volSymmTensorField S
    (
        "S2PK", 2.0*mu_*biotStrain_ + lambda_*tr(biotStrain_)*I
    );

    // Calculate the Cauchy stress
    sigma = J*symm(F() & S & F().T());
}


void Foam::rotationInvariantLinearElastic::correct(surfaceSymmTensorField& sigma)
{
    notImplemented
    (
        "void Foam::rotationInvariantLinearElastic::correct"
        "(surfaceSymmTensorField& sigma)"
    );

    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(Ff()));

    // Calculate the right Cauchy Green tensor
    const surfaceSymmTensorField C(symm(Ff().T() & Ff()));

    // Eigen value field of C
    surfaceVectorField lambda
    (
        IOobject
        (
            "eigenVal(" + C.name() + ")",
            C.time().timeName(),
            C.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C.mesh(),
        dimensionedVector("zero", C.dimensions(), vector::zero)
    );

    // Eigen vectors will be store in the rows i.e. the first eigen vector
    // is (eigenVec.xx() eigenVec.xy() eigenVec.xz())
    surfaceTensorField eigVec
    (
        IOobject
        (
            "eigenVec(" + C.name() + ")",
            C.time().timeName(),
            C.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C.mesh(),
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    // Calculate eigen values and eigen vectors of C
    eig3Field(C, eigVec, lambda);

    // Create a surfaceVectorField for each principal direction
    const vector i(1, 0, 0);
    const vector j(0, 1, 0);
    const vector k(0, 0, 1);
    const surfaceVectorField eigVec1
    (
        "eigVec1",
        eigVec.component(tensor::XX)*i
      + eigVec.component(tensor::XY)*j
      + eigVec.component(tensor::XZ)*k
    );
    const surfaceVectorField eigVec2
    (
        "eigVec2",
        eigVec.component(tensor::YX)*i
      + eigVec.component(tensor::YY)*j
      + eigVec.component(tensor::YZ)*k
    );
    const surfaceVectorField eigVec3
    (
        "eigVec3",
        eigVec.component(tensor::ZX)*i
      + eigVec.component(tensor::ZY)*j
      + eigVec.component(tensor::ZZ)*k
    );

    // Calculate the Biot strain
    const surfaceSymmTensorField biotStrain
    (
        "biotStrain",
        sqrt(lambda.component(vector::X))*sqr(eigVec1)
      + sqrt(lambda.component(vector::Y))*sqr(eigVec2)
      + sqrt(lambda.component(vector::Z))*sqr(eigVec3)
      - I
    );

    // Calculate the 2nd Piola-Kirchhoff stress using Hooke's law
    const surfaceSymmTensorField S
    (
        "S2PK", 2.0*mu_*biotStrain + lambda_*tr(biotStrain)*I
    );

    // Calculate the Cauchy stress
    sigma = J*symm(Ff() & S & Ff().T());
}


void Foam::rotationInvariantLinearElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}

// ************************************************************************* //
