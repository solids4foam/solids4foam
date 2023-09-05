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

#include "StVenantKirchhoffOrthotropicElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "fvc.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(StVenantKirchhoffOrthotropicElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, StVenantKirchhoffOrthotropicElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


void Foam::StVenantKirchhoffOrthotropicElastic::makeElasticC() const
{
    if (elasticCPtr_)
    {
        FatalErrorIn
        (
            "void Foam::StVenantKirchhoffOrthotropicElastic::makeElasticC()"
        )   << "pointer already set" << abort(FatalError);
    }

    // Set local fourth order tensor stiffness

    const dimensionedScalar J =
        (
            1.0 - nu12_*nu21_ - nu23_*nu32_ - nu31_*nu13_ - 2*nu21_*nu32_*nu13_
        )/(E1_*E2_*E3_);

    const dimensionedScalar A11 = (1.0 - nu23_*nu32_)/(J*E2_*E3_);
    const dimensionedScalar A22 = (1.0 - nu13_*nu31_)/(J*E1_*E3_);
    const dimensionedScalar A33 = (1.0 - nu21_*nu12_)/(J*E2_*E1_);
    const dimensionedScalar A12 = (nu12_ + nu32_*nu13_)/(J*E1_*E3_);
    const dimensionedScalar A31 = (nu31_ + nu21_*nu32_)/(J*E2_*E3_);
    const dimensionedScalar A23 = (nu23_ + nu21_*nu13_)/(J*E1_*E2_);
    const dimensionedScalar A44 =  2.0*G12_;
    const dimensionedScalar A55 =  2.0*G23_;
    const dimensionedScalar A66 =  2.0*G31_;

    // Set elasticC in local coordinate system and then we will rotate it to the
    // global coordinate system

    elasticCPtr_ =
        new volSymmTensor4thOrderField
        (
            IOobject
            (
                "elasticC",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor4thOrder
            (
                "localElasticC",
                dimPressure,
                symmTensor4thOrder
                (
                    A11.value(), A12.value(), A31.value(),
                    A22.value(), A23.value(),
                    A33.value(),
                    A44.value(),
                    A55.value(),
                    A66.value()
                )
            )
        );

    volSymmTensor4thOrderField& C = *elasticCPtr_;

    C = transform(matDir_, C);
}


const Foam::volSymmTensor4thOrderField&
Foam::StVenantKirchhoffOrthotropicElastic::elasticC() const
{
    if (!elasticCPtr_)
    {
        makeElasticC();
    }

    return *elasticCPtr_;
}


void Foam::StVenantKirchhoffOrthotropicElastic::makeElasticCf() const
{
    if (elasticCfPtr_)
    {
        FatalErrorIn
        (
            "void Foam::StVenantKirchhoffOrthotropicElastic::makeElasticCf()"
        )   << "pointer already set" << abort(FatalError);
    }

    // Set local fourth order tensor stiffness

    const dimensionedScalar J =
        (
            1.0 - nu12_*nu21_ - nu23_*nu32_ - nu31_*nu13_ - 2*nu21_*nu32_*nu13_
        )/(E1_*E2_*E3_);

    const dimensionedScalar A11 = (1.0 - nu23_*nu32_)/(J*E2_*E3_);
    const dimensionedScalar A22 = (1.0 - nu13_*nu31_)/(J*E1_*E3_);
    const dimensionedScalar A33 = (1.0 - nu21_*nu12_)/(J*E2_*E1_);
    const dimensionedScalar A12 = (nu12_ + nu32_*nu13_)/(J*E1_*E3_);
    const dimensionedScalar A31 = (nu31_ + nu21_*nu32_)/(J*E2_*E3_);
    const dimensionedScalar A23 = (nu23_ + nu21_*nu13_)/(J*E1_*E2_);
    const dimensionedScalar A44 =  2.0*G12_;
    const dimensionedScalar A55 =  2.0*G23_;
    const dimensionedScalar A66 =  2.0*G31_;

    // Set elasticC in local coordinate system and then we will rotate it to the
    // global coordinate system

    elasticCfPtr_ =
        new surfaceSymmTensor4thOrderField
        (
            IOobject
            (
                "elasticCf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor4thOrder
            (
                "localElasticC",
                dimPressure,
                symmTensor4thOrder
                (
                    A11.value(), A12.value(), A31.value(),
                    A22.value(), A23.value(),
                    A33.value(),
                    A44.value(),
                    A55.value(),
                    A66.value()
                )
            )
        );

    surfaceSymmTensor4thOrderField& C = *elasticCfPtr_;

    C = transform(fvc::interpolate(matDir_), C);
}


const Foam::surfaceSymmTensor4thOrderField&
Foam::StVenantKirchhoffOrthotropicElastic::elasticCf() const
{
    if (!elasticCfPtr_)
    {
        makeElasticCf();
    }

    return *elasticCfPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::StVenantKirchhoffOrthotropicElastic::StVenantKirchhoffOrthotropicElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(typeName, name, mesh, dict, nonLinGeom),
    E1_(dict.lookup("E1")),
    E2_(dict.lookup("E2")),
    E3_(dict.lookup("E3")),
    nu12_(dict.lookup("nu12")),
    nu23_(dict.lookup("nu23")),
    nu31_(dict.lookup("nu31")),
    nu21_(nu12_*E2_/E1_),
    nu32_(nu23_*E3_/E2_),
    nu13_(nu31_*E1_/E3_),
    G12_(dict.lookup("G12")),
    G23_(dict.lookup("G23")),
    G31_(dict.lookup("G31")),
    mu_("mu", dimPressure, 0.0),
    K_("K", dimPressure, 0.0),
    elasticCPtr_(NULL),
    elasticCfPtr_(NULL),
    matDir_
    (
        IOobject
        (
            "materialDirections",
             mesh.time().timeName(),
             mesh,
             IOobject::READ_IF_PRESENT,
             IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor
        (
            "localMaterialDirections",
            dimless,
            tensor
            (
                vector(1,0,0)
               *vector
                (
                    dict.lookupOrAddDefault<vector>
                    (
                        "materialDirection1", vector(1,0,0)
                    )
                )
              + vector(0,1,0)
               *vector
                (
                    dict.lookupOrAddDefault<vector>
                    (
                        "materialDirection2", vector(0,1,0)
                    )
                )
              + vector(0,0,1)
               *vector
                (
                    dict.lookupOrAddDefault<vector>
                    (
                        "materialDirection2", vector(0,0,1)
                    )
                )
            )
        )
    )
{
    if (planeStress())
    {
        FatalErrorIn
        (
            "Foam::StVenantKirchhoffOrthotropicElastic::"
            "StVenantKirchhoffOrthotropicElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Material law not implemented for planeStress!"
            << abort(FatalError);
    }

    // Check material properties are physical
    Info<< "Checking physical constraints on the orthotropic material "
        << "properties" << endl;

    // E and G should be greater than zero
    if
    (
        E1_.value() < 0.0 || E2_.value() < 0.0 || E3_.value() < 0.0
     || G12_.value() < 0.0 || G23_.value() < 0.0 || G31_.value() < 0.0
    )
    {
        FatalErrorIn
        (
            "Foam::StVenantKirchhoffOrthotropicElastic::"
            "StVenantKirchhoffOrthotropicElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "E1, E2, E3, G12, G23, G31 should all be greater than zero!"
            << abort(FatalError);
    }

    // Restriction on Poisson's ratio
    if
    (
        mag(nu12_.value()) >= sqrt(E1_.value()/E2_.value())
     || mag(nu23_.value()) >= sqrt(E2_.value()/E3_.value())
     || mag(nu31_.value()) >= sqrt(E3_.value()/E1_.value())
    )
    {
        FatalErrorIn
        (
            "Foam::StVenantKirchhoffOrthotropicElastic::"
            "StVenantKirchhoffOrthotropicElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Unphysical Poisson's ratio!"
            << " mag(nu_ij) should be less sqrt(E_i/E_j)"
            << abort(FatalError);
    }

    if
    (
        1.0 - nu12_.value()*nu21_.value() - nu23_.value()*nu32_.value()
      - nu31_.value()*nu13_.value()
      - 2.0*nu21_.value()*nu32_.value()*nu13_.value()
     <= 0.0
    )
    {
        FatalErrorIn
        (
            "Foam::StVenantKirchhoffOrthotropicElastic::"
            "StVenantKirchhoffOrthotropicElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Unphysical Poisson's ratio!"
            << " (1 - nu12*nu21 - nu23*nu32 "
            << "- nu31*nu13 - 2*nu21*nu32*nu13) should be greater than zero!"
            << abort(FatalError);
    }

    // Check the direction vectors
    if (min(mag(matDir_)).value() < SMALL)
    {
        FatalErrorIn
        (
            "Foam::StVenantKirchhoffOrthotropicElastic::"
            "StVenantKirchhoffOrthotropicElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "The direction vectors must not have zero length!"
            << abort(FatalError);
    }

    // Normalise the direction vectors
    matDir_ /= mag(matDir_);

    // Initialise mu and K
    mu_ = (G12_ + G23_ + G31_)/3.0;
    const dimensionedScalar averageE = (E1_ + E2_ + E3_)/3.0;
    K_ = averageE*mu_/(3*(3*mu_ - averageE));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::StVenantKirchhoffOrthotropicElastic::
~StVenantKirchhoffOrthotropicElastic()
{
    deleteDemandDrivenData(elasticCPtr_);
    deleteDemandDrivenData(elasticCfPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::
StVenantKirchhoffOrthotropicElastic::impK() const
{
    // We could employ a diagTensor as the implicit stiffness; however, the
    // average of this diagTensor should achieve similar convergence
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
            (
                elasticC().component(symmTensor4thOrder::XXXX)
              + elasticC().component(symmTensor4thOrder::YYYY)
              + elasticC().component(symmTensor4thOrder::ZZZZ)
            )/3.0
        )
    );
}


void Foam::StVenantKirchhoffOrthotropicElastic::correct
(
    volSymmTensorField& sigma
)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the right Cauchy–Green deformation tensor
    const volSymmTensorField c = symm(F().T() & F());

    // Calculate the Green strain tensor
    const volSymmTensorField E = 0.5*(c - I);

    // Calculate the 2nd Piola Kirchhoff stress
    const volSymmTensorField S = elasticC() && E;

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F()));

    // Convert the 2nd Piola Kirchhoff stress to the Cauchy stress
    // sigma = (1.0/J)*symm(F() & S & F().T());
    sigma = (1.0/J)*transform(F(), S);
}


void Foam::StVenantKirchhoffOrthotropicElastic::correct
(
    surfaceSymmTensorField& sigma
)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the right Cauchy–Green deformation tensor
    const surfaceSymmTensorField c = symm(Ff().T() & Ff());

    // Calculate the Green strain tensor
    const surfaceSymmTensorField E = 0.5*(c - I);

    // Calculate the 2nd Piola Kirchhoff stress
    const surfaceSymmTensorField S = elasticCf() && E;

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(Ff()));

    // Convert the 2nd Piola Kirchhoff stress to the Cauchy stress
    // sigma = (1.0/J)*symm(Ff() & S & Ff().T());
    sigma = (1.0/J)*transform(Ff(), S);
}


void Foam::StVenantKirchhoffOrthotropicElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}

// ************************************************************************* //
