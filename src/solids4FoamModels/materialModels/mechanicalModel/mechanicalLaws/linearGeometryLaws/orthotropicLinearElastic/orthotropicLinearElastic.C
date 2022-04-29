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

#include "orthotropicLinearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "fvc.H"
#include "doubleDotProduct.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(orthotropicLinearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, orthotropicLinearElastic, linGeomMechLaw
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::orthotropicLinearElastic::makeElasticC() const
{
    if (elasticCPtr_.valid())
    {
        FatalErrorIn
        (
            "void Foam::orthotropicLinearElastic::makeElasticC()"
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

#ifdef FOAMEXTEND
    elasticCPtr_.set
    (
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
        )
    );
    volSymmTensor4thOrderField& C = elasticCPtr_();

    // Calculate rotating matrix from local directions to global directions
    const volTensorField matDir =
        matDirX_*matDirX_ + matDirY_*matDirY_ + matDirZ_*matDirZ_;

    // Rotate C from local directions to global directions
    C = transform(matDir, C);
#else
    elasticCPtr_.set
    (
        new volTensorField
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
            dimensionedTensor
            (
                "localElasticC",
                dimPressure,
                tensor
                (
                    A11.value(), A12.value(), A31.value(),
                    A22.value(), A23.value(),
                    A33.value(),
                    A44.value(),
                    A55.value(),
                    A66.value()
                )
            )
        )
    );

    // Rotation not performed for OpenFOAM versions
#endif
}


#ifdef FOAMEXTEND
const Foam::volSymmTensor4thOrderField&
#else
const Foam::volTensorField&
#endif
Foam::orthotropicLinearElastic::elasticC() const
{
    if (elasticCPtr_.empty())
    {
        makeElasticC();
    }

    return elasticCPtr_();
}


void Foam::orthotropicLinearElastic::makeElasticCf() const
{
    if (elasticCfPtr_.valid())
    {
        FatalErrorIn
        (
            "void Foam::orthotropicLinearElastic::makeElasticCf()"
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

#ifdef FOAMEXTEND
    elasticCfPtr_.set
    (
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
        )
    );

    surfaceSymmTensor4thOrderField& C = elasticCfPtr_();

    // Calculate rotating matrix from local directions to global directions
    const volTensorField matDir =
        matDirX_*matDirX_ + matDirY_*matDirY_ + matDirZ_*matDirZ_;

    // Rotate C from local directions to global directions
    C = transform(fvc::interpolate(matDir), C);
#else
    elasticCfPtr_.set
    (
        new surfaceTensorField
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
            dimensionedTensor
            (
                "localElasticC",
                dimPressure,
                tensor
                (
                    A11.value(), A12.value(), A31.value(),
                    A22.value(), A23.value(),
                    A33.value(),
                    A44.value(),
                    A55.value(),
                    A66.value()
                )
            )
        )
    );

    // Rotation not performed for OpenFOAM versions
#endif
}


#ifdef FOAMEXTEND
const Foam::surfaceSymmTensor4thOrderField&
#else
const Foam::surfaceTensorField&
#endif
Foam::orthotropicLinearElastic::elasticCf() const
{
    if (elasticCfPtr_.empty())
    {
        makeElasticCf();
    }

    return elasticCfPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::orthotropicLinearElastic::orthotropicLinearElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
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
    elasticCPtr_(),
    elasticCfPtr_(),
    matDirX_
    (
        IOobject
        (
            "materialDirectionsX",
            mesh.time().timeName(),
            mesh,
#ifdef FOAMEXTEND
            IOobject::READ_IF_PRESENT,
#else
            IOobject::NO_READ, // Only global directions allowed
#endif
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
#ifdef FOAMEXTEND
            dict.lookupOrDefault<vector>("materialDirectionX", vector(1,0,0))
#else
            vector(1,0,0) // Only global directions allowed
#endif
        )
    ),
    matDirY_
    (
        IOobject
        (
           "materialDirectionsY",
            mesh.time().timeName(),
            mesh,
#ifdef FOAMEXTEND
            IOobject::READ_IF_PRESENT,
#else
            IOobject::NO_READ, // Only global directions allowed
#endif
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
#ifdef FOAMEXTEND
            dict.lookupOrDefault<vector>("materialDirectionX", vector(0,1,0))
#else
            vector(0,1,0) // Only global directions allowed
#endif
        )
    ),
    matDirZ_
    (
        IOobject
        (
           "materialDirectionsZ",
            mesh.time().timeName(),
            mesh,
#ifdef FOAMEXTEND
            IOobject::READ_IF_PRESENT,
#else
            IOobject::NO_READ, // Only global directions allowed
#endif
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
#ifdef FOAMEXTEND
            dict.lookupOrDefault<vector>("materialDirectionX", vector(0,0,1))
#else
            vector(0,0,1) // Only global directions allowed
#endif
        )
    )
{
    // Store old-time epsilon fields
    epsilon().storeOldTime();
    epsilonf().storeOldTime();

    if (planeStress())
    {
        FatalErrorIn
        (
            "Foam::orthotropicLinearElastic::orthotropicLinearElastic\n"
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
            "Foam::orthotropicLinearElastic::orthotropicLinearElastic\n"
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
            "Foam::orthotropicLinearElastic::orthotropicLinearElastic\n"
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
            "Foam::orthotropicLinearElastic::orthotropicLinearElastic\n"
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
    if
    (
        min(mag(matDirX_)).value() < SMALL
     || min(mag(matDirY_)).value() < SMALL
     || min(mag(matDirZ_)).value() < SMALL
    )
    {
        FatalErrorIn
        (
            "Foam::orthotropicLinearElastic::orthotropicLinearElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "The direction vectors must not have zero length!"
            << abort(FatalError);
    }

    // Normalise the direction vectors
    matDirX_ /= mag(matDirX_);
    matDirY_ /= mag(matDirY_);
    matDirZ_ /= mag(matDirZ_);

    // Check the material direcction vectors are locally orthogonal
    if
    (
        max(matDirX_ & matDirY_).value() > SMALL
     || max(matDirX_ & matDirZ_).value() > SMALL
     || max(matDirY_ & matDirZ_).value() > SMALL
    )
    {
        FatalErrorIn
        (
            "Foam::orthotropicLinearElastic::orthotropicLinearElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "The direction vectors should be locally orthogonal!"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::orthotropicLinearElastic::~orthotropicLinearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::orthotropicLinearElastic::impK() const
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
                elasticC().component(0) // symmTensor4thOrder::XXXX
              + elasticC().component(3) //symmTensor4thOrder::YYYY
              + elasticC().component(5) // symmTensor4thOrder::ZZZZ
            )/3.0
        )
    );
}


void Foam::orthotropicLinearElastic::correct(volSymmTensorField& sigma)
{
    // Update epsilon field
    updateEpsilon();

    // Calculate stress
#ifdef FOAMEXTEND
    sigma = elasticC() && epsilon();
#else
    doubleDotProduct(sigma, elasticC(), epsilon());
#endif
}


void Foam::orthotropicLinearElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update epsilonf field
    updateEpsilonf();

    // Calculate stress
#ifdef FOAMEXTEND
    sigma = elasticCf() && epsilonf();
#else
    doubleDotProduct(sigma, elasticCf(), epsilonf());
#endif
}


// ************************************************************************* //
