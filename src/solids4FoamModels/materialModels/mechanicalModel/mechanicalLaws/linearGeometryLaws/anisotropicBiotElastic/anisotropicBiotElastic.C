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

#include "anisotropicBiotElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(anisotropicBiotElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, anisotropicBiotElastic, linGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::anisotropicBiotElastic::anisotropicBiotElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    model2d_(bool(mesh.solutionD()[vector::Z] > 0)),
    rho_(dict.lookup("rho")),
    A11_(0.0),
    A22_(0.0),
    A33_(0.0),
    A44_(0.0),
    A55_(0.0),
    A66_(0.0),
    A12_(0.0),
    A21_(0.0),
    A31_(0.0),
    A23_(0.0),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    )
{
    // Set elastic stiffness parameters
    if (model2d_)
    {
        // Only the z direction is allow for 2-D models
        if (mesh.solutionD()[vector::X] < 0 || mesh.solutionD()[vector::Y] < 0)
        {
            FatalErrorIn(type() + "::" + type())
                << "For 2-D models, z must be the empty direction"
                << abort(FatalError);
        }

        const scalar Ex = readScalar(dict.lookup("Ex"));
        const scalar Ey = readScalar(dict.lookup("Ey"));
        const scalar vxy = readScalar(dict.lookup("nuxy"));
        const scalar Gxy = readScalar(dict.lookup("Gxy"));

        // Material contraints
        const scalar vyx = vxy*Ey/Ex;

        Info<< type() << "2D model"
            << nl
            << "Linear elastic orthotropic properties are:\n"
            << "Ex " << Ex/1e9 << " GPa" << nl
            << "Ey " << Ey/1e9 << " GPa" << nl
            << "nuxy " << vxy << nl
            << "Gxy " << Gxy/1e9 << " GPa" << endl;

        const scalar J = 1/(1 - vxy*vyx);
        A11_ = J*Ex;
        A22_ = J*Ey;
        A12_ = J*vyx*Ex;
        A21_ = J*vxy*Ey;
        A44_ = 2*Gxy;

        Info<< type() << ": 2D stiffness coefficients are:" << nl
            << "A11 " << A11_ << nl
            << "A22 " << A22_ << nl
            << "A12 " << A12_ << nl
            << "A21 " << A21_ << nl
            << "A44 " << A44_ << endl;
    }
    else // 3D
    {
        const scalar Ex = readScalar(dict.lookup("Ex"));
        const scalar Ey = readScalar(dict.lookup("Ey"));
        const scalar Ez = readScalar(dict.lookup("Ez"));
        const scalar vxy = readScalar(dict.lookup("nuxy"));
        const scalar vyz = readScalar(dict.lookup("nuyz"));
        const scalar vzx = readScalar(dict.lookup("nuzx"));
        const scalar Gxy = readScalar(dict.lookup("Gxy"));
        const scalar Gyz = readScalar(dict.lookup("Gyz"));
        const scalar Gzx = readScalar(dict.lookup("Gzx"));

        // Material contraints
        const scalar vyx = vxy*Ey/Ex;
        const scalar vxz = vzx*Ex/Ez;
        const scalar vzy = vyz*Ez/Ey;

        Info<< type() << ": 3D model" << nl
            << "Linear elastic orthotropic properties are:\n"
            << "Ex " << Ex/1e9 << " GPa" << nl
            << "Ey " << Ey/1e9 << " GPa" << nl
            << "Ez " << Ez/1e9 << " GPa" << nl
            << "nuxy " << vxy << nl
            << "nuyz " << vyz << nl
            << "nuzx " << vzx << nl
            << "Gxy " << Gxy/1e9 << " GPa" << nl
            << "Gyz " << Gyz/1e9 << " GPa" << nl
            << "Gzx " << Gzx/1e9 << " GPa" << endl;

        const scalar J =
            (1.0 - vxy*vyx - vyz*vzy - vzx*vxz - 2*vyx*vzy*vxz)/(Ex*Ey*Ez);
        A11_ = (1.0 - vyz*vzy)/(J*Ey*Ez);
        A22_ = (1.0 - vxz*vzx)/(J*Ex*Ez);
        A33_ = (1.0 - vyx*vxy)/(J*Ey*Ex);
        A12_ = (vxy + vzy*vxz)/(J*Ex*Ez);
        A31_ = (vzx + vyx*vzy)/(J*Ey*Ez);
        A23_ = (vyz + vyx*vxz)/(J*Ex*Ey);
        A44_ = 2*Gxy;
        A55_ = 2*Gyz;
        A66_ = 2*Gzx;

        Info<< type() << ": 3D stiffness coefficients:" << nl
            << "A11 " << A11_ << nl
            << "A22 " << A22_ << nl
            << "A33 " << A33_ << nl
            << "A12 " << A12_ << nl
            << "A31 " << A31_ << nl
            << "A23 " << A23_ << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::anisotropicBiotElastic::~anisotropicBiotElastic()
{}


// * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::anisotropicBiotElastic::rho() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
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


Foam::tmp<Foam::volScalarField> Foam::anisotropicBiotElastic::impK() const
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
            dimensionedScalar("impK", dimPressure, max(A11_, max(A22_, A33_)))
            // (A11_ + A22_ + A33_)/3.0 // could be better?
        )
    );
}

void Foam::anisotropicBiotElastic::correct(volSymmTensorField& sigma)
{
    // Calculate total strain
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        epsilon_ = epsilon_.oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        epsilon_ = symm(gradD);
    }

    // Calculate the effective stress
    // sigma = C:epsilon

    // Take references for convenience and efficiency
#ifdef OPENFOAMESIORFOUNDATION
    symmTensorField& sigmaI = sigma.primitiveFieldRef();
#else
    symmTensorField& sigmaI = sigma.internalField();
#endif
    const symmTensorField& epsilonI = epsilon_.internalField();

    // Internal field
    forAll(sigmaI, celli)
    {
        const scalar& e11 = epsilonI[celli][symmTensor::XX];
        const scalar& e22 = epsilonI[celli][symmTensor::YY];
        const scalar& e33 = epsilonI[celli][symmTensor::ZZ];
        const scalar& e12 = epsilonI[celli][symmTensor::XY];
        const scalar& e23 = epsilonI[celli][symmTensor::YZ];
        const scalar& e31 = epsilonI[celli][symmTensor::XZ];

        if (model2d_)
        {
            sigmaI[celli][symmTensor::XX] = A11_*e11 + A12_*e22;
            sigmaI[celli][symmTensor::YY] = A21_*e11 + A22_*e22;
            sigmaI[celli][symmTensor::XY] = A44_*e12;
        }
        else
        {
            sigmaI[celli][symmTensor::XX] = A11_*e11 + A12_*e22 + A31_*e33;
            sigmaI[celli][symmTensor::YY] = A12_*e11 + A22_*e22 + A23_*e33;
            sigmaI[celli][symmTensor::ZZ] = A31_*e11 + A23_*e22 + A33_*e33;
            sigmaI[celli][symmTensor::XY] = A44_*e12;
            sigmaI[celli][symmTensor::YZ] = A55_*e23;
            sigmaI[celli][symmTensor::XZ] = A66_*e31;
        }
    }

    // BOundary patches
    forAll(sigma.boundaryField(), patchI)
    {
        // Take references for convenience and efficiency
#ifdef OPENFOAMESIORFOUNDATION
        symmTensorField& sigmaP = sigma.boundaryFieldRef()[patchI];
#else
        symmTensorField& sigmaP = sigma.boundaryField()[patchI];
#endif
        const symmTensorField& epsilonP = epsilon_.boundaryField()[patchI];

        forAll(sigmaP, faceI)
        {
            const scalar& e11 = epsilonP[faceI][symmTensor::XX];
            const scalar& e22 = epsilonP[faceI][symmTensor::YY];
            const scalar& e33 = epsilonP[faceI][symmTensor::ZZ];
            const scalar& e12 = epsilonP[faceI][symmTensor::XY];
            const scalar& e23 = epsilonP[faceI][symmTensor::YZ];
            const scalar& e31 = epsilonP[faceI][symmTensor::XZ];
      
            if(model2d_)
            {
                sigmaP[faceI][symmTensor::XX] = A11_*e11 + A12_*e22;
                sigmaP[faceI][symmTensor::YY] = A21_*e11 + A22_*e22;
                sigmaP[faceI][symmTensor::XY] = A44_*e12;
            }
            else
            {
                sigmaP[faceI][symmTensor::XX] =
                    A11_*e11 + A12_*e22 + A31_*e33;
                sigmaP[faceI][symmTensor::YY] =
                    A12_*e11 + A22_*e22 + A23_*e33;
                sigmaP[faceI][symmTensor::ZZ] =
                    A31_*e11 + A23_*e22 + A33_*e33;
                sigmaP[faceI][symmTensor::XY] = A44_*e12;
                sigmaP[faceI][symmTensor::YZ] = A55_*e23;
                sigmaP[faceI][symmTensor::XZ] = A66_*e31;
            }
        }
    }
}


void Foam::anisotropicBiotElastic::correct(surfaceSymmTensorField& sigma)
{
    notImplemented
    (
        "void Foam::anisotropicBiotElastic::correct(surfaceSymmTensorField&)"
    );
}


// ************************************************************************* //
