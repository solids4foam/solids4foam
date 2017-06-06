/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
    
\*---------------------------------------------------------------------------*/

#include "EulerDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<surfaceScalarField> EulerDdtScheme<vector>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const surfaceScalarField& phi
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    surfaceScalarField ddtPhiCoeff
    (
        IOobject
        (
            "ddtPhiCoeff",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensioned<scalar>("1", dimless, 1.0)
    );

    forAll(U.boundaryField(), patchI)
    {
        if
        (
            U.boundaryField()[patchI].fixesValue()
         || isA<symmetryFvPatchVectorField>(U.boundaryField()[patchI])
         || isA<slipFvPatchVectorField>(U.boundaryField()[patchI])
        )
        {
            ddtPhiCoeff.boundaryField()[patchI] = 0.0;
        }
    }

    if (mesh().moving())
    {
        volScalarField V0oV
        (
            IOobject
            (
                "V0oV",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );

        V0oV.internalField() = mesh().V0()/mesh().V();
        V0oV.correctBoundaryConditions();

        const surfaceVectorField& Sf = 
            mesh().objectRegistry::lookupObject<surfaceVectorField>("Sf");

        // Non-conservative cell-face velocity
        surfaceVectorField U0 = fvc::interpolate(V0oV*U.oldTime());
        forAll(U0.boundaryField(), patchI)
        {
            if (!U.boundaryField()[patchI].coupled())
            {
                U0.boundaryField()[patchI] = 
                    U.oldTime().boundaryField()[patchI]
                   .patchInternalField()
                   *V0oV.boundaryField()[patchI];
            }
        }

        // Conservataive cell-face velocity
        surfaceVectorField U0c = fvc::interpolate(U.oldTime());
        U0c -= (Sf.oldTime()&U0c)*Sf.oldTime()/magSqr(Sf.oldTime());
        U0c += phi.oldTime()*Sf.oldTime()/magSqr(Sf.oldTime());

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                ddtIOobject,
                rDeltaT*ddtPhiCoeff
               *((fvc::interpolate(V0oV)*U0c - U0) & mesh().Sf())
               /fvc::interpolate(1/rA)
            )
        );
    }
    else
    {
        // Non-conservative cell-face velocity
        surfaceVectorField U0 = fvc::interpolate(U.oldTime());
        forAll(U0.boundaryField(), patchI)
        {
            if (!U.boundaryField()[patchI].coupled())
            {
                U0.boundaryField()[patchI] = 
                    U.oldTime().boundaryField()[patchI]
                   .patchInternalField();
            }
        }

        // Conservataive cell-face velocity
        surfaceVectorField U0c = fvc::interpolate(U.oldTime());
        U0c -= (Sf.oldTime()&U0)*Sf.oldTime()/magSqr(Sf.oldTime());
        U0c += phi.oldTime()*Sf.oldTime()/magSqr(Sf.oldTime());

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                ddtIOobject,
                rDeltaT*ddtPhiCoeff
               *((U0c - U0) & mesh().Sf())
               /fvc::interpolate(1/rA)
            )
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
