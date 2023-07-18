/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "backwardDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
#include "slipFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<typename backwardDdtScheme<vector>::fluxFieldType>
backwardDdtScheme<vector>::fvcDdtConsistentPhiCorr
(
    const GeometricField<vector, fvsPatchField, surfaceMesh>& Uf,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const surfaceScalarField& rAUf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(U);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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
        dimensionedScalar("1", dimless, 1),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField V00oV
    (
        IOobject
        (
            "V00oV",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("1", dimless, 1),
        zeroGradientFvPatchScalarField::typeName
    );

    if (mesh().moving())
    {
        V0oV.internalField() = mesh().V0()/mesh().V();
        V0oV.correctBoundaryConditions();
    
        V00oV.internalField() = mesh().V00()/mesh().V();
        V00oV.correctBoundaryConditions();
    }

    surfaceVectorField Uof = fvc::interpolate(U.oldTime()*V0oV);
    surfaceVectorField Uoof = fvc::interpolate(U.oldTime().oldTime()*V00oV);

    fluxFieldType phiCorr
    (
        mesh().Sf()
      & (
            (
                coefft0*Uf.oldTime()*fvc::interpolate(V0oV)
              - coefft00*Uf.oldTime().oldTime()*fvc::interpolate(V00oV)
            )
          - (
                coefft0*Uof
              - coefft00*Uoof
            )
        )
    );

    forAll(phiCorr.boundaryField(), patchI)
    {
        if
        (
            U.boundaryField()[patchI].fixesValue()
         || isA<slipFvPatchVectorField>(U.boundaryField()[patchI])
        )
        {
            phiCorr.boundaryField()[patchI] = 0;
        }
    }

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
                mesh().time().timeName(),
                mesh()
            ),
            // this->fvcDdtPhiCoeff(U.oldTime(), (mesh().Sf() & Uf.oldTime()))
            rDeltaT*phiCorr*rAUf
        )
    );
}


template<>
tmp<surfaceScalarField> backwardDdtScheme<vector>::meshPhi
(
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Coefficient for t-3/2 (between times 0 and 00)
    scalar coefft0_00 = deltaT/(deltaT + deltaT0);

    // Coefficient for t-1/2 (between times n and 0)
    scalar coefftn_0 = 1 + coefft0_00;

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                mesh().phi().name(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            coefftn_0*mesh().phi() - coefft0_00*mesh().phi().oldTime()
        )
    );
    
    // ZT: This bugfix is incorrect. ESI version is correct!
    // Bugfix: missing possibility of having the variable time step
    // Reported by Sopheak Seng, Bureau Veritas, 6/Sep/2018.
    // const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    // const scalar coefft  = 1 + deltaT/(deltaT + deltaT0);

    // return coefft*mesh().phi() - coefft00*mesh().phi().oldTime();
}

    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
