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

#include "EulerDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
#include "zeroGradientFvPatchFields.H"
#include "slipFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<typename EulerDdtScheme<vector>::fluxFieldType>
EulerDdtScheme<vector>::fvcDdtConsistentPhiCorr
(
    const GeometricField<vector, fvsPatchField, surfaceMesh>& Uf,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const surfaceScalarField& rAUf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

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

    if (mesh().moving())
    {
        V0oV.internalField() = mesh().V0()/mesh().V();
        V0oV.correctBoundaryConditions();
    }
    
    surfaceVectorField Uof = fvc::interpolate(U.oldTime()*V0oV);

    fluxFieldType phiUf0
    (
        (mesh().Sf() & Uf.oldTime())*
        fvc::interpolate(V0oV)
    );

    fluxFieldType phiCorr
    (
        phiUf0 - (mesh().Sf() & Uof)
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
            rDeltaT*phiCorr*rAUf
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
