/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "wedgeFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<GeometricField<vector, fvPatchField, volMesh> >
EulerDdtScheme<vector>::fvcDdt
(
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        if
        (
            mesh().objectRegistry::found("grad(" + vf.name() + ")")
         && mesh().objectRegistry::found("meshU")
        )
        {
            const volTensorField& gradVf =
                mesh().objectRegistry::lookupObject<volTensorField>
                (
                    "grad(" + vf.name() + ")"
                );

            const volVectorField& meshU =
                mesh().objectRegistry::lookupObject<volVectorField>
                (
                    "meshU"
                );

            return tmp<GeometricField<vector, fvPatchField, volMesh> >
            (
                new GeometricField<vector, fvPatchField, volMesh>
                (
                    ddtIOobject,
                    rDeltaT*(vf - vf.oldTime()) - (meshU&gradVf.oldTime())
                )
            );
        }
        else
        {
            return tmp<GeometricField<vector, fvPatchField, volMesh> >
            (
                new GeometricField<vector, fvPatchField, volMesh>
                (
                    ddtIOobject,
                    rDeltaT*(vf - vf.oldTime())
                )
            );
        }

//         return tmp<GeometricField<vector, fvPatchField, volMesh> >
//         (
//             new GeometricField<vector, fvPatchField, volMesh>
//             (
//                 ddtIOobject,
//                 mesh(),
//                 rDeltaT.dimensions()*vf.dimensions(),
//                 rDeltaT.value()*
//                 (
//                     vf.internalField()
//                   - vf.oldTime().internalField()
//                 ),
//                 rDeltaT.value()*
//                 (
//                     vf.boundaryField() - vf.oldTime().boundaryField()
//                 )
//             )
//         );
    }
    else
    {
        return tmp<GeometricField<vector, fvPatchField, volMesh> >
        (
            new GeometricField<vector, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


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
         || isA<wedgeFvPatchVectorField>(U.boundaryField()[patchI])
        )
        {
            ddtPhiCoeff.boundaryField()[patchI] = 0.0;
        }
    }

    if (mesh().moving())
    {
        Info << "Moving mesh ddtPhiCorr " << U.name() << endl;

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
        surfaceVectorField U0 =
            fvc::interpolate(V0oV*U.oldTime(), "interpolate(U)");
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
        surfaceVectorField U0c =
            fvc::interpolate(U.oldTime(), "interpolate(U)");
        U0c -= (Sf.oldTime()&U0c)*Sf.oldTime()/magSqr(Sf.oldTime());
        U0c += phi.oldTime()*Sf.oldTime()/magSqr(Sf.oldTime());

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                ddtIOobject,
                rDeltaT*ddtPhiCoeff
                *((fvc::interpolate(V0oV)*U0c - U0) & mesh().Sf())
               /fvc::interpolate(1/rA, "interpolate(U)")
            )
        );
    }
    else
    {
        Info << "Fixed mesh ddtPhiCorr " << U.name() << endl;

        // Non-conservative cell-face velocity
        surfaceVectorField U0 =
            fvc::interpolate(U.oldTime(), "interpolate(U)");
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
        surfaceVectorField U0c = U0;
//             fvc::interpolate(U.oldTime(), "interpolate(U)");
        U0c -= (mesh().Sf() & U0c)*mesh().Sf()/magSqr(mesh().Sf());
        U0c += phi.oldTime()*mesh().Sf()/magSqr(mesh().Sf());

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                ddtIOobject,
                rDeltaT*ddtPhiCoeff
               *((U0c - U0) & mesh().Sf())
               /fvc::interpolate(1/rA, "interpolate(U)")
            )
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
