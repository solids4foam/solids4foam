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

#include "backwardDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
#include "symmetryFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "wedgeFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// template<>
// tmp<fvMatrix<vector> >
// backwardDdtScheme<vector>::fvmDdt
// (
//     GeometricField<vector, fvPatchField, volMesh>& vf
// )
// {
//     tmp<fvMatrix<vector> > tfvm
//     (
//         new fvMatrix<vector>
//         (
//             vf,
//             vf.dimensions()*dimVol/dimTime
//         )
//     );

//     fvMatrix<vector>& fvm = tfvm();

//     scalar rDeltaT = 1.0/deltaT_();

//     scalar deltaT = deltaT_();
//     scalar deltaT0 = deltaT0_(vf);

//     scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
//     scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
//     scalar coefft0  = coefft + coefft00;

//     fvm.diag() = rDeltaT*mesh().V();
// //     fvm.diag() = (coefft*rDeltaT)*mesh().V();

//     if (mesh().moving())
//     {
//         Info << "Corrected backward ddt" << endl;

//         fvm.source() = rDeltaT*
//         (
//             vf.oldTime().internalField()*mesh().V0()
//           - (
//                 // backward
//                 (
//                     coefft*vf.internalField()*mesh().V()
//                   - coefft0*vf.oldTime().internalField()*mesh().V0()
//                   + coefft00*vf.oldTime().oldTime().internalField()
//                    *mesh().V00()
//                 )
//                 // Euler
//               - (
//                     vf.internalField()*mesh().V()
//                   - vf.oldTime().internalField()*mesh().V0()
//                 )
//             )
//         );
//     }
//     else
//     {
//         fvm.source() = rDeltaT*mesh().V()*
//         (
//             coefft0*vf.oldTime().internalField()
//           - coefft00*vf.oldTime().oldTime().internalField()
//         );
//     }

//     return tfvm;
// }



template<>
tmp<surfaceScalarField>
backwardDdtScheme<vector>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const surfaceScalarField& phi
)
{
    Info << "Consistent backwardDdtPhiCorr" << endl;

    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(U);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

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
//         if (!U.boundaryField()[patchI].coupled())
//         {
//             ddtPhiCoeff.boundaryField()[patchI] = 0.0;
//         }
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
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );

        V00oV.internalField() = mesh().V00()/mesh().V();
        V00oV.correctBoundaryConditions();

        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            surfaceVectorField U0 =
                fvc::interpolate(U.oldTime(), "interpolate(U)");

            surfaceVectorField U00 =
                fvc::interpolate(U.oldTime().oldTime(), "interpolate(U)");


//             surfaceVectorField dU0 =
//                 fvc::interpolate(U.oldTime());
//             forAll(dU0.boundaryField(), patchI)
//             {
//                 if (!U.boundaryField()[patchI].coupled())
//                 {
//                     dU0.boundaryField()[patchI] =
//                         U.oldTime().boundaryField()[patchI]
//                        .patchInternalField();
//                 }
//             }

//             surfaceVectorField dU00 =
//                 fvc::interpolate(U.oldTime().oldTime());
//             forAll(dU00.boundaryField(), patchI)
//             {
//                 if (!U.boundaryField()[patchI].coupled())
//                 {
//                     dU00.boundaryField()[patchI] =
//                         U.oldTime().oldTime().boundaryField()[patchI]
//                        .patchInternalField();
//                 }
//             }

            const surfaceVectorField& Sf =
                mesh().objectRegistry::lookupObject<surfaceVectorField>
                (
                    "Sf"
                );


            U0 += Sf.oldTime()
               *(phi.oldTime() - (Sf.oldTime()&U0))
               /(
                    magSqr(Sf.oldTime())
                  + dimensionedScalar("SMALL", sqr(dimArea), SMALL)
                );

            U00 += Sf.oldTime().oldTime()
               *(phi.oldTime().oldTime() - (Sf.oldTime().oldTime()&U00))
               /(
                    magSqr(Sf.oldTime().oldTime())
                  + dimensionedScalar("SMALL", sqr(dimArea), SMALL)
                );

//             dU0 = Sf.oldTime()
//                *(phi.oldTime() - (Sf.oldTime()&dU0))
//                /(
//                     magSqr(Sf.oldTime())
//                   + dimensionedScalar("SMALL", sqr(dimArea), SMALL)
//                 );

//             dU00 = Sf.oldTime().oldTime()
//                *(phi.oldTime().oldTime() - (Sf.oldTime().oldTime()&dU00))
//                /(
//                     magSqr(Sf.oldTime().oldTime())
//                   + dimensionedScalar("SMALL", sqr(dimArea), SMALL)
//                 );

            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    ddtIOobject,
                    rDeltaT*ddtPhiCoeff
                   *(
                        (
                            (
                                coefft0*fvc::interpolate(V0oV)*U0
                              - coefft00*fvc::interpolate(V00oV)*U00
                            ) & mesh().Sf()
                        )
                      - (
                            fvc::interpolate
                            (
                                coefft0*V0oV*U.oldTime()
                              - coefft00*V00oV*U.oldTime().oldTime(),
                                "interpolate(U)"
                            ) & mesh().Sf()
                        )
                    )/fvc::interpolate(1/rA, "interpolate(U)")
                )
            );

//             return tmp<surfaceScalarField>
//             (
//                 new surfaceScalarField
//                 (
//                     ddtIOobject,
//                     rDeltaT*ddtPhiCoeff
//                    *(
//                         coefft0*fvc::interpolate(V0oV)
//                        *(mesh().Sf()&dU0)
//                       - coefft00
//                        *fvc::interpolate(V00oV)
//                        *(mesh().Sf()&dU00)
//                     )
//                    /fvc::interpolate(1.0/rA)
//                 )
//             );
        }
        else
        {
            FatalErrorIn
            (
                "backwardDdtScheme<vector>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return surfaceScalarField::null();
        }
    }
    else
    {
        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            surfaceVectorField dU0 =
                fvc::interpolate(U.oldTime(), "interpolate(U)");
            forAll(dU0.boundaryField(), patchI)
            {
                if (!U.boundaryField()[patchI].coupled())
                {
                    dU0.boundaryField()[patchI] =
                        U.oldTime().boundaryField()[patchI]
                       .patchInternalField();
                }
            }

            surfaceVectorField dU00 =
                fvc::interpolate(U.oldTime().oldTime(), "interpolate(U)");
            forAll(dU00.boundaryField(), patchI)
            {
                if (!U.boundaryField()[patchI].coupled())
                {
                    dU00.boundaryField()[patchI] =
                        U.oldTime().oldTime().boundaryField()[patchI]
                       .patchInternalField();
                }
            }

            dU0 = mesh().Sf()
               *(phi.oldTime() - (mesh().Sf()&dU0))
               /(
                    magSqr(mesh().Sf())
                  + dimensionedScalar("SMALL", sqr(dimArea), SMALL)
                );

            dU00 = mesh().Sf()
               *(phi.oldTime().oldTime() - (mesh().Sf()&dU00))
               /(
                    magSqr(mesh().Sf())
                  + dimensionedScalar("SMALL", sqr(dimArea), SMALL)
                );

            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    ddtIOobject,
                    rDeltaT*ddtPhiCoeff
                   *(
                        coefft0*(mesh().Sf()&dU0)
                      - coefft00*(mesh().Sf()&dU00)
                    )
                   /fvc::interpolate(1.0/rA, "interpolate(U)")
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "backwardDdtScheme<vector>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return surfaceScalarField::null();
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
