/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "robinInterfaceHelper.H"
#include "fvc.H"
#include "EulerDdtScheme.H"
#include "backwardDdtScheme.H"
#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "elasticWallVelocityFvPatchVectorField.H"
#include "elasticWallPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(robinInterfaceHelper, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::robinInterfaceHelper::robinInterfaceHelper
(
    const volVectorField& U,
    const volScalarField& p
)
:
    patchIDs_(0),
    ddtScheme_()
{
    // Check for Robin boundary conditions

    labelHashSet robinPatches;

    forAll(p.boundaryField(), patchI)
    {
        if
        (
            (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                )
             && isA<elasticSlipWallVelocityFvPatchVectorField>
                (
                    U.boundaryField()[patchI]
                )
            )
         || (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                )
             && isA<elasticWallVelocityFvPatchVectorField>
                (
                    U.boundaryField()[patchI]
                )
            )
        )
        {
            robinPatches.insert(patchI);
        }
    }

    patchIDs_ = robinPatches.toc();


    // Lookup U ddt scheme
    ddtScheme_ = word
    (
        U.mesh().schemesDict().ddtScheme("ddt(" + U.name() +')')
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::robinInterfaceHelper::setInterfaceFluxToZero
(
    surfaceScalarField& phi
) const
{
    forAll(patchIDs_, pI)
    {
        const label patchID = patchIDs_[pI];
        phi.boundaryField()[patchID] = 0;

        Info<< type() << ": setting phi to 0 on patch " << patchID << endl;
    }
}


void Foam::robinInterfaceHelper::updateInterface
(
    surfaceScalarField& phi,
    volScalarField& rAU
) const
{
    forAll(patchIDs_, pI)
    {
        const label patchID = patchIDs_[pI];
        const fvMesh& mesh = phi.mesh();

        if (ddtScheme_ == fv::EulerDdtScheme<vector>::typeName)
        {
            phi.boundaryField()[patchID] = phi.oldTime().boundaryField()[patchID];

            rAU.boundaryField()[patchID] = mesh.time().deltaT().value();
        }
        else if
        (
            ddtScheme_ == fv::backwardDdtScheme<vector>::typeName
        )
        {
            if (mesh.time().timeIndex() == 1)
            {
                phi.boundaryField()[patchID] =
                    phi.oldTime().boundaryField()[patchID];

                rAU.boundaryField()[patchID] = mesh.time().deltaT().value();

                phi.oldTime().oldTime();
            }
            else
            {
                const scalar deltaT = mesh.time().deltaT().value();
                const scalar deltaT0 = mesh.time().deltaT0().value();

                const scalar Cn = 1 + deltaT/(deltaT + deltaT0);
                const scalar Coo =
                    deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
                const scalar Co = Cn + Coo;

                phi.boundaryField()[patchID] =
                    (Co/Cn)*phi.oldTime().boundaryField()
                    [
                        patchID
                    ]
                  - (Coo/Cn)*phi.oldTime().oldTime().boundaryField()
                    [
                        patchID
                    ];

                rAU.boundaryField()[patchID] = deltaT/Cn;
            }
        }
    }
}

void Foam::robinInterfaceHelper::correctInterface
(
    const volScalarField& p,
    surfaceScalarField& phi,
    volScalarField& rAU
) const
{
    forAll(patchIDs_, pI)
    {
        const label patchID = patchIDs_[pI];
        const fvMesh& mesh = phi.mesh();

        if (ddtScheme_ == fv::EulerDdtScheme<vector>::typeName)
        {
            phi.boundaryField()[patchID] =
                phi.oldTime().boundaryField()[patchID]
              - p.boundaryField()
                [
                    patchID
                ].snGrad()*mesh.time().deltaT().value()
               *mesh.magSf().boundaryField()
                [
                    patchID
                ];
        }
        else if
        (
            ddtScheme_ == fv::backwardDdtScheme<vector>::typeName
        )
        {
            if (mesh.time().timeIndex() == 1)
            {
                phi.boundaryField()[patchID] =
                    phi.oldTime().boundaryField()[patchID];

                rAU.boundaryField()[patchID] = mesh.time().deltaT().value();

                phi.oldTime().oldTime();
            }
            else
            {
                const scalar deltaT = mesh.time().deltaT().value();
                const scalar deltaT0 = mesh.time().deltaT0().value();

                const scalar Cn = 1 + deltaT/(deltaT + deltaT0);
                const scalar Coo = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
                const scalar Co = Cn + Coo;

                phi.boundaryField()[patchID] =
                    (Co/Cn)*phi.oldTime().boundaryField()
                    [
                        patchID
                    ]
                  - (Coo/Cn)*phi.oldTime().oldTime().boundaryField()
                    [
                        patchID
                    ];

                rAU.boundaryField()[patchID] = deltaT/Cn;
            }
        }
    }
}


// ************************************************************************* //
