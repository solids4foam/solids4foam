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

#include "AitkenCouplingInterface.H"
// #include "volFields.H"
// #include "fvm.H"
// #include "fvc.H"
// #include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
//#include "adjustPhi.H"
//#include "findRefCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(AitkenCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface, AitkenCouplingInterface, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AitkenCouplingInterface::AitkenCouplingInterface
(
    dynamicFvMesh& fluidMesh,
    fvMesh& solidMesh
)
:
    fluidSolidInterface(typeName, fluidMesh, solidMesh),
    relaxationFactor_(readScalar(lookup("relaxationFactor"))),
    aitkenRelaxationFactor_(relaxationFactor_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AitkenCouplingInterface::updateDisplacement()
{
    Info<< nl << "Time = " << fluid().runTime().timeName()
        << ", iteration: " << outerCorr() << endl;

    if (outerCorr() < 3)
    {
        Info<< "Current fsi under-relaxation factor: "
            << relaxationFactor_ << endl;

        fluidZonePointsDisplPrev() = fluidZonePointsDispl();

        fluidZonePointsDispl() += relaxationFactor_*residual();
    }
    else
    {
        aitkenRelaxationFactor_ =
            -aitkenRelaxationFactor_
           *(
                sum
                (
                    residualPrev()
                  & (residual() - residualPrev())
                )
               /(
                    sum
                    (
                        (residual() - residualPrev())
                      & (residual() - residualPrev())
                    )
                )
            );

        if (Pstream::parRun())
        {
            if (!Pstream::master())
            {
                aitkenRelaxationFactor_ = 0.0;
            }

            // Pass to all procs
            reduce(aitkenRelaxationFactor_, sumOp<scalar>());
        }

        aitkenRelaxationFactor_ = mag(aitkenRelaxationFactor_);

        if (aitkenRelaxationFactor_ > 1)
        {
            aitkenRelaxationFactor_ = relaxationFactor_;
        }

        Info<< "Current fsi under-relaxation factor (Aitken): "
            << aitkenRelaxationFactor_ << endl;

        fluidZonePointsDisplPrev() = fluidZonePointsDispl();

        fluidZonePointsDispl() += aitkenRelaxationFactor_*residual();
    }


    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    if (Pstream::parRun())
    {
        if (!Pstream::master())
        {
            fluidZonePointsDispl() = vector::zero;
        }

        //- pass to all procs
        reduce(fluidZonePointsDispl(), sumOp<vectorField>());

        label globalFluidZoneIndex =
            findIndex(fluid().globalFaceZones(), fluidZoneIndex());

        if (globalFluidZoneIndex == -1)
        {
            FatalErrorIn
            (
                "fluidSolidInterface::updateDisplacement()"
            )   << "global zone point map is not availabel"
                << abort(FatalError);
        }

        const labelList& map =
            fluid().globalToLocalFaceZonePointMap()[globalFluidZoneIndex];

        if (!Pstream::master())
        {
            vectorField fluidZonePointsDisplGlobal =
                fluidZonePointsDispl();

            forAll(fluidZonePointsDisplGlobal, globalPointI)
            {
                label localPoint = map[globalPointI];

                fluidZonePointsDispl()[localPoint] =
                    fluidZonePointsDisplGlobal[globalPointI];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
