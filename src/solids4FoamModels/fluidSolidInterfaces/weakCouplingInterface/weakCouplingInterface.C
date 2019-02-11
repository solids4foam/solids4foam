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

#include "weakCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(weakCouplingInterface, 0);
addToRunTimeSelectionTable
(
    physicsModel, weakCouplingInterface, fluidSolidInteraction
);
addToRunTimeSelectionTable
(
    fluidSolidInterface, weakCouplingInterface, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

weakCouplingInterface::weakCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    solidZoneTraction_
    (
        IOobject
        (
            "solidZoneTraction",
            runTime.timeName(),
            fluidMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        vectorField()
    ),
    solidZoneTractionPrev_(),
    predictedSolidZoneTraction_
    (
        IOobject
        (
            "predictedSolidZoneTraction",
            runTime.timeName(),
            fluidMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        vectorField()
    ),
    relaxationFactor_
    (
        fsiProperties().lookupOrDefault<scalar>("relaxationFactor", 0.01)
    )
{
    // Initialize solid zone traction fields
    if (solidZoneTraction_.size() == 0)
    {
        solidZoneTraction_ =
            vectorField
            (
                solid().globalPatch().globalPatch().size(),
                vector::zero
            );
    }

    solidZoneTractionPrev_ =
        vectorField
        (
            solid().globalPatch().globalPatch().size(),
            vector::zero
        );

    if (predictedSolidZoneTraction_.size() == 0)
    {
        predictedSolidZoneTraction_ =
            vectorField
            (
                solid().globalPatch().globalPatch().size(),
                vector::zero
            );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool weakCouplingInterface::evolve()
{
    initializeFields();

    updateInterpolatorAndGlobalPatches();

    solid().evolve();

    updateWeakDisplacement();

    moveFluidMesh();

    fluid().evolve();

    updateWeakTraction();

    solid().updateTotalFields();

    return 0;
}


void weakCouplingInterface::initializeFields()
{
    predictedSolidZoneTraction_ =
        vectorField
        (
            solid().globalPatch().globalPatch().size(),
            vector::zero
        );

    fluidSolidInterface::initializeFields();
}


bool weakCouplingInterface::updateWeakDisplacement()
{
    const vectorField solidZonePointsDisplAtSolid =
        solid().faceZonePointDisplacementIncrement();

    solidZonePointsDispl() =
        ggiInterpolator().slaveToMasterPointInterpolate
        (
            solidZonePointsDisplAtSolid
        );

    const vectorField solidZonePointsTotDisplAtSolid =
        solid().faceZonePointDisplacementOld();

    const vectorField solidZonePointsTotDispl =
        ggiInterpolator().slaveToMasterPointInterpolate
        (
            solidZonePointsTotDisplAtSolid
        );


    residualPrev() = residual();

    residual() = solidZonePointsDispl() - fluidZonePointsDispl();

    fluidZonePointsDisplPrev() = fluidZonePointsDispl();

    fluidZonePointsDispl() += residual();

    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    fluidSolidInterface::syncFluidZonePointsDispl(fluidZonePointsDispl());


    // Update elasticWallPressure boundary conditions, if found
    fluidSolidInterface::updateElasticWallPressureAcceleration();


    // Calculate residual norm

    scalar residualNorm = ::sqrt(gSum(magSqr(residual())));
    scalar residualNorm_2 = residualNorm;

    if (residualNorm > maxResidualNorm())
    {
        maxResidualNorm() = residualNorm;
    }

    residualNorm /= maxResidualNorm() + SMALL;

    Info<< "Current fsi relative residual norm: " << residualNorm << endl;

    interfacePointsDisplPrev() = interfacePointsDispl();

    interfacePointsDispl() = solidZonePointsDispl();

    const vectorField intTotDispl =
        interfacePointsDispl() + solidZonePointsTotDispl;
    const scalar intTotDisplNorm = ::sqrt(gSum(magSqr(intTotDispl)));
    if (intTotDisplNorm > maxIntDisplNorm())
    {
        maxIntDisplNorm() = intTotDisplNorm;
    }

    residualNorm_2 /= maxIntDisplNorm() + SMALL;

    Info<< "Alternative fsi residual: " << residualNorm_2 << endl;

    return min(residualNorm_2, residualNorm);
}


void weakCouplingInterface::updateWeakTraction()
{
    Info<< "Update weak traction on solid patch" << endl;

    // Calc fluid traction

    const vectorField& p =
        fluid().globalPatch().globalPatch().localPoints();
    const faceList& f =
        fluid().globalPatch().globalPatch().localFaces();

    vectorField n(f.size(), vector::zero);
    forAll(n, faceI)
    {
        n[faceI] = f[faceI].normal(p);
        n[faceI] /= mag(n[faceI]);
    }

    const vectorField fluidZoneTraction =
        fluid().faceZoneViscousForce() - fluid().faceZonePressureForce()*n;

    const vectorField fluidZoneTractionAtSolid =
        ggiInterpolator().masterToSlave
        (
            -fluidZoneTraction
        );

    const scalar beta_ = relaxationFactor_;

    solidZoneTractionPrev_ = solidZoneTraction_;

    solidZoneTraction_ =
        beta_*fluidZoneTractionAtSolid
      + (1.0 - beta_)*predictedSolidZoneTraction_;

    FatalError
        << "Check beta here!" << abort(FatalError);

    predictedSolidZoneTraction_ =
        2.0*solidZoneTraction_ - solidZoneTractionPrev_;

    if (coupled())
    {
        Info<< "Setting weak traction on solid patch" << endl;

        solid().setTraction
        (
            solidPatchIndex(),
            predictedSolidZoneTraction_
        );
    }

    // Total force at the fluid side of the interface
    {
        const vectorField& p =
            fluid().globalPatch().globalPatch().localPoints();
        const faceList& f =
            fluid().globalPatch().globalPatch().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        const vector totalTractionForce = sum(fluidZoneTraction*mag(S));

        Info<< "Total force (fluid) = " << totalTractionForce << endl;
    }

    // Total force at the solid side of the interface
    {
        const vectorField& p =
            solid().globalPatch().globalPatch().localPoints();
        const faceList& f =
            solid().globalPatch().globalPatch().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        const vector totalTractionForce =
            sum(fluidZoneTractionAtSolid*mag(S));

        Info<< "Total force (solid) = " << totalTractionForce << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
