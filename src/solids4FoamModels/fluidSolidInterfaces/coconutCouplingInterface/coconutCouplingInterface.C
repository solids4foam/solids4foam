/*---------------------------------------------------------------------------*\
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

#include "coconutCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coconutCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface, coconutCouplingInterface, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar coconutCouplingInterface::computeResidual()
{
    // Placeholder for residual calculation logic
    // This should include the logic from iqni.py's compute step
    return 0.0; // Example value, needs to be replaced with actual computation
}


tmp<vectorField> coconutCouplingInterface::predictorPredict
(
    const vectorField& x
) const
{
    tmp<vectorField> tresult(new vectorField(x.size(), vector::zero));

    // CAREFUL: we store xPrev in the opposite order to what is used in coconut
    // xPrev[0]: old-old-old
    // xPrev[1]: old-old
    // xPrev[2]: old

    if (xPrev_.size() == 3)
    {
        // Quadratic predictor
        tresult = 3.0*xPrev_[2] - 3.0*xPrev_[1] + 1.0*xPrev_[0];
    }
    else if (xPrev_.size() == 2)
    {
        // Linear predictor
        tresult = 2.0*xPrev_[1] - xPrev_[0];
    }
    else if (xPrev_.size() == 1)
    {
        // Constant predictor
        tresult = xPrev_[0];
    }

    return tresult;
}


tmp<vectorField> coconutCouplingInterface::solver0SolveSolutionStep
(
    const vectorField& x
)
{
    // Fluid solver
    // Return interface traction given the interface motion

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(x.size(), vector::zero));

    // For now, only implemented for one interface
    if (fluid().globalPatches().size() > 1)
    {
        FatalErrorInFunction
            << "Only implemented for one interface" << abort(FatalError);
    }
    const int interfaceI = 0;

    // Interpolate x from solid faces to fluid points

    // Take references to zones
    const standAlonePatch& fluidZone =
        fluid().globalPatches()[interfaceI].globalPatch();
    const standAlonePatch& solidZone =
        solid().globalPatches()[interfaceI].globalPatch();

    // Interpolate x from solid faces to solid points
    PrimitivePatchInterpolation<standAlonePatch> interp
    (
        solid().globalPatches()[interfaceI].globalPatch()
    );
    const vectorField xPoints(interp.faceToPointInterpolate(x));

    // Initialise x at fluid points
    vectorField xFluidPoints(fluidZone.nPoints(), vector::zero);

    // Transfer displacement field from the solid to the fluid
    interfaceToInterfaceList()[interfaceI].transferPointsZoneToZone
    (
        solidZone,              // from zone
        fluidZone,              // to zone
        xPoints,                // from field
        xFluidPoints            // to field
    );


    // Apply the interface motion
    fluidZonesPointsDisplsPrev()[interfaceI] =
        fluidZonesPointsDispls()[interfaceI];
    fluidZonesPointsDispls()[interfaceI] = xFluidPoints;


    // Solve the fluid
    fluid().evolve();


    // Retrieve the interface traction

    // Calculate total traction of fluid zone
    vectorField fluidZoneTotalTraction
    (
        fluid().faceZoneViscousForce(0)
      - fluid().faceZonePressureForce(0)*fluidZone.faceNormals()
    );

    // Initialise the solid zone traction field that is to be interpolated
    // from the fluid zone
    vectorField solidZoneTotalTraction(solidZone.size(), vector::zero);

    // Transfer the field frm the fluid interface to the solid interface
    interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
    (
        fluidZone,                 // from zone
        solidZone,                 // to zone
        fluidZoneTotalTraction,    // from field
        solidZoneTotalTraction     // to field
    );

    // Flip traction sign after transferring from fluid to solid
    solidZoneTotalTraction = -solidZoneTotalTraction;

    tresult.ref() = solidZoneTotalTraction;

    return tresult;
}


tmp<vectorField> coconutCouplingInterface::solver1SolveSolutionStep
(
    const vectorField& x
)
{
    // Solid solver
    // Return interface motion given the interface traction

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(x.size(), vector::zero));

    // For now, only implemented for one interface
    if (fluid().globalPatches().size() > 1)
    {
        FatalErrorInFunction
            << "Only implemented for one interface" << abort(FatalError);
    }
    const int interfaceI = 0;

    if (coupled())
    {
        // Apply the traction to the solid
        solid().setTraction
        (
            interfaceI,
            solidPatchIndices()[interfaceI],
            x // solidZoneTotalTraction
        );


        // Solve the solid
        solid().evolve();


        // Retrieve the solid interface motion

        // Take references to zones
        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        const standAlonePatch& solidZone =
            solid().globalPatches()[interfaceI].globalPatch();

        // Calculate the point displacements of the solid interface
        const vectorField solidZonePointsDisplsAtSolid
        (
            solid().faceZonePointDisplacementIncrement(interfaceI)
        );
        // const vectorField solidZonePointsTotDisplsAtSolid
        // (
        //     solid().faceZonePointDisplacementOld(interfaceI)
        // );

        // Initialise point displacement field at fluid interface
        vectorField solidZonePointsTotDispl
        (
            solidZonesPointsDispls()[interfaceI].size(), vector::zero
        );

        // Transfer displacement field from the solid to the fluid
        interfaceToInterfaceList()[interfaceI].transferPointsZoneToZone
        (
            solidZone,                              // from zone
            fluidZone,                              // to zone
            solidZonePointsDisplsAtSolid,           // from field
            solidZonesPointsDispls()[interfaceI]    // to field
        );
        // interfaceToInterfaceList()[interfaceI].transferPointsZoneToZone
        // (
        //     solidZone,                              // from zone
        //     fluidZone,                              // to zone
        //     solidZonePointsTotDisplsAtSolid,        // from field
        //     solidZonePointsTotDispl                 // to field
        // );

        // Update interface residuals
        residualsPrev()[interfaceI] = residuals()[interfaceI];

        // Residual calculated based on the displacement increments of the
        // fluid and solid interfaces
        residuals()[interfaceI] =
            solidZonesPointsDispls()[interfaceI]
          - fluidZonesPointsDispls()[interfaceI];

        // Interface displacement
        // CHECK: is this correct?
        // tresult.ref() = residuals()[interfaceI];
        tresult.ref() = solidZonesPointsDispls()[interfaceI];
    }


    return tresult;
}


void coconutCouplingInterface::finalizeIteration(const vectorField& r)
{
    // Update the convergence criterion with the current residual
    convergenceCriteria_.update(r);

    // // Print iteration information
    // printIterationInfo(r);

    // // Output iteration data
    // outputIteration(r);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coconutCouplingInterface::coconutCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    residual_(0.0),
    relaxationFactor_(readScalar(fsiProperties().lookup("relaxationFactor"))),
    model_(readInt(fsiProperties().lookup("q"))),
    convergenceCriteria_(outerCorrTolerance(), nOuterCorr()),
    x_(solid().globalPatches()[0].globalPatch().size(), vector::zero),
    xPrev_(0)
{
    if (solid().globalPatches().size() > 1)
    {
        FatalErrorInFunction
            << "Only implemented for one interface" << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coconutCouplingInterface::initialize()
{
    Info<< "Initializing coconutCouplingInterface..." << endl;
    residual_ = 0.0;
}


void coconutCouplingInterface::solveTimeStep()
{
    Info << "Solving time step..." << endl;

    // Initial value using predictor
    vectorField x(predictorPredict(x_));

    // First coupling iteration
    vectorField y(solver0SolveSolutionStep(x));
    vectorField xt(solver1SolveSolutionStep(y));

    // Compute residual: r = xt - x
    vectorField r(xt - x);

    // Add the residual and solution to the model
    model_.add(r, xt);

    // Finalize the iteration
    finalizeIteration(r);

    // Coupling iteration loop
    while (!convergenceCriteria_.isSatisfied())
    {
        vectorField dx(r.size(), vector::zero);

        if (!model_.isReady())
        {
            dx = relaxationFactor_*r;
        }
        else
        {
            const vectorField dr(-r);
            dx = model_.predict(dr) - dr;
        }

        // Update x
        x += dx;

        // Solve with updated x
        y = solver0SolveSolutionStep(x);
        xt = solver1SolveSolutionStep(y);

        // Compute residual: r = xt - x
        r = xt - x;

        // Add the residual and solution to the model
        model_.add(r, xt);

        // Finalize the iteration
        finalizeIteration(r);
    }
}


//void coconutCouplingInterface::restart(int timeStep)
//{}


void coconutCouplingInterface::finalize()
{
    Info<< "Finalizing coconutCouplingInterface..." << endl;
    // Finalization logic can be added here if needed
}


bool coconutCouplingInterface::evolve()
{
    // Check if coupling switch needs to be updated
    if (!coupled())
    {
        updateCoupled();
    }

    // Initialise solution step
    model_.initializeSolutionStep();
    convergenceCriteria_.initializeSolutionStep();

    // Call the coconut solution procedure for a time step
    solveTimeStep();

    // Finalise solution step
    model_.finaliseSolutionStep();
    convergenceCriteria_.finaliseSolutionStep();

    // Update xPrev, which is used by predictorPredict
    if (xPrev_.size() < 3)
    {
        // Most recent value at the end (opposite order used in coconut)
        xPrev_.append(new vectorField(x_));
    }
    else
    {
        xPrev_[0] = xPrev_[1];  // Old-old
        xPrev_[1] = xPrev_[2];  // Old
        xPrev_[2] = x_;         // Current
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
