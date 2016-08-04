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

Application
    fsiFoam

Description
    Finite volume fluid solid interaction solver based on partitioned approach
    and strong coupling, where the fluid, solid and interface models are run-
    time selectable.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createSolidMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Create fluid-solid interface coupling class
    autoPtr<fluidSolidInterface> fsi =
        fluidSolidInterface::New(mesh, solidMesh);

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        fsi().initializeFields();

        fsi().updateInterpolator();

        scalar residualNorm = 0;

        do
        {
            fsi().outerCorr()++;

            // Transfer the displacement from the solid to the fluid
            fsi().updateDisplacement();

            // Move the fluid mesh
            fsi().moveFluidMesh();

            // Solve fluid
            fsi().fluid().evolve();

            // Transfer the force from the fluid to the solid
            fsi().updateForce();

            // Solve solid
            fsi().solid().evolve();

            // Calculate the FSI residual
            residualNorm = fsi().updateResidual();
        }
        while
        (
            (residualNorm > fsi().outerCorrTolerance())
         && (fsi().outerCorr() < fsi().nOuterCorr())
        );

        fsi().solid().updateTotalFields();

        if (runTime.outputTime())
        {
            fsi().solid().writeFields(runTime);
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
