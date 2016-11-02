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

File
    setSolidDeltaT

Description
    Update the time-step based on the solid model.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep && runTime.timeIndex() > 0)
{
    Info<< "Adjusting the time step" << endl;

    // Calculate the time-step scaling factor and limit it so as it does not
    // grow or shrink too quickly
    const scalar deltaTFact =
        max
        (
            min
            (
                solid().newDeltaT()/runTime.deltaTValue(),
                1.2
            ),
            0.8
        );

    // Increase the time-step at most every increaseTimeStepFreq time-steps
    if (deltaTFact > 1.0)
    {
        increaseTimeStep++;
    }
    else
    {
        increaseTimeStep = 0;
    }

    if (deltaTFact <= 1.0 || increaseTimeStep == increaseTimeStepFreq)
    {
        increaseTimeStep = 0;

        // Set the time-step and limit it to the specified max and min values

        scalar newDeltaT = deltaTFact*runTime.deltaTValue();

        if (newDeltaT > maxDeltaT)
        {
            Info<< nl << "    Setting deltaT to the maxDeltaT" << endl;
            newDeltaT = maxDeltaT;
        }
        else if (newDeltaT < minDeltaT)
        {
            Info<< nl << "    Setting deltaT to the minDeltaT" << endl;
            newDeltaT = minDeltaT;
        }
        else if (deltaTFact > 1.0)
        {
            Info<< nl << "    Increasing the time-step" << endl;
        }
        else
        {
            Info<< nl << "    Decreasing the time-step" << endl;
        }

        runTime.setDeltaT(newDeltaT);
    }
    else
    {
        Info<< nl << "    The time-step is unchanged" << endl;
    }

    Info<< nl << "    Time step = " <<  runTime.deltaT().value() << nl
        << "    Minimum time step = "<< minDeltaT << nl
        << "    Maximum time step = "<< maxDeltaT << nl
        << "    Desired time step scaling factor = "<< deltaTFact
        << nl << nl << endl;
}

// ************************************************************************* //
