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
    readSolidTimeControls

Description
    Read the time controls for adaptive time-stepping.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

// Adjustable time-step
const Switch adjustTimeStep
(
    runTime.controlDict().lookupOrDefault<Switch>("adjustTimeStep", false)
);

// Maximum allowed time-step
const scalar maxDeltaT
(
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT)
);

// Minimum allowed time-step
const scalar minDeltaT
(
    runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", 0.0)
);

// Increase the time-step at most every increaseTimeStepFreq time-steps
int increaseTimeStep = 0;
const int increaseTimeStepFreq = 4;
(
    runTime.controlDict().lookupOrDefault<int>("increaseTimeStepFreq", 4)
);

if (adjustTimeStep)
{
    Info<< nl << "Adjustable time-step" << nl
        << "    Max deltaT " << maxDeltaT << nl
        << "    Min deltaT " << minDeltaT << nl << endl;
}

// ************************************************************************* //
