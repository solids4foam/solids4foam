/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fv::backwardDdtScheme

Group
    grpFvDdtSchemes

Description
    Second-order backward-differencing ddt using the current and
    two previous time-step values.

SourceFiles
    backwardDdtScheme.H

\*---------------------------------------------------------------------------*/

// There are significant differences between OpenFOAM/FOAM versions so to
// improve readability and maintainability, the ported versions will be given
// in separate files
#ifdef OPENFOAM_ORG
    #include "backwardDdtScheme.foundation.H"
#elif OPENFOAM_COM
    #include "backwardDdtScheme.esi.H"
#endif


// ************************************************************************* //
