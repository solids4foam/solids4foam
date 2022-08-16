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

Class
    pimpleFluid

Description
    pimpleFoam fluid model.

Author
    Based on pimpleFoam.
    Ported by Philip Cardiff, UCD.

SourceFiles
    pimpleFluid.C

\*---------------------------------------------------------------------------*/

// There are significant differences between OpenFOAM/FOAM versions so to
// improve readability and maintainability, the ported versions will be given
// in separate files
#ifdef OPENFOAMFOUNDATION
    #include "pimpleFluid.foundation.H"
#elif OPENFOAMESI
    #include "pimpleFluid.esi.H"
#else
    #include "pimpleFluid.foamextend.H"
#endif

// ************************************************************************* //
