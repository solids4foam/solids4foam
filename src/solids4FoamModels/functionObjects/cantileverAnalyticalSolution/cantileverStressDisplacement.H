/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
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

Class
    cantileverStressDisplacement

Description
    Functions that calculate the analytical stress and displacement for a
    bending slender cantilever problem.

    The solution is taken from:
    C.E. Augarde, A.J. Deeks, The use of Timoshenko’s exact solution for a
    cantilever beam in adaptive analysis. Finite Elements in Analysis and
    Design, 44, 2008, 595–601, 10.1016/j.finel.2008.01.010.

    It is assumed that:
    - the beam is fixed on the left
    - the load is applied on the right
    - the beam is L long (left to right) and D deep/thick (top to bottom)
    - the beam is 2-D (1 m in out of plane direction)
    - the origin is at the centre of the fixed (left) patch

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    cantileverStressDisplacement.C

\*---------------------------------------------------------------------------*/

#ifndef cantileverStressDisplacement_H
#define cantileverStressDisplacement_H

#include "symmTensor.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    //- Analytical stress solution
    symmTensor cantileverStress
    (
        const vector& C, // Position vector
        const scalar P, // Load applied in the minus y direction
        const scalar E, // Young's modulus
        const scalar nu, // Poisson's ratio
        const scalar L, // Length of the beam
        const scalar D, // Depth of the beam
        const scalar I // Second moment of area
    );

    //- Analytical displacement solution
    vector cantileverDisplacement
    (
        const vector& C, // position vector
        const scalar P, // Load applied in the minus y direction
        const scalar E, // Young's modulus
        const scalar nu, // Poisson's ratio
        const scalar L, // Length of the beam
        const scalar D, // Depth of the beam
        const scalar I // Second moment of area
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
