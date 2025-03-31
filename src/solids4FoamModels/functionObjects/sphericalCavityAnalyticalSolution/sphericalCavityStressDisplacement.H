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

Functions
    sphericalCavityStressDisplacement

Description
    Functions that generates the analytical displacement and stress fields
    for the spherical cavity (radius = a) in an infinite linear elastic solid
    (Young's modulus of E and Poisson's ration of nu) subjected to a distant
    uniaxial stress T0 applied in the positive Z direction.

    The analytical solution requires the following inputs:
        - farFieldTractionZ : far-field traction in the positive Z direction
        - cavityRadius : radius of the cavity centred on the origin
        - E : Young's modulus
        - nu : Poisson's ratio

    The analytical solution for stress is taken from:

    R.V. Southwell M.A. F.R.S. & H.J. Gough M.B.E. B.Sc. (1926) VI. On the
    concentration of stress in the neighbourhood of a small spherical flaw;
    and on the propagation of fatigue fractures in “Statistically Isotropic”
    materials , The London, Edinburgh, and Dublin Philosophical Magazine and
    Journal of Science: Series 7, 1:1, 71-97, DOI: 10.1080/14786442608633614

    The analytical solution for displacement is taken from:

    Allan F. Bower (2009) Applied Mechanics of Solids, CRC Press, Taylor and
    Francis Group.

    Note: Bower also gives the stress solutions but they seem to have errors
    as the axial (zz) stress formula does not provide the correct stress
    concentration factor at the edge of the cavity.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    sphericalCavityStressDisplacement.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalCavityStressDisplacement_H
#define sphericalCavityStressDisplacement_H

#include "vector.H"
#include "symmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    //- Kronecker delta function
    scalar KroneckerDelta
    (
       const int i, const int j
    );

    //- Calculate the displacement component i at a coordinate x
    //  From Bower (2009)
    scalar sphericalCavityDisplacementComponent
    (
       const int i,
       const scalar nu,
       const scalar sigma_0,
       const scalar E,
       const scalar a,
       const vector& x
    );

    //- Calculate the displacement vector at a coordinate x
    //  From Bower (2009)
    vector sphericalCavityDisplacement
    (
       const scalar nu,
       const scalar sigma_0,
       const scalar E,
       const scalar a,
       const vector& x
    );

    //- Calculate the stress at point x
    //  From Southwell (1926)
    symmTensor sphericalCavityStress
    (
       const scalar T,     // distant ZZ stress
       const scalar nu,    // Poisson's ratio
       const scalar a,     // Sphere radius
       const vector& x     // Position
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
