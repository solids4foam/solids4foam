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

Class
    plateHoleAnalyticalFields

Description
    Shared Kirsch analytical stress, displacement, pressure and traction
    helpers for the plate-hole tutorials.

SourceFiles
    plateHoleAnalyticalFields.C

\*---------------------------------------------------------------------------*/

#ifndef plateHoleAnalyticalFields_H
#define plateHoleAnalyticalFields_H

#include "scalar.H"
#include "symmTensor.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class plateHoleAnalyticalFields Declaration
\*---------------------------------------------------------------------------*/

class plateHoleAnalyticalFields
{
public:

    //- Analytical stress field
    static symmTensor stress
    (
        const vector& C,
        const scalar farFieldTractionX,
        const scalar holeRadius
    );

    //- Analytical displacement field
    static vector displacement
    (
        const vector& C,
        const scalar farFieldTractionX,
        const scalar holeRadius,
        const scalar mu,
        const scalar kappa
    );

    //- Analytical displacement field from E and nu
    static vector displacement
    (
        const vector& C,
        const scalar farFieldTractionX,
        const scalar holeRadius,
        const scalar E,
        const scalar nu,
        const bool planeStress
    );

    //- Hydrostatic pressure from the analytical stress state
    static scalar hydPressure
    (
        const vector& C,
        const scalar farFieldTractionX,
        const scalar holeRadius,
        const scalar nu
    );

    //- Analytical traction
    static vector traction
    (
        const vector& C,
        const vector& n,
        const scalar farFieldTractionX,
        const scalar holeRadius
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
