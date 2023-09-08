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
    nonLinearGeometry

Description
    Enumerator inidacting whether the mathematical is employs a linear geometry,
    total Lagrangian nonlinear geometry or updated Lagrangian nonlinear geometry
    approach.

Author
    Hrvoje Jasak
    Philip Cardiff, UCD.

SourceFiles
    nonLinearGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinearGeometry_H
#define nonLinearGeometry_H

#ifdef OPENFOAM_COM
#include "Enum.H"
#else
#include "NamedEnum.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class nonLinearGeometry Declaration
\*---------------------------------------------------------------------------*/

class nonLinearGeometry
{
public:

        //- Non-linear solver options
        enum nonLinearType
        {
            LINEAR_GEOMETRY,
            UPDATED_LAGRANGIAN,
            TOTAL_LAGRANGIAN
        };

#ifdef OPENFOAM_COM
        static const Enum<nonLinearType> nonLinearNames_;
#else
        static const NamedEnum<nonLinearType, 3> nonLinearNames_;
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
