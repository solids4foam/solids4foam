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
    RhieChowStab

Description
    Backwards-compatibility alias for `diffStencilLaplacianStab`.

SourceFiles
    RhieChowStab.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef RhieChowStab_H
#define RhieChowStab_H

#include "diffStencilLaplacianStab.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RhieChowStab Declaration
\*---------------------------------------------------------------------------*/

class RhieChowStab
:
    public diffStencilLaplacianStab
{

public:

    //- Runtime type information
    TypeName("RhieChow");

    // Constructors

        //- Construct from dictionary
        RhieChowStab
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
