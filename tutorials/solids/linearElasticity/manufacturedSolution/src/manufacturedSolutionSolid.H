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
    Foam::solidModels::manufacturedSolutionSolid

Description
    Linear solid with the manufactured body force for OpenFOAM.org and
    foam-extend. Source coefficients use the same constant/fvOptions dictionary
    as the OpenFOAM.com tutorial.

SourceFiles
    manufacturedSolutionSolid.C

\*----------------------------------------------------------------------------*/

#ifndef manufacturedSolutionSolid_H
#define manufacturedSolutionSolid_H

#ifndef OPENFOAM_COM

#include "linGeomTotalDispSolid.H"
#include "manufacturedSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                     Class manufacturedSolutionSolid Declaration
\*---------------------------------------------------------------------------*/

class manufacturedSolutionSolid
:
    public linGeomTotalDispSolid
{
    // Private Data

        //- Source coefficients, shared with the OpenFOAM.com fvOption
        const IOdictionary sourceDict_;

        //- Manufactured solution supplying the body force
        manufacturedSolution mms_;


protected:

    //- Return the manufactured source per unit volume
    virtual tmp<volVectorField> fvOptionsSource() const;


public:

    //- Runtime type information
    TypeName("manufacturedSolutionSolid");

    // Constructors

        //- Construct from time and region
        manufacturedSolutionSolid(Time& runTime, const word& region);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
#endif

// ************************************************************************* //
