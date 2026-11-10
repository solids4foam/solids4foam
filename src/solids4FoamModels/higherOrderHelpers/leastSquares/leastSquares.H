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
    Foam::leastSquares

Description
    Point-value least-squares reconstruction with one polynomial per cell.
    The stored cell value is the field value at the cell centre, not a cell
    average. For each cell a polynomial of order N, centred at the cell
    centre and with its constant fixed to the cell value, is fitted by
    weighted least squares to the point values of the stencil cells. The cell
    polynomial is then evaluated at the cell centre, the face centres, the
    face quadrature points and arbitrary points for values and derivatives.

    This differs from movingLeastSquares, which fits a new polynomial at each
    evaluation point, and from kExactLeastSquares, whose stored values are
    cell averages and whose basis is corrected with the cell moments. The
    scheme is the kExactLeastSquares formulation with every positive-order
    cell moment set to zero, so no cell moments are computed or exchanged.

    Selected with "type leastSquares;" in the reconstruction dictionary. This
    is unrelated to the OpenFOAM gradient scheme of the same name, which lives
    in a different run-time selection table.

    The per-cell polynomial machinery is implemented in cellLeastSquares.

SourceFiles
    leastSquares.C

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef leastSquares_H
#define leastSquares_H

#include "cellLeastSquares.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class leastSquares Declaration
\*---------------------------------------------------------------------------*/

class leastSquares
:
    public cellLeastSquares
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        leastSquares(const leastSquares&);

        //- Disallow default bitwise assignment
        void operator=(const leastSquares&);


public:

    //- Runtime type information
    TypeName("leastSquares");


    // Constructors

        //- Construct from components
        leastSquares
        (
            const fvMesh& mesh,
            const boolList& includePatchInStencils,
            const dictionary& dict
        );


    //- Destructor
    virtual ~leastSquares();


    // Member Functions

        //- The stored cell values are point values at the cell centres
        virtual bool cellAverageUnknowns() const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
