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

\*---------------------------------------------------------------------------*/

#include "kExactLeastSquares.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kExactLeastSquares, 0);
addToRunTimeSelectionTable
(
    leastSquaresScheme,
    kExactLeastSquares,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kExactLeastSquares::kExactLeastSquares
(
    const fvMesh& mesh,
    const boolList& includePatchInStencils,
    const dictionary& dict
)
:
    cellLeastSquares(mesh, includePatchInStencils, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

kExactLeastSquares::~kExactLeastSquares()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

const List<symmTensor>& kExactLeastSquares::secondOrderCellMoments() const
{
    return quadrature().secondOrderCellMoments();
}


const List<symmTensor3rdOrder>&
kExactLeastSquares::thirdOrderCellMoments() const
{
    return quadrature().thirdOrderCellMoments();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
