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

Description
     Block matrix member static data members

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "blockLduMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Note: Matrix debug level defaults to 1.  HJ, 4/Jun/2015
defineNamedTemplateTypeNameAndDebug(blockScalarMatrix, 1);
defineNamedTemplateTypeNameAndDebug(blockVectorMatrix, 1);
defineNamedTemplateTypeNameAndDebug(blockSphericalTensorMatrix, 1);
defineNamedTemplateTypeNameAndDebug(blockSymmTensorMatrix, 1);
defineNamedTemplateTypeNameAndDebug(blockTensorMatrix, 1);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
