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
    Function to calculate eigen values and eigen vectors of volSymmTensorField
    Using the main procedure/code from here:
    http://barnesc.blogspot.ie/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
    and code here:
    http://www.connellybarnes.com/code/c/eig3-1.0.0.zip
    Note: built-in OpenFOAM functions mess-up on a number of different tensors.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef eig3Field_H
#define eig3Field_H

#include "eig3.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void eig3Field
(
    const volTensorField& A, volTensorField& V, volVectorField& d
);

void eig3Field
(
    const volSymmTensorField& A, volTensorField& V, volVectorField& d
);

#ifndef OPENFOAM_NOT_EXTEND
void eig3Field
(
    const volSymmTensorField& A, volTensorField& V, volDiagTensorField& d
);
#endif

void eig3Field
(
    const surfaceSymmTensorField& A,
    surfaceTensorField& V,
    surfaceVectorField& d
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
