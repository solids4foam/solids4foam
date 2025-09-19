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

Namespace
    Foam::PetscUtils

Description
    Utilities for OpenFOAM / PETSc interaction

    This file was copied and modified from the petsc4Foam library.

SourceFiles
    petscUtils.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_module_petscUtils_H
#define Foam_module_petscUtils_H

// OpenFOAM
#include "scalarField.H"

// PETSc
#include "petscvec.h"
#include "petscksp.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace PetscUtils
{

/*---------------------------------------------------------------------------*\
                             Namespace PetscUtils
\*---------------------------------------------------------------------------*/

//- Set PETSc flags.
//  The prefix does not contain a leading '-'
void setFlags
(
    const word& prefix,
    const dictionary& dict,
    const bool verbose = false
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PetscUtils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
