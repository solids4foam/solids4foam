/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Stefano Zampini
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Namespace
    Foam::petscErrorHandling

Description
    PETSc error handling for OpenFOAM

    This file was copied from the petsc4Foam library.

SourceFiles
    petscErrorHandling.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_module_petscErrorHandling_H
#define Foam_module_petscErrorHandling_H

#include "petscsys.h"  // For PetscErrorCode

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PetscErrorHandling
{

void Handler
(
    PetscErrorCode ierr,
    const char* invocation,
    const char* filename,
    int lineno
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PetscErrorHandling
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Macro to wrap PETSc calls inside OpenFOAM code.
#define AssertPETSc(code)                                                     \
    do                                                                        \
    {                                                                         \
        PetscErrorCode ierr = (code);                                         \
        if (ierr)                                                             \
        {                                                                     \
            Foam::PetscErrorHandling::Handler(ierr,#code,__FILE__,__LINE__);  \
        }                                                                     \
    } while (false)

// Macro to wrap PETSc inside callbacks.
#ifndef PetscCall
#define PetscCall(code)                                                       \
    do                                                                        \
    {                                                                         \
        PetscErrorCode ierr = (code); CHKERRQ(ierr);                          \
    } while (false)
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
