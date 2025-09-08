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

\*---------------------------------------------------------------------------*/

#include "petscErrorHandling.H"
#include "error.H"

#include <unordered_map>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

static
const std::unordered_map<unsigned int, const char*> human_readable_error =
{
    { PETSC_ERR_MEM              , "unable to allocate requested memory" },
    { PETSC_ERR_SUP              , "no support for requested operation" },
    { PETSC_ERR_SUP_SYS          , "no support for requested operation on this computer system" },
    { PETSC_ERR_ORDER            , "operation done in wrong order" },
    { PETSC_ERR_SIG              , "signal received" },
    { PETSC_ERR_FP               , "floating point exception" },
    { PETSC_ERR_COR              , "corrupted PETSc object" },
    { PETSC_ERR_LIB              , "error in library called by PETSc" },
    { PETSC_ERR_PLIB             , "PETSc library generated inconsistent data" },
    { PETSC_ERR_MEMC             , "memory corruption" },
    { PETSC_ERR_CONV_FAILED      , "iterative method (KSP or SNES) failed" },
    { PETSC_ERR_USER             , "user has not provided needed function" },
    { PETSC_ERR_SYS              , "error in system call" },
    { PETSC_ERR_POINTER          , "pointer does not point to valid address" },
    { PETSC_ERR_MPI_LIB_INCOMP   , "MPI library at runtime is not compatible with MPI user compiled with" },
    { PETSC_ERR_ARG_SIZ          , "nonconforming object sizes used in operation" },
    { PETSC_ERR_ARG_IDN          , "two arguments not allowed to be the same" },
    { PETSC_ERR_ARG_WRONG        , "wrong argument (but object probably ok)" },
    { PETSC_ERR_ARG_CORRUPT      , "null or corrupted PETSc object as argument" },
    { PETSC_ERR_ARG_OUTOFRANGE   , "input argument, out of range" },
    { PETSC_ERR_ARG_BADPTR       , "invalid pointer argument" },
    { PETSC_ERR_ARG_NOTSAMETYPE  , "two args must be same object type" },
    { PETSC_ERR_ARG_NOTSAMECOMM  , "two args must be same communicators" },
    { PETSC_ERR_ARG_WRONGSTATE   , "object in argument is in wrong state, e.g. unassembled mat" },
    { PETSC_ERR_ARG_TYPENOTSET   , "the type of the object has not yet been set" },
    { PETSC_ERR_ARG_INCOMP       , "two arguments are incompatible" },
    { PETSC_ERR_ARG_NULL         , "argument is null that should not be" },
    { PETSC_ERR_ARG_UNKNOWN_TYPE , "type name doesn't match any registered type" },
    { PETSC_ERR_FILE_OPEN        , "unable to open file" },
    { PETSC_ERR_FILE_READ        , "unable to read from file" },
    { PETSC_ERR_FILE_WRITE       , "unable to write to file" },
    { PETSC_ERR_FILE_UNEXPECTED  , "unexpected data in file" },
    { PETSC_ERR_MAT_LU_ZRPVT     , "detected a zero pivot during LU factorization" },
    { PETSC_ERR_MAT_CH_ZRPVT     , "detected a zero pivot during Cholesky factorization" },
    { PETSC_ERR_INT_OVERFLOW     , "integer overflow" },
    { PETSC_ERR_NOT_CONVERGED    , "solver did not converge" },
    { PETSC_ERR_MISSING_FACTOR   , "MatGetFactor() failed" },
    { PETSC_ERR_OPT_OVERWRITE    , "attempted to over write options which should not be changed" },
    { PETSC_ERR_WRONG_MPI_SIZE   , "application run with number of MPI ranks it does not support" },
    { PETSC_ERR_USER_INPUT       , "missing or incorrect user input" },
    { PETSC_ERR_GPU_RESOURCE     , "unable to load a GPU resource, for example cuBLAS" },
    { PETSC_ERR_GPU              , "An error from a GPU call, this may be due to lack of resources on the GPU or a true error in the call" },
    { PETSC_ERR_MPI              , "general MPI error" }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::PetscErrorHandling::Handler
(
    PetscErrorCode ierr,
    const char* invocation,
    const char* filename,
    int lineno
)
{
    const auto iter = human_readable_error.find(ierr);
    const char* message =
    (
        (iter != human_readable_error.end())
      ? iter->second
      : "Unknown error to OpenFOAM"
    );

    FatalError(invocation, filename, lineno)
        << "Error in PETSc. See stacktrace above." << nl
        << "PETSc error code: " << ierr << nl
        << "Error message: \"" << message << '"' << nl
        << abort(FatalError);
}


// ************************************************************************* //
