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
    foamPetscSnesHelper

Description
    Helper class and functions for interacting with the PETSc nonlinear
    solver (SNES).

    Solid models can use this SNES interface by deriving from this class and
    defining two virtual functions:

    1. residualMomentum
       This should return the momentum equation imbalance (residual) for
       the given displacement field
    2. JacobianMomentum
       This should return the Jacobian (or an approximation of it) of the
       residual with respect to the displacement field.
       A Newton-Krylov method can work well with a compact stencil
       approximation of the Jacobian, e.g. fvm::laplacian.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    foamPetscSnesHelper.C
    foamPetscSnesHelperTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamPetscSnesHelper_H
#define foamPetscSnesHelper_H

#ifdef USE_PETSC

#include <petscsnes.h>
#include <petscksp.h>

#include "globalIndex.H"
#include "volFields.H"
#include "sparseMatrix.H"
#include "fvMatrices.H"


// * * * * * * * * * * * * * * External Structs  * * * * * * * * * * * * * * //

// Forward declaration of the Foam::foamPetscSnesHelper class
namespace Foam
{
    class foamPetscSnesHelper;
}

// User data "context" for PETSc functions
typedef struct appCtxfoamPetscSnesHelper
{
    // Reference to the solid model object
    Foam::foamPetscSnesHelper& solMod_;

    // Constructor
    appCtxfoamPetscSnesHelper
    (
        Foam::foamPetscSnesHelper& solMod
    )
    :
        solMod_(solMod)
    {}
} appCtxfoamPetscSnesHelper;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class foamPetscSnesHelper Declaration
\*---------------------------------------------------------------------------*/

class foamPetscSnesHelper
{
    // Private data

        //- Was PETSc initialised?
        const Switch initialised_;

        //- Reference to the solution field
        volVectorField& solution_;

        //- Helper object for converting between local (on processor) and global
        //  cell indices
        globalIndex globalCells_;

        //- Is this a 2-D case? If so, it is assumed that Z is the empty
        //  direction
        const Switch twoD_;

        //- Flag to set whether the solver stops when PETSs SNES returns an
        //  error
        const Switch stopOnPetscError_;

        //- Pointer to the PETSc SNES object
        autoPtr<SNES> snesPtr_;

        //- Pointer to the PETSc solution vector
        autoPtr<Vec> xPtr_;

        //- Pointer to the PETSc matrix
        autoPtr<Mat> APtr_;

        //- PETSc user context
        autoPtr<appCtxfoamPetscSnesHelper> snesUserPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        foamPetscSnesHelper(const foamPetscSnesHelper&);

        //- Disallow default bitwise assignment
        void operator=(const foamPetscSnesHelper&);


public:

    //- Runtime type information
    TypeName("foamPetscSnesHelper");

    // Constructors

        //- Construct from components
        foamPetscSnesHelper
        (
            fileName optionsFile,       // PETSc options file address
            volVectorField& solution,   // The solution field
            const Switch twoD,          // Is this a 2-D model?
            const Switch stopOnPetscError = true,
            const Switch initialise = true // Will PETSc be used?
        );

    // Destructor

        virtual ~foamPetscSnesHelper();


    // Member Functions

        // Access

            //- Non-const access to the solution field
            volVectorField& solution()
            {
                return solution_;
            }

            //- Const access to the global index object
            const globalIndex& globalCells() const
            {
                return globalCells_;
            }

            //- Const access to the twoD flag
            const Switch& twoD() const
            {
                return twoD_;
            }


        // Edit

            //- Calculate the residual of the governing equations
            //  A non-zero return value indicates an error
            virtual label formResidual
            (
                PetscScalar *f,        // Residual
                const PetscScalar *x   // Solution
            ) = 0;

            //- Form the Jacobian
            //  This Jacobian is the derivative of residual with respect
            //  to the solution vector. This can also be an approximation of
            //  the true Jacobian, e.g. a compact molecule Laplacian. A non-zero
            //  return value indicates an error
            virtual label formJacobian
            (
                Mat jac,              // Jacobian
                const PetscScalar *x  // Solution
            ) = 0;

            //- Insert an fvMatrix into a PETSc matrix
            //  Currently only implemented for fvVectorMatrix
            //  A non-zero return value indicates an error
            label insertFvMatrixIntoPETScMatrix
            (
                const fvVectorMatrix& fvM,
                Mat jac
            ) const;

            //- Extracts a subset of components from a raw scalar array 'x' into
            //  the OpenFOAM field 'vf'. 'x' is organized in blocks of
            //  'xBlockSize' scalars per element, and extraction starts at
            //  'offset'. For each element in 'vf', the function copies the
            //  first pTraits<Type>::nComponents scalars, assuming that 'vf' is
            //  stored contiguously.
            //  Parameters:
            //    x          - Pointer to the raw scalar array.
            //    vf         - The OpenFOAM field (e.g., Field<scalar>,
            //                 Field<vector>, Field<tensor>) to populate.
            //    xBlockSize - Block size (number of scalars per element) in the
            //                 array 'x'.
            //    offset     - Starting index within each block in 'x' to begin
            //                 copying.
            template <class Type>
            void extractFieldComponents
            (
                const PetscScalar *x,
                Field<Type>& vf,
                const label xBlockSize,
                const label offset
            ) const;

            //- Inserts a subset of components from the OpenFOAM field 'vf' into
            //  a raw scalar array 'x'. The array 'x' is organized in blocks of
            //  'xBlockSize' scalars per element, and insertion begins at the
            //  specified 'offset'. For each element in 'vf', the function
            //  copies the first pTraits<Type>::nComponents scalars (e.g., for a
            //  vector, the x, y, z components), assuming that 'vf' is stored
            //  contiguously.
            //  Parameters:
            //    vf         - The OpenFOAM field (e.g., Field<scalar>,
            //                 Field<vector>, Field<tensor>) from which the
            //                 components are extracted.
            //    x          - Pointer to the raw scalar array where the
            //                 components will be inserted.
            //    xBlockSize - Block size (number of scalars per element) in the
            //                 array 'x'.
            //    offset     - Starting index within each block in 'x' where
            //                 insertion of the components begins.
            template <class Type>
            void insertFieldComponents
            (
                const Field<Type>& vf,
                PetscScalar *x,
                const label xBlockSize,
                const label offset
            ) const;

            //- Map the Foam solution volVectorField to the PETSc solution Vec
            //  Note that PETSc defines a Vec as a pointer to a struct type
            //  It is assumed that the size of x is 3*f.size() for 3-D and
            //  2*f.size() for 2-D
            void mapSolutionFoamToPetsc();

            //- Map the PETSc solution Vec to the Foam solution volVectorField
            //  Note that PETSc defines a Vec as a pointer to a struct type
            //  It is assumed that the size of x is 3*f.size() for 3-D and
            //  2*f.size() for 2-D
            void mapSolutionPetscToFoam();

            //- Call the PETSc SNES solve function and check its convergence
            //  Give a warning or fatal error if convergence was not reached
            //  If returnOnSnesError = true, a PETSc error code is returned;
            //  otherwise, a fatal error is thrown within this function when
            //  the SNES solver fails
            int solve(const bool returnOnSnesError = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "foamPetscSnesHelperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // #ifdef USE_PETSC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
