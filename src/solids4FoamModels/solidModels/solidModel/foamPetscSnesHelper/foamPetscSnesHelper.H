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
#include "leastSquaresS4fVectors.H"


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
        //volVectorField& solution_;
        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Helper object for converting between local (on processor) and global
        //  cell indices
        globalIndex globalCells_;

        //- Block size of the solution
        //  The solution 'x' is organized in blocks of 'blockSize' scalars per
        //  element
        const label blockSize_;

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

        //- Global cell IDs for the neighbouring processor patch cells
        mutable PtrList<labelList> neiProcGlobalIDs_;

        //- Cell volumes for the neighbouring processor patch cells
        mutable PtrList<scalarField> neiProcVolumes_;

        //- Least squares vectors for the neighbouring processor patch cells
        mutable PtrList<vectorField> neiProcLs_;

        //- Least squares vectors
        //  Note: these are different than the ones used for the displacement
        //  field as different boundary patch face values may be included
        mutable autoPtr<leastSquaresS4fVectors> lsVectorsPtr_;


    // Private Member Functions

        //- Make the neighbour processor patch fields
        void makeNeiProcFields() const;

        //- Return a reference to neiProcGlobalIDs
        const PtrList<labelList>& neiProcGlobalIDs() const;

        //- Return a reference to neiProcGlobalIDs
        const PtrList<scalarField>& neiProcVolumes() const;

        //- Return a reference to neiProcLs_
        const PtrList<vectorField>& neiProcLs() const;

        //- Make least squares vectors
        void makeLsVectors(const volScalarField& p) const;

        //- Return a reference to the least squares vectors
        const leastSquaresS4fVectors& lsVectors(const volScalarField& p) const;

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
            const fvMesh& mesh,         // Mesh
            const label blockSize,      // Number of scalar unknowns per cell
            const labelListList& fieldDefs, // Indicate which rows are different fields
            const Switch stopOnPetscError = true,
            const Switch initialise = true // Will PETSc be used?
        );

    // Destructor

        virtual ~foamPetscSnesHelper();


    // Member Functions

        // Access

            //- Non-const access to the mesh
            const fvMesh& fmesh()
            {
                return mesh_;
            }

            //- Non-const access to the PETSc solution vector
            Vec solution()
            {
                return xPtr_.ref();
            }

            //- Const access to the global index object
            const globalIndex& globalCells() const
            {
                return globalCells_;
            }

            //- Return the blockSize
            label blockSize() const
            {
                return blockSize_;
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

            //- Inserts the OpenFOAM fvMatrix 'fvM' into the PETSc matrix 'jac'.
            //  Converts the finite volume matrix into block coefficients using
            //  the specified row and column offsets (rowOffset, colOffset)
            //  within each block.
            //  nScalarEqns indicates the number of scalar equations to add, e.g.
            //  3 for a 3-D vector equation, 2 for 2-D vector equation, 1 for a
            //  scalar equation
            //  A non-zero return value indicates an error.
            template <class Type>
            label InsertFvMatrixIntoPETScMatrix
            (
                const fvMatrix<Type>& fvM,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns
            ) const;

            //- Inserts least squares gradient into a PETSc matrix at location
            //  (rowOffset, colOffset) within the block
            //  Only implemented for the gradient of a scalar, e.g. grad(p)
            //  Note that the -grad(p) is insert, not grad(p)
            label InsertFvmGradIntoPETScMatrix
            (
                const volScalarField& p,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns
            ) const;

            //- Extracts a subset of components from a raw scalar array 'x' into
            //  the OpenFOAM field 'vf'. 'x' is organized in blocks of
            //  'blockSize_' scalars per element, and extraction starts at
            //  'offset'. For each element in 'vf', the function copies the
            //  first pTraits<Type>::nComponents scalars, assuming that 'vf' is
            //  stored contiguously.
            //  Parameters:
            //    x          - Pointer to the raw scalar array.
            //    vf         - The OpenFOAM field (e.g., Field<scalar>,
            //                 Field<vector>, Field<tensor>) to populate.
            //    offset     - Starting index within each block in 'x' to begin
            //                 copying.
            //    nScalarEqns- Number of components to extract, e.g. 3 for a 3-D
            //                 vector
            template <class Type>
            void ExtractFieldComponents
            (
                const PetscScalar *x,
                Field<Type>& vf,
                const label offset,
                const label nScalarEqns
            ) const;

            //- Same as above, but takes a PETSc "Vec x" instead of a raw
            //  PetscScalar array pointer
            template <class Type>
            void ExtractFieldComponents
            (
                const Vec x,
                Field<Type>& vf,
                const label offset,
                const label nScalarEqns
            ) const;

            //- Inserts a subset of components from the OpenFOAM field 'vf' into
            //  a raw scalar array 'x'. The array 'x' is organized in blocks of
            //  'blockSize_' scalars per element, and insertion begins at the
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
            //    offset     - Starting index within each block in 'x' where
            //                 insertion of the components begins.
            //    nScalarEqns- Number of components to extract, e.g. 3 for a 3-D
            //                 vector
            template <class Type>
            void InsertFieldComponents
            (
                const Field<Type>& vf,
                PetscScalar *x,
                const label offset,
                const label nScalarEqns
            ) const;

            //- Same as above, but takes a PETSc "Vec x" instead of a raw
            //  PetscScalar array pointer
            template <class Type>
            void InsertFieldComponents
            (
                const Field<Type>& vf,
                Vec x,
                const label offset,
                const label nScalarEqns
            ) const;

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
