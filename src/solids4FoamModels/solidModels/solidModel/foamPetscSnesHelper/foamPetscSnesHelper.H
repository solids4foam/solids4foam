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


// * * * * * * * * *  External Helper Functions  * * * * * * * * * * * * * * //

namespace Foam
{

    //- Initialise Joacobian based on compact stencil of the given mesh with the
    //  given blockSize (number of scalar equations per cell)
    label initialiseJacobian
    (
        Mat& jac,
        const fvMesh& mesh,
        const label blockSize,
        const bool createMat = true  // MatCreate or similar should be called
    );

    //- Initialise the solution vector based on the given mesh and
    //  blockSize (number of scalar equations per cell)
    label initialiseSolution
    (
        Vec& x,
        const fvMesh& mesh,
        const label blockSize,
        const bool createVec = true  // VecCreate or similar should be called
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                     Class foamPetscSnesHelper Declaration
\*---------------------------------------------------------------------------*/

class foamPetscSnesHelper
{
    // Private data

        //- Was PETSc initialised?
        const Switch initialised_;

        //- Flag to set whether the solver stops when PETSs SNES returns an
        //  error
        const Switch stopOnPetscError_;

        //- Pointer to the PETSc SNES object
        //  Note that SNES is a pointer type
        SNES snes_;

        //- Pointer to the PETSc solution vector
        //  Note that Vec is a pointer type
        Vec x_;

        //- Pointer to the PETSc matrix
        //  Note that Mat is a pointer type
        Mat A_;

        //- PETSc user context
        autoPtr<appCtxfoamPetscSnesHelper> snesUserPtr_;

        //- Helper object for converting between local (on processor) and global
        //  cell indices
        autoPtr<globalIndex> globalCellsPtr_;

        //- Global cell IDs for the neighbouring processor patch cells
        mutable PtrList<labelList> neiProcGlobalIDs_;

        //- Cell volumes for the neighbouring processor patch cells
        mutable PtrList<scalarField> neiProcVolumes_;


    // Private Member Functions

        //- Make the neighbour processor patch fields for the given mesh
        void makeNeiProcFields(const fvMesh& mesh) const;

        //- Return a reference to neiProcGlobalIDs for the given mesh
        const PtrList<labelList>& neiProcGlobalIDs(const fvMesh& mesh) const;

        //- Return a reference to neiProcGlobalIDs for the given mesh
        const PtrList<scalarField>& neiProcVolumes(const fvMesh& mesh) const;

        //- Return a reference to the least squares vectors
        const leastSquaresS4fVectors& lsVectors(const volScalarField& p) const;

        //- Initialise the snes object
        virtual label initialiseSnes();

        //- Disallow default bitwise copy construct
        foamPetscSnesHelper(const foamPetscSnesHelper&);

        //- Disallow default bitwise assignment
        void operator=(const foamPetscSnesHelper&);


public:

    //- Runtime type information
    TypeName("foamPetscSnesHelper");

    // Constructors

        //- Construct from nBlocks and blockSize
        //
        //  Parameters:
        //    optionsFile         - PETSc options file address
        //    nLocalBlocks        - Optional: Number of blocks (e.g., cells) on
        //                          this processor
        //    stopOnPetscError    - Optional: Give a FatalError when PETSc
        //                          returns an error
        //    initialise          - Optional: Initialise PETSc
        foamPetscSnesHelper
        (
            fileName optionsFile,
            const label nLocalBlocks = -1,
            const Switch stopOnPetscError = true,
            const Switch initialise = true
        );


    // Destructor

        virtual ~foamPetscSnesHelper();


    // Member Functions

        // Access

            //- Non-const access to the PETSc solution vector
            Vec solution()
            {
                if (!x_)
                {
                    // Initialise snes object
                    if (initialiseSnes() != 0)
                    {
                        FatalErrorInFunction
                        << "initialiseSnes failed" << abort(FatalError);
                    }
                }

                return x_;
            }

            //- Const access to the global index object
            const globalIndex& globalCells() const
            {
                return globalCellsPtr_();
            }


        // Edit

            //- Initialise the Jacobian matrix
            //  This procedure should set the size and non-zero structure of
            //  the matrix
            //  We pass a reference ot the matrix pointer so the pointer can be
            //  assigned
            virtual label initialiseJacobian(Mat& jac) = 0;

            //- Initialise the solution vector
            //  This procedure should set the size of the solution vector
            virtual label initialiseSolution(Vec& x) = 0;

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
            //  If "flipSign" is true, -grad(p) is inserted rather than grad(p)
            label InsertFvmGradIntoPETScMatrix
            (
                const volScalarField& p,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns,
                const bool flipSign = true
            ) const;

            //- Insert linearisation of div(phi,U) into a PETSc matrix
            //  where phi = Sf & Uf
            //  An upwind discretisation is assumed
            //  If "flipSign" is true, -div(phi,U) is inserted rather than
            //  div(phi,U)
            label InsertFvmDivPhiUIntoPETScMatrix
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns,
                const bool flipSign = true
            ) const;

            //- Insert div(U) into the scalar p equation of a a PETSc matrix
            //  If "flipSign" is true, -div(phi,U) is inserted rather than
            //  div(phi,U)
            label InsertFvmDivUIntoPETScMatrix
            (
                const volScalarField& p,
                const volVectorField& U,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns,
                const bool flipSign = false
            ) const;

            //- Extract a subset of components a raw scalar array 'x' into the
            //  OpenFOAM field 'vf'. The array 'x' is organized in blocks of
            //  'xBlockSize' scalars per block, and extraction begins at the
            //  specified 'offset'. For each element in 'vf', the function
            //  inserts  the components specified in 'compIndices'. If
            //  'compIndices' is empty, it copies all available components
            //  (i.e., the first pTraits<Type>::nComponents scalars), assuming
            //  that 'vf' is stored contiguously.
            //
            //  Parameters:
            //    vf          - The OpenFOAM field (e.g., Field<scalar>,
            //                  Field<vector>, Field<tensor>) into which the
            //                  components will be inserted.
            //    x           - Pointer to the raw scalar array from where the
            //                  components will be extracted.
            //    offset      - Starting index within each block in 'x' where
            //                  extraction begins.
            //    xBlockSize  - The total number of scalar values stored per
            //                  block in 'x' (i.e., the stride).
            //    compIndices - (Optional) List of component indices specifying
            //                  which scalars from each element of 'vf' should
            //                  be inserted. If empty, all components are copied.
            template <class Type>
            void ExtractFieldComponents
            (
                const PetscScalar *x,
                Field<Type>& vf,
                const label offset,
                const label xBlockSize,
                const labelList& compIndices = labelList()
            ) const;

            //- Same as above, but takes a PETSc "Vec x" instead of a raw
            //  PetscScalar array pointer. The xBlockSize is looked up
            //  from "Vec x"
            template <class Type>
            void ExtractFieldComponents
            (
                const Vec x,
                Field<Type>& vf,
                const label offset,
                const labelList& compIndices = labelList()
            ) const;

            //- Inserts a subset of components from the OpenFOAM field 'vf' into
            //  a raw scalar array 'x'. The array 'x' is organized in blocks of
            //  'blockSize' scalars per cell, and insertion begins at the
            //  specified 'offset'. For each element in 'vf', the function
            //  copies the components specified in 'compIndices'. If
            //  'compIndices' is empty, it copies all available components
            //  (i.e., the first pTraits<Type>::nComponents scalars), assuming
            //  that 'vf' is stored contiguously.
            //
            //  Parameters:
            //    vf          - The OpenFOAM field (e.g., Field<scalar>,
            //                  Field<vector>, Field<tensor>) from which the
            //                  components are extracted.
            //    x           - Pointer to the raw scalar array where the
            //                  components will be inserted.
            //    offset      - Starting index within each block in 'x' where
            //                  insertion begins.
            //    xBlockSize  - The total number of scalar values stored per
            //                  block in 'x' (i.e., the stride).
            //    compIndices - (Optional) List of component indices specifying
            //                  which scalars from each element of 'vf' should
            //                  be copied. If empty, all components are copied.
            template <class Type>
            void InsertFieldComponents
            (
                const Field<Type>& vf,
                PetscScalar *x,
                const label offset,
                const label xBlockSize,
                const labelList& compIndices = labelList()
            ) const;

            //- Same as above, but takes a PETSc "Vec x" instead of a raw
            //  PetscScalar array pointer. The xBlockSize is looked up
            //  from "Vec x"
            template <class Type>
            void InsertFieldComponents
            (
                const Field<Type>& vf,
                Vec x,
                const label offset,
                const labelList& compIndices = labelList()
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
