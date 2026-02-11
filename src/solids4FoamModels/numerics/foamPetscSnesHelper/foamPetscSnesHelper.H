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

#include "PetscHandles.H"
#include "globalIndex.H"
#include "globalPointIndices.H"
#include "volFields.H"
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


namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------* \
                     Class foamPetscSnesHelper Declaration
\*---------------------------------------------------------------------------*/

class foamPetscSnesHelper
{
public:

        //- Location of the solution unknowns
        enum solutionLocation
        {
            CELLS,
            POINTS,
            NONE
        };

        //- Solution location names
#ifdef OPENFOAM_COM
        static const Enum<solutionLocation> solutionLocationNames_;
#else
        static const NamedEnum<solutionLocation, 3> solutionLocationNames_;
#endif

private:

    // Private data

        //- Location of the solution unknowns, e.g. cells or points
        const solutionLocation location_;

        //- Was PETSc initialised?
        const Switch initialised_;

        //- Store the options file
        //  We re-set the options file before each call to solve
        PetscOptions options_;

        //- Flag to set whether the solver stops when PETSc SNES returns an
        //  error
        Switch stopOnPetscError_;

        //- Flag to track if SNES diverged
        bool diverged_;

        //- Pointer to the PETSc SNES object
        //  Note that SNES is a pointer type
        SnesHandle snes_;

        //- Pointer to the PETSc solution vector
        //  Note that Vec is a pointer type
        VecHandle x_;

        //- Pointer to a backup of the PETSc solution vector
        VecHandle xBackup_;

        //- Pointer to the PETSc matrix
        //  Note that Mat is a pointer type
        MatHandle A_;

        //- PETSc user context
        autoPtr<appCtxfoamPetscSnesHelper> snesUserPtr_;

        //- Helper object for converting between local (on processor) and global
        //  cell indices
        autoPtr<globalIndex> globalCellsPtr_;

        //- Helper object for converting between local (on processor) and global
        //  point indices
        autoPtr<globalPointIndices> globalPointsPtr_;

        //- Global cell IDs for the neighbouring processor patch cells
        mutable PtrList<labelList> neiProcGlobalIDs_;

        //- Cell volumes for the neighbouring processor patch cells
        mutable PtrList<scalarField> neiProcVolumes_;

        //- A flag indicating if the solve function has been called
        bool snesHasRun_;


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
        //    solution field name - Name of the solution field used to look up
        //                          PETSc options from fvSolution
        //    optionsFile         - PETSc options file address
        //    mesh                - Used to access the number of cells or points
        //                          and global syncing activities for points,
        //                          and fvSolution solvers dictionary
        //    location            - Optional: Location of the solution unknowns
        //                          This is used to determine the local-to-global
        //                          mapping
        //    stopOnPetscError    - Optional: Give a FatalError when PETSc
        //                          returns an error
        //    initialise          - Optional: Initialise PETSc
        foamPetscSnesHelper
        (
            const word& fieldName,
            fileName optionsFile,
            const fvMesh& mesh,
            const solutionLocation& location = solutionLocation::CELLS,
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
                if (!x_.v)
                {
                    // Initialise snes object
                    if (initialiseSnes() != 0)
                    {
                        FatalErrorInFunction
                        << "initialiseSnes failed" << abort(FatalError);
                    }
                }

                return x_.v;
            }

            //- Non-const access to the PETSc solution vector
            Vec solutionBackup()
            {
                if (!xBackup_.v)
                {
                    // Initialise based on x
                    VecDuplicate(x_.v, xBackup_.put());
                }

                return xBackup_.v;
            }

            //- Return a reference to the snes object
            SNES& snes()
            {
                return snes_.s;
            }

            //- Return the stopOnPetscError flag
            const Switch stopOnPetscError() const
            {
                return stopOnPetscError_;
            }

            //- Return non-const access to the diverged flag
            bool& diverged()
            {
                return diverged_;
            }

            //- Const access to the global index object
            const globalIndex& globalCells() const
            {
                if (globalCellsPtr_.empty())
                {
                    FatalErrorInFunction
                        << "Pointer not set" << exit(FatalError);
                }

                return globalCellsPtr_();
            }

            //- Const access to the global point object
            const globalPointIndices& globalPoints() const
            {
                if (globalPointsPtr_.empty())
                {
                    FatalErrorInFunction
                        << "Pointer not set" << exit(FatalError);
                }

                return globalPointsPtr_();
            }


        // Edit

            //- Store the solution to the solution backup
            void storeSolutionBackup()
            {
                VecCopy(solution(), solutionBackup());
            }

            //- Reset the solution and snes object
            void resetSnes();

            //- Initialise the Jacobian matrix
            //  This procedure should set the size and non-zero structure of
            //  the matrix
            //  We pass a reference ot the matrix pointer so the pointer can be
            //  assigned
            virtual label initialiseJacobian(Mat& jac) = 0;

            //- Initialise Joacobian based on compact stencil of the given mesh with the
            //  given blockSize (number of scalar equations per cell)
            label initialiseJacobian
            (
                Mat& jac,
                const fvMesh& mesh,
                const label blockSize,
                const bool createMat = true  // MatCreate or similar should be called
            );

            //- Initialise the solution vector
            //  This procedure should set the size of the solution vector
            virtual label initialiseSolution(Vec& x) = 0;

            //- Initialise the solution vector based on the given mesh and
            //  blockSize
            label initialiseSolution
            (
                Vec& x,
                const fvMesh& mesh,
                const label blockSize,
                const bool createVec = true  // VecCreate or similar should be called
            );

            //- Calculate the residual of the governing equations
            //  A non-zero return value indicates an error
            virtual label formResidual
            (
                Vec f,         // Residual
                const Vec x    // Solution
            ) = 0;

            //- Form the Jacobian
            //  This Jacobian is the derivative of residual with respect
            //  to the solution vector. This can also be an approximation of
            //  the true Jacobian, e.g. a compact molecule Laplacian. A non-zero
            //  return value indicates an error
            virtual label formJacobian
            (
                Mat jac,       // Jacobian
                const Vec x    // Solution
            ) = 0;

            //- Apply a physics-based preconditioner
            //  This should (approximately) solve y = M^{-1} x, where M is an
            //  approximation of the solid model's Jacobian. This is used inside
            //  the `-pc_type physics` preconditioner.
            //  For example, we may use an approximate segregated solver as a
            //  preconditioner for a Jacobian-free Newton-Krylov solver
            //  A non-zero return value indicates an error
            virtual label precondition
            (
                Vec y,         // Output: y = M^{-1} x
                const Vec x    // Input: vector to be preconditioned
            )
            {
                notImplemented
                (
                    "`label precondition(...)` is not implemented for this "
                    "solid model. The `physics` preconditioner can only be used"
                    " if the chosen solidModel implements this function"
                );

                return 0;
            }

            //- Inserts the OpenFOAM fvMatrix 'fvM' into the PETSc matrix 'jac'.
            //  Converts the finite volume matrix into block coefficients using
            //  the specified row and column offsets (rowOffset, colOffset)
            //  within each block.
            //  nScalarEqns indicates the number of scalar equations to add, e.g.
            //  3 for a 3-D vector equation, 2 for 2-D vector equation, 1 for a
            //  scalar equation
            //  A non-zero return value indicates an error.
            //  Note: fvM should be const but OF9 does not have required const
            //  access functions
            template <class Type>
            label InsertFvMatrixIntoPETScMatrix
            (
                fvMatrix<Type>& fvM,
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
            //  Do NOT use this function if x has a non-constant block size
            //  Use the version of this function above, where the block size
            //  is explicitly passed
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
            //  Do NOT use this function if x has a non-constant block size
            //  Use the version of this function above, where the block size
            //  is explicitly passed
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

            //- Return true if the solve function has been called
            bool snesHasRun() const
            {
                return snesHasRun_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "foamPetscSnesHelperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#else  // !USE_PETSC

// Minimal opaque PETSc types so signatures compile without PETSc headers
typedef struct _p_Mat*  Mat;
typedef struct _p_Vec*  Vec;
typedef struct _p_SNES* SNES;

#include "Switch.H"
#include "word.H"
#include "fileName.H"
#include "fvMesh.H"
#include "List.H"
#include "label.H"

// Keep it in the Foam namespace like the real one
namespace Foam
{

class foamPetscSnesHelper
{
public:
    // Same enum so downstream code compiles
    enum solutionLocation { CELLS, POINTS, NONE };

#ifdef OPENFOAM_COM
    // If you rely on solutionLocationNames_, you can still declare it.
    // Provide a definition in a .C for the no-PETSc build if needed,
    // or omit it if not referenced.
    static const Enum<solutionLocation> solutionLocationNames_;
#else
    static const NamedEnum<solutionLocation, 3> solutionLocationNames_;
#endif

    // Constructor with identical signature: does nothing in no-PETSc builds
    foamPetscSnesHelper
    (
        const word& /*fieldName*/,
        fileName /*optionsFile*/,
        const fvMesh& /*mesh*/,
        const solutionLocation& /*location*/ = solutionLocation::CELLS,
        const Switch /*stopOnPetscError*/ = true,
        const Switch /*initialise*/ = true
    )
    {}

    virtual ~foamPetscSnesHelper() {}

    // --- Minimal API surface you might call even when PETSc is off

    // If something accidentally calls solve() without PETSc, you can either:
    //  a) return an error code (0=success in your code? pick a non-zero),
    //  b) or abort loudly.
    int solve(const bool /*returnOnSnesError*/ = false)
    {
        FatalErrorInFunction
            << "PETSc support disabled at compile time (USE_PETSC not defined)"
            << abort(FatalError);
        return -1;
    }

    // If you ever take the SNES by reference, provide a stub accessor
    // to keep overload resolution happy (don’t call it!)
    SNES& snes()
    {
        FatalErrorInFunction
            << "PETSc support disabled at compile time"
            << abort(FatalError);
        static SNES dummy = nullptr; // never reached
        return dummy;
    }

    // --- Keep the pure virtuals so derived classes still compile & define them
    virtual label initialiseJacobian(Mat& /*jac*/) { return 0; }
    virtual label initialiseSolution(Vec& /*x*/) { return 0; }
    virtual label formResidual(Vec /*f*/, const Vec /*x*/) { return 0; }
    virtual label formJacobian(Mat /*jac*/, const Vec /*x*/) { return 0; }

    // If you have template helpers that don’t depend on PETSc directly,
    // you can keep them here as well (or omit to keep the stub minimal).
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // namespace Foam

#endif // #ifndef USE_PETSC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
