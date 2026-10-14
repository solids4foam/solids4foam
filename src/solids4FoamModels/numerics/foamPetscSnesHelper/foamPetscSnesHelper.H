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
    Helper class for models that are solved with the PETSc nonlinear solver
    (SNES).

    A model uses this SNES interface by deriving from this class and defining
    the functions of the petscSnesSystem interface, principally:
    1. formResidual
       This should return the residual (imbalance) of the governing equations
       for the given solution vector
    2. formJacobian
       This should return the Jacobian (or an approximation of it) of the
       residual with respect to the solution vector.
       A Newton-Krylov method can work well with a compact stencil
       approximation of the Jacobian, e.g. fvm::laplacian.

    The helper combines two components, which can also be used on their own:
    - a foamPetscMeshAssembly, which numbers the unknowns of the model mesh
      and assembles fvMatrix and other operator coefficients into PETSc
      matrices and vectors;
    - a foamPetscSnesSolver, which holds the PETSc SNES objects and solves the
      system. It is created on first use, so a model whose residual and
      Jacobian are only used as blocks of a larger system, e.g. by the
      monolithic fluid-solid interaction coupling, holds no SNES objects.

Author
    Philip Cardiff, UCD.

SourceFiles
    foamPetscSnesHelper.C
    foamPetscSnesHelperEnums.C

\*---------------------------------------------------------------------------*/

#ifndef foamPetscSnesHelper_H
#define foamPetscSnesHelper_H

#ifdef USE_PETSC

#include "petscSnesSystem.H"
#include "foamPetscMeshAssembly.H"
#include "foamPetscSnesSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class foamPetscSnesHelper Declaration
\*---------------------------------------------------------------------------*/

class foamPetscSnesHelper
:
    public petscSnesSystem
{
public:

        //- Location of the solution unknowns
        typedef foamPetscMeshAssembly::solutionLocation solutionLocation;


private:

    // Private data

        //- Was PETSc initialised, i.e. can this system be solved?
        const Switch initialised_;

        //- Name of the solution field used to look up PETSc options
        const word fieldName_;

        //- Default PETSc options file used when fvSolution does not define
        //  this solver
        const fileName optionsFile_;

        //- Optional explicit location of the fvSolution dictionary to use
        //  when looking up PETSc options. If empty, use the mesh fvSolution
        const fileName fvSolutionLocation_;

        //- Flag to set whether the solver stops when PETSc SNES returns an
        //  error
        const Switch stopOnPetscError_;

        //- Assembly of the mesh discretisation into PETSc objects
        foamPetscMeshAssembly assembly_;

        //- Solver, created on first use
        mutable autoPtr<foamPetscSnesSolver> solverPtr_;


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
        //    initialise          - Optional: Initialise PETSc so that the
        //                          system can be solved with solve()
        //    fvSolutionLocation  - Optional: explicit location of the
        //                          fvSolution dictionary used to look up PETSc
        //                          options. If unset, use the mesh fvSolution
        foamPetscSnesHelper
        (
            const word& fieldName,
            fileName optionsFile,
            const fvMesh& mesh,
            const solutionLocation& location = solutionLocation::CELLS,
            const Switch stopOnPetscError = true,
            const Switch initialise = true,
            const fileName& fvSolutionLocation = fileName()
        );


    //- Destructor
    virtual ~foamPetscSnesHelper();


    // Member Functions

        // Access

            //- Return the assembly component
            const foamPetscMeshAssembly& assembly() const
            {
                return assembly_;
            }

            //- Return the solver component, creating it on first use
            foamPetscSnesSolver& solver();

            //- Non-const access to the PETSc solution vector
            Vec solution()
            {
                return solver().solution();
            }

            //- Non-const access to the PETSc solution backup vector
            Vec solutionBackup()
            {
                return solver().solutionBackup();
            }

            //- Return a reference to the snes object
            SNES& snes()
            {
                return solver().snes();
            }

            //- Return the stopOnPetscError flag
            Switch stopOnPetscError() const
            {
                return stopOnPetscError_;
            }

            //- Return non-const access to the diverged flag
            bool& diverged()
            {
                return solver().diverged();
            }

            //- Const access to the global index object
            const globalIndex& globalCells() const
            {
                return assembly_.globalCells();
            }

            //- Const access to the global point object
            const globalPointIndices& globalPoints() const
            {
                return assembly_.globalPoints();
            }

            //- Return true if the solve function has been called
            bool snesHasRun() const
            {
                return solverPtr_.valid() && solverPtr_->snesHasRun();
            }


        // Edit

            //- Store the solution to the solution backup
            void storeSolutionBackup()
            {
                solver().storeSolutionBackup();
            }

            //- Reset the solution and snes object
            void resetSnes()
            {
                solver().resetSnes();
            }

            //- Reset retry-sensitive SNES/KSP state while preserving the
            //  solution vector
            void resetSnesSolverState()
            {
                if (solverPtr_.valid())
                {
                    solverPtr_->resetSnesSolverState();
                }
            }

            //- Call the PETSc SNES solve function and check its convergence
            //  See foamPetscSnesSolver::solve
            int solve(const bool returnOnSnesError = false)
            {
                return solver().solve(returnOnSnesError);
            }

            // The petscSnesSystem functions initialiseJacobian(Mat&),
            // initialiseSolution(Vec&), formResidual and formJacobian are
            // defined by the derived class
            using petscSnesSystem::initialiseJacobian;
            using petscSnesSystem::initialiseSolution;

            //- Initialise a Jacobian based on the compact stencil of the mesh
            //  with the given blockSize (number of scalar equations per cell)
            //  The mesh must be the mesh of this helper
            label initialiseJacobian
            (
                Mat& jac,
                const fvMesh& mesh,
                const label blockSize,
                const bool createMat = true  // MatCreate or similar should be called
            ) const;

            //- Initialise the solution vector based on the mesh and blockSize
            //  The mesh must be the mesh of this helper
            label initialiseSolution
            (
                Vec& x,
                const fvMesh& mesh,
                const label blockSize,
                const bool createVec = true  // VecCreate or similar should be called
            ) const;


        // Assembly
        // See foamPetscMeshAssembly for the descriptions of these functions

            template<class Type>
            label InsertFvMatrixIntoPETScMatrix
            (
                fvMatrix<Type>& fvM,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns
            ) const
            {
                return assembly_.InsertFvMatrixIntoPETScMatrix
                (
                    fvM, jac, rowOffset, colOffset, nScalarEqns
                );
            }

            label InsertFvmGradIntoPETScMatrix
            (
                const volScalarField& p,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns,
                const scalar scale = -1.0
            ) const
            {
                return assembly_.InsertFvmGradIntoPETScMatrix
                (
                    p, jac, rowOffset, colOffset, nScalarEqns, scale
                );
            }

            label InsertFvcDivGradInterpolateIntoPETScMatrix
            (
                const volScalarField& p,
                const surfaceScalarField& gamma,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const scalar scale = 1.0
            ) const
            {
                return assembly_.InsertFvcDivGradInterpolateIntoPETScMatrix
                (
                    p, gamma, jac, rowOffset, colOffset, scale
                );
            }

            label InsertFvmDivPhiUIntoPETScMatrix
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns,
                const scalar scale = -1.0
            ) const
            {
                return assembly_.InsertFvmDivPhiUIntoPETScMatrix
                (
                    U, phi, jac, rowOffset, colOffset, nScalarEqns, scale
                );
            }

            label InsertFvmDivUIntoPETScMatrix
            (
                const volScalarField& p,
                const volVectorField& U,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns,
                const scalar scale = 1.0
            ) const
            {
                return assembly_.InsertFvmDivUIntoPETScMatrix
                (
                    p, U, jac, rowOffset, colOffset, nScalarEqns, scale
                );
            }

            label InsertFvmGradPGaussIntoPETScMatrix
            (
                const volScalarField& p,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns,
                const scalar scale = 1.0
            ) const
            {
                return assembly_.InsertFvmGradPGaussIntoPETScMatrix
                (
                    p, jac, rowOffset, colOffset, nScalarEqns, scale
                );
            }

            template<class Type>
            void ExtractFieldComponents
            (
                const PetscScalar *x,
                Field<Type>& vf,
                const label offset,
                const label xBlockSize,
                const labelList& compIndices = labelList()
            ) const
            {
                assembly_.ExtractFieldComponents
                (
                    x, vf, offset, xBlockSize, compIndices
                );
            }

            template<class Type>
            void ExtractFieldComponents
            (
                const Vec x,
                Field<Type>& vf,
                const label offset,
                const labelList& compIndices = labelList()
            ) const
            {
                assembly_.ExtractFieldComponents(x, vf, offset, compIndices);
            }

            template<class Type>
            void InsertFieldComponents
            (
                const Field<Type>& vf,
                PetscScalar *x,
                const label offset,
                const label xBlockSize,
                const labelList& compIndices = labelList()
            ) const
            {
                assembly_.InsertFieldComponents
                (
                    vf, x, offset, xBlockSize, compIndices
                );
            }

            template<class Type>
            void InsertFieldComponents
            (
                const Field<Type>& vf,
                Vec x,
                const label offset,
                const labelList& compIndices = labelList()
            ) const
            {
                assembly_.InsertFieldComponents(vf, x, offset, compIndices);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

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
        const Switch /*initialise*/ = true,
        const fileName& /*fvSolutionLocation*/ = fileName()
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
