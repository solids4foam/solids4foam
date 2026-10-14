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
    foamPetscSnesSolver

Description
    Solves a petscSnesSystem with the PETSc nonlinear solver (SNES).

    The solver owns the PETSc objects: the SNES object, the solution vector
    and its backup, the Jacobian matrix, and the options database. The options
    are read from the solvers/<fieldName> sub-dictionary of fvSolution when it
    specifies `solver petsc;`, and otherwise from an options file.

    The callbacks given to PETSc call the residual, Jacobian and
    preconditioner functions of the system.

Author
    Philip Cardiff, UCD.

SourceFiles
    foamPetscSnesSolver.C

\*---------------------------------------------------------------------------*/

#ifndef foamPetscSnesSolver_H
#define foamPetscSnesSolver_H

#ifdef USE_PETSC

#include "petscSnesSystem.H"
#include "PetscHandles.H"
#include "fvMesh.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class foamPetscSnesSolver;
class IOdictionary;

//- User data "context" given to the PETSc callback functions
struct foamPetscSnesSolverContext
{
    //- The system being solved
    petscSnesSystem& system_;

    //- The solver
    foamPetscSnesSolver& solver_;

    //- Construct from components
    foamPetscSnesSolverContext
    (
        petscSnesSystem& system,
        foamPetscSnesSolver& solver
    )
    :
        system_(system),
        solver_(solver)
    {}
};


/*---------------------------------------------------------------------------*\
                      Class foamPetscSnesSolver Declaration
\*---------------------------------------------------------------------------*/

class foamPetscSnesSolver
{
    // Private data

        //- The system to be solved
        petscSnesSystem& system_;

        //- Name of the solution field used to look up PETSc options
        const word fieldName_;

        //- Default PETSc options file used when fvSolution does not define
        //  this solver
        const fileName optionsFile_;

        //- Mesh used for the default fvSolution lookup
        const fvMesh& mesh_;

        //- Optional explicit location of the fvSolution dictionary to use
        //  when looking up PETSc options. If empty, use the mesh fvSolution
        const fileName fvSolutionLocation_;

        //- Store the options database
        //  We re-set the options database before each call to solve
        PetscOptions options_;

        //- Flag indicating if the PETSc options database has been populated
        bool optionsLoaded_;

        //- Flag to set whether the solver stops when PETSc SNES returns an
        //  error
        const Switch stopOnPetscError_;

        //- Flag to track if SNES diverged
        bool diverged_;

        //- Force stale lagged Jacobian/preconditioner state to rebuild
        bool forceSnesSolverStateRebuild_;

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
        autoPtr<foamPetscSnesSolverContext> contextPtr_;

        //- A flag indicating if the solve function has been called
        bool snesHasRun_;


    // Private Member Functions

        //- Return an explicit fvSolution dictionary override if requested
        autoPtr<IOdictionary> makeSolutionDictOverride() const;

        //- Populate the PETSc options database from the given fvSolution dict
        bool populateOptionsFromDict
        (
            const dictionary& solutionDict,
            const fileName& solutionDictPath
        );

        //- Populate the PETSc options database on first use
        void loadOptions();

        //- Initialise the snes object
        label initialiseSnes();

        //- Disallow default bitwise copy construct
        foamPetscSnesSolver(const foamPetscSnesSolver&);

        //- Disallow default bitwise assignment
        void operator=(const foamPetscSnesSolver&);


public:

    //- Runtime type information
    TypeName("foamPetscSnesSolver");


    // Static Member Functions

        //- Initialise PETSc, if it is not already initialised
        static void initialisePetsc();


    // Constructors

        //- Construct from components
        //
        //  Parameters:
        //    system              - The system to be solved, which must
        //                          outlive the solver
        //    fieldName           - Name of the solution field used to look up
        //                          PETSc options from fvSolution
        //    optionsFile         - PETSc options file address
        //    mesh                - Used for the fvSolution solvers dictionary
        //    stopOnPetscError    - Optional: Give a FatalError when PETSc
        //                          returns an error
        //    fvSolutionLocation  - Optional: explicit location of the
        //                          fvSolution dictionary used to look up PETSc
        //                          options. If unset, use the mesh fvSolution
        foamPetscSnesSolver
        (
            petscSnesSystem& system,
            const word& fieldName,
            const fileName& optionsFile,
            const fvMesh& mesh,
            const Switch stopOnPetscError = true,
            const fileName& fvSolutionLocation = fileName()
        );


    //- Destructor
    virtual ~foamPetscSnesSolver();


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

            //- Non-const access to the PETSc solution backup vector
            Vec solutionBackup()
            {
                if (!xBackup_.v)
                {
                    // Initialise based on x
                    // Note: we call solution() rather than accessing x_.v
                    // directly, as x_ is lazily created by solution(); this
                    // way, the backup vector can be created before the first
                    // call to solve()
                    VecDuplicate(solution(), xBackup_.put());
                }

                return xBackup_.v;
            }

            //- Return the Jacobian matrix assembled by the system
            //  This may differ from the operator used by the SNES object,
            //  e.g. with -snes_mf_operator, or if the system replaces the
            //  preconditioning matrix
            Mat jacobian() const
            {
                return A_.m;
            }

            //- Return a reference to the snes object
            SNES& snes()
            {
                return snes_.s;
            }

            //- Return the stopOnPetscError flag
            Switch stopOnPetscError() const
            {
                return stopOnPetscError_;
            }

            //- Return non-const access to the diverged flag
            bool& diverged()
            {
                return diverged_;
            }

            //- Return true if the solve function has been called
            bool snesHasRun() const
            {
                return snesHasRun_;
            }


        // Edit

            //- Store the solution to the solution backup
            void storeSolutionBackup()
            {
                // Note: the solution vector is lazily created, so we take a
                // copy of it before calling solutionBackup(); the evaluation
                // order of function arguments is unspecified in C++, so
                // passing solution() and solutionBackup() directly to VecCopy
                // would not guarantee that x is created first
                Vec x = solution();

                VecCopy(x, solutionBackup());
            }

            //- Rebuild the SNES object, the solution vector and the Jacobian
            //  The system is told through resetCustomSolverState()
            void resetSnes();

            //- Reset retry-sensitive SNES/KSP state while preserving the
            //  solution vector
            void resetSnesSolverState();

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

#endif // #ifdef USE_PETSC

#endif

// ************************************************************************* //
