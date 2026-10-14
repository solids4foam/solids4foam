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
    petscSnesSystem

Description
    Abstract interface of a nonlinear system solved with PETSc SNES.

    A class implements this interface to define its system: the sizes and
    sparsity of the Jacobian and solution vector, the residual, the Jacobian
    (or an approximation of it, e.g. a compact stencil Laplacian, for a
    Jacobian-free Newton-Krylov method), and optionally a physics-based
    preconditioner. The system is solved by a foamPetscSnesSolver, which holds
    the PETSc objects and calls these functions.

    The system holds no PETSc state, so a class can define a system without
    solving it itself, e.g. a fluid or solid model whose residual and Jacobian
    are blocks of a monolithic fluid-solid system.

Author
    Philip Cardiff, UCD.

SourceFiles
    petscSnesSystem.H

\*---------------------------------------------------------------------------*/

#ifndef petscSnesSystem_H
#define petscSnesSystem_H

#include "error.H"
#include "label.H"

#ifdef USE_PETSC
    #include <petscsnes.h>
#else
    // Minimal opaque PETSc types so signatures compile without PETSc headers
    typedef struct _p_Mat*  Mat;
    typedef struct _p_Vec*  Vec;
    typedef struct _p_SNES* SNES;
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class petscSnesSystem Declaration
\*---------------------------------------------------------------------------*/

class petscSnesSystem
{
public:

    //- Destructor
    virtual ~petscSnesSystem()
    {}


    // Member Functions

        //- Initialise the Jacobian matrix
        //  This procedure should set the size and non-zero structure of the
        //  matrix. We pass a reference to the matrix pointer so the pointer
        //  can be assigned
        virtual label initialiseJacobian(Mat& jac) = 0;

        //- Initialise the solution vector
        //  This procedure should set the size of the solution vector
        virtual label initialiseSolution(Vec& x) = 0;

        //- Calculate the residual of the governing equations
        //  A non-zero return value indicates an error
        virtual label formResidual
        (
            Vec f,         // Residual
            const Vec x    // Solution
        ) = 0;

        //- Form the Jacobian
        //  This Jacobian is the derivative of residual with respect to the
        //  solution vector. This can also be an approximation of the true
        //  Jacobian, e.g. a compact molecule Laplacian. A non-zero return
        //  value indicates an error
        virtual label formJacobian
        (
            Mat jac,       // Jacobian
            const Vec x    // Solution
        ) = 0;

        //- Apply a physics-based preconditioner
        //  This should (approximately) solve y = M^{-1} x, where M is an
        //  approximation of the system Jacobian. This is used inside the
        //  `-pc_type physics` preconditioner.
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
                "system. The `physics` preconditioner can only be used"
                " if the system implements this function"
            );
            return 0;
        }

        //- Called by the solver after the PETSc options are pushed, before
        //  SNESSolve(), e.g. to adjust the KSP or PC
        virtual void customiseSolver(SNES snes)
        {}

        //- Called by the solver after it has rebuilt the SNES object, so
        //  that state derived from the previous SNES object, e.g. cached
        //  sub-matrices or index sets, can be cleared
        virtual void resetCustomSolverState()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
