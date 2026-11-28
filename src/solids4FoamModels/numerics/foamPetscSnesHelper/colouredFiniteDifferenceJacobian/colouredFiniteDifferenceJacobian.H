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
    Foam::colouredFiniteDifferenceJacobian

Description
    Coloured forward-difference Jacobian of the residual of a
    foamPetscSnesHelper system, for use as (part of) a preconditioning matrix.

    The sparsity pattern is that of the square of a given compact
    (first-neighbour) matrix of the same system, i.e. the two-layer stencil of
    a cell-centred finite-volume residual whose face values use cell
    gradients. PETSc colours the pattern (distance-two colouring, so that no
    two columns of a colour share a row) and differences the residual once per
    colour (MatFDColoring). The PETSc objects and the colouring are created on
    the first build and reused; the values can be rebuilt or reused at each
    request. Works in serial and in parallel. PETSc options of the colouring
    and of the differencing take the prefix "solid_fd_" (e.g.
    -solid_fd_mat_fd_coloring_err).

    The residual must be a function of the unknowns only (as required by
    -snes_mf_operator); the residual is re-evaluated at the unperturbed state
    after differencing, so that the fields of the system are left as they were.

SourceFiles
    colouredFiniteDifferenceJacobian.C

\*---------------------------------------------------------------------------*/

#ifndef colouredFiniteDifferenceJacobian_H
#define colouredFiniteDifferenceJacobian_H

#ifdef USE_PETSC

#include "foamPetscSnesHelper.H"
#include <petscmat.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class colouredFiniteDifferenceJacobian Declaration
\*---------------------------------------------------------------------------*/

class colouredFiniteDifferenceJacobian
{
    // Private Data

        //- Name used in log messages
        const word name_;

        //- The system whose residual is differenced
        foamPetscSnesHelper& system_;

        //- Finite-difference Jacobian with the two-layer pattern
        Mat J_;

        //- PETSc colouring context
        MatFDColoring fdc_;

        //- Number of colours
        label nColours_;

        //- Number of finite-difference builds performed
        label nBuilds_;

        //- Wall-clock time spent building, summed over the builds
        scalar buildTime_;


    // Private Member Functions

        //- Residual callback for MatFDColoring
        static PetscErrorCode residual(void*, Vec x, Vec f, void* ctx);

        //- Create the pattern, colouring and differencing context
        void setUp(Mat compact);

        //- No copy construct
        colouredFiniteDifferenceJacobian
        (
            const colouredFiniteDifferenceJacobian&
        ) = delete;

        //- No copy assignment
        void operator=(const colouredFiniteDifferenceJacobian&) = delete;


public:

    // Constructors

        //- Construct from a name and the system to difference
        colouredFiniteDifferenceJacobian
        (
            const word& name,
            foamPetscSnesHelper& system
        );


    //- Destructor
    ~colouredFiniteDifferenceJacobian();


    // Member Functions

        //- Overwrite the assembled compact matrix A with the finite-difference
        //  Jacobian of the residual at x, building the Jacobian values if
        //  rebuild is true or if they have never been built, and reusing the
        //  stored values otherwise. A must be the (first-neighbour) matrix of
        //  the same system with entries added (ADD_VALUES) or assembled; it
        //  is assembled here and its pattern is enlarged on the first call.
        void apply(Mat A, const Vec x, const bool rebuild);

        //- Relative difference between J v and a forward difference of the
        //  residual along a random direction v at x (a consistency check of
        //  the pattern and values; restores the system state)
        scalar check(const Vec x, const scalar relStep = 1e-7);

        //- Has the Jacobian been built
        bool built() const
        {
            return nBuilds_ > 0;
        }

        //- Number of builds performed
        label nBuilds() const
        {
            return nBuilds_;
        }

        //- Number of colours (residual evaluations per build)
        label nColours() const
        {
            return nColours_;
        }

        //- Wall-clock time spent building, summed over the builds
        scalar buildTime() const
        {
            return buildTime_;
        }

        //- Destroy the PETSc objects so that the next apply sets them up again
        //  (e.g. after a topology or decomposition change)
        void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif // USE_PETSC

#endif

// ************************************************************************* //
