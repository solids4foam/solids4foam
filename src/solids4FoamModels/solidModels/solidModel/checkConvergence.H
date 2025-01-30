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

Description
    Function for checking convergence.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef checkConvergence_H
#define checkConvergence_H

#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Convergence check
//  Determines whether the nonlinear solver has converged based on
//  various criteria, including residual norms, step norms, and
//  iteration limits.
//  Similar approach to PETSc SNES solver.
//  Arguments:
//  - currentResidualNorm: The 2-norm of the current residual,
//    indicating the deviation of the current solution from equilibrium.
//  - initialResidualNorm: The 2-norm of the initial residual, used to
//    evaluate relative convergence.
//  - deltaXNorm: The 2-norm of the last step (increment in the
//    solution vector) taken by the solver.
//  - xNorm: The 2-norm of the current solution vector, used to scale
//    the step tolerance for relative checks.
//  - iteration: The current iteration count of the nonlinear solver.
//  - maxIterations: The maximum number of allowed iterations before
//    the solver terminates.
//  - rtol: The relative residual norm tolerance. Convergence is
//    achieved if the current residual is less than
//    rtol*initialResidualNorm.
//  - atol: The absolute residual norm tolerance. Convergence is
//    achieved if the current residual is less than atol.
//  - stol: The relative step norm tolerance. Convergence is achieved
//    if deltaXNorm is less than stol*xNorm.
//  - divtol: The divergence tolerance. If the current residual norm
//    grows beyond divtol*initialResidualNorm, the solver is
//    considered diverged.
//  - writeResidualFrequency: Frequency at which the residuals are
//    printed to the console. Set to 0 or a negative number to disable
//    output.
bool checkConvergence
(
    const scalar currentResidualNorm,
    const scalar initialResidualNorm,
    const scalar deltaXNorm,
    const scalar xNorm,
    const int iteration,
    const int maxIterations,
    const scalar rtol,
    const scalar atol,
    const scalar stol,
    const scalar divtol,
    const int writeResidualFrequency,
    const bool writeConvergedReason
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
