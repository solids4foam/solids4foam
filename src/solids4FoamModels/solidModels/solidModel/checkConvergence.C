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

\*---------------------------------------------------------------------------*/

#include "checkConvergence.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


convergenceParameters readConvergenceParameters(const dictionary& dict)
{
    convergenceParameters p;
    p.rTol_ = dict.lookupOrDefault<scalar>("rTol", 1e-6);
    p.aTol_ = dict.lookupOrDefault<scalar>("aTol", 1e-50);
    p.sTol_ = dict.lookupOrDefault<scalar>("sTol", 1e-10);
    p.divTol_ = dict.lookupOrDefault<scalar>("divTol", 1e4);
    p.maxIterations_ = dict.lookupOrDefault<label>("nCorrectors", 10000);
    p.writeResidualFrequency_ =
        dict.lookupOrDefault<label>("infoFrequency", 100);
    p.writeConvergedReason_ =
        dict.lookupOrDefault<Switch>("writeConvergedReason", true);

    return p;
}


bool checkConvergence
(
    const scalar currentResidualNorm,
    const scalar initialResidualNorm,
    const scalar deltaXNorm,
    const scalar xNorm,
    const label iteration,
    const convergenceParameters& param
)
{
    return checkConvergence
    (
        currentResidualNorm,
        initialResidualNorm,
        deltaXNorm,
        xNorm,
        iteration,
        param.maxIterations_,
        param.rTol_,
        param.aTol_,
        param.sTol_,
        param.divTol_,
        param.writeResidualFrequency_,
        param.writeConvergedReason_
    );
}


bool checkConvergence
(
    const scalar currentResidualNorm,
    const scalar initialResidualNorm,
    const scalar deltaXNorm,
    const scalar xNorm,
    const label iteration,
    const label maxIterations,
    const scalar rtol,
    const scalar atol,
    const scalar stol,
    const scalar divtol,
    const label writeResidualFrequency,
    const bool writeConvergedReason
)
{
    // Precompute tolerances
    const scalar relativeResidualTol = rtol*initialResidualNorm;
    const scalar stepTolerance = stol*xNorm;

    // Log residuals if enabled
    if (writeResidualFrequency > 0)
    {
        if (iteration == 1)
        {
            // Print the header with fixed widths
            Info<< setw(10) << "Iteration"
                << setw(20) << "Residual Norm"
                << setw(20) << "Step Norm" << endl;
        }

        // Print each iteration's data with aligned fields
        if (iteration % writeResidualFrequency == 0)
        {
            Info<< setw(10) << iteration
                << setw(20) << currentResidualNorm
                << setw(20) << deltaXNorm << endl;
        }
    }

    // 1. Check Residual Norm against absolute tolerance
    if (currentResidualNorm <= atol)
    {
        if (writeConvergedReason)
        {
            Info<< setw(10) << iteration
                << setw(20) << ": Converged - Absolute residual tolerance met."
                << endl;
        }
        return true;
    }

    // 2. Check Residual Norm Convergence
    if (currentResidualNorm <= relativeResidualTol)
    {
        if (writeConvergedReason)
        {
            Info<< setw(10) << iteration
                << setw(20) << ": Converged - Relative residual tolerance met."
                << endl;
        }
        return true;
    }

    // 3. Check Step Norm Convergence
    if (deltaXNorm <= stepTolerance)
    {
        if (writeConvergedReason)
        {
            Info<< setw(10) << iteration
                << setw(20) << ": Converged - Step norm relative tolerance met."
                << endl;
        }
        return true;
    }

    // 4. Check Divergence
    if (currentResidualNorm >= divtol*initialResidualNorm)
    {
        FatalErrorInFunction
            << "Iteration " << iteration
            << setw(20) << ": Diverged - Residual grew excessively."
            << abort(FatalError);
        return false;
    }

    // 5. Check Maximum Iterations
    if (iteration >= maxIterations)
    {
        FatalErrorInFunction
            << "Iteration " << iteration
            << setw(20) << ": Failed - Maximum iterations reached."
            << abort(FatalError);
        return false;
    }

    // 6. Not Converged Yet
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
