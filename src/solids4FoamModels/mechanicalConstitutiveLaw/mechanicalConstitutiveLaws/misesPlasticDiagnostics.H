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
    The end-of-step diagnostics the von Mises plastic laws share: how many
    integration points yielded over the step, out of how many, and the
    largest increment of equivalent plastic strain. The laws differ in how
    they return to the yield surface, not in what they report.

    Reported, not reduced: see mechanicalConstitutiveLaw::endTimeStep.

SourceFiles
    None: header only

\*---------------------------------------------------------------------------*/

#ifndef misesPlasticDiagnostics_H
#define misesPlasticDiagnostics_H

#include "mechanicalConstitutiveLawState.H"
#include "mechanicalConstitutiveLawDiagnostic.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Append a Mises law's diagnostics for one state, from its equivalent
//  plastic strain at the two times
inline void appendMisesPlasticDiagnostics
(
    const mechanicalConstitutiveLawState& state,
    DynamicList<mechanicalConstitutiveLawDiagnostic>& diagnostics
)
{
    const Field<scalar>& epsilonPEq = state.scalarField("epsilonPEq");
    const Field<scalar>& epsilonPEq0 = state.getScalarField0("epsilonPEq");

    label nYielding = 0;
    scalar curDEpsilonPEq = 0.0;
    scalar maxDEpsilonPEq = 0.0;
    forAll(epsilonPEq, i)
    {
        curDEpsilonPEq = epsilonPEq[i] - epsilonPEq0[i];

        if (curDEpsilonPEq > SMALL)
        {
            ++nYielding;

            maxDEpsilonPEq = max(maxDEpsilonPEq, curDEpsilonPEq);
        }
    }

    typedef mechanicalConstitutiveLawDiagnostic diagnostic;

    diagnostics.append
    (
        diagnostic
        (
            "yielding integration points",
            scalar(nYielding),
            diagnostic::combineOperation::sum
        )
    );

    diagnostics.append
    (
        diagnostic
        (
            "integration points",
            scalar(epsilonPEq.size()),
            diagnostic::combineOperation::sum
        )
    );

    diagnostics.append
    (
        diagnostic
        (
            "max DEpsilonPEq",
            maxDEpsilonPEq,
            diagnostic::combineOperation::maximum
        )
    );
}


//- Say what was found, once the manager has reduced it, under the law's own
//  debug switch. The boundary integration points are counted too, which
//  they were not while each law reduced for itself
inline void reportMisesPlasticDiagnostics
(
    const bool report,
    const UList<mechanicalConstitutiveLawDiagnostic>& diagnostics
)
{
    if (!report || diagnostics.size() != 3)
    {
        return;
    }

    Info<< nl << "Max DEpsilonPEq is " << diagnostics[2].value() << nl
        << "Number of yielding integration points = "
        << label(diagnostics[0].value())
        << "/" << label(diagnostics[1].value())
        << nl << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
