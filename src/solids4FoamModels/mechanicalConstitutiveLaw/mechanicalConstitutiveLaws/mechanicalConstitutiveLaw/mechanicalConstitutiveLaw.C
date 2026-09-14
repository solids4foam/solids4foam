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

InClass
    Foam::mechanicalConstitutiveLaw

\*---------------------------------------------------------------------------*/

#include "mechanicalConstitutiveLaw.H"
#include "mat66.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mechanicalConstitutiveLaw, 0);
    defineRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw, mechanicalConstitutiveLaw
    );
}

// * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * * //

void Foam::mechanicalConstitutiveLaw::finiteDifferenceFourthOrder
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    // Not available in this increment.
    //
    // A finite-difference tangent perturbs the kinematics and re-evaluates the
    // law once per Voigt component. For a history-dependent law each of those
    // evaluations must start from the same history and must not leave its
    // outputs in the current-time state fields, so it needs a shadow of the
    // constitutive state that aliases the old-time fields. That shadow does
    // not exist yet, and without it the perturbed evaluations would corrupt
    // the history rather than measure the tangent.
    //
    // Until it does, tangentRequest::fourthOrderFiniteDifference is rejected
    // rather than silently answered with a wrong tangent

    // Silence unused-parameter warnings
    (void)kin;
    (void)state;

    FatalErrorInFunction
        << "A " << tangentRequestName(response.tangentReq())
        << " tangent was requested, but the finite-difference tangent is not "
        << "yet implemented." << nl
        << "It needs a copy-on-write shadow of the constitutive state so that "
        << "each perturbed evaluation starts from the same history; that is "
        << "added with the laws that need it."
        << exit(FatalError);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLaw::mechanicalConstitutiveLaw
(
    const dictionary& dict
)
{
    // Suppress unused-parameter warning
    (void)dict;
}


// ************************************************************************* //
