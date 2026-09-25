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
    evaluateResponse: build the response a mechanicalConstitutiveLaw writes
    into and evaluate the law, as every evaluation path of the
    mechanicalConstitutiveLawManager does. A header of its own so that the
    manager's evaluation can be shared between translation units.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawEvaluateResponse_H
#define mechanicalConstitutiveLawEvaluateResponse_H

#include "mechanicalConstitutiveLaw.H"
#include "mechanicalConstitutiveLawInputs.H"
#include "mechanicalConstitutiveLawState.H"
#include "mechanicalConstitutiveLawResponse.H"
#include "mechanicalConstitutiveLawTangentRequest.H"
#include "mat66.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Detaches the standing stress from the inputs however evaluateResponse()
// exits. The volumetric branches there return early, so that the view they
// build outlives the evaluation that uses it, and a detach written at the end
// of the function would be stepped over on exactly those paths - leaving
// inputs holding a pointer to a list that has gone out of scope, for the next
// law to read
class incomingStressGuard
{
    const mechanicalConstitutiveLawInputs& inputs_;

public:

    explicit incomingStressGuard
    (
        const mechanicalConstitutiveLawInputs& inputs
    )
    :
        inputs_(inputs)
    {}

    ~incomingStressGuard()
    {
        inputs_.clearIncomingStress();
    }
};


// True if the law, or any law it wraps, is fully incompressible. Asked of
// the whole tree because a wrapper such as electroMechanicalLaw evaluates
// its passive law directly: an incompressible law inside one would otherwise
// reach a total-stress evaluation unseen
inline bool incompressibleLawTree(const mechanicalConstitutiveLaw& law)
{
    if (law.incompressible())
    {
        return true;
    }

    const wordList childNames(law.childStateNames());

    forAll(childNames, i)
    {
        if (incompressibleLawTree(law.childLaw(childNames[i])))
        {
            return true;
        }
    }

    return false;
}


//- Build the response a law writes into, and evaluate the law.
//
//  Every evaluation in the manager ends in the same three lines: work out
//  which tangent storage was supplied, wrap it and the stress in a response,
//  and call the law. This is those lines, written once.
//
//  What the caller keeps is everything that actually differs between those
//  places: which state to evaluate against - real, shadow or scratch - which
//  points to address, whether to evaluate at all, and what happens to the
//  stress afterwards. Those are not incidental, and folding any of them in
//  here would change results.
//
//  Three things are deliberate. The tangent storage is passed rather than a
//  view, so that a view is built only when one is wanted and a null pointer is
//  never dereferenced: some callers guard on the pointer and some on the
//  request, and both stay correct. The request passed is the *effective* one,
//  which is not always the caller's own: the surface boundary path asks for no
//  tangent whatever was requested elsewhere, because it computes none.
//
//  And the storage is const, which looks wrong for something a law writes
//  into. UIndirectList takes a const reference and casts it away itself, so
//  this is what the callers were already doing when they built the view
//  inline. Taking non-const here would push callers onto boundaryFieldRef and
//  primitiveFieldRef, and those are not the same thing: they call
//  setUpToDate() and storeOldTimes(), so merely reaching for the pointer would
//  snapshot an old time that nothing asked for.
template<class KinematicsType>
void evaluateResponse
(
    const mechanicalConstitutiveLaw& law,
    const KinematicsType& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    UIndirectList<symmTensor>& stressView,
    const labelUList& tangentIDs,
    const UList<scalar>* scalarTangentStore,
    const UList<mat66>* fourthOrderTangentStore,
    const tangentRequest tangentReq,
    const UList<scalar>* volumetricStore = nullptr,
    const bool stressReturned = true
)
{
    // A fully incompressible law has no volumetric response of its own, so
    // only a caller that replaces it - one asking for the split - can use it,
    // and only a tangent that leaves the bulk stiffness out means anything.
    // Checked here because every evaluation, stress or tangent, comes this way
    if (incompressibleLawTree(law))
    {
        // Whatever tangent comes with it: a total stress handed back to the
        // caller is undefined for this law. Only a tangent query, which
        // evaluates into a shadow state and discards the stress, may leave
        // the split out
        const bool totalStress = stressReturned && !volumetricStore;

        const bool bulkTangent =
            tangentReq == tangentRequest::scalar
         || needsFourthOrderTangent
            (
                tangentReq
            );

        if (totalStress || bulkTangent)
        {
            FatalErrorInFunction
                << "The mechanical constitutive law " << law.type()
                << " is fully incompressible (nu = 0.5), and was asked for "
                << (bulkTangent ? "a tangent that includes" : "a total stress")
                << (bulkTangent ? " the bulk stiffness." : ".") << nl << nl
                << "    Its bulk modulus is infinite, so neither exists. An "
                << "incompressible material needs a mixed displacement-"
                << "pressure formulation, which solves for the pressure: for "
                << "example coupledPressureDisplacementSolid, or a solid "
                << "model with solvePressure. Otherwise set nu below 0.5."
                << exit(FatalError);
        }
    }

    // One place, because every law is reached through here.
    //
    // A law that needs the stress already standing at its points is handed it
    // as a declared input. Every other law has its output cleared first, so
    // reading it yields zero rather than whatever the last caller left - which
    // turns a silent dependence on buffer contents into an obvious one
    List<symmTensor> incoming;

    // Armed before anything can return, and for every law: detaching a stress
    // that was never attached is a no-op
    const incomingStressGuard guard(inputs);

    if (law.requiresIncomingStress())
    {
        // Copied in the law's own loop order, so that it indexes the same way
        // as the stress it is writing
        incoming.setSize(stressView.size());

        forAll(stressView, i)
        {
            incoming[i] = stressView[i];
        }

        inputs.setIncomingStress(incoming);
    }

    forAll(stressView, i)
    {
        stressView[i] = symmTensor::zero;
    }

    if
    (
        scalarTangentStore
     && needsScalarTangent(tangentReq)
    )
    {
        UIndirectList<scalar> tangentView(*scalarTangentStore, tangentIDs);

        mechanicalConstitutiveLawResponse response
        (
            stressView, tangentView, tangentReq
        );

        if (volumetricStore)
        {
            // Declared inside, so that the view outlives the evaluation using
            // it and is not built when it is not wanted
            UIndirectList<scalar> volumetricView(*volumetricStore, tangentIDs);
            response.requestVolumetricSplit(volumetricView);

            law.evaluate(kin, inputs, state, response);

            return;
        }

        law.evaluate(kin, inputs, state, response);
    }
    else if
    (
        fourthOrderTangentStore
     && needsFourthOrderTangent(tangentReq)
    )
    {
        UIndirectList<mat66> tangentView(*fourthOrderTangentStore, tangentIDs);

        mechanicalConstitutiveLawResponse response
        (
            stressView, tangentView, tangentReq
        );

        if (volumetricStore)
        {
            UIndirectList<scalar> volumetricView(*volumetricStore, tangentIDs);
            response.requestVolumetricSplit(volumetricView);

            law.evaluate(kin, inputs, state, response);

            return;
        }

        law.evaluate(kin, inputs, state, response);
    }
    else
    {
        mechanicalConstitutiveLawResponse response(stressView, tangentReq);

        if (volumetricStore)
        {
            UIndirectList<scalar> volumetricView(*volumetricStore, tangentIDs);
            response.requestVolumetricSplit(volumetricView);

            law.evaluate(kin, inputs, state, response);

            return;
        }

        law.evaluate(kin, inputs, state, response);
    }

    // The detach is the guard's, above: it happens on every path out of here,
    // including the early returns in the volumetric branches
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
