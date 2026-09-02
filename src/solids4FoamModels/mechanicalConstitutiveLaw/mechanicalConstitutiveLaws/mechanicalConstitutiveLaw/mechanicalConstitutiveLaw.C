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
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    if (!response.hasFourthOrderTangent())
    {
        FatalErrorInFunction
            << "Finite difference tangent requested but response has no storage"
            << exit(FatalError);
    }

    // Relative perturbation floor. An absolute step is not usable across the
    // range of strains encountered in practice: it is not a small perturbation
    // below strains of order 1e-6, and it approaches round-off in the stress
    // above strains of order 1e-2
    const scalar hMin = 1e-8;
    const scalar hRel = 1e-6;

    UIndirectList<symmTensor>& stress = response.stress();
    UIndirectList<mat66>& Cfield = response.fourthOrderTangent();

    const label nIP = stress.size();

    if (nIP == 0)
    {
        return;
    }

    // The unperturbed stress, which the caller has already evaluated into the
    // response. Taking a copy leaves the caller's storage untouched below.
    // Copied element-wise because Field has no constructor from a
    // UIndirectList on every supported fork
    Field<symmTensor> sigma0(nIP);

    // Per-integration-point perturbation, scaled with the local gradient
    Field<scalar> h(nIP);
    forAll(h, ipI)
    {
        sigma0[ipI] = stress[ipI];
        h[ipI] = max(hMin, hRel*mag(kin.gradD()[ipI]));
    }

    // Evaluate the perturbed states into a shadow of the caller's state. The
    // shadow aliases the old-time fields, so each perturbed evaluation starts
    // from the same history - which is what a consistent tangent means - and
    // writes its outputs where they are discarded. Without this the last
    // perturbed evaluation would be left in the current-time fields, to be
    // read by endTimeStep() and committed by the next storeOldTime()
    mechanicalConstitutiveLawState shadow
    (
        state, mechanicalConstitutiveLawState::SHADOW
    );

    // Work storage. The perturbation is applied to every integration point at
    // once and the law evaluated once per Voigt component, rather than six
    // times per integration point: same arithmetic, one sixth of the calls
    Field<tensor> gradDPert(nIP);
    Field<symmTensor> sigmaPert(nIP, symmTensor::zero);

    // Identity addressing, so the work fields can be presented as the
    // UIndirectList views the kinematics and response wrappers expect
    labelList addr(nIP);
    forAll(addr, ipI)
    {
        addr[ipI] = ipI;
    }

    const UIndirectList<tensor> gradDPertView(gradDPert, addr);
    UIndirectList<symmTensor> sigmaPertView(sigmaPert, addr);

    forAll(Cfield, ipI)
    {
        // mat66 is not zero-initialised
        Cfield[ipI].clear();
    }

    // Loop over the Voigt strain components
    for (label j = 0; j < 6; ++j)
    {
        forAll(gradDPert, ipI)
        {
            gradDPert[ipI] = kin.gradDPerturbed(ipI, j, h[ipI]);
        }

        smallStrainMechanicalConstitutiveLawKinematics kinPert
        (
            gradDPertView, kin.gradD0()
        );

        mechanicalConstitutiveLawResponse respPert
        (
            sigmaPertView, tangentRequest::none
        );

        evaluate(kinPert, inputs, shadow, respPert);

        // Column j of the tangent
        forAll(Cfield, ipI)
        {
            mat66& C = Cfield[ipI];
            const symmTensor ds((sigmaPert[ipI] - sigma0[ipI])/h[ipI]);

            for (label i = 0; i < 6; ++i)
            {
                C(i, j) = ds[i];
            }
        }
    }
}


void Foam::mechanicalConstitutiveLaw::finiteDifferenceFourthOrder
(
    const finiteStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    if (!response.hasFourthOrderTangent())
    {
        FatalErrorInFunction
            << "Finite difference tangent requested but response has no storage"
            << exit(FatalError);
    }

    // Relative perturbation floor, as for the small-strain form. The scale
    // here is the displacement gradient F - I rather than F itself, since it
    // is the deformation that sets the useful step, not the identity
    const scalar hMin = 1e-8;
    const scalar hRel = 1e-6;

    UIndirectList<symmTensor>& stress = response.stress();
    UIndirectList<mat66>& Cfield = response.fourthOrderTangent();

    const label nIP = stress.size();

    if (nIP == 0)
    {
        return;
    }

    // The unperturbed stress, which the caller has already evaluated into the
    // response. Copied element-wise because Field has no constructor from a
    // UIndirectList on every supported fork
    Field<symmTensor> sigma0(nIP);

    Field<scalar> h(nIP);
    forAll(h, ipI)
    {
        sigma0[ipI] = stress[ipI];
        h[ipI] = max(hMin, hRel*mag(kin.F()[ipI] - I));
    }

    // Evaluate the perturbed states against a shadow, so each column starts
    // from the same history and the perturbed results are discarded
    mechanicalConstitutiveLawState shadow
    (
        state, mechanicalConstitutiveLawState::SHADOW
    );

    // Work storage. Every integration point is perturbed at once and the law
    // evaluated once per Voigt component, as in the small-strain form
    Field<tensor> FPert(nIP);
    Field<tensor> FinvPert(nIP);
    Field<scalar> JPert(nIP);
    Field<symmTensor> sigmaPert(nIP, symmTensor::zero);

    labelList addr(nIP);
    forAll(addr, ipI)
    {
        addr[ipI] = ipI;
    }

    const UIndirectList<tensor> FPertView(FPert, addr);
    const UIndirectList<tensor> FinvPertView(FinvPert, addr);
    const UIndirectList<scalar> JPertView(JPert, addr);
    UIndirectList<symmTensor> sigmaPertView(sigmaPert, addr);

    forAll(Cfield, ipI)
    {
        // mat66 is not zero-initialised
        Cfield[ipI].clear();
    }

    for (label j = 0; j < 6; ++j)
    {
        forAll(FPert, ipI)
        {
            FPert[ipI] = kin.FPerturbed(ipI, j, h[ipI]);

            // The inverse and determinant must follow the perturbation, or the
            // kinematics handed to the law would be inconsistent
            FinvPert[ipI] = inv(FPert[ipI]);
            JPert[ipI] = det(FPert[ipI]);
        }

        finiteStrainMechanicalConstitutiveLawKinematics kinPert
        (
            FPertView,
            kin.F0(),
            JPertView,
            kin.J0(),
            FinvPertView,
            kin.Finv0()
        );

        mechanicalConstitutiveLawResponse respPert
        (
            sigmaPertView, tangentRequest::none
        );

        evaluate(kinPert, inputs, shadow, respPert);

        forAll(Cfield, ipI)
        {
            mat66& C = Cfield[ipI];
            const symmTensor ds((sigmaPert[ipI] - sigma0[ipI])/h[ipI]);

            for (label i = 0; i < 6; ++i)
            {
                C(i, j) = ds[i];
            }
        }
    }
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
