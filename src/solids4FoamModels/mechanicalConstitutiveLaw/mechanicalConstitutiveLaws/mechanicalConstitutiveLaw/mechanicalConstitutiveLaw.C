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
    if (!response.hasFourthOrderTangent())
    {
        FatalErrorInFunction
            << "Finite difference tangent requested but response has no storage"
            << exit(FatalError);
    }

    const scalar h = 1e-8;  // configurable later

    auto& stress = response.stress();
    auto& Cfield = response.fourthOrderTangent();

    forAll(stress, ipI)
    {
        mat66& C = Cfield[ipI];
        C.clear();   // mat66 is not zero-initialised

        // Copy of the base stress
        symmTensor sigma0 = stress[ipI];

        // Single-entry addressing
        labelList addr(1, 0);

        // Loop over Voigt strain components
        for (label j = 0; j < 6; ++j)
        {
            // Perturb gradD (copy)
            tensor gPert = kin.gradDPerturbed(ipI, j, h);

            UList<tensor> gradDBuf(&gPert, 1);
            UIndirectList<tensor> gradDView(gradDBuf, addr);

            // Unperturbed gradD0
            // The const_cast is safe as we do not modify g0
            const tensor& g0 = kin.gradD0()[ipI];
            UList<tensor> gradD0Buf(const_cast<tensor*>(&g0), 1);
            UIndirectList<tensor> gradD0View(gradD0Buf, addr);

            // Output stress
            symmTensor sigmaPert = symmTensor::zero;
            UList<symmTensor> stressBuf(&sigmaPert, 1);
            UIndirectList<symmTensor> stressView(stressBuf, addr);

            // Kinematics (standard constructor)
            smallStrainMechanicalConstitutiveLawKinematics kinPert
            (
                gradDView, gradD0View, kin.dt()
            );

            mechanicalConstitutiveLawResponse respPert
            (
                stressView, tangentRequest::none
            );

            // FIX!!! WE MUST MAKE A COPY OF THE STATE!
            FatalErrorInFunction
                << "In FD material tangent, make copy of state"
                << exit(FatalError);
            evaluate(kinPert, state, respPert);

            // Finite difference column
            const symmTensor ds((sigmaPert - sigma0)/h);
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
