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

#include "diffusionHyperElasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffusionHyperElasticMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        diffusionHyperElasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dictionary
Foam::diffusionHyperElasticMechanicalConstitutiveLaw::neoHookeanElasticDict
(
    const dictionary& dict
)
{
    if (dict.found("diffusivity"))
    {
        WarningInFunction
            << "The diffusivity entry of the diffusionHyperElastic law is "
            << "ignored: the stiffness scale factor is now calculated by the "
            << "solids4foamSolidModel mesh motion solver, from the "
            << "diffusivity entry in its dynamicMeshDict" << endl;
    }

    dictionary lawDict(dict);

    if (!lawDict.found("K") && lawDict.found("kappa"))
    {
        lawDict.add("K", dimensionedScalar(lawDict.lookup("kappa")));
    }

    return lawDict;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionHyperElasticMechanicalConstitutiveLaw::
diffusionHyperElasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    neoHookeanElasticMechanicalConstitutiveLaw(neoHookeanElasticDict(dict))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList
Foam::diffusionHyperElasticMechanicalConstitutiveLaw::
requiredScalarInputs() const
{
    return wordList(1, "stiffnessScaleFactor");
}


void Foam::diffusionHyperElasticMechanicalConstitutiveLaw::evaluate
(
    const finiteStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    const scalarField& d = inputs.getScalar("stiffnessScaleFactor");

    if (response.wantsFourthOrderTangent())
    {
        // The neo-Hookean law has only a finite difference fourth-order
        // tangent, which must be that of the scaled law, so it is calculated
        // below, from the scaled stress
        mechanicalConstitutiveLawResponse stressResponse(response.stress());

        neoHookeanElasticMechanicalConstitutiveLaw::evaluate
        (
            kin, inputs, state, stressResponse
        );
    }
    else
    {
        neoHookeanElasticMechanicalConstitutiveLaw::evaluate
        (
            kin, inputs, state, response
        );

        if (response.wantsScalarTangent())
        {
            UIndirectList<scalar>& K = response.scalarTangent();

            forAll(K, i)
            {
                K[i] *= d[i];
            }
        }
    }

    UIndirectList<symmTensor>& sigma = response.stress();

    forAll(sigma, i)
    {
        sigma[i] *= d[i];
    }

    fourthOrderByFiniteDifferenceOnly(kin, inputs, state, response);
}


// ************************************************************************* //
