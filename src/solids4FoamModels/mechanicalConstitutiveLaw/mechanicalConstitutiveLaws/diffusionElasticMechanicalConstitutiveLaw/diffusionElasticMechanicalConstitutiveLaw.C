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

#include "diffusionElasticMechanicalConstitutiveLaw.H"
#include "mat66.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffusionElasticMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        diffusionElasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dictionary
Foam::diffusionElasticMechanicalConstitutiveLaw::linearElasticDict
(
    const dictionary& dict
)
{
    if (dict.found("sigma0"))
    {
        FatalIOErrorInFunction(dict)
            << "The diffusionElastic law does not take an initial stress, "
            << "sigma0" << exit(FatalIOError);
    }

    if (dict.found("diffusivity"))
    {
        WarningInFunction
            << "The diffusivity entry of the diffusionElastic law is "
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

Foam::diffusionElasticMechanicalConstitutiveLaw::
diffusionElasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    linearElasticMechanicalConstitutiveLaw(linearElasticDict(dict))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList
Foam::diffusionElasticMechanicalConstitutiveLaw::requiredScalarInputs() const
{
    return wordList(1, "stiffnessScaleFactor");
}


void Foam::diffusionElasticMechanicalConstitutiveLaw::evaluate
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    const scalarField& d = inputs.getScalar("stiffnessScaleFactor");

    const bool finiteDifferenceTangent =
        response.tangentReq()
     == tangentRequest::fourthOrderFiniteDifference;

    if (finiteDifferenceTangent)
    {
        // The finite difference tangent must be that of the scaled law, so it
        // is calculated below, from the scaled stress
        mechanicalConstitutiveLawResponse stressResponse(response.stress());

        linearElasticMechanicalConstitutiveLaw::evaluate
        (
            kin, inputs, state, stressResponse
        );
    }
    else
    {
        // The tangents of the linear elastic law are linear in the stiffness,
        // so they are scaled with the stress
        linearElasticMechanicalConstitutiveLaw::evaluate
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

        if (response.wantsFourthOrderTangent())
        {
            UIndirectList<mat66>& C = response.fourthOrderTangent();

            forAll(C, i)
            {
                for (label k = 0; k < 36; ++k)
                {
                    C[i].a_[k] *= d[i];
                }
            }
        }
    }

    UIndirectList<symmTensor>& sigma = response.stress();

    forAll(sigma, i)
    {
        sigma[i] *= d[i];
    }

    if (finiteDifferenceTangent)
    {
        finiteDifferenceFourthOrder(kin, inputs, state, response);
    }
}


// ************************************************************************* //
