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

#include "backflowStabilisedOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::backflowStabilisedOutletVelocityFvPatchVectorField::nu()
const
{
    if (nuFromDict_)
    {
        return nu_;
    }

    if (!db().foundObject<IOdictionary>("transportProperties"))
    {
        FatalErrorIn
        (
            "backflowStabilisedOutletVelocityFvPatchVectorField::nu()"
        )   << "Specify nu for patch " << patch().name()
            << ": transportProperties was not found" << abort(FatalError);
    }

    const dictionary& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");

    const dimensionedScalar nu(transportProperties.lookup("nu"));

    return nu.value();
}


void Foam::backflowStabilisedOutletVelocityFvPatchVectorField::
updateValueFraction()
{
    // Velocity of the patch-adjacent cells, from the velocity field rather
    // than from the field this condition belongs to, so that a field with
    // the same boundary types, e.g. a Newton correction, gets the same
    // coefficients
    const volVectorField& U = db().lookupObject<volVectorField>(UName_);
    const label patchID = patch().index();
    const vectorField UP(U.boundaryField()[patchID].patchInternalField());
    const vectorField n(patch().nf());
    const scalarField Un(n & UP);

    // Smoothed min(Un, 0) = 0.5*(Un - |Un|), with |Un| replaced by
    // sqrt(Un^2 + eps^2), so the condition is smooth as faces change
    // between outflow and backflow. The smoothing velocity, eps, is a
    // fraction of the mean speed of the patch-adjacent cells at the
    // previous time step, so it is fixed during the solution of a time step
    const vectorField UPold
    (
        U.oldTime().boundaryField()[patchID].patchInternalField()
    );
    const scalar eps = smoothing_*max(gAverage(mag(UPold)), SMALL);

    // For Un > 0, the algebraically equal form -0.5*eps^2/(Un + sqrt(...))
    // avoids the cancellation of Un - sqrt(...), whose rounding error
    // would spoil finite-difference (matrix-free) Jacobians
    scalarField UnMinus(Un.size(), 0);
    forAll(Un, faceI)
    {
        const scalar r = sqrt(sqr(Un[faceI]) + sqr(eps));

        if (Un[faceI] > 0)
        {
            UnMinus[faceI] = -0.5*sqr(eps)/(Un[faceI] + r);
        }
        else
        {
            UnMinus[faceI] = 0.5*(Un[faceI] - r);
        }
    }

    // Robin coefficient, dU/dn = c*U, non-zero only on backflow faces (up
    // to the smoothing)
    const scalarField c(beta_*UnMinus/nu());

    // Implicit Robin condition as a mixed condition with zero reference
    // value and gradient: U_b = U_P*deltaCoeffs/(deltaCoeffs - c)
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = -c/(deltaCoeffs - c);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::backflowStabilisedOutletVelocityFvPatchVectorField::
backflowStabilisedOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    UName_("U"),
    beta_(0.2),
    smoothing_(0.01),
    nu_(0),
    nuFromDict_(false)
{
    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = 0;
}


Foam::backflowStabilisedOutletVelocityFvPatchVectorField::
backflowStabilisedOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    beta_(dict.lookupOrDefault<scalar>("beta", 0.2)),
    smoothing_(dict.lookupOrDefault<scalar>("smoothing", 0.01)),
    nu_(dict.lookupOrDefault<scalar>("nu", 0)),
    nuFromDict_(dict.found("nu"))
{
    if (beta_ < 0)
    {
        FatalIOErrorIn
        (
            "backflowStabilisedOutletVelocityFvPatchVectorField::"
            "backflowStabilisedOutletVelocityFvPatchVectorField(...)",
            dict
        )   << "beta must be non-negative on patch " << p.name()
            << exit(FatalIOError);
    }

    if (smoothing_ < 0)
    {
        FatalIOErrorIn
        (
            "backflowStabilisedOutletVelocityFvPatchVectorField::"
            "backflowStabilisedOutletVelocityFvPatchVectorField(...)",
            dict
        )   << "smoothing must be non-negative on patch " << p.name()
            << exit(FatalIOError);
    }

    if (nuFromDict_ && nu_ <= 0)
    {
        FatalIOErrorIn
        (
            "backflowStabilisedOutletVelocityFvPatchVectorField::"
            "backflowStabilisedOutletVelocityFvPatchVectorField(...)",
            dict
        )   << "nu must be positive on patch " << p.name()
            << exit(FatalIOError);
    }

    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = 0;

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }
}


Foam::backflowStabilisedOutletVelocityFvPatchVectorField::
backflowStabilisedOutletVelocityFvPatchVectorField
(
    const backflowStabilisedOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    beta_(ptf.beta_),
    smoothing_(ptf.smoothing_),
    nu_(ptf.nu_),
    nuFromDict_(ptf.nuFromDict_)
{}


#ifndef OPENFOAM_ORG
Foam::backflowStabilisedOutletVelocityFvPatchVectorField::
backflowStabilisedOutletVelocityFvPatchVectorField
(
    const backflowStabilisedOutletVelocityFvPatchVectorField& ptf
)
:
    mixedFvPatchVectorField(ptf),
    UName_(ptf.UName_),
    beta_(ptf.beta_),
    smoothing_(ptf.smoothing_),
    nu_(ptf.nu_),
    nuFromDict_(ptf.nuFromDict_)
{}
#endif


Foam::backflowStabilisedOutletVelocityFvPatchVectorField::
backflowStabilisedOutletVelocityFvPatchVectorField
(
    const backflowStabilisedOutletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(ptf, iF),
    UName_(ptf.UName_),
    beta_(ptf.beta_),
    smoothing_(ptf.smoothing_),
    nu_(ptf.nu_),
    nuFromDict_(ptf.nuFromDict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::backflowStabilisedOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    updateValueFraction();

    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::backflowStabilisedOutletVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // The value fraction depends on the velocity, so it is updated at every
    // evaluation: updateCoeffs is skipped while the patch is flagged as
    // updated, e.g. after an fvMatrix of the velocity is assembled without
    // a subsequent evaluation, and a stale value fraction would spoil
    // finite-difference (matrix-free) Jacobians of a Newton solver
    updateValueFraction();

    mixedFvPatchVectorField::evaluate(commsType);
}


void Foam::backflowStabilisedOutletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    mixedFvPatchVectorField::write(os);

    if (UName_ != "U")
    {
        os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    }

    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("smoothing") << smoothing_ << token::END_STATEMENT << nl;

    if (nuFromDict_)
    {
        os.writeKeyword("nu") << nu_ << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        backflowStabilisedOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
