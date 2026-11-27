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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar backflowStabilisedOutletVelocityFvPatchVectorField::nu() const
{
    if (nu_ > 0)
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

    if (nu.value() <= 0)
    {
        FatalErrorIn
        (
            "backflowStabilisedOutletVelocityFvPatchVectorField::nu()"
        )   << "nu in transportProperties must be positive for patch "
            << patch().name() << abort(FatalError);
    }

    return nu.value();
}


void backflowStabilisedOutletVelocityFvPatchVectorField::checkLaminar()
{
    if (laminarChecked_)
    {
        return;
    }

    laminarChecked_ = true;

    // The momentum transport model is read from turbulenceProperties
    // (OpenFOAM.com, foam-extend) or momentumTransport (OpenFOAM.org)
    const label nDicts = 2;
    const word dictNames[nDicts] =
    {
        "turbulenceProperties",
        "momentumTransport"
    };

    for (label i = 0; i < nDicts; ++i)
    {
        IOobject io
        (
            dictNames[i],
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

#ifdef OPENFOAM_NOT_EXTEND
        if (!io.typeHeaderOk<IOdictionary>(true))
#else
        if (!io.headerOk())
#endif
        {
            continue;
        }

        const IOdictionary dict(io);

        const word simulationType(dict.lookup("simulationType"));

        if (simulationType != "laminar")
        {
            FatalErrorIn
            (
                "backflowStabilisedOutletVelocityFvPatchVectorField::"
                "checkLaminar()"
            )   << "The condition on patch " << patch().name()
                << " assumes laminar flow with a constant viscosity, but "
                << "simulationType is " << simulationType << " in "
                << dictNames[i] << abort(FatalError);
        }

        // Laminar sub-model, e.g. generalisedNewtonian or Maxwell
        if (dict.found("laminar"))
        {
            const dictionary& laminarDict = dict.subDict("laminar");

            word model("Stokes");

            if (laminarDict.found("model"))
            {
                model = word(laminarDict.lookup("model"));
            }
            else if (laminarDict.found("laminarModel"))
            {
                model = word(laminarDict.lookup("laminarModel"));
            }

            if (model != "Stokes")
            {
                FatalErrorIn
                (
                    "backflowStabilisedOutletVelocityFvPatchVectorField::"
                    "checkLaminar()"
                )   << "The condition on patch " << patch().name()
                    << " assumes a Newtonian fluid with a constant "
                    << "viscosity, but the laminar model is " << model
                    << " in " << dictNames[i] << abort(FatalError);
            }
        }
    }

    // Non-Newtonian viscosity models of transportProperties
    if (db().foundObject<IOdictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
            db().lookupObject<IOdictionary>("transportProperties");

        if (transportProperties.found("transportModel"))
        {
            const word transportModel
            (
                transportProperties.lookup("transportModel")
            );

            if (transportModel != "Newtonian")
            {
                FatalErrorIn
                (
                    "backflowStabilisedOutletVelocityFvPatchVectorField::"
                    "checkLaminar()"
                )   << "The condition on patch " << patch().name()
                    << " assumes a Newtonian fluid with a constant "
                    << "viscosity, but transportModel is " << transportModel
                    << " in transportProperties" << abort(FatalError);
            }
        }
    }
}


scalar backflowStabilisedOutletVelocityFvPatchVectorField::eps
(
    const volVectorField& U
)
{
    if (smoothing_ <= 0)
    {
        return 0;
    }

    // eps is computed on the first evaluation of each time step, which every
    // processor makes, as the time index is the same on all of them, so the
    // global average below is called on all processors
    const label timeIndex = db().time().timeIndex();

    if (epsTimeIndex_ != timeIndex)
    {
        const vectorField UPold
        (
            U.oldTime().boundaryField()[patch().index()].patchInternalField()
        );

        eps_ = smoothing_*max(gAverage(mag(UPold)), SMALL);
        epsTimeIndex_ = timeIndex;
    }

    return eps_;
}


void backflowStabilisedOutletVelocityFvPatchVectorField::updateValueFraction()
{
    // On the first call, once the momentum transport model has been read
    checkLaminar();

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
    const scalar eps = this->eps(U);

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

backflowStabilisedOutletVelocityFvPatchVectorField::
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
    eps_(0),
    epsTimeIndex_(-1),
    laminarChecked_(false)
{
    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = 0;
}


backflowStabilisedOutletVelocityFvPatchVectorField::
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
    eps_(0),
    epsTimeIndex_(-1),
    laminarChecked_(false)
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

    if (dict.found("nu") && nu_ <= 0)
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


backflowStabilisedOutletVelocityFvPatchVectorField::
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
    eps_(0),
    epsTimeIndex_(-1),
    laminarChecked_(false)
{}


#ifndef OPENFOAM_ORG
backflowStabilisedOutletVelocityFvPatchVectorField::
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
    eps_(ptf.eps_),
    epsTimeIndex_(ptf.epsTimeIndex_),
    laminarChecked_(ptf.laminarChecked_)
{}
#endif


backflowStabilisedOutletVelocityFvPatchVectorField::
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
    eps_(ptf.eps_),
    epsTimeIndex_(ptf.epsTimeIndex_),
    laminarChecked_(ptf.laminarChecked_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void backflowStabilisedOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    updateValueFraction();

    mixedFvPatchVectorField::updateCoeffs();
}


void backflowStabilisedOutletVelocityFvPatchVectorField::evaluate
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


void backflowStabilisedOutletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);

    if (UName_ != "U")
    {
        os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    }

    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("smoothing") << smoothing_ << token::END_STATEMENT << nl;

    if (nu_ > 0)
    {
        os.writeKeyword("nu") << nu_ << token::END_STATEMENT << nl;
    }

#ifdef OPENFOAM_ORG
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    backflowStabilisedOutletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
