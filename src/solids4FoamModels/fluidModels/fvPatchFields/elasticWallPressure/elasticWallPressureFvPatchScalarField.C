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

#include "elasticWallPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidSolidInterface.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * hsControls  * * * * * * * * * * * * * * * * //

elasticWallPressureFvPatchScalarField::hsControls::hsControls()
:
    model("secant"),
    modelSpecified(false),
    scale(1.0),
    thicknessBlend("tanh"),
    twoSidedHalving(true),
    waveSpeed(-1.0),
    seedModel("thicknessLimited"),
    secantUpdate("iteration"),
    secantFit("minimax"),
    secantSpatial(false),
    secantSmoothing(5),
    secantMinFactor(0.1),
    secantMaxFactor(100.0),
    secantMaxChange(2.0),
    secantMaxIncrease(1.25),
    secantMemory(5),
    divergenceSafeguard(false),
    writeDiagnostics(false)
{}


elasticWallPressureFvPatchScalarField::hsControls::hsControls
(
    const dictionary& dict,
    const scalar constantHs
)
:
    hsControls()
{
    modelSpecified = dict.found("hsModel");

    if (modelSpecified)
    {
        model = word(dict.lookup("hsModel"));
    }
    else if (constantHs >= SMALL || dict.found("hs"))
    {
        model = "constant";
    }

    scale = dict.lookupOrDefault<scalar>("hsScale", scale);
    thicknessBlend =
        dict.lookupOrDefault<word>("thicknessBlend", thicknessBlend);
    twoSidedHalving =
        dict.lookupOrDefault<Switch>("twoSidedHalving", twoSidedHalving);
    waveSpeed = dict.lookupOrDefault<scalar>("waveSpeed", waveSpeed);
    seedModel = dict.lookupOrDefault<word>("seedModel", seedModel);
    secantUpdate = dict.lookupOrDefault<word>("secantUpdate", secantUpdate);
    secantFit = dict.lookupOrDefault<word>("secantFit", secantFit);
    secantSpatial =
        dict.lookupOrDefault<Switch>("secantSpatial", secantSpatial);
    secantSmoothing =
        dict.lookupOrDefault<label>("secantSmoothing", secantSmoothing);
    secantMinFactor =
        dict.lookupOrDefault<scalar>("secantMinFactor", secantMinFactor);
    secantMaxFactor =
        dict.lookupOrDefault<scalar>("secantMaxFactor", secantMaxFactor);
    secantMaxChange =
        dict.lookupOrDefault<scalar>("secantMaxChange", secantMaxChange);
    secantMaxIncrease =
        dict.lookupOrDefault<scalar>("secantMaxIncrease", secantMaxIncrease);
    secantMemory = dict.lookupOrDefault<label>("secantMemory", secantMemory);
    divergenceSafeguard =
        dict.lookupOrDefault<Switch>
        (
            "divergenceSafeguard", divergenceSafeguard
        );
    writeDiagnostics =
        dict.lookupOrDefault<Switch>("writeDiagnostics", writeDiagnostics);

    const wordList models
    (
        IStringStream("(pWaveSpeed constant thicknessLimited secant)")()
    );

    if (findIndex(models, model) == -1)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown hsModel " << model << ". Valid models are "
            << models << exit(FatalIOError);
    }

    if (seedModel == "secant" || findIndex(models, seedModel) == -1)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown seedModel " << seedModel << ". Valid models are "
            << "(pWaveSpeed constant thicknessLimited)" << exit(FatalIOError);
    }

    if (thicknessBlend != "tanh" && thicknessBlend != "min")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown thicknessBlend " << thicknessBlend
            << ". Valid options are (tanh min)" << exit(FatalIOError);
    }

    if (secantUpdate != "timeStep" && secantUpdate != "iteration")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown secantUpdate " << secantUpdate
            << ". Valid options are (timeStep iteration)"
            << exit(FatalIOError);
    }

    if (secantFit != "median" && secantFit != "minimax")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown secantFit " << secantFit
            << ". Valid options are (median minimax)"
            << exit(FatalIOError);
    }
}


void elasticWallPressureFvPatchScalarField::hsControls::write
(
    Ostream& os
) const
{
    const hsControls defaults;

    if (modelSpecified)
    {
        os.writeKeyword("hsModel") << model << token::END_STATEMENT << nl;
    }

    if (scale != defaults.scale)
    {
        os.writeKeyword("hsScale") << scale << token::END_STATEMENT << nl;
    }

    if (thicknessBlend != defaults.thicknessBlend)
    {
        os.writeKeyword("thicknessBlend")
            << thicknessBlend << token::END_STATEMENT << nl;
    }

    if (twoSidedHalving != defaults.twoSidedHalving)
    {
        os.writeKeyword("twoSidedHalving")
            << twoSidedHalving << token::END_STATEMENT << nl;
    }

    if (waveSpeed != defaults.waveSpeed)
    {
        os.writeKeyword("waveSpeed")
            << waveSpeed << token::END_STATEMENT << nl;
    }

    if (model == "secant")
    {
        os.writeKeyword("seedModel")
            << seedModel << token::END_STATEMENT << nl;
        os.writeKeyword("secantUpdate")
            << secantUpdate << token::END_STATEMENT << nl;
        os.writeKeyword("secantFit")
            << secantFit << token::END_STATEMENT << nl;
        os.writeKeyword("secantSpatial")
            << secantSpatial << token::END_STATEMENT << nl;
        os.writeKeyword("secantSmoothing")
            << secantSmoothing << token::END_STATEMENT << nl;
        os.writeKeyword("secantMinFactor")
            << secantMinFactor << token::END_STATEMENT << nl;
        os.writeKeyword("secantMaxFactor")
            << secantMaxFactor << token::END_STATEMENT << nl;
        os.writeKeyword("secantMaxChange")
            << secantMaxChange << token::END_STATEMENT << nl;
        os.writeKeyword("secantMaxIncrease")
            << secantMaxIncrease << token::END_STATEMENT << nl;
        os.writeKeyword("secantMemory")
            << secantMemory << token::END_STATEMENT << nl;
    }

    if (divergenceSafeguard != defaults.divergenceSafeguard)
    {
        os.writeKeyword("divergenceSafeguard")
            << divergenceSafeguard << token::END_STATEMENT << nl;
    }

    if (writeDiagnostics != defaults.writeDiagnostics)
    {
        os.writeKeyword("writeDiagnostics")
            << writeDiagnostics << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * * //

const scalarField& elasticWallPressureFvPatchScalarField::rhoSolidHs() const
{
    const scalarField& seed = seedCoeff();

    if (secantScale_.size() != seed.size())
    {
        FatalErrorInFunction
            << "Inconsistent secant scale size on patch " << patch().name()
            << abort(FatalError);
    }

    rhoSolidHsPtr_.reset
    (
        new scalarField(hsControls_.scale*safeguardScale_*secantScale_*seed)
    );

    return rhoSolidHsPtr_();
}


void elasticWallPressureFvPatchScalarField::resetRobinState()
{
    kinematicConsistency_ = false;
    fluidMeshFollowsSolid_ = false;
    solidThicknessPtr_.clear();
    seedCoeffPtr_.clear();
    seedCoeffTimeIndex_ = -1;
    // The learned scales are kept (restart, mapping); only their size is
    // made consistent with the patch
    if (secantScale_.size() != patch().size())
    {
        secantScale_.setSize(patch().size());
        secantScale_ = 1.0;
    }
    secantTimeIndex_ = -1;
    secantIter_ = 0;
    iterPressure_.clear();
    iterSolidAccN_.clear();
    iterFluidAccN_.clear();
    secantLambdas_.clear();
    secantSamples_.clear();
    secantSamplesHist_.clear();
    secantSamplesHistStep_.clear();
    secantNumField_.setSize(patch().size(), 0.0);
    secantNumField_ = 0.0;
    secantDenField_.setSize(patch().size(), 0.0);
    secantDenField_ = 0.0;
    prevDeltaPNorm_ = -1;
    nDiverging_ = 0;
}


void elasticWallPressureFvPatchScalarField::reportRobinNumber
(
    const scalarField& robinNumber
)
{
    #ifdef OPENFOAM_NOT_EXTEND
    const fvMesh& mesh = internalField().mesh();
    #else
    const fvMesh& mesh = dimensionedInternalField().mesh();
    #endif

    if (robinNumberTimeIndex_ != mesh.time().timeIndex())
    {
        Info<< type() << " " << patch().name() << ": Robin number: min = "
            << gMin(robinNumber) << ", max = " << gMax(robinNumber)
            << ", mean = " << gAverage(robinNumber) << endl;

        robinNumberTimeIndex_ = mesh.time().timeIndex();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    robinFvPatchScalarField(p, iF),
    prevPressure_(p.patch().size(), 0),
    prevAcceleration_(p.patch().size(), vector::zero),
    rhoSolidHsPtr_(),
    robinNumberTimeIndex_(-1),
    constantHs_(-1.0),
    hsControls_(),
    hsUserPtr_()
{
    safeguardScale_ = 1.0;
    resetRobinState();
}


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const elasticWallPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    robinFvPatchScalarField(ptf, p, iF, mapper),
#ifdef OPENFOAM_ORG
    prevPressure_(mapper(ptf.prevPressure_)),
    prevAcceleration_(mapper(ptf.prevAcceleration_)),
#else
    prevPressure_(ptf.prevPressure_, mapper),
    prevAcceleration_(ptf.prevAcceleration_, mapper),
#endif
    rhoSolidHsPtr_(),
    robinNumberTimeIndex_(-1),
    constantHs_(ptf.constantHs_),
    hsControls_(ptf.hsControls_),
    hsUserPtr_()
{
    if (ptf.hsUserPtr_.valid())
    {
#ifdef OPENFOAM_ORG
        hsUserPtr_.reset(new scalarField(mapper(ptf.hsUserPtr_())));
#else
        hsUserPtr_.reset(new scalarField(ptf.hsUserPtr_(), mapper));
#endif
    }

    // Map the learned Robin coefficient scales
#ifdef OPENFOAM_ORG
    secantScale_ = mapper(ptf.secantScale_);
#else
    secantScale_ = scalarField(ptf.secantScale_, mapper);
#endif
    safeguardScale_ = ptf.safeguardScale_;

    resetRobinState();
}


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    robinFvPatchScalarField(p, iF),
    prevPressure_(p.patch().size(), 0),
    prevAcceleration_(p.patch().size(), vector::zero),
    rhoSolidHsPtr_(),
    robinNumberTimeIndex_(-1),
    constantHs_(dict.lookupOrDefault<scalar>("constantHs", -1.0)),
    hsControls_(dict, constantHs_),
    hsUserPtr_()
{
    if (dict.found("value"))
    {
        Field<scalar>::operator=(scalarField("value", dict, p.size()));
    }

    if (dict.found("prevPressure"))
    {
        prevPressure_ = scalarField("prevPressure", dict, p.size());
    }

    if (dict.found("hs"))
    {
        hsUserPtr_.reset(new scalarField("hs", dict, p.size()));
    }

    // Constant models (and constant secant seeds) need a positive thickness
    if
    (
        (
            hsControls_.model == "constant"
         || (
                hsControls_.model == "secant"
             && hsControls_.seedModel == "constant"
            )
        )
     && constantHs_ < SMALL
     && hsUserPtr_.empty()
    )
    {
        FatalIOErrorInFunction(dict)
            << "A constant Robin coefficient model requires constantHs > 0 or"
            << " an hs field" << exit(FatalIOError);
    }

    if (hsUserPtr_.valid() && gMin(hsUserPtr_()) <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "The hs field must be positive" << exit(FatalIOError);
    }

    safeguardScale_ = 1.0;
    resetRobinState();

    // Restart: restore the learned Robin coefficient scales
    if (dict.found("secantScale"))
    {
        secantScale_ = scalarField("secantScale", dict, p.size());
    }
    safeguardScale_ =
        dict.lookupOrDefault<scalar>("safeguardScale", safeguardScale_);

    Info<< type() << " " << patch().name() << ": hsModel = "
        << hsControls_.model;
    if (hsControls_.model == "constant" && constantHs_ >= SMALL)
    {
        Info<< ", constantHs = " << constantHs_;
    }
    if (hsControls_.model == "secant")
    {
        Info<< ", seedModel = " << hsControls_.seedModel
            << ", secantUpdate = " << hsControls_.secantUpdate
            << ", secantFit = " << hsControls_.secantFit
            << ", secantSpatial = " << hsControls_.secantSpatial;
    }
    if (hsControls_.scale != 1.0)
    {
        Info<< ", hsScale = " << hsControls_.scale;
    }
    Info<< endl;

    if (hsControls_.divergenceSafeguard && hsControls_.model == "secant")
    {
        Info<< type() << " " << patch().name() << ": divergenceSafeguard "
            << "is not used with hsModel secant" << endl;
    }

    this->coeff0() = 1.0;
    this->coeff1() = 1.0;
}


#ifndef OPENFOAM_ORG
elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const elasticWallPressureFvPatchScalarField& pivpvf
)
:
    robinFvPatchScalarField(pivpvf),
    prevPressure_(pivpvf.prevPressure_),
    prevAcceleration_(pivpvf.prevAcceleration_),
    rhoSolidHsPtr_(),
    robinNumberTimeIndex_(-1),
    constantHs_(pivpvf.constantHs_),
    hsControls_(pivpvf.hsControls_),
    hsUserPtr_()
{
    if (pivpvf.hsUserPtr_.valid())
    {
        hsUserPtr_.reset(new scalarField(pivpvf.hsUserPtr_()));
    }

    secantScale_ = pivpvf.secantScale_;
    safeguardScale_ = pivpvf.safeguardScale_;
    resetRobinState();
}
#endif


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const elasticWallPressureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    robinFvPatchScalarField(pivpvf, iF),
    prevPressure_(pivpvf.prevPressure_),
    prevAcceleration_(pivpvf.prevAcceleration_),
    rhoSolidHsPtr_(),
    robinNumberTimeIndex_(-1),
    constantHs_(pivpvf.constantHs_),
    hsControls_(pivpvf.hsControls_),
    hsUserPtr_()
{
    if (pivpvf.hsUserPtr_.valid())
    {
        hsUserPtr_.reset(new scalarField(pivpvf.hsUserPtr_()));
    }

    secantScale_ = pivpvf.secantScale_;
    safeguardScale_ = pivpvf.safeguardScale_;
    resetRobinState();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void elasticWallPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    robinFvPatchScalarField::autoMap(m);

#ifdef OPENFOAM_ORG
    m(prevPressure_, prevPressure_);
    m(prevAcceleration_, prevAcceleration_);
    if (hsUserPtr_.valid())
    {
        m(hsUserPtr_(), hsUserPtr_());
    }
    m(secantScale_, secantScale_);
#else
    prevPressure_.autoMap(m);
    prevAcceleration_.autoMap(m);
    if (hsUserPtr_.valid())
    {
        hsUserPtr_().autoMap(m);
    }
    secantScale_.autoMap(m);
#endif

    rhoSolidHsPtr_.clear();
    resetRobinState();
}


void elasticWallPressureFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    robinFvPatchScalarField::rmap(ptf, addr);

    const elasticWallPressureFvPatchScalarField& mptf =
        refCast<const elasticWallPressureFvPatchScalarField>(ptf);

    prevPressure_.rmap(mptf.prevPressure_, addr);
    prevAcceleration_.rmap(mptf.prevAcceleration_, addr);

    if (hsUserPtr_.valid() && mptf.hsUserPtr_.valid())
    {
        hsUserPtr_().rmap(mptf.hsUserPtr_(), addr);
    }

    secantScale_.rmap(mptf.secantScale_, addr);

    rhoSolidHsPtr_.clear();
    resetRobinState();
}

void elasticWallPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

#ifdef OPENFOAM_NOT_EXTEND
    const fvMesh& mesh = internalField().mesh();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
#endif

    // Map the solid density times the solid virtual thickness mapped to the
    // current (fluid) patch
    const scalarField& rhoSolidHs = this->rhoSolidHs();

    // Lookup the pressure dimensions
#ifdef OPENFOAM_NOT_EXTEND
    const word fieldName = internalField().name();
#else
    const word fieldName = dimensionedInternalField().name();
#endif

    const dimensionSet& pDims =
        mesh.lookupObject<volScalarField>(fieldName).dimensions();;

    // The previous acceleration is updated at the end of each
    // FSI iteration in the fluidSolidInterface
    const scalarField prevDdtUn(patch().nf() & prevAcceleration_);

    // Check if a density field is present
    if (fsi().fluid().mesh().foundObject<volScalarField>("rho"))
    {
        const scalarField& rhoFluid =
            patch().lookupPatchField<volScalarField, scalar>("rho");
        // Gravity and surface-tension source of the wall-normal momentum
        // balance, snGrad(p_rgh) = gSrc - rhoFluid*a_n, in pressure-gradient
        // units: phig is the corresponding flux, (gSrc)*rAUf*|Sf|
        const scalarField& phigFlux =
            patch().lookupPatchField<surfaceScalarField, scalar>("phig");
        scalarField phig(patch().size(), 0.0);
        if (db().foundObject<surfaceScalarField>("rAUf"))
        {
            const scalarField& rAUf =
                patch().lookupPatchField<surfaceScalarField, scalar>("rAUf");
            phig = phigFlux/(rAUf*patch().magSf());
        }

        const scalarField c1(rhoSolidHs/rhoFluid);
        reportRobinNumber(c1*patch().deltaCoeffs());

        if (pDims == dimPressure/dimDensity)
        {
            // p/rho
            // Divide RHS by rhoFluid
            this->coeff0() = 1.0;
            this->coeff1() = c1;
            this->rhs() =
                (prevPressure_ - rhoSolidHs*prevDdtUn + c1*phig)/rhoFluid;
        }
        else
        {
            // p
            this->coeff0() = 1.0;
            this->coeff1() = c1;
            this->rhs() = prevPressure_ - rhoSolidHs*prevDdtUn + c1*phig;
        }
    }
    else
    {
        if (debug)
        {
            Info<< "Did not find rho: looking up from transportProperties"
                << endl;
        }

        // Fluid properties
        const dictionary& transportProperties =
            db().lookupObject<IOdictionary>("transportProperties");

        // Lookup the density from the transport properties
        const dimensionedScalar rhoFluid
        (
            transportProperties.lookup("rho")
        );

        const scalarField c1(rhoSolidHs/rhoFluid.value());
        reportRobinNumber(c1*patch().deltaCoeffs());

        if (debug)
        {
            Info<< "rhoSolidHs = " << max(rhoSolidHs)
                << ", rhoFluid = " << rhoFluid.value()
                << endl;
        }

        if (pDims == dimPressure/dimDensity)
        {
            // p/rho
            this->coeff0() = 1.0;
            this->coeff1() = c1;
            this->rhs() =
                prevPressure_/rhoFluid.value()
              - rhoSolidHs*prevDdtUn/rhoFluid.value();
        }
        else
        {
            // p
            this->coeff0() = 1.0;
            this->coeff1() = c1;
            this->rhs() = prevPressure_ - rhoSolidHs*prevDdtUn;
        }
    }

    robinFvPatchField<scalar>::updateCoeffs();
}


void elasticWallPressureFvPatchScalarField::patchFlux
(
    GeometricField<scalar, fvsPatchField, surfaceMesh>& flux,
    const fvMatrix<scalar>& matrix
) const
{
    scalarField rAU(patch().size(), 0.0);
    if (db().foundObject<volScalarField>("rAU"))
    {
        rAU = patch().lookupPatchField<volScalarField, scalar>("rAU");
    }
    else
    {
        rAU = patch().lookupPatchField<surfaceScalarField, scalar>("rAUf");
    }

    boundaryFieldRef(flux)[patch().index()] = rAU*snGrad()*patch().magSf();
}


void elasticWallPressureFvPatchScalarField::write(Ostream& os) const
{
    robinFvPatchScalarField::write(os);
#ifdef OPENFOAM_ORG
    writeEntry(os, "prevPressure", prevPressure_);
#else
    prevPressure_.writeEntry("prevPressure", os);
#endif

    if (constantHs_ >= SMALL)
    {
        os.writeKeyword("constantHs")
            << constantHs_ << token::END_STATEMENT << nl;
    }

    if (hsUserPtr_.valid())
    {
#ifdef OPENFOAM_ORG
        writeEntry(os, "hs", hsUserPtr_());
#else
        hsUserPtr_().writeEntry("hs", os);
#endif
    }

    hsControls_.write(os);

    // Learned Robin coefficient scales, restored on restart
    if (hsControls_.model == "secant")
    {
#ifdef OPENFOAM_ORG
        writeEntry(os, "secantScale", secantScale_);
#else
        secantScale_.writeEntry("secantScale", os);
#endif
    }

    if (safeguardScale_ != 1.0)
    {
        os.writeKeyword("safeguardScale")
            << safeguardScale_ << token::END_STATEMENT << nl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    elasticWallPressureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
