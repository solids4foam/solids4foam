/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "pulsedParabolicVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pulsedParabolicVelocityFvPatchVectorField::
pulsedParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    steadyValue_(0),
    n_(1, 0, 0),
    y_(0, 1, 0),
    t1_(SMALL),
    t2_(SMALL),
    boundBoxMin_(0, 0, 0),
    boundBoxMax_(0, 0, 0)
{}


pulsedParabolicVelocityFvPatchVectorField::
pulsedParabolicVelocityFvPatchVectorField
(
    const pulsedParabolicVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    steadyValue_(ptf.steadyValue_),
    n_(ptf.n_),
    y_(ptf.y_),
    t1_(ptf.t1_),
    t2_(ptf.t2_),
    boundBoxMin_(ptf.boundBoxMin_),
    boundBoxMax_(ptf.boundBoxMax_)
{}


pulsedParabolicVelocityFvPatchVectorField::
pulsedParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    steadyValue_(readScalar(dict.lookup("steadyValue"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y")),
    t1_(readScalar(dict.lookup("t1"))),
    t2_(readScalar(dict.lookup("t2"))),
    boundBoxMin_(dict.lookup("boundBoxMin")),
    boundBoxMax_(dict.lookup("boundBoxMax"))
{
    if (mag(n_) < SMALL || mag(y_) < SMALL)
    {
        FatalErrorIn("pulsedParabolicVelocityFvPatchVectorField(dict)")
            << "n or y given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);

    evaluate();
}


pulsedParabolicVelocityFvPatchVectorField::
pulsedParabolicVelocityFvPatchVectorField
(
    const pulsedParabolicVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    steadyValue_(fcvpvf.steadyValue_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_),
    t1_(fcvpvf.t1_),
    t2_(fcvpvf.t2_),
    boundBoxMin_(fcvpvf.boundBoxMin_),
    boundBoxMax_(fcvpvf.boundBoxMax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pulsedParabolicVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Scale the max value in time

    const scalar t = this->db().time().value();

    const scalar curMaxValue =
        steadyValue_*pow(t, 2)/::sqrt(pow(t1_ - pow(t, 2), 2) + pow(t2_*t, 2));

    // Calculate the parabolic profile
    const vector ctr = 0.5*(boundBoxMax_ + boundBoxMin_);

    const vectorField& c = patch().Cf();

    // Calculate local 1-D coordinate for the parabolic profile
    const scalarField coord(2*((c - ctr) & y_)/((boundBoxMax_ - boundBoxMin_) & y_));

    vectorField::operator=(n_*curMaxValue*(1.0 - sqr(coord)));
}


// Write
void pulsedParabolicVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("steadyValue")
        << steadyValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    os.writeKeyword("t1")
        << t1_ << token::END_STATEMENT << nl;
    os.writeKeyword("t2")
        << t2_ << token::END_STATEMENT << nl;
    os.writeKeyword("boundBoxMin")
        << boundBoxMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("boundBoxMax")
        << boundBoxMax_ << token::END_STATEMENT << nl;
#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, pulsedParabolicVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
