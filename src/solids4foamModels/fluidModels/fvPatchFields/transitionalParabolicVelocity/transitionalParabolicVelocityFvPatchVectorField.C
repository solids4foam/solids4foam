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

#include "transitionalParabolicVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

transitionalParabolicVelocityFvPatchVectorField::transitionalParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(0),
    n_(1, 0, 0),
    y_(0, 1, 0),
    transitionPeriod_(SMALL),
    boundBoxMin_(0, 0, 0),
    boundBoxMax_(0, 0, 0)
{}


transitionalParabolicVelocityFvPatchVectorField::transitionalParabolicVelocityFvPatchVectorField
(
    const transitionalParabolicVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    maxValue_(ptf.maxValue_),
    n_(ptf.n_),
    y_(ptf.y_),
    transitionPeriod_(ptf.transitionPeriod_),
    boundBoxMin_(ptf.boundBoxMin_),
    boundBoxMax_(ptf.boundBoxMax_)
{}


transitionalParabolicVelocityFvPatchVectorField::transitionalParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(readScalar(dict.lookup("maxValue"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y")),
    transitionPeriod_(readScalar(dict.lookup("transitionPeriod"))),
    boundBoxMin_(dict.lookup("boundBoxMin")),
    boundBoxMax_(dict.lookup("boundBoxMax"))
{
    if (mag(n_) < SMALL || mag(y_) < SMALL)
    {
        FatalErrorIn("transitionalParabolicVelocityFvPatchVectorField(dict)")
            << "n or y given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);

    if (transitionPeriod_ < SMALL)
    {
        transitionPeriod_ = SMALL;
    }

    evaluate();
}


transitionalParabolicVelocityFvPatchVectorField::transitionalParabolicVelocityFvPatchVectorField
(
    const transitionalParabolicVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    maxValue_(fcvpvf.maxValue_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_),
    transitionPeriod_(fcvpvf.transitionPeriod_),
    boundBoxMin_(fcvpvf.boundBoxMin_),
    boundBoxMax_(fcvpvf.boundBoxMax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void transitionalParabolicVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar curMaxValue = maxValue_;

    if (this->db().time().value() < transitionPeriod_)
    {
        scalar t = this->db().time().value();
        curMaxValue *= 0.5*(1.0-cos(M_PI*t/transitionPeriod_));
    }

//     // Get range and orientation
//     boundBox bb(patch().patch().localPoints(), true);
    
//     boundBox bb(vector(0, 0, -0.025334), vector(0, 0.41, 0.025334));

    vector ctr = 0.5*(boundBoxMax_ + boundBoxMin_);

    const vectorField& c = patch().Cf();

    // Calculate local 1-D coordinate for the parabolic profile
    scalarField coord = 2*((c - ctr) & y_)/((boundBoxMax_ - boundBoxMin_) & y_);

    vectorField::operator=(n_*curMaxValue*(1.0 - sqr(coord)));
}


// Write
void transitionalParabolicVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("maxValue")
        << maxValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    os.writeKeyword("transitionPeriod")
        << transitionPeriod_ << token::END_STATEMENT << nl;
    os.writeKeyword("boundBoxMin")
        << boundBoxMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("boundBoxMax")
        << boundBoxMax_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, transitionalParabolicVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
