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

#include "mokCavityVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mokCavityVelocityFvPatchVectorField::mokCavityVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    period_(1),
    linearProfile_(false),
    yBottom_(0),
    yTop_(1)
{}


Foam::mokCavityVelocityFvPatchVectorField::mokCavityVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    period_(readScalar(dict.lookup("period"))),
    linearProfile_(dict.found("yBottom")),
    yBottom_(linearProfile_ ? readScalar(dict.lookup("yBottom")) : 0),
    yTop_(linearProfile_ ? readScalar(dict.lookup("yTop")) : 1)
{
    if (period_ <= 0)
    {
        FatalIOErrorIn
        (
            "mokCavityVelocityFvPatchVectorField::"
            "mokCavityVelocityFvPatchVectorField(...)",
            dict
        )   << "period must be positive, not " << period_
            << exit(FatalIOError);
    }

    if (linearProfile_ && yTop_ <= yBottom_)
    {
        FatalIOErrorIn
        (
            "mokCavityVelocityFvPatchVectorField::"
            "mokCavityVelocityFvPatchVectorField(...)",
            dict
        )   << "yTop (" << yTop_ << ") must be greater than yBottom ("
            << yBottom_ << ")" << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        updateCoeffs();
    }
}


Foam::mokCavityVelocityFvPatchVectorField::mokCavityVelocityFvPatchVectorField
(
    const mokCavityVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    period_(ptf.period_),
    linearProfile_(ptf.linearProfile_),
    yBottom_(ptf.yBottom_),
    yTop_(ptf.yTop_)
{}


#ifndef OPENFOAM_ORG
Foam::mokCavityVelocityFvPatchVectorField::mokCavityVelocityFvPatchVectorField
(
    const mokCavityVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    period_(ptf.period_),
    linearProfile_(ptf.linearProfile_),
    yBottom_(ptf.yBottom_),
    yTop_(ptf.yTop_)
{}
#endif


Foam::mokCavityVelocityFvPatchVectorField::mokCavityVelocityFvPatchVectorField
(
    const mokCavityVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    period_(ptf.period_),
    linearProfile_(ptf.linearProfile_),
    yBottom_(ptf.yBottom_),
    yTop_(ptf.yTop_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mokCavityVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

#ifdef OPENFOAM_NOT_EXTEND
    const scalar pi = constant::mathematical::pi;
#else
    const scalar pi = mathematicalConstant::pi;
#endif

    const scalar t = db().time().value();
    const scalar vBar = 1.0 - Foam::cos(2.0*pi*t/period_);

    if (linearProfile_)
    {
        const scalarField y(patch().Cf().component(vector::Y));

        fvPatchVectorField::operator==
        (
            vBar*(y - yBottom_)/(yTop_ - yBottom_)*vector(1, 0, 0)
        );
    }
    else
    {
        fvPatchVectorField::operator==(vector(vBar, 0, 0));
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::mokCavityVelocityFvPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("period") << period_ << token::END_STATEMENT << nl;

    if (linearProfile_)
    {
        os.writeKeyword("yBottom") << yBottom_ << token::END_STATEMENT << nl;
        os.writeKeyword("yTop") << yTop_ << token::END_STATEMENT << nl;
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        mokCavityVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
