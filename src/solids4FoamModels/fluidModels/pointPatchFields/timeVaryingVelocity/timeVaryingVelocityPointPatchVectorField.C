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

#include "timeVaryingVelocityPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "Time.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingVelocityPointPatchVectorField::
timeVaryingVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    timeSeries_(),
    curTimeIndex_(-1)
{}


timeVaryingVelocityPointPatchVectorField::
timeVaryingVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF, dict),
    timeSeries_(dict),
    curTimeIndex_(-1)
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


timeVaryingVelocityPointPatchVectorField::
timeVaryingVelocityPointPatchVectorField
(
    const timeVaryingVelocityPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
#ifdef OPENFOAMESIORFOUNDATION
    const pointPatchFieldMapper& mapper
#else
    const PointPatchFieldMapper& mapper
#endif
)
:
    fixedValuePointPatchVectorField(ptf, p, iF, mapper),
    timeSeries_(ptf.timeSeries_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


timeVaryingVelocityPointPatchVectorField::
timeVaryingVelocityPointPatchVectorField
(
    const timeVaryingVelocityPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    timeSeries_(ptf.timeSeries_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeVaryingVelocityPointPatchVectorField::autoMap
(
#ifdef OPENFOAMESIORFOUNDATION
    const pointPatchFieldMapper& m
#else
    const PointPatchFieldMapper& m
#endif
)
{
    fixedValuePointPatchVectorField::autoMap(m);
}


void timeVaryingVelocityPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const timeVaryingVelocityPointPatchVectorField& oVptf =
        refCast<const timeVaryingVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchVectorField::rmap(oVptf, addr);
}


void timeVaryingVelocityPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vector velocity = vector::zero;

    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Lookup velocity from interpolation table
        velocity = timeSeries_(db().time().timeOutputValue());
    }

    Field<vector>::operator=(velocity);

    fixedValuePointPatchVectorField::updateCoeffs();
}


void timeVaryingVelocityPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    timeSeries_.write(os);
#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    timeVaryingVelocityPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
