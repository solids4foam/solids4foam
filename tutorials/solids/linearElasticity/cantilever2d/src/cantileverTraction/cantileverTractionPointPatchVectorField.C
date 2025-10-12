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

#include "cantileverTractionPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "pointPatchFields.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "Time.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cantileverTractionPointPatchVectorField::cantileverTractionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    solidTractionPointPatchVectorField(p, iF),
    analyticalSol_(),
    curTimeIndex_(-1)
{}


cantileverTractionPointPatchVectorField::cantileverTractionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    solidTractionPointPatchVectorField(p, iF),
    analyticalSol_(dict),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        solidTractionPointPatchVectorField::operator==
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        solidTractionPointPatchVectorField::operator==(vector::zero);
    }
}


cantileverTractionPointPatchVectorField::cantileverTractionPointPatchVectorField
(
    const cantileverTractionPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    solidTractionPointPatchVectorField(p, iF),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


#ifndef OPENFOAM_ORG
cantileverTractionPointPatchVectorField::cantileverTractionPointPatchVectorField
(
    const cantileverTractionPointPatchVectorField& ptf
)
:
    solidTractionPointPatchVectorField(ptf),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}
#endif


cantileverTractionPointPatchVectorField::cantileverTractionPointPatchVectorField
(
    const cantileverTractionPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    solidTractionPointPatchVectorField(ptf, iF),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
void cantileverTractionPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    solidTractionPointPatchVectorField::autoMap(m);
}


// Grab the values using rmap
void cantileverTractionPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidTractionPointPatchVectorField::rmap(ptf, addr);
}


void cantileverTractionPointPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        curTimeIndex_ = this->db().time().timeIndex();

        // Patch point coordinates
        const vectorField& p = patch().localPoints();

        // Calculate point traction field
        const vectorField trac(analyticalSol_.traction(p));

        // Set patch point traction field
        traction() = trac;
        pressure() = 0.0;
    }

    solidTractionPointPatchVectorField::initEvaluate(commsType);
}


void cantileverTractionPointPatchVectorField::write(Ostream& os) const
{
    solidTractionPointPatchVectorField::write(os);

    analyticalSol_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    cantileverTractionPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
