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

#include "cantileverDisplacementPointPatchVectorField.H"
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

cantileverDisplacementPointPatchVectorField::cantileverDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    analyticalSol_(),
    curTimeIndex_(-1)
{}


cantileverDisplacementPointPatchVectorField::cantileverDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF),
    analyticalSol_(dict),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fixedValuePointPatchVectorField::operator==
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        fixedValuePointPatchVectorField::operator==(vector::zero);
    }
}


cantileverDisplacementPointPatchVectorField::cantileverDisplacementPointPatchVectorField
(
    const cantileverDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(p, iF),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


#ifndef OPENFOAM_ORG
cantileverDisplacementPointPatchVectorField::cantileverDisplacementPointPatchVectorField
(
    const cantileverDisplacementPointPatchVectorField& ptf
)
:
    fixedValuePointPatchVectorField(ptf),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}
#endif

cantileverDisplacementPointPatchVectorField::cantileverDisplacementPointPatchVectorField
(
    const cantileverDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
void cantileverDisplacementPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    fixedValuePointPatchVectorField::autoMap(m);
}


// Grab the values using rmap
void cantileverDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValuePointPatchVectorField::rmap(ptf, addr);
}


void cantileverDisplacementPointPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        curTimeIndex_ = this->db().time().timeIndex();

        // Patch point coordinates
        const vectorField& p = patch().localPoints();

        // Calculate point displacement field
        const vectorField disp(analyticalSol_.displacement(p));

        // Set patch displacement field
        fixedValuePointPatchVectorField::operator==(disp);
    }

    fixedValuePointPatchVectorField::initEvaluate(commsType);
}


void cantileverDisplacementPointPatchVectorField::write(Ostream& os) const
{
    fixedValuePointPatchVectorField::write(os);

    analyticalSol_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    cantileverDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
