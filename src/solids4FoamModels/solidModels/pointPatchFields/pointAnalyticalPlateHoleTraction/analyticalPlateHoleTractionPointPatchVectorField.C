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

#include "analyticalPlateHoleTractionPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "pointPatchFields.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#include "plateHoleAnalyticalFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

symmTensor analyticalPlateHoleTractionPointPatchVectorField::plateHoleSolution
(
    const vector& C
)
{
    return plateHoleAnalyticalFields::stress(C, T_, holeR_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticalPlateHoleTractionPointPatchVectorField::
analyticalPlateHoleTractionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    solidTractionPointPatchVectorField(p, iF),
    T_(0.0),
    holeR_(0.0)
{}


analyticalPlateHoleTractionPointPatchVectorField::analyticalPlateHoleTractionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    solidTractionPointPatchVectorField(p, iF),
    T_(readScalar(dict.lookup("farFieldTractionX"))),
    holeR_(readScalar(dict.lookup("holeRadius")))
{
    traction() = vector::zero;
    pressure() = 0.0;

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


analyticalPlateHoleTractionPointPatchVectorField::analyticalPlateHoleTractionPointPatchVectorField
(
    const analyticalPlateHoleTractionPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    solidTractionPointPatchVectorField(p, iF),
    T_(ptf.T_),
    holeR_(ptf.holeR_)
{}


#ifndef OPENFOAM_ORG
analyticalPlateHoleTractionPointPatchVectorField::analyticalPlateHoleTractionPointPatchVectorField
(
    const analyticalPlateHoleTractionPointPatchVectorField& ptf
)
:
    solidTractionPointPatchVectorField(ptf),
    T_(ptf.T_),
    holeR_(ptf.holeR_)
{}
#endif


analyticalPlateHoleTractionPointPatchVectorField::analyticalPlateHoleTractionPointPatchVectorField
(
    const analyticalPlateHoleTractionPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    solidTractionPointPatchVectorField(ptf, iF),
    T_(ptf.T_),
    holeR_(ptf.holeR_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
void analyticalPlateHoleTractionPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    solidTractionPointPatchVectorField::autoMap(m);
}


// Grab the values using rmap
void analyticalPlateHoleTractionPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidTractionPointPatchVectorField::rmap(ptf, addr);
}


void analyticalPlateHoleTractionPointPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    // Patch point normals
    const vectorField& n = patch().pointNormals();

    // Patch points
    const vectorField& p = patch().localPoints();

    // Set the patch point traction

    vectorField& trac = traction();

    forAll(trac, pointI)
    {
        vector curP(p[pointI].x(), p[pointI].y(), 0);
        vector curN = n[pointI];

        if (patch().name() == "hole")
        {
            curP /= mag(curP);
            curP *= holeR_;

            curN = -curP/mag(curP);
        }

        trac[pointI] =
            plateHoleAnalyticalFields::traction(curP, n[pointI], T_, holeR_);

        // Info<< "p = " << curP << ", n = " << curN
        //     << ", trac " << trac[pointI]
        //     << ", s =  " << plateHoleSolution(curP) << endl;
    }

    solidTractionPointPatchVectorField::initEvaluate(commsType);
}


// Write
void analyticalPlateHoleTractionPointPatchVectorField::write(Ostream& os) const
{
    solidTractionPointPatchVectorField::write(os);

    os.writeKeyword("farFieldTractionX")
        << T_ << token::END_STATEMENT << nl;

    os.writeKeyword("holeRadius")
        << holeR_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    analyticalPlateHoleTractionPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
