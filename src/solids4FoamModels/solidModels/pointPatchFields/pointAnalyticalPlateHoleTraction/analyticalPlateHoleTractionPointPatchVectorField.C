/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

symmTensor analyticalPlateHoleTractionPointPatchVectorField::plateHoleSolution
(
    const vector& C
)
{
    tensor sigma = tensor::zero;

    // Calculate radial coordinate
    const scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));

    // Calculate circumferential coordinate
    const scalar theta = Foam::atan2(C.y(), C.x());

    const coordinateSystem cs("polarCS", C, vector(0, 0, 1), C/mag(C));

    sigma.xx() =
        T_*(1 - sqr(holeR_)/sqr(r))/2
      + T_
       *(1 + 3*pow(holeR_,4)/pow(r,4) - 4*sqr(holeR_)/sqr(r))*::cos(2*theta)/2;

    sigma.xy() =
      - T_
       *(1 - 3*pow(holeR_,4)/pow(r,4) + 2*sqr(holeR_)/sqr(r))*::sin(2*theta)/2;

    sigma.yx() = sigma.xy();

    sigma.yy() =
        T_*(1 + sqr(holeR_)/sqr(r))/2
      - T_*(1 + 3*pow(holeR_,4)/pow(r,4))*::cos(2*theta)/2;


    // Transformation to Cartesian coordinate system
#ifdef OPENFOAMFOUNDATION
    sigma = ((cs.R().R() & sigma) & cs.R().R().T());
#else
    sigma = ((cs.R() & sigma) & cs.R().T());
#endif

    symmTensor S = symmTensor::zero;

    S.xx() = sigma.xx();
    S.xy() = sigma.xy();
    S.yy() = sigma.yy();

    return S;
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


#ifndef OPENFOAMFOUNDATION
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

        trac[pointI] = (n[pointI] & plateHoleSolution(curP));

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
