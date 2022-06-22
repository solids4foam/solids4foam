/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solidTractionPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "PointPatchFieldMapper.H"
#include "pointPatchFields.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidTractionPointPatchVectorField::solidTractionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    calculatedPointPatchVectorField(p, iF),
    traction_(p.size()),
    pressure_(p.size()),
    tractionSeries_(),
    pressureSeries_(),
    curTimeIndex_(-1)
{}


solidTractionPointPatchVectorField::solidTractionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    calculatedPointPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_(),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        calculatedPointPatchVectorField::operator==
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        calculatedPointPatchVectorField::operator==(vector::zero);
    }

    // Check how traction is defined
    if (dict.found("traction") && dict.found("tractionSeries"))
    {
        FatalErrorIn
        (
            "solidTractionPointPatchVectorField::"
            "solidTractionPointPatchVectorField"
        )   << "Only traction or tractionSeries can be specified!"
            << abort(FatalError);
    }
    else if (dict.found("tractionSeries"))
    {
        Info<< "    traction is time-varying" << endl;
        tractionSeries_ =
            interpolationTable<vector>(dict.subDict("tractionSeries"));
    }
    else
    {
        traction_ = vectorField("traction", dict, p.size());
    }

    // Check how pressure is defined
    if (dict.found("pressure") && dict.found("pressureSeries"))
    {
        FatalErrorIn
        (
            "solidTractionPointPatchVectorField::"
            "solidTractionPointPatchVectorField"
        )   << "Only pressure or pressureSeries can be specified!"
            << abort(FatalError);
    }
    else if (dict.found("pressureSeries"))
    {
        Info<< "    pressure is time-varying" << endl;
        pressureSeries_ =
            interpolationTable<scalar>(dict.subDict("pressureSeries"));
    }
    else
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }
}


solidTractionPointPatchVectorField::solidTractionPointPatchVectorField
(
    const solidTractionPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    calculatedPointPatchVectorField(p, iF),
    traction_(ptf.traction_, mapper),
    pressure_(ptf.pressure_, mapper),
    tractionSeries_(),
    pressureSeries_(),
    curTimeIndex_(-1)
{}


solidTractionPointPatchVectorField::solidTractionPointPatchVectorField
(
    const solidTractionPointPatchVectorField& ptf
)
:
    calculatedPointPatchVectorField(ptf),
    traction_(ptf.traction_),
    pressure_(ptf.pressure_),
    tractionSeries_(),
    pressureSeries_(),
    curTimeIndex_(-1)
{}


solidTractionPointPatchVectorField::solidTractionPointPatchVectorField
(
    const solidTractionPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    calculatedPointPatchVectorField(ptf),
    traction_(ptf.traction_),
    pressure_(ptf.pressure_),
    tractionSeries_(),
    pressureSeries_(),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
void solidTractionPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    //Field<vector>::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Grab the values using rmap
void solidTractionPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    calculatedPointPatchVectorField::rmap(ptf, addr);

    const solidTractionPointPatchVectorField& tiptf =
      refCast<const solidTractionPointPatchVectorField>(ptf);

    traction_.rmap(tiptf.traction_, addr);
    pressure_.rmap(tiptf.pressure_, addr);
}


void solidTractionPointPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        // If time-varying, update the traction field
        if (tractionSeries_.size())
        {
            traction_ = tractionSeries_(this->db().time().timeOutputValue());
        }

        // If time-varying, update the pressure field
        if (pressureSeries_.size())
        {
            pressure_ = pressureSeries_(this->db().time().timeOutputValue());
        }

        curTimeIndex_ = this->db().time().timeIndex();
    }

    calculatedPointPatchVectorField::initEvaluate(commsType);
}


// Write
void solidTractionPointPatchVectorField::write(Ostream& os) const
{
    calculatedPointPatchVectorField::write(os);

    if (tractionSeries_.size())
    {
        os.writeKeyword("tractionSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        tractionSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "traction", traction_);
#else
        traction_.writeEntry("traction", os);
#endif
    }

    if (pressureSeries_.size())
    {
        os.writeKeyword("pressureSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        pressureSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "pressure", pressure_);
#else
        pressure_.writeEntry("pressure", os);
#endif
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    solidTractionPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
