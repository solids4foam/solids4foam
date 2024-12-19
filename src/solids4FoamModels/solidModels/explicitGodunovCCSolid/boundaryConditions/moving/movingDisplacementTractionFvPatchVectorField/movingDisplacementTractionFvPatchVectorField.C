/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "movingDisplacementTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingDisplacementTractionFvPatchVectorField::
movingDisplacementTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    movingTractionFvPatchVectorField(p, iF),
    rho_(0.0),
    uMax_(vector::zero),
    tEnd_(0.0)
{}


movingDisplacementTractionFvPatchVectorField::
movingDisplacementTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    movingTractionFvPatchVectorField(p, iF),
    rho_(dict.lookup("density")),
    uMax_(dict.lookup("displacement")),
    tEnd_(readScalar(dict.lookup("endTime")))
{}


movingDisplacementTractionFvPatchVectorField::
movingDisplacementTractionFvPatchVectorField
(
    const movingDisplacementTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    movingTractionFvPatchVectorField(ptf, p, iF, mapper),
    rho_(ptf.rho_),
    uMax_(ptf.uMax_),
    tEnd_(ptf.tEnd_)
{}


movingDisplacementTractionFvPatchVectorField::
movingDisplacementTractionFvPatchVectorField
(
    const movingDisplacementTractionFvPatchVectorField& rifvpvf
)
:
    movingTractionFvPatchVectorField(rifvpvf),
    rho_(rifvpvf.rho_),
    uMax_(rifvpvf.uMax_),
    tEnd_(rifvpvf.tEnd_)
{}


movingDisplacementTractionFvPatchVectorField::
movingDisplacementTractionFvPatchVectorField
(
    const movingDisplacementTractionFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    movingTractionFvPatchVectorField(rifvpvf, iF),
    rho_(rifvpvf.rho_),
    uMax_(rifvpvf.uMax_),
    tEnd_(rifvpvf.tEnd_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingDisplacementTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    movingTractionFvPatchVectorField::autoMap(m);
}


void movingDisplacementTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    movingTractionFvPatchVectorField::rmap(ptf, addr);
}


void movingDisplacementTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar& t = this->db().time().value();
    const scalar A = 3*10;
    const scalar B = 4*15/tEnd_;
    const scalar C = 5*6/pow(tEnd_,2);

    linearMomentum() =
        rho_.value()*(uMax_/pow(tEnd_,3))*(A - B*t + C*t*t)*pow(t,2);

    movingTractionFvPatchVectorField::updateCoeffs();
}


void movingDisplacementTractionFvPatchVectorField::write(Ostream& os) const
{
    movingTractionFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    movingDisplacementTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //