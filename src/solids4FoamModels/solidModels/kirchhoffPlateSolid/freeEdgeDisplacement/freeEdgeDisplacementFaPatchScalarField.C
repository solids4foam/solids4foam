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

#ifndef OPENFOAM_ORG

#include "freeEdgeDisplacementFaPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "freeEdgeStencil.H"
#include "kirchhoffPlateSolid.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedGradientFaPatchField<scalar>(p, iF),
    relaxFac_(1.0)
{}


Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFaPatchField<scalar>(p, iF),
    relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0))

{
   if (dict.found("value"))
   {
       faPatchField<scalar>::operator==(Field<scalar>("value", dict, p.size()));
   }
   else
   {
       faPatchField<scalar>::operator==(0.0);
   }

   gradient() = 0.0;
}


Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const freeEdgeDisplacementFaPatchScalarField& ptf,
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedGradientFaPatchField<scalar>(ptf, p, iF, mapper),
    relaxFac_(ptf.relaxFac_)
{}


#if (OPENFOAM < 2606)
Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const freeEdgeDisplacementFaPatchScalarField& ptf
)
:
    fixedGradientFaPatchField<scalar>(ptf),
    relaxFac_(ptf.relaxFac_)
{}
#endif


Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const freeEdgeDisplacementFaPatchScalarField& ptf,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedGradientFaPatchField<scalar>(ptf, iF),
    relaxFac_(ptf.relaxFac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeEdgeDisplacementFaPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Lookup angle of rotation field
    const faPatchField<vector>& theta =
        patch().lookupPatchField<areaVectorField, vector>("theta");

    // Lookup the gradient of rotation field
    const faPatchField<tensor>& gradTheta =
        patch().lookupPatchField<areaTensorField, tensor>("grad(theta)");

    // Calculate correction vectors
    const vectorField n(patch().edgeNormals());
    const vectorField delta(patch().delta());
    const vectorField k((I - sqr(n)) & delta);

    // Normal slope at the foot of the normal from the cell centre, including
    // the non-orthogonal correction, which matters on curved edges
    const scalarField nDotThetaPif
    (
        n & (theta.patchInternalField() + (k & gradTheta.patchInternalField()))
    );

    // Lookup moment sum field
    const faPatchField<scalar>& M =
        patch().lookupPatchField<areaScalarField, scalar>("M");

    // Lookup fvMesh
    const fvMesh& vmesh =
        refCast<const fvMesh>(patch().boundaryMesh().mesh().mesh());

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(vmesh);

    // Cast the solid model to a Kirchhoff plate solid model
    const solidModels::kirchhoffPlateSolid& plateSolid =
        refCast<const solidModels::kirchhoffPlateSolid>(solMod);

    // Lookup flexural stiffness from the solver
    const scalar D = plateSolid.bendingStiffness().value();
    const scalar nu = plateSolid.nu().value();

    // Normal slope at the edge extrapolated from the cell centre: the normal
    // curvature on a free edge is w_nn = nu/(1 - nu) M/D, where M is taken
    // from the adjacent cell
    const scalarField extrapolated
    (
       -nDotThetaPif
      + (nu/(1.0 - nu))*(M.patchInternalField()/D)/patch().deltaCoeffs()
    );
    scalarField newGradient(extrapolated);

    // At a corner shared with another free edge the twisting moment vanishes,
    // so the normal slope does not vary along the edge there: take the slope
    // of the neighbouring edge, since the extrapolation above leaves the
    // corner cell undetermined
    const labelListList nbrs(freeEdgeStencil::edgeNeighbours(patch()));
    const labelListList ends(freeEdgeStencil::edgeEndPoints(patch()));
    const areaScalarField& w = db().lookupObject<areaScalarField>("w");

    forAll(newGradient, edgeI)
    {
        if
        (
            nbrs[edgeI].size() == 1
         && ends[edgeI].size() == 1
         && !freeEdgeStencil::fixesValueAtPoint(w, patch(), ends[edgeI][0])
        )
        {
            newGradient[edgeI] = extrapolated[nbrs[edgeI][0]];
        }
    }

    gradient() = relaxFac_*newGradient + (1.0 - relaxFac_)*gradient();

    fixedGradientFaPatchField<scalar>::updateCoeffs();
}


void Foam::freeEdgeDisplacementFaPatchScalarField::write
(
    Ostream& os
) const
{
    faPatchField<scalar>::write(os);

    os.writeKeyword("relaxationFactor")
        << relaxFac_ << token::END_STATEMENT << endl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeFaPatchTypeField
    (
        faPatchScalarField,
        freeEdgeDisplacementFaPatchScalarField
    );
}

#endif // OPENFOAM_ORG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
