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

#include "freeEdgeMomentFaPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "freeEdgeStencil.H"
#include "kirchhoffPlateSolid.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeEdgeMomentFaPatchScalarField::
freeEdgeMomentFaPatchScalarField
(
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedGradientFaPatchField<scalar>(p, iF),
    relaxFac_(1.0)
{}


Foam::freeEdgeMomentFaPatchScalarField::
freeEdgeMomentFaPatchScalarField
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
        faPatchField<scalar>::operator==
        (
            Field<scalar>("value", dict, p.size())
        );
    }
    else
    {
        faPatchField<scalar>::operator==(0.0);
    }

    gradient() = 0.0;
}


Foam::freeEdgeMomentFaPatchScalarField::
freeEdgeMomentFaPatchScalarField
(
    const freeEdgeMomentFaPatchScalarField& ptf,
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedGradientFaPatchField<scalar>(ptf, p, iF, mapper),
    relaxFac_(ptf.relaxFac_)
{}


#if (OPENFOAM < 2606)
Foam::freeEdgeMomentFaPatchScalarField::
freeEdgeMomentFaPatchScalarField
(
    const freeEdgeMomentFaPatchScalarField& ptf
)
:
    fixedGradientFaPatchField<scalar>(ptf),
    relaxFac_(ptf.relaxFac_)
{}
#endif


Foam::freeEdgeMomentFaPatchScalarField::
freeEdgeMomentFaPatchScalarField
(
    const freeEdgeMomentFaPatchScalarField& ptf,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedGradientFaPatchField<scalar>(ptf, iF),
    relaxFac_(ptf.relaxFac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeEdgeMomentFaPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Lookup the transverse displacement and its normal slope on this patch
    const areaScalarField& w = db().lookupObject<areaScalarField>("w");
    const scalarField dwdn(w.boundaryField()[patch().index()].snGrad());

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

    // Edge connectivity along the patch
    const labelListList nbrs(freeEdgeStencil::edgeNeighbours(patch()));
    const labelListList ends(freeEdgeStencil::edgeEndPoints(patch()));
    const vectorField t(freeEdgeStencil::edgeTangents(patch(), nbrs));
    const vectorField n(patch().edgeNormals());
    const vectorField& Ce = patch().edgeCentres();
    const pointField& points = patch().boundaryMesh().mesh().points();

    // Boundary values of the transverse displacement
    const scalarField& wb = w.boundaryField()[patch().index()];

    // Curvature of the edge times the tangential slope, kappa*dw/ds, where
    // kappa = t & dn/ds. This vanishes on a straight edge
    scalarField kappaDwds(dwdn.size(), 0.0);

    forAll(kappaDwds, edgeI)
    {
        const labelList& nbr = nbrs[edgeI];

        if (nbr.size() == 0)
        {
            continue;
        }

        // Edges behind and ahead along the tangent; one-sided at an end
        label m = edgeI;
        label p = edgeI;

        forAll(nbr, i)
        {
            if (((Ce[nbr[i]] - Ce[edgeI]) & t[edgeI]) > 0)
            {
                p = nbr[i];
            }
            else
            {
                m = nbr[i];
            }
        }

        const scalar h = mag(Ce[p] - Ce[m]);

        kappaDwds[edgeI] = (t[edgeI] & (n[p] - n[m]))*(wb[p] - wb[m])/sqr(h);
    }

    // Second derivative of the normal slope along the edge and first
    // derivative of the curvature term
    scalarField d2dwdn(dwdn.size(), 0.0);
    scalarField dKappaDwds(dwdn.size(), 0.0);

    forAll(d2dwdn, edgeI)
    {
        const labelList& nbr = nbrs[edgeI];

        // Distances, slopes and curvature terms behind (m) and ahead (p) of
        // the edge centre along the tangent
        scalar hm = 0.0;
        scalar hp = 0.0;
        scalar gm = 0.0;
        scalar gp = 0.0;
        scalar qm = 0.0;
        scalar qp = 0.0;

        if (nbr.size() == 2)
        {
            label m = nbr[0];
            label p = nbr[1];

            if (((Ce[p] - Ce[m]) & t[edgeI]) < 0)
            {
                Swap(m, p);
            }

            hm = mag(Ce[m] - Ce[edgeI]);
            gm = dwdn[m];
            qm = kappaDwds[m];
            hp = mag(Ce[p] - Ce[edgeI]);
            gp = dwdn[p];
            qp = kappaDwds[p];
        }
        else if (nbr.size() == 1 && ends[edgeI].size() == 1)
        {
            const label nbrI = nbr[0];
            const label pointI = ends[edgeI][0];
            const scalar h0 = mag(points[pointI] - Ce[edgeI]);

            // Virtual node at the corner
            scalar hc = 0.0;
            scalar gc = 0.0;
            scalar qc = 0.0;

            if (freeEdgeStencil::fixesValueAtPoint(w, patch(), pointI))
            {
                // Corner with a clamped or supported edge: the slope along
                // that edge, which is the normal slope here, vanishes, as
                // does the tangential slope
                hc = h0;
            }
            else
            {
                // Corner with another free edge: the twisting moment
                // vanishes so the slope is symmetric about the corner
                hc = 2.0*h0;
                gc = dwdn[edgeI];
                qc = kappaDwds[edgeI];
            }

            if (((Ce[nbrI] - Ce[edgeI]) & t[edgeI]) > 0)
            {
                hm = hc;
                gm = gc;
                qm = qc;
                hp = mag(Ce[nbrI] - Ce[edgeI]);
                gp = dwdn[nbrI];
                qp = kappaDwds[nbrI];
            }
            else
            {
                hm = mag(Ce[nbrI] - Ce[edgeI]);
                gm = dwdn[nbrI];
                qm = kappaDwds[nbrI];
                hp = hc;
                gp = gc;
                qp = qc;
            }
        }
        else
        {
            continue;
        }

        const scalar denom = hm*hp*(hm + hp);

        d2dwdn[edgeI] = 2.0*(hm*gp - (hm + hp)*dwdn[edgeI] + hp*gm)/denom;

        dKappaDwds[edgeI] =
            (
                sqr(hm)*qp - sqr(hp)*qm + (sqr(hp) - sqr(hm))*kappaDwds[edgeI]
            )/denom;
    }

    // Update the gradient
    gradient() =
        relaxFac_*D*(1.0 - nu)*(d2dwdn - dKappaDwds)
      + (1.0 - relaxFac_)*gradient();

    fixedGradientFaPatchField<scalar>::updateCoeffs();
}


void Foam::freeEdgeMomentFaPatchScalarField::write
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
        freeEdgeMomentFaPatchScalarField
    );
}

#endif // OPENFOAM_ORG

// ************************************************************************* //
