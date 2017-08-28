/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "wedgeFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePatchFields(wedge);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<vectorField> wedgeFvPatchField<vector>::snGrad() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    Field<vector> pif = this->patchInternalField();

    if 
    (
     
        mesh.foundObject<volTensorField>
        (
            "grad(" + this->dimensionedInternalField().name() + ")"
        )
     && true
    )
    {
        // Philip Cardiff, UCD
        // Zeljko Tukovic, FSB Zagreb

//         Info << "wedgeFvPatchField<vector>::snGrad(): " 
//             << "Using stored gradient for field: " 
//             << this->dimensionedInternalField().name()<< endl;

        // Method:
        // We project the patch normal back to the wedge centrePlane and 
        // find the intersection point. We then extrapolate U from 
        // the Cn to this intersection point (projC) using gradU 
        // looked up from the solver.
        // We can then calculate snGrad by getting the difference 
        // between projU and transformed projU (i.e. across the wedge), 
        // and then dividing by the magnitude of the distance between them.

        const wedgePolyPatch& wedgePatch =
            refCast<const wedgePolyPatch>(patch().patch());

        const vectorField& patchC = patch().patch().faceCentres();
        vectorField nHat = this->patch().nf();
        const vector& centreN = wedgePatch.centreNormal();
        scalarField d = ((patch().Cn() - patchC) & centreN)/(nHat & centreN);
        vectorField projC = d*nHat + patchC;


        // Calculate correction vector which connects actual cell centre to the
        // transformed cell centre
        const vectorField k = projC - patch().Cn();

        const fvPatchField<tensor>& gradField =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + this->dimensionedInternalField().name() + ")"
            );

        Field<vector> projU =
            pif + (k & gradField.patchInternalField());

        // Calculate delta coeffs from proj position on centre plane 
        // to transformed projected position
        scalarField projDeltaCoeff =
            1.0/mag(transform(wedgePatch.cellT(), projC) - projC);

        return
        (
            transform(wedgePatch.cellT(), projU) - projU
        )*projDeltaCoeff;
    }

    return
    (
        transform(refCast<const wedgeFvPatch>(this->patch()).cellT(), pif) 
      - pif
    )*(0.5*this->patch().deltaCoeffs());
}


template<>
void wedgeFvPatchField<vector>::evaluate(const Pstream::commsTypes)
{
    if (!updated())
    {
        updateCoeffs();
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    
    if 
    (
        mesh.foundObject<volTensorField>
        (
            "grad(" + this->dimensionedInternalField().name() + ")"
        )
     && true
    )
    {
        // Philip Cardiff, UCD
        // Zeljko Tukovic, FSB Zagreb

//         Info << "wedgeFvPatchField<vector>::evaluate(): " 
//             << "Using stored gradient for field: "
//             << this->dimensionedInternalField().name()<< endl;

        const wedgeFvPatch& wedgePatch =
            refCast<const wedgeFvPatch>(this->patch());

        // Rotate patchC field back to centre plane to find 
        // transformed cell centres
        const vectorField& patchC = patch().patch().faceCentres();
        vectorField transC = wedgePatch.faceT().T() & patchC;

        // Calculate correction vector which connects actual cell centre 
        // to the transformed cell centre
        const vectorField k = transC - patch().Cn();

        const fvPatchField<tensor>& gradField =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + this->dimensionedInternalField().name() + ")"
            );

        Field<vector> pif = this->patchInternalField();
        pif += (k & gradField.patchInternalField());

        fvPatchVectorField::operator==
        (
            transform(wedgePatch.faceT(), pif)
        );
    }
    else
    {
        fvPatchField<vector>::operator==
        (
            transform
            (
                refCast<const wedgeFvPatch>(this->patch()).faceT(),
                this->patchInternalField()
            )
        );
    }
}



} // End namespace Foam

// ************************************************************************* //
