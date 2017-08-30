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

#include "wedgeFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
tmp<scalarField> wedgeFvPatchField<scalar>::snGrad() const
{
    return tmp<scalarField >(new scalarField(size(), 0.0));
}


template<>
void wedgeFvPatchField<scalar>::evaluate(const Pstream::commsTypes)
{
    if (!updated())
    {
        updateCoeffs();
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // word objName("grad(" + this->dimensionedInternalField().name() + ")");
    
    if 
    (
        mesh.objectRegistry::foundObject<volVectorField>
        (
            "grad(" + this->dimensionedInternalField().name() + ")"
        )
     && true
    )
    {
        // Philip Cardiff, UCD
        // Zeljko Tukovic, FSB Zagreb

//         Info << "wedgeFvPatchField<scalar>::evaluate(): "
//             << "Using stored gradient" << endl;

        const wedgeFvPatch& wedgePatch =
            refCast<const wedgeFvPatch>(this->patch());

        // Rotate patchC field back to centre plane to find 
        // transformed cell centres
        const vectorField& patchC = patch().patch().faceCentres();
        vectorField transC = wedgePatch.faceT().T() & patchC;

        // Calculate correction vector which connects actual cell centre 
        // to the transformed cell centre
        const vectorField k = transC - patch().Cn();

        const fvPatchField<vector>& gradField =
            patch().lookupPatchField<volVectorField, vector>
            (
                "grad(" + this->dimensionedInternalField().name() + ")"
            );

        Field<scalar> pif = this->patchInternalField();
        pif += (k & gradField.patchInternalField());

        fvPatchScalarField::operator==
        (
            transform(wedgePatch.faceT(), pif)
        );
    }
    else
    {
        fvPatchScalarField::operator==(patchInternalField());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
