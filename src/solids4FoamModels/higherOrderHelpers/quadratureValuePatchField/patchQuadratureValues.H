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

InNamespace
    Foam

Description
    Values of a fixed-value patch field at its face quadrature points, as
    used by the high-order least-squares schemes for the boundary data of
    their reconstructions.

    A patch field that inherits quadratureValuePatchField evaluates the values
    itself. Any other patch field must have a uniform value, so that copying
    it to the quadrature points is exact.

Author
    Amirhossein Taran, UCD.

\*---------------------------------------------------------------------------*/

#ifndef patchQuadratureValues_H
#define patchQuadratureValues_H

#include "quadratureValuePatchField.H"
#include "fvMeshQuadrature.H"
#include "fvPatchField.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Return the values of patch field pf at the face quadrature points of quad
template<class Type>
autoPtr<CompactListList<Type>> patchQuadratureValues
(
    const fvPatchField<Type>& pf,
    const fvMeshQuadrature& quad
)
{
    const quadratureValuePatchField<Type>* qvpfPtr =
        dynamic_cast<const quadratureValuePatchField<Type>*>(&pf);

    if (qvpfPtr)
    {
        return qvpfPtr->quadratureValues(quad);
    }

    // Copying the face-centre value to the quadrature points of a non-uniform
    // patch would limit the boundary data to second-order accuracy
    forAll(pf, faceI)
    {
        if (pf[faceI] != pf[0])
        {
            FatalErrorInFunction
                << "The value of field " << internalFieldName(pf)
                << " on patch " << pf.patch().name() << " is not uniform, "
                << "but patch field type " << pf.type() << " does not "
                << "evaluate it at the face quadrature points." << nl
                << "Derive the patch field from quadratureValuePatchField to "
                << "use it with a high-order least-squares scheme"
                << abort(FatalError);
        }
    }

    // faceQuadPoints is a list for the whole mesh
    auto& faceQuadPoints = compactListListCRef(quad.faceQuadPoints());

    labelList nQpPerFace(pf.size(), 0);
    const label start = pf.patch().patch().start();
    forAll(nQpPerFace, faceI)
    {
        nQpPerFace[faceI] = faceQuadPoints[faceI + start].size();
    }

    autoPtr<CompactListList<Type>> tQuadPointsValue
    (
        new CompactListList<Type>(nQpPerFace)
    );

    // Get a reference to the actual data for easier access
    CompactListList<Type>& quadPointsValue = autoPtrRef(tQuadPointsValue);

    forAll(pf, faceI)
    {
        forAll(quadPointsValue[faceI], pointI)
        {
            quadPointsValue[faceI][pointI] = pf[faceI];
        }
    }

    return tQuadPointsValue;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
