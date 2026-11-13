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
    itself. A fixedValue patch field must have a uniform value, so that copying
    it to the quadrature points is exact. Any other patch field type that fixes
    the value (e.g. mixed and directionMixed) gives an error.

Author
    Amirhossein Taran, UCD.

\*---------------------------------------------------------------------------*/

#ifndef patchQuadratureValues_H
#define patchQuadratureValues_H

#include "quadratureValuePatchField.H"
#include "fvMeshQuadrature.H"
#include "fvPatchField.H"
#include "fixedValueFvPatchField.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Return true if patch field pf has the same value on all faces
//  Values that differ only by round-off (e.g. after mapping or a
//  transformation) are taken as uniform, using a tolerance of 1e-12 times the
//  largest value magnitude on the patch
template<class Type>
bool uniformPatchValue(const fvPatchField<Type>& pf)
{
    scalar maxMag = 0;
    forAll(pf, faceI)
    {
        maxMag = max(maxMag, mag(pf[faceI]));
    }

    const scalar tol = 1e-12*maxMag;

    forAll(pf, faceI)
    {
        if (mag(pf[faceI] - pf[0]) > tol)
        {
            return false;
        }
    }

    return true;
}


//- Return the face-centre values of patch field pf copied to the face
//  quadrature points of quad
//  This is exact only where the value is constant over each face
template<class Type>
autoPtr<CompactListList<Type>> uniformQuadratureValues
(
    const fvPatchField<Type>& pf,
    const fvMeshQuadrature& quad
)
{
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


//- Stop with an error if the face quadrature points of quad on the patch of
//  pf differ from those of ownQuad, the quadrature at which pf evaluates its
//  values
//  Two quadratures of the same mesh with the same settings (e.g. those of two
//  least-squares schemes) have the same points and are accepted
template<class Type>
void checkSameQuadraturePoints
(
    const fvPatchField<Type>& pf,
    const fvMeshQuadrature& quad,
    const fvMeshQuadrature& ownQuad
)
{
    if (&quad == &ownQuad)
    {
        return;
    }

    auto& faceQuadPoints = compactListListCRef(quad.faceQuadPoints());
    auto& ownFaceQuadPoints = compactListListCRef(ownQuad.faceQuadPoints());
    const label start = pf.patch().patch().start();
    const scalarField& magSf = pf.patch().magSf();

    bool samePoints = true;
    for (label faceI = 0; samePoints && faceI < pf.size(); ++faceI)
    {
        const label globalFaceID = faceI + start;

        samePoints =
        (
            faceQuadPoints[globalFaceID].size()
         == ownFaceQuadPoints[globalFaceID].size()
        );

        // The tolerance is relative to the face size
        const scalar tol = 1e-10*Foam::sqrt(magSf[faceI]);

        forAll(faceQuadPoints[globalFaceID], pointI)
        {
            if (!samePoints)
            {
                break;
            }

            samePoints =
            (
                mag
                (
                    faceQuadPoints[globalFaceID][pointI]
                  - ownFaceQuadPoints[globalFaceID][pointI]
                ) <= tol
            );
        }
    }

    if (!samePoints)
    {
        FatalErrorInFunction
            << "Patch field type " << pf.type() << " of field "
            << internalFieldName(pf) << " on patch " << pf.patch().name()
            << " evaluates its values at the face quadrature points of its "
            << "own least-squares scheme, but the values were requested at "
            << "different quadrature points" << exit(FatalError);
    }
}


//- Return the values of patch field pf at the face quadrature points of quad
template<class Type>
autoPtr<CompactListList<Type>> patchQuadratureValues
(
    const fvPatchField<Type>& pf,
    const fvMeshQuadrature& quad
)
{
    // faceQuadPoints is a list for the whole mesh
    auto& faceQuadPoints = compactListListCRef(quad.faceQuadPoints());
    const label start = pf.patch().patch().start();

    const quadratureValuePatchField<Type>* qvpfPtr =
        dynamic_cast<const quadratureValuePatchField<Type>*>(&pf);

    if (qvpfPtr)
    {
        autoPtr<CompactListList<Type>> tQuadPointsValue =
            qvpfPtr->quadratureValues(quad);

        // The values must match the quadrature points of quad on every face
        auto& quadPointsValue =
            compactListListCRef(autoPtrRef(tQuadPointsValue));

        bool sizesMatch = (quadPointsValue.size() == pf.size());
        for (label faceI = 0; sizesMatch && faceI < pf.size(); ++faceI)
        {
            sizesMatch =
            (
                quadPointsValue[faceI].size()
             == faceQuadPoints[faceI + start].size()
            );
        }

        if (!sizesMatch)
        {
            FatalErrorInFunction
                << "Patch field type " << pf.type() << " of field "
                << internalFieldName(pf) << " on patch " << pf.patch().name()
                << " returned a number of values per face that differs from "
                << "the number of face quadrature points" << exit(FatalError);
        }

        return tQuadPointsValue;
    }

    // Only a plain prescribed value can be copied to the quadrature points.
    // Mixed and directionMixed patch fields also fix the value, but their
    // value depends on the solution.
    if (!isA<fixedValueFvPatchField<Type>>(pf))
    {
        FatalErrorInFunction
            << "Patch field type " << pf.type() << " of field "
            << internalFieldName(pf) << " on patch " << pf.patch().name()
            << " fixes the value, but it is neither a fixedValue patch field "
            << "nor derived from quadratureValuePatchField." << nl
            << "The high-order least-squares schemes need the boundary "
            << "values at the face quadrature points, which this patch field "
            << "does not provide" << exit(FatalError);
    }

    // Copying the face-centre value to the quadrature points of a non-uniform
    // patch would limit the boundary data to second-order accuracy
    if (!uniformPatchValue(pf))
    {
        FatalErrorInFunction
            << "The value of field " << internalFieldName(pf)
            << " on patch " << pf.patch().name() << " is not uniform (to a "
            << "relative tolerance of 1e-12), but patch field type "
            << pf.type() << " does not evaluate it at the face quadrature "
            << "points." << nl
            << "Derive the patch field from quadratureValuePatchField to "
            << "use it with a high-order least-squares scheme"
            << exit(FatalError);
    }

    return uniformQuadratureValues(pf, quad);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
