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

#include "extrapolatedFvPatchField.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "emptyPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "scalarMatrices.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
extrapolatedFvPatchField<Type>::extrapolatedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    zeroGradient_(true)
{}


template<class Type>
extrapolatedFvPatchField<Type>::extrapolatedFvPatchField
(
    const extrapolatedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper),
    zeroGradient_(ptf.zeroGradient_)
{}


template<class Type>
extrapolatedFvPatchField<Type>::extrapolatedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    zeroGradient_(true)
{
    if (dict.found("zeroGradient"))
    {
        zeroGradient_ = Switch(dict.lookup("zeroGradient"));
    }

    fixedGradientFvPatchField<Type>::evaluate();
}


template<class Type>
extrapolatedFvPatchField<Type>::extrapolatedFvPatchField
(
    const extrapolatedFvPatchField<Type>& ptf
)
:
    fixedGradientFvPatchField<Type>(ptf),
    zeroGradient_(ptf.zeroGradient_)
{}


template<class Type>
extrapolatedFvPatchField<Type>::extrapolatedFvPatchField
(
    const extrapolatedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(ptf, iF),
    zeroGradient_(ptf.zeroGradient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void extrapolatedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchField<Type>::autoMap(m);
//     gradient_.autoMap(m);
}


template<class Type>
void extrapolatedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchField<Type>::rmap(ptf, addr);

//     const extrapolatedFvPatchField<Type>& fgptf =
//         refCast<const extrapolatedFvPatchField<Type> >(ptf);

//     gradient_.rmap(fgptf.gradient_, addr);
}


template<class Type>
void extrapolatedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const cellList& cells = mesh.cells();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const unallocLabelList& patchCells = this->patch().faceCells();

    const surfaceScalarField& weights = mesh.weights();
    const vectorField& faceCentres = mesh.faceCentres();
    const vectorField& cellCentres = mesh.cellCentres();

    const Field<Type>& phiI = this->internalField();

    word fieldName =
        this->dimensionedInternalField().name();

    const GeometricField<Type, fvPatchField, volMesh>& phi =
        mesh.lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            fieldName
        );

    vectorField n = this->patch().nf();
    vectorField C = this->patch().Cf();

    Field<Type> patchPhi(this->patch().size(), pTraits<Type>::zero);

    forAll(patchCells, faceI)
    {
        label curCell = patchCells[faceI];

        const labelList& curCellFaces = cells[curCell];

        DynamicList<Type> iPhi;
        DynamicList<vector> iPoint;

        // Add first cell centre point
        iPhi.append(phiI[curCell]);
        iPoint.append(cellCentres[curCell]);

        // Add face centre points
        forAll(curCellFaces, fI)
        {
            label curFace = curCellFaces[fI];

            if(mesh.isInternalFace(curFace))
            {
                iPhi.append
                (
                    weights[curFace]
                   *(
                       phiI[owner[curFace]]
                     - phiI[neighbour[curFace]]
                    )
                  + phiI[neighbour[curFace]]
                );

                vector curFaceIntersection =
                    weights[curFace]
                   *(
                       cellCentres[owner[curFace]]
                     - cellCentres[neighbour[curFace]]
                    )
                  + cellCentres[neighbour[curFace]];

                iPoint.append(curFaceIntersection);
            }
            else
            {
                label patchID = mesh.boundaryMesh().whichPatch(curFace);
                if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == cyclicPolyPatch::typeName
                )
                {
                    label start = mesh.boundaryMesh()[patchID].start();
                    label localFaceID = curFace - start;

                    const unallocLabelList& cycPatchCells =
                        mesh.boundaryMesh()[patchID].faceCells();

                    label sizeby2 = cycPatchCells.size()/2;

                    if (localFaceID < sizeby2)
                    {
                        iPhi.append
                        (
                            weights.boundaryField()[patchID][localFaceID]
                           *(
                               phiI[cycPatchCells[localFaceID]]
                             - phiI[cycPatchCells[localFaceID + sizeby2]]
                            )
                          + phiI[cycPatchCells[localFaceID + sizeby2]]
                        );

                        vector curFaceIntersection =
                            weights[curFace]
                           *(
                                cellCentres[cycPatchCells[localFaceID]]
                              - cellCentres
                                [
                                    cycPatchCells[localFaceID + sizeby2]
                                ]
                            )
                          + cellCentres
                            [
                                cycPatchCells[localFaceID + sizeby2]
                            ];

                        iPoint.append(curFaceIntersection);
                    }
                    else
                    {
                        iPhi.append
                        (
                            weights.boundaryField()[patchID][localFaceID]
                           *(
                               phiI[cycPatchCells[localFaceID]]
                             - phiI[cycPatchCells[localFaceID - sizeby2]]
                            )
                          + phiI[cycPatchCells[localFaceID - sizeby2]]
                        );

                        vector curFaceIntersection =
                            weights[curFace]
                           *(
                                cellCentres[cycPatchCells[localFaceID]]
                              - cellCentres
                                [
                                    cycPatchCells[localFaceID - sizeby2]
                                ]
                            )
                          + cellCentres
                            [
                                cycPatchCells[localFaceID - sizeby2]
                            ];

                        iPoint.append(curFaceIntersection);
                    }
                }
                else if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == processorPolyPatch::typeName
                )
                {
                    label start = mesh.boundaryMesh()[patchID].start();
                    label localFaceID = curFace - start;

                    const unallocLabelList& procPatchCells =
                        mesh.boundaryMesh()[patchID].faceCells();

                    iPhi.append
                    (
                        weights.boundaryField()[patchID][localFaceID]
                       *(
                            phiI[procPatchCells[localFaceID]]
                          - phi.boundaryField()[patchID][localFaceID]
                        )
                      + phi.boundaryField()[patchID][localFaceID]
                    );

                    vector curFaceIntersection =
                        weights[curFace]
                       *(
                            cellCentres[procPatchCells[localFaceID]]
                          - mesh.C().boundaryField()[patchID][localFaceID]
                        )
                      + mesh.C().boundaryField()[patchID][localFaceID];

                    iPoint.append(curFaceIntersection);
                }
                else if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == emptyPolyPatch::typeName
                )
                {
                    iPhi.append(phiI[curCell]);
                    iPoint.append(faceCentres[curFace]);
                }
            }
        }

        Type avgPhi = average(iPhi);
        vector avgPoint = average(iPoint);

        // Weights
        scalarField W(iPoint.size(), 1.0);

        label nCoeffs = 3;
        scalarRectangularMatrix M
        (
            iPoint.size(),
            nCoeffs,
            0.0
        );

        for (label i=0; i<iPoint.size(); i++)
        {
            scalar X = iPoint[i].x() - avgPoint.x();
            scalar Y = iPoint[i].y() - avgPoint.y();
            scalar Z = iPoint[i].z() - avgPoint.z();

            M[i][0] = X;
            M[i][1] = Y;
            M[i][2] = Z;
        }

        // Applying weights
        for (label i=0; i<M.n(); i++)
        {
            for (label j=0; j<M.m(); j++)
            {
                M[i][j] *= W[i];
            }
        }

        tensor lsM = tensor::zero;

        for (label i=0; i<3; i++)
        {
            for (label j=0; j<3; j++)
            {
                for (label k=0; k<M.n(); k++)
                {
                    lsM(i,j) += M[k][i]*M[k][j];
                }
            }
        }

        // Calculate inverse
        tensor invLsM = inv(lsM);

        scalarRectangularMatrix curInvMatrix
        (
            nCoeffs,
            iPoint.size(),
            0.0
        );

        for (label i=0; i<3; i++)
        {
            for (label j=0; j<M.n(); j++)
            {
                for (label k=0; k<3; k++)
                {
                    curInvMatrix[i][j] += invLsM(i,k)*M[j][k]*W[j];
                }
            }
        }

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source(iPoint.size(), pTraits<Type>::zero);

        for (label i=0; i<iPoint.size(); i++)
        {
            source[i] = iPhi[i] - avgPhi;
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector dr = C[faceI] - avgPoint;

        patchPhi[faceI] =
            avgPhi
          + coeffs[0]*dr.x()
          + coeffs[1]*dr.y()
          + coeffs[2]*dr.z();
    }

    Field<Type>::operator=(patchPhi);

    if (zeroGradient_)
    {
        this->gradient() = pTraits<Type>::zero;
    }
    else
    {
        this->gradient() =
            (patchPhi - this->patchInternalField())
           *this->patch().deltaCoeffs();
    }

    fvPatchField<Type>::evaluate();
}


template<class Type>
void extrapolatedFvPatchField<Type>::write(Ostream& os) const
{
    fixedGradientFvPatchField<Type>::write(os);

    os.writeKeyword("zeroGradient")
        << zeroGradient_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
