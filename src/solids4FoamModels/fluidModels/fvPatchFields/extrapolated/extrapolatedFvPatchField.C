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

#include "extrapolatedFvPatchField.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "emptyPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "scalarMatrices.H"
#include "volFields.H"
#include "skewCorrectionVectors.H"

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
    iPoints_(p.size()),
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
    iPoints_(p.size()),
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
    iPoints_(p.size()),
    zeroGradient_(true)
{
    if (dict.found("zeroGradient"))
    {
        zeroGradient_ = Switch(dict.lookup("zeroGradient"));
    }

    if (dict.found("value"))
    {
        Field<Type>::operator=(Field<Type>("value", dict, p.size()));
    }
    else
    {
        fixedGradientFvPatchField<Type>::evaluate();
    }
}


template<class Type>
extrapolatedFvPatchField<Type>::extrapolatedFvPatchField
(
    const extrapolatedFvPatchField<Type>& ptf
)
:
    fixedGradientFvPatchField<Type>(ptf),
    iPoints_(ptf.size()),
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
    iPoints_(ptf.size()),
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
//      gradient_.autoMap(m);
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
void extrapolatedFvPatchField<Type>::updateCoeffs()
{

    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const cellList& cells = mesh.cells();

#ifdef OPENFOAMESI
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const labelList& patchCells = this->patch().faceCells();
#else
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    const unallocLabelList& patchCells = this->patch().faceCells();
#endif

    const surfaceScalarField& weights = mesh.weights();
    const vectorField& faceCentres = mesh.faceCentres();
    const vectorField& cellCentres = mesh.cellCentres();

    skewCorrectionVectors scv(mesh);

    vectorField n = this->patch().nf();

    iPoints_.clear();
    iPoints_.setSize(this->patch().size());

    forAll(patchCells, faceI)
    {
        label curCell = patchCells[faceI];

        const labelList& curCellFaces = cells[curCell];

        iPoints_.set(faceI, new DynamicList<vector>());

        // Add face centre points
        forAll(curCellFaces, fI)
        {
            label curFace = curCellFaces[fI];

            label patchID = mesh.boundaryMesh().whichPatch(curFace);

            if(mesh.isInternalFace(curFace))
            {
                vector curFaceIntersection =
                    weights[curFace]
                   *(
                       cellCentres[owner[curFace]]
                     - cellCentres[neighbour[curFace]]
                    )
                  + cellCentres[neighbour[curFace]];

                iPoints_[faceI].append(curFaceIntersection);
            }
            else if (patchID != this->patch().index())
            {
                if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == cyclicPolyPatch::typeName
                )
                {
                    label start = mesh.boundaryMesh()[patchID].start();
                    label localFaceID = curFace - start;

#ifdef OPENFOAMESI
                    const labelList& cycPatchCells =
                        mesh.boundaryMesh()[patchID].faceCells();
#else
                    const unallocLabelList& cycPatchCells =
                        mesh.boundaryMesh()[patchID].faceCells();
#endif
                    label sizeby2 = cycPatchCells.size()/2;

                    if (localFaceID < sizeby2)
                    {
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

                        iPoints_[faceI].append(curFaceIntersection);
                    }
                    else
                    {
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

                        iPoints_[faceI].append(curFaceIntersection);
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

#ifdef OPENFOAMESI
                    const labelList& procPatchCells =
                        mesh.boundaryMesh()[patchID].faceCells();
#else
                    const unallocLabelList& procPatchCells =
                        mesh.boundaryMesh()[patchID].faceCells();
#endif

                    iPoints_[faceI].append
                    (
                        weights.boundaryField()[patchID][localFaceID]
                       *(
                            cellCentres[procPatchCells[localFaceID]]
                          - mesh.C().boundaryField()[patchID][localFaceID]
                        )
                      + mesh.C().boundaryField()[patchID][localFaceID]
                    );
                }
                else if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == emptyPolyPatch::typeName
                )
                {
                    iPoints_[faceI].append(faceCentres[curFace]);
                }
                else
                {
                    // Normal patches

                    iPoints_[faceI].append(faceCentres[curFace]);
                }
            }
        }
    }

    fixedGradientFvPatchField<Type>::updateCoeffs();
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

#ifdef OPENFOAMESI
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const labelList& patchCells = this->patch().faceCells();
#else
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    const unallocLabelList& patchCells = this->patch().faceCells();
#endif

    const surfaceScalarField& weights = mesh.weights();
//     const vectorField& faceCentres = mesh.faceCentres();
    const vectorField& cellCentres = mesh.cellCentres();

#ifdef OPENFOAMESIORFOUNDATION
    const Field<Type>& phiI = this->primitiveField();

    word fieldName =
        this->internalField().name();
#else
    const Field<Type>& phiI = this->internalField();

    word fieldName =
        this->dimensionedInternalField().name();
#endif

    const GeometricField<Type, fvPatchField, volMesh>& phi =
        mesh.lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            fieldName
        );

    vectorField n = this->patch().nf();
    const vectorField& C = this->patch().Cf();

    Field<Type> patchPhi(this->patch().size(), pTraits<Type>::zero);

    forAll(patchCells, faceI)
    {
        label curCell = patchCells[faceI];

        const labelList& curCellFaces = cells[curCell];

        DynamicList<Type> iPhi;
//         DynamicList<vector> iPoint;

//         // Add first cell centre point
//         iPhi.append(phiI[curCell]);
//         iPoint.append(cellCentres[curCell]);

        // Add face centre points
        forAll(curCellFaces, fI)
        {
            label curFace = curCellFaces[fI];

            label patchID = mesh.boundaryMesh().whichPatch(curFace);

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

//                 vector curFaceIntersection =
//                     weights[curFace]
//                    *(
//                        cellCentres[owner[curFace]]
//                      - cellCentres[neighbour[curFace]]
//                     )
//                   + cellCentres[neighbour[curFace]];

//                 iPoint.append(curFaceIntersection);
            }
            else if (patchID != this->patch().index())
            {
                if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == cyclicPolyPatch::typeName
                )
                {
                    label start = mesh.boundaryMesh()[patchID].start();
                    label localFaceID = curFace - start;

#ifdef OPENFOAMESI
                    const labelList& cycPatchCells =
                        mesh.boundaryMesh()[patchID].faceCells();
#else
                    const unallocLabelList& cycPatchCells =
                        mesh.boundaryMesh()[patchID].faceCells();
#endif

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

//                         vector curFaceIntersection =
//                             weights[curFace]
//                            *(
//                                 cellCentres[cycPatchCells[localFaceID]]
//                               - cellCentres
//                                 [
//                                     cycPatchCells[localFaceID + sizeby2]
//                                 ]
//                             )
//                           + cellCentres
//                             [
//                                 cycPatchCells[localFaceID + sizeby2]
//                             ];

//                         iPoint.append(curFaceIntersection);
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

//                         vector curFaceIntersection =
//                             weights[curFace]
//                            *(
//                                 cellCentres[cycPatchCells[localFaceID]]
//                               - cellCentres
//                                 [
//                                     cycPatchCells[localFaceID - sizeby2]
//                                 ]
//                             )
//                           + cellCentres
//                             [
//                                 cycPatchCells[localFaceID - sizeby2]
//                             ];

//                         iPoint.append(curFaceIntersection);
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

#ifdef OPENFOAMESI
                    const labelList& procPatchCells =
                        mesh.boundaryMesh()[patchID].faceCells();
#else
                    const unallocLabelList& procPatchCells =
                        mesh.boundaryMesh()[patchID].faceCells();
#endif

                    iPhi.append
                    (
                        weights.boundaryField()[patchID][localFaceID]
                       *(
                            phiI[procPatchCells[localFaceID]]
                          - phi.boundaryField()[patchID][localFaceID]
                        )
                      + phi.boundaryField()[patchID][localFaceID]
                    );


//                     const processorPolyPatch& procPatch =
//                         refCast<const processorPolyPatch>
//                         (
//                             mesh.boundaryMesh()[patchID]
//                         );

//                     vector ngbCellCentre =
//                         procPatch.neighbFaceCellCentres()[localFaceID];

//                     vector curFaceIntersection =
//                         faceCentres[curFace];

//                     iPoint.append(curFaceIntersection);
                }
                else if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == emptyPolyPatch::typeName
                )
                {
                    iPhi.append(phiI[curCell]);
//                     iPoint.append(faceCentres[curFace]);
                }
                else
                {
                    // Normal patches

                    label start = mesh.boundaryMesh()[patchID].start();
                    label localFaceID = curFace - start;

                    iPhi.append
                    (
                        phi.boundaryField()[patchID][localFaceID]
                    );

//                     iPoint.append(faceCentres[curFace]);
                }
            }
        }

        Type avgPhi = phiI[curCell];
        vector avgPoint = cellCentres[curCell];

//         Type avgPhi = average(iPhi);
//         vector avgPoint = average(iPoint);

        // Weights
        scalarField W(iPoints_[faceI].size(), 1.0);

        label nCoeffs = 3;
        scalarRectangularMatrix M
        (
            iPoints_[faceI].size(),
            nCoeffs,
            0.0
        );

//         vector origin = C[faceI];

        scalar L = max(mag(iPoints_[faceI]-avgPoint));

        for (label i=0; i<iPoints_[faceI].size(); i++)
        {
            scalar X = (iPoints_[faceI][i].x() - avgPoint.x())/L;
            scalar Y = (iPoints_[faceI][i].y() - avgPoint.y())/L;
            scalar Z = (iPoints_[faceI][i].z() - avgPoint.z())/L;

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
//         Info << M << endl;
//         Info << lsM << endl;
//         Info << det(lsM) << endl;
        tensor invLsM = inv(lsM);

        scalarRectangularMatrix curInvMatrix
        (
            nCoeffs,
            iPoints_[faceI].size(),
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
        Field<Type> source(iPoints_[faceI].size(), pTraits<Type>::zero);

        for (label i=0; i<iPoints_[faceI].size(); i++)
        {
//             source[i] = iPhi[i];
            source[i] = iPhi[i] - avgPhi;
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

//         vector dr = C[faceI] - origin;
        vector dr = (C[faceI] - avgPoint)/L;

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

#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
#else
    this->writeEntry("value", os);
#endif

    os.writeKeyword("zeroGradient")
        << zeroGradient_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
