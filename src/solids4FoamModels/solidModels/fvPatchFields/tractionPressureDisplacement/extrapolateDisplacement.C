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

#include "tractionPressureDisplacementFvPatchVectorField.H"
#include "volFields.H"
#include "cyclicPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionPressureDisplacementFvPatchVectorField::
makeInterpolationPoints()
{
    // Info << "makeInterpolationPoints()" << endl;

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const cellList& cells = mesh.cells();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const unallocLabelList& patchCells = this->patch().faceCells();

    const surfaceScalarField& weights = mesh.weights();
    const vectorField& faceCentres = mesh.faceCentres();
    const vectorField& cellCentres = mesh.cellCentres();

    vectorField n = this->patch().nf();

    word UName =
        this->dimensionedInternalField().name();

    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
            UName
        );

    iPointsIndices_.clear();
    iPointsIndices_.setSize(this->patch().size());

    iPoints_.clear();
    iPoints_.setSize(this->patch().size());

    iMatrices_.clear();
    iMatrices_.setSize(this->patch().size());

    forAll(patchCells, faceI)
    {
        label curCell = patchCells[faceI];

        const labelList& curCellFaces = cells[curCell];

        iPointsIndices_.set(faceI, new DynamicList<label>());
        iPoints_.set(faceI, new DynamicList<vector>());

        // Add first cell centre point
        iPointsIndices_[faceI].append(curCell);
        iPoints_[faceI].append(cellCentres[curCell]);

        // Add face centre points
        forAll(curCellFaces, fI)
        {
            label curFace = curCellFaces[fI];

            label patchID = mesh.boundaryMesh().whichPatch(curFace);

            if(mesh.isInternalFace(curFace))
            {
                if (curCell == owner[curFace])
                {
                    iPointsIndices_[faceI].append(neighbour[curFace]);
                    iPoints_[faceI].append(cellCentres[neighbour[curFace]]);
                }
                else
                {
                    iPointsIndices_[faceI].append(owner[curFace]);
                    iPoints_[faceI].append(cellCentres[owner[curFace]]);
                }
            }
            else if
            (
               !isA<tractionPressureDisplacementFvPatchVectorField>
                (
                    U.boundaryField()[patchID]
                )
            )
            {
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

                    const processorPolyPatch& ppp =
                        refCast<const processorPolyPatch>
                        (
                            mesh.boundaryMesh()[patchID]
                        );

                    const vector& neiCellCentre =
                        ppp.neighbFaceCellCentres()[localFaceID];

                    iPoints_[faceI].append(neiCellCentre);
                }
                else if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == wedgePolyPatch::typeName
                )
                {
                    const wedgePolyPatch& wpp =
                        refCast<const wedgePolyPatch>
                        (
                            mesh.boundaryMesh()[patchID]
                        );

                    vector curCellCentre =
                        transform
                        (
                            wpp.cellT(),
                            cellCentres[curCell]
                        );

                    iPointsIndices_[faceI].append(curCell);
                    iPoints_[faceI].append(curCellCentre);
                }
                else if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == emptyPolyPatch::typeName
                )
                {
                    iPointsIndices_[faceI].append(curCell);
                    iPoints_[faceI].append(faceCentres[curFace]);
                }
                else
                {
                    // // Normal patches
                    iPointsIndices_[faceI].append(curCell);
                    iPoints_[faceI].append(faceCentres[curFace]);
                }
            }
        }

        // Weights
        scalarField W(iPoints_[faceI].size(), 1.0);

        // Weighted averate
        vector avgPoint = average(iPoints_[faceI]);

        label nCoeffs = 3;
        scalarRectangularMatrix M
        (
            iPoints_[faceI].size(),
            nCoeffs,
            0.0
        );

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

        // Check determinant
        scalar detLsMat = det(lsM);

        if (mag(detLsMat)>SMALL)
        {
            // Calculate inverse
            tensor invLsM = hinv(lsM);

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

            // Save matrix
            iMatrices_.set(faceI, new scalarRectangularMatrix(curInvMatrix));
        }
        else
        {
            Info << this->patch().name() << ", "<< faceI
                 << ", detLsMat: " << detLsMat << endl;

            scalarRectangularMatrix curInvMatrix
            (
                nCoeffs,
                iPoints_[faceI].size(),
                0.0
            );

            // Save matrix
            iMatrices_.set(faceI, new scalarRectangularMatrix(curInvMatrix));
        }
    }
}


// XXX: never used?
tmp<vectorField> tractionPressureDisplacementFvPatchVectorField::
extrapolatedDerivative
(
    const volTensorField& gradU,
    const vectorField& dir
)
{
    tmp<vectorField> tExDeriv
    (
        new vectorField(this->size(), vector::zero)
    );
    vectorField& patchDirGradU = tExDeriv();

    if (iPoints_.size() == 0)
    {
        makeInterpolationPoints();
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const cellList& cells = mesh.cells();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const unallocLabelList& patchCells = this->patch().faceCells();

    const surfaceScalarField& weights = mesh.weights();

    const tensorField& gradUI = gradU.internalField();

    word UName =
        this->dimensionedInternalField().name();

    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
            UName
        );

    vectorField n = this->patch().nf();
    const vectorField& C = this->patch().Cf();

    forAll(patchCells, faceI)
    {
        label curCell = patchCells[faceI];

        const labelList& curCellFaces = cells[curCell];
        DynamicList<vector> iDirGradU;

        // Add first cell centre point
        iDirGradU.append(dir[faceI] & gradUI[curCell]);

        // Add face centre points
        forAll(curCellFaces, fI)
        {
            label curFace = curCellFaces[fI];

            label patchID = mesh.boundaryMesh().whichPatch(curFace);

            if(mesh.isInternalFace(curFace))
            {
                if (curCell == owner[curFace])
                {
                    iDirGradU.append(dir[faceI] & gradUI[neighbour[curFace]]);
                }
                else
                {
                    iDirGradU.append(dir[faceI] & gradUI[owner[curFace]]);
                }
            }
            else if
            (
               !isA<tractionPressureDisplacementFvPatchVectorField>
                (
                    U.boundaryField()[patchID]
                )
            )
            {
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
                        iDirGradU.append
                        (
                            weights.boundaryField()[patchID][localFaceID]
                           *(
                                (
                                    dir[faceI]
                                  & gradUI[cycPatchCells[localFaceID]]
                                )
                              - (
                                    dir[faceI]
                                  & gradUI[cycPatchCells[localFaceID + sizeby2]]
                                )
                            )
                          + (
                                dir[faceI]
                              & gradUI[cycPatchCells[localFaceID + sizeby2]]
                            )
                        );
                    }
                    else
                    {
                        iDirGradU.append
                        (
                            weights.boundaryField()[patchID][localFaceID]
                           *(
                                (
                                    dir[faceI]
                                  & gradUI[cycPatchCells[localFaceID]]
                                )
                              - (
                                    dir[faceI]
                                  & gradUI[cycPatchCells[localFaceID - sizeby2]]
                                )
                            )
                          + (
                                dir[faceI]
                              & gradUI[cycPatchCells[localFaceID - sizeby2]]
                            )
                        );
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

                    iDirGradU.append
                    (
                        weights.boundaryField()[patchID][localFaceID]
                       *(
                            (
                                dir[faceI]
                              & gradUI[procPatchCells[localFaceID]]
                            )
                          - (
                                dir[faceI]
                              & gradU.boundaryField()[patchID][localFaceID]
                            )
                        )
                      + (
                            dir[faceI]
                          & gradU.boundaryField()[patchID][localFaceID]
                        )
                    );
                }
                else if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == emptyPolyPatch::typeName
                )
                {
                    iDirGradU.append(dir[faceI] & gradUI[curCell]);
                }
                else
                {
                    // Normal patches
                }
            }
        }

        // Weights
        scalarField W(iPoints_[faceI].size(), 1.0);

        // Weighted averate
        vector avgDirGradU = average(iDirGradU);
        vector avgPoint = average(iPoints_[faceI]);

        label nCoeffs = 3;

        scalar L = max(mag(iPoints_[faceI]-avgPoint));

        Field<vector> coeffs(nCoeffs, vector::zero);
        Field<vector> source(iPoints_[faceI].size(), vector::zero);

        for (label i=0; i<iPoints_[faceI].size(); i++)
        {
            source[i] = iDirGradU[i] - avgDirGradU;
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += iMatrices_[faceI][i][j]*source[j];
                // coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector dr = (C[faceI] - avgPoint)/L;

        patchDirGradU[faceI] =
            avgDirGradU
          + coeffs[0]*dr.x()
          + coeffs[1]*dr.y()
          + coeffs[2]*dr.z();
    }

    if (iPoints_.size())
    {
        iPoints_.clear();
        iPointsIndices_.clear();
        iMatrices_.clear();
    }

    return tExDeriv;
}


tmp<tensorField> tractionPressureDisplacementFvPatchVectorField::
extrapolatedGradient
(
    const volTensorField& gradU
)
{
    tmp<tensorField> tExDeriv
    (
        new tensorField(this->size(), tensor::zero)
    );
    tensorField& patchDirGradU = tExDeriv();

    if (iPoints_.size() == 0)
    {
        makeInterpolationPoints();
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const cellList& cells = mesh.cells();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const unallocLabelList& patchCells = this->patch().faceCells();

    const surfaceScalarField& weights = mesh.weights();

    const tensorField& gradUI = gradU.internalField();

    word UName =
        this->dimensionedInternalField().name();

    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
            UName
        );

    vectorField n = this->patch().nf();
    const vectorField& C = this->patch().Cf();

    forAll(patchCells, faceI)
    {
        label curCell = patchCells[faceI];

        const labelList& curCellFaces = cells[curCell];
        DynamicList<tensor> iDirGradU;

        // // Add first cell centre point
        iDirGradU.append(gradUI[curCell]);

        // Add face centre points
        forAll(curCellFaces, fI)
        {
            label curFace = curCellFaces[fI];

            label patchID = mesh.boundaryMesh().whichPatch(curFace);

            if(mesh.isInternalFace(curFace))
            {
                if (curCell == owner[curFace])
                {
                    iDirGradU.append(gradUI[neighbour[curFace]]);
                }
                else
                {
                    iDirGradU.append(gradUI[owner[curFace]]);
                }
            }
            else if
            (
               !isA<tractionPressureDisplacementFvPatchVectorField>
                (
                    U.boundaryField()[patchID]
                )
            )
            {
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
                        iDirGradU.append
                        (
                            weights.boundaryField()[patchID][localFaceID]
                           *(
                                gradUI[cycPatchCells[localFaceID]]
                              - gradUI[cycPatchCells[localFaceID + sizeby2]]
                            )
                          + gradUI[cycPatchCells[localFaceID + sizeby2]]
                        );
                    }
                    else
                    {
                        iDirGradU.append
                        (
                            weights.boundaryField()[patchID][localFaceID]
                           *(
                                gradUI[cycPatchCells[localFaceID]]
                              - gradUI[cycPatchCells[localFaceID - sizeby2]]
                            )
                          + gradUI[cycPatchCells[localFaceID - sizeby2]]
                        );
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

                    iDirGradU.append
                    (
                        weights.boundaryField()[patchID][localFaceID]
                       *(
                            gradUI[procPatchCells[localFaceID]]
                          - gradU.boundaryField()[patchID][localFaceID]
                        )
                      + gradU.boundaryField()[patchID][localFaceID]
                    );
                }
                else if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == emptyPolyPatch::typeName
                )
                {
                    iDirGradU.append(gradUI[curCell]);
                }
                else
                {
                }
            }
        }

        // Weights
        scalarField W(iPoints_[faceI].size(), 1.0);

        tensor avgDirGradU = average(iDirGradU);
        vector avgPoint = average(iPoints_[faceI]);

        label nCoeffs = 3;

        scalar L = max(mag(iPoints_[faceI]-avgPoint));
        Field<tensor> coeffs(nCoeffs, tensor::zero);
        Field<tensor> source(iPoints_[faceI].size(), tensor::zero);

        for (label i=0; i<iPoints_[faceI].size(); i++)
        {
            source[i] = iDirGradU[i] - avgDirGradU;
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += iMatrices_[faceI][i][j]*source[j];
            }
        }

        vector dr = (C[faceI] - avgPoint)/L;

        patchDirGradU[faceI] =
            avgDirGradU
          + coeffs[0]*dr.x()
          + coeffs[1]*dr.y()
          + coeffs[2]*dr.z();
    }

    if (iPoints_.size())
    {
        iPoints_.clear();
        iPointsIndices_.clear();
        iMatrices_.clear();
    }

    return tExDeriv;
}


tmp<vectorField> tractionPressureDisplacementFvPatchVectorField::
extrapolatedDisplacement()
{
    tmp<vectorField> tExDisp
    (
        new vectorField(this->size(), vector::zero)
    );
    vectorField& patchU = tExDisp();

    if (iPoints_.size() == 0)
    {
        makeInterpolationPoints();
    }

    const fvMesh& mesh =
        this->patch().boundaryMesh().mesh();

    word UName =
        this->dimensionedInternalField().name();

    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
            UName
        );

    const cellList& cells = mesh.cells();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const unallocLabelList& patchCells = this->patch().faceCells();

    const surfaceScalarField& weights = mesh.weights();

    const vectorField& UI = U.internalField();

    vectorField n = this->patch().nf();
    const vectorField& C = this->patch().Cf();

    forAll(patchCells, faceI)
    {
        label curCell = patchCells[faceI];

        const labelList& curCellFaces = cells[curCell];
        DynamicList<vector> iU;

        // Add first cell centre point
        iU.append(UI[curCell]);

        // Add face centre points
        forAll(curCellFaces, fI)
        {
            label curFace = curCellFaces[fI];

            label patchID = mesh.boundaryMesh().whichPatch(curFace);

            if(mesh.isInternalFace(curFace))
            {
                if (curCell == owner[curFace])
                {
                    iU.append(UI[neighbour[curFace]]);
                }
                else
                {
                    iU.append(UI[owner[curFace]]);
                }
            }
            else if
            (
               !isA<tractionPressureDisplacementFvPatchVectorField>
                (
                    U.boundaryField()[patchID]
                )
            )
            {
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
                        iU.append
                        (
                            weights.boundaryField()[patchID][localFaceID]
                           *(
                                UI[cycPatchCells[localFaceID]]
                              - UI[cycPatchCells[localFaceID + sizeby2]]
                            )
                          + UI[cycPatchCells[localFaceID + sizeby2]]
                        );
                    }
                    else
                    {
                        iU.append
                        (
                            weights.boundaryField()[patchID][localFaceID]
                           *(
                                UI[cycPatchCells[localFaceID]]
                              - UI[cycPatchCells[localFaceID - sizeby2]]
                            )
                          + UI[cycPatchCells[localFaceID - sizeby2]]
                        );
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

                    iU.append(U.boundaryField()[patchID][localFaceID]);
                }
                else if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == wedgePolyPatch::typeName
                )
                {
                    const wedgePolyPatch& wpp =
                        refCast<const wedgePolyPatch>
                        (
                            mesh.boundaryMesh()[patchID]
                        );

                    vector curU =
                        transform
                        (
                            wpp.cellT(),
                            UI[curCell]
                        );

                    iU.append(curU);
                }
                else if
                (
                    mesh.boundaryMesh()[patchID].type()
                 == emptyPolyPatch::typeName
                )
                {
                    iU.append(UI[curCell]);
                }
                else
                {
                    // Normal patches
                    label start = mesh.boundaryMesh()[patchID].start();
                    label localFaceID = curFace - start;

                    iU.append(U.boundaryField()[patchID][localFaceID]);
                }
            }
        }

        // Weights
        scalarField W(iPoints_[faceI].size(), 1.0);

        // Weighted averate
        vector avgU = average(iU);
        vector avgPoint = average(iPoints_[faceI]);

        label nCoeffs = 3;

        scalar L = max(mag(iPoints_[faceI]-avgPoint));

        Field<vector> coeffs(nCoeffs, vector::zero);
        Field<vector> source(iPoints_[faceI].size(), vector::zero);

        for (label i=0; i<iPoints_[faceI].size(); i++)
        {
            source[i] = iU[i] - avgU;
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += iMatrices_[faceI][i][j]*source[j];
            }
        }

        vector dr = (C[faceI] - avgPoint)/L;

        patchU[faceI] =
            avgU
          + coeffs[0]*dr.x()
          + coeffs[1]*dr.y()
          + coeffs[2]*dr.z();
    }

    if (iPoints_.size())
    {
        iPoints_.clear();
        iPointsIndices_.clear();
        iMatrices_.clear();
    }

    return tExDisp;
}

// XXX: never used?
tmp<vectorField> tractionPressureDisplacementFvPatchVectorField::
patchValue() const
{
    tmp<vectorField> tPatchValue
    (
        new vectorField(this->size(), vector::zero)
    );
    vectorField& pValue = tPatchValue();

    const fvMesh& mesh =
        this->patch().boundaryMesh().mesh();

    word UName =
        this->dimensionedInternalField().name();

    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
            UName
        );

    const volTensorField& gradU =
        mesh.lookupObject<volTensorField>
        (
            "grad(" + UName + ")"
        );

    const unallocLabelList& patchCells =
        this->patch().faceCells();

    const vectorField n = this->patch().nf();
    const vectorField& C = mesh.cellCentres();
    const vectorField& Cf = this->patch().Cf();

    forAll(pValue, faceI)
    {
        const label P = patchCells[faceI];
        const label PP = PP_[faceI];

        const scalar xPP = mag(C[PP] - C[P]);
        const scalar xf = -mag(Cf[faceI] - C[P]);

        const vector UP = U[P];
        const vector UPP = U[PP];

        const vector nGradUP = -(n[faceI] & gradU[P]);
        const vector nGradUPP = -(n[faceI] & gradU[PP]);

        const vector a =
            3*(UPP - UP)/sqr(xPP)
          - (nGradUPP + 2*nGradUP)/xPP;

        const vector b =
            (nGradUPP + nGradUP)/sqr(xPP)
          - 2*(UPP - UP)/pow(xPP, 3);

        pValue[faceI] = UP + nGradUP*xf + a*sqr(xf) + b*pow(xf, 3);
    }

    return tPatchValue;
}

// XXX: never used?
tmp<vectorField> tractionPressureDisplacementFvPatchVectorField::
patchGradient() const
{
    tmp<vectorField> tPatchGradient
    (
        new vectorField(this->size(), vector::zero)
    );
    vectorField& pGradient = tPatchGradient();

    const fvMesh& mesh =
        this->patch().boundaryMesh().mesh();

    word UName =
        this->dimensionedInternalField().name();

    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
            UName
        );

    const unallocLabelList& patchCells =
        this->patch().faceCells();

    const vectorField n = this->patch().nf();
    const vectorField& C = mesh.cellCentres();

    forAll(pGradient, faceI)
    {
        const label P = patchCells[faceI];

        const label PP = PP_[faceI];

        const scalar xPP = mag(C[PP] - C[P]);

        const vector UP = U[P];
        const vector UPP = U[PP];

        pGradient[faceI] = (UP - UPP)/xPP;
    }

    return tPatchGradient;
}

// XXX: never used?
tmp<vectorField> tractionPressureDisplacementFvPatchVectorField::
patchGradient2() const
{
    tmp<vectorField> tPatchGradient
    (
        new vectorField(this->size(), vector::zero)
    );
    vectorField& pGradient = tPatchGradient();

    const fvMesh& mesh =
        this->patch().boundaryMesh().mesh();

    word UName =
        this->dimensionedInternalField().name();

    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
            UName
        );

    const fvsPatchVectorField& Uf =
        patch().lookupPatchField<surfaceVectorField, vector>("Uf");

    const unallocLabelList& patchCells =
        this->patch().faceCells();

    const vectorField n = this->patch().nf();
    const vectorField& C = mesh.cellCentres();
    const vectorField& Cf = this->patch().Cf();

    forAll(pGradient, faceI)
    {
        const label P = patchCells[faceI];

        const label PP = PP_[faceI];


        const scalar xP = mag(C[P] - Cf[faceI]);
        const scalar xPP = xP + mag(C[PP] - C[P]);

        const vector Ub = Uf[faceI];
        const vector UP = U[P];
        const vector UPP = U[PP];

        const vector a =
            ((UPP-Ub)*sqr(xP) - (UP-Ub)*sqr(xPP))/(xP*xPP*(xP-xPP));

        const vector b = (UP-Ub)/sqr(xP) - a/xP;

        pGradient[faceI] = -b;
    }

    return tPatchGradient;
}

// XXX: never used?
tmp<vectorField> tractionPressureDisplacementFvPatchVectorField::
valueBndCoeffs() const
{
    tmp<vectorField> tValBndCoeffs
    (
        new vectorField(this->size(), vector::zero)
    );
    vectorField& valBndCoeffs = tValBndCoeffs();

    const fvMesh& mesh =
        this->patch().boundaryMesh().mesh();

    word UName =
        this->dimensionedInternalField().name();

    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
            UName
        );

    const volTensorField& gradU =
        mesh.lookupObject<volTensorField>
        (
            "grad(" + UName + ")"
        );

    const unallocLabelList& patchCells =
        this->patch().faceCells();

    const vectorField n = this->patch().nf();
    const vectorField& C = mesh.cellCentres();
    const vectorField& Cf = this->patch().Cf();

    forAll(valBndCoeffs, faceI)
    {
        const label P = patchCells[faceI];
        const label PP = PP_[faceI];

        const scalar xPP = mag(C[PP] - C[P]);
        const scalar xf = -mag(Cf[faceI] - C[P]);

        const vector UP = U[P];
        const vector UPP = U[PP];

        const vector nGradUP = -(n[faceI] & gradU[P]);
        const vector nGradUPP = -(n[faceI] & gradU[PP]);

        const vector a =
            3*(UPP - UP)/sqr(xPP)
          - (nGradUPP + 2*nGradUP)/xPP;

        const vector b =
            (nGradUPP + nGradUP)/sqr(xPP)
          - 2*(UPP - UP)/pow(xPP, 3);

        valBndCoeffs[faceI] = nGradUP*xf + a*sqr(xf) + b*pow(xf, 3);
    }

    return tValBndCoeffs;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
