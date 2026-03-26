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

#include "alphaStab.H"
#include "fvc.H"
#include "compatibilityFunctions.H"
#include "deltaVectors.H"
#include "symmTensor3rdOrder.H"
#include "processorPolyPatch.H"
#include "processorFvPatch.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
inline void Foam::alphaStab::exchangeProcessorPatchField
(
    const processorFvPatch& procPatch,
    const Field<Type>& sendField,
    Field<Type>& recvField
) const
{
    if (Pstream::myProcNo() < procPatch.neighbProcNo())
    {
        procPatch.send(Pstream::commsTypes::blocking, sendField);
        procPatch.receive(Pstream::commsTypes::blocking, recvField);
    }
    else
    {
        procPatch.receive(Pstream::commsTypes::blocking, recvField);
        procPatch.send(Pstream::commsTypes::blocking, sendField);
    }
}


inline void Foam::alphaStab::computeDiffStencilHighOrderScalar
(
    const volScalarField& field,
    const volVectorField& gradField,
    surfaceScalarField& faceFlux,
    scalar scaleFactor
) const
{
    const_cast<volScalarField&>(field).correctBoundaryConditions();
    const_cast<volVectorField&>(gradField).correctBoundaryConditions();

    // Take references for brevity and efficiency
    const fvMesh& mesh = field.mesh();
    const vectorField& C = mesh.C();
    const vectorField& Cf = mesh.Cf();
    const scalarField& fieldI = field.internalField();
    const vectorField& gradFieldI = gradField.internalField();
    scalarField& faceFluxI = Foam::primitiveFieldRef(faceFlux);
    const solidModel& solMod = lookupSolidModel(field.mesh());
    const label order = solMod.MLS().order();
    const vectorField deltaI(deltaVectors(mesh));
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());
    const vectorField& nI = n.internalField();

    // owner/neighbour addressing
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    // Initialise storage for second gradint
    tmp<volSymmTensorField> tSecondGrad;
    const symmTensorField* secondGradI = nullptr;

    // Initialise storage for third gradient
    autoPtr<List<symmTensor3rdOrder>> tThirdGrad;
    const List<symmTensor3rdOrder>* thirdGrad = nullptr;

    // Second gradient is computed on demand by MLS
    if (order >= 2)
    {
        tSecondGrad = solMod.MLS().secondGrad(field);
        tmpRef(tSecondGrad).correctBoundaryConditions();
        secondGradI = &tSecondGrad().internalField();
    }

    // Third gradient is computed on demand by MLS
    if (order >= 3)
    {
        tThirdGrad = solMod.MLS().thirdGrad(field);
        thirdGrad = &tThirdGrad();
    }

    forAll(owner, facei)
    {
        // Delta vector and face normal
        const vector& d = deltaI[facei];
        const vector& nf = nI[facei];

        // Owner and neighbour cell IDs
        const label ownCellID = owner[facei];
        const label neiCellID = neighbour[facei];

        // Vectors from cell centres to the face
        const vector dL = Cf[facei] - C[ownCellID];
        const vector dR = Cf[facei] - C[neiCellID];

        // Linear extrapolation from the owner to the face
        scalar uL = fieldI[ownCellID] + (dL & gradFieldI[ownCellID]);

        // Linear extrapolation from the neighbour to the face
        scalar uR = fieldI[neiCellID] + (dR & gradFieldI[neiCellID]);

        // Quadratic extrapolation term
        if (order >= 2)
        {
            uL += 0.5*(dL & ((*secondGradI)[ownCellID] & dL));
            uR += 0.5*(dR & ((*secondGradI)[neiCellID] & dR));
        }

        // Cubic extrapolation term
        if (order >= 3)
        {
            const scalar oneOverSix = 1.0/6.0;
            uL += oneOverSix*cubicForm((*thirdGrad)[ownCellID], dL);
            uR += oneOverSix*cubicForm((*thirdGrad)[neiCellID], dR);
        }

        const scalar denom = max(mag(nf & d), VSMALL);

        // Stabilisation based on the jump across the face
        faceFluxI[facei] = scaleFactor*(uR - uL)/denom;
    }

    auto& faceFluxfBf = Foam::boundaryFieldRef(faceFlux);

    forAll(faceFlux.boundaryField(), patchi)
    {
        scalarField& faceFluxfBfI = faceFluxfBf[patchi];
        const fvPatchScalarField& pvf = field.boundaryField()[patchi];
        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchCf = mesh.boundary()[patchi].Cf();
        tmp<vectorField> tPatchN = mesh.boundary()[patchi].nf();
        const vectorField& patchN = tPatchN();

        if (pvf.coupled())
        {
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(mesh.boundary()[patchi]);
            const vectorField& neiProcC =
                procPatch.procPolyPatch().neighbFaceCellCentres();

            const scalarField neiProcField
            (
                field.boundaryField()[patchi].patchNeighbourField()
            );
            const vectorField neiProcGrad
            (
                gradField.boundaryField()[patchi].patchNeighbourField()
            );
            const symmTensorField neiProcSecondGrad
            (
                tSecondGrad().boundaryField()[patchi].patchNeighbourField()
            );

            Field<symmTensor3rdOrder> sendThirdGrad(faceCells.size());
            Field<symmTensor3rdOrder> neiProcThirdGrad(faceCells.size());
            if (order >= 3)
            {
                forAll(faceCells, facei)
                {
                    const label ownCellID = faceCells[facei];
                    sendThirdGrad[facei] = (*thirdGrad)[ownCellID];
                }

                // patchNeighbourField() for symmTensor3rdOrder field
                exchangeProcessorPatchField
                (
                    procPatch,
                    sendThirdGrad,
                    neiProcThirdGrad
                );
            }

            forAll(pvf, facei)
            {
                // Owner cell IDs
                const label ownCellID = faceCells[facei];
                const vector d = neiProcC[facei] - C[ownCellID];
                const vector normal = patchN[facei];

                // Vectors from cell centres to the face
                const vector dL = patchCf[facei] - C[ownCellID];
                const vector dR = patchCf[facei] - neiProcC[facei];

                // Linear extrapolate from the owner to the face
                scalar uL = fieldI[ownCellID] + (dL & gradFieldI[ownCellID]);

                // Linear extrapolate from the neighbour to the face
                scalar uR = neiProcField[facei] + (dR & neiProcGrad[facei]);

                // Quadratic extrapolation term
                if (order >= 2)
                {
                    uL += 0.5*(dL & ((*secondGradI)[ownCellID] & dL));
                    uR += 0.5*(dR & (neiProcSecondGrad[facei] & dR));
                }

                // Cubic extrapolation term
                if (order >= 3)
                {
                    const scalar oneOverSix = 1.0/6.0;
                    uL += oneOverSix*cubicForm((*thirdGrad)[ownCellID], dL);
                    uR += oneOverSix*cubicForm(neiProcThirdGrad[facei], dR);
                }

                const scalar denom = max(mag(normal & d), VSMALL);

                // Stabilisation based on the jump across the face
                faceFluxfBfI[facei] = scaleFactor*(uR - uL)/denom;
            }
        }
        else if (pvf.fixesValue())
        {
            fvPatchScalarField& fieldPatch =
                const_cast<volScalarField&>(field).boundaryFieldRef()[patchi];
            fieldPatch.updateCoeffs();
            scalarField patchField(fieldPatch);

            forAll(pvf, facei)
            {
                // Owner cell IDs
                const label ownCellID = faceCells[facei];

                // Vectors from cell centres to the face
                const vector dL = patchCf[facei] - C[ownCellID];
                const vector normal = patchN[facei];

                // Linear extrapolate from the owner to the face
                scalar uL = fieldI[ownCellID] + (dL & gradFieldI[ownCellID]);

                // Use the actual boundary value
                const scalar uR = patchField[facei];

                // Quadratic extrapolation term
                if (order >= 2)
                {
                    uL += 0.5*(dL & ((*secondGradI)[ownCellID] & dL));
                }

                // Cubic extrapolation term
                if (order >= 3)
                {
                    const scalar oneOverSix = 1.0/6.0;
                    uL += oneOverSix*cubicForm((*thirdGrad)[ownCellID], dL);
                }

                const scalar denom = max(mag(normal & dL), VSMALL);

                // Stabilisation based on the jump across the face
                faceFluxfBfI[facei] = scaleFactor*(uR - uL)/denom;
            }
        }
        else
        {
            faceFluxfBfI = 0.0;
        }
    }
}

inline void Foam::alphaStab::computeDiffStencilHighOrderVector
(
    const volVectorField& field,
    const volTensorField& gradField,
    surfaceVectorField& faceFlux,
    scalar scaleFactor
) const
{
    const_cast<volVectorField&>(field).correctBoundaryConditions();
    const_cast<volTensorField&>(gradField).correctBoundaryConditions();

    // Take references for brevity and efficiency
    const fvMesh& mesh = field.mesh();
    const vectorField& C = mesh.C();
    const vectorField& Cf = mesh.Cf();
    const vectorField& fieldI = field.internalField();
    const tensorField& gradFieldI = gradField.internalField();
    vectorField& faceFluxI = Foam::primitiveFieldRef(faceFlux);
    const solidModel& solMod = lookupSolidModel(field.mesh());
    const label order = solMod.MLS().order();
    const vectorField deltaI(deltaVectors(mesh));
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());
    const vectorField& nI = n.internalField();

    // owner/neighbour addressing
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    // Decompose vector field into scalar components
    volScalarField fieldX(field.component(vector::X));
    volScalarField fieldY(field.component(vector::Y));
    volScalarField fieldZ(field.component(vector::Z));

    // Initialise storage for second gradint
    tmp<volSymmTensorField> tSecondGradX;
    tmp<volSymmTensorField> tSecondGradY;
    tmp<volSymmTensorField> tSecondGradZ;
    const symmTensorField* secondGradXI = nullptr;
    const symmTensorField* secondGradYI = nullptr;
    const symmTensorField* secondGradZI = nullptr;

    // Initialise storage for third gradient
    autoPtr<List<symmTensor3rdOrder>> tThirdGradX;
    autoPtr<List<symmTensor3rdOrder>> tThirdGradY;
    autoPtr<List<symmTensor3rdOrder>> tThirdGradZ;
    const List<symmTensor3rdOrder>* thirdGradX = nullptr;
    const List<symmTensor3rdOrder>* thirdGradY = nullptr;
    const List<symmTensor3rdOrder>* thirdGradZ = nullptr;

    // Second gradient is computed on demand by MLS
    if (order >= 2)
    {
        tSecondGradX = solMod.MLS().secondGrad(fieldX);
        tSecondGradY = solMod.MLS().secondGrad(fieldY);
        tSecondGradZ = solMod.MLS().secondGrad(fieldZ);
        tmpRef(tSecondGradX).correctBoundaryConditions();
        tmpRef(tSecondGradY).correctBoundaryConditions();
        tmpRef(tSecondGradZ).correctBoundaryConditions();

        secondGradXI = &tSecondGradX().internalField();
        secondGradYI = &tSecondGradY().internalField();
        secondGradZI = &tSecondGradZ().internalField();
    }

    // Third gradient is computed on demand by MLS
    if (order >= 3)
    {
        tThirdGradX = solMod.MLS().thirdGrad(fieldX);
        tThirdGradY = solMod.MLS().thirdGrad(fieldY);
        tThirdGradZ = solMod.MLS().thirdGrad(fieldZ);

        thirdGradX = &tThirdGradX();
        thirdGradY = &tThirdGradY();
        thirdGradZ = &tThirdGradZ();
    }

    forAll(owner, facei)
    {
        // Delta vector and face normal
        const vector& d = deltaI[facei];
        const vector& nf = nI[facei];

        // Owner and neighbour cell IDs
        const label ownCellID = owner[facei];
        const label neiCellID = neighbour[facei];

        // Vectors from cell centres to the face
        const vector dL = Cf[facei] - C[ownCellID];
        const vector dR = Cf[facei] - C[neiCellID];

        // Linear extrapolation from the owner to the face
        vector uL = fieldI[ownCellID] + (dL & gradFieldI[ownCellID]);

        // Linear extrapolation from the neighbour to the face
        vector uR = fieldI[neiCellID] + (dR & gradFieldI[neiCellID]);

        // Quadratic extrapolation term
        if (order >= 2)
        {
            uL.x() += 0.5*(dL & ((*secondGradXI)[ownCellID] & dL));
            uL.y() += 0.5*(dL & ((*secondGradYI)[ownCellID] & dL));
            uL.z() += 0.5*(dL & ((*secondGradZI)[ownCellID] & dL));

            uR.x() += 0.5*(dR & ((*secondGradXI)[neiCellID] & dR));
            uR.y() += 0.5*(dR & ((*secondGradYI)[neiCellID] & dR));
            uR.z() += 0.5*(dR & ((*secondGradZI)[neiCellID] & dR));
        }

        // Cubic extrapolation term
        if (order >= 3)
        {
            const scalar oneOverSix = 1.0/6.0;
            uL.x() += oneOverSix*cubicForm((*thirdGradX)[ownCellID], dL);
            uL.y() += oneOverSix*cubicForm((*thirdGradY)[ownCellID], dL);
            uL.z() += oneOverSix*cubicForm((*thirdGradZ)[ownCellID], dL);

            uR.x() += oneOverSix*cubicForm((*thirdGradX)[neiCellID], dR);
            uR.y() += oneOverSix*cubicForm((*thirdGradY)[neiCellID], dR);
            uR.z() += oneOverSix*cubicForm((*thirdGradZ)[neiCellID], dR);
        }

        const scalar denom = max(mag(nf & d), VSMALL);

        // Stabilisation based on the jump across the face
        faceFluxI[facei] = scaleFactor*(uR - uL)/denom;
    }

    auto& faceFluxfBf = Foam::boundaryFieldRef(faceFlux);

    forAll(faceFlux.boundaryField(), patchi)
    {
        vectorField& faceFluxfBfI = faceFluxfBf[patchi];
        const fvPatchVectorField& pvf = field.boundaryField()[patchi];
        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchCf = mesh.boundary()[patchi].Cf();
        tmp<vectorField> tPatchN = mesh.boundary()[patchi].nf();
        const vectorField& patchN = tPatchN();

        if (pvf.coupled())
        {
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(mesh.boundary()[patchi]);
            const vectorField& neiProcC =
                procPatch.procPolyPatch().neighbFaceCellCentres();

            const vectorField neiProcField
            (
                field.boundaryField()[patchi].patchNeighbourField()
            );
            const tensorField neiProcGrad
            (
                gradField.boundaryField()[patchi].patchNeighbourField()
            );
            const symmTensorField neiProcSecondGradX
            (
                tSecondGradX().boundaryField()[patchi].patchNeighbourField()
            );
            const symmTensorField neiProcSecondGradY
            (
                tSecondGradY().boundaryField()[patchi].patchNeighbourField()
            );
            const symmTensorField neiProcSecondGradZ
            (
                tSecondGradZ().boundaryField()[patchi].patchNeighbourField()
            );
            Field<symmTensor3rdOrder> sendThirdGradX(faceCells.size());
            Field<symmTensor3rdOrder> sendThirdGradY(faceCells.size());
            Field<symmTensor3rdOrder> sendThirdGradZ(faceCells.size());

            forAll(faceCells, facei)
            {
                const label ownCellID = faceCells[facei];
                if (order >= 3)
                {
                    sendThirdGradX[facei] = (*thirdGradX)[ownCellID];
                    sendThirdGradY[facei] = (*thirdGradY)[ownCellID];
                    sendThirdGradZ[facei] = (*thirdGradZ)[ownCellID];
                }
            }

            Field<symmTensor3rdOrder> neiProcThirdGradX(faceCells.size());
            Field<symmTensor3rdOrder> neiProcThirdGradY(faceCells.size());
            Field<symmTensor3rdOrder> neiProcThirdGradZ(faceCells.size());

            // patchNeighbourField()x for symmTensor3rdOrder field
            if (order >= 3)
            {
                exchangeProcessorPatchField
                (
                    procPatch,
                    sendThirdGradX,
                    neiProcThirdGradX
                );
                exchangeProcessorPatchField
                (
                    procPatch,
                    sendThirdGradY,
                    neiProcThirdGradY
                );
                exchangeProcessorPatchField
                (
                    procPatch,
                    sendThirdGradZ,
                    neiProcThirdGradZ
                );
            }
            forAll(pvf, facei)
            {
                // Owner cell IDs
                const label ownCellID = faceCells[facei];
                const vector d = neiProcC[facei] - C[ownCellID];
                const vector normal = patchN[facei];

                // Vectors from cell centres to the face
                const vector dL = patchCf[facei] - C[ownCellID];
                const vector dR = patchCf[facei] - neiProcC[facei];

                // Linear extrapolate from the owner to the face
                vector uL = fieldI[ownCellID] + (dL & gradFieldI[ownCellID]);

                // Linear extrapolate from the neighbour to the face
                vector uR = neiProcField[facei] + (dR & neiProcGrad[facei]);

                // Quadratic extrapolation term
                if (order >= 2)
                {
                    uL.x() += 0.5*(dL & ((*secondGradXI)[ownCellID] & dL));
                    uL.y() += 0.5*(dL & ((*secondGradYI)[ownCellID] & dL));
                    uL.z() += 0.5*(dL & ((*secondGradZI)[ownCellID] & dL));

                    uR.x() += 0.5*(dR & (neiProcSecondGradX[facei] & dR));
                    uR.y() += 0.5*(dR & (neiProcSecondGradY[facei] & dR));
                    uR.z() += 0.5*(dR & (neiProcSecondGradZ[facei] & dR));
                }

                // Cubic extrapolation term
                if (order >= 3)
                {
                    const scalar oneOverSix = 1.0/6.0;
                    uL.x() += oneOverSix*cubicForm((*thirdGradX)[ownCellID], dL);
                    uL.y() += oneOverSix*cubicForm((*thirdGradY)[ownCellID], dL);
                    uL.z() += oneOverSix*cubicForm((*thirdGradZ)[ownCellID], dL);

                    uR.x() += oneOverSix*cubicForm(neiProcThirdGradX[facei], dR);
                    uR.y() += oneOverSix*cubicForm(neiProcThirdGradY[facei], dR);
                    uR.z() += oneOverSix*cubicForm(neiProcThirdGradZ[facei], dR);
                }

                const scalar denom = max(mag(normal & d), VSMALL);

                // Stabilisation based on the jump across the face
                faceFluxfBfI[facei] = scaleFactor*(uR - uL)/denom;
            }
        }
        else if (pvf.fixesValue())
        {
            fvPatchVectorField& fieldPatch =
                const_cast<volVectorField&>(field).boundaryFieldRef()[patchi];
            fieldPatch.updateCoeffs();
            vectorField patchField(fieldPatch);

            forAll(pvf, facei)
            {
                // Owner cell IDs
                const label ownCellID = faceCells[facei];

                // Vectors from cell centres to the face
                const vector dL = patchCf[facei] - C[ownCellID];
                const vector normal = patchN[facei];

                // Linear extrapolate from the owner to the face
                vector uL = fieldI[ownCellID] + (dL & gradFieldI[ownCellID]);

                // Use the actual boundary value
                const vector uR = patchField[facei];

                // Quadratic extrapolation term
                if (order >= 2)
                {
                    uL.x() += 0.5*(dL & ((*secondGradXI)[ownCellID] & dL));
                    uL.y() += 0.5*(dL & ((*secondGradYI)[ownCellID] & dL));
                    uL.z() += 0.5*(dL & ((*secondGradZI)[ownCellID] & dL));
                }

                // Cubic extrapolation term
                if (order >= 3)
                {
                    const scalar oneOverSix = 1.0/6.0;
                    uL.x() += oneOverSix*cubicForm((*thirdGradX)[ownCellID], dL);
                    uL.y() += oneOverSix*cubicForm((*thirdGradY)[ownCellID], dL);
                    uL.z() += oneOverSix*cubicForm((*thirdGradZ)[ownCellID], dL);
                }

                const scalar denom = max(mag(normal & dL), VSMALL);

                // Stabilisation based on the jump across the face
                faceFluxfBfI[facei] = scaleFactor*(uR - uL)/denom;
            }
        }
        else
        {
            faceFluxfBfI = vector::zero;
        }
    }
}


template<class GeoField, class GradField, class SurfaceField>
void Foam::alphaStab::computeDiffStencil
(
    const GeoField& field,
    const GradField& gradField,
    SurfaceField& faceFlux,
    scalar scaleFactor
) const
{
    // Value and gradient types, e.g. scalar, vector
    using Type = typename GeoField::value_type;
    using GradType = typename GradField::value_type;

    // Take references for brevity and efficiency
    const fvMesh& mesh = field.mesh();
    const vectorField& C = mesh.C();
    const vectorField& Cf = mesh.Cf();
    const Field<Type>& fieldI = field.internalField();
    const Field<GradType>& gradFieldI = gradField.internalField();
    Field<Type>& faceFluxI = Foam::primitiveFieldRef(faceFlux);

    // Unit normals
    //const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // set reference to difference factors array
    const scalarField& deltaCoeffs = mesh.deltaCoeffs();

    // owner/neighbour addressing
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        // Owner and neighbour cell IDs
        const label ownCellID = owner[facei];
        const label neiCellID = neighbour[facei];

        // Vectors from cell centres to the face
        const vector dL = Cf[facei] - C[ownCellID];
        const vector dR = Cf[facei] - C[neiCellID];

        // Extrapolate from the owner to the face
        const Type uL = fieldI[ownCellID] + (dL & gradFieldI[ownCellID]);

        // Extrapolate from the neighbour to the face
        const Type uR = fieldI[neiCellID] + (dR & gradFieldI[neiCellID]);

        // Stabilisation based on the jump across the face
        faceFluxI[facei] = scaleFactor_*deltaCoeffs[facei]*(uR - uL);
    }

    auto& faceFluxfBf = Foam::boundaryFieldRef(faceFlux);

    forAll(faceFlux.boundaryField(), patchi)
    {
        Field<Type>& faceFluxfBfI = faceFluxfBf[patchi];
        const fvPatchField<Type>& pvf = field.boundaryField()[patchi];
        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchCf = mesh.boundary()[patchi].Cf();
        const scalarField& patchDeltaCoeffs =
            mesh.boundary()[patchi].deltaCoeffs();

        if (pvf.coupled())
        {
            const Field<Type> neiProcField(pvf.patchNeighbourField());
            const Field<GradType> neiProcGrad
            (
                gradField.boundaryField()[patchi].patchNeighbourField()
            );
            const vectorField neiProcC
            (
                mesh.C().boundaryField()[patchi].patchNeighbourField()
            );

            forAll(pvf, facei)
            {
                // Owner cell IDs
                const label ownCellID = faceCells[facei];

                // Vectors from cell centres to the face
                const vector dL = patchCf[facei] - C[ownCellID];
                const vector dR = patchCf[facei] - neiProcC[facei];

                // Extrapolate from the owner to the face
                const Type uL =
                    fieldI[ownCellID] + (dL & gradFieldI[ownCellID]);

                // Extrapolate from the neighbour (other proc) to the face
                const Type uR =
                    neiProcField[facei] + (dR & neiProcGrad[facei]);

                // Stabilisation based on the jump across the face
                faceFluxfBfI[facei] =
                    scaleFactor_*patchDeltaCoeffs[facei]*(uR - uL);
            }
        }
        else
        {
            forAll(pvf, facei)
            {
                // Owner cell IDs
                const label ownCellID = faceCells[facei];

                // Vectors from cell centres to the face
                const vector dL = patchCf[facei] - C[ownCellID];

                // Extrapolate from the owner to the face
                const Type uL = fieldI[ownCellID] + (dL & gradFieldI[ownCellID]);

                // Use the actual boundary value on uncoupled patches
                const Type& uR = pvf[facei];

                // Stabilisation based on the jump across the face
                faceFluxfBfI[facei] =
                    scaleFactor_*patchDeltaCoeffs[facei]*(uR - uL);
            }
        }
    }

#ifdef OPENFOAM_COM
    faceFlux.correctBoundaryConditions();
#endif
}


// ************************************************************************* //
