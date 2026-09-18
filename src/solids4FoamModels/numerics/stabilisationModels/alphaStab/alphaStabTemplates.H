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
#include "processorFvPatch.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class GeoField, class SurfaceField>
void Foam::alphaStab::computeDiffStencilHighOrder
(
    const GeoField& field,
    SurfaceField& faceFlux,
    const scalar scaleFactor,
    const leastSquaresScheme& reconstruction
) const
{
    using Type = typename GeoField::value_type;

    const_cast<GeoField&>(field).correctBoundaryConditions();

    const fvMesh& mesh = field.mesh();
    const vectorField& C = mesh.C();
    Field<Type>& faceFluxI = Foam::primitiveFieldRef(faceFlux);
    const vectorField deltaI(deltaVectors(mesh));
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());
    const vectorField& nI = n.internalField();
    const labelUList& owner = mesh.owner();

    SurfaceField ownerFaceValue
    (
        IOobject
        (
            "ownerFaceCentreValue(" + field.name() + ")",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<Type>
        (
            "zero",
            field.dimensions(),
            pTraits<Type>::zero
        )
    );
    SurfaceField neighbourFaceValue
    (
        IOobject
        (
            "neighbourFaceCentreValue(" + field.name() + ")",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<Type>
        (
            "zero",
            field.dimensions(),
            pTraits<Type>::zero
        )
    );

    reconstruction.faceCentreValues
    (
        field,
        ownerFaceValue,
        neighbourFaceValue
    );

    const Field<Type>& ownerFaceValueI = ownerFaceValue.internalField();
    const Field<Type>& neighbourFaceValueI =
        neighbourFaceValue.internalField();

    forAll(owner, facei)
    {
        const scalar denom =
            max(mag(nI[facei] & deltaI[facei]), VSMALL);

        faceFluxI[facei] =
            scaleFactor
           *(neighbourFaceValueI[facei] - ownerFaceValueI[facei])
           /denom;
    }

    auto& faceFluxBf = Foam::boundaryFieldRef(faceFlux);
    const auto& ownerFaceValueBf = ownerFaceValue.boundaryField();
    const auto& neighbourFaceValueBf = neighbourFaceValue.boundaryField();

    forAll(faceFlux.boundaryField(), patchi)
    {
        Field<Type>& patchFaceFlux = faceFluxBf[patchi];
        const fvPatchField<Type>& pvf = field.boundaryField()[patchi];
        const fvPatch& patch = mesh.boundary()[patchi];
        const Field<Type>& patchOwnerValue = ownerFaceValueBf[patchi];
        tmp<vectorField> tPatchN = patch.nf();
        const vectorField& patchN = tPatchN();

        if (isA<processorFvPatch>(patch))
        {
#ifdef FOAMEXTEND
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patch.patch());
            const vectorField& neighbourCellCentres =
                procPatch.neighbFaceCellCentres();
#else
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(patch);
            const vectorField& neighbourCellCentres =
                procPatch.procPolyPatch().neighbFaceCellCentres();
#endif
            const labelUList& faceCells = patch.faceCells();
            const Field<Type>& patchNeighbourValue =
                neighbourFaceValueBf[patchi];

            forAll(patch, facei)
            {
                const vector d =
                    neighbourCellCentres[facei] - C[faceCells[facei]];
                const scalar denom =
                    max(mag(patchN[facei] & d), VSMALL);

                patchFaceFlux[facei] =
                    scaleFactor
                   *(patchNeighbourValue[facei] - patchOwnerValue[facei])
                   /denom;
            }
        }
        else if (pvf.fixesValue())
        {
            const vectorField& patchCf = patch.Cf();
            const labelUList& faceCells = patch.faceCells();

            forAll(patch, facei)
            {
                const vector dL =
                    patchCf[facei] - C[faceCells[facei]];
                const scalar denom =
                    max(mag(patchN[facei] & dL), VSMALL);

                patchFaceFlux[facei] =
                    scaleFactor
                   *(pvf[facei] - patchOwnerValue[facei])
                   /denom;
            }
        }
        else
        {
            patchFaceFlux = pTraits<Type>::zero;
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
