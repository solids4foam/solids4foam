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

#include "newtonCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "directMapInterfaceToInterfaceMapping.H"
#include "fixedValueFvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"
#include "newtonIcoFluid.H"
#include "dynamicMotionSolverFvMesh.H"
#include "motionSolver.H"
#include "meshMotionSolidModelFvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(newtonCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface, newtonCouplingInterface, dictionary
);


// * * * * * * * * * * * Private Member Functions* * * * * * * * * * * * * * //


foamPetscSnesHelper& newtonCouplingInterface::motion()
{
    // meshMotionSolidModelFvMotionSolver is the only currently supported
    // mesh motion solver

    if (!isA<dynamicMotionSolverFvMesh>(fluidMesh()))
    {
        FatalErrorInFunction
            << "meshMotionSolidModelFvMotionSolver is the only currently "
            << "supported mesh motion solver" << exit(FatalError);
    }

    if
    (
        !isA<meshMotionSolidModelFvMotionSolver>
        (
            refCast<dynamicMotionSolverFvMesh>(fluidMesh()).motion()
        )
    )
    {
        FatalErrorInFunction
            << "meshMotionSolidModelFvMotionSolver is the only currently "
            << "supported mesh motion solver" << exit(FatalError);
    }

    // Cast the fluid mesh to a dynamicMotionSolverFvMesh then access its motion
    // solver and cast it to a foamPetscSnesHelper
    // This will only suceed if using a motion solver based on a
    // foamPetscSnesHelper solver
    const foamPetscSnesHelper& motion =
        refCast<const foamPetscSnesHelper>
        (
            refCast<const meshMotionSolidModelFvMotionSolver>
            (
                refCast<dynamicMotionSolverFvMesh>(fluidMesh()).motion()
            ).model()
        );

    // Cast away the const-ness (sorry)
    return const_cast<foamPetscSnesHelper&>(motion);
}


solidModel& newtonCouplingInterface::motionSolid()
{
    return refCast<solidModel>(motion());
}


label newtonCouplingInterface::initialiseAfm
(
    Mat Afm,
    const fvMesh& fluidMesh,
    const label fluidBlockSize,
    const label motionBlockSize,
    const bool twoD
) const
{
    if (debug)
    {
        InfoInFunction
            << "Initialising" << endl;
    }

    // The motion appears as the mesh flux in the momentum convection and
    // pressure equation flux divergence terms.
    // This means each fluid equation is coupled to the mesh motion in all
    // cells in the compact stencil

    // CAREFUL: we are setting non-zeros here based on the scalar rows, not
    // the block rows

    // Set matrix type to AIJ (since BAIJ does not support non-square
    // blocks)
    CHKERRQ(MatSetType(Afm, MATMPIAIJ));

    // Total number of scalar rows in Amf (same as in the fluid region)
    const label scalarRowN = fluidMesh.nCells()*fluidBlockSize;

    // Allocate per-scalar-row nonzeros
    // d is initialised to 1 to count the diagonal
    // while o is initialised to 0
    std::vector<int> d_nnz(scalarRowN, fluidBlockSize);
    std::vector<int> o_nnz(scalarRowN, 0);

    // Count neighbours sharing an internal face
    const Foam::labelUList& own = fluidMesh.owner();
    const Foam::labelUList& nei = fluidMesh.neighbour();
    forAll(own, faceI)
    {
        const Foam::label ownCellID = own[faceI];
        const Foam::label neiCellID = nei[faceI];

        // Increase count for all scalar rows
        label ownRowID = ownCellID*fluidBlockSize;
        d_nnz[ownRowID++] += motionBlockSize;
        d_nnz[ownRowID++] += motionBlockSize;
        d_nnz[ownRowID++] += motionBlockSize;

        label neiRowID = neiCellID*fluidBlockSize;
        d_nnz[neiRowID++] += motionBlockSize;
        d_nnz[neiRowID++] += motionBlockSize;
        d_nnz[neiRowID++] += motionBlockSize;

        if (!twoD)
        {
            d_nnz[ownRowID++] += motionBlockSize;
            d_nnz[neiRowID] += motionBlockSize;
        }
    }

    // Count off-processor neighbour cells
    forAll(fluidMesh.boundary(), patchI)
    {
        if (fluidMesh.boundary()[patchI].type() == "processor")
        {
            const Foam::labelUList& faceCells =
                fluidMesh.boundary()[patchI].faceCells();

            forAll(faceCells, fcI)
            {
                const Foam::label cellID = faceCells[fcI];
                label rowID = cellID*fluidBlockSize;
                o_nnz[rowID] += motionBlockSize;
                o_nnz[rowID] += motionBlockSize;
                o_nnz[rowID] += motionBlockSize;

                if (!twoD)
                {
                    o_nnz[rowID] += motionBlockSize;
                }
            }
        }
        else if (fluidMesh.boundary()[patchI].coupled())
        {
            // Other coupled boundaries are not implemented
            FatalErrorInFunction
                << "Coupled boundary are not implemented, except for"
                << " processor boundaries" << exit(FatalError);
        }
    }

    // Allocate parallel matrix using AIJ
    CHKERRQ
    (
        MatMPIAIJSetPreallocation(Afm, 0, d_nnz.data(), 0, o_nnz.data())
    );

    // Raise an error if mallocs are required during matrix assembly
    CHKERRQ(MatSetOption(Afm, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));

    return 0;
}


label newtonCouplingInterface::initialiseAfs
(
    Mat Afs,
    const fvMesh& fluidMesh,
    const label fluidBlockSize,
    const label solidBlockSize,
    const bool twoD
) const
{
    // Initialise Afs solid-in-fluid-coupling matrix
    if (debug)
    {
        InfoInFunction
            << "Initialising" << endl;
    }

    // Initially we assume a conformal FSI interface, where each fluid cell
    // shares a face with a solid cell. So we assume the number of blocks in
    // the Afs (and Asf) matrix is equal to the number of cells at the
    // interface

    // Initially, we only allow one FSI interface
    if (interfaceToInterfaceList().size() != 1)
    {
        FatalErrorInFunction
            << "Currently, only one interface is allowed when using "
            << typeName << abort(FatalError);
    }

    // Allow only a direct map (conformal interface)
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // CAREFUL: we are setting non-zeros here based on the scalar rows, not
    // the block rows

    // Set matrix type to AIJ (since BAIJ does not support non-square
    // blocks)
    CHKERRQ(MatSetType(Afs, MATMPIAIJ));

    // Total number of scalar rows in the fluid region
    const label scalarRowN = fluidMesh.nCells()*fluidBlockSize;

    // Allocate per-scalar-row nonzeros, initialised to 0
    std::vector<int> d_nnz(scalarRowN, 0);
    std::vector<int> o_nnz(scalarRowN, 0);

    // Set non-zeros for each interface fluid cells
    const labelList& fluidFaceMap = interfaceMap.zoneAToZoneBFaceMap();
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const labelList& fluidFaceCells =
        fluidMesh.boundary()[fluidPatchID].faceCells();
    forAll(fluidFaceMap, fluidFaceI)
    {
        const label fluidCellID = fluidFaceCells[fluidFaceI];

        // Calculate the row index for this cells first scalar equation
        label rowI = fluidCellID*fluidBlockSize;

        // Set the number of non-zeros to be the number of solid equations
        // e.g., The x-momentum equation could have a coefficient for the
        // solid x/y/z displacement
        d_nnz[rowI++] = solidBlockSize;
        d_nnz[rowI++] = solidBlockSize;
        d_nnz[rowI++] = solidBlockSize;

        if (!twoD)
        {
            d_nnz[rowI++] = solidBlockSize;
        }
    }

    // Parallel: we need to decide how to deal with this, e.g. do we allow
    // the same general decompositions as in the partitioned approach
    // For now, only allow serial
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Currently, serial runs are allowed in "
            << typeName << abort(FatalError);
    }

    // Allocate parallel matrix using AIJ
    CHKERRQ
    (
        MatMPIAIJSetPreallocation(Afs, 0, d_nnz.data(), 0, o_nnz.data())
    );

    // Raise an error if mallocs are required during matrix assembly
    CHKERRQ(MatSetOption(Afs, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));

    return 0;
}


label newtonCouplingInterface::initialiseAms
(
    Mat Ams,
    const fvMesh& fluidMesh,
    const label motionBlockSize,
    const label solidBlockSize,
    const bool twoD
) const
{
    if (debug)
    {
        InfoInFunction
            << "Initialising" << endl;
    }

    // Initially we assume a conformal FSI interface, where each fluid cell
    // shares a face with a solid cell. So we assume the number of blocks in
    // the Ams matrix is equal to the number of cells at the
    // interface

    // Initially, we only allow one FSI interface
    if (interfaceToInterfaceList().size() != 1)
    {
        FatalErrorInFunction
            << "Currently, only one interface is allowed when using "
            << typeName << abort(FatalError);
    }

    // Allow only a direct map (conformal interface)
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // CAREFUL: we are setting non-zeros here based on the scalar rows, not
    // the block rows

    // Set matrix type to AIJ (since BAIJ does not support non-square
    // blocks)
    CHKERRQ(MatSetType(Ams, MATMPIAIJ));

    // Total number of scalar rows in the mesh motion region
    const label scalarRowN = fluidMesh.nCells()*motionBlockSize;

    // Allocate per-scalar-row nonzeros, initialised to 0
    std::vector<int> d_nnz(scalarRowN, 0);
    std::vector<int> o_nnz(scalarRowN, 0);

    // Set non-zeros for each interface fluid cells
    const labelList& fluidFaceMap = interfaceMap.zoneAToZoneBFaceMap();
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const labelList& fluidFaceCells =
        fluidMesh.boundary()[fluidPatchID].faceCells();
    forAll(fluidFaceMap, fluidFaceI)
    {
        const label fluidCellID = fluidFaceCells[fluidFaceI];

        // Calculate the row index for this cells first scalar equation
        label rowI = fluidCellID*motionBlockSize;

        // Set the number of non-zeros to be the number of solid equations
        d_nnz[rowI++] = solidBlockSize;
        d_nnz[rowI++] = solidBlockSize;

        if (!twoD)
        {
            d_nnz[rowI++] = solidBlockSize;
        }
    }

    // Parallel: we need to decide how to deal with this, e.g. do we allow
    // the same general decompositions as in the partitioned approach
    // For now, only allow serial
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Currently, serial runs are allowed in "
            << typeName << abort(FatalError);
    }

    // Allocate parallel matrix using AIJ
    CHKERRQ
    (
        MatMPIAIJSetPreallocation(Ams, 0, d_nnz.data(), 0, o_nnz.data())
    );

    // Raise an error if mallocs are required during matrix assembly
    CHKERRQ(MatSetOption(Ams, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));

    return 0;
}


label newtonCouplingInterface::initialiseAsf
(
    Mat Asf,
    const fvMesh& solidMesh,
    const label fluidBlockSize,
    const label solidBlockSize,
    const bool twoD
) const
{
    if (debug)
    {
        InfoInFunction
            << "Initialising" << endl;
    }

    // Allow only a direct map (conformal interface)
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // CAREFUL: we are setting non-zeros her based on the scalar rows, not
    // the block rows

    // Set matrix type to AIJ (since BAIJ does not support non-square
    // blocks)
    CHKERRQ(MatSetType(Asf, MATMPIAIJ));

    // Total number of scalar rows in the solid region
    const label scalarRowN = solidMesh.nCells()*solidBlockSize;

    // Allocate per-scalar-row nonzeros, initialised to 0
    std::vector<int> d_nnz(scalarRowN, 0);
    std::vector<int> o_nnz(scalarRowN, 0);

    // Set non-zeros for each interface solid cells
    const labelList& solidFaceMap = interfaceMap.zoneAToZoneBFaceMap();
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const labelList& solidFaceCells =
        solidMesh.boundary()[solidPatchID].faceCells();
    forAll(solidFaceMap, solidFaceI)
    {
        const label solidCellID = solidFaceCells[solidFaceI];

        // Calculate the row index for this cells first scalar equation
        label rowI = solidCellID*solidBlockSize;

        // Set the number of non-zeros to be the number of fluid equations
        // e.g., The x-momentum equation could have a coefficient for the
        // fluid x/y/z velocity and fluid pressure
        d_nnz[rowI++] = fluidBlockSize;
        d_nnz[rowI++] = fluidBlockSize;

        if (!twoD)
        {
            d_nnz[rowI++] = fluidBlockSize;
        }
    }

    // Parallel: we need to decide how to deal with this, e.g. do we allow
    // the same general decompositions as in the partitioned approach
    // For now, only allow serial
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Currently, serial runs are allowed in "
            << typeName << abort(FatalError);
    }

    // Allocate parallel matrix using AIJ
    CHKERRQ
    (
        MatMPIAIJSetPreallocation(Asf, 0, d_nnz.data(), 0, o_nnz.data())
    );

    // Raise an error if mallocs are required during matrix assembly
    CHKERRQ(MatSetOption(Asf, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));

    return 0;
}


label newtonCouplingInterface::formAfm
(
    Mat Afm,
    const label fluidBlockSize,
    const label motionBlockSize,
    const bool twoD
)
{
    // Next, we will manually insert to mesh-in-fluid coupling
    // This represents how the fluid equations depend on the mesh motion
    // The coupling in the momentum equation comes from the convection term
    // where the mesh flux appears. The coupling in the continuity equation
    // appears in the divergence of mesh motion

    // Take references
    const scalarField& phiI = fluid().phi();
    const vectorField& UI = fluid().U();
    const vectorField& SfI = fluidMesh().Sf();
    const scalarField& wI = fluidMesh().weights();
    const labelList& own = fluidMesh().owner();
    const labelList& nei = fluidMesh().neighbour();

    const globalIndex& globalCells =
        refCast<foamPetscSnesHelper>(fluid()).globalCells();

    // We insert the momentum coupling first
    // For the momentum, we will assume the momentum is 1st order upwind
    // discretised
    // The term is
    //      -div(phi, U) == -div(Sf & (U - meshU)*U)
    // Discretising the meshU term (minuses cancel):
    //     div(Sf & meshU*U)
    //       = sum_faces Sf*meshUf*Uf
    // where meshUf = w*meshUown + (1 - w)*meshUnei
    //
    // Assuming an upwind discretisation,
    // if (Sf & Uf) > 0:
    //     Uf = Uown
    // So, the term for the face becomes
    //     Sf*meshUf*Uown
    //       = Sf*[w*meshUown + (1 - w)*meshUnei]*Uown
    // The coefficient contributions become
    //     for the mesh own cell (differentiating wrt meshUown):
    //         coeff += Sf*w*Uown
    //     for the mesh nei cell (differentiating wrt meshUnei):
    //         coeff += Sf*(1 - w)*Uown
    // We must also add contributions considering face f from the
    // perspective of the fluid nei cell (i.e. Sf in opposite direction)
    //
    // Similarly
    // if (Sf & Uf) <= 0:
    //     Uf = Unei
    // So, the term for the face becomes
    //     Sf*meshUf*Unei
    //       = Sf*[w*meshUown + (1 - w)*meshUnei]*Unei
    // The coefficient contributions become
    //     for the mesh own cell (differentiating wrt meshUown):
    //         coeff += Sf*w*Unei
    //     for the mesh nei cell (differentiating wrt meshUnei):
    //         coeff += Sf*(1 - w)*Unei
    // We must also add contributions considering face f from the
    // perspective of the fluid nei cell (i.e. Sf in opposite direction)

    tensor coeff = tensor::zero;
    label rowI = -1;
    label colI = -1;
    forAll(own, faceI)
    {
        // Local row ID
        const label ownCellID = own[faceI];

        // Local column ID
        const label neiCellID = nei[faceI];

        // Global block row ID
        const label globalBlockRowI =
            globalCells.toGlobal(ownCellID);

        // Global block column ID
        const label globalBlockColI =
            globalCells.toGlobal(neiCellID);

        // Note: we cannot insert the coefficients as blocks since they are
        // not square (a requirement of MatSetValuesBlocked). Instead we
        // insert scalar coefficients (using MatSetValues)
        if (phiI[faceI] > 0)
        {
            // Add Sf*w*Uown to (fluid own, motion own)
            coeff = SfI[faceI]*wI[faceI]*UI[ownCellID];
            rowI = globalBlockRowI*fluidBlockSize;
            colI = globalBlockRowI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add Sf*(1 - w)*Uown to (fluid own, motion nei)
            coeff = SfI[faceI]*(1.0 - wI[faceI])*UI[ownCellID];
            rowI = globalBlockRowI*fluidBlockSize;
            colI = globalBlockColI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add -Sf*w*Uown to (fluid nei, motion nei)
            coeff = -SfI[faceI]*wI[faceI]*UI[ownCellID];
            rowI = globalBlockColI*fluidBlockSize;
            colI = globalBlockColI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add -Sf*(1 - w)*Uown to (fluid nei, motion own)
            coeff = -SfI[faceI]*(1.0 - wI[faceI])*UI[ownCellID];
            rowI = globalBlockColI*fluidBlockSize;
            colI = globalBlockRowI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }
        }
        else
        {
            // Add Sf*w*Unei to (fluid own, motion own)
            coeff = SfI[faceI]*wI[faceI]*UI[neiCellID];
            rowI = globalBlockRowI*fluidBlockSize;
            colI = globalBlockRowI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add Sf*(1 - w)*Unei to (fluid own, motion nei)
            coeff = SfI[faceI]*(1.0 - wI[faceI])*UI[neiCellID];
            rowI = globalBlockRowI*fluidBlockSize;
            colI = globalBlockColI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add -Sf*w*Unei to (fluid nei, motion nei)
            coeff = -SfI[faceI]*wI[faceI]*UI[neiCellID];
            rowI = globalBlockColI*fluidBlockSize;
            colI = globalBlockColI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add -Sf*(1 - w)*Unei to (fluid nei, motion own)
            coeff = -SfI[faceI]*(1.0 - wI[faceI])*UI[neiCellID];
            rowI = globalBlockColI*fluidBlockSize;
            colI = globalBlockRowI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }
        }
    }

    // Next we insert the continuity coupling
    // -div(phi - meshPhi) + pressure terms
    // Hence, the coefficients are (d/dD)(-div(meshPhi)), where D are the
    // cell centred mesh motion (displacements)
    // div(meshPhi) = sum_faces Sf & (w*Down + (1 - w)*Dnei)
    // So coefficients contributions per face are:
    //     own cell: Sf*w
    //     nei cell: Sf*(1 - w)
    vector coeff2 = vector::zero;
    forAll(own, faceI)
    {
        // Local row ID
        const label ownCellID = own[faceI];

        // Local column ID
        const label neiCellID = nei[faceI];

        // Global block row ID
        const label globalBlockRowI = globalCells.toGlobal(ownCellID);

        // Global block column ID
        const label globalBlockColI = globalCells.toGlobal(neiCellID);
        Warning
            << "CHECK WE ARE USING CORRECT GLOBAL CELLS" << endl;

        // The scalar row index for the pressure equation
        const label rowI =
            globalBlockRowI*fluidBlockSize + (fluidBlockSize - 1);

        // The scalar column index for the first component of the mesh
        // motion equation
        label colI = globalBlockRowI*motionBlockSize;

        // Add w*Sf to owner eqn
        coeff2 = wI[faceI]*SfI[faceI];
        for (int cmptI = 0; cmptI < motionBlockSize; ++cmptI)
        {
            CHKERRQ
            (
                MatSetValues
                (
                    Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                )
            );

            colI++;
        }

        // Scalar column index for neighbour
        colI = globalBlockColI*motionBlockSize;

        // Add (1 - w)*Sf to owner eqn
        coeff2 = (1.0 - wI[faceI])*SfI[faceI];
        for (int cmptI = 0; cmptI < motionBlockSize; ++cmptI)
        {
            CHKERRQ
            (
                MatSetValues
                (
                    Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                )
            );

            colI++;
        }
    }

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "div(meshPhi) coeffs for parallel not yet implemented"
            << exit(FatalError);
    }

    return 0;
}


label newtonCouplingInterface::formAfs
(
    Mat Afs,
    const label fluidBlockSize,
    const label solidBlockSize,
    const bool twoD
)
{
    if (debug)
    {
        InfoInFunction
            << "Start" << endl;
    }

    // The fluid interface is a prescribed velocity (fixedValue) condition
    // where we assume the fluid wall velocity is equal to the velocity of
    // the adjacent solid cell centre. This approximatation is sufficiently
    // accurate as a preconditioner for the matrix and will not affected the
    // converged solution (which is entirely governed by formResidual)

    if (fluidSolidInterface::fluidPatchIndices().size() != 1)
    {
        FatalError
            << "Only one interface patch is currently allowed"
            << abort(FatalError);
    }

    // Lookup the interface map from the fluid faces to the solid faces
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // The face map gives the solid face ID for each fluid face
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

    // Lookup the fluid interface patch
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const fvPatch& fluidPatch = fluidMesh().boundary()[fluidPatchID];
    const labelList& fluidFaceCells = fluidPatch.faceCells();

    // Lookup the solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const fvPatch& solidPatch = solidMesh().boundary()[solidPatchID];
    const labelList& solidFaceCells = solidPatch.faceCells();

    // Lookup the fluid patch
    const fvPatchVectorField& fluidPatchU =
        fluid().U().boundaryField()[fluidPatchID];
    if (!isA<fixedValueFvPatchVectorField>(fluidPatchU))
    {
        FatalErrorInFunction
            << "The fluid interface patch must be of type 'fixedValue'"
            << abort(FatalError);
    }

    if (!isA<fluidModels::newtonIcoFluid>(fluid()))
    {
        FatalErrorInFunction
            << "Currently, the fluid model must be of type 'newtonIcoFluid'"
            << abort(FatalError);
    }

    // Patch viscosity
    const scalarField fluidPatchNuEff
    (
        refCast<fluidModels::newtonIcoFluid>
        (
            fluid()
        ).turbulence().nuEff(fluidPatchID)
    );

    // Fluid interface area vectors
    const vectorField& fluidPatchSf = fluidPatch.Sf();

    // First we will insert the contribution to the fluid momentum equation
    // coming from the diffusion term
    // The known fluid boundary face value is now replaced by the adjacent
    // solid cell velocity
    // Diffusion coefficient is nu*|Sf|/(|n & d|*dt), where dt comes
    // from converting the solid displacement to velocity
    const scalar deltaT = solid().time().deltaTValue();
    const scalarField fluidPatchCoeffs
    (
        fluidPatch.magSf()*fluidPatch.deltaCoeffs()*fluidPatchNuEff/deltaT
    );

    // Second we will insert the contribution to the fluid continuity
    // (pressure) equation, where the div(U) term should use the adjacent
    // solid cell velocity instead of the known boundary face velocity

    forAll(fluidPatch, fluidFaceI)
    {
        const label solidFaceID = fluidFaceMap[fluidFaceI];

        // Fluid and solid cells, which are coupled
        const label fluidCellID = fluidFaceCells[fluidFaceI];
        const label solidCellID = solidFaceCells[solidFaceID];

        // We will add coefficients at block row "fluidCellID" and block
        // column "solidCellID"
        // Note that the nested monolithic matrix takes care of offsetting
        // the rows/columns when the submatrices are inserted into the
        // nested matrix

        // Global block row ID of fluid matrix
        const label globalBlockRowI =
            refCast<foamPetscSnesHelper>(fluid()).globalCells().toGlobal
            (
                fluidCellID
            );

        // Global block column ID of solid matrix
        const label globalBlockColI =
            refCast<foamPetscSnesHelper>(solid()).globalCells().toGlobal
            (
                solidCellID
            );

        // CAREFUL: we cannot use MatSetValuesBlocked as it only works with
        // square block coefficients, so we will insert the scalar
        // coefficients manually

        // Calculate the scalar global row ID (not the block row ID)
        label globalRowI = globalBlockRowI*fluidBlockSize;
        label globalColI = globalBlockColI*solidBlockSize;

        // Momentum coefficient for this face
        PetscScalar value = fluidPatchCoeffs[fluidFaceI];

        // Manually insert the 3 scalar coefficients (2 in 2-D)
        CHKERRQ
        (
            MatSetValues
            (
                Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        globalRowI++;
        globalColI++;
        CHKERRQ
        (
            MatSetValues
            (
                Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        if (!twoD)
        {
            globalRowI++;
            globalColI++;
            CHKERRQ
            (
                MatSetValues
                (
                    Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );
        }

        // Secondly we will insert the contributions for the pressure
        // equation

        // Manually insert the 3 scalar coefficients (2 in 2-D)
        value = -fluidPatchSf[fluidFaceI][vector::X]/deltaT;
        globalRowI++; // pressure equation
        globalColI = globalBlockColI*solidBlockSize; // x solid displacement
        CHKERRQ
        (
            MatSetValues
            (
                Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        value = -fluidPatchSf[fluidFaceI][vector::Y]/deltaT;
        globalColI++; // y solid displacement
        CHKERRQ
        (
            MatSetValues
            (
                Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        if (!twoD)
        {
            value = -fluidPatchSf[fluidFaceI][vector::Z]/deltaT;
            globalColI++; // y solid displacement
            CHKERRQ
            (
                MatSetValues
                (
                    Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );
        }
    }

    return 0;
}


label newtonCouplingInterface::formAms
(
    Mat Ams,
    const label solidBlockSize,
    const label motionBlockSize,
    const bool twoD
)
{
    // The mesh motion interface is a prescribed displacement (fixedValue)
    // condition where we assume the mesh motion interface displacement is
    // equal to the displacement of the adjacent solid cell centre. This
    // approximatation is sufficiently accurate as a preconditioner for the
    // matrix and will not affected the converged solution (which is entirely
    // governed by formResidual)

    if (fluidSolidInterface::fluidPatchIndices().size() != 1)
    {
        FatalErrorInFunction
            << "Only one interface patch is currently allowed"
            << abort(FatalError);
    }

    // Note: the fluid mesh is the same as the motion mesh, but we will use the
    // notation "motion" below rather than fluid to avoid confusion and
    // emphasise that we are referring to the mesh motion equations

    // Lookup the interface map from the motion faces to the solid faces
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // The face map gives the solid face ID for each motion face
    const labelList& motionFaceMap = interfaceMap.zoneBToZoneAFaceMap();

    // Lookup the motion interface patch
    const label motionPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const fvPatch& motionPatch = fluidMesh().boundary()[motionPatchID];
    const labelList& motionFaceCells = motionPatch.faceCells();

    // Lookup the solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const fvPatch& solidPatch = solidMesh().boundary()[solidPatchID];
    const labelList& solidFaceCells = solidPatch.faceCells();

    // Lookup the motion patch
    const fvPatchVectorField& motionPatchD =
        motionSolid().D().boundaryField()[motionPatchID];
    if (!isA<fixedValueFvPatchVectorField>(motionPatchD))
    {
        FatalErrorInFunction
            << "The motion interface patch must be of type 'fixedValue'"
            << abort(FatalError);
    }

    // Patch impK
    const scalarField motionPatchImpK
    (
        motionSolid().mechanical().impK()().boundaryField()[motionPatchID]
    );

    // For the motion momentum equation, the known interface displacement is now
    // replaced by the adjacent solid cell displacement
    // Diffusion coefficient is impK*|Sf|/(|n & d|)
    const scalarField motionPatchCoeffs
    (
        motionPatch.magSf()*motionPatch.deltaCoeffs()*motionPatchImpK
    );

    forAll(motionPatch, motionFaceI)
    {
        const label solidFaceID = motionFaceMap[motionFaceI];

        // Motion and solid cells, which are coupled
        const label motionCellID = motionFaceCells[motionFaceI];
        const label solidCellID = solidFaceCells[solidFaceID];

        // We will add coefficients at block row "motionCellID" and block
        // column "solidCellID"
        // Note that the nested monolithic matrix takes care of offsetting
        // the rows/columns when the submatrices are inserted into the
        // nested matrix

        // Global block row ID of motion matrix
        const label globalBlockRowI =
            refCast<foamPetscSnesHelper>(motionSolid()).globalCells().toGlobal
            (
                motionCellID
            );

        // Global block column ID of solid matrix
        const label globalBlockColI =
            refCast<foamPetscSnesHelper>(solid()).globalCells().toGlobal
            (
                solidCellID
            );

        // CAREFUL: we cannot use MatSetValuesBlocked as it only works with
        // square block coefficients, so we will insert the scalar
        // coefficients manually

        // Calculate the scalar global row ID (not the block row ID)
        label globalRowI = globalBlockRowI*motionBlockSize;
        label globalColI = globalBlockColI*solidBlockSize;

        // Momentum coefficient for this face
        PetscScalar value = motionPatchCoeffs[motionFaceI];

        // Manually insert the 3 scalar coefficients (2 in 2-D)
        CHKERRQ
        (
            MatSetValues
            (
                Ams, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        globalRowI++;
        globalColI++;
        CHKERRQ
        (
            MatSetValues
            (
                Ams, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        if (!twoD)
        {
            globalRowI++;
            globalColI++;
            CHKERRQ
            (
                MatSetValues
                (
                    Ams, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );
        }
    }

    return 0;
}


label newtonCouplingInterface::formAsf
(
    Mat Asf,
    const label fluidBlockSize,
    const label solidBlockSize,
    const bool twoD
)
{
    if (debug)
    {
        Info<< "Forming Asf" << endl;
    }

    // The solid interface is a prescribed traction condition, where we
    // approximate the traction on the fluid side of the interface using a
    // compact stencil. For this approximate Jacobian, we assume the traction
    // on a fluid interface face is equal to the pressure at the centre of the
    // adjacent fluid cell. This approximatation is sufficiently
    // accurate as a preconditioner for the matrix and will not affected the
    // converged solution (which is entirely governed by formResidual)

    if (fluidSolidInterface::solidPatchIndices().size() != 1)
    {
        FatalError
            << "Only one interface patch is currently allowed"
            << abort(FatalError);
    }

    // Lookup the interface map from the solid faces to the fluid faces
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // The face map gives the fluid face ID for each solid face
    const labelList& solidFaceMap = interfaceMap.zoneAToZoneBFaceMap();

    // Lookup the fluid interface patch
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const fvPatch& fluidPatch = fluidMesh().boundary()[fluidPatchID];
    const labelList& fluidFaceCells = fluidPatch.faceCells();

    // Lookup the solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const fvPatch& solidPatch = solidMesh().boundary()[solidPatchID];
    const labelList& solidFaceCells = solidPatch.faceCells();

    // The approximate force on the solid interface is the solid face area
    // vector multiplied by the adjacent fluid cell centre pressure. We also
    // need to multiply by the fluid density, if the kinematic pressure is
    // used. So the coefficient is the solid face area vector times the
    // fluid density
    // To-do: determine on the fly whether kinematic or dynamic pressure is
    // used
    const fvPatchVectorField& solidPatchD =
        solid().D().boundaryField()[solidPatchID];
    if (!isA<solidTractionFvPatchVectorField>(solidPatchD))
    {
        FatalErrorInFunction
            << "The solid interface patch must be of type 'solidTraction'"
            << abort(FatalError);
    }
    // Todo: add rho() to the fluidModel base class
    const vectorField patchCoeffs
    (
        solidPatch.Sf()
       *refCast<fluidModels::newtonIcoFluid>(fluid()).rho().value()
    );

    forAll(solidPatch, solidFaceI)
    {
        const label fluidFaceID = solidFaceMap[solidFaceI];

        // Fluid and solid cells, which are coupled
        const label solidCellID = solidFaceCells[solidFaceI];
        const label fluidCellID = fluidFaceCells[fluidFaceID];

        // We will add a coefficient at block row "solidCellID" and block
        // column "fluidCellID"
        // Note that the nested monolithic matrix takes care of offsetting
        // the rows/columns when the submatrices are inserted into the
        // nested matrix

        // Global block row ID of solid matrix
        const label globalBlockRowI =
            refCast<foamPetscSnesHelper>(solid()).globalCells().toGlobal
            (
                solidCellID
            );

        // Global block column ID of fluid matrix
        const label globalBlockColI =
            refCast<foamPetscSnesHelper>(fluid()).globalCells().toGlobal
            (
                fluidCellID
            );

        // CAREFUL: we cannot use MatSetValuesBlocked as it only works with
        // square block coefficients, so we will insert the scalar
        // coefficients manually

        // Calculate the scalar global row ID (not the block row ID)
        // The column corresponds to the pressure equation in the fluid cell
        label globalRowI = globalBlockRowI*solidBlockSize;
        const label globalColI =
            globalBlockColI*fluidBlockSize + (fluidBlockSize - 1);

        // Manually insert the 3 scalar coefficients (2 in 2-D)
        PetscScalar value = -patchCoeffs[solidFaceI][vector::X];
        CHKERRQ
        (
            MatSetValues
            (
                Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        globalRowI++;
        //globalColI++;
        value = -patchCoeffs[solidFaceI][vector::Y];
        CHKERRQ
        (
            MatSetValues
            (
                Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        if (!twoD)
        {
            globalRowI++;
            //globalColI++;
            value = -patchCoeffs[solidFaceI][vector::Z];
            CHKERRQ
            (
                MatSetValues
                (
                    Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );
        }
    }

    return 0;
}


void newtonCouplingInterface::updateMeshMotionInterfaceDisplacement()
{
    // Map solid interface motion to the mesh motion interface

    // Lookup the interface map from the fluid faces to the solid faces
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // The face map gives the solid face ID for each fluid face
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

    // Lookup the fluid mesh interface patch
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];

    // Lookup the mesh motion displacement field
    fvPatchVectorField& motionPatchD =
        motionSolid().D().boundaryFieldRef()[fluidPatchID];
    if (!isA<fixedValueFvPatchVectorField>(motionPatchD))
    {
        FatalErrorInFunction
            << "The meshMotionFluid interface patch must be of type "
            << "'fixedValue'" << abort(FatalError);
    }

    // Lookup the solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const fvPatchVectorField& solidPatchD =
        solid().D().boundaryField()[solidPatchID];

    // Map the solid interface displacement to the motion interface
    forAll(motionPatchD, fluidFaceI)
    {
        const label solidFaceID = fluidFaceMap[fluidFaceI];

        // Extrapolated patch value (larger stencil)
        motionPatchD[fluidFaceI] = solidPatchD[solidFaceID];
    }

    if (interfaceToInterfaceList().size() != 1)
    {
        FatalErrorInFunction
            << "Only one interface allowed" << abort(FatalError);
    }

    // Take references to zones
    const standAlonePatch& fluidZone =
        fluid().globalPatches()[0].globalPatch();
    const standAlonePatch& solidZone =
        solid().globalPatches()[0].globalPatch();

    // Get the solid interface patch pointD field
    const vectorField solidPatchPointD
    (
        solid().pointD().boundaryField()[solidPatchID].patchInternalField()
    );

    // Solid point zone pointD
    const vectorField solidZonePointD
    (
        solid().globalPatches()[0].patchPointToGlobal(solidPatchPointD)
    );

    // Prepare the fluid mesh interface pointD field
    vectorField fluidZonePointD(fluidZone.nPoints(), vector::zero);

    // Transfer the point displacements
    interfaceToInterfaceList()[0].transferPointsZoneToZone
    (
        solidZone,
        fluidZone,
        solidZonePointD,
        fluidZonePointD
    );

    // Map the zone field to the patch
    const vectorField meshPatchPointD
    (
        fluid().globalPatches()[0].globalPointToPatch(fluidZonePointD)
    );

    // Set the mesh interface pointD
    // Use "==" to reassign fixedValue
    motionSolid().pointD().boundaryFieldRef()[fluidPatchID] ==
        meshPatchPointD;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

newtonCouplingInterface::newtonCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    foamPetscSnesHelper
    (
        fileName
        (
            fsiProperties().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        2*fluid().mesh().nCells() + solid().mesh().nCells(),
        fsiProperties().lookupOrDefault<Switch>("stopOnPetscError", true),
        true // Will PETSc be used
    ),
    motionSystemScaleFactor_
    (
        fsiProperties().lookupOrDefault<scalar>
        (
            "motionSystemScaleFactor",
            gAverage
            (
                1.0/motionSolid().mechanical().shearModulus()().primitiveField()
            )
        )
    ),
    solidSystemScaleFactor_
    (
        fsiProperties().lookupOrDefault<scalar>
        (
            "solidSystemScaleFactor",
            gAverage(1.0/solid().mechanical().shearModulus()().primitiveField())
        )
    ),
    fluidToSolidCoupling_(fsiProperties().lookup("fluidToSolidCoupling")),
    solidToFluidCoupling_(fsiProperties().lookup("solidToFluidCoupling")),
    meshToFluidCoupling_(fsiProperties().lookup("meshToFluidCoupling")),
    solidToMeshCoupling_(fsiProperties().lookup("solidToMeshCoupling")),
    nRegions_(3),
    subMatsPtr_(nullptr)
{
    if (!isA<dynamicMotionSolverFvMesh>(fluidMesh()))
    {
        FatalErrorInFunction
            << "The fluid dynamic mesh must be dynamicMotionSolverFvMesh"
            << exit(FatalError);
    }

    if
    (
        !isA<meshMotionSolidModelFvMotionSolver>
        (
            refCast<dynamicMotionSolverFvMesh>(fluidMesh()).motion()
        )
    )
    {
        FatalErrorInFunction
            << "The fluid mesh motion solver must be "
            << meshMotionSolidModelFvMotionSolver::typeName
            << exit(FatalError);
    }

    // Not required as we manually move the mesh
    // // Tell the motion solver not to solve its own equations and the motion
    // // will be calculated here
    // refCast<const meshMotionSolidModelFvMotionSolver>
    // (
    //     refCast<dynamicMotionSolverFvMesh>(fluidMesh()).motion()
    // ).solveEquations() = false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool newtonCouplingInterface::evolve()
{
    // Steps
    // 1. Optional: predict solution field using old-time fields
    // 2. Map foam solution fields to PETSc
    // 3. Solve PETSc system
    // 4. Map PETSc solution to foam fields
    // 5. Update an secondary foam fields

    // Take references
    volVectorField& U = fluid().U();
    volScalarField& p = fluid().p();
    surfaceScalarField& phi = fluid().phi();
    volVectorField& D = solid().D();

    const label fluidBlockSize = fluid().twoD() ? 3 : 4;
    const label solidBlockSize = solid().twoD() ? 2 : 3;

    // WarningInFunction
    //     << "Todo: Add linear predictor" << endl;

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Access the raw solution data
    const PetscScalar *xx;
    VecGetArrayRead(foamPetscSnesHelper::solution(), &xx);

    // Retrieve the solution
    // Map the PETSc solution to the U field
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        //foamPetscSnesHelper::solution(),
        xx,
        U.primitiveFieldRef(),
        0, // Location of U
        fluidBlockSize,
        fluid().twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    U.correctBoundaryConditions();

    // Update the flux
    // CHECK
    // Make relative?
    phi = fvc::interpolate(U) & fluidMesh().Sf();

    // Map the PETSc solution to the p field
    // p is located in the final component
    foamPetscSnesHelper::ExtractFieldComponents<scalar>
    (
        //foamPetscSnesHelper::solution(),
        xx,
        p.primitiveFieldRef(),
        fluidBlockSize - 1, // Location of p component
        fluidBlockSize
    );

    p.correctBoundaryConditions();

    // Extract the displacement, which starts at row fluid.nCells*fluidBlockSize
    const label solidFirstEqnID = fluidMesh().nCells()*fluidBlockSize;
    ExtractFieldComponents
    (
        &xx[solidFirstEqnID],
        D.primitiveFieldRef(),
        0, // Location of first component
        solidBlockSize,
        solid().twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    // Restore the x vector
    VecRestoreArrayRead(foamPetscSnesHelper::solution(), &xx);

    // Interpolate cell displacements to vertices
    solid().mechanical().interpolate(solid().D(), solid().gradD(), solid().pointD());

    // Increment of displacement
    solid().DD() = solid().D() - solid().D().oldTime();

    // Increment of point displacement
    solid().pointDD() = solid().pointD() - solid().pointD().oldTime();

    // Velocity
    solid().U() = fvc::ddt(solid().D());

    // Update the mesh motion interface
    Info<< nl << "Solving mesh motion" << endl;
    updateMeshMotionInterfaceDisplacement();
    // meshMotion().evolve();
    Warning
        << "FIX EVOLVE: extract mesh motion" << endl;

    return 0;
}


label newtonCouplingInterface::initialiseJacobian(Mat& jac)
{
    // Create a nest matrix.

    // We will initialise sub-matrices:
    //
    //     *-----------------*
    //     | Aff | Afm | Afs |
    //     |-----------------|
    //     |  0  | Amm | Ams |
    //     |-----------------|
    //     | Asf |  0  | Ass |
    //     *-----------------*
    //
    // where the diagonal submatrices are
    //  - Aff: fluid equations (momentum and pressure)
    //  - Amm: mesh motion equations in the fluid domain
    //  - Ass: solid equations
    //
    // and the off diagonal submatrices, which represent coupling between
    // equations, are
    //  - Afm: mesh motion terms in the fluid equations
    //  - Afs: solid terms in the fluid equations
    //  - Ams: solid terms (interface motion) in the mesh motion equation
    //  - Asf: fluid terms (pressure on interface) in solid equations
    //
    // Note
    //  - Amf is empty as the the mesh motion does not directly
    //    depend on the fluid solution
    //  - Asm is empty as the solid equations do not depend on the fluid
    //    mesh motion

    if (nRegions_ != 3)
    {
        FatalErrorInFunction
            << "nRegions must be 3!" << abort(FatalError);
    }

    // Set twoD flag
    if (solid().twoD() != fluid().twoD())
    {
        FatalErrorInFunction
            << "Either the solid and fluid are both 2-D or both not 2-D"
            << exit(FatalError);
    }
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;
    const label motionBlockSize = solidBlockSize;

    // For brevity and convenience, we will store the size and blockSize of the
    // regions in a labelPairList
    const labelPairList nBlocksAndBlockSize
    (
        {
            {fluidMesh().nCells(), fluidBlockSize},
            {fluidMesh().nCells(), motionBlockSize},
            {solidMesh().nCells(), solidBlockSize}
        }
    );

    // Create arrays (vectors) of IS objects for rows and columns.
    // These will indicate where in the matrix the different regions are
    // located
    std::vector<IS> isRow(nRegions_), isCol(nRegions_);
    label globalOffset = 0;
    for (label r = 0; r < nRegions_; ++r)
    {
        const label nBlocks = nBlocksAndBlockSize[r].first();
        const label blockSize = nBlocksAndBlockSize[r].second();
        const label regionSize = nBlocks*blockSize;

        // Create an IS that covers the indices for region r
        ISCreateStride
        (
            PETSC_COMM_WORLD, regionSize, globalOffset, 1, &isRow[r]
        );
        ISCreateStride
        (
            PETSC_COMM_WORLD, regionSize, globalOffset, 1, &isCol[r]
        );
        globalOffset += regionSize;
    }

    // Create an array of submatrices.
    // Each submatrix represents the coupling between region i and
    // region j.
    // For example, subMat[i*nRegions_ + i] might be the matrix for
    // region i, while subMat[i*nRegions_ + j] (i != j) are the coupling
    // matrices.
    subMatsPtr_ = new Mat[nRegions_*nRegions_];
    for (label i = 0; i < nRegions_; ++i)
    {
        for (label j = 0; j < nRegions_; ++j)
        {
            // Create the submatrix for regions i and j.
            Mat subA = nullptr;

            // Skip Amf and Asm
            if
            (
                (i == 1 && j == 0) // Amf
             || (i == 2 && j == 1) // Asm                
            )
            {
                subMatsPtr_[i*nRegions_ + j] = subA;
                continue;
            }

            MatCreate(PETSC_COMM_WORLD, &subA);

            // Number of rows
            const label nRowBlocks = nBlocksAndBlockSize[i].first();
            const label rowBlockSize = nBlocksAndBlockSize[i].second();
            const label nRowsLocal = nRowBlocks*rowBlockSize;
            const label nRowsGlobal =
                returnReduce(nRowsLocal, sumOp<label>());

            // Number of columns
            const label nColBlocks = nBlocksAndBlockSize[j].first();
            const label colBlockSize = nBlocksAndBlockSize[j].second();
            const label nColsLocal = nColBlocks*colBlockSize;
            const label nColsGlobal =
                returnReduce(nColsLocal, sumOp<label>());

            if (debug)
            {
                Info<< "subMat(" << i << "," << j << ") " << nl
                    << "nRowsLocal = " << nRowsLocal << nl
                    << "nColsLocal = " << nColsLocal << endl;
            }

            MatSetSizes
            (
                subA,
                nRowsLocal,
                nColsLocal,
                nRowsGlobal,
                nColsGlobal
            );
            MatSetFromOptions(subA);
            MatSetType(subA, MATMPIAIJ);

            // Store subA in row-major order at index [i*nRegions_ + j]
            subMatsPtr_[i*nRegions_ + j] = subA;
        }
    }

    // Create the nest matrix using the index sets and submatrices.
    jac = Mat();
    MatCreateNest
    (
        PETSC_COMM_WORLD,
        nRegions_,
        isRow.data(),
        nRegions_,
        isCol.data(),
        (Mat*)subMatsPtr_,
        &jac
    );

    // Cleanup IS objects.
    for (label r = 0; r < nRegions_; ++r)
    {
        ISDestroy(&isRow[r]);
        ISDestroy(&isCol[r]);
    }


    // Get access to the sub-matrices
    PetscInt nr, nc;
    Mat **subMats;
    CHKERRQ(MatNestGetSubMats(jac, &nr, &nc, &subMats));

    // Initialise the diagonal submatrices:
    //  - Aff: fluid equations (momentum and pressure)
    //  - Amm: mesh motion equations in the fluid domain
    //  - Ass: solid equations

    // Aff
    Foam::initialiseJacobian
    (
        subMats[0][0], fluid().mesh(), fluidBlockSize, false
    );

    // Amm
    Foam::initialiseJacobian
    (
        subMats[1][1], fluid().mesh(), motionBlockSize, false
    );

    // Ass
    Foam::initialiseJacobian
    (
        subMats[2][2], solid().mesh(), solidBlockSize, false
    );

    // Initialise the four off diagonal submatrices:
    //  - Afm: mesh motion terms in the fluid equations
    //  - Afs: solid terms in the fluid equations
    //  - Ams: solid terms (interface motion) in the mesh motion equation
    //  - Asf: fluid terms (pressure on interface) in solid equations

    // Afm
    initialiseAfm
    (
        subMats[0][1], fluidMesh(), fluidBlockSize, motionBlockSize, twoD
    );

    // Afs
    initialiseAfs
    (
        subMats[0][2], fluidMesh(), fluidBlockSize, solidBlockSize, twoD
    );

    // Ams
    initialiseAms
    (
        subMats[1][2], fluidMesh(), motionBlockSize, solidBlockSize, twoD
    );

    // Asf
    initialiseAsf
    (
        subMats[2][0], solidMesh(), fluidBlockSize, solidBlockSize, twoD
    );

    return 0;
}


label newtonCouplingInterface::initialiseSolution(Vec& x)
{
    // Set twoD flag
    if (solid().twoD() != fluid().twoD())
    {
        FatalErrorInFunction
            << "Either the solid and fluid are both 2-D or both not 2-D"
            << exit(FatalError);
    }
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;
    const label motionBlockSize = solidBlockSize;

    // For brevity and convenience, we will store the size and blockSize of the
    // regions in a labelPairList
    const labelPairList nBlocksAndBlockSize
    (
        {
            {fluidMesh().nCells(), fluidBlockSize},
            {fluidMesh().nCells(), motionBlockSize},
            {solidMesh().nCells(), solidBlockSize}
        }
    );

    // Count number of local blocks and local scalar equations
    label n = 0;
    forAll(nBlocksAndBlockSize, regionI)
    {
        const label nBlocksRegionI = nBlocksAndBlockSize[regionI].first();
        const label blockSizeRegionI = nBlocksAndBlockSize[regionI].second();
        n += nBlocksRegionI*blockSizeRegionI;
    }

    // Global system size: total number of scalar equation across all
    // processors
    const label N = returnReduce(n, sumOp<label>());

    // Create solution vector
    x = Vec();
    CHKERRQ(VecCreate(PETSC_COMM_WORLD, &x));
    CHKERRQ(VecSetSizes(x, n, N));
    CHKERRQ(VecSetType(x, VECMPI));
    CHKERRQ(PetscObjectSetName((PetscObject) x, "Solution"));
    CHKERRQ(VecSetFromOptions(x));
    CHKERRQ(VecZeroEntries(x));

    return 0;
}


label newtonCouplingInterface::formResidual
(
    PetscScalar *f,
    const PetscScalar *x
)
{
    // Considerations on the order of updating the fluid, solid, and mesh motion
    // residuals
    //  - solid depends on the fluid interface traction
    //  - mesh motion depends on the solid interface displacement
    //  - fluid depends on the solid interface velocity
    //  - fluid depends on the entire mesh motion flux field
    //
    // Approach
    // 1. Update the fluid velocity and pressure fields and calculate the
    //    traction at the fluid interface
    // 2. Map the fluid interface traction to the solid interface
    // 3. Update the solid residual, which now has the correct interface
    //    traction
    // 4. Map the solid interface velocity to the fluid interface
    // 5. Map the solid interface displacement to the mesh motion interface
    // 6. Update the mesh motion residual, which now has the correct interface
    //    displacement
    // 7. Move the fluid mesh using the mesh motion field
    // 8. Update the fluid residual, which now has the correct interface
    //    velocity and mesh motion
    // 9. Apply scaling factor to solid and motion equations to preserve the
    //    condition number of the monolithic system

    // Set twoD flag
    if (solid().twoD() != fluid().twoD())
    {
        FatalErrorInFunction
            << "Either the solid and fluid are both 2-D or both not 2-D"
            << exit(FatalError);
    }
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;
    const label motionBlockSize = solidBlockSize;

    // The scalar row at which the motion equations start
    const label motionFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

    // The scalar row at which the solid equations start
    const label solidFirstEqnID =
        motionFirstEqnID + fluidMesh().nCells()*motionBlockSize;

    // Currently limited to one interface: it should be straight-forward to add
    // a loop over multiple interface => todo
    if (fluidSolidInterface::fluidPatchIndices().size() != 1)
    {
        FatalErrorInFunction
            << "Only one interface patch is currently allowed"
            << abort(FatalError);
    }

    // 1. Update the fluid velocity and pressure fields and calculate the
    //    traction at the fluid interface
    {
        // Take references
        volVectorField& U = fluid().U();
        volScalarField& p = fluid().p();

        // Retrieve the solution
        // Map the PETSc solution to the U field
        foamPetscSnesHelper::ExtractFieldComponents<vector>
        (
            x,
            U.primitiveFieldRef(),
            0, // Location of U
            fluidBlockSize,
            fluid().twoD() ? labelList({0,1}) : labelList({0,1,2})
        );

        U.correctBoundaryConditions();

        const_cast<volTensorField&>
        (
            fluidMesh().lookupObject<volTensorField>("grad(U)")
        ) = fvc::grad(U);

        U.correctBoundaryConditions();

        // Map the PETSc solution to the p field
        // p is located in the final component
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x,
            p.primitiveFieldRef(),
            fluidBlockSize - 1, // Location of p component
            fluidBlockSize
        );

        p.correctBoundaryConditions();
    }

    // 2. Map the fluid interface traction to the solid interface
    if (fluidToSolidCoupling_)
    {
        // Fluid interface traction
        const label fluidPatchID =
            fluidSolidInterface::fluidPatchIndices()[0];
        const vectorField fluidNf
        (
            fluidMesh().boundary()[fluidPatchID].nf()
        );
        const vectorField fluidTraction
        (
            fluid().patchViscousForce(fluidPatchID)
          - fluid().patchPressureForce(fluidPatchID)*fluidNf
        );

        // Lookup the interface map from the fluid faces to the solid faces
        const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping& interfaceMap =
            refCast
            <
                const interfaceToInterfaceMappings::
                directMapInterfaceToInterfaceMapping
            >
            (
                interfaceToInterfaceList()[0]
            );

        // The face map gives the solid face ID for each fluid face
        const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

        // Calculate the solid interface traction
        // Flip the sign as the solid normals point in the opposite direction to
        // the fluid normals
        const label solidPatchID =
            fluidSolidInterface::solidPatchIndices()[0];
        vectorField solidTraction
        (
            solidMesh().boundary()[solidPatchID].size()
        );
        forAll(fluidTraction, fluidFaceI)
        {
            solidTraction[fluidFaceMap[fluidFaceI]] =
                -fluidTraction[fluidFaceI];
        }

        // Lookup the displacement interface traction patch and set the traction
        fvPatchVectorField& solidPatchD =
            solid().D().boundaryFieldRef()[solidPatchID];
        if (!isA<solidTractionFvPatchVectorField>(solidPatchD))
        {
            FatalErrorInFunction
                << "The solidinterface patch must be of type "
                << "'solidTraction'"
                << abort(FatalError);
        }

        solidTractionFvPatchVectorField& solidTractionPatch =
            refCast<solidTractionFvPatchVectorField>(solidPatchD);

        solidTractionPatch.traction() = solidTraction;
    }


    // 3. Update the solid residual, which now has the correct interface
    //    traction
    refCast<foamPetscSnesHelper>(solid()).formResidual
    (
        &f[solidFirstEqnID], &x[solidFirstEqnID]
    );

    // 4. Map the solid interface velocity to the fluid interface
    // Note: it is assumed the solid.U() is dD/dt even for steady-state cases
    if (solidToFluidCoupling_)
    {
        // Lookup the interface map from the fluid faces to the solid faces
        const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping& interfaceMap =
            refCast
            <
                const interfaceToInterfaceMappings::
                directMapInterfaceToInterfaceMapping
            >
            (
                interfaceToInterfaceList()[0]
            );

        // The face map gives the solid face ID for each fluid face
        const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

        // Lookup the fluid interface patch
        const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
        fvPatchVectorField& fluidPatchU =
            fluid().U().boundaryFieldRef()[fluidPatchID];
        if (!isA<fixedValueFvPatchVectorField>(fluidPatchU))
        {
            FatalErrorInFunction
                << "The fluid interface patch must be of type 'fixedValue'"
                << abort(FatalError);
        }

        // Lookup the solid interface velocity and displacement
        const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
        const fvPatchVectorField& solidPatchU =
            solid().U().boundaryField()[solidPatchID];

        // Map the solid interface velocity to the fluid interface and the solid
        // displacement to the motion interface
        // const labelList& solidFaceCells =
        //     solidMesh().boundary()[solidPatchID].faceCells();
        // const vectorField& solidUI = solid().U();
        forAll(fluidPatchU, fluidFaceI)
        {
            const label solidFaceID = fluidFaceMap[fluidFaceI];

            // Extrapolated patch value (larger stencil)
            fluidPatchU[fluidFaceI] = solidPatchU[solidFaceID];

            // Adjacent cell value
            // const label solidCellID = solidFaceCells[solidFaceID];
            // fluidPatchU[fluidFaceI] = solidUI[solidCellID];
        }

        fluid().phi() = fvc::interpolate(fluid().U()) & fluidMesh().Sf();
    }

    // 5. Map the solid interface displacement to the mesh motion interface
    if (solidToMeshCoupling_)
    {
        // Lookup the interface map from the fluid faces to the solid faces
        const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping& interfaceMap =
            refCast
            <
                const interfaceToInterfaceMappings::
                directMapInterfaceToInterfaceMapping
            >
            (
                interfaceToInterfaceList()[0]
            );

        // The face map gives the solid face ID for each motion face
        const labelList& motionFaceMap = interfaceMap.zoneBToZoneAFaceMap();

        // Lookup the motion interface patch
        const label motionPatchID = fluidSolidInterface::fluidPatchIndices()[0];
        fvPatchVectorField& motionPatchD =
            motionSolid().D().boundaryFieldRef()[motionPatchID];
        if (!isA<fixedValueFvPatchVectorField>(motionPatchD))
        {
            FatalErrorInFunction
                << "The motion interface patch must be of type 'fixedValue'"
                << abort(FatalError);
        }

        // Lookup the solid interface velocity and displacement
        const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
        const fvPatchVectorField& solidPatchD =
            solid().D().boundaryField()[solidPatchID];

        // Map the solid interface displacement to the motion interface
        // const labelList& solidFaceCells =
        //     solidMesh().boundary()[solidPatchID].faceCells();
        // const vectorField& solidDI = solid().D();
        forAll(motionPatchD, motionFaceI)
        {
            const label solidFaceID = motionFaceMap[motionFaceI];

            // Extrapolated patch value (larger stencil)
            motionPatchD[motionFaceI] = solidPatchD[solidFaceID];

            // Adjacent cell value
            // const label solidCellID = solidFaceCells[solidFaceID];
            // motionPatchD[motionFaceI] = solidDI[solidCellID];
        }
    }

    // 6. Update the mesh motion residual, which now has the correct interface
    //    displacement
    refCast<foamPetscSnesHelper>(motionSolid()).formResidual
    (
        &f[motionFirstEqnID], &x[motionFirstEqnID]
    );

    // 7. Move the fluid mesh using the mesh motion field
    if (meshToFluidCoupling_)
    {
        fluidMesh().movePoints(motionSolid().pointD());
        fluid().U().correctBoundaryConditions();
    }

    // 8. Update the fluid residual, which now has the correct interface
    //    velocity and mesh motion
    // Note that the fluid equations are first in the f (residual) and x
    // (solution) lists
    refCast<foamPetscSnesHelper>(fluid()).formResidual(f, x);

    // 9. Apply scaling factor to solid and motion equations to preserve the
    // condition number of the monolithic system

    const label solidSystemEnd =
        solidFirstEqnID + solidMesh().nCells()*solidBlockSize;
    for (int i = solidFirstEqnID; i < solidSystemEnd; ++i)
    {
        f[i] *= solidSystemScaleFactor_;
    }
    
    const label motionSystemEnd =
        motionFirstEqnID + fluidMesh().nCells()*motionBlockSize;
    for (int i = motionFirstEqnID; i < motionSystemEnd; ++i)
    {
        f[i] *= motionSystemScaleFactor_;
    }

    return 0;
}


label newtonCouplingInterface::formJacobian
(
    Mat jac,
    const PetscScalar *x
)
{
    // We will assembly the seven sub-matrices:

    // Diagonal submatrices:
    //  - Aff: fluid equations (momentum and pressure)
    //  - Amm: mesh motion equations in the fluid domain
    //  - Ass: solid equations
    //
    // Off-diagonal submatrices:
    //  - Afm: mesh motion terms in the fluid equations
    //  - Afs: solid terms in the fluid equations
    //  - Ams: solid terms (interface motion) in the mesh motion equation
    //  - Asf: fluid terms (pressure on interface) in solid equations

    // Zero entries
    MatZeroEntries(jac);

    // Get access to the two sub-matrices
    PetscInt nr, nc;
    Mat **subMats;
    CHKERRQ(MatNestGetSubMats(jac, &nr, &nc, &subMats));

    if (nr != 3 || nc != 3)
    {
        FatalErrorInFunction
            << "The matrix has the wrong number of sub matrices: "
            << "nr = " << nr << ", nc = " << nc << abort(FatalError);
    }

    // Set twoD flag
    if (solid().twoD() != fluid().twoD())
    {
        FatalErrorInFunction
            << "Either the solid and fluid are both 2-D or both not 2-D"
            << exit(FatalError);
    }
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;
    const label motionBlockSize = solidBlockSize;

    // The scalar row at which the motion equations start
    const label motionFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

    // The scalar row at which the solid equations start
    const label solidFirstEqnID =
        motionFirstEqnID + fluidMesh().nCells()*motionBlockSize;

    // Check the solid and fluid derive from foamPetscSnesHelper
    if
    (
        !isA<foamPetscSnesHelper>(fluid()) || !isA<foamPetscSnesHelper>(solid())
    )
    {
        FatalErrorInFunction
            << "You must use solid and fluid models derived from the "
            << "foamPetscSnesHelper class" << exit(FatalError);
    }

    // Form diagonal submatrices
    //  - Aff: fluid equations (momentum and pressure)
    //  - Amm: mesh motion equations in the fluid domain
    //  - Ass: solid equations
    //

    // Aff
    refCast<foamPetscSnesHelper>(fluid()).formJacobian
    (
        subMats[0][0], x
    );

    // Amm
    refCast<foamPetscSnesHelper>(motionSolid()).formJacobian
    (
        subMats[1][1], &x[motionFirstEqnID]
    );

    // Ass
    refCast<foamPetscSnesHelper>(solid()).formJacobian
    (
        subMats[2][2], &x[solidFirstEqnID]
    );

    // Form off-diagonal submatrices:
    //  - Afm: mesh motion terms in the fluid equations
    //  - Afs: solid terms in the fluid equations
    //  - Ams: solid terms (interface motion) in the mesh motion equation
    //  - Asf: fluid terms (pressure on interface) in solid equations

    // Afm
    if (meshToFluidCoupling_)
    {
        formAfm(subMats[0][1], fluidBlockSize, motionBlockSize, twoD);
    }

    // Afs
    if (solidToFluidCoupling_)
    {
        formAfs(subMats[0][2], fluidBlockSize, solidBlockSize, twoD);
    }

    // Ams
    if (solidToMeshCoupling_)
    {
        formAms(subMats[1][2], solidBlockSize, motionBlockSize, twoD);
    }

    // Asf
    if (fluidToSolidCoupling_)
    {
        formAsf(subMats[2][0], fluidBlockSize, solidBlockSize, twoD);
    }


    // Scale the solid matrix to preserve the condition number of the
    // monolithic system
    // We must assembly the matrix before we can use MatScale
    // Complete matrix assembly
    CHKERRQ(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));

    // Amm
    MatScale(subMats[1][1], motionSystemScaleFactor_);

    // Ams
    MatScale(subMats[1][2], motionSystemScaleFactor_);

    // Asf
    MatScale(subMats[2][0], solidSystemScaleFactor_);

    // Ass
    MatScale(subMats[2][2], solidSystemScaleFactor_);

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
