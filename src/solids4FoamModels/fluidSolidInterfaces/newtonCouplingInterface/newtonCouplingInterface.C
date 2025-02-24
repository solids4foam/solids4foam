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
        labelPairList
        (
            {
                {
                    fluid().mesh().nCells(),
                    fluid().twoD() ? 3 : 4 // Number of unknowns per block
                },
                {
                    solid().mesh().nCells(),
                    solid().twoD() ? 2 : 3 // Number of unknowns per block
                }
            }
        ),
        labelListList(),
        fsiProperties().lookupOrDefault<Switch>("stopOnPetscError", true),
        true // Will PETSc be used
    )
{}


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
    const label solidBlockSize = solid().twoD() ? 2 : 4;

    WarningInFunction
        << "Todo: map U/p/D to x" << endl;

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

    return 0;
}


label newtonCouplingInterface::initialiseJacobian(Mat jac)
{
    // Get access to the two sub-matrices
    PetscInt nr, nc;
    Mat **subMats;
    CHKERRQ(MatNestGetSubMats(jac, &nr, &nc, &subMats));

    if (nr != 2 || nc != 2)
    {
        FatalErrorInFunction
            << "The matrix has the wrong number of sub matrices: "
            << "nr = " << nr << ", nc = " << nc << abort(FatalError);
    }

    // Set fluid and solid block sizes
    const label fluidBlockSize = fluid().twoD() ? 3 : 4;
    const label solidBlockSize = solid().twoD() ? 2 : 4;

    // We will initialise the four sub-matrices:
    //  - Afluid: fluid
    //  - Asolid: solid
    //  - Asif: solid-in-fluid-coupling
    //  - Afis: fluid-in-solid coupling

    // Initialise Afluid based on compact stencil fvMesh
    {
        Mat Afluid = subMats[0][0];
        Foam::initialiseJacobian(Afluid, fluid().mesh(), fluidBlockSize);
    }

    // Initialise Asolid based on compact stencil fvMesh
    {
        Mat Asolid = subMats[1][1];
        Foam::initialiseJacobian(Asolid, solid().mesh(), solidBlockSize);
    }

    // Initially we assume a conformal FSI interface, where each fluid cell
    // shares a face with a solid cell. So we assume the number of blocks in
    // the Asif (and Afis) matrix is equal to the number of cells at the
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

    // Initialise solid-in-fluid-coupling matrix
    {
        Mat Asif = subMats[0][1];

        // Set the block size: this corresponds to the row size for a
        // rectangular matrix
        CHKERRQ(MatSetBlockSizes(jac, fluidBlockSize, solidBlockSize));

        // Number of blocks: number of faces at the interface
        const label blockn = interfaceMap.zoneBToZoneAFaceMap().size();

        // Number of on-processor non-zeros blocks per row
        int* d_nnz = (int*)malloc(blockn*sizeof(int));

        // Number of off-processor non-zeros per row
        int* o_nnz = (int*)malloc(blockn*sizeof(int));

        // Initialise d_nnz to 1 and o_nnz to zero
        for (int i = 0; i < blockn; ++i)
        {
            d_nnz[i] = 1;
            o_nnz[i] = 0;
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
        // // Count off-processor neighbour cells: todo
        // forAll(mesh.boundary(), patchI)
        // {
        //     if (mesh.boundary()[patchI].type() == "processor")
        //     {
        //         // Todo
        //     }
        // }

        // Allocate parallel matrix
        // Use row block size (assumes square matrix)
        CHKERRQ
        (
            MatMPIBAIJSetPreallocation(Asif, fluidBlockSize, 0, d_nnz, 0, o_nnz)
        );

        // Raise an error if mallocs are required during matrix assembly
        MatSetOption(Asif, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    }

    // Initialise fluid-in-solid-coupling matrix
    {
        Mat Afis = subMats[1][0];

        // Set the block size: this corresponds to the row size for a
        // rectangular matrix
        CHKERRQ(MatSetBlockSizes(jac, solidBlockSize, fluidBlockSize));

        // Number of blocks: number of faces at the interface
        const label blockn = interfaceMap.zoneAToZoneBFaceMap().size();

        // Number of on-processor non-zeros blocks per row
        int* d_nnz = (int*)malloc(blockn*sizeof(int));

        // Number of off-processor non-zeros per row
        int* o_nnz = (int*)malloc(blockn*sizeof(int));

        // Initialise d_nnz to 1 and o_nnz to zero
        for (int i = 0; i < blockn; ++i)
        {
            d_nnz[i] = 1;
            o_nnz[i] = 0;
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
        // // Count off-processor neighbour cells: todo
        // forAll(mesh.boundary(), patchI)
        // {
        //     if (mesh.boundary()[patchI].type() == "processor")
        //     {
        //         // Todo
        //     }
        // }

        // Allocate parallel matrix
        // Use row block size (assumes square matrix)
        CHKERRQ
        (
            MatMPIBAIJSetPreallocation(Afis, solidBlockSize, 0, d_nnz, 0, o_nnz)
        );

        // Raise an error if mallocs are required during matrix assembly
        MatSetOption(Afis, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    }

    return 0;
}


label newtonCouplingInterface::formResidual
(
    PetscScalar *f,
    const PetscScalar *x
)
{
    // Approach
    // 1. Update the fluid velocity and pressure fields and calculate the
    //     traction at the fluid interface
    // 2. Map the fluid interface traction to the solid interface
    // 3. Update the solid residual, which now has the correct interface
    //    traction
    // 4. Map the solid interface velocity to the fluid interface
    // 5. Update the fluid residual, which now has the correct interface
    //    velocity

    // Block sizes
    const label fluidBlockSize = fluid().twoD() ? 3 : 4;
    // const label solidBlockSize = solid().twoD() ? 2 : 3;

    // Currently limited to one interface: it should be straight-forward to add
    // a loop over multiple interface => todo
    if (fluidSolidInterface::fluidPatchIndices().size() != 1)
    {
        FatalError
            << "Only one interface patch is currently allowed"
            << abort(FatalError);
    }

    // 1. Update the fluid velocity and pressure fields and calculate the
    //     traction at the fluid interface
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

        // Fluid interface traction
        const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
        const vectorField fluidNf(fluidMesh().boundary()[fluidPatchID].nf());
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
        const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
        vectorField solidTraction(solidMesh().boundary()[solidPatchID].size());
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
                << "The solidinterface patch must be of type 'solidTraction'"
                << abort(FatalError);
        }

        solidTractionFvPatchVectorField& solidTractionPatch =
            refCast<solidTractionFvPatchVectorField>(solidPatchD);

        solidTractionPatch.traction() = solidTraction;
    }


    // Calculate the solid residual first
    // The first solid scalar equation is at row fluidMesh.nCells*fluidBlockSize
    const label solidFirstEqnID = fluidMesh().nCells()*fluidBlockSize;
    refCast<foamPetscSnesHelper>(solid()).formResidual
    (
        &f[solidFirstEqnID], &x[solidFirstEqnID]
    );

    // We update the velocity at the fluid interface with the solid interface
    // velocity
    // Note: it is assumed the solid.U is dD/dt even for steady-state cases
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

        // Lookup the solid interface patch
        const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
        const fvPatchVectorField& solidPatchU =
            solid().U().boundaryField()[solidPatchID];

        // Map the solid interface velocity to the fluid interface
        forAll(fluidPatchU, fluidFaceI)
        {
            fluidPatchU[fluidFaceI] = solidPatchU[fluidFaceMap[fluidFaceI]];
        }
    }

    // Then, calculate the fluid residual
    // Note that the fluid equations are first in the f (residual) and x
    // (solution) lists
    refCast<foamPetscSnesHelper>(fluid()).formResidual(f, x);

    return 0;
}


label newtonCouplingInterface::formJacobian
(
    Mat jac,
    const PetscScalar *x
)
{
    // We will assembly the four sub-matrices:
    //  - Afluid: fluid
    //  - Asolid: solid
    //  - Asif: solid-in-fluid-coupling
    //  - Afis: fluid-in-solid coupling

    // Get access to the two sub-matrices
    PetscInt nr, nc;
    Mat **subMats;
    CHKERRQ(MatNestGetSubMats(jac, &nr, &nc, &subMats));

    if (nr != 2 || nc != 2)
    {
        FatalErrorInFunction
            << "The matrix has the wrong number of sub matrices: "
            << "nr = " << nr << ", nc = " << nc << abort(FatalError);
    }

    // Assemble the fluid matrix
    {
        // Mat Afluid;
        // MatNestGetSubMat(jac, 0, 0, Afluid);
        Mat Afluid = subMats[0][0];

        if (!isA<foamPetscSnesHelper>(fluid()))
        {
            FatalErrorInFunction
                << "A fluid model derived from foamPetscSnesHelper must be used!"
                << abort(FatalError);
        }

        // The FSI interface is a prescribed velocity condition
        refCast<foamPetscSnesHelper>(fluid()).formJacobian(Afluid, x);
    }

    // Assemble the solid matrix
    {
        Mat Asolid = subMats[1][1];

        if (!isA<foamPetscSnesHelper>(solid()))
        {
            FatalErrorInFunction
                << "A solid model derived from foamPetscSnesHelper must be used!"
                << abort(FatalError);
        }

        // The FSI interface is a prescribed velocity condition
        refCast<foamPetscSnesHelper>(solid()).formJacobian(Asolid, x);
    }

    // Set fluid and solid block sizes
    const label fluidBlockSize = fluid().twoD() ? 3 : 4;
    const label solidBlockSize = solid().twoD() ? 2 : 4;

    // const label nFluidCells = fluidMesh().nCells();

    // Assemble the solid-in-fluid coupling matrix
    {
        Mat Asif = subMats[0][1];

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

        // Lookup the fluid patch boundary gradient coefficients
        // These are the coefficients to be multiplied by the prescribed patch
        // velocity
        const fvPatchVectorField& fluidPatchU =
            fluid().U().boundaryField()[fluidPatchID];
        if (!isA<fixedValueFvPatchVectorField>(fluidPatchU))
        {
            FatalErrorInFunction
                << "The fluid interface patch must be of type 'fixedValue'"
                << abort(FatalError);
        }
        //const vectorField fluidPatchCoeffs(fluidPatchU.gradientBoundaryCoeffs());
        if (!isA<fluidModels::newtonIcoFluid>(fluid()))
        {
            FatalErrorInFunction
                << "Currently, the fluid model must be of type 'newtonIcoFluid'"
                << abort(FatalError);
        }
        const scalarField fluidPatchNuEff
        (
            refCast<fluidModels::newtonIcoFluid>
            (
                fluid()
            ).turbulence().nuEff(fluidPatchID)
        );
        const vectorField fluidPatchCoeffs
        (
            fluidPatch.Sf()*fluidPatch.deltaCoeffs()*fluidPatchNuEff
        );

        // Initialise the block coefficient
        // const label nCoeffCmpts = fluidBlockSize*solidBlockSize;
        // PetscScalar values[nCoeffCmpts];
        // std::memset(values, 0, sizeof(values));

        forAll(fluidPatch, fluidFaceI)
        {
            const label solidFaceID = fluidFaceMap[fluidFaceI];

            // Fluid and solid cells, which are coupled
            const label fluidCellID = fluidFaceCells[fluidFaceI];
            const label solidCellID = solidFaceCells[solidFaceID];

            // We will add a coefficient at block row "fluidCellID" and block
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

            // Manually insert the 3 scalar coefficients (2 in 2-D)
            PetscScalar value = fluidPatchCoeffs[fluidFaceI][vector::X];
            CHKERRQ
            (
                MatSetValues
                (
                    Asif, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );

            globalRowI++;
            globalColI++;
            value = fluidPatchCoeffs[fluidFaceI][vector::Y];
            CHKERRQ
            (
                MatSetValues
                (
                    Asif, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );

            if (fluidBlockSize > 3)
            {
                globalRowI++;
                globalColI++;
                value = fluidPatchCoeffs[fluidFaceI][vector::Z];
                CHKERRQ
                (
                    MatSetValues
                    (
                        Asif, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                    )
                );
            }
        }
    }

    // Assemble the fluid-in-solid coupling matrix
    {
        Mat Afis = subMats[1][0];

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
        // vector multiplied by the adjacent fluid cell centre pressure. So the
        // coefficient is the solid face area vector
        const fvPatchVectorField& solidPatchD =
            solid().D().boundaryField()[solidPatchID];
        if (!isA<solidTractionFvPatchVectorField>(solidPatchD))
        {
            FatalErrorInFunction
                << "The solid interface patch must be of type 'solidTraction'"
                << abort(FatalError);
        }
        const vectorField& solidPatchSf = solidPatch.Sf();

        // Initialise the block coefficient
        // const label nCoeffCmpts = solidBlockSize*fluidBlockSize;
        // PetscScalar values[nCoeffCmpts];
        // std::memset(values, 0, sizeof(values));

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
            PetscScalar value = -solidPatchSf[solidFaceI][vector::X];
            CHKERRQ
            (
                MatSetValues
                (
                    Afis, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );

            globalRowI++;
            //globalColI++;
            value = -solidPatchSf[solidFaceI][vector::Y];
            CHKERRQ
            (
                MatSetValues
                (
                    Afis, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );

            if (solidBlockSize > 2)
            {
                globalRowI++;
                //globalColI++;
                value = -solidPatchSf[solidFaceI][vector::Z];
                CHKERRQ
                (
                    MatSetValues
                    (
                        Afis, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                    )
                );
            }
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
