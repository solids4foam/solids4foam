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

#include "vertexCentredLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "SparseMatrixTemplate.H"
#include "vfvcCellPoint.H"
#include "vfvmCellPoint.H"
#include "fvcDiv.H"
#include "fixedValuePointPatchFields.H"
#include "zeroGradientPointPatchFields.H"
#include "sparseMatrixTools.H"
#include "symmetryPointPatchFields.H"
#ifdef USE_PETSC
    #include <petscksp.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(vertexCentredLaplacian, 0);
addToRunTimeSelectionTable(solidModel, vertexCentredLaplacian, dictionary);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void vertexCentredLaplacian::updateSource
(
    scalarField& source,
    const pointScalarField& pointT,
    const surfaceVectorField& dualGradTf,
    const dimensionedScalar& diffusivity,
    const labelList& dualCellToPoint
)
{
    if (debug)
    {
        Info<< "void vertexCentredLaplacian::updateSource(...): start"
            << endl;
    }

    // Reset to zero
    source = 0.0;

    // The source vector is -F, where:
    // F = div(diffusivity*grad(T))

    // Point volume field
    const scalarField& pointVolI = pointVol_.internalField();

    // Calculate the flux (n & gradT) at the dual faces
    surfaceScalarField dualFlux
    (
        (dualMesh().Sf()/dualMesh().magSf()) & dualGradTf*diffusivity
    );

    // Enforce flux on flux boundaries
    enforceGradientBoundaries
    (
        pointT, dualFlux, mesh(), dualMeshMap().pointToDualFaces()
    );

    // Set coupled boundary (e.g. processor) flux fields to zero: this
    // ensures their global contribution is zero
    forAll(dualFlux.boundaryField(), patchI)
    {
        if (dualFlux.boundaryField()[patchI].coupled())
        {
#ifdef OPENFOAM_NOT_EXTEND
            dualFlux.boundaryFieldRef()[patchI] = 0.0;
#else
            dualFlux.boundaryField()[patchI] = 0.0;
#endif
        }
    }

    // Calculate divergence of diffusivity*gradT for the dual cells
    const scalarField dualLaplacian = fvc::div(dualFlux*dualMesh().magSf());

    // Map dual cell field to primary mesh point field
    scalarField pointLaplacian(mesh().nPoints(), 0.0);
    forAll(dualLaplacian, dualCellI)
    {
        const label pointID = dualCellToPoint[dualCellI];
        pointLaplacian[pointID] = dualLaplacian[dualCellI];
    }

    // Add to the source
    source -= pointLaplacian*pointVolI;

    // Add source terms
    // Can be added

    // Add transient term
    // Can be added

    if (debug)
    {
        Info<< "void vertexCentredLaplacian::updateSource(...): end"
            << endl;
    }
}


void vertexCentredLaplacian::setFixedDofs
(
    const pointScalarField& pointT,
    boolList& fixedDofs,
    scalarField& fixedDofValues
) const
{
    // Flag all fixed DOFs
    forAll(pointT.boundaryField(), patchI)
    {
        if
        (
            isA<fixedValuePointPatchScalarField>
            (
                pointT.boundaryField()[patchI]
            )
        )
        {
            const labelList& meshPoints =
                pointT.mesh().mesh().boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];
                const scalar val = pointT[pointID];

                // Check if this point has already been fixed
                if (fixedDofs[pointID])
                {
                    // Check if the new value is consistent with the previous
                    // fixed value
                    if (mag(fixedDofValues[pointID] - val) > SMALL)
                    {
                        WarningIn
                        (
                            "void vertexCentredLaplacian::setFixedDofs(...)"
                        )   << "Point " << pointID << " "
                            << pointT.mesh().mesh().points()[pointID]
                            << " has a prescribed value from more than one "
                            << "patch: the value from the patch with the lowest"
                            << " index will be used." << endl;
                    }
                }
                else
                {
                    fixedDofs[pointID] = true;
                    fixedDofValues[pointID] = val;
                }
            }
        }
    }

    if (debug > 1)
    {
        Info<< "fixedDofs:" << endl;
        forAll(fixedDofs, pI)
        {
            if (fixedDofs[pI])
            {
                Info<< pI << ": " << fixedDofValues[pI] << endl;
            }
        }
    }
}


void vertexCentredLaplacian::enforceGradientBoundaries
(
    const pointScalarField& pointT,
    surfaceScalarField& dualFlux,
    const fvMesh& mesh,
    const labelListList& pointToDualFaces
) const
{
    forAll(pointT.boundaryField(), patchI)
    {
        // Optionally we could add a fixedGradient condition
        if
        (
            isA<zeroGradientPointPatchScalarField>
            (
                pointT.boundaryField()[patchI]
            )
         || isA<symmetryPointPatchScalarField>(pointT.boundaryField()[patchI])
        )
        {
            // Set the dual patch face flux to zero
#ifdef OPENFOAM_NOT_EXTEND
            dualFlux.boundaryFieldRef()[patchI] = 0.0;
#else
            dualFlux.boundaryField()[patchI] = 0.0;
#endif
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vertexCentredLaplacian::vertexCentredLaplacian
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    pointT_
    (
        IOobject
        (
            "pointT",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimTemperature, 0.0)
    ),
    diffusivity_(solidModelDict().lookup("diffusivity")),
    twoD_(sparseMatrixTools::checkTwoD(mesh())),
    fixedDofs_(mesh().nPoints(), false),
    fixedDofValues_(fixedDofs_.size(), 0.0),
    fixedDofScale_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "fixedDofScale",
            diffusivity_.value()*Foam::sqrt(gAverage(mesh().magSf()))
        )
    ),
    pointVol_
    (
        IOobject
        (
            "pointVolumes",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimVolume, 0.0)
    ),
    gradT_
    (
        IOobject
        (
            "grad(T)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", pointT_.dimensions()/dimLength, vector::zero) //,
        //"zeroGradient"
    ),
    dualGradTf_
    (
        IOobject
        (
            "grad(T)f",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedVector("zero", pointT_.dimensions()/dimLength, vector::zero),
        "calculated"
    ),
    globalPointIndices_(mesh())
#ifdef OPENFOAM_COM
    ,
    pointVolInterp_(pMesh(), mesh())
#endif
{
    // Create dual mesh and set write option
    dualMesh().objectRegistry::writeOpt() = IOobject::NO_WRITE;

    // Set fixed degree of freedom list
    setFixedDofs(pointT_, fixedDofs_, fixedDofValues_);

    // Set the pointVol field
    // Map dualMesh cell volumes to the primary mesh points
#ifdef OPENFOAM_NOT_EXTEND
    scalarField& pointVolI = pointVol_.primitiveFieldRef();
#else
    scalarField& pointVolI = pointVol_.internalField();
#endif
    const scalarField& dualCellVol = dualMesh().V();
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualCellToPoint, dualCellI)
    {
        // Find point which maps to this dual cell
        const label pointID = dualCellToPoint[dualCellI];

        // Map the cell volume
        pointVolI[pointID] = dualCellVol[dualCellI];
    }

    // Write fixed degree of freedom equation scale
    Info<< "fixedDofScale: " << fixedDofScale_ << endl;

    // Disable the writing of the unused fields
    D().writeOpt() = IOobject::NO_WRITE;
    D().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    DD().writeOpt() = IOobject::NO_WRITE;
    DD().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    U().writeOpt() = IOobject::NO_WRITE;
    pointDD().writeOpt() = IOobject::NO_WRITE;
}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

vertexCentredLaplacian::~vertexCentredLaplacian()
{
#ifdef USE_PETSC
    if (Switch(solidModelDict().lookup("usePETSc")))
    {
        PetscFinalize();
    }
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool vertexCentredLaplacian::evolve()
{
    Info<< "Evolving vertex-centred Laplacian solver" << endl;

    // Initialise matrix
    sparseScalarMatrix matrix(sum(globalPointIndices_.stencilSize()));

    // Global point index lists
    // const boolList& ownedByThisProc = globalPointIndices_.ownedByThisProc();
    // const labelList& localToGlobalPointMap =
    //     globalPointIndices_.localToGlobalPointMap();

    // Lookup compact edge gradient factor
    const scalar zeta(solidModelDict().lookupOrDefault<scalar>("zeta", 0.2));
    if (debug)
    {
        Info<< "zeta: " << zeta << endl;
    }

    // Assemble matrix once per time-step
    Info<< "    Assembling the matrix" << endl;

    // Add Laplacian coefficients
    vfvm::laplacian
    (
        matrix,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        scalarField(mesh().nCells(), diffusivity_.value()),
        debug
    );

    // Solution field: pointT correction
    scalarField pointTcorr(pointT_.internalField().size(), 0.0);

    // Calculate grad at dual faces
    dualGradTf_ = vfvc::fGrad
    (
        pointT_,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta,
        debug
    );

    // Update the source vector
    scalarField source(mesh().nPoints(), 0.0);
    pointT_.correctBoundaryConditions();
    updateSource
    (
        source, 
        pointT_,
        dualGradTf_,
        diffusivity_,
        dualMeshMap().dualCellToPoint()
    );

    if (debug > 1)
    {
        // Print the matrix
        matrix.print();
    }

    // Enforce fixed DOF on the linear system
    sparseMatrixTools::enforceFixedDof
    (
        matrix,
        source,
        fixedDofs_,
        fixedDofValues_,
        fixedDofValues_,
        fixedDofScale_,
        debug
    );

    if (debug > 1)
    {
        // Print the matrix
        matrix.print();
    }

    // Solve linear system for displacement correction
    if (debug)
    {
        Info<< "bool vertexCentredLaplacian::evolve(): "
            << " solving linear system: start" << endl;
    }
    else
    {
        Info<< "    Solving" << endl;
    }

    if (Switch(solidModelDict().lookup("usePETSc")))
    {
#ifdef USE_PETSC
        notImplemented
        (
            "sparseMatrixTools::solveLinearSystemPETSc(...scalar...) "
            "to be implemented!"
        );
        // fileName optionsFile(solidModelDict().lookup("optionsFile"));
        // solverPerf = sparseMatrixTools::solveLinearSystemPETSc
        // (
        //     matrix,
        //     source,
        //     pointDcorr,
        //     twoD_,
        //     optionsFile,
        //     mesh().points(),
        //     ownedByThisProc,
        //     localToGlobalPointMap,
        //     globalPointIndices_.stencilSizeOwned(),
        //     globalPointIndices_.stencilSizeNotOwned(),
        //     solidModelDict().lookupOrDefault<bool>("debugPETSc", false)
        // );
#else
        FatalErrorIn("vertexCentredLaplacian::evolve()")
            << "PETSc not available. Please set the PETSC_DIR environment "
            << "variable and re-compile solids4foam" << abort(FatalError);
#endif
    }
    else
    {
        // Lookup exportToMatlab flag
        const Switch writeMatlabMatrix
        (
            solidModelDict().lookup("writeMatlabMatrix")
        );

        // Use Eigen SparseLU direct solver
        sparseMatrixTools::solveLinearSystemEigen
        (
            matrix, source, pointTcorr, writeMatlabMatrix, debug
        );
    }

    if (debug)
    {
        Info<< "bool vertexCentredLaplacian::evolve(): "
            << " solving linear system: end" << endl;
    }


    // Update solution field
#ifdef OPENFOAM_NOT_EXTEND
    pointT_.primitiveFieldRef() += pointTcorr;
#else
    pointT_.internalField() += pointTcorr;
#endif

    // Calculate grad at dual faces
    dualGradTf_ = vfvc::fGrad
    (
        pointT_,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta,
        debug
    );

    // Calculate cell gradient
    // This assumes a constant gradient within each primary mesh cell
    // This is a first-order approximation
    gradT_ = vfvc::grad(pointT_, mesh());

    return true;
}


void vertexCentredLaplacian::writeFields(const Time& runTime)
{
    // Optionally we could also calculate the flux field in the cells for
    // visualisation

    solidModel::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
