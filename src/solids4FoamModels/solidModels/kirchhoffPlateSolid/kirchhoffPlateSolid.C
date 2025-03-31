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

#include "kirchhoffPlateSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "faCFD.H"
#include "linearElastic.H"
#include "BlockLduSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kirchhoffPlateSolid, 0);
addToRunTimeSelectionTable(solidModel, kirchhoffPlateSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool kirchhoffPlateSolid::converged
(
    const int iCorr,
    const lduSolverPerformance& solverPerfM,
    const lduSolverPerformance& solverPerfw,
    const areaScalarField& M,
    const areaScalarField& w
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate relative residuals
    const scalar residualM =
        gMax
        (
            mag(M.internalField() - M.prevIter().internalField())
           /max
            (
                gMax(mag(M.internalField() - M.oldTime().internalField())),
                SMALL
            )
        );

    const scalar residualw =
        gMax
        (
            mag(w.internalField() - w.prevIter().internalField())
           /max
            (
                gMax(mag(w.internalField() - w.oldTime().internalField())),
                SMALL
            )
        );

    // Calculate material residual
    const scalar materialResidual = mechanical().residual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at leaast 1 outer iteration and the material law must be converged
    if (iCorr > 1 && materialResidual < materialTol())
    {
        bool convergedM = false;
        bool convergedw = false;

        if
        (
            (
                solverPerfM.initialResidual() < solutionTol()
             && residualM < solutionTol()
            )
         || solverPerfM.initialResidual() < alternativeTol()
         || residualM < alternativeTol()
        )
        {
            convergedM = true;
        }

        if
        (
            (
                solverPerfw.initialResidual() < solutionTol()
             && residualw < solutionTol()
            )
         || solverPerfw.initialResidual() < alternativeTol()
         || residualw < alternativeTol()
        )
        {
            convergedw = true;
        }


        if (convergedM && convergedw)
        {
            Info<< "    The residuals have converged" << endl;
            converged = true;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res (M & w), relRes (M & w), matRes, iters (M & w)"
            << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfM.initialResidual()
            << ", " << solverPerfw.initialResidual()
            << ", " << residualM
            << ", " << residualw
            << ", " << materialResidual
            << ", " << solverPerfM.nIterations()
            << ", " << solverPerfw.nIterations() << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr() - 1)
    {
        maxIterReached()++;
        Warning
            << "Max iterations reached within the M-w loop" << endl;
    }

    return converged;
}


const fvPatch& kirchhoffPlateSolid::areaPatch() const
{
    if (areaPatchID_ == -1)
    {
        calcAreaPatches();
    }

    return mesh().boundary()[areaPatchID_];
}


const fvPatch& kirchhoffPlateSolid::areaShadowPatch() const
{
    if (areaShadowPatchID_ == -1)
    {
        calcAreaPatches();
    }

    return mesh().boundary()[areaShadowPatchID_];
}


void kirchhoffPlateSolid::calcAreaPatches() const
{
    // Note: face0PatchID may be -1 if this processor has no faces on the
    // finiteArea patch

    // Check that all areaMesh faces map to the same patch

    const polyMesh& pMesh = mesh();
    const polyBoundaryMesh& bm = pMesh.boundaryMesh();
    const labelList& faceLabels = aMesh_.faceLabels();
    const label pMeshNFaces = pMesh.nFaces();

    if (faceLabels.size() > 0)
    {
        const label face0ID = faceLabels[0];

        if (face0ID < pMeshNFaces)
        {
            areaPatchID_ = bm.whichPatch(face0ID);

            // Check all faces map to the same fvMesh patch

            forAll(faceLabels, aFaceI)
            {
                const label faceID = faceLabels[aFaceI];

                // Escape if face is beyond active faces, eg belongs to a face
                // zone
                if (faceID < pMeshNFaces)
                {
                    const label curPatchID = bm.whichPatch(face0ID);

                    if (curPatchID != areaPatchID_)
                    {
                        FatalErrorIn
                        (
                            "void kirchhoffPlateSolid::calcAreaPatches() const"
                        )   << "The finiteArea patch should correspond to a "
                            << "patch on the boundary of the polyMesh!"
                            << abort(FatalError);
                    }
                }
            }
        }
    }


    // We will now check if the polyMesh has the same number of cells as the
    // number of faces on the areaPatch, as we are assuming the polyMesh to be
    // one cell thick
    if (pMesh.nCells() != pMesh.boundaryMesh()[areaPatchID_].size())
    {
        FatalErrorIn
        (
            "void kirchhoffPlateSolid::calcAreaPatches() const"
        )   << "The solid polyMesh should be one cell thick, where there is "
            << "the same number of cells as the number of faces on the "
            << "areaPatch" << endl
            << "areaPatchID: " << areaPatchID_
            << abort(FatalError);
    }


    // To find the areaShadowPatch, we will ...
    const unallocLabelList& faceCells =
        pMesh.boundaryMesh()[areaPatchID_].faceCells();

    if (faceCells.size())
    {
        const cellList& cells = pMesh.cells();

        const label face0ID = bm[areaPatchID_].start();
        const label cell0ID = faceCells[0];
        const vector& face0N = bm[areaPatchID_].faceNormals()[0];
        const labelList& curCellFaces = cells[cell0ID];

        scalar mostNegativeDotProduct = GREAT;

        forAll(curCellFaces, fI)
        {
            const label curFaceID = curCellFaces[fI];

            if (curFaceID != face0ID)
            {
                if (!pMesh.isInternalFace(curFaceID))
                {
                    const label otherPatchID = bm.whichPatch(curFaceID);
                    const label curLocalFaceID =
                        curFaceID - bm[otherPatchID].start();

                    const vector& curFaceN =
                        bm[otherPatchID].faceNormals()[curLocalFaceID];

                    const scalar dotProduct = face0N & curFaceN;

                    if (dotProduct < mostNegativeDotProduct)
                    {
                        mostNegativeDotProduct = dotProduct;
                        areaShadowPatchID_ = otherPatchID;
                    }
                }
            }
        }
    }


    // Check if the areaPatch and areaShadowPatch have the same number of faces
    if
    (
        pMesh.boundaryMesh()[areaShadowPatchID_].size()
     != pMesh.boundaryMesh()[areaPatchID_].size()
    )
    {
        FatalErrorIn
        (
            "void kirchhoffPlateSolid::calcAreaPatches() const"
        )   << "The polyMesh should be one cell thick, where there should be "
            << "two patches opposite each other that have the same number of "
            << "faces" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kirchhoffPlateSolid::kirchhoffPlateSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    aMesh_(mesh()),
    w_
    (
        IOobject
        (
            "w",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    wVf_
    (
        IOobject
        (
            "wVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    M_
    (
        IOobject
        (
            "M",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    MVf_
    (
        IOobject
        (
            "MVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimPressure/dimArea, 0.0)
    ),
    p_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    pVf_
    (
        IOobject
        (
            "pVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", p_.dimensions(), 0.0)
    ),
    theta_
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_,
        dimensionedVector("zero", dimless, vector::zero)
    ),
    thetaVf_
    (
        IOobject
        (
            "thetaVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimless, vector::zero)
    ),
    gradTheta_(fac::grad(theta_)),
    rho_("zero", dimDensity, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    h_(solidModelDict().lookup("plateThickness")),
    bendingStiffness_("zero", dimPressure*dimVolume, 0.0),
    areaPatchID_(-1),
    areaShadowPatchID_(-1)
{
    const PtrList<mechanicalLaw>& mechLaws = mechanical();

    // Only the linearElastic mechanicalLaw is allow and one material
    if (mechLaws.size() != 1)
    {
        FatalErrorIn(type() + "::" + type())
            << " can currently only be used with a single material"
            << abort(FatalError);
    }
    else if (!isA<linearElastic>(mechLaws[0]))
    {
        FatalErrorIn(type() + "::" + type())
            << " can only be used with the linearElastic "
            << "mechanicalLaw" << nl
            << abort(FatalError);
    }

    // Cast the mechanical law to a linearElastic mechanicalLaw
    const linearElastic& mech = refCast<const linearElastic>(mechLaws[0]);

    // Set plate properties
    rho_ = mech.rhoScalar();
    E_ = mech.E();
    nu_ = mech.nu();
    bendingStiffness_ = E_*pow(h_, 3)/(12*(1 - pow(nu_, 2)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool kirchhoffPlateSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Create volume-to surface mapping object
    volSurfaceMapping vsm(aMesh_);

    // Mesh update loop
    do
    {
        int iCorr = 0;
        lduSolverPerformance solverPerfM;
        lduSolverPerformance solverPerfw;
        blockLduMatrix::debug = 0;

        Info<< "Solving the Kirchhoff plate equation for w and M" << endl;

        // w and M equation loop
        do
        {
            // Algorithm
            // Solve M equation
            // Solve w equation
            // where
            // M is the moment sum
            // w is the transvere (out of plane) displacement
            // The M equation is:
            //     rho*h*fac::d2dt2(w) = fam::laplacian(M) + p
            // and the w equation is:
            //     fam::laplacian(D, w) + M
            // where
            // rho is the density
            // h is the plate thickness
            // p is the net transverse pressure
            // D is the bending stiffness
            // THESE ARE BOTH SCALARS: IDEA: USE BLOCK COUPLED!

            // Store fields for under-relaxation and residual calculation
            M_.storePrevIter();

            // Solve M equation
            // d2dt2 not implemented so calculate it manually using ddt
            // Also, "==" complains so we will move all terms to left
            faScalarMatrix MEqn
            (
                rho_*h_
               *(
                   fac::ddt(w_) - fac::ddt(w_.oldTime())
                )/runTime().deltaT()
              - fam::laplacian(M_) - p_
            );

            // Relax the linear system
            MEqn.relax();

            // Solve the linear system
            solverPerfM = MEqn.solve();

            // Relax the field
            M_.relax();

            // Store fields for under-relaxation and residual calculation
            w_.storePrevIter();

            // Solve w equation
            faScalarMatrix wEqn
            (
                fam::laplacian(bendingStiffness_, w_) + M_
            );

            // Relax the linear system
            wEqn.relax();

            // Solve the linear system
            solverPerfw = wEqn.solve();

            // Relax the field
            w_.relax();

            // Update the angle of rotation
            theta_ = -fac::grad(w_);

            // Update the gradient of rotation field, used for non-orthogonal
            // correction in clamped boundary conditions
            gradTheta_ = fac::grad(theta_);

            // Work in progress
            // Use block-coupled framework to implicitly couple the equations
            //if (coupled_)
            // {
            //     // Prepare block system
            //     BlockLduSystem<vector2, vector2> blockM(aMesh_);

            //     // Grab block diagonal and set it to zero
            //     Field<tensor2>& d = blockM.diag().asSquare();
            //     d = tensor2::zero;

            //     // Grab linear off-diagonal and set it to zero
            //     Field<vector2>& l = blockM.lower().asLinear();
            //     Field<vector2>& u = blockM.upper().asLinear();
            //     u = vector2::zero;
            //     l = vector2::zero;

            //     // Insert MEqn and wEqn (without coupling terms)
            //     to-do

            //     // Insert coupling terms implicitly
            //     to-do

            //     // Solve the linear system
            //     to-do

            //     // Retrieve solution
            //     to-do

            //     // Correct the boundary conditions for M and w
            //     M.correctBoundaryConditions();
            //     w.correctBoundaryConditions();
            // }
        }
        while
        (
            !converged(iCorr, solverPerfM, solverPerfw, M_, w_)
         && ++iCorr < nCorr()
        );

        // Map area fields to vol fields
        mapAreaFieldToSingleLayerVolumeField(M_, MVf_);
        mapAreaFieldToSingleLayerVolumeField(w_, wVf_);
        mapAreaFieldToSingleLayerVolumeField(theta_, thetaVf_);
        mapAreaFieldToSingleLayerVolumeField(p_, pVf_);
        {
            const areaVectorField Ds = w_*aMesh_.faceAreaNormals();
            mapAreaFieldToSingleLayerVolumeField(Ds, D());
        }

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());

    // Disable writing of the stress fields, as it is not calculated
    sigma().writeOpt() = IOobject::NO_WRITE;

    return true;
}


tmp<vectorField> kirchhoffPlateSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    notImplemented(type() + "::tractionBoundarySnGrad(...)");

    // Keep compiler happy
    return tmp<vectorField>();
}


void kirchhoffPlateSolid::setTraction
(
    const label patchID,
    const vectorField& traction
)
{
    if (traction.size() != p_.size())
    {
        FatalErrorIn("void kirchhoffPlateSolid::setTraction(...)")
            << "Something is wrong with the length of the traction field passed"
            << "to the solid!"
            << abort(FatalError);
    }

    // Take normal component of the traction field
    // Note: p is the net pressure on the plate (from both sides)
    p_.internalField() = aMesh_.faceAreaNormals().internalField() & traction;
}


void kirchhoffPlateSolid::writeFields(const Time& runTime)
{
    // Do not call solidModel::writeFields() as we do not want to write the
    // stress and strain fields

    physicsModel::writeFields(runTime);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
