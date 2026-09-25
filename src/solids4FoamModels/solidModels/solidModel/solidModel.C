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

#include "solidModel.H"
#include "mechanicalConstitutiveLawManager.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "symmetryPolyPatch.H"
#include "twoDPointCorrector.H"
#include "solidTractionFvPatchVectorField.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "primitivePatchInterpolation.H"
#else
    #include "blockSolidTractionFvPatchVectorField.H"
#endif
#include "fvcGradf.H"
#include "wedgePolyPatch.H"
#include "meshTools.H"
#include "addToRunTimeSelectionTable.H"
#include "compatibilityFunctions.H"
#include "hofvc.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "zeroGradientFvPatchFields.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "enhancedVolPointInterpolation.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidModel, 0);
    defineRunTimeSelectionTable(solidModel, dictionary);
    addToRunTimeSelectionTable(physicsModel, solidModel, physicsModel);

#ifdef OPENFOAM_COM
    const Enum<solidModel::solutionAlgorithm>
    solidModel::solutionAlgorithmNames_
    ({
        {
            solidModel::solutionAlgorithm::PETSC_SNES,
            "PETScSNES"
        },
        {
            solidModel::solutionAlgorithm::IMPLICIT_COUPLED,
            "implicitCoupled"
        },
        {
            solidModel::solutionAlgorithm::IMPLICIT_SEGREGATED,
            "implicitSegregated"
        },
        {
            solidModel::solutionAlgorithm::EXPLICIT,
            "explicit"
        },
    });
#else
    template<>
    const char* NamedEnum<solidModel::solutionAlgorithm, 4>::names[] =
    {
        "PETScSNES",
        "implicitCoupled",
        "implicitSegregated",
        "explicit"
    };
#endif

#ifdef OPENFOAM_ORG
    typedef meshFaceZones faceZoneMesh;
#endif
}


#ifndef OPENFOAM_COM
const Foam::NamedEnum<Foam::solidModel::solutionAlgorithm, 4>
    Foam::solidModel::solutionAlgorithmNames_;
#endif


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidModel::makeDualMesh() const
{
    if (dualMeshPtr_.valid())
    {
        FatalErrorIn("void Foam::solidModel::makeDualMesh() const")
            << "Pointer already set!" << abort(FatalError);
    }

    Info<< "Creating dualMesh" << endl;

    dualMeshPtr_.set(new meshDual(mesh(), solidModelDict()));
}


void Foam::solidModel::makeRAUf() const
{
    if (rAUfPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    rAUfPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "rAUf",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("0", dimPressure, 0.0)
        )
    );
}


const Foam::surfaceScalarField& Foam::solidModel::rAUf() const
{
    if (rAUfPtr_.empty())
    {
        makeRAUf();
    }

    return autoPtrRef(rAUfPtr_);
}


Foam::surfaceScalarField& Foam::solidModel::rAUf()
{
    if (rAUfPtr_.empty())
    {
        makeRAUf();
    }

    return autoPtrRef(rAUfPtr_);
}


void Foam::solidModel::checkWedges() const
{
    const fvMesh& mesh = this->mesh();

    label nWedgePatches = 0;
    vector wedgeDirVec = vector::zero;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>
            (
                mesh.boundaryMesh()[patchI]
            );

            nWedgePatches++;
            wedgeDirVec += cmptMag(wpp.centreNormal());

            // Make sure that solidWedge is used instead of wedge
            if
            (
                DD_.boundaryField()[patchI].type() == "wedge"
             && D_.boundaryField()[patchI].type() == "wedge"
            )
            {
                FatalErrorIn("void Foam::solidModel::checkWedges() const")
                    << "solidWedge should be used on displacement solution "
                    << "field wedge patches as non-orthogonal corrections "
                    << "are important!"
                    << abort(FatalError);
            }
        }
    }

    reduce(nWedgePatches, maxOp<label>());

    if (nWedgePatches)
    {
        if (nWedgePatches != 2)
        {
            FatalErrorIn("void Foam::solidModel::checkWedges() const")
                << "For axisymmetric cases, there should be exactly two wedge "
                << "patches!" << abort(FatalError);
        }

        Info<< nl << "Axisymmetric case: disabling the solution in the "
            << "out-of-plane direction" << endl;

        // We will use const_cast to disable the out-of-lane direction
        Vector<label>& solD = const_cast<Vector<label>&>(mesh.solutionD());

        reduce(wedgeDirVec, sumOp<vector>());

        wedgeDirVec /= mag(wedgeDirVec);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (wedgeDirVec[cmpt] > 1e-6)
            {
                solD[cmpt] = -1;

                wordList dirs(3);
                dirs[0] = "x";
                dirs[1] = "y";
                dirs[2] = "z";
                Info<< "    out-of-plane direction: " << dirs[cmpt] << nl
                    << endl;
            }
            else
            {
                solD[cmpt] = 1;
            }
        }
    }


    // Check all the face normals are in the same direction on the wedge patches
    // This is to avoid the case where a wedge patch is composed of two
    // disconnected regions with one on the front and one on the back
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            // Unit face normals on processor
            const vectorField nf = mesh.boundaryMesh()[patchI].faceNormals();

            if (nf.size() == 0)
            {
                FatalErrorIn("void Foam::solidModel::checkWedges() const")
                    << "There are no faces on the wedge patch "
                    << mesh.boundaryMesh()[patchI].name() << " on this processor:"
                    << nl << "Every processor should have at least one face on "
                    << "each wedge patch"
                    << abort(FatalError);
            }

            // Check that all the wedge face normals point in the same direction

            vector firstFaceNOnMasterProc = vector::zero;

            if (Pstream::master())
            {
                firstFaceNOnMasterProc = nf[0];
            }

            // Sync in parallel so that all processors have the master vector
            reduce(firstFaceNOnMasterProc, sumOp<vector>());

            forAll(nf, faceI)
            {
                if ((nf[faceI] & firstFaceNOnMasterProc) < 0)
                {
                    FatalErrorIn("void Foam::solidModel::checkWedges() const")
                        << "On wedge patch "
                        << mesh.boundaryMesh()[patchI].name()
                        << " there are at "
                        << "least two faces with unit normals in the opposite "
                        << "directions" << nl
                        << "Please check that the wedge patches are correctly "
                        << "defined"
                        << abort(FatalError);
                }
            }
        }
    }
}


void Foam::solidModel::makeThermalModel() const
{
    if (!thermalPtr_.empty())
    {
        FatalErrorIn("void Foam::solidModel::makeThermalModel() const")
            << "pointer already set!" << abort(FatalError);
    }

    thermalPtr_.set
    (
        new thermalModel(mesh())
    );
}


void Foam::solidModel::makeMechanicalProperties() const
{
    if (!mechanicalPropertiesPtr_.empty())
    {
        FatalErrorInFunction
            << "pointer already set!" << abort(FatalError);
    }

    mechanicalPropertiesPtr_.set
    (
        new IOdictionary
        (
            IOobject
            (
                "mechanicalProperties",
                mesh().time().constant(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


void Foam::solidModel::checkRemovedMechanicalModelEntries() const
{
    // The legacy laws under-relaxed their plastic strain increment, so the
    // constitutive update lagged the displacement between outer iterations and
    // its own convergence was tested against materialTolerance. A
    // mechanicalConstitutiveLaw is a pure function of the kinematics and the
    // old-time state, so there is nothing for a material residual to measure,
    // and convergence is governed by the displacement residuals alone
    if (solidModelDict().found("materialTolerance"))
    {
        Info<< "    'materialTolerance' is ignored and can be removed: the "
            << "mechanicalConstitutiveLaw framework has no material residual"
            << endl;
    }

    // The switch chose between the mechanicalConstitutiveLaw framework and the
    // legacy mechanicalModel. The legacy model has been removed, so a case
    // that asks for it must stop rather than silently run on the framework
    // and give different results; one that asks for the framework gets what
    // it asked for, and is only told that the entry is no longer needed
    const word key("useMechanicalConstitutiveLawManager");

    if (!solidModelDict().found(key))
    {
        return;
    }

    if (Switch(solidModelDict().lookup(key)))
    {
        Info<< "    '" << key << " yes' is obsolete and can be removed: the "
            << "mechanicalConstitutiveLaw framework is the only mechanical "
            << "model" << endl;

        return;
    }

    FatalIOErrorInFunction(solidModelDict())
        << "'" << key << " no' selects the legacy mechanicalModel, which has "
        << "been removed from solids4foam." << nl << nl
        << "    The mechanicalConstitutiveLaw framework is now the only "
        << "mechanical model. It reads the same constant/mechanicalProperties, "
        << "and a law it does not provide stops the run at construction. "
        << "Remove '" << key
        << "' from " << type_ << "Coeffs in constant/solidProperties to run "
        << "on it, and check the results: they are not guaranteed to match "
        << "the legacy model's." << nl
        << "    To reproduce a result from the legacy model, use a solids4foam "
        << "release that still contains it."
        << exit(FatalIOError);
}


void Foam::solidModel::makeRho() const
{
    if (!rhoPtr_.empty())
    {
        FatalErrorIn("void Foam::solidModel::makeRho() const")
            << "pointer already set!" << abort(FatalError);
    }

    // Built from a tmp, so the field takes over the tmp's registration as
    // "rho"
    rhoPtr_.set(new volScalarField(initialRho()));
}


void Foam::solidModel::makeU() const
{
    if (!UPtr_.empty())
    {
        FatalErrorIn("void Foam::solidModel::makeU() const")
            << "pointer already set!" << abort(FatalError);
    }

    UPtr_.set
    (
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector("0", dimLength/dimTime, vector::zero)
        )
    );
}


void Foam::solidModel::makeP() const
{
    if (!pPtr_.empty())
    {
        FatalErrorIn("void Foam::solidModel::makep() const")
            << "pointer already set!" << abort(FatalError);
    }

    pPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "p",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimPressure, 0.0),
            "zeroGradient"
        )
    );
}


namespace
{
    template<class Type>
    Foam::boolList fixedValuePatchMask
    (
        const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&
            field
    )
    {
        Foam::boolList includePatchInStencils
        (
            field.boundaryField().size(),
            false
        );

        forAll(includePatchInStencils, patchI)
        {
            includePatchInStencils[patchI] =
                field.boundaryField()[patchI].fixesValue();
        }

        return includePatchInStencils;
    }
}


const Foam::dictionary& Foam::solidModel::highOrderCoeffsDict() const
{
    if (!solidModelDict().found("highOrderCoeffs"))
    {
        FatalErrorInFunction
            << "The solid model dictionary does not contain 'highOrderCoeffs'."
            << abort(FatalError);
    }

    return solidModelDict().subDict("highOrderCoeffs");
}


const Foam::dictionary& Foam::solidModel::displacementHighOrderCoeffs() const
{
    const dictionary& hoDict = highOrderCoeffsDict();

    if (!hoDict.found("displacement"))
    {
        FatalErrorInFunction
            << "Expected 'highOrderCoeffs.displacement' to be defined."
            << abort(FatalError);
    }

    return hoDict.subDict("displacement");
}


const Foam::dictionary& Foam::solidModel::pressureHighOrderCoeffs() const
{
    const dictionary& hoDict = highOrderCoeffsDict();

    if (!hoDict.found("pressure"))
    {
        FatalErrorInFunction
            << "Expected 'highOrderCoeffs.pressure' to be defined."
            << abort(FatalError);
    }

    return hoDict.subDict("pressure");
}


const Foam::leastSquaresScheme&
Foam::solidModel::displacementLeastSquares() const
{
    // Note: the mask is taken from the primary solution field: for incremental
    // solid models this is DD, as it is that field which carries the boundary
    // conditions
    const leastSquaresReconstruction& reconstructions =
        leastSquaresReconstruction::New(mesh());

    return reconstructions.scheme
    (
        "displacement",
        fixedValuePatchMask(incremental() ? DD_ : D_),
        displacementHighOrderCoeffs()
    );
}


const Foam::leastSquaresScheme&
Foam::solidModel::pressureLeastSquares() const
{
    const leastSquaresReconstruction& reconstructions =
        leastSquaresReconstruction::New(mesh());

    return reconstructions.scheme
    (
        "pressure",
        fixedValuePatchMask(p()),
        pressureHighOrderCoeffs()
    );
}


void Foam::solidModel::makeSigmaQuad() const
{
    if (!sigmaQuadPtr_.empty())
    {
        FatalErrorInFunction
            << "pointer already set!" << abort(FatalError);
    }

    auto& faceQuadPts = compactListListCRef
    (
        displacementLeastSquares().quadrature().faceQuadPoints()
    );

    labelList rowSizes(faceQuadPts.size(), 0);
    forAll(faceQuadPts, faceI)
    {
        rowSizes[faceI] = faceQuadPts[faceI].size();
    }

    sigmaQuadPtr_.set(new CompactListList<symmTensor>(rowSizes));

    CompactListList<symmTensor>& sigmaQuad = sigmaQuadPtr_();
    forAll(sigmaQuad, faceI)
    {
        forAll(sigmaQuad[faceI], qpI)
        {
            sigmaQuad[faceI][qpI] = symmTensor::zero;
        }
    }
}


void Foam::solidModel::makeGradDQuad() const
{
    if (!gradDQuadPtr_.empty())
    {
        FatalErrorInFunction
            << "pointer already set!" << abort(FatalError);
    }

    auto& faceQuadPts = compactListListCRef
    (
        displacementLeastSquares().quadrature().faceQuadPoints()
    );

    labelList rowSizes(faceQuadPts.size(), 0);
    forAll(faceQuadPts, faceI)
    {
        rowSizes[faceI] = faceQuadPts[faceI].size();
    }

    gradDQuadPtr_.set(new CompactListList<tensor>(rowSizes));

    CompactListList<tensor>& gradDQ = gradDQuadPtr_();
    forAll(gradDQ, faceI)
    {
        forAll(gradDQ[faceI], qpI)
        {
            gradDQ[faceI][qpI] = tensor::zero;
        }
    }
}


const Foam::CompactListList<Foam::tensor>&
Foam::solidModel::gradDQuad0() const
{
    if (gradDQuad0Ptr_.empty())
    {
#ifdef OPENFOAM_NOT_EXTEND
        if (gradDQuadPtr_.empty())
        {
            makeGradDQuad();
        }

        // Build the previous-time store from the previous-time displacement,
        // rather than copying the current gradient into it: on a restart the
        // current gradient is not the previous step's
        gradDQuad0Ptr_.set
        (
            new CompactListList<tensor>(gradDQuadPtr_().sizes())
        );
        hofvc::fGrad(D_.oldTime(), gradDQuad0Ptr_());
#else
        gradDQuad0Ptr_.set(new CompactListList<tensor>());
        copyQuadGradient(gradDQuad(), gradDQuad0Ptr_());
#endif
    }

    return autoPtrRef(gradDQuad0Ptr_);
}


void Foam::solidModel::correctPointDisplacement
(
    pointVectorField& pointD
) const
{
    pointD.correctBoundaryConditions();

    const polyMesh& pMesh = pointD.mesh().mesh();
    vectorField& pointDI = pointD;

    forAll(pMesh.boundaryMesh(), patchI)
    {
        if (isA<symmetryPolyPatch>(pMesh.boundaryMesh()[patchI]))
        {
            const polyPatch& patch = pMesh.boundaryMesh()[patchI];

            if (returnReduce(patch.size(), sumOp<int>()) == 0)
            {
                continue;
            }

            const labelList& meshPoints = patch.meshPoints();
            const vector avgN = gAverage(patch.pointNormals());

            forAll(meshPoints, pointI)
            {
                vector& pointValue = pointDI[meshPoints[pointI]];

                if (mag(avgN.x()) > 0.95)
                {
                    pointValue.x() = 0;
                }
                else if (mag(avgN.y()) > 0.95)
                {
                    pointValue.y() = 0;
                }
                else if (mag(avgN.z()) > 0.95)
                {
                    pointValue.z() = 0;
                }
            }
        }
    }

    twoDCorrector_.correctPoints(pointDI);
}


const Foam::pointVectorField& Foam::solidModel::pointDorPointDD() const
{
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        // Updated Lagrangian approaches move the mesh at the end of each
        // time-step so we use the increment of displacement field to calculate
        // the current deformed face zone points
        return pointDD();
    }
    else
    {
        // As linearGeometry and total Lagrangian approaches do not move the
        // mesh, we use the total displacement field to calculate the current
        // deformed face zone points
        return pointD();
    }
}


void Foam::solidModel::makeSetCellDisps() const
{
    if (setCellDispsPtr_.valid())
    {
        FatalErrorIn(type() + "::makeSetCellDisps() const")
            << "pointer already set!" << abort(FatalError);
    }

    if (solidModelDict().found("cellDisplacements"))
    {
        setCellDispsPtr_.set
        (
            new setCellDisplacements
            (
                mesh(), solidModelDict().subDict("cellDisplacements")
            )
        );
    }
    else
    {
        dictionary dict;
        setCellDispsPtr_.set(new setCellDisplacements(mesh(), dict));
    }
}


const Foam::setCellDisplacements& Foam::solidModel::setCellDisps() const
{
    if (setCellDispsPtr_.empty())
    {
        makeSetCellDisps();
    }

    return setCellDispsPtr_();
}


// * * * * * * * * * * Protected Member Function * * * * * * * * * * * * * * //

const Foam::meshDual& Foam::solidModel::dualMesh() const
{
    if (dualMeshPtr_.empty())
    {
        makeDualMesh();
    }

    return dualMeshPtr_();
}


Foam::meshDual& Foam::solidModel::dualMesh()
{
    if (dualMeshPtr_.empty())
    {
        makeDualMesh();
    }

    return dualMeshPtr_();
}


Foam::thermalModel& Foam::solidModel::thermal()
{
    if (thermalPtr_.empty())
    {
        makeThermalModel();
    }

    return thermalPtr_();
}


bool Foam::solidModel::newTimeStep() const
{
    if (curTimeIndex_ != runTime().timeIndex())
    {
        curTimeIndex_ = runTime().timeIndex();
        return true;
    }

    return false;
}


Foam::volScalarField& Foam::solidModel::rho()
{
    if (rhoPtr_.empty())
    {
        makeRho();
    }

    return rhoPtr_();
}


void Foam::solidModel::setCellDisps(fvVectorMatrix& DEqn)
{
    if (setCellDisps().cellIDs().size() == 0)
    {
        return;
    }

    if (incremental())
    {
        // Prepare the list of incremental displacements
        const vectorField& Dold = D().oldTime().internalField();
        const vectorField cellDisps = setCellDisps().cellDisps();
        vectorField cellIncrDisps(cellDisps.size(), vector::zero);
        const labelList cellIDs = setCellDisps().cellIDs();
        forAll(cellIncrDisps, cI)
        {
            cellIncrDisps[cI] = cellDisps[cI] - Dold[cellIDs[cI]];
        }

        DEqn.setValues(cellIDs, cellIncrDisps);
    }
    else
    {
        DEqn.setValues(setCellDisps().cellIDs(), setCellDisps().cellDisps());
    }
}


void Foam::solidModel::relaxField(volVectorField& D, int iCorr)
{
    // Hack to avoid expensive copy of residuals
#ifdef OPENFOAM_COM
    #if (OPENFOAM >= 2312)
        const_cast<dictionary&>
        (
            D.mesh().data().solverPerformanceDict()
        ).clear();
    #else
        const_cast<dictionary&>(D.mesh().solverPerformanceDict()).clear();
    #endif
#endif

    if (relaxationMethod_ == "fixed")
    {
        // Fixed under-relaxation
        D.relax();
    }
    else if (relaxationMethod_ == "Aitken")
    {
        // See Aitken method at:
        // http://empire-multiphysics.com/projects/empire/wiki/Aitken_Relaxation
        // and
        // A partitioned solution approach for electro-thermo-
        // problems, Patrick Erbts, Stefan Hartmann, Alexander Duster.

        // Store aitkenResidual previous iteration
        aitkenResidual_.storePrevIter();

        // Calculate new aitkenResidual
        aitkenResidual_ = D.prevIter() - D;

        if (iCorr == 0)
        {
            // Fixed under-relaxation is applied on the first iteration
            aitkenAlpha_ = 1.0;

#ifdef OPENFOAM_NOT_EXTEND
            if (mesh().relaxField(D.name()))
            {
                aitkenAlpha_ =
                    mesh().fieldRelaxationFactor(D.name());
            }
#else
            if (mesh().solutionDict().relaxField(D.name()))
            {
                aitkenAlpha_ =
                    mesh().solutionDict().fieldRelaxationFactor(D.name());
            }
#endif
        }
        else
        {
            const volVectorField aitkenResidualDelta
            (
                aitkenResidual_.prevIter() - aitkenResidual_
            );

            // Update the relaxation factor field
            aitkenAlpha_ =
                aitkenAlpha_*(aitkenResidual_.prevIter() & aitkenResidualDelta)
               /(
                    magSqr(aitkenResidualDelta)
                  + dimensionedScalar("SMALL", dimLength*dimLength, SMALL)
                );

            // Bound alpha between 0.0 and 2.0
            // This may not be necessary but it seems to help convergence
            aitkenAlpha_ = max(0.0, min(2.0, aitkenAlpha_));
        }

        // Relax the field
        D -= aitkenAlpha_*aitkenResidual_;
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::relaxField(volVectorField& D, int iCorr)"
        )   << "relaxationMethod '" << relaxationMethod_ << "' unknown!"
            << " Options are fixed or Aitken" << abort(FatalError);
    }
}


Foam::dictionary& Foam::solidModel::solidModelDict()
{
    return solidProperties_.subDict(type_ + "Coeffs");
}


void Foam::solidModel::makeRhoD2dt2D() const
{
    if (rhoD2dt2DPtr_.valid())
    {
        FatalErrorIn("void Foam::solidModel::makeRhoD2dt2D() const")
            << "Pointer already set" << abort(FatalError);
    }

    rhoD2dt2DPtr_.set
    (
        new volVectorField
        (
            IOobject
            (
                "rhoD2dt2D",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector("zero", dimForce/dimVolume, vector::zero)
        )
    );
}


Foam::volVectorField& Foam::solidModel::rhoD2dt2D() const
{
    if (rhoD2dt2DPtr_.empty())
    {
        makeRhoD2dt2D();
    }

    return rhoD2dt2DPtr_();
}


Foam::volScalarField& Foam::solidModel::p()
{
    if (pPtr_.empty())
    {
        makeP();
    }

    return autoPtrRef(pPtr_);
}


const Foam::volScalarField& Foam::solidModel::p() const
{
    if (pPtr_.empty())
    {
        makeP();
    }

    return pPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModel::solidModel
(
    const word& type,
    Time& runTime,
    const word& region
)
:
    physicsModel(type, runTime, region),
    regIOobject // ZT, Jul18: allow for multiple solid regions
    (
        IOobject
        (
            "solidModel_" + region,
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    meshOwnerPtr_
    (
        runTime.foundObject<dynamicFvMesh>(region)
      ? Foam::autoPtr<Foam::dynamicFvMesh>()
      : dynamicFvMesh::New
        (
            IOobject
            (
                region,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    ),
    meshPtr_
    (
        meshOwnerPtr_.valid()
      ? meshOwnerPtr_.ptr()
      : &const_cast<dynamicFvMesh&>(runTime.lookupObject<dynamicFvMesh>(region))
    ),
    dualMeshPtr_(),
    solidProperties_
    (
        // If region == "region0" then read from the main case
        // Otherwise, read from the region/sub-mesh directory e.g.
        // constant/fluid or constant/solid
        bool(region == dynamicFvMesh::defaultRegion)
      ? IOobject
        (
            "solidProperties",
            runTime.caseConstant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
      : IOobject
        (
            "solidProperties",
            runTime.caseConstant(),
            region, // using 'local' property of IOobject
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    type_(type),
    solutionAlgorithm_
    (
        solidModelDict().found("solutionAlgorithm")
#ifdef OPENFOAM_COM
      ? solutionAlgorithmNames_.get("solutionAlgorithm", solidModelDict())
#else
      ? solutionAlgorithmNames_.read(solidModelDict().lookup("solutionAlgorithm"))
#endif
      : solutionAlgorithm::IMPLICIT_SEGREGATED
    ),
    thermalPtr_(),
    useBoundaryFaceValuesD_
    (
        IOobject
        (
            "useBoundaryFaceValues_D",
            runTime.constant(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        boolList(mesh().boundary().size(), false)
    ),
    useBoundaryFaceValuesDD_
    (
        IOobject
        (
            "useBoundaryFaceValues_DD",
            runTime.constant(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        boolList(mesh().boundary().size(), false)
    ),
    useBoundaryFaceValuesp_
    (
        IOobject
        (
            "useBoundaryFaceValues_p",
            runTime.constant(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        boolList(mesh().boundary().size(), false)
    ),
    Dheader_("D", runTime.timeName(), mesh(), IOobject::MUST_READ),
    DDheader_("DD", runTime.timeName(), mesh(), IOobject::MUST_READ),
    pointDheader_("pointD", runTime.timeName(), mesh(), IOobject::MUST_READ),
    D_
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    DD_
    (
        IOobject
        (
            "DD",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    UPtr_(),
    pMesh_(pointMesh::New(*meshPtr_)),
    pointD_
    (
        IOobject
        (
            "pointD",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    pointDD_
    (
        IOobject
        (
            "pointDD",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    gradD_
    (
        IOobject
        (
            "grad(" + D_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradDQuadPtr_(),
    gradDD_
    (
        IOobject
        (
            "grad(" + DD_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    sigmaQuadPtr_(),
    curTimeIndex_(-1),
    rhoPtr_(),
    g_
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    dampingCoeff_
    (
        solidModelDict().lookupOrAddDefault<dimensionedScalar>
        (
            "dampingCoeff",
            dimensionedScalar("dampingCoeff", dimless/dimTime, 0)
        )
    ),
    momentumStabilisationPtr_(),
    pressureStabilisationPtr_(),
    rAUfPtr_(),
    solutionTol_
    (
        solidModelDict().lookupOrAddDefault<scalar>("solutionTolerance", 1e-06)
    ),
    alternativeTol_
    (
        solidModelDict().lookupOrAddDefault<scalar>
        (
            "alternativeTolerance", 1e-07
        )
    ),
    infoFrequency_
    (
        solidModelDict().lookupOrAddDefault<int>("infoFrequency", 100)
    ),
    nCorr_(solidModelDict().lookupOrAddDefault<int>("nCorrectors", 10000)),
    maxIterReached_(0),
    residualFilePtr_(),
    writeResidualField_
    (
        solidModelDict().lookupOrAddDefault<Switch>("writeResidualField", false)
    ),
    enforceLinear_(false),
    relaxationMethod_
    (
        solidModelDict().lookupOrAddDefault<word>("relaxationMethod", "fixed")
    ),
    aitkenAlpha_
    (
        IOobject
        (
            "aitkenAlpha",
            runTime.constant(),
            *meshPtr_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *meshPtr_,
        dimensionedScalar("one", dimless, 1.0)
    ),
    aitkenResidual_
    (
        IOobject
        (
            "aitkenResidual",
            runTime.constant(),
            *meshPtr_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *meshPtr_,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    globalPatchesPtrList_(),
    setCellDispsPtr_(),
    mechanicalManagerPtr_(),
    jacobianTangentCached_(false),
    jacobianTangent_(tangentRequest::none),
    restartSpecified_(solidModelDict().found("restart")),
    restart_
    (
        solidModelDict().lookupOrAddDefault<Switch>("restart", false)
    ),
    restartKinematicsAvailable_(true),
    rhoD2dt2DPtr_(),
    twoDCorrector_(mesh()),
    twoD_(mesh().nGeometricD() == 2),
    solvePressure_
    (
        solidModelDict().lookupOrDefault<Switch>("solvePressure", false)
    ),
    highOrderJacobian_(false),
    highOrderResidual_(false),
    pPtr_(),
    hydrostaticSmoothingRead_(false),
    hydrostaticSmoothingRequested_(false),
    hydrostaticSmoothingChecked_(false),
    pressureSmoothingScaleFactor_(100.0),
    sigmaHydPtr_(),
    useBoundaryFaceValuesSigmaHydPtr_(),
    gradSigmaHydPtr_(),
    smoothingVolumetricResponsePtr_()
#ifdef OPENFOAM_COM
    ,
    fvOptions_(fv::options::New(*meshPtr_))
#endif
{
    // As far as the finite volume discretisation is concerned, the solid mesh
    // is not a "moving mesh": the solid models which do update the mesh, e.g.
    // the updated Lagrangian ones, explicitly clear the moving flag after each
    // mesh motion, as the mesh motion is not a flow of material through the
    // faces.
    // On a restart, however, the fvMesh constructor marks the mesh as moving if
    // it finds a 'meshPhi' or 'V0' field in the start time directory, and then,
    // for example, backwardD2dt2Scheme stops with a "not implemented for a
    // moving mesh" error (issue #184). So the flag is cleared here
    mesh().moving(false);

    checkRemovedMechanicalModelEntries();

    // Set the useBoundaryFaceValues fields
    forAll(useBoundaryFaceValuesD_, patchI)
    {
        if
        (
            isA<solidTractionFvPatchVectorField>
            (
                D_.boundaryField()[patchI]
            )
         || isA<fixedDisplacementZeroShearFvPatchVectorField>
            (
                D_.boundaryField()[patchI]
            )
        )
        {
            useBoundaryFaceValuesD_[patchI] = false;
        }
        else
        {
            useBoundaryFaceValuesD_[patchI] = true;
        }
    }
    forAll(useBoundaryFaceValuesDD_, patchI)
    {
        if
        (
            isA<solidTractionFvPatchVectorField>
            (
                DD_.boundaryField()[patchI]
            )
         || isA<fixedDisplacementZeroShearFvPatchVectorField>
            (
                DD_.boundaryField()[patchI]
            )
        )
        {
            useBoundaryFaceValuesDD_[patchI] = false;
        }
        else
        {
            useBoundaryFaceValuesDD_[patchI] = true;
        }
    }
    if (pPtr_.valid())
    {
        forAll(useBoundaryFaceValuesp_, patchI)
        {
            if (pPtr_->boundaryField()[patchI].fixesValue())
            {
                useBoundaryFaceValuesp_[patchI] = true;
            }
        }
    }

    // Force old time fields to be stored
    D_.oldTime().oldTime();
    DD_.oldTime().oldTime();
    pointD_.oldTime();
    pointDD_.oldTime();
    gradD_.oldTime();
    gradDD_.oldTime();
    sigma_.oldTime();

    // Whether a continued run has the old-time gradient a constitutive
    // history is measured against. Looked for rather than assumed from
    // 'restart yes', since the run that wrote this time may not have asked
    // for it
    if (runTime.startTimeIndex() > 0)
    {
        IOobject gradD0IO
        (
            gradD_.oldTime().name(),
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        );

#ifdef OPENFOAM_NOT_EXTEND
        restartKinematicsAvailable_ =
            gradD0IO.typeHeaderOk<volTensorField>(false);
#else
        restartKinematicsAvailable_ = gradD0IO.headerOk();
#endif
    }

    if (restart_)
    {
        // Enable writing of fields which are needed for restart
        D_.oldTime().writeOpt() = IOobject::AUTO_WRITE;
        D_.oldTime().oldTime().writeOpt() = IOobject::AUTO_WRITE;
        DD_.writeOpt() = IOobject::AUTO_WRITE;
        DD_.oldTime().writeOpt() = IOobject::AUTO_WRITE;
        DD_.oldTime().oldTime().writeOpt() = IOobject::AUTO_WRITE;
        pointD_.writeOpt() = IOobject::AUTO_WRITE;
        pointD_.oldTime().writeOpt() = IOobject::AUTO_WRITE;
        pointDD_.writeOpt() = IOobject::AUTO_WRITE;
        gradD_.writeOpt() = IOobject::AUTO_WRITE;
        gradD_.oldTime().writeOpt() = IOobject::AUTO_WRITE;
        gradDD_.writeOpt() = IOobject::AUTO_WRITE;
    }
    else
    {
        // Starting from a time that is not the first, without having been
        // asked to write what a restart needs. The fields below are switched
        // off in this branch, so the run is about to continue from a state it
        // only partly has: the displacement comes back, the increment and the
        // gradient it is measured against do not.
        //
        // Whether that matters depends on the material. A law written in total
        // strain will not notice; an incremental one reads the whole run's
        // strain as a single step's and is wrong by tens of percent while
        // running happily to the end. That is too quiet a way to be wrong to
        // leave unsaid, and too common a mistake to assume: the flag defaults
        // to off and most cases never set it
        if (runTime.startTimeIndex() > 0)
        {
            // Continuing from a time that is not the first, without having
            // been asked to write what a restart needs. Whether that matters
            // depends on the material: one written in total strain will not
            // notice, while an incremental one reads the whole run's strain as
            // a single step's and is wrong by tens of percent while running
            // happily to the end.
            //
            // The fields may still be there, if the run that produced this
            // time directory did ask for them, so look before complaining
            const bool present = restartKinematicsAvailable_;

            if (!present && !restartSpecified_)
            {
                // The case has not said anything about restarting, and is
                // restarting. Refuse: this is a mistake far more often than it
                // is a choice, and the cost of being wrong is a plausible
                // answer rather than an obvious failure
                FatalErrorInFunction
                    << "Continuing from time " << runTime.timeName()
                    << ", but the fields a consistent restart needs were "
                    << "never written." << nl << nl
                    << "    The displacement increment, the old-time "
                    << "displacement gradient and the point fields are only "
                    << "written when the case asks for them." << nl << nl
                    << "    Either" << nl << nl
                    << "        restart yes;" << nl << nl
                    << "    in the solidModel's coefficients dictionary, and "
                    << "run again from the start; or" << nl << nl
                    << "        restart no;" << nl << nl
                    << "    to say that this material does not need them and "
                    << "continue. A material written in total strain does not; "
                    << "an incremental one does, and without them reads the "
                    << "whole run's strain as one step's."
                    << exit(FatalError);
            }
            else if (!present)
            {
                // Said 'no' deliberately. Their call, said once
                WarningInFunction
                    << "Continuing from time " << runTime.timeName()
                    << " with 'restart no': the displacement increment and "
                    << "old-time gradient were not written, so an incremental "
                    << "material would continue from the wrong strain."
                    << endl;
            }
        }

        D_.oldTime().writeOpt() = IOobject::NO_WRITE;
        D_.oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
        DD_.writeOpt() = IOobject::NO_WRITE;
        DD_.oldTime().writeOpt() = IOobject::NO_WRITE;
        DD_.oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
        pointD_.writeOpt() = IOobject::AUTO_WRITE;
        pointD_.oldTime().writeOpt() = IOobject::NO_WRITE;
        pointDD_.writeOpt() = IOobject::NO_WRITE;
        gradD_.writeOpt() = IOobject::NO_WRITE;
        gradD_.oldTime().writeOpt() = IOobject::NO_WRITE;
        gradDD_.writeOpt() = IOobject::NO_WRITE;
    }

    if (solidModelDict().found("highOrderCoeffs"))
    {
        const dictionary& hoDict = highOrderCoeffsDict();

        if (!hoDict.found("displacement"))
        {
            FatalErrorInFunction
                << "Expected 'highOrderCoeffs.displacement' to be defined."
                << abort(FatalError);
        }

        highOrderJacobian_ =
            hoDict.lookupOrDefault<Switch>("highOrderJacobian", false);

        highOrderResidual_ =
            hoDict.lookupOrDefault<Switch>("highOrderResidual", false);

        if
        (
            (highOrderJacobian_ || highOrderResidual_)
         && solutionAlg() != solutionAlgorithm::PETSC_SNES
        )
        {
            FatalErrorInFunction
                << "highOrderResidual/highOrderJacobian are only supported "
                << "with "
                << solidModel::solutionAlgorithmNames_
                   [solidModel::solutionAlgorithm::PETSC_SNES]
                << abort(FatalError);
        }
    }

    // Print out the relaxation factor
    Info<< "    under-relaxation method: " << relaxationMethod_ << endl;

    // If requested, create the residual file
    if (solidModelDict().lookupOrAddDefault<Switch>("residualFile", false))
    {
        if (Pstream::master())
        {
            Info<< "Creating residual.dat" << endl;
            residualFilePtr_.set
            (
                new OFstream(runTime.path()/"residual.dat")
            );
        }
    }

    // Create momentum stabilisation

    dictionary defaultStabSubDict;
    defaultStabSubDict.add("type", "diffStencilLaplacian");
    defaultStabSubDict.add("scaleFactor", 0.1);

    if (!solidModelDict().found("stabilisation"))
    {
        // If the stabilisation sub-dict is not found, we will add it with
        // default settings
        dictionary stabDict;
        stabDict.add("momentum", defaultStabSubDict);
        stabDict.add("pressure", defaultStabSubDict);

        // Add stabilisation dict
        solidModelDict().add("stabilisation", stabDict);
    }

    dictionary& stabDict = solidModelDict().subDict("stabilisation");

    // Check for previous stabilisation definition
    if (stabDict.found("type") || stabDict.found("scaleFactor"))
    {
        FatalErrorInFunction
            << "Found 'type' or 'scaleFactor' in stabilisation subDict of "
            << "solidProperties: this is the old format. Instead, define a "
            << "stabilisation/momentum sub-dict" << exit(FatalError);
    }

    if (!stabDict.found("momentum"))
    {
        stabDict.add("momentum", defaultStabSubDict);
    }

    if (!stabDict.found("pressure"))
    {
        stabDict.add("pressure", defaultStabSubDict);
    }

    momentumStabilisationPtr_ =
        stabilisationModel::New
        (
            mesh(),
            stabDict.subDict("momentum"),
            dimless
        );

    // Only stabilisation models that support high-order residual/Jacobian
    // calculation are allowed when high-order is enabled
    if
    (
        (highOrderResidual() || highOrderJacobian())
     && !momentumStabilisationPtr_->supportsHighOrderResidual()
    )
    {
        FatalErrorInFunction
            << "Momentum stabilisation type "
            << momentumStabilisationPtr_->type()
            << " does not support high-order residual or Jacobian "
            << "calculation" << abort(FatalError);
    }

    pressureStabilisationPtr_ =
        stabilisationModel::New
        (
            mesh(),
            stabDict.subDict("pressure"),
            dimPressure/dimLength
        );

#ifdef OPENFOAM_COM
    if (!fvOptions_.optionList::size())
    {
        Info<< "No finite volume options present" << endl;
    }
#endif

    // If the case is axisymmetric, we will disable solving in the out-of-plane
    // direction
    // PC, 12-Nov-18: disabling the 3rd direction slows down convergence a lot
    // in some elastic cases: disabled for now
    //checkWedges();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidModel::~solidModel()
{
    thermalPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::solidModel::rho() const
{
    if (rhoPtr_.empty())
    {
        makeRho();
    }

    return rhoPtr_();
}

const Foam::thermalModel& Foam::solidModel::thermal() const
{
    if (thermalPtr_.empty())
    {
        makeThermalModel();
    }

    return thermalPtr_();
}


const Foam::IOdictionary& Foam::solidModel::mechanicalProperties() const
{
    if (mechanicalPropertiesPtr_.empty())
    {
        makeMechanicalProperties();
    }

    return mechanicalPropertiesPtr_();
}


namespace Foam
{
    // Search a law's dictionary and every sub-dictionary for a smoothing
    // request. A wrapper law - thermoMechanicalLaw, poroMechanicalLaw,
    // electroMechanicalLaw - keeps its inner law in a sub-dictionary, and a
    // case could set solvePressureEqn there, where the inner law read it
    static bool findSmoothingRequest
    (
        const dictionary& dict,
        scalar& scaleFactor
    )
    {
        if (dict.lookupOrDefault<Switch>("solvePressureEqn", false))
        {
            scaleFactor =
                dict.lookupOrDefault<scalar>
                (
                    "pressureSmoothingScaleFactor", 100.0
                );

            return true;
        }

        forAllConstIter(dictionary, dict, iter)
        {
            if
            (
                iter().isDict()
             && findSmoothingRequest(iter().dict(), scaleFactor)
            )
            {
                return true;
            }
        }

        return false;
    }
}


void Foam::solidModel::readHydrostaticSmoothing() const
{
    hydrostaticSmoothingRead_ = true;

    // Read from the law entries, where cases have always set it. It is a
    // solid-model setting in all but location: one equation is solved, with
    // the solid model's own momentum diagonal, so one answer is needed for
    // the whole mesh
    const PtrList<entry> lawEntries(mechanicalProperties().lookup("mechanical"));

    label nRequested = 0;
    word requestingLaw;

    forAll(lawEntries, lawI)
    {
        scalar scaleFactor = 100.0;

        if (findSmoothingRequest(lawEntries[lawI].dict(), scaleFactor))
        {
            nRequested++;
            requestingLaw = lawEntries[lawI].keyword();
            pressureSmoothingScaleFactor_ = scaleFactor;
        }
    }

    if (nRequested == 0)
    {
        return;
    }

    // Not supported for more than one material. One equation over the whole
    // mesh would diffuse the hydrostatic stress across a material interface,
    // where it genuinely jumps
    if (lawEntries.size() > 1)
    {
        FatalIOErrorInFunction(mechanicalProperties())
            << "solvePressureEqn is set for material '" << requestingLaw
            << "', and there are " << lawEntries.size() << " materials."
            << nl << nl
            << "    The hydrostatic stress smoothing is done by the solid "
            << "model, as one equation over the whole mesh, and that would "
            << "smear the hydrostatic stress across the material interfaces, "
            << "where it jumps. It is therefore supported for a single "
            << "material only." << nl
            << "    Remove solvePressureEqn."
            << exit(FatalIOError);
    }

    hydrostaticSmoothingRequested_ = true;
}


bool Foam::solidModel::hydrostaticSmoothingRequested() const
{
    if (!hydrostaticSmoothingRead_)
    {
        readHydrostaticSmoothing();
    }

    return hydrostaticSmoothingRequested_;
}


bool Foam::solidModel::smoothHydrostaticStress() const
{
    if (!hydrostaticSmoothingRequested())
    {
        return false;
    }

    if (hydrostaticSmoothingChecked_)
    {
        return true;
    }

    if (solvePressure())
    {
        FatalErrorInFunction
            << "solvePressureEqn is set in mechanicalProperties and "
            << "solvePressure in the " << type() << " coefficients." << nl
            << "    They are alternatives: the mixed formulation solves for "
            << "the pressure itself, and the smoothing would be applied to a "
            << "volumetric response it has already replaced. Remove one."
            << exit(FatalError);
    }

    if (solutionAlg() != solutionAlgorithm::IMPLICIT_SEGREGATED)
    {
        FatalErrorInFunction
            << "solvePressureEqn is set in mechanicalProperties, and the "
            << type() << " solution algorithm is "
            << solutionAlgorithmNames_[solutionAlg()] << "." << nl
            << "    The hydrostatic stress smoothing is scaled by the diagonal "
            << "of the segregated momentum equation, so it is available with "
            << "the implicitSegregated algorithm only. For another algorithm "
            << "use the mixed formulation (solvePressure) instead."
            << exit(FatalError);
    }

    // Asked by name here rather than left to the split update, whose refusal
    // is written for the mixed formulation
    if (!mechanicalManager().allLawsProvideVolumetricSplit())
    {
        const PtrList<entry> lawEntries
        (
            mechanicalProperties().lookup("mechanical")
        );

        FatalErrorInFunction
            << "solvePressureEqn is set for material '"
            << lawEntries[0].keyword() << "', of type "
            << word(lawEntries[0].dict().lookup("type"))
            << ", and that law cannot separate its isochoric stress from its "
            << "volumetric response." << nl
            << "    The smoothing replaces the volumetric response with a "
            << "smoothed one, so it needs the two apart: taking the trace of "
            << "a total stress would also smooth any spherical stress that "
            << "is not a volumetric response." << nl
            << "    Remove solvePressureEqn, or use a law that separates them."
            << exit(FatalError);
    }

    Info<< type() << ": smoothing the hydrostatic stress (solvePressureEqn), "
        << "with pressureSmoothingScaleFactor " << pressureSmoothingScaleFactor_
        << endl;

    hydrostaticSmoothingChecked_ = true;

    return true;
}


Foam::volScalarField& Foam::solidModel::smoothingVolumetricResponse()
{
    if (smoothingVolumetricResponsePtr_.empty())
    {
        smoothingVolumetricResponsePtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "smoothingVolumetricResponse",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimPressure, 0.0)
            )
        );
    }

    return smoothingVolumetricResponsePtr_();
}


void Foam::solidModel::addSmoothedHydrostaticStress
(
    volSymmTensorField& sigma,
    const volScalarField& impK,
    const volScalarField* JPtr
)
{
    if (sigmaHydPtr_.empty())
    {
        // Not read, zero-gradient, and written
        sigmaHydPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "sigmaHyd",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimPressure, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        // Looked up by name by the leastSquaresS4f gradient scheme
        useBoundaryFaceValuesSigmaHydPtr_.set
        (
            new boolIOList
            (
                IOobject
                (
                    "useBoundaryFaceValues_sigmaHyd",
                    mesh().time().constant(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                boolList(mesh().boundary().size(), false)
            )
        );

        gradSigmaHydPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "grad(sigmaHyd)",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedVector("zero", dimPressure/dimLength, vector::zero)
            )
        );
    }

    volScalarField& sigmaHyd = sigmaHydPtr_();
    volVectorField& gradSigmaHyd = gradSigmaHydPtr_();

    // The solid model registers its momentum diagonal under this name for the
    // duration of each outer iteration
    if (!mesh().foundObject<volScalarField>("DEqnA"))
    {
        FatalErrorInFunction
            << "The hydrostatic stress smoothing (solvePressureEqn) needs the "
            << "momentum equation diagonal, DEqnA, and " << type()
            << " has not registered one at this stress update." << nl
            << "    It exists only inside the segregated momentum loop, so "
            << "the smoothing is not available from a stress update outside "
            << "it, such as the linear predictor."
            << exit(FatalError);
    }

    const volScalarField& AD = mesh().lookupObject<volScalarField>("DEqnA");

    // The explicit hydrostatic Kirchhoff stress, J*dU/dJ.
    //
    // Kirchhoff for every law, deliberately. Most of the legacy laws passed
    // their hydrostatic Kirchhoff stress to the smoothing, but GuccioneElastic
    // passed its Cauchy one, and since multiplying by a varying J does not
    // commute with the smoothing, a Guccione case smoothed this way is
    // stabilised a little differently from before. No case combines the two;
    // one measure for all laws is the simpler contract
    const volScalarField& volumetricResponse = smoothingVolumetricResponse();
    const volScalarField sigmaHydExplicit
    (
        "sigmaHydExplicit",
        JPtr ? (*JPtr)*volumetricResponse : 1.0*volumetricResponse
    );

#ifdef OPENFOAM_NOT_EXTEND
    const int oldDebug = SolverPerformance<scalar>::debug;
    SolverPerformance<scalar>::debug = 0;
#endif

    // Store previous iteration to allow relaxation, if needed
    sigmaHyd.storePrevIter();

    // Pressure diffusivity field
    const surfaceScalarField rDAf
    (
        "rDAf",
        pressureSmoothingScaleFactor_*fvc::interpolate
        (
            impK/AD, "interpolate(" + gradSigmaHyd.name() + ")"
        )
    );

    const dimensionedScalar one("one", dimless, 1.0);

    // The fvm and fvc Laplacians agree for a smooth field, and their
    // difference damps the oscillations of one that is not
    fvScalarMatrix sigmaHydEqn
    (
        fvm::Sp(one, sigmaHyd)
      - fvm::laplacian(rDAf, sigmaHyd, "laplacian(rDA,sigmaHyd)")
     ==
        sigmaHydExplicit
      - fvc::div(rDAf*fvc::interpolate(gradSigmaHyd) & mesh().Sf())
    );

    sigmaHydEqn.solve();

    sigmaHyd.relax();

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<scalar>::debug = oldDebug;
#endif

    gradSigmaHyd = fvc::grad(sigmaHyd);

    // Recombine, including on the boundary, where this takes the
    // zero-gradient value of sigmaHyd
    if (JPtr)
    {
        sigma = sigma + (sigmaHyd/(*JPtr))*I;
    }
    else
    {
        sigma = sigma + sigmaHyd*I;
    }
}


#ifdef OPENFOAM_NOT_EXTEND
const Foam::enhancedVolPointInterpolation& Foam::solidModel::volToPoint() const
{
    return enhancedVolPointInterpolation::New(mesh());
}
#else
const Foam::newLeastSquaresVolPointInterpolation&
Foam::solidModel::volToPoint() const
{
    return newLeastSquaresVolPointInterpolation::New(mesh());
}
#endif


void Foam::solidModel::DisRequired()
{
#ifdef OPENFOAM_NOT_EXTEND
    if (!Dheader_.typeHeaderOk<volVectorField>(true))
#else
    if (!Dheader_.headerOk())
#endif
    {
        FatalErrorIn(type() + "::DisRequired()")
            << "This solidModel requires the 'D' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::solidModel::DDisRequired()
{
#ifdef OPENFOAM_NOT_EXTEND
    if (!DDheader_.typeHeaderOk<volVectorField>(true))
#else
    if (!DDheader_.headerOk())
#endif
    {
        FatalErrorIn(type() + "::DDisRequired()")
            << "This solidModel requires the 'DD' field to be specified!"
            << abort(FatalError);
    }

#ifdef OPENFOAM_NOT_EXTEND
    if (incremental() && !restart() && Dheader_.typeHeaderOk<volVectorField>(true))
#else
    if (incremental() && !restart() && Dheader_.headerOk())
#endif
    {
        FatalErrorIn(type() + "::DDisRequired()")
            << "This solidModel solves for the displacement increment 'DD', "
            << "but a 'D' field was found at the start time." << nl
            << "Remove 'D' from the initial time directory, or set "
            << "'restart true' for a consistent restart." << abort(FatalError);
    }
}


void Foam::solidModel::pointDisRequired()
{
#ifdef OPENFOAM_NOT_EXTEND
    if (!pointDheader_.typeHeaderOk<pointVectorField>(true))
#else
    if (!pointDheader_.headerOk())
#endif
    {
        FatalErrorIn(type() + "::pointDisRequired()")
            << "This solidModel requires the 'pointD' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::solidModel::makeGlobalPatches
(
    const wordList& patchNames,
    const bool currentConfiguration
) const
{
    globalPatchesPtrList_.setSize(patchNames.size());

    forAll(patchNames, i)
    {
        if (globalPatchesPtrList_.set(i))
        {
            FatalErrorIn
            (
                type() + "::makeGlobalPatches(const wordList&) const"
            )   << "Pointer already set for global patch: "
                << patchNames[i] << "!"
                << abort(FatalError);
        }

        // Lookup patch index
        if (mesh().boundaryMesh().findPatchID(patchNames[i]) == -1)
        {
            FatalErrorIn("void Foam::solidModel::makeGlobalPatches(...)")
                << "Patch " << patchNames[i] << " not found!"
                << abort(FatalError);
        }

        // Create the global patch based on the undeformed mesh
        globalPatchesPtrList_.set
        (
            i,
            new globalPolyPatch(patchNames[i], mesh())
        );

        if (currentConfiguration)
        {
            // Force creation of standAlonePatch, so that its points can be set
            // to the current configuration by syncGlobalPatches() below
            globalPatchesPtrList_[i].globalPatch();
        }
    }

    if (currentConfiguration)
    {
        // The global patches are required in the current (deformed)
        // configuration, so displace their points by pointD/pointDD
        // Note: the global patches are always constructed on the undeformed
        // mesh above and then moved here; previously, the mesh itself was
        // temporarily moved to the deformed configuration while the global
        // patches were constructed, and then moved back. That approach had two
        // undesirable side-effects (issue #184):
        //   - fvMesh::movePoints() creates the mesh motion flux field (meshPhi)
        //     and marks the mesh points as AUTO_WRITE; consequently, meshPhi
        //     and polyMesh/points were written to every time directory, even
        //     though the solid mesh is not moved for linear geometry and total
        //     Lagrangian approaches. On restart, the presence of meshPhi makes
        //     OpenFOAM consider the solid mesh to be moving, e.g. causing
        //     backwardD2dt2Scheme to stop with a "not implemented for a moving
        //     mesh" error, and the reduced-precision points written to the time
        //     directory can break the processor patch face matching checks in
        //     parallel.
        //   - globalPolyPatch merges coincident points across processor
        //     boundaries using exact point coordinate comparisons; the
        //     undeformed point coordinates are bit-identical on either side of
        //     a processor boundary, whereas the deformed ones need not be.
        syncGlobalPatches();
    }
}


const Foam::PtrList<Foam::globalPolyPatch>&
Foam::solidModel::globalPatches() const
{
    if (globalPatchesPtrList_.empty())
    {
        FatalErrorIn(type() + "::globalPatches() const")
            << "makeGlobalPatches(const wordList&) must be called "
            << "before globalPatch can be called!"
            << abort(FatalError);
    }

    return globalPatchesPtrList_;
}


void Foam::solidModel::clearGlobalPatches() const
{
    globalPatchesPtrList_.clear();
}


void Foam::solidModel::syncGlobalPatches() const
{
    forAll(globalPatchesPtrList_, i)
    {
        const polyPatch& ppatch = globalPatchesPtrList_[i].patch();

        const vectorField patchPointDisplacement
        (
            pointDorPointDD().internalField(), ppatch.meshPoints()
        );

        const pointField patchPoints
        (
            ppatch.localPoints() + patchPointDisplacement
        );

        globalPatchesPtrList_[i].syncPoints(patchPoints);
    }
}


Foam::vector Foam::solidModel::pointU(const label pointID) const
{
    pointVectorField pointU
    (
        IOobject
        (
            "pointU",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimVelocity, vector::zero)
    );

    volToPoint().interpolate(U(), pointU);

    return pointU.internalField()[pointID];
}


Foam::tmp<Foam::vectorField>
Foam::solidModel::faceZonePointDisplacementIncrement
(
    const label interfaceI
) const
{
    // Create patch point field
    const vectorField patchPointDispIncr
    (
        pointDD().internalField(),
        globalPatches()[interfaceI].patch().meshPoints()
    );

    // Return the global patch field
    return globalPatches()[interfaceI].patchPointToGlobal(patchPointDispIncr);
}


Foam::tmp<Foam::vectorField>
Foam::solidModel::faceZonePointDisplacementOld
(
    const label interfaceI
) const
{
    // Create patch point field
    const vectorField patchPointDispOld
    (
        pointD().oldTime().internalField(),
        globalPatches()[interfaceI].patch().meshPoints()
    );

    // Return the global patch field
    return globalPatches()[interfaceI].patchPointToGlobal(patchPointDispOld);
}


Foam::tmp<Foam::vectorField> Foam::solidModel::faceZoneAcceleration
(
    const label interfaceI
) const
{
    const volVectorField a(fvc::d2dt2(D()));

    return globalPatches()[interfaceI].patchFaceToGlobal
    (
        a.boundaryField()[globalPatches()[interfaceI].patch().index()]
    );
}


void Foam::solidModel::rollOverQuadratureHistory()
{
    // The high-order discretisation, and hence gradDQuad(), does not exist on
    // foam-extend
#ifndef FOAMEXTEND
    if (gradDQuadPtr_.valid())
    {
        if (gradDQuad0Ptr_.empty())
        {
            gradDQuad0Ptr_.set(new CompactListList<tensor>());
        }

        copyQuadGradient(gradDQuad(), gradDQuad0Ptr_());
    }
#endif
}


void Foam::solidModel::updateTotalFields()
{
    // Here rather than in each solid model because every model needs it
    mechanicalManager().endTimeStep();

    // A model overriding this must call rollOverQuadratureHistory() itself
    rollOverQuadratureHistory();
}


// The high-order face quadrature does not run on foam-extend, but these
// helpers are portable and the models that call them are compiled there
void Foam::solidModel::quadDeformationGradient
(
    const CompactListList<tensor>& gradD,
    autoPtr<CompactListList<tensor>>& FPtr
)
{
    if (FPtr.empty())
    {
        FPtr.set(new CompactListList<tensor>());
    }

    FPtr().setSize(gradD.sizes());

    const List<tensor>& gradDv = gradD.m();
    List<tensor>& F = FPtr().m();

    forAll(gradDv, i)
    {
        F[i] = I + gradDv[i].T();
    }
}


void Foam::solidModel::quadInverseAndJacobian
(
    const CompactListList<tensor>& F,
    autoPtr<CompactListList<tensor>>& FinvPtr,
    autoPtr<CompactListList<scalar>>& JPtr
)
{
    if (FinvPtr.empty())
    {
        FinvPtr.set(new CompactListList<tensor>());
        JPtr.set(new CompactListList<scalar>());
    }

    FinvPtr().setSize(F.sizes());
    JPtr().setSize(F.sizes());

    const List<tensor>& Fv = F.m();
    List<tensor>& Finv = FinvPtr().m();
    List<scalar>& J = JPtr().m();

    forAll(Fv, i)
    {
        Finv[i] = inv(Fv[i]);
        J[i] = det(Fv[i]);
    }
}


Foam::mechanicalConstitutiveLawManager&
Foam::solidModel::mechanicalManager() const
{
    if (mechanicalManagerPtr_.empty())
    {
        Info<< "Creating the mechanicalConstitutiveLawManager" << endl;

        mechanicalManagerPtr_.set
        (
            new mechanicalConstitutiveLawManager(mesh(), mechanicalProperties())
        );

        mechanicalManagerPtr_->setRestartKinematicsAvailable
        (
            restartKinematicsAvailable_
        );
    }

    return mechanicalManagerPtr_();
}


void Foam::solidModel::frameworkInterpolate
(
    const volVectorField& D,
    const volTensorField& gradD,
    pointVectorField& pointD
)
{
#ifdef OPENFOAM_NOT_EXTEND
    enhancedVolPointInterpolation::New(mesh()).interpolate(D, gradD, pointD);
#else
    if (mechanicalManager().nLaws() > 1)
    {
        // The least squares fit below would straddle a material interface,
        // where the displacement is continuous but its gradient jumps, and
        // smear it. Instead each cell's value is extrapolated with its own
        // gradient, which the material-aware leastSquaresS4f scheme keeps to
        // one material, as the other forks do for any number of materials
        volToPoint().interpolate(D, gradD, pointD);
    }
    else
    {
        volToPoint().interpolate(D, pointD);
    }
#endif

    correctPointDisplacement(pointD);
}


Foam::tmp<Foam::volScalarField> Foam::solidModel::initialRho() const
{
    // A registered copy, so that a field built from it takes over the
    // registration. The manager's own cache is not registered, and a plain
    // copy of it would not be either; an unregistered field is not mapped
    // when the mesh topology changes
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mechanicalManager().rho()
        )
    );
}


void Foam::solidModel::gradQuad
(
    const volVectorField& D,
    CompactListList<tensor>& gradDQuad
) const
{
    // The framework has no quadrature-point gradient for more than one
    // material
    if (mechanicalManager().nLaws() > 1)
    {
        FatalErrorInFunction
            << "The face quadrature gradient is not implemented for more than "
            << "one material." << exit(FatalError);
    }

    hofvc::fGrad(D, gradDQuad);
}


void Foam::solidModel::checkFrameworkGradScheme(const word& fieldName) const
{
    if (mechanicalManager().nLaws() < 2)
    {
        return;
    }

    const word gradScheme
    (
#ifdef OPENFOAM_NOT_EXTEND
        mesh().gradScheme("grad(" + fieldName + ')')
#else
        mesh().schemesDict().gradScheme("grad(" + fieldName + ')')
#endif
    );

    if (gradScheme != "leastSquaresS4f")
    {
        FatalErrorInFunction
            << "More than one material needs a material-aware gradient for "
            << "grad(" << fieldName << "), and `" << gradScheme << "` is not "
            << "one." << nl << nl
            << "    One gradient is computed on the whole mesh, which only "
            << "works if the scheme keeps a cell's stencil within its own "
            << "material. Set "
            << "`grad(" << fieldName << ") leastSquaresS4f;` in fvSchemes."
            << abort(FatalError);
    }
}


Foam::tmp<Foam::volScalarField> Foam::solidModel::frameworkImpK
(
    mechanicalConstitutiveLawManager& manager,
    const tangentRequest req
) const
{
    tmp<volScalarField> tImpK
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimForce/dimArea, 0),
            calculatedFvPatchScalarField::typeName
        )
    );

#ifdef OPENFOAM_NOT_EXTEND
    volScalarField& impK = tImpK.ref();
#else
    volScalarField& impK = tImpK();
#endif

    // Evaluate the tangent at zero gradient against a state with no history.
    // impK is formed once and kept, so this gives the same elastic
    // preconditioner on a cold start and on a restart without disturbing the
    // stored constitutive state
    const volTensorField zeroGradD
    (
        IOobject
        (
            "zeroGradD",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("zero", gradD().dimensions(), tensor::zero)
    );

    manager.updateScalarTangent
    (
        zeroGradD,
        zeroGradD,
        mesh().time().deltaTValue(),
        impK,
        req,
        true        // evaluate against a state with no history
    );

    return tImpK;
}


void Foam::solidModel::end()
{
    solidProperties_.IOobject::rename
    (
        solidProperties().IOobject::name() + ".withDefaultValues"
    );
    solidProperties_.regIOobject::write();

    if (!thermalPtr_.empty())
    {
        thermal().IOobject::rename
        (
            thermal().IOobject::name() + ".withDefaultValues"
        );
        static_cast<const IOdictionary>(thermal()).regIOobject::write();
    }

    if (maxIterReached_ > 0)
    {
        WarningIn(type() + "::end()")
            << "The maximum momentum correctors were reached in "
            << maxIterReached_ << " time-steps" << nl << endl;
    }
    else
    {
        Info<< "The momentum equation converged in all time-steps"
            << nl << endl;
    }

    physicsModel::end();
}


Foam::tmp<Foam::volTensorField> Foam::solidModel::couplingField
(
    const word& fieldName
) const
{
    FatalErrorInFunction
        << "Coupling field \"" << fieldName
        << "\" requested from solidModel type " << type() << nl
        << "This solidModel does not support coupling fields." << nl
        << "Override couplingField() in the derived class to enable coupling."
        << abort(FatalError);

    // Keep compiler happy
    return tmp<volTensorField>(nullptr);
}


Foam::autoPtr<Foam::solidModel> Foam::solidModel::New
(
    Time& runTime,
    const word& region
)
{

    // It is possible to run a single region of a multi-region case (e.g., to
    // check the convergence of that single region) by setting it in
    // physicsProperties and in the controlDict under the subDict 'solid'.
    // This also allows the region name not to be 'region0'
    // or 'solid' but user defined.
    // See https://github.com/solids4foam/solids4foam/pull/83
    const word runRegion
    (
        runTime.controlDict().subOrEmptyDict("solid").lookupOrDefault<word>
        (
            "region", region
        )
    );

    // NB: dictionary must be unregistered to avoid adding to the database

    IOdictionary props
    (
        IOobject
        (
            "solidProperties",
            bool(runRegion == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/runRegion),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false  // Do not register
        )
    );

    const word modelType(props.lookup("solidModel"));

    Info<< "Selecting solidModel " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            props,
            "solidModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "solidModel::New(Time&, const word&)"
        )   << "Unknown solidModel type " << modelType
            << endl << endl
            << "Valid solidModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    autoPtr<solidModel> modelPtr(ctorPtr(runTime, runRegion));

    // Asked here, once the model is fully constructed and its override can be
    // seen, so that a solid model which cannot smooth the hydrostatic stress
    // refuses a case that asks for it rather than silently ignoring it.
    //
    // Only of a model that has built the framework manager, and so has read
    // mechanicalProperties already
    if
    (
        modelPtr->mechanicalManagerPtr_.valid()
     && !modelPtr->supportsHydrostaticSmoothing()
     && modelPtr->hydrostaticSmoothingRequested()
    )
    {
        FatalErrorIn("solidModel::New(Time&, const word&)")
            << "solvePressureEqn is set in mechanicalProperties, and solid "
            << "model " << modelType << " cannot smooth the hydrostatic "
            << "stress." << nl
            << "    It is supported by the updated Lagrangian, total "
            << "Lagrangian and linear geometry total displacement models. "
            << "Remove solvePressureEqn, or use one of those."
            << exit(FatalError);
    }

    // And of a model that can, the request is validated here too, rather than
    // at the first stress update: a model running the mixed formulation
    // returns from its stress update before reaching the smoothing, so
    // solvePressure with solvePressureEqn would otherwise be accepted and the
    // smoothing silently ignored
    if
    (
        modelPtr->mechanicalManagerPtr_.valid()
     && modelPtr->supportsHydrostaticSmoothing()
     && modelPtr->hydrostaticSmoothingRequested()
    )
    {
        modelPtr->smoothHydrostaticStress();
    }

    return modelPtr;
}


void Foam::solidModel::setTraction
(
    fvPatchVectorField& tractionPatch,
    const vectorField& traction
)
{
    if (tractionPatch.type() == solidTractionFvPatchVectorField::typeName)
    {
        solidTractionFvPatchVectorField& patchD =
            refCast<solidTractionFvPatchVectorField>(tractionPatch);

        patchD.setTraction(traction);
    }
#ifdef FOAMEXTEND
    else if
    (
        tractionPatch.type() == blockSolidTractionFvPatchVectorField::typeName
    )
    {
        blockSolidTractionFvPatchVectorField& patchD =
            refCast<blockSolidTractionFvPatchVectorField>(tractionPatch);

        patchD.traction() = traction;
    }
#endif
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::setTraction\n"
            "(\n"
            "    fvPatchVectorField& tractionPatch,\n"
            "    const vectorField& traction\n"
            ")"
        )   << "Boundary condition "
            << tractionPatch.type()
            << " for patch " << tractionPatch.patch().name()
            << " should instead be type "
            << solidTractionFvPatchVectorField::typeName
#ifdef FOAMEXTEND
            << " or "
            << blockSolidTractionFvPatchVectorField::typeName
#endif
            << abort(FatalError);
    }
}

void Foam::solidModel::setTractionQuadrature
(
    fvPatchVectorField& tractionPatch,
    const CompactListList<vector>& traction
)
{
    if (tractionPatch.type() == solidTractionFvPatchVectorField::typeName)
    {
        solidTractionFvPatchVectorField& patchD =
            refCast<solidTractionFvPatchVectorField>(tractionPatch);

        patchD.setTractionQuadrature(traction);
    }
    else
    {
        FatalErrorInFunction
            << "Boundary condition " << tractionPatch.type() << " for patch "
            << tractionPatch.patch().name() << " should instead be type "
            << solidTractionFvPatchVectorField::typeName << abort(FatalError);
    }
}


void Foam::solidModel::setTraction
(
    const label interfaceI,
    const label patchID,
    const vectorField& faceZoneTraction
)
{
    const vectorField patchTraction
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneTraction)
    );

    setTraction(boundaryFieldRef(solutionD())[patchID], patchTraction);
}


void Foam::solidModel::setPressure
(
    fvPatchVectorField& tractionPatch,
    const scalarField& pressure
)
{
    if (tractionPatch.type() == solidTractionFvPatchVectorField::typeName)
    {
        solidTractionFvPatchVectorField& patchD =
            refCast<solidTractionFvPatchVectorField>(tractionPatch);

        patchD.pressure() = pressure;
    }
#ifdef FOAMEXTEND
    else if
    (
        tractionPatch.type() == blockSolidTractionFvPatchVectorField::typeName
    )
    {
        blockSolidTractionFvPatchVectorField& patchD =
            refCast<blockSolidTractionFvPatchVectorField>(tractionPatch);

        patchD.pressure() = pressure;
    }
#endif
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::setTraction\n"
            "(\n"
            "    fvPatchVectorField& tractionPatch,\n"
            "    const vectorField& traction\n"
            ")"
        )   << "Boundary condition "
            << tractionPatch.type()
            << " for patch " << tractionPatch.patch().name()
            << " should instead be type "
            << solidTractionFvPatchVectorField::typeName
#ifdef FOAMEXTEND
            << " or "
            << blockSolidTractionFvPatchVectorField::typeName
#endif
            << abort(FatalError);
    }
}


void Foam::solidModel::setPressure
(
    const label interfaceI,
    const label patchID,
    const scalarField& faceZonePressure
)
{
    const scalarField patchPressure
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZonePressure)
    );

    setPressure(boundaryFieldRef(solutionD())[patchID], patchPressure);
}


void Foam::solidModel::recalculateRho()
{
    rhoPtr_.clear();
    makeRho();
}


void Foam::solidModel::clearLeastSquaresData()
{
    gradDQuadPtr_.clear();
    sigmaQuadPtr_.clear();
    leastSquaresReconstruction::New(mesh()).clear();
}


Foam::Switch& Foam::solidModel::checkEnforceLinear(const volScalarField& J)
{
    scalar minJ = min(J).value();
    reduce(minJ, minOp<scalar>());

    scalar maxJ = max(J).value();
    reduce(maxJ, maxOp<scalar>());

    if ((minJ < 0.01) || (maxJ > 100))
    {
        Info<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear() = true;
    }

    return enforceLinear();
}


Foam::Switch& Foam::solidModel::checkEnforceLinear(const surfaceScalarField& J)
{
    scalar minJ = min(J).value();
    reduce(minJ, minOp<scalar>());

    scalar maxJ = max(J).value();
    reduce(maxJ, maxOp<scalar>());

    if ((minJ < 0.01) || (maxJ > 100))
    {
        Info<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear() = true;
    }

    return enforceLinear();
}


void Foam::solidModel::writeFields(const Time& runTime)
{
    // Write strain fields
    // Currently only defined for linear geometry
    if (nonLinGeom() == nonLinearGeometry::LINEAR_GEOMETRY)
    {
        // Total strain
        volSymmTensorField epsilon("epsilon", symm(gradD()));
        epsilon.write();

        // Equivalent strain
        volScalarField epsilonEq
        (
            "epsilonEq", sqrt((2.0/3.0)*magSqr(dev(epsilon)))
        );
        epsilonEq.write();

        Info<< "Max epsilonEq = " << gMax(epsilonEq) << endl;
    }

    // Calculate equivalent (von Mises) stress
    volScalarField sigmaEq
    (
        "sigmaEq", sqrt((3.0/2.0)*magSqr(dev(sigma())))
    );
    sigmaEq.write();

    Info<< "Max sigmaEq (von Mises stress) = " << gMax(sigmaEq) << endl;

    // The constitutive history, such as the plastic strain, as fields that
    // can be viewed. The restart files hold it in a form ParaView cannot read
    if (mechanicalManagerPtr_.valid())
    {
        mechanicalManagerPtr_->writeStateFields();
    }

    // If asked, write the residual field
    if (writeResidualField_)
    {
        const volVectorField& D = solutionD();
#ifdef OPENFOAM_NOT_EXTEND
        scalar denom =
            gMax(mag(D.primitiveField() - D.oldTime().primitiveField()));
        if (denom < SMALL)
        {
            denom = max(gMax(mag(D.primitiveField())), SMALL);
        }
#else
        scalar denom =
            gMax(mag(D.internalField() - D.oldTime().internalField()));
        if (denom < SMALL)
        {
            denom = max(gMax(mag(D.internalField())), SMALL);
        }
#endif

        const volVectorField residualD
        (
            "residualD",
            (D - D.prevIter())/denom
        );

        Info<< "Writing residualD field" << endl;
        residualD.write();
    }

    physicsModel::writeFields(runTime);
}


void Foam::solidModel::moveMesh
(
    const pointField& oldPoints,
    const pointVectorField& pointDD
)
{
    Info<< "Moving the mesh to the deformed configuration" << nl << endl;

    const vectorField& pointDDI = pointDD;

    // Calculate the new points and apply 2-D corrections
    vectorField newPoints(oldPoints + pointDDI);
    twoDCorrector_.correctPoints(newPoints);

    // Move the mesh
    mesh().movePoints(newPoints);
    mesh().V00();
    mesh().moving(false);
#ifdef FOAMEXTEND
    mesh().changing(false);
#endif
#if (OPENFOAM >= 2206)
    {
        auto tmeshPhi(mesh().setPhi());
        if (tmeshPhi)
        {
            tmeshPhi.ref().writeOpt(IOobject::NO_WRITE);
        }
    }
#else
    mesh().setPhi().writeOpt() = IOobject::NO_WRITE;
#endif
}


const Foam::dictionary& Foam::solidModel::solidModelDict() const
{
    return solidProperties_.subDict(type_ + "Coeffs");
}


Foam::tangentRequest Foam::solidModel::jacobianTangent
(
    const tangentRequest deflt
) const
{
    if (jacobianTangentCached_)
    {
        return jacobianTangent_;
    }

    const dictionary& dict = solidModelDict();

    if (dict.found("approximateJacobian"))
    {
        if (dict.found("jacobianTangent"))
        {
            FatalIOErrorInFunction(dict)
                << "Both 'approximateJacobian' and 'jacobianTangent' are set."
                << nl
                << "'approximateJacobian' is deprecated: use "
                << "'jacobianTangent' only."
                << exit(FatalIOError);
        }

        const Switch approximate(dict.lookup("approximateJacobian"));

        jacobianTangent_ =
            approximate
          ? tangentRequest::scalar
          : tangentRequest::fourthOrder;

        WarningInFunction
            << "'approximateJacobian' is deprecated. Replace it with "
            << "'jacobianTangent " << tangentRequestName(jacobianTangent_)
            << ";'" << endl;
    }
    else if (dict.found("jacobianTangent"))
    {
        jacobianTangent_ =
            tangentRequestNamed(word(dict.lookup("jacobianTangent")));
    }
    else
    {
        jacobianTangent_ = deflt;
    }

    jacobianTangentCached_ = true;

    return jacobianTangent_;
}

// ************************************************************************* //
