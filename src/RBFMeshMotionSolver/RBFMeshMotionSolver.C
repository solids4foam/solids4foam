
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFMeshMotionSolver.H"
#include <unordered_map>
#include <memory>
#include <vector>
#include "FieldSumOp.H"

using namespace Foam;

defineTypeNameAndDebug(RBFMeshMotionSolver, 0);

addToRunTimeSelectionTable
(
    motionSolver,
    RBFMeshMotionSolver,
    dictionary
);

RBFMeshMotionSolver::RBFMeshMotionSolver
(
    const polyMesh& mesh,
#if defined(OPENFOAM_ORG)
    const dictionary& msData
#elif defined(OPENFOAM_COM)
    const IOdictionary& msData
#else
    Istream & msData
#endif
)
    :
#ifdef OPENFOAM_ORG
    motionSolver(mesh, msData, typeName),
#else
    motionSolver(mesh),
#endif
    motionCenters(mesh.boundaryMesh().size(), vectorField(0)),
    accumulatedMotionCentersField_
    (
        IOobject
        (
            "rbfMotionCentersField",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh.parent().lookupObject<fvMesh>(mesh.name()),
        dimensionedVector("0", dimLength, vector::zero)
    ),
#ifdef OPENFOAM_NOT_EXTEND
    staticPatches(msData.subDict("RBFMeshMotionSolverCoeffs").lookup("staticPatches")),
    staticPatchIDs(staticPatches.size()),
    movingPatches(msData.subDict("RBFMeshMotionSolverCoeffs").lookup("movingPatches")),
    movingPatchIDs(movingPatches.size()),
    fixedPatches(msData.subDict("RBFMeshMotionSolverCoeffs").lookup("fixedPatches")),
    fixedPatchIDs(fixedPatches.size()),
#else
    staticPatches(lookup("staticPatches")),
    staticPatchIDs(staticPatches.size()),
    movingPatches(lookup("movingPatches")),
    movingPatchIDs(movingPatches.size()),
    fixedPatches(lookup("fixedPatches")),
    fixedPatchIDs(fixedPatches.size()),
#endif
    newPoints(mesh.points().size(), vector::zero),
    rbf(),
    nbGlobalFaceCenters(Pstream::nProcs(), 0),
    nbGlobalMovingFaceCenters(Pstream::nProcs(), 0),
    nbGlobalStaticFaceCenters(Pstream::nProcs(), 0),
    nbGlobalFixedFaceCenters(Pstream::nProcs(), 0),
    globalMovingPointsLabelList(mesh.boundaryMesh().size(), labelList(0)),
    twoDCorrector(mesh),
    nbPoints(0),
    faceCellCenters(true),
    cpu(false),
    corrector(false),
    k(0)
{
#ifdef OPENFOAM_NOT_EXTEND
    const dictionary& coeffDict = msData.subDict("RBFMeshMotionSolverCoeffs");
#else
    const dictionary& coeffDict = *this;
#endif

    // Create previous iteration field
    accumulatedMotionCentersField_.storePrevIter();

    // Find IDs of staticPatches
    forAll(staticPatches, patchI)
    {
        const label patchIndex = mesh.boundaryMesh().findPatchID(staticPatches[patchI]);

        assert(patchIndex >= 0);

        staticPatchIDs[patchI] = patchIndex;
    }

    // Find IDs of movingPatches
    forAll(movingPatches, patchI)
    {
        const label patchIndex = mesh.boundaryMesh().findPatchID(movingPatches[patchI]);

        assert(patchIndex >= 0);

        movingPatchIDs[patchI] = patchIndex;
    }

    // Find IDs of fixedPatches
    forAll(fixedPatches, patchI)
    {
        const label patchIndex = mesh.boundaryMesh().findPatchID(fixedPatches[patchI]);

        assert(patchIndex >= 0);

        fixedPatchIDs[patchI] = patchIndex;
    }

    // Verify that a patch is not defined as a static and a moving patch

    forAll(staticPatchIDs, staticPatchI)
    {
        // Search the moving patches for static patchI
        forAll(movingPatchIDs, movingPatchI)
        {
            assert(movingPatchIDs[movingPatchI] != staticPatchIDs[staticPatchI]);
        }

        // Search the fixed patches for static patchI
        forAll(fixedPatchIDs, fixedPatchI)
        {
            assert(fixedPatchIDs[fixedPatchI] != staticPatchIDs[staticPatchI]);
        }
    }

    forAll(fixedPatchIDs, fixedPatchI)
    {
        // Search the moving patches for fixed patchI
        forAll(movingPatchIDs, movingPatchI)
        {
            assert(movingPatchIDs[movingPatchI] != fixedPatchIDs[fixedPatchI]);
        }
    }

    // Initialize RBF interpolator

    const dictionary & dict = coeffDict.subDict("interpolation");

    word function = word(dict.lookup("function"));

    assert
    (
        function == "TPS"
     || function == "WendlandC0"
     || function == "WendlandC2"
     || function == "WendlandC4"
     || function == "WendlandC6"
   );

    std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction;

    Info << "Radial Basis Function interpolation: Selecting RBF function: " << function << endl;

    if (function == "TPS")
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> (new rbf::TPSFunction());

    if (function == "WendlandC0")
    {
        scalar radius = readScalar(dict.lookup("radius"));
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> (new rbf::WendlandC0Function(radius));
    }

    if (function == "WendlandC2")
    {
        scalar radius = readScalar(dict.lookup("radius"));
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> (new rbf::WendlandC2Function(radius));
    }

    if (function == "WendlandC4")
    {
        scalar radius = readScalar(dict.lookup("radius"));
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> (new rbf::WendlandC4Function(radius));
    }

    if (function == "WendlandC6")
    {
        scalar radius = readScalar(dict.lookup("radius"));
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> (new rbf::WendlandC6Function(radius));
    }

    assert(rbfFunction);

    bool polynomialTerm = dict.lookupOrDefault("polynomial", false);
    bool cpu = dict.lookupOrDefault("cpu", false);
    this->cpu = dict.lookupOrDefault("fullCPU", false);
    std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator(new rbf::RBFInterpolation(rbfFunction, polynomialTerm, cpu));

    if (this->cpu == true)
        assert(cpu == true);

    const bool coarsening =
        readBool(coeffDict.subDict("coarsening").lookup("enabled"));
    scalar tol = 0.1;
    scalar tolLivePointSelection = 0.1;
    bool livePointSelection = false;
    bool exportSelectedPoints = false;
    int coarseningMinPoints = 1;
    int coarseningMaxPoints = 2;
    bool twoPointSelection = false;
    bool surfaceCorrection = false;
    scalar ratioRadiusError = 10.0;

    if (coarsening)
    {
        tol = readScalar(coeffDict.subDict("coarsening").lookup("tol"));
        coarseningMinPoints =
            readLabel(coeffDict.subDict("coarsening").lookup("minPoints"));
        coarseningMaxPoints =
            readLabel(coeffDict.subDict("coarsening").lookup("maxPoints"));
        livePointSelection =
            readBool
            (
                coeffDict.subDict("coarsening").lookup("livePointSelection")
            );
        exportSelectedPoints =
            readBool
            (
                coeffDict.subDict("coarsening").lookup("exportSelectedPoints")
            );
        twoPointSelection =
            coeffDict.subDict("coarsening").lookupOrDefault
            (
                "twoPointSelection", false
            );
    }

    if (livePointSelection)
    {
        tolLivePointSelection = readScalar
        (
            coeffDict.subDict("coarsening").lookup("tolLivePointSelection")
        );
        surfaceCorrection =
            coeffDict.subDict("coarsening").lookupOrDefault
            (
                "surfaceCorrection", false
            );

        if (surfaceCorrection)
        {
            ratioRadiusError =
                coeffDict.subDict("coarsening").lookupOrDefault
                (
                    "ratioRadiusError", 10.0
                );
        }
    }

    rbf = std::shared_ptr<rbf::RBFCoarsening>
    (
        new rbf::RBFCoarsening
        (
            rbfInterpolator, coarsening, livePointSelection, true, tol,
            tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints,
            twoPointSelection, surfaceCorrection, ratioRadiusError,
            exportSelectedPoints
        )
    );

    faceCellCenters = coeffDict.lookupOrDefault("faceCellCenters", true);

    Info << "RBF mesh deformation settings:" << endl;
    Info << "    interpolation function = " << function << endl;
    Info << "    interpolation polynomial term = " << polynomialTerm << endl;
    Info << "    interpolation cpu formulation = " << cpu << endl;
    Info << "    coarsening = " << coarsening << endl;
    Info << "        coarsening tolerance = " << tol << endl;
    Info << "        coarsening reselection tolerance = " << tolLivePointSelection << endl;
    Info << "        coarsening two-point selection = " << twoPointSelection << endl;
}

RBFMeshMotionSolver::~RBFMeshMotionSolver()
{}

tmp<pointField> RBFMeshMotionSolver::curPoints() const
{
    // Prepare new points: same as old point
    tmp<pointField> tnewPoints
    (
        new vectorField(mesh().nPoints(), vector::zero)
   );

#ifdef OPENFOAM_NOT_EXTEND
    pointField & newPoints = tnewPoints.ref();
#else
    pointField & newPoints = tnewPoints();
#endif

    newPoints = this->newPoints;

    // Add old point positions
    newPoints += mesh().points();

    return tnewPoints;
}

// As a first step, the motion is defined in the
void RBFMeshMotionSolver::setMotion(const Field<vectorField> & motion)
{
    // Input checking

    assert(motion.size() == mesh().boundaryMesh().size());

    forAll(motion, ipatch)
    {
        const vectorField & mpatch = motion[ipatch];

        // Check whether the size of patch motion is equal to number of face centers in patch
        if (faceCellCenters && mpatch.size() > 0)
            assert(mpatch.size() == mesh().boundaryMesh()[ipatch].faceCentres().size());

        if (not faceCellCenters && mpatch.size() > 0)
            assert(mpatch.size() == mesh().boundaryMesh()[ipatch].meshPoints().size());

        // Check whether the size of a moving patch is equal to the number of face centers in the patch
        // First check if patchid is a moving patch
        bool movingPatch = false;
        forAll(movingPatchIDs, movingPatchI)
        {
            if (movingPatchIDs[movingPatchI] == ipatch)
                movingPatch = true;
        }

        if (faceCellCenters && movingPatch)
            assert(mpatch.size() == mesh().boundaryMesh()[ipatch].faceCentres().size());

        if (not faceCellCenters && movingPatch)
            assert(mpatch.size() == mesh().boundaryMesh()[ipatch].meshPoints().size());

        // Set values on motionCentersField boundary
        if (mesh().boundaryMesh()[ipatch].type() == "empty")
        {
            // Store data in internal field as field empty patches have zero
            // size
        #ifdef OPENFOAM_NOT_EXTEND
            const labelUList& faceCells = mesh().boundaryMesh()[ipatch].faceCells();
            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];
                accumulatedMotionCentersField_[cellID] = motion[ipatch][faceI];
            }
        #else
            const labelList& faceCells = mesh().boundaryMesh()[ipatch].faceCells();
            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];
                accumulatedMotionCentersField_[cellID] = motion[ipatch][faceI];
            }
        #endif
        }
        else
        {
        #ifdef OPENFOAM_NOT_EXTEND
            accumulatedMotionCentersField_.boundaryFieldRef()[ipatch] = motion[ipatch];
        #else
            accumulatedMotionCentersField_.boundaryField()[ipatch] = motion[ipatch];
        #endif
        }
    }

    motionCenters = motion;
}

void RBFMeshMotionSolver::updateMesh(const mapPolyMesh &)
{
    assert(false);
}

void RBFMeshMotionSolver::solve()
{
    assert(motionCenters.size() == mesh().boundaryMesh().size());

    // Copy motionCentersField to motionCenters
    forAll(accumulatedMotionCentersField_.boundaryField(), patchI)
    {
        if (mesh().boundaryMesh()[patchI].type() == "empty")
        {
        #ifdef OPENFOAM_NOT_EXTEND
            const labelUList& faceCells = mesh().boundaryMesh()[patchI].faceCells();
            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];
                motionCenters[patchI][faceI] =
                    accumulatedMotionCentersField_[cellID]
                  - accumulatedMotionCentersField_.prevIter()[cellID];
            }
        #else
            const labelList& faceCells = mesh().boundaryMesh()[patchI].faceCells();
            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];
                motionCenters[patchI][faceI] =
                    accumulatedMotionCentersField_[cellID]
                  - accumulatedMotionCentersField_.prevIter()[cellID];
            }
        #endif
        }
        else
        {
            motionCenters[patchI] =
                accumulatedMotionCentersField_.boundaryField()[patchI]
              - accumulatedMotionCentersField_.prevIter().boundaryField()[patchI];
        }
    }

    // Update previous field
    accumulatedMotionCentersField_.storePrevIter();


    /*
     * RBF interpolator from face centers to local complete mesh vertices
     * The interpolation consists of the following steps:
     * 1. Build a matrix with the face center positions of the static patches and the moving patches
     * 2. Build a matrix with the positions of every vertex in the local mesh
     * 3. Build a matrix with the displacement/motion of the face center positions of the static patches and the moving patches
     * 4. Perform the interpolation from the face centers to the complete mesh
     * 5. Correct the mesh vertices of the static patches. Set these displacement to zero.
     * 6. Set the motion of the mesh vertices
     */

    /*
     * Step 1: Build a matrix with the face center positions of the static patches and the moving patches
     * The order of the matrix is defined as first a list of the moving patch face centers,
     * thereafter the static patch face centers. These are the control points used by the
     * radial basis function interpolation.
     * The control points should be exactly the same at each processor, and are therefore communicated
     * to each process. As only an absolute RBF interpolation is implemented, this communication step is only
     * performed once per simulation.
     * The global ordering of the data is first the information of the moving patches,
     * thereafter the static patches.
     */

    unsigned int nbFaceCenters = 0;
    unsigned int nbMovingFaceCenters = 0;
    unsigned int nbStaticFaceCenters = 0;
    unsigned int nbFixedFaceCenters = 0;

    std::unordered_map<unsigned int, unsigned int> staticControlPointLabels;
    std::unordered_map<unsigned int, unsigned int> fixedControlPointLabels;
    std::vector<unsigned int> movingControlPointLabelsVector;
    std::unordered_map<unsigned int, unsigned int> movingControlPointLabelsMap;

    std::vector<unsigned int> movingControlPointPatchIds;
    std::vector<unsigned int> movingControlPointIndices;

    std::unordered_map<unsigned int, unsigned int> staticControlGlobalPointLabels;
    std::unordered_map<unsigned int, unsigned int> fixedControlGlobalPointLabels;
    std::unordered_map<unsigned int, unsigned int> movingControlGlobalPointLabelsMap;

    labelList globalStaticPointsListEnabled(nbStaticFaceCenters, 0);
    labelList globalFixedPointsListEnabled(nbFixedFaceCenters, 0);
    labelList globalMovingPointsListEnabled(nbMovingFaceCenters, 0);
    unsigned int globalStaticOffsetNonUnique = 0;
    unsigned int globalFixedOffsetNonUnique = 0;
    unsigned int globalMovingOffsetNonUnique = 0;

    if (sum(nbGlobalFaceCenters) == 0)
    {
        // Determine the number of face centers
        // The total number of face centers is simply the sum of the face centers
        // on each processor.
        nbFaceCenters = 0;

        // First add the static patches, thereafter the fixed patches, and
        // the moving patches as last.

        forAll(staticPatchIDs, i)
        {
            const labelList & meshPoints = mesh().boundaryMesh()[staticPatchIDs[i]].meshPoints();

            forAll(meshPoints, j)
            {
                if (twoDCorrector.marker()[meshPoints[j]] != 0)
                {
                    continue;
                }

                if (staticControlPointLabels.find(meshPoints[j]) == staticControlPointLabels.end())
                {
                    // staticControlPointLabels[meshPoints[j]] = staticControlPointLabels.size() - 1;
                    staticControlPointLabels[meshPoints[j]] = staticControlPointLabels.size();
                }
            }
        }

        nbStaticFaceCenters = staticControlPointLabels.size();

        forAll(fixedPatchIDs, i)
        {
            const labelList & meshPoints = mesh().boundaryMesh()[fixedPatchIDs[i]].meshPoints();

            forAll(meshPoints, j)
            {
                if (twoDCorrector.marker()[meshPoints[j]] != 0)
                    continue;

                if (staticControlPointLabels.find(meshPoints[j]) == staticControlPointLabels.end()
                    && fixedControlPointLabels.find(meshPoints[j]) == fixedControlPointLabels.end())
                    fixedControlPointLabels[meshPoints[j]] = fixedControlPointLabels.size();
                    // fixedControlPointLabels[meshPoints[j]] = fixedControlPointLabels.size() - 1;
            }
        }

        nbFixedFaceCenters = fixedControlPointLabels.size();

        if (faceCellCenters)
        {
            forAll(movingPatchIDs, i)
            {
                nbMovingFaceCenters += mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres().size();
            }
        }

        if (not faceCellCenters)
        {
            forAll(movingPatchIDs, patchI)
            {
                const labelList & meshPoints = mesh().boundaryMesh()[movingPatchIDs[patchI]].meshPoints();
                globalMovingPointsLabelList[movingPatchIDs[patchI]] = labelList(meshPoints.size(), 0);

                forAll(meshPoints, j)
                {
                    if (twoDCorrector.marker()[meshPoints[j]] != 0)
                        continue;

                    if (staticControlPointLabels.find(meshPoints[j]) == staticControlPointLabels.end()
                        && fixedControlPointLabels.find(meshPoints[j]) == fixedControlPointLabels.end()
                        && movingControlPointLabelsMap.find(meshPoints[j]) == movingControlPointLabelsMap.end())
                    {
                        movingControlPointLabelsMap[meshPoints[j]] = movingControlPointLabelsMap.size() - 1;
                        movingControlPointLabelsVector.push_back(meshPoints[j]);
                        movingControlPointPatchIds.push_back(movingPatchIDs[patchI]);
                        movingControlPointIndices.push_back(j);
                        globalMovingPointsLabelList[movingPatchIDs[patchI]][j] = 1;
                    }
                }
            }

            nbMovingFaceCenters = movingControlPointLabelsVector.size();
        }

        if (Pstream::nProcs() == 1)
        {
            globalStaticPointsListEnabled.resize(nbStaticFaceCenters);
            globalStaticPointsListEnabled = 1;
            globalFixedPointsListEnabled.resize(nbFixedFaceCenters);
            globalFixedPointsListEnabled = 1;
            globalMovingPointsListEnabled.resize(nbMovingFaceCenters);
            globalMovingPointsListEnabled = 1;
        }

        if (Pstream::nProcs() > 1)
        {
            IOobject addrHeader
            (
                "pointProcAddressing",
                mesh().facesInstance(),
                mesh().meshSubDir,
                mesh(),
                IOobject::MUST_READ
           );

#ifdef OPENFOAM_NOT_EXTEND
            assert(addrHeader.typeHeaderOk<labelIOList>(true));
#else
            assert(addrHeader.headerOk());
#endif
            labelIOList pointProcAddressing(addrHeader);

            assert(pointProcAddressing.size() == mesh().points().size());

            // Count the number of global static points including scalar points
            nbGlobalStaticFaceCenters[Pstream::myProcNo()] = nbStaticFaceCenters;
            nbGlobalFixedFaceCenters[Pstream::myProcNo()] = nbFixedFaceCenters;
            nbGlobalMovingFaceCenters[Pstream::myProcNo()] = nbMovingFaceCenters;
            reduce(nbGlobalStaticFaceCenters, sumOp<labelList>());
            reduce(nbGlobalFixedFaceCenters, sumOp<labelList>());
            reduce(nbGlobalMovingFaceCenters, sumOp<labelList>());
            nbStaticFaceCenters = sum(nbGlobalStaticFaceCenters);
            nbFixedFaceCenters = sum(nbGlobalFixedFaceCenters);

            if (not faceCellCenters)
                nbMovingFaceCenters = sum(nbGlobalMovingFaceCenters);

            // Construct a list with all the global point labels, thus including
            // also scalar points. Thereafter, construct a list of static control
            // list which indicates whether the point is already included or not.
            // Use this later to build a list of the unique static control points.
            labelList globalStaticPointsList(nbStaticFaceCenters, 0);
            labelList globalFixedPointsList(nbFixedFaceCenters, 0);
            labelList globalMovingPointsList(nbMovingFaceCenters, 0);
            labelList globalMovingPointsPatchIds(nbMovingFaceCenters, 0);
            labelList globalMovingPointsIndices(nbMovingFaceCenters, 0);

            globalStaticOffsetNonUnique = 0;

            for (int i = 0; i < Pstream::myProcNo(); i++)
                globalStaticOffsetNonUnique += nbGlobalStaticFaceCenters[i];

            globalFixedOffsetNonUnique = 0;

            for (int i = 0; i < Pstream::myProcNo(); i++)
                globalFixedOffsetNonUnique += nbGlobalFixedFaceCenters[i];

            globalMovingOffsetNonUnique = 0;

            for (int i = 0; i < Pstream::myProcNo(); i++)
                globalMovingOffsetNonUnique += nbGlobalMovingFaceCenters[i];

            for (auto label : staticControlPointLabels)
            {
                globalStaticPointsList[label.second + globalStaticOffsetNonUnique] = pointProcAddressing[label.first];
            }

            for (auto label : fixedControlPointLabels)
            {
                globalFixedPointsList[label.second + globalFixedOffsetNonUnique] = pointProcAddressing[label.first];
            }

            for (unsigned int i = 0; i < movingControlPointLabelsVector.size(); ++i)
            {
                globalMovingPointsList[i + globalMovingOffsetNonUnique] = pointProcAddressing[movingControlPointLabelsVector[i]];
                globalMovingPointsPatchIds[i + globalMovingOffsetNonUnique] = movingControlPointPatchIds[i];
                globalMovingPointsIndices[i + globalMovingOffsetNonUnique] = movingControlPointIndices[i];
            }

            reduce(globalStaticPointsList, sumOp<labelList>());
            reduce(globalFixedPointsList, sumOp<labelList>());

            if (not faceCellCenters)
            {
                reduce(globalMovingPointsList, sumOp<labelList>());
                reduce(globalMovingPointsPatchIds, sumOp<labelList>());
                reduce(globalMovingPointsIndices, sumOp<labelList>());
            }

            // Construct a list of static control points which indicate whether
            // should be included or not.

            globalStaticPointsListEnabled.resize(nbStaticFaceCenters);
            globalStaticPointsListEnabled = 0;
            globalFixedPointsListEnabled.resize(nbFixedFaceCenters);
            globalFixedPointsListEnabled = 0;
            globalMovingPointsListEnabled.resize(nbMovingFaceCenters);
            globalMovingPointsListEnabled = 0;
            forAll(globalStaticPointsList, i)
            {
                if (staticControlGlobalPointLabels.find(globalStaticPointsList[i]) == staticControlGlobalPointLabels.end())
                {
                    staticControlGlobalPointLabels[globalStaticPointsList[i]] = staticControlGlobalPointLabels.size() - 1;
                    globalStaticPointsListEnabled[i] = 1;
                }
            }

            forAll(globalFixedPointsList, i)
            {
                if (staticControlGlobalPointLabels.find(globalFixedPointsList[i]) == staticControlGlobalPointLabels.end()
                    && fixedControlGlobalPointLabels.find(globalFixedPointsList[i]) == fixedControlGlobalPointLabels.end())
                {
                    fixedControlGlobalPointLabels[globalFixedPointsList[i]] = fixedControlGlobalPointLabels.size() - 1;
                    globalFixedPointsListEnabled[i] = 1;
                }
            }

            if (not faceCellCenters)
            {
                forAll(movingPatchIDs, patchI)
                {
                    const labelList & meshPoints = mesh().boundaryMesh()[movingPatchIDs[patchI]].meshPoints();
                    globalMovingPointsLabelList[movingPatchIDs[patchI]] = labelList(meshPoints.size(), 0);
                }

                forAll(globalMovingPointsList, i)
                {
                    if (staticControlGlobalPointLabels.find(globalMovingPointsList[i]) == staticControlGlobalPointLabels.end()
                        && fixedControlGlobalPointLabels.find(globalMovingPointsList[i]) == fixedControlGlobalPointLabels.end()
                        && movingControlGlobalPointLabelsMap.find(globalMovingPointsList[i]) == movingControlGlobalPointLabelsMap.end())
                    {
                        movingControlGlobalPointLabelsMap[globalMovingPointsList[i]] = movingControlGlobalPointLabelsMap.size() - 1;
                        globalMovingPointsListEnabled[i] = 1;

                        if (static_cast<unsigned int>(i) < movingControlPointLabelsVector.size() + globalMovingOffsetNonUnique
                            && static_cast<unsigned int>(i) >= globalMovingOffsetNonUnique)
                        {
                            label patchId = globalMovingPointsPatchIds[i];
                            label index = globalMovingPointsIndices[i];
                            globalMovingPointsLabelList[patchId][index] = 1;
                        }
                    }
                }
            }

            // Count the number of local unique static points
            nbStaticFaceCenters = 0;

            for (auto label : staticControlPointLabels)
            {
                if (globalStaticPointsListEnabled[label.second + globalStaticOffsetNonUnique] == 1)
                    nbStaticFaceCenters++;
            }

            nbFixedFaceCenters = 0;

            for (auto label : fixedControlPointLabels)
            {
                if (globalFixedPointsListEnabled[label.second + globalFixedOffsetNonUnique] == 1)
                    nbFixedFaceCenters++;
            }

            if (not faceCellCenters)
            {
                nbMovingFaceCenters = 0;

                for (unsigned int i = 0; i < movingControlPointLabelsVector.size(); ++i)
                {
                    if (globalMovingPointsListEnabled[i + globalMovingOffsetNonUnique] == 1)
                        nbMovingFaceCenters++;
                }
            }
        }

        // Calculate sum of all faces on each processor
        nbGlobalStaticFaceCenters = 0;
        nbGlobalFixedFaceCenters = 0;
        nbGlobalMovingFaceCenters = 0;
        nbGlobalMovingFaceCenters[Pstream::myProcNo()] = nbMovingFaceCenters;
        nbGlobalStaticFaceCenters[Pstream::myProcNo()] = nbStaticFaceCenters;
        nbGlobalFixedFaceCenters[Pstream::myProcNo()] = nbFixedFaceCenters;
        nbGlobalFaceCenters[Pstream::myProcNo()] = nbMovingFaceCenters + nbStaticFaceCenters + nbFixedFaceCenters;

        reduce(nbGlobalMovingFaceCenters, sumOp<labelList>());
        reduce(nbGlobalStaticFaceCenters, sumOp<labelList>());
        reduce(nbGlobalFixedFaceCenters, sumOp<labelList>());
        reduce(nbGlobalFaceCenters, sumOp<labelList>());
    }

    nbMovingFaceCenters = sum(nbGlobalMovingFaceCenters);
    nbStaticFaceCenters = sum(nbGlobalStaticFaceCenters);
    nbFixedFaceCenters = sum(nbGlobalFixedFaceCenters);
    nbFaceCenters = sum(nbGlobalFaceCenters);


    // Determine the offset taking into account multiple processors

    int globalMovingOffset = 0;

    for (int i = 0; i < Pstream::myProcNo(); i++)
        globalMovingOffset += nbGlobalMovingFaceCenters[i];

    int globalStaticOffset = nbMovingFaceCenters;

    for (int i = 0; i < Pstream::myProcNo(); i++)
        globalStaticOffset += nbGlobalStaticFaceCenters[i];

    int globalFixedOffset = nbMovingFaceCenters + nbStaticFaceCenters;

    for (int i = 0; i < Pstream::myProcNo(); i++)
        globalFixedOffset += nbGlobalFixedFaceCenters[i];

    if (!rbf->rbf->computed)
    {
        rbf::matrix positions(nbFaceCenters, mesh().nGeometricD());
        positions.setZero();

        const Foam::pointField & points = mesh().points();

        vectorField positionsField(positions.rows(), vector::zero);

        if (faceCellCenters)
        {
            int offset = 0;

            forAll(movingPatchIDs, i)
            {
                const Foam::vectorField::subField faceCentres = mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres();

                // Set the positions for patch i
                forAll(faceCentres, j)
                {
                    positionsField[j + offset + globalMovingOffset] = faceCentres[j];
                }

                offset += faceCentres.size();
            }
        }

        int index = 0;

        if (not faceCellCenters)
        {
            for (unsigned int i = 0; i < movingControlPointLabelsVector.size(); ++i)
            {
                if (globalMovingPointsListEnabled[i + globalMovingOffsetNonUnique] == 1)
                {
                    assert(index + globalMovingOffset < positionsField.size());
                    positionsField[index + globalMovingOffset] = points[movingControlPointLabelsVector[i]];
                    index++;
                }
            }

            assert(index == nbGlobalMovingFaceCenters[Pstream::myProcNo()]);
        }

        index = 0;

        for (auto label : staticControlPointLabels)
        {
            if (globalStaticPointsListEnabled[label.second + globalStaticOffsetNonUnique] == 1)
            {
                positionsField[index + globalStaticOffset] = points[label.first];
                index++;
            }
        }

        assert(index == nbGlobalStaticFaceCenters[Pstream::myProcNo()]);

        index = 0;

        for (auto label : fixedControlPointLabels)
        {
            if (globalFixedPointsListEnabled[label.second + globalFixedOffsetNonUnique] == 1)
            {
                assert(index + globalFixedOffset < positionsField.size());
                positionsField[index + globalFixedOffset] = points[label.first];
                index++;
            }
        }

        assert(index == nbGlobalFixedFaceCenters[Pstream::myProcNo()]);

        reduce(positionsField, FieldSumOp<vector>());

        // Copy the FOAM vector field to an Eigen matrix
        for (int i = 0; i < positions.rows(); i++)
            for (int j = 0; j < positions.cols(); j++)
                positions(i, j) = positionsField[i][j];

        /*
         * Step 2: Build a matrix with the positions of every vertex in the local mesh.
         * This is only local information and does not need to be communicated to other
         * processors.
         */

        // Determine the number of points by using the 2d corrector
        nbPoints = 0;
        forAll(points, i)
        {
            if (twoDCorrector.marker()[i] == 0)
                nbPoints++;
        }

        rbf::matrix positionsInterpolation(nbPoints, positions.cols());

        index = 0;
        forAll(points, i)
        {
            if (twoDCorrector.marker()[i] == 0)
            {
                for (int j = 0; j < positionsInterpolation.cols(); j++)
                    positionsInterpolation(index, j) = points[i][j];

                index++;
            }
        }

        rbf->compute(positions, positionsInterpolation);

        rbf->setNbMovingAndStaticFaceCenters(nbMovingFaceCenters, nbStaticFaceCenters + nbFixedFaceCenters);
    }

    /*
     * Step 3: Build a matrix with the displacement/motion of the face center
     * positions of the static patches and the moving patches.
     * The motion needs to be communicated to every process at every mesh deformation.
     * This is considered to be the most expensive step with regards to parallel
     * scalability of the overall algorithm.
     */

    rbf::matrix values(nbFaceCenters, mesh().nGeometricD());
    values.setZero();

    vectorField valuesField(values.rows(), vector::zero);

    if (faceCellCenters)
    {
        int offset = 0;

        forAll(movingPatchIDs, i)
        {
            const Foam::vectorField::subField faceCentres = mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres();

            forAll(motionCenters[movingPatchIDs[i]], j)
            {
                valuesField[j + offset + globalMovingOffset] = motionCenters[movingPatchIDs[i]][j];
            }

            offset += faceCentres.size();
        }
    }

    if (not faceCellCenters)
    {
        int index = 0;

        forAll(movingPatchIDs, patchI)
        {
            forAll(globalMovingPointsLabelList[movingPatchIDs[patchI]], j)
            {
                if (globalMovingPointsLabelList[movingPatchIDs[patchI]][j] == 1)
                {
                    valuesField[index + globalMovingOffset] = motionCenters[movingPatchIDs[patchI]][j];
                    index++;
                }
            }
        }

        assert(index == nbGlobalMovingFaceCenters[Pstream::myProcNo()]);
    }

    reduce(valuesField, FieldSumOp<vector>());

    // Copy the FOAM vector field to an Eigen matrix
    for (int i = 0; i < values.rows(); i++)
        for (int j = 0; j < values.cols(); j++)
            values(i, j) = valuesField[i][j];

    /*
     * Step 4: Perform the interpolation from the face centers to the complete mesh
     */

    rbf::matrix valuesInterpolation(nbPoints, values.cols());
    valuesInterpolation.setZero();

    if (cpu)
        rbf->rbf->computed = false;

    rbf->interpolate(values, valuesInterpolation);

    // Apply the 2d correction

    vectorField valuesInterpolationField(mesh().points().size(), Foam::vector::zero);
    int index = 0;
    forAll(valuesInterpolationField, i)
    {
        if (twoDCorrector.marker()[i] == 0)
        {
            for (int j = 0; j < valuesInterpolation.cols(); j++)
                valuesInterpolationField[i][j] = valuesInterpolation(index, j);

            index++;
        }
    }

    twoDCorrector.setShadowSide(valuesInterpolationField);

    /*
     * Step 5: Correct the mesh vertices of the fixed patches. Set these displacements to zero.
     */

    // Loop over all the patches, and set the fixed patches to zero.

    forAll(mesh().boundaryMesh(), i)
    {
        const labelList & meshPoints = mesh().boundaryMesh()[i].meshPoints();

        bool isFixedPatch = false;
        forAll(fixedPatchIDs, j)
        {
            if (i == fixedPatchIDs[j])
                isFixedPatch = true;
        }

        if (isFixedPatch)
        {
            for (int j = 0; j < meshPoints.size(); j++)
                valuesInterpolationField[meshPoints[j]] = Foam::vector::zero;
        }
    }

    /*
     * Step 6: Set the motion of the mesh vertices
     */

    assert(newPoints.size() == valuesInterpolationField.size());

    newPoints = valuesInterpolationField;
}
