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

#include "mechanicalConstitutiveLawManager.H"
#include "compatibilityFunctions.H"
#include "integrationPointTopologies.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mechanicalConstitutiveLawManager, 0);
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const Foam::integrationPointTopology&
Foam::mechanicalConstitutiveLawManager::topologyFor
(
    const word& topologyTypeName
) const
{
    // Already constructed?
    if (topologyCache_.found(topologyTypeName))
    {
        return topology(*topologyCache_[topologyTypeName]).topology_;
    }

    // Lazily construct via OpenFOAM runtime selection
    autoPtr<integrationPointTopology> topoPtr
    (
        integrationPointTopology::New(topologyTypeName, mesh_)
    );

    if (!topoPtr.valid())
    {
        FatalErrorInFunction
            << "Failed to construct integrationPointTopology of type "
            << topologyTypeName
            << exit(FatalError);
    }

    // Cache and return
    topologyCache_.insert(topologyTypeName, topoPtr);

    return topology(*topologyCache_[topologyTypeName]).topology_;
}


const Foam::integrationPointTopology&
Foam::mechanicalConstitutiveLawManager::compactCellTopologyFor
(
    const CompactListList<tensor>& layout
) const
{
    // Unique key per layout instance
    const word key =
        "compactCell:" + Foam::name(reinterpret_cast<std::uintptr_t>(&layout));

    // Already constructed?
    if (topologyCache_.found(key))
    {
        return topology(*topologyCache_[key]).topology_;
    }

    // Construct topology lazily

    // We know this is cell-based compact storage:
    //  - one sub-list per cell
    //  - integration-point counts encoded in sub-list sizes

    // Build cell → IP addressing
    CompactListList<label> cellToIP(layout.size(), layout.totalSize());

    for (label cellI = 0; cellI < layout.size(); ++cellI)
    {
        const label n = layout.localSize(cellI);
        for (label j = 0; j < n; ++j)
        {
            cellToIP(cellI, j) = layout.toGlobal(cellI, j);
        }
    }

    autoPtr<integrationPointTopology> topoPtr
    (
        new compactCellIntegrationPointTopology(mesh_, std::move(cellToIP))
    );

    // Cache and return
    topologyCache_.insert(key, topoPtr);

    return topology(*topologyCache_[key]).topology_;
}


const Foam::integrationPointTopology&
Foam::mechanicalConstitutiveLawManager::compactFaceTopologyFor
(
    const CompactListList<tensor>& layout
) const
{
    FatalErrorInFunction
        << "compactFaceTopologyFor is not yet implemented"
        << exit(FatalError);

    // Dummy return to silence compiler warnings
    return topologyFor(cellCentredIntegrationPointTopology::typeName);
}


Foam::mechanicalConstitutiveLawManager::topologyEntry&
Foam::mechanicalConstitutiveLawManager::topology
(
    const integrationPointTopology& topo
) const
{
    // Use the address of the topology object as a unique key
    const word key = Foam::name(reinterpret_cast<std::uintptr_t>(&topo));

    // Return existing entry if already constructed
    if (topologyEntries_.found(key))
    {
        return topologyEntries_[key];
    }

    // ---------------------------------------------------------------------
    // Construct new topology entry (lazy initialisation)
    // ---------------------------------------------------------------------

    DebugInfo
        << "Creating topologyEntry for " << key << endl;

    topologyEntries_.insert(key, topologyEntry(topo));
    topologyEntry& entry = topologyEntries_[key];

    const label nLaws = laws_.size();

    entry.lawIntegrationPointIDs_.setSize(nLaws);
    entry.states_.setSize(nLaws);
    entry.boundaryStates_.setSize(nLaws);

    // Detect whether this topology supports boundary integration points
    entry.boundaryAware_ = topo.boundaryAware();

    // ---------------------------------------------------------------------
    // Build integration-point addressing per law
    // ---------------------------------------------------------------------

    forAll(laws_, lawI)
    {
        DynamicList<label> ipIDs;
        labelHashSet seen;   // only used when needed

        const labelList& cells = lawCells_[lawI];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const labelUList cellIPs = topo.cellIntegrationPointIDs(cellI);

            forAll(cellIPs, j)
            {
                const label ip = cellIPs[j];

                if (topo.requiresUniqueIntegrationPointsPerMaterial())
                {
                    if (seen.insert(ip))
                    {
                        ipIDs.append(ip);
                    }
                }
                else
                {
                    ipIDs.append(ip);
                }
            }
        }

        entry.lawIntegrationPointIDs_[lawI].transfer(ipIDs);
    }


    // ---------------------------------------------------------------------
    // Allocate and initialise constitutive states (per law)
    // ---------------------------------------------------------------------

    forAll(laws_, lawI)
    {
        entry.states_.set
        (
            lawI,
            new mechanicalConstitutiveLawState
            (
                entry.lawIntegrationPointIDs_[lawI].size()
            )
        );

        // Law-specific state initialisation
        laws_[lawI].initialiseState(entry.states_[lawI]);
    }

    // ---------------------------------------------------------------------
    // Allocate and initialise boundary states (if applicable)
    // ---------------------------------------------------------------------

    if (entry.boundaryAware_)
    {
        forAll(laws_, lawI)
        {
            entry.boundaryStates_[lawI].setSize(mesh_.boundary().size());

            forAll(mesh_.boundary(), patchI)
            {
                const label nFaces =
                    lawBoundaryFaces_[lawI][patchI].size();

                entry.boundaryStates_[lawI].set
                (
                    patchI,
                    new mechanicalConstitutiveLawState(nFaces)
                );

                laws_[lawI].initialiseState
                (
                    entry.boundaryStates_[lawI][patchI]
                );
            }
        }
    }

    return entry;
}


void Foam::mechanicalConstitutiveLawManager::updateOldTimeIfNeeded()
{
    const label timeIndex = mesh_.time().timeIndex();

    if (timeIndex != curTimeIndex_)
    {
        if (debug)
        {
            InfoInFunction
                << "Updating old-time states for all topology entries"
                << endl;
        }

        // Loop over all topology entries
        forAllIter
        (
            HashTable<topologyEntry>, topologyEntries_, topoIter
        )
        {
            topologyEntry& entry = topoIter();

            // Internal states
            forAll(entry.states_, lawI)
            {
                entry.states_[lawI].storeOldTime();
            }

            // Boundary states (if present)
            if (entry.boundaryAware_)
            {
                forAll(entry.boundaryStates_, lawI)
                {
                    forAll(entry.boundaryStates_[lawI], patchI)
                    {
                        entry.boundaryStates_[lawI][patchI].storeOldTime();
                    }
                }
            }
        }

        curTimeIndex_ = timeIndex;
    }
}


Foam::surfaceSymmTensorField&
Foam::mechanicalConstitutiveLawManager::surfaceStressSum() const
{
    if (!surfaceStressSumPtr_.valid())
    {
        surfaceStressSumPtr_.reset
        (
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "surfaceStressSum",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            )
        );
    }

    return autoPtrRef(surfaceStressSumPtr_);
}


Foam::surfaceScalarField&
Foam::mechanicalConstitutiveLawManager::surfaceStressWeight() const
{
    if (!surfaceStressWeightPtr_.valid())
    {
        surfaceStressWeightPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "surfaceStressWeight",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }

    return autoPtrRef(surfaceStressWeightPtr_);
}


Foam::surfaceScalarField&
Foam::mechanicalConstitutiveLawManager::surfaceTangentWeight() const
{
    if (!surfaceTangentWeightPtr_.valid())
    {
        surfaceTangentWeightPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "surfaceTangentWeight",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimPressure, 0.0)
            )
        );
    }

    return autoPtrRef(surfaceTangentWeightPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawManager::mechanicalConstitutiveLawManager
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    curTimeIndex_(-1),
    laws_(),
    lawCells_(),
    lawBoundaryFaces_(),
    surfaceStressSumPtr_(),
    surfaceStressWeightPtr_(),
    surfaceTangentWeightPtr_(),
    rhoPtr_(),
    kappaPtr_(),
    topologyCache_(),
    topologyEntries_()
{
    // Read the mechanical laws
    const PtrList<entry> lawEntries(dict.lookup("mechanical"));

    laws_.setSize(lawEntries.size());
    lawCells_.setSize(lawEntries.size());
    lawBoundaryFaces_.setSize(lawEntries.size());

    // Look list of law names
    wordList lawNames(lawEntries.size());
    forAll(lawNames, lawI)
    {
        lawNames[lawI] = lawEntries[lawI].keyword();
    }

    // Create a map for each cell to its mechanical law
    labelList cellToLaw(mesh_.nCells(), -1);

    forAll(lawNames, lawI)
    {
        const word& lawName = lawNames[lawI];
        const dictionary& lawDict = lawEntries[lawI].dict();

        // Construct law
        laws_.set
        (
            lawI,
            mechanicalConstitutiveLaw::New(lawDict)
        );

        if (laws_.size() == 1)
        {
            lawCells_[lawI] = Foam::identity(mesh_.nCells());
        }
        else // more than one material law
        {
            // Look up cell zone of the same name as the law
            const label zoneID = mesh_.cellZones().findZoneID(lawName);

            if (zoneID < 0)
            {
                FatalErrorInFunction
                    << "CellZone " << lawName << " not found"
                    << "When defining more than one mechanical constitutive "
                    << "law, each cell must belong to exactly one cellZone "
                    << "with the same name as the law."
                    << exit(FatalError);
            }

            lawCells_[lawI] = mesh_.cellZones()[zoneID];
        }

        // Check that each cell appears in only one cell zone
        forAll(lawCells_[lawI], i)
        {
            const label cellI = lawCells_[lawI][i];

            if (cellToLaw[cellI] != -1)
            {
                FatalErrorInFunction
                    << "Cell " << cellI
                    << " appears in more than one mechanical constitutive law: "
                    << lawNames[cellToLaw[cellI]] << " and " << lawName
                    << exit(FatalError);
            }

            cellToLaw[cellI] = lawI;
        }
    }

    // Check for any cells without a material law
    forAll(cellToLaw, cellI)
    {
        if (cellToLaw[cellI] == -1)
        {
            FatalErrorInFunction
                << "Cell " << cellI
                << " is not assigned to any mechanical constitutive law. "
                << "When defining more than one material, every cell must "
                << "belong to exactly one cellZone."
                << exit(FatalError);
        }
    }

    // Set lawBoundaryFaces
    forAll(lawNames, lawI)
    {
        lawBoundaryFaces_[lawI].resize(mesh.boundary().size());

        forAll(lawBoundaryFaces_[lawI], patchI)
        {
            const labelList& faceCells = mesh.boundary()[patchI].faceCells();

            labelHashSet curFaceSet;

            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];

                if (cellToLaw[cellID] == lawI)
                {
                    curFaceSet.insert(faceI);
                }
            }

            lawBoundaryFaces_[lawI][patchI] = curFaceSet.toc();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawManager::~mechanicalConstitutiveLawManager()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const Foam::volScalarField& Foam::mechanicalConstitutiveLawManager::rho() const
{
    if (!rhoPtr_.valid())
    {
        rhoPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("rho", dimDensity, 0.0),
                "zeroGradient"
            )
        );

        volScalarField& rhoField = autoPtrRef(rhoPtr_);

        forAll(laws_, lawI)
        {
            const dimensionedScalar rhoLaw = laws_[lawI].rho();
            const labelList& cells = lawCells_[lawI];

            forAll(cells, i)
            {
                rhoField[cells[i]] = rhoLaw.value();
            }
        }

        rhoField.correctBoundaryConditions();
    }

    return rhoPtr_();
}


const Foam::volScalarField&
Foam::mechanicalConstitutiveLawManager::kappa() const
{
    if (!kappaPtr_.valid())
    {
        kappaPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "kappa",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("kappa", dimDensity, 0.0),
                "zeroGradient"
            )
        );

        volScalarField& kappaField = autoPtrRef(kappaPtr_);

        forAll(laws_, lawI)
        {
            const dimensionedScalar kappaLaw = laws_[lawI].kappa();
            const labelList& cells = lawCells_[lawI];

            forAll(cells, i)
            {
                kappaField[cells[i]] = kappaLaw.value();
            }
        }

        kappaField.correctBoundaryConditions();
    }

    return kappaPtr_();
}


void Foam::mechanicalConstitutiveLawManager::resetMaterialPropertyFields()
{
    rhoPtr_.clear();
    kappaPtr_.clear();
}


void Foam::mechanicalConstitutiveLawManager::updateStressSmallStrain
(
    const volTensorField& gradD,
    const volTensorField& gradD0,
    const scalar dt,
    volSymmTensorField& stress,
    volScalarField* scalarTangentPtr,
    const tangentRequest tangentReq
)
{
    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Look up the map and state for cell-based topologies
    const integrationPointTopology& topo =
        topologyFor(cellCentredIntegrationPointTopology::typeName);

    topologyEntry& tp = topology(topo);

    forAll(laws_, lawI)
    {
        const labelList& ipIDs = tp.lawIntegrationPointIDs_[lawI];

        // Update the internal field
        {
            // "View" into the gradD for this material => does not copy data
            const UIndirectList<tensor> gradDView
            (
                gradD.internalField(), ipIDs
            );

            // "View" into the gradD0 for this material => does not copy data
            const UIndirectList<tensor> gradD0View
            (
                gradD0.internalField(), ipIDs
            );

            // "View" into the stress for this material => does not copy data
            UIndirectList<symmTensor> stressView
            (
                stress.internalField(), ipIDs
            );

            // Create wrapper for kinematic data: input to material law
            // This does not copy data
            smallStrainMechanicalConstitutiveLawKinematics kin
            (
                gradDView, gradD0View, dt
            );

            // Create wrapper for output
            if (scalarTangentPtr && needsScalarTangent(tangentReq))
            {
                // "View" into the tangent for this material => does not copy
                // data
                UIndirectList<scalar> tangentView
                (
                    scalarTangentPtr->internalField(),
                    ipIDs
                );

                // Create wrapper for material: output from material law
                // This does not copy data
                mechanicalConstitutiveLawResponse response
                (
                    stressView,
                    tangentView,
                    tangentReq
                );

                // Update the material response, e.g. update the stress
                laws_[lawI].evaluate(kin, tp.states_[lawI], response);
            }
            else
            {
                // Create wrapper for material: output from material law
                // This does not copy data
                mechanicalConstitutiveLawResponse response
                (
                    stressView,
                    tangentReq
                );

                // Update the material response, e.g. update the stress
                laws_[lawI].evaluate(kin, tp.states_[lawI], response);
            }
        }

        // Optionally, update the boundary field
        // Boundary constitutive response uses independent state objects,
        // allowing history-dependent laws to operate correctly on boundary faces.
        if (tp.boundaryAware_)
        {
            forAll(gradD.boundaryField(), patchI)
            {
                if (!gradD.boundaryField()[patchI].coupled())
                {
                    // Select all faces on the patch for which the adjacent
                    // cell is in this material
                    const labelList& faces = lawBoundaryFaces_[lawI][patchI];

                    if (faces.empty())
                    {
                        continue;
                    }

                    // "View" into the kinematic and stress fields for this
                    // material => does not copy data
                    const UIndirectList<tensor> gradDView
                    (
                        gradD.boundaryField()[patchI], faces
                    );
                    const UIndirectList<tensor> gradD0View
                    (
                        gradD0.boundaryField()[patchI], faces
                    );
                    UIndirectList<symmTensor> stressView
                    (
                        stress.boundaryFieldRef()[patchI], faces
                    );

                    // Create wrapper for kinematic data: input to material law
                    // This does not copy data
                    smallStrainMechanicalConstitutiveLawKinematics kin
                    (
                        gradDView, gradD0View, dt
                    );

                    // Create wrapper for output
                    if (scalarTangentPtr && needsScalarTangent(tangentReq))
                    {
                        // "View" into the tangent for this material => does not copy
                        // data
                        UIndirectList<scalar> tangentView
                        (
                            scalarTangentPtr->boundaryField()[patchI],
                            faces
                        );

                        // Create wrapper for material: output from material law
                        // This does not copy data
                        mechanicalConstitutiveLawResponse response
                        (
                            stressView, tangentView, tangentReq
                        );

                        // Update the material response, e.g. update the stress
                        laws_[lawI].evaluate
                        (
                            kin, tp.boundaryStates_[lawI][patchI], response
                        );
                    }
                    else
                    {
                        // Create wrapper for material: output from material law
                        // This does not copy data
                        mechanicalConstitutiveLawResponse response
                        (
                            stressView, tangentReq
                        );

                        // Update the material response, e.g. update the stress
                        laws_[lawI].evaluate
                        (
                            kin, tp.boundaryStates_[lawI][patchI], response
                        );
                    }
                }
            }
        }
    }

    // Update boundaries including syncing coupled boundaries
    stress.correctBoundaryConditions();

    if (scalarTangentPtr && needsScalarTangent(tangentReq))
    {
        scalarTangentPtr->correctBoundaryConditions();
    }
}


void Foam::mechanicalConstitutiveLawManager::updateStressSmallStrain
(
    const surfaceTensorField& gradD,
    const surfaceTensorField& gradD0,
    const scalar dt,
    surfaceSymmTensorField& stress,
    const stressCollapseRule collapseRule,
    surfaceScalarField* scalarTangentPtr,
    const tangentRequest tangentReq
)
{
    FatalError
        << "face integration points not yet available" << exit(FatalError);

    // Look up the map and state for face-based topologies
    // const integrationPointTopology& topo =
    //     topologyFor(faceCentredIntegrationPointTopology::typeName);

    // topologyEntry& tp = topology(topo);

    /*
    if (!isA<faceCentredIntegrationPointTopology>(ipTopology_))
    {
        FatalErrorInFunction
            << "updateStressSmallStrain(surfaceField) requires "
            << "face-centred integration points, but topology is "
            << ipTopology_.type()
            << exit(FatalError);
    }

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Collapse rule: none → direct write, no accumulation allowed
    if (collapseRule == stressCollapseRule::none)
    {
        stress = symmTensor::zero;

        forAll(laws_, lawI)
        {
            const labelList& ipIDs = lawIntegrationPointIDs_[lawI];

            const UIndirectList<tensor> gradDView
            (
                gradD.internalField(),
                ipIDs
            );

            const UIndirectList<tensor> gradD0View
            (
                gradD0.internalField(),
                ipIDs
            );

            UIndirectList<symmTensor> stressView
            (
                stress.internalField(),
                ipIDs
            );

            smallStrainMechanicalConstitutiveLawKinematics kin
            (
                gradDView,
                gradD0View,
                dt
            );

            mechanicalConstitutiveLawResponse response
            (
                stressView,
                tangentReq
            );

            laws_[lawI].evaluate
            (
                kin,
                states_[lawI],
                response
            );
        }

        stress.correctBoundaryConditions();
        return;
    }

    // Accumulation path with stress "collapse" given by collapse rule, e.g.
    // average, harmonic

    surfaceSymmTensorField& stressSum = surfaceStressSum();
    surfaceScalarField& weightSum = surfaceStressWeight();

    stressSum = symmTensor::zero;
    weightSum = 0.0;

    surfaceScalarField* tangentWeightPtr = nullptr;

    if (scalarTangentPtr && needsScalarTangent(tangentReq))
    {
        surfaceScalarField& tangentWeight = surfaceTangentWeight();
        tangentWeight = 0.0;
        tangentWeightPtr = &tangentWeight;
    }

    // Loop over constitutive laws
    forAll(laws_, lawI)
    {
        const labelList& ipIDs = lawIntegrationPointIDs_[lawI];

        const UIndirectList<tensor> gradDView
        (
            gradD.internalField(),
            ipIDs
        );

        const UIndirectList<tensor> gradD0View
        (
            gradD0.internalField(),
            ipIDs
        );

        UIndirectList<symmTensor> stressView
        (
            stressSum.internalField(),
            ipIDs
        );

        smallStrainMechanicalConstitutiveLawKinematics kin
        (
            gradDView,
            gradD0View,
            dt
        );

        if (scalarTangentPtr && needsScalarTangent(tangentReq))
        {
            UIndirectList<scalar> tangentView
            (
                *tangentWeightPtr,
                ipIDs
            );

            mechanicalConstitutiveLawResponse response
            (
                stressView,
                tangentView,
                tangentReq
            );

            laws_[lawI].evaluate
            (
                kin,
                states_[lawI],
                response
            );
        }
        else
        {
            mechanicalConstitutiveLawResponse response
            (
                stressView,
                tangentReq
            );

            laws_[lawI].evaluate
            (
                kin,
                states_[lawI],
                response
            );
        }

        // Accumulate weights
        forAll(ipIDs, i)
        {
            weightSum[ipIDs[i]] += 1.0;
        }
    }

    // Collapse accumulated stress

    forAll(stress.internalField(), faceI)
    {
        const scalar w = weightSum[faceI];

        if (w > SMALL)
        {
            if (collapseRule == stressCollapseRule::average)
            {
                stress[faceI] = stressSum[faceI]/w;
            }
            else if
            (
                collapseRule == stressCollapseRule::harmonic
             && scalarTangentPtr
             && needsScalarTangent(tangentReq)
            )
            {
                stress[faceI] = stressSum[faceI]/max(w, SMALL);
                (*scalarTangentPtr)[faceI] =
                    w/max((*tangentWeightPtr)[faceI], SMALL);
            }
            else
            {
                FatalErrorInFunction
                    << "Invalid stress collapse rule combination:\n"
                    << "collapseRule = " << stressCollapseRuleName(collapseRule)
                    << nl
                    << "tangentReq   = " << tangentRequestName(tangentReq) << nl
                    << "scalarTangentPtr = "
                    << (scalarTangentPtr ? "set" : "null")
                    << exit(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction
                << "At least one mesh face does not have a material law!\n"
                << "faceI = " << faceI << nl
                << "Check cellZones and integrationPointTopology."
                << exit(FatalError);
        }
    }

    stress.correctBoundaryConditions();

    if (scalarTangentPtr && needsScalarTangent(tangentReq))
    {
        scalarTangentPtr->correctBoundaryConditions();
    }
        */
}


void Foam::mechanicalConstitutiveLawManager::updateStressSmallStrain
(
    const CompactListList<tensor>& gradD,
    const CompactListList<tensor>& gradD0,
    const scalar dt,
    CompactListList<symmTensor>& stress,
    List<scalar>* scalarTangentPtr,
    const tangentRequest tangentReq
)
{
    // Look up the map and state for compact list cell-based topologies
    const integrationPointTopology& topo = compactCellTopologyFor(gradD);

    topologyEntry& tp = topology(topo);

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Access packed integration-point storage
    const List<tensor>& gradDVals  = gradD.values();
    const List<tensor>& gradD0Vals = gradD0.values();
    List<symmTensor>& stressVals   = stress.values();

    // Consistency check => encapusate in a function: should be called in all
    // update stress functions!
    FatalError
        << "Consistency check => encapusate" << exit(FatalError);
    if
    (
        gradD.totalSize()  != stress.totalSize()
     || gradD0.totalSize() != stress.totalSize()
    )
    {
        FatalErrorInFunction
            << "Inconsistent CompactListList sizes supplied to "
            << "updateStressSmallStrain." << nl
            << "gradD size:  " << gradD.totalSize()  << nl
            << "gradD0 size: " << gradD0.totalSize() << nl
            << "stress size:" << stress.totalSize()
            << exit(FatalError);
    }

    if
    (
        scalarTangentPtr
     && scalarTangentPtr->size() != stress.totalSize()
    )
    {
        FatalErrorInFunction
            << "Scalar tangent list has incorrect size. Expected "
            << stress.totalSize() << " but got "
            << scalarTangentPtr->size()
            << exit(FatalError);
    }

    // Loop over mechanical constitutive laws
    forAll(laws_, lawI)
    {
        const labelList& ipIDs = tp.lawIntegrationPointIDs_[lawI];

        // Views into integration-point data (no copies)
        const UIndirectList<tensor> gradDView(gradDVals, ipIDs);

        const UIndirectList<tensor> gradD0View(gradD0Vals, ipIDs);

        UIndirectList<symmTensor> stressView(stressVals, ipIDs);

        // Kinematics wrapper
        smallStrainMechanicalConstitutiveLawKinematics kin
        (
            gradDView, gradD0View, dt
        );

        // Constitutive response
        if (scalarTangentPtr && needsScalarTangent(tangentReq))
        {
            UIndirectList<scalar> tangentView
            (
                *scalarTangentPtr, ipIDs
            );

            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentView, tangentReq
            );

            laws_[lawI].evaluate
            (
                kin, tp.states_[lawI], response
            );
        }
        else
        {
            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentReq
            );

            laws_[lawI].evaluate
            (
                kin, tp.states_[lawI], response
            );
        }
    }
}


void Foam::mechanicalConstitutiveLawManager::updateStressFiniteStrain
(
    const volTensorField& F,
    const volTensorField& F0,
    const volScalarField& J,
    const volScalarField& J0,
    const volTensorField& Finv,
    const volTensorField& Finv0,
    const scalar dt,
    volSymmTensorField& stress,
    volScalarField* scalarTangentPtr,
    const tangentRequest tangentReq
)
{
    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Look up the map and state for cell-based topologies
    const integrationPointTopology& topo =
        topologyFor(cellCentredIntegrationPointTopology::typeName);

    topologyEntry& tp = topology(topo);

    forAll(laws_, lawI)
    {
        const labelList& ipIDs = tp.lawIntegrationPointIDs_[lawI];

        // Update the internal field
        {
            // "View" into the F for this material => does not copy data
            const UIndirectList<tensor> FView(F.internalField(), ipIDs);

            // "View" into the F0 for this material => does not copy data
            const UIndirectList<tensor> F0View(F0.internalField(), ipIDs);

            // "View" into the J for this material => does not copy data
            const UIndirectList<scalar> JView(J.internalField(), ipIDs);

            // "View" into the J0 for this material => does not copy data
            const UIndirectList<scalar> J0View(J0.internalField(), ipIDs);

            // "View" into the Finv for this material => does not copy data
            const UIndirectList<tensor> FinvView(Finv.internalField(), ipIDs);

            // "View" into the Finv0 for this material => does not copy data
            const UIndirectList<tensor> Finv0View(Finv0.internalField(), ipIDs);

            // "View" into the stress for this material => does not copy data
            UIndirectList<symmTensor> stressView(stress.internalField(), ipIDs);

            // Create wrapper for kinematic data: input to material law
            // This does not copy data
            finiteStrainMechanicalConstitutiveLawKinematics kin
            (
                FView, F0View, JView, J0View, FinvView, Finv0View, dt
            );

            // Create wrapper for output
            if (scalarTangentPtr && needsScalarTangent(tangentReq))
            {
                // "View" into the tangent for this material => does not copy
                // data
                UIndirectList<scalar> tangentView
                (
                    scalarTangentPtr->internalField(), ipIDs
                );

                // Create wrapper for material: output from material law
                // This does not copy data
                mechanicalConstitutiveLawResponse response
                (
                    stressView, tangentView, tangentReq
                );

                // Update the material response, e.g. update the stress
                laws_[lawI].evaluate(kin, tp.states_[lawI], response);
            }
            else
            {
                // Create wrapper for material: output from material law
                // This does not copy data
                mechanicalConstitutiveLawResponse response
                (
                    stressView, tangentReq
                );

                // Update the material response, e.g. update the stress
                laws_[lawI].evaluate(kin, tp.states_[lawI], response);
            }
        }

        // Update the boundary field
        // Boundary constitutive response uses independent state objects,
        // allowing history-dependent laws to operate correctly on boundary faces.
        forAll(F.boundaryField(), patchI)
        {
            if (!F.boundaryField()[patchI].coupled())
            {
                // Select all faces on the patch for which the adjacent cell is
                // in this material
                const labelList& faces = lawBoundaryFaces_[lawI][patchI];

                // "View" into the J for this material => does not copy data
                const UIndirectList<scalar> JView
                (
                    J.boundaryField()[patchI], faces
                );

                // "View" into the J0 for this material => does not copy data
                const UIndirectList<scalar> J0View
                (
                    J0.boundaryField()[patchI], faces
                );

                // "View" into the F for this material => does not copy data
                const UIndirectList<tensor> FView
                (
                    F.boundaryField()[patchI], faces
                );

                // "View" into the F0 for this material => does not copy data
                const UIndirectList<tensor> F0View
                (
                    F0.boundaryField()[patchI], faces
                );

                // "View" into the Finv for this material => does not copy data
                const UIndirectList<tensor> FinvView
                (
                    Finv.boundaryField()[patchI], faces
                );

                // "View" into the Finv0 for this material => does not copy data
                const UIndirectList<tensor> Finv0View
                (
                    Finv0.boundaryField()[patchI], faces
                );

                // "View" into the stress for this material => does not copy data
                UIndirectList<symmTensor> stressView
                (
                    stress.boundaryFieldRef()[patchI], faces
                );

                // Create wrapper for kinematic data: input to material law
                // This does not copy data
                finiteStrainMechanicalConstitutiveLawKinematics kin
                (
                    FView, F0View, JView, J0View, FinvView, Finv0View, dt
                );

                // Create wrapper for output
                if (scalarTangentPtr && needsScalarTangent(tangentReq))
                {
                    // "View" into the tangent for this material => does not copy
                    // data
                    UIndirectList<scalar> tangentView
                    (
                        scalarTangentPtr->boundaryField()[patchI],
                        faces
                    );

                    // Create wrapper for material: output from material law
                    // This does not copy data
                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tangentView, tangentReq
                    );

                    // Update the material response, e.g. update the stress
                    laws_[lawI].evaluate
                    (
                        kin, tp.boundaryStates_[lawI][patchI], response
                    );
                }
                else
                {
                    // Create wrapper for material: output from material law
                    // This does not copy data
                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tangentReq
                    );

                    // Update the material response, e.g. update the stress
                    laws_[lawI].evaluate
                    (
                        kin, tp.boundaryStates_[lawI][patchI], response
                    );
                }
            }
        }
    }

    // Update boundaries including syncing coupled boundaries
    stress.correctBoundaryConditions();

    if (scalarTangentPtr && needsScalarTangent(tangentReq))
    {
        scalarTangentPtr->correctBoundaryConditions();
    }
}


void Foam::mechanicalConstitutiveLawManager::updateStressFiniteStrain
(
    const CompactListList<tensor>& F,
    const CompactListList<tensor>& F0,
    const CompactListList<tensor>& Finv,
    const CompactListList<tensor>& Finv0,
    const CompactListList<scalar>& J,
    const CompactListList<scalar>& J0,
    const scalar dt,
    CompactListList<symmTensor>& stress,
    List<scalar>* scalarTangentPtr,
    const tangentRequest tangentReq
)
{
    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Look up the map and state for compact list cell-based topologies
    const integrationPointTopology& topo = compactCellTopologyFor(F);

    topologyEntry& tp = topology(topo);

    // Access packed integration-point storage
    const List<tensor>& FVals = F.values();
    const List<tensor>& F0Vals = F0.values();
    const List<tensor>& FinvVals = Finv.values();
    const List<tensor>& Finv0Vals = Finv0.values();
    const List<scalar>& JVals = J.values();
    const List<scalar>& J0Vals = J0.values();
    List<symmTensor>& stressVals = stress.values();

    // Consistency checks
    const label nIP = stress.totalSize();

    if
    (
        F.totalSize() != nIP
     || F0.totalSize() != nIP
     || Finv.totalSize() != nIP
     || Finv0.totalSize() != nIP
     || J.totalSize() != nIP
     || J0.totalSize() != nIP
    )
    {
        FatalErrorInFunction
            << "Inconsistent CompactListList sizes supplied to "
            << "updateStressFiniteStrain." << nl
            << "Expected totalSize = " << nIP << exit(FatalError);
    }

    if (scalarTangentPtr && scalarTangentPtr->size() != nIP)
    {
        FatalErrorInFunction
            << "Scalar tangent list has incorrect size. Expected "
            << nIP << " but got " << scalarTangentPtr->size()
            << exit(FatalError);
    }

    // Loop over mechanical constitutive laws
    forAll(laws_, lawI)
    {
        const labelList& ipIDs = tp.lawIntegrationPointIDs_[lawI];

        // Views into integration-point data (no copies)
        const UIndirectList<tensor> FView(FVals, ipIDs);
        const UIndirectList<tensor> F0View(F0Vals, ipIDs);
        const UIndirectList<tensor> FinvView(FinvVals, ipIDs);
        const UIndirectList<tensor> Finv0View(Finv0Vals, ipIDs);
        const UIndirectList<scalar> JView(JVals, ipIDs);
        const UIndirectList<scalar> J0View(J0Vals, ipIDs);

        UIndirectList<symmTensor> stressView(stressVals, ipIDs);

        // Kinematics wrapper
        finiteStrainMechanicalConstitutiveLawKinematics kin
        (
            FView,
            F0View,
            JView,
            J0View,
            FinvView,
            Finv0View,
            dt
        );

        // Constitutive response
        if (scalarTangentPtr && needsScalarTangent(tangentReq))
        {
            UIndirectList<scalar> tangentView
            (
                *scalarTangentPtr, ipIDs
            );

            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentView, tangentReq
            );

            laws_[lawI].evaluate
            (
                kin, tp.states_[lawI], response
            );
        }
        else
        {
            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentReq
            );

            laws_[lawI].evaluate
            (
                kin, tp.states_[lawI], response
            );
        }
    }
}


void Foam::mechanicalConstitutiveLawManager::endTimeStep()
{
    const scalar time = mesh_.time().value();
    const label timeIndex = mesh_.time().timeIndex();

    // Loop over all topology entries
    forAllIter
    (
        HashTable<topologyEntry>, topologyEntries_, topoIter
    )
    {
        topologyEntry& tp = topoIter();

        forAll(laws_, lawI)
        {
            // Call endTimeStep with each topo type for each law
            laws_[lawI].endTimeStep(tp.states_[lawI], time, timeIndex);
        }
    }
}


// ************************************************************************* //
