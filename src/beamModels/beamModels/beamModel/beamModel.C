/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "beamModel.H"
#include "volFields.H"
#include "polyPatchID.H"

// #include "beamContactModel.H"
#include "IFstream.H"
#include "beamHelperFunctions.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(beamModel, 0);
    defineRunTimeSelectionTable(beamModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::beamModel::makeUpperLowerNeiCellFaces() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::makePointFaces() : "
            << "constructing upper and lower neighbour cell faces"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if
    (
        upperNeiCellFacesPtr_ || lowerNeiCellFacesPtr_
    )
    {
        FatalErrorIn
        (
            "beamModel::makeUpperLowerNeiCellFaces() const"
        )
            << "upper and lower neighbour cells and faces already exist"
            << abort(FatalError);
    }

    // Allocate storage for addressing
    upperNeiCellFacesPtr_ = new labelList(mesh().nCells(), -1);
    labelList& upperNeiCellFaces = *upperNeiCellFacesPtr_;

    // Allocate storage for addressing
    lowerNeiCellFacesPtr_ = new labelList(mesh().nCells(), -1);
    labelList& lowerNeiCellFaces = *lowerNeiCellFacesPtr_;

    const cellList& cells = mesh().cells();
    const labelListList& cc = mesh().cellCells();

    forAll(cells, cellI)
    {
        label upperNei = -1;
        label lowerNei = -1;

        forAll(cc[cellI], cI)
        {
            label curNei = cc[cellI][cI];

            if (curNei < cellI)
            {
                lowerNei = curNei;
            }
            else
            {
                upperNei = curNei;
            }
        }

        if (lowerNei != -1)
        {
            lowerNeiCellFaces[cellI] =
                sharedFace
                (
                    cells[cellI],
                    cells[lowerNei]
                );
        }
        else
        {
            lowerNeiCellFaces[cellI] = //-(startPatchIndex() + 1);
               -(neiPatchIndex(cells[cellI], mesh()) + 1);
        }

        if (upperNei != -1)
        {
            upperNeiCellFaces[cellI] =
                sharedFace
                (
                    cells[cellI],
                    cells[upperNei]
                );
        }
        else
        {
            upperNeiCellFaces[cellI] = //-(endPatchIndex() + 1);
                -(neiPatchIndex(cells[cellI], mesh()) + 1);
        }
    }
}

Foam::tmp<Foam::surfaceScalarField>
Foam::beamModel::indicator(const Foam::label bI) const
{
    tmp<surfaceScalarField> tIndicator
    (
        new surfaceScalarField
        (
            IOobject
            (
                "beamIndicator",
                runTime().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0)
        )
    );

    const labelList& beamCells = mesh().cellZones()[bI];

    volScalarField cellIndicator
    (
        IOobject
        (
            "cellIndicator",
            runTime().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(beamCells, cellI)
    {
        cellIndicator[beamCells[cellI]] = 1.0;
    }
    cellIndicator.correctBoundaryConditions();

    tIndicator.ref() = fvc::interpolate(cellIndicator);

    return tIndicator;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::beamModel::beamModel
(
    const word& type,
    Time& runTime,
    const word& region
)
:
    IOdictionary
    (
        IOobject
        (
            "beamProperties",
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(runTime),
    meshOwnerPtr_
    (
        // Use the region mesh if another model has already registered it
        runTime.foundObject<dynamicFvMesh>(region)
      ? autoPtr<dynamicFvMesh>()
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
      ? &meshOwnerPtr_()
      : &runTime.lookupObject<dynamicFvMesh>(region)
    ),
    localToGlobalCellAddressing_(),
    globalToLocalCellAddressing_(),
    localToGlobalBeamPointsAddressing_(),
    globalToLocalBeamPointsAddressing_(),
    upperNeiCellFacesPtr_(NULL),
    lowerNeiCellFacesPtr_(NULL),
    beamProperties_(subDict(type + "Coeffs")),
    startPatchIndex_(),
    endPatchIndex_(),
    startCells_(),
    nBeamCells_(),
    crossSections_(),
    R_(),
    // U_(),
    // E_("E",beamProperties()),
    // G_("G",beamProperties()),
    E_(),
    G_(),
    rho_
    (
        // beamProperties().getOrDefault<dimensionedScalar>
        // (
        //     "rho",
        //     dimensionedScalar("rho", dimDensity, 0)
        // )
    ),
    rhof_
    (
        beamProperties().getOrDefault<dimensionedScalar>
        (
            "rhoFluid",
            dimensionedScalar("rhoFluid", dimDensity, 0)
        )
    ),
    gPtr_(),
    L_
    (
        IOobject
        (
            "L",
            runTime.timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("L", dimLength, 0)
    ),
    q_
    (
        IOobject
        (
            "q",
            runTime.timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedVector("q", dimForce/dimLength, vector::zero)
    ),
    m_
    (
        IOobject
        (
            "m",
            runTime.timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedVector("m", dimForce, vector::zero)
    ),
    pointForces_(),
    iOuterCorr_(0),
    kCI_(beamProperties_.getOrDefault<scalar>("scalingInertiaTensor", 1)),
    // contactPtr_(),
    deltaTseries_()
{
    // if (beamProperties().found("rho"))
    // {
    //     rho_ = dimensionedScalar("rho",beamProperties());
    // }

     // Header for gravitational acceleration
    IOobject gHeader
    (
        "g",
        runTime.constant(),
        mesh(),
        IOobject::MUST_READ
    );

    if (gHeader.typeHeaderOk<uniformDimensionedVectorField>(true))
    {
        gPtr_.reset
        (
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.caseConstant(),
                    this->mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
           )
        );

        Info<< "g is set for gravitational body force: "
            << g()
            << endl;

        // To check whether density of fluid is specified to calculate
        // buoyant body force
        if (beamProperties().found("rhoFluid"))
        {
            rhof_ = dimensionedScalar("rhoFluid",beamProperties());
        }
        else
        {
            WarningIn("Constructor of beamModel.C")
                << "density of fluid rhoFluid is not set in constant/beamProperties "
                << "--> buoyant body force = 0"
                << nl << endl;
        }
    }
    else
    {
        gPtr_.reset
        (
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.caseConstant(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
            dimensionedVector("g", dimVelocity/dimTime, vector::zero)
            )
        );
    }

    // Read cross-sections of beams
    // if (found("beams"))
    // {
    Info<< "Reading beam properties from cross-section model" << endl;

    const PtrList<entry> entries(lookup("beams"));

    const label nBeams = entries.size();

    Info<< "nBeams: " << nBeams << endl;

    crossSections_.setSize(nBeams);

    forAll(entries, beamI)
    {
        crossSections_.set
        (
            beamI,
            crossSectionModel::New
            (
                word(entries[beamI].dict().lookup("crossSectionModel")),
                entries[beamI].dict()
            )
        );
    }

    R_.setSize(nBeams);
    // U_.setSize(nBeams);
    forAll(crossSections_, beamI)
    {
        R_[beamI] = crossSections_[beamI].R();
        // U_[beamI] = crossSections_[beamI].U();
    }
    // }
    // else if (this->found("nBeams")) // Read beam radius
    // {
    //     label nBeams = readInt(this->lookup("nBeams"));
    //     R_.setSize
    //     (
    //         nBeams,
    //         dimensionedScalar(this->lookup("R")).value()
    //      );
    //     if (this->found("U"))
    //     {
    //         U_.setSize
    //         (
    //             nBeams,
    //             dimensionedScalar(this->lookup("U")).value()
    //         );
    //     }
    //     else
    //     {
    //         U_.setSize(nBeams, 0);
    //     }
    // }
    // else
    // {
    //     R_ = scalarField(this->lookup("R"));
    //     U_ = scalarField(this->lookup("U"));
    // }

    Info<< "Attempting to read global mechanical variables "
        << "from constant/beamProperties "
        << "coupledTotalLagNetwonRaphsonBeamCoeffs() dictionary!"
        << endl;

    dimensionedScalar globalE("E", dimPressure, 0.0);
    dimensionedScalar globalG("G", dimPressure, 0.0);
    dimensionedScalar globalRho("rho", dimDensity, 0.0);

    bool foundGlobalE = false;
    bool foundGlobalG = false;
    bool foundGlobalRho = false;

    if (beamProperties().found("E"))
    {
        globalE = dimensionedScalar("E", beamProperties());
        foundGlobalE = true;
    }

    if (beamProperties().found("G"))
    {
        globalG = dimensionedScalar("G", beamProperties());
        foundGlobalG = true;
    }

    if (beamProperties().found("rho"))
    {
        globalRho = dimensionedScalar("rho", beamProperties());
        foundGlobalRho = true;
    }

    E_.setSize(nBeams, 0.0);
    G_.setSize(nBeams, 0.0);
    rho_.setSize(nBeams, 0.0);

    forAll(crossSections_, beamI)
    {
        const dictionary& bDict = entries[beamI].dict();

        if (
               bDict.found("E")
            || bDict.found("G")
            || bDict.found("rho")
           )
        {
            if (foundGlobalE || foundGlobalG || foundGlobalRho)
            {
                FatalErrorInFunction
                    << "Both local and global inputs for "
                    << "beam mechanical properties found!\n"
                    << "Either specify the variables (E, G, rho) of type "
                    << "*scalar* locally inside crossSectionModel dictionary, OR "
                    << "specify ALL of them as global variables of type "
                    << "*dimensionedScalar*, i.e., specify with dimensions.\n"
                    << "The global values are common "
                    << "to all beams!\n"
                    << abort(FatalError);
            }
            else
            {
                Info<< "Local mechanical properties specific to beam "
                    << beamI << " found!"  << endl;

                E_[beamI] = readScalar(bDict.lookup("E"));
                G_[beamI] = readScalar(bDict.lookup("G"));
                rho_[beamI] = bDict.getOrDefault<scalar>("rho", 0.0);
            }
        }
        else if (foundGlobalE || foundGlobalG || foundGlobalRho)
        {
            E_[beamI] = globalE.value();
            G_[beamI] = globalG.value();
            rho_[beamI] = globalRho.value();
        }
        else
        {
            FatalErrorInFunction
                << "Missing E, G: neither per-beam 'E'/'G' for beam " << beamI
                << " nor global 'E'/'G' in " << beamProperties().name() << "."
                << endl;
        }

    }

    const bool gravityBodyForce =
        mag(g().component(0).value()) > SMALL
     || mag(g().component(1).value()) > SMALL
     || mag(g().component(2).value()) > SMALL;

    if (gravityBodyForce)
    {
        forAll(rho_, beamI)
        {
            if (mag(rho_[beamI]) < SMALL)
            {
                WarningInFunction
                    << "rho for beam " << beamI << " = zero"
                    << " --> gravitational body force = 0 for this beam"
                    << nl << endl;
            }
        }
    }

    // Write beam cross-section properties
    for (label bI = 0; bI < nBeams; ++bI)
    {
        Info<< "\n\n Writing mechanical properties of beam no: " << bI << endl;
        Info<< "E " << E(bI) << endl;
        Info<< "G " << G(bI) << endl;
        Info<< "rho " << rho(bI) << endl;

        Info<< "A: " << A(bI) << endl;

        Info<< "Iyy: " << Iyy(bI) << endl;
        Info<< "Izz: " << Izz(bI) << endl;

        Info<< "I: " << I(bI) << endl;
        Info<< "J: " << J(bI) << endl;

        Info<< "EI: " << EI(bI) << endl;
        Info<< "GJ: " << GJ(bI) << endl;

        Info<< "EA: " << EA(bI) << endl;
        Info<< "GA: " << GA(bI) << endl;
    }

    // word startPatchName(beamProperties_.lookup("startPatchName"));
    // word endPatchName(beamProperties_.lookup("endPatchName"));
    word startPatchName = "left";
    word endPatchName = "right";

    label nCellZones = mesh().cellZones().size();

    // Check that there is at least one cell zone
    if (nCellZones == 0)
    {
        FatalErrorIn("beamModel::beamModel(...)")
            << "No cell zones are defined: one cell zone per beam is required"
            << abort(FatalError);
    }
    else if (nCellZones < 2)
    {
        polyPatchID startPatch(startPatchName, mesh().boundaryMesh());

        if (!startPatch.active())
        {
            FatalErrorIn("beamModel::beamModel(...)")
              << "Patch name " << startPatchName << " not found."
              << abort(FatalError);
        }

        startPatchIndex_ = labelList(1, startPatch.index());

        polyPatchID endPatch(endPatchName, mesh().boundaryMesh());

        if (!endPatch.active())
        {
            FatalErrorIn("beamModel::beamModel(...)")
              << "Patch name " << endPatchName << " not found."
              << abort(FatalError);
        }

        endPatchIndex_ = labelList(1, endPatch.index());
    }
    else
    {
        // Start patches
        startPatchIndex_ = labelList(nCellZones, -1);

        polyPatchID startPatch(startPatchName, mesh().boundaryMesh());

        if (startPatch.active())
        {
            startPatchIndex_ = startPatch.index();
        }
        else
        {
            forAll(startPatchIndex_, zI)
            {
                OStringStream StartPatchName;
                StartPatchName() << startPatchName << '_' << zI;

                word zoneStartPatchName(StartPatchName.str());

                polyPatchID startPatch
                (
                    zoneStartPatchName,
                    mesh().boundaryMesh()
                );

                if (!startPatch.active())
                {
                    FatalErrorIn("beamModel::beamModel(...)")
                        << "Patch name " << zoneStartPatchName
                        << " not found."
                        << abort(FatalError);
                }

                startPatchIndex_[zI] = startPatch.index();
            }
        }

        // End patches
        endPatchIndex_ = labelList(nCellZones, -1);

        polyPatchID endPatch(endPatchName, mesh().boundaryMesh());

        if (endPatch.active())
        {
            endPatchIndex_ = endPatch.index();
        }
        else
        {
            forAll(endPatchIndex_, zI)
            {
                OStringStream EndPatchName;
                EndPatchName() << endPatchName << '_' << zI;

                word zoneEndPatchName(EndPatchName.str());

                polyPatchID endPatch
                (
                    zoneEndPatchName,
                    mesh().boundaryMesh()
                );

                if (!endPatch.active())
                {
                    FatalErrorIn("beamModel::beamModel(...)")
                        << "Patch name " << zoneEndPatchName
                        << " not found."
                        << abort(FatalError);
                }

                endPatchIndex_[zI] = endPatch.index();
            }
        }
    }

    // Calculate start cells
    nBeamCells_.setSize(nCellZones, 0);
    startCells_.setSize(nCellZones, -1);
    startCells_[0] = 0;
    for(label i=1; i<nCellZones; i++)
    {
        const cellZone& cz = mesh().cellZones()[i-1];

        nBeamCells_[i-1] = cz.size();

        // Global zone size in parallel
        if (Pstream::parRun())
        {
            reduce(nBeamCells_[i-1], sumOp<scalar>());
        }

        startCells_[i] = startCells_[i-1] + nBeamCells_[i-1];
    }

    nBeamCells_[nCellZones-1] =
        mesh().cellZones()[nCellZones-1].size();

    // Global zone size in parallel
    if (Pstream::parRun())
    {
        reduce(nBeamCells_[nCellZones-1], sumOp<scalar>());
    }

    // Create global-to-local and local-to-global cell addressing
    // if (Pstream::parRun())
    // {
    //     localToGlobalCellAddressing_ =
    //         labelIOList
    //         (
    //             IOobject
    //             (
    //                 "cellProcAddressing",
    //                 mesh().facesInstance(),
    //                 mesh().meshSubDir,
    //                 mesh(),
    //                 IOobject::MUST_READ,
    //                 IOobject::NO_WRITE
    //             )
    //         );

    //     label nGlobalCells = sum(nBeamCells_);

    //     globalToLocalCellAddressing_.setSize(Pstream::nProcs());

    //     forAll(globalToLocalCellAddressing_, procI)
    //     {
    //         globalToLocalCellAddressing_[procI] =
    //             labelList(nGlobalCells, -1);
    //     }

    //     forAll(globalToLocalCellAddressing_[Pstream::myProcNo()], gcI)
    //     {
    //         label locCellIndex =
    //             findIndex(localToGlobalCellAddressing_, gcI);

    //         if (locCellIndex != -1)
    //         {
    //             locCellIndex += csrAddr().globalNCellsOffset();
    //         }

    //         globalToLocalCellAddressing_[Pstream::myProcNo()][gcI] =
    //             locCellIndex;
    //     }

    //     // Send data
    //     for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
    //     {
    //         if (domainI != Pstream::myProcNo())
    //         {
    //             // if (allPoints[Pstream::myProcNo()].size())
    //             {
    //                 // Parallel data exchange
    //                 OPstream::write
    //                 (
    //                     Pstream::blocking,
    //                     domainI,
    //                     reinterpret_cast<const char*>
    //                     (
    //                         globalToLocalCellAddressing_
    //                         [Pstream::myProcNo()].begin()
    //                     ),
    //                     globalToLocalCellAddressing_
    //                     [Pstream::myProcNo()].byteSize()
    //                 );
    //             }
    //         }
    //     }

    //     // Receive data
    //     for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
    //     {
    //         if (domainI != Pstream::myProcNo())
    //         {
    //             // if (procNPoints[domainI])
    //             {
    //                 // Parallel data exchange
    //                 IPstream::read
    //                 (
    //                     Pstream::blocking,
    //                     domainI,
    //                     reinterpret_cast<char*>
    //                     (
    //                         globalToLocalCellAddressing_[domainI].begin()
    //                     ),
    //                     globalToLocalCellAddressing_[domainI].byteSize()
    //                 );
    //             }
    //         }
    //     }

    //     // Pout << globalToLocalCellAddressing_ << endl;
    //     // sleep(5);
    // }

    // Calculate beam points local-to-global and global-to-local addressing
    // {
    //     label nBeams = mesh().cellZones().size();

    //     localToGlobalBeamPointsAddressing_.setSize(nBeams);

    //     forAll(localToGlobalBeamPointsAddressing_, bI)
    //     {
    //         localToGlobalBeamPointsAddressing_.set
    //         (
    //             bI,
    //             new labelListList(Pstream::nProcs())
    //         );
    //     }

    //     forAll(localToGlobalBeamPointsAddressing_, bI)
    //     {
    //         forAll(localToGlobalBeamPointsAddressing_[bI], procI)
    //         {
    //             if (procI == Pstream::myProcNo())
    //             {
    //                 localToGlobalBeamPointsAddressing_[bI][procI] =
    //                     globalPointsIndices(bI);
    //             }
    //         }

    //         Pstream::gatherList(localToGlobalBeamPointsAddressing_[bI]);
    //         Pstream::scatterList(localToGlobalBeamPointsAddressing_[bI]);
    //     }

    //     globalToLocalBeamPointsAddressing_.setSize(nBeams);

    //     forAll(globalToLocalBeamPointsAddressing_, bI)
    //     {
    //         globalToLocalBeamPointsAddressing_[bI] =
    //             localPointsIndices(bI);
    //     }
    // }

    // Pout << "nBeamCells: " << nBeamCells_ << endl;
    // sleep(5);

    // Read point forces
    IFstream pointForcesStream(runTime.caseConstant()/word("pointForces"));
    if (pointForcesStream.good())
    {
        pointForcesStream >> pointForces_;
    }
    else
    {
        Info<< "Point forces do not exist" << endl;
    }

    // Read conical pulleys
    // if (found("conicalPulleys"))
    // {
    //     // Info<< "Found conical pulleys" << endl;

    //     const PtrList<entry> entries(lookup("conicalPulleys"));

    //     label nPulleys = entries.size();

    //     // Info<< "nPulleys: " << nPulleys << endl;

    //     conicalPulleys_.setSize(nPulleys);

    //     forAll(entries, pulleyI)
    //     {
    //         conicalPulleys_.set
    //         (
    //             pulleyI,
    //             new conicalPulley(entries[pulleyI].dict())
    //         );

    //         // Info<< entries[pulleyI].dict() << endl;
    //     }
    // }

    // Read toroidal pulleys
    // if (found("toroidalPulleys"))
    // {
    //     const PtrList<entry> entries(lookup("toroidalPulleys"));

    //     label nPulleys = entries.size();

    //     toroidalPulleys_.setSize(nPulleys);

    //     forAll(entries, pulleyI)
    //     {
    //         toroidalPulleys_.set
    //         (
    //             pulleyI,
    //             new toroidalPulley(entries[pulleyI].dict())
    //         );
    //     }
    // }

    // Read density if present
    //    if (beamProperties().found("rho"))
    //    {
    //        rho_ = dimensionedScalar(beamProperties().lookup("rho"));
    //    }

    // Check if deltaT is time-varying
    if (beamProperties().found("deltaTseries"))
    {
        Info<< "deltaT is time-varying" << endl;
        deltaTseries_ =
            interpolationTable<scalar>
            (
                beamProperties().subDict("deltaTseries")
            );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::beamModel::~beamModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::beamModel::evolve()
{
    // if (conicalPulleys_.size())
    // {
    //     if (runTime().value() > deleteConicalPulleysAt_)
    //     {
    //         Info<< "Delete conical pulleys ";
    //         conicalPulleys_.clear();
    //         Info<< conicalPulleys_.size() << endl;
    //     }
    // }

    return 0;
}

// const Foam::beamContactModel& Foam::beamModel::contact() const
// {
//     if (contactPtr_.empty())
//     {
//         FatalErrorIn("beamModel::contact() const")
//           << "Contact is not updated"
//           << abort(FatalError);
//     }

//     return contactPtr_();
// }

// Foam::beamContactModel& Foam::beamModel::contact()
// {
//     if (contactPtr_.empty())
//     {
//         contactPtr_.set(new beamContactModel(*this));
//     }

//     return contactPtr_();
// }

bool Foam::beamModel::read()
{
    if (regIOobject::read())
    {
        beamProperties_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}

void Foam::beamModel::writeVTK() const
{
    // Create directory if does not exist.
    fileName vtkDir(runTime().path()/"VTK");
    mkDir(vtkDir);

    OStringStream FileName;
    FileName() << mesh().name() << "_" << runTime().timeIndex() << ".vtk";

    fileName vtkFileName(word(FileName.str()));

    OFstream vtkFile(vtkDir/vtkFileName);

    // Write header
    vtkFile << "# vtk DataFile Version 2.0" << endl;
    vtkFile << "\nASCII" << endl;
    vtkFile << "\nDATASET UNSTRUCTURED_GRID" << endl;

    // Write points
    vectorField curPoints(currentBeamPoints());
    vtkFile << "\nPOINTS " << curPoints.size() << " float" << endl;
    for (label i=0; i<curPoints.size(); i++)
    {
        vtkFile << curPoints[i].x() << " "
                << curPoints[i].y() << " "
                << curPoints[i].z() << endl;
    }

    // Write cells
    label nCells = curPoints.size()-1;
    vtkFile << "\nCELLS " << nCells << " " << 3*nCells << endl;
    for (label i=0; i<nCells; i++)
    {
        vtkFile << 2 << " " << i << " " << i+1 << endl;
    }

    // Write cell types
    vtkFile << "\nCELL_TYPES " << nCells << endl;
    for (label i=0; i<nCells; i++)
    {
        vtkFile << 3 << endl;
    }

    //Write data
    vtkFile << "\nPOINT_DATA " << curPoints.size() << endl;
    vtkFile << "SCALARS scalars float 1" << endl;
    vtkFile << "LOOKUP_TABLE default" << endl;
    for (label i=0; i<curPoints.size(); i++)
    {
        vtkFile << i+1 << endl;
    }
}

Foam::tmp<Foam::vectorField> Foam::beamModel::points(const label bI) const
{
    const fvMesh& mesh = this->mesh();

    label nCellZones = mesh.cellZones().size();

    if (nCellZones < 2)
    {
        label nPoints = mesh.nCells() + 1;

        tmp<vectorField> tCurrentPoints
        (
            new vectorField(nPoints, vector::zero)
        );
        vectorField& curPoints = tCurrentPoints.ref();

        const surfaceVectorField& curCf = mesh.Cf();
        const vectorField& curCfI = curCf.internalField();

        curPoints[0] = curCf.boundaryField()[startPatchIndex()][0];
        curPoints[nPoints-1] = curCf.boundaryField()[endPatchIndex()][0];
        for (label i=0; i<curCfI.size(); i++)
        {
            curPoints[i+1] = curCfI[i];
        }

        return tCurrentPoints;
    }
    else
    {
        const cellZone& cz = mesh.cellZones()[bI];

        label nPoints = cz.size() + 1;

        tmp<vectorField> tCurrentPoints
        (
            new vectorField(nPoints, vector::zero)
        );
        vectorField& curPoints = tCurrentPoints.ref();
        const surfaceVectorField& curCf = mesh.Cf();
        const vectorField& curCfI = curCf.internalField();

        const labelList startPatchCells =
            mesh.boundary()[startPatchIndex(bI)].faceCells();

        forAll(curCf.boundaryField()[startPatchIndex(bI)], faceI)
        {
            if (cz.whichCell(startPatchCells[faceI]) != -1)
            {
                curPoints[0] =
                    curCf.boundaryField()[startPatchIndex(bI)][faceI];
            }
        }

        const labelList endPatchCells =
            mesh.boundary()[endPatchIndex(bI)].faceCells();

        forAll(curCf.boundaryField()[startPatchIndex(bI)], faceI)
        {
            if (cz.whichCell(endPatchCells[faceI]) != -1)
            {
                curPoints[nPoints-1] =
                    curCf.boundaryField()[endPatchIndex(bI)][faceI];
            }
        }

        const labelList& nei = mesh.neighbour();

        label index = 1;
        for (label i=0; i<curCfI.size(); i++)
        {
            if (cz.whichCell(nei[i]) != -1)
            {
                curPoints[index] = curCfI[i];
                index++;
            }
        }

        return tCurrentPoints;
    }
}

void Foam::beamModel::writeFields()
{
    runTime().write();

    // if (runTime().outputTime() && contactActive() && Pstream::master())
    // {
    //     // Write contact force and contact gap for first beam
    //     label nBeams = contact().splines().size();
    //     for (label bI=0; bI<nBeams; bI++)
    //     {
    //         // label bI = 0;
    //         // cubicSpline spline
    //         // (
    //         //     points(bI),
    //         //     cubicSpline::CLAMPED_CALC,
    //         //     vector(0, 0, 0),
    //         //     cubicSpline::CLAMPED_CALC,
    //         //     vector(0, 0, 0)
    //         // );

    //         scalarField t = contact().splines()[bI].midPointParameters();
    //         // scalarField t = spline.midPointParameters();
    //         // Info<< "First beam length: "
    //         //      << sum(spline.segLengths()) << endl;

    //         if (false)
    //         {
    //             Info<< "\nTotal forces for beam " << bI << endl;
    //             for (label i=0; i<nBeams; i++)
    //             {
    //                 if (i != bI)
    //                 {
    //                     Info<< ' ' << sum(contact().contactForces()[bI][i]);
    //                 }
    //                 else
    //                 {
    //                     Info<< ' ' << 0;
    //                 }
    //             }
    //             Info<< endl;
    //         }


    //         {
    //             OStringStream FileName;
    //             FileName() << "beam-" << bI << "_force.txt";

    //             // OFstream forceFile(runTime().timePath()/"beam-0_force.txt");
    //             OFstream forceFile
    //             (
    //                 runTime().timePath()/word(FileName.str())
    //             );
    //             forAll(t, segI)
    //             {
    //                 forceFile << t[segI];
    //                 for (label i=0; i<nBeams; i++)
    //                 {
    //                     if (i != bI)
    //                     {
    //                         forceFile << ' ' <<
    //                  mag
    //                  (
    //                                 contact().lineContacts()[bI][segI][i]
    //                     .normalContactForce()
    //                  );
    //                     }
    //                     else
    //                     {
    //                         forceFile << ' ' << 0;
    //                     }
    //                 }
    //                 forceFile << endl;
    //             }
    //         }


    //         {
    //             OStringStream FileName;
    //             FileName() << "beam-" << bI << "_gap.txt";

    //             OFstream gapFile
    //             (
    //                 runTime().timePath()/word(FileName.str())
    //             );

    //             // OFstream gapFile(runTime().timePath()/"beam-0_gap.txt");
    //             forAll(t, segI)
    //             {
    //                 gapFile << t[segI];
    //                 for (label i=0; i<nBeams; i++)
    //                 {
    //                     if (i != bI)
    //                     {
    //                         gapFile << ' ' <<
    //                  contact().lineContacts()[bI][segI][i]
    //                 .normalGap();

    //                         // gapFile << ' '
    //                         //     << contact().contactGaps()[bI][i][segI];
    //                     }
    //                     else
    //                     {
    //                         gapFile << ' ' << 0;
    //                     }
    //                 }
    //                 gapFile << endl;
    //             }
    //         }
    //     }
    // }


    // if (runTime().outputTime() && conicalPulleys_.size())
    // {
    //     // Write moving pulleys
    //     forAll(conicalPulleys_, pulleyI)
    //     {
    //         if (conicalPulleys_[pulleyI].moving())
    //         {
    //             // Create directory if does not exist.
    //             fileName vtkDir(runTime().path()/"VTK");
    //             mkDir(vtkDir);

    //             if (conicalPulleys_[pulleyI].stlModel().size())
    //             {
    //                 fileName conicalPulleysDir
    //                 (
    //                     vtkDir/"conicalPulleys"
    //                 );
    //                 mkDir(conicalPulleysDir);

    //                 OStringStream conicalPulleyNvtk;
    //                 conicalPulleyNvtk()
    //                     << "conicalPulley-" << pulleyI
    //                     << "_" << runTime().timeIndex() << ".vtk";

    //                 fileName vtkFileName
    //                 (
    //                     // "conicalPulleys"/
    //                     conicalPulleysDir/word(conicalPulleyNvtk.str())
    //                 );

    //                 OStringStream conicalPulleyNvtp;
    //                 conicalPulleyNvtp()
    //                     << "conicalPulley-" << pulleyI
    //                     << "_" << runTime().timeIndex() << ".vtp";

    //                 fileName vtpFileName
    //                 (
    //                     word(conicalPulleyNvtp.str())
    //                 );

    //                 vectorField oldPoints =
    //                     conicalPulleys_[pulleyI].stlModel().points();

    //                 vectorField newPoints =
    //                     oldPoints
    //                   + conicalPulleys_[pulleyI].temporalOriginDisplacement();

    //                 conicalPulleys_[pulleyI].stlModel().movePoints(newPoints);
    //                 conicalPulleys_[pulleyI].stlModel().write(vtkFileName);
    //                 conicalPulleys_[pulleyI].stlModel().movePoints(oldPoints);

    //                 // Write temporal collection data
    //                 {
    //                     OStringStream conicalPulleyNpvd;
    //                     conicalPulleyNpvd() << "conicalPulley-"
    //                                         << pulleyI << ".pvd";

    //                     fileName collectionFileName
    //                     (
    //                         conicalPulleysDir/conicalPulleyNpvd.str()
    //                     );
    //                     ifstream collectionFile(collectionFileName);
    //                     // IFstream collectionFile(collectionFileName);

    //                     fileName newCollectionFileName
    //                     (
    //                         conicalPulleysDir/"new.pvd"
    //                     );

    //                     if (collectionFile.good())
    //                     {
    //                         ofstream newCollectionFile
    //                         (
    //                             newCollectionFileName
    //                         );

    //                         // Add to existing collectin file
    //                         label lineIndex = 0;
    //                         do
    //                         {
    //                             lineIndex++;

    //                             std::string line;
    //                             std::getline(collectionFile, line);

    //                             // collectionFile.getLine(line);
    //                             // newCollectionFile << line << '\n';
    //                             // Info<< line << endl;

    //                             if
    //                             (
    //                                 (lineIndex < 3)
    //                              || (
    //                                     line.find("DataSet")
    //                                  != std::string::npos
    //                                 )
    //                             )
    //                             {
    //                                 newCollectionFile << line << '\n';
    //                             }
    //                         }
    //                         while(!collectionFile.eof());

    //                         // Add current pulley data
    //                         newCollectionFile
    //                             << "    <DataSet timestep=\""
    //                             << runTime().value()
    //                             << "\" group=\"\" part=\"0\" file=\""
    //                             << vtpFileName << "\"/>" << '\n';

    //                         // Add last two lines
    //                         newCollectionFile << "  </Collection>" << '\n';
    //                         newCollectionFile << "</VTKFile>";
    //                     }
    //                     else
    //                     {
    //                         OFstream newCollectionFile
    //                         (
    //                             newCollectionFileName
    //                         );

    //                         // Add first two lines
    //                         newCollectionFile
    //                             << "<VTKFile type=\"Collection\" version=\"0.1\" "
    //                             << "byte_order=\"LittleEndian\">" << endl;
    //                         newCollectionFile << "  <Collection>" << endl;

    //                         // Add current pulley data
    //                         newCollectionFile
    //                             << "    <DataSet timestep=\""
    //                             << runTime().value()
    //                             << "\" group=\"\" part=\"0\" file="
    //                             << vtpFileName << "/>" << endl;

    //                         // Add last two lines
    //                         newCollectionFile << "  </Collection>" << endl;
    //                         newCollectionFile << "</VTKFile>";
    //                     }

    //                     mv(newCollectionFileName, collectionFileName);

    //                     // collectionFile.close();
    //                     // newCollectionFile.close();

    //                     // if
    //                     // (
    //                     //     std::rename
    //                     //     (
    //                     //         collectionFile.name().c_str(),
    //                     //         newCollectionFile.name().c_str()
    //                     //     )
    //                     // )
    //                     // {
    //                     //     FatalErrorIn("void Foam::beamModel::writeFields()")
    //                     //         << "Error renaming"
    //                     //         << abort(FatalError);
    //                     // }
    //                 }
    //             }
    //         }
    //     }
    // }


    // if (runTime().outputTime() && toroidalPulleys_.size())
    // {
    //     // Write moving pulleys
    //     forAll(toroidalPulleys_, pulleyI)
    //     {
    //         if (toroidalPulleys_[pulleyI].moving())
    //         {
    //             // Create directory if does not exist.
    //             fileName vtkDir(runTime().path()/"VTK");
    //             mkDir(vtkDir);

    //             if (toroidalPulleys_[pulleyI].stlModel().size())
    //             {
    //                 fileName toroidalPulleysDir
    //                 (
    //                     vtkDir/"toroidalPulleys"
    //                 );
    //                 mkDir(toroidalPulleysDir);

    //                 OStringStream toroidalPulleyNvtk;
    //                 toroidalPulleyNvtk()
    //                     << "toroidalPulley-" << pulleyI
    //                     << "_" << runTime().timeIndex() << ".vtk";

    //                 fileName vtkFileName
    //                 (
    //                     // "toroidalPulleys"/
    //                     toroidalPulleysDir/
    //                     word(toroidalPulleyNvtk.str())
    //                 );

    //                 OStringStream toroidalPulleyNvtp;
    //                 toroidalPulleyNvtp()
    //                     << "toroidalPulley-" << pulleyI
    //                     << "_" << runTime().timeIndex() << ".vtp";

    //                 fileName vtpFileName
    //                 (
    //                     word(toroidalPulleyNvtp.str())
    //                 );

    //                 vectorField oldPoints =
    //                     toroidalPulleys_[pulleyI].stlModel().points();

    //                 vectorField newPoints =
    //                     oldPoints
    //                   + toroidalPulleys_[pulleyI].temporalOriginDisplacement();

    //                 toroidalPulleys_[pulleyI].stlModel().movePoints(newPoints);
    //                 toroidalPulleys_[pulleyI].stlModel().write(vtkFileName);
    //                 toroidalPulleys_[pulleyI].stlModel().movePoints(oldPoints);

    //                 // Write temporal collection data
    //                 {
    //                     OStringStream toroidalPulleyNpvd;
    //                     toroidalPulleyNpvd() << "toroidalPulley-"
    //                                         << pulleyI << ".pvd";

    //                     fileName collectionFileName
    //                     (
    //                         toroidalPulleysDir/toroidalPulleyNpvd.str()
    //                     );
    //                     ifstream collectionFile(collectionFileName);
    //                     // IFstream collectionFile(collectionFileName);

    //                     fileName newCollectionFileName
    //                     (
    //                         toroidalPulleysDir/"new.pvd"
    //                     );

    //                     if (collectionFile.good())
    //                     {
    //                         ofstream newCollectionFile
    //                         (
    //                             newCollectionFileName
    //                         );

    //                         // Add to existing collectin file
    //                         label lineIndex = 0;
    //                         do
    //                         {
    //                             lineIndex++;

    //                             std::string line;
    //                             std::getline(collectionFile, line);

    //                             // collectionFile.getLine(line);
    //                             // newCollectionFile << line << '\n';
    //                             // Info<< line << endl;

    //                             if
    //                             (
    //                                 (lineIndex < 3)
    //                              || (
    //                                     line.find("DataSet")
    //                                  != std::string::npos
    //                                 )
    //                             )
    //                             {
    //                                 newCollectionFile << line << '\n';
    //                             }
    //                         }
    //                         while(!collectionFile.eof());

    //                         // Add current pulley data
    //                         newCollectionFile
    //                             << "    <DataSet timestep=\""
    //                             << runTime().value()
    //                             << "\" group=\"\" part=\"0\" file=\""
    //                             << vtpFileName << "\"/>" << '\n';

    //                         // Add last two lines
    //                         newCollectionFile << "  </Collection>" << '\n';
    //                         newCollectionFile << "</VTKFile>";
    //                     }
    //                     else
    //                     {
    //                         OFstream newCollectionFile
    //                         (
    //                             newCollectionFileName
    //                         );

    //                         // Add first two lines
    //                         newCollectionFile
    //                             << "<VTKFile type=\"Collection\" version=\"0.1\" "
    //                             << "byte_order=\"LittleEndian\">" << endl;
    //                         newCollectionFile << "  <Collection>" << endl;

    //                         // Add current pulley data
    //                         newCollectionFile
    //                             << "    <DataSet timestep=\""
    //                             << runTime().value()
    //                             << "\" group=\"\" part=\"0\" file="
    //                             << vtpFileName << "/>" << endl;

    //                         // Add last two lines
    //                         newCollectionFile << "  </Collection>" << endl;
    //                         newCollectionFile << "</VTKFile>";
    //                     }

    //                     mv(newCollectionFileName, collectionFileName);

    //                     // collectionFile.close();
    //                     // newCollectionFile.close();

    //                     // if
    //                     // (
    //                     //     std::rename
    //                     //     (
    //                     //         collectionFile.name().c_str(),
    //                     //         newCollectionFile.name().c_str()
    //                     //     )
    //                     // )
    //                     // {
    //                     //     FatalErrorIn("void Foam::beamModel::writeFields()")
    //                     //         << "Error renaming"
    //                     //         << abort(FatalError);
    //                     // }
    //                 }
    //             }
    //         }
    //     }
    // }
}

const Foam::labelList& Foam::beamModel::upperNeiCellFaces() const
{
    if (!upperNeiCellFacesPtr_)
    {
        makeUpperLowerNeiCellFaces();
    }

    return *upperNeiCellFacesPtr_;
}

const Foam::labelList& Foam::beamModel::lowerNeiCellFaces() const
{
    if (!lowerNeiCellFacesPtr_)
    {
        makeUpperLowerNeiCellFaces();
    }

    return *lowerNeiCellFacesPtr_;
}


Foam::label Foam::beamModel::whichBeam(const label cellIndex) const
{
    // cellIndex is global cell index in parallel simulations

    // label nBeams = nBeamCells_.size();
    // Pout << "nBeams: " << nBeams << endl;

    label nCells = sum(nBeamCells_);
    // label nCells = (startCells_[nBeams-1] + nBeamCells_[nBeams-1]);

    if (cellIndex >= nCells)
    // if (cellIndex >= mesh().nCells())
    {
        FatalErrorIn
        (
            "beamModel::whichBeam(const label cellIndex) const"
        )   << "given label greater than the number of cells, "
            << cellIndex << ">=" << nCells
            << abort(FatalError);
    }

    if (cellIndex < 0)
    {
        FatalErrorIn
        (
            "beamModel::whichBeam(const label cellIndex) const"
        )   << "negativ global cell index "
            << cellIndex << abort(FatalError);

        // return -1;
    }

    forAll (startCells_, bI)
    {
        // const cellZone& cz = mesh().cellZones()[bI];

        if
        (
            cellIndex >= startCells_[bI]
         && cellIndex < startCells_[bI] + nBeamCells_[bI]
         // && cellIndex < startCells_[bI] + cz.size()
        )
        {
            return bI;
        }
    }

    FatalErrorIn
    (
        "beamModel::whichBeam(const label cellIndex) const"
    )  << "Cannot find cell " << cellIndex << " in any of the beams"
       << abort(FatalError);

    return -1;
}


Foam::label Foam::beamModel::whichSegment(const label cellIndex) const
{
    label bI = whichBeam(cellIndex);

    return cellIndex - startCells_[bI];
}


Foam::label Foam::beamModel::whichCell
(
    const label beamIndex,
    const label segIndex
) const
{
    // if (Pstream::parRun())
    // {
    //     // Return global cell index
    // }

    return startCells_[beamIndex] + segIndex;
}


// Foam::label Foam::beamModel::localCellIndex
// (
//     const label globalCellIndex
// ) const
// {
//     if (Pstream::parRun())
//     {
//         label lci =
//             globalToLocalCellAddressing_[Pstream::myProcNo()][globalCellIndex];

//         if (lci != -1)
//         {
//             lci -= csrAddr().globalNCellsOffset();
//         }

//         return lci;
//     }

//     return globalCellIndex;
// }


Foam::labelPair Foam::beamModel::procLocalCellIndex
(
    const label globalCellIndex
) const
{
    if (Pstream::parRun())
    {
        forAll(globalToLocalCellAddressing_, procI)
        {
            if (globalToLocalCellAddressing_[procI][globalCellIndex] != -1)
            {
                return
                    labelPair
                    (
                        procI,
                        globalToLocalCellAddressing_[procI][globalCellIndex]
                    );
            }
        }
    }

    return labelPair(0, globalCellIndex);
}


Foam::label Foam::beamModel::globalCellIndex
(
    const label localCellIndex
) const
{
    if (Pstream::parRun())
    {
        return localToGlobalCellAddressing_[localCellIndex];
    }

    return localCellIndex;
}


Foam::scalar Foam::beamModel::deltaT() const
{
    scalar newDeltaT = runTime().deltaT().value();

    if (deltaTseries_.size())
    {
        newDeltaT = deltaTseries_(runTime().value());
    }

    return newDeltaT;
}


Foam::tmp<Foam::labelField> Foam::beamModel::globalPointsIndices
(
    const label bI
) const
{
    const fvMesh& mesh = this->mesh();

    label nCellZones = mesh.cellZones().size();

    if (nCellZones < 2)
    {
        label nPoints = mesh.nCells() + 1;

        tmp<labelField> tGlPtIndices
        (
            new labelField(nPoints, -1)
        );

        forAll(tGlPtIndices(), pI)
        {
            tGlPtIndices.ref()[pI] = pI;
        }

        if (Pstream::parRun())
        {
            // Apply offset
            tGlPtIndices.ref() += globalCellIndex(0);
        }

        return tGlPtIndices;
    }
    else
    {
        const cellZone& cz = mesh.cellZones()[bI];
        label nPoints = cz.size() + 1;

        tmp<labelField> tGlPtIndices
        (
            new labelField(nPoints, -1)
        );

        for (label i=0; i<cz.size(); i++)
        {
            label locCellIndex = cz[i];
            label glCellIndex = globalCellIndex(locCellIndex);

            label glSegIndex = whichSegment(glCellIndex);
            tGlPtIndices.ref()[i] = glSegIndex;
        }
        // add last point
        tGlPtIndices.ref()[nPoints-1] =
            tGlPtIndices.ref()[nPoints-2] + 1;

        return tGlPtIndices;
    }
}


Foam::tmp<Foam::labelField> Foam::beamModel::localPointsIndices
(
    const label bI
) const
{
    const fvMesh& mesh = this->mesh();

    label nCellZones = mesh.cellZones().size();

    if (nCellZones < 2)
    {
        label nCells = mesh.nCells();
        reduce(nCells, sumOp<label>());

        label nPoints = nCells + 1;

        tmp<labelField> tLocPtIndices
        (
            new labelField(nPoints, -1)
        );

        label globalCellIndexOffset = 0;
        if (Pstream::parRun())
        {
            // Get global cell index offset
            globalCellIndexOffset = globalCellIndex(0);
        }

        label curGlPtIndex = -1;
        for (label i=0; i<mesh.nCells(); i++)
        {
            curGlPtIndex = globalCellIndexOffset + i;
            tLocPtIndices.ref()[curGlPtIndex] = i;
        }
        curGlPtIndex++;
        tLocPtIndices.ref()[curGlPtIndex] = mesh.nCells();

        return tLocPtIndices;
    }
    else
    {
        const cellZone& cz = mesh.cellZones()[bI];

        label nCells = cz.size();
        reduce(nCells, sumOp<label>());

        label nPoints = nCells + 1;

        tmp<labelField> tLocPtIndices
        (
            new labelField(nPoints, -1)
        );

        label glSegIndex = -1;
        for (label i=0; i<cz.size(); i++)
        {
            label locCellIndex = cz[i];
            label glCellIndex = globalCellIndex(locCellIndex);

            glSegIndex = whichSegment(glCellIndex);
            tLocPtIndices.ref()[glSegIndex] = i;
        }

        // add last point
        glSegIndex++;
        tLocPtIndices.ref()[glSegIndex] = cz.size();


        // if (Pstream::myProcNo() == 0)
        // {
        //     if (bI == 0)
        //     {
        //         forAll(tLocPtIndices(), gpI)
        //         {
        //             Pout << gpI << ", " << tLocPtIndices()[gpI] << endl;
        //         }
        //     }
        // }

        return tLocPtIndices;
    }
}

// ************************************************************************* //
