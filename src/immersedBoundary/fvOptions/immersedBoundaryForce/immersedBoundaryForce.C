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

#include "immersedBoundaryForce.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "pimpleControl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(immersedBoundaryForce, 0);
    addToRunTimeSelectionTable(option, immersedBoundaryForce, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::pimpleControl& Foam::fv::immersedBoundaryForce::pimple() const
{
    // The PIMPLE and PISO controls are registered as "solutionControl"
    const pimpleControl* pimplePtr =
        mesh_.findObject<pimpleControl>("solutionControl");

    if (!pimplePtr)
    {
        FatalErrorInFunction
            << "The " << typeName << " option " << name_
            << " requires a PIMPLE or PISO solution control on mesh "
            << mesh_.name() << abort(FatalError);
    }

    return *pimplePtr;
}


Foam::scalar Foam::fv::immersedBoundaryForce::rho()
{
    if (!rhoSet_)
    {
        rhoSet_ = true;

        const IOdictionary* transportPropertiesPtr =
            mesh_.findObject<IOdictionary>("transportProperties");

        if (transportPropertiesPtr && transportPropertiesPtr->found("rho"))
        {
            rho_ =
                dimensionedScalar
                (
                    "rho",
                    dimDensity,
                    *transportPropertiesPtr
                ).value();

            Info<< "    Immersed boundary " << name_ << ": density for the "
                << "forces " << rho_ << " from transportProperties" << endl;
        }
        else
        {
            Info<< "    Immersed boundary " << name_ << ": no density given, "
                << "the forces are per unit density" << endl;
        }
    }

    return rho_;
}


void Foam::fv::immersedBoundaryForce::updateBodies(const bool move)
{
    lambda_.primitiveFieldRef() = 0;

    forAll(bodies_, bodyi)
    {
        immersedBody& body = bodies_[bodyi];

        if (move)
        {
            body.move(mesh_.time().value());
        }

        body.addOccupancy(lambda_, surfaceThreshold_);

        Info<< "    Immersed body " << body.name() << ": "
            << returnReduce(body.internalCells().size(), sumOp<label>())
            << " internal cells, "
            << returnReduce(body.surfaceCells().size(), sumOp<label>())
            << " surface cells" << endl;
    }

    lambda_.correctBoundaryConditions();

    // Record the mesh points for which lambda was calculated
    mesh_.setUpToDatePoints(pointsStamp_);
}


void Foam::fv::immersedBoundaryForce::calcForces(const volVectorField& U)
{
    const label timeIndex = mesh_.time().timeIndex();
    const scalar deltaT = mesh_.time().deltaTValue();
    const scalar rho = this->rho();

    // Store the fluid momentum of the previous time step, when the forces are
    // first calculated in this time step
    if (momentumTimeIndex_ != timeIndex)
    {
        momentum0Valid_ =
            momentumTimeIndex_ >= 0 && momentumTimeIndex_ == timeIndex - 1;
        momentum0_ = momentum_;
        momentumTimeIndex_ = timeIndex;
    }

    forAll(bodies_, bodyi)
    {
        const immersedBody& body = bodies_[bodyi];

        force_[bodyi] = body.force(f_, rho);
        torque_[bodyi] = body.torque(f_, rho);

        momentum_[bodyi] = body.momentum(U, lambda_);

        inertia_[bodyi] = Zero;
        if (momentum0Valid_)
        {
            inertia_[bodyi] =
                rho*(momentum_[bodyi] - momentum0_[bodyi])/deltaT;
        }

        Info<< "    Immersed body " << body.name() << ": force "
            << force_[bodyi] << ", torque " << torque_[bodyi]
            << ", inertia of the fluid inside " << inertia_[bodyi] << endl;
    }

    // In a partitioned fluid-solid interaction, the fluid may be solved
    // several times per time step: only the last solution is written
    forceTimeName_ = mesh_.time().timeName();
    forcesPending_ = true;
}


void Foam::fv::immersedBoundaryForce::writeForces()
{
    if (!forcesPending_)
    {
        return;
    }

    forcesPending_ = false;

    if (Pstream::master())
    {
        forAll(bodies_, bodyi)
        {
            OFstream& os = forceFiles_[bodyi];
            const vector& F = force_[bodyi];
            const vector& T = torque_[bodyi];
            const vector& Fi = inertia_[bodyi];

            os  << forceTimeName_ << token::TAB
                << F.x() << token::SPACE << F.y() << token::SPACE << F.z()
                << token::TAB
                << T.x() << token::SPACE << T.y() << token::SPACE << T.z()
                << token::TAB
                << Fi.x() << token::SPACE << Fi.y() << token::SPACE << Fi.z()
                << endl;
        }
    }
}


void Foam::fv::immersedBoundaryForce::setPenalty(const volVectorField& U)
{
    const scalar deltaT = mesh_.time().deltaTValue();
    const bool volumeFraction = (weighting_ == "volumeFraction");
    const scalarField& lambdaI = lambda_.primitiveField();

    kappa_ = 0;

    for (const immersedBody& body : bodies_)
    {
        if (weighting_ == "sharp")
        {
            for (const label celli : body.insideCells())
            {
                kappa_[celli] = penaltyCoeff_/deltaT;
            }
            continue;
        }

        for (const labelList* cellsPtr :
            {&body.internalCells(), &body.surfaceCells()})
        {
            for (const label celli : *cellsPtr)
            {
                const scalar lambda = lambdaI[celli];

                if (volumeFraction)
                {
                    // The penalised velocity is lambda*Ui + (1 - lambda)*U
                    kappa_[celli] =
                        min
                        (
                            penaltyCoeff_,
                            lambda/max(1 - lambda, 1/penaltyCoeff_)
                        )/deltaT;
                }
                else
                {
                    kappa_[celli] = penaltyCoeff_*lambda/deltaT;
                }
            }
        }
    }

    Ui_ = U;

    for (const immersedBody& body : bodies_)
    {
        body.setVelocity(Ui_);
    }
}


void Foam::fv::immersedBoundaryForce::findGhostCells()
{
    // polyMesh::findCell uses the tet base points, whose construction is
    // collective: construct them on every processor before searching
    (void)mesh_.tetBasePtIs();
    (void)mesh_.cellTree();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const pointField& meshPoints = mesh_.points();
    const labelListList& cellPoints = mesh_.cellPoints();

    // Penalised cells: those covered by a body, or whose centre is inside a
    // body for the sharp weighting
    volScalarField inside
    (
        IOobject
        (
            IOobject::scopedName(name_, "inside"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    DynamicList<label> insideCells;
    boolList isPartial(mesh_.nCells(), false);
    for (const immersedBody& body : bodies_)
    {
        if (weighting_ == "sharp")
        {
            for (const label celli : body.insideCells())
            {
                if (inside[celli] < 0.5)
                {
                    inside[celli] = 1;
                    insideCells.append(celli);
                }
            }
        }
        else
        {
            for (const labelList* cellsPtr :
                {&body.internalCells(), &body.surfaceCells()})
            {
                for (const label celli : *cellsPtr)
                {
                    if (inside[celli] < 0.5)
                    {
                        inside[celli] = 1;
                        insideCells.append(celli);
                    }
                }
            }

            for (const label celli : body.surfaceCells())
            {
                isPartial[celli] = true;
            }
        }
    }
    inside.correctBoundaryConditions();

    penalisedOld_.setSize(mesh_.nCells(), false);
    penalisedOld_ = false;
    for (const label celli : insideCells_)
    {
        penalisedOld_[celli] = true;
    }
    insideCells_.transfer(insideCells);

    // Inside cells with a fluid neighbour
    boolList isGhost(mesh_.nCells(), false);
    forAll(nei, facei)
    {
        const bool ownIn = inside[own[facei]] > 0.5;
        const bool neiIn = inside[nei[facei]] > 0.5;
        if (ownIn && !neiIn)
        {
            isGhost[own[facei]] = true;
        }
        else if (neiIn && !ownIn)
        {
            isGhost[nei[facei]] = true;
        }
    }
    forAll(inside.boundaryField(), patchi)
    {
        const fvPatchScalarField& pf = inside.boundaryField()[patchi];
        if (pf.coupled())
        {
            const scalarField nbr(pf.patchNeighbourField());
            const labelUList& faceCells = pf.patch().faceCells();
            forAll(faceCells, i)
            {
                if (inside[faceCells[i]] > 0.5 && nbr[i] < 0.5)
                {
                    isGhost[faceCells[i]] = true;
                }
            }
        }
    }

    // Reconstructed cells: the penalised cells next to a fluid cell, and the
    // partially covered cells
    DynamicList<label> ghostCells;
    for (const label celli : insideCells_)
    {
        if (isGhost[celli] || isPartial[celli])
        {
            ghostCells.append(celli);
        }
    }
    ghostCells_.transfer(ghostCells);

    const label nGhosts = ghostCells_.size();
    ghostDonors_.setSize(nGhosts);
    ghostDonors_ = -1;
    ghostWallPoints_.setSize(nGhosts, Zero);
    ghostNormals_.setSize(nGhosts, Zero);
    ghostWallVelocities_.setSize(nGhosts, Zero);
    ghostDistances_.setSize(nGhosts, Zero);

    const pointField ghostCentres(mesh_.C(), ghostCells_);

    vectorField spans(nGhosts);
    scalarField searchDistSqr(nGhosts);
    forAll(ghostCells_, gi)
    {
        spans[gi] =
            boundBox(meshPoints, cellPoints[ghostCells_[gi]], false).span();
        searchDistSqr[gi] = 4*magSqr(spans[gi]);
    }

    // Nearest surface point of any body
    scalarField bestDistSqr(nGhosts, GREAT);
    forAll(bodies_, bodyi)
    {
        List<pointIndexHit> hits;
        vectorField normals;
        bodies_[bodyi].nearest(ghostCentres, searchDistSqr, hits, normals);

        forAll(hits, gi)
        {
            if (hits[gi].hit())
            {
                const scalar dSqr =
                    magSqr(hits[gi].hitPoint() - ghostCentres[gi]);

                if (dSqr < bestDistSqr[gi])
                {
                    bestDistSqr[gi] = dSqr;
                    ghostWallPoints_[gi] = hits[gi].hitPoint();
                    ghostNormals_[gi] = normals[gi];
                    // Signed distance from the surface, along the outward
                    // normal: negative inside the body
                    ghostDistances_[gi] =
                        (ghostCentres[gi] - hits[gi].hitPoint()) & normals[gi];
                    ghostWallVelocities_[gi] =
                        bodies_[bodyi].velocity
                        (
                            pointField(1, hits[gi].hitPoint())
                        )()[0];
                }
            }
        }
    }

    // Donor fluid cell of each ghost cell: the fluid cell containing an image
    // point along the normal, about a cell width beyond the surface; if it is
    // not on this processor, it is requested from the other processors
    scalarField widths(nGhosts, Zero);
    DynamicList<label> remoteGhosts;
    forAll(ghostCells_, gi)
    {
        const vector& n = ghostNormals_[gi];
        widths[gi] = cmptSum(cmptMultiply(cmptMag(n), spans[gi]));

        if (magSqr(n) > SMALL)
        {
            ghostDonors_[gi] =
                findDonor(ghostWallPoints_[gi], n, widths[gi], inside);

            if (ghostDonors_[gi] < 0 && Pstream::parRun())
            {
                remoteGhosts.append(gi);
            }
        }
    }
    remoteGhosts_.transfer(remoteGhosts);

    const label nProcs = Pstream::nProcs();
    const label myProc = Pstream::myProcNo();
    sendDonors_.clear();
    sendDonors_.setSize(nProcs);
    recvGhosts_.clear();
    recvGhosts_.setSize(nProcs);
    remoteDonorProcs_.setSize(remoteGhosts_.size(), -1);

    if (Pstream::parRun())
    {
        // Requests of every processor
        List<pointField> reqPoints(nProcs);
        List<vectorField> reqNormals(nProcs);
        List<scalarField> reqWidths(nProcs);
        reqPoints[myProc] = pointField(ghostWallPoints_, remoteGhosts_);
        reqNormals[myProc] = vectorField(ghostNormals_, remoteGhosts_);
        reqWidths[myProc] = scalarField(widths, remoteGhosts_);
        Pstream::allGatherList(reqPoints);
        Pstream::allGatherList(reqNormals);
        Pstream::allGatherList(reqWidths);

        // Requests this processor can serve
        List<labelList> canServe(nProcs);
        List<labelList> canServeCells(nProcs);
        forAll(reqPoints, proci)
        {
            if (proci == myProc)
            {
                continue;
            }

            DynamicList<label> reqs;
            DynamicList<label> cells;
            forAll(reqPoints[proci], k)
            {
                const label donori = findDonor
                (
                    reqPoints[proci][k],
                    reqNormals[proci][k],
                    reqWidths[proci][k],
                    inside
                );

                if (donori >= 0)
                {
                    reqs.append(k);
                    cells.append(donori);
                }
            }
            canServe[proci].transfer(reqs);
            canServeCells[proci].transfer(cells);
        }

        // Tell the requesters which requests can be served
        PstreamBuffers offerBufs(UPstream::commsTypes::nonBlocking);
        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci != myProc)
            {
                UOPstream os(proci, offerBufs);
                os << canServe[proci];
            }
        }
        offerBufs.finishedSends();

        // Choose the lowest processor that can serve each request
        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci != myProc)
            {
                UIPstream is(proci, offerBufs);
                labelList offers(is);
                for (const label k : offers)
                {
                    if (remoteDonorProcs_[k] < 0)
                    {
                        remoteDonorProcs_[k] = proci;
                    }
                }
            }
        }

        List<DynamicList<label>> chosen(nProcs);
        forAll(remoteDonorProcs_, k)
        {
            if (remoteDonorProcs_[k] >= 0)
            {
                chosen[remoteDonorProcs_[k]].append(k);
            }
        }

        // Tell the donors which of their offers are taken
        PstreamBuffers chosenBufs(UPstream::commsTypes::nonBlocking);
        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci != myProc)
            {
                recvGhosts_[proci] = chosen[proci];
                UOPstream os(proci, chosenBufs);
                os << recvGhosts_[proci];
            }
        }
        chosenBufs.finishedSends();

        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci != myProc)
            {
                UIPstream is(proci, chosenBufs);
                labelList taken(is);

                // Map the requests to the donor cells found above
                Map<label> requestToCell(2*canServe[proci].size());
                forAll(canServe[proci], i)
                {
                    requestToCell.insert
                    (
                        canServe[proci][i],
                        canServeCells[proci][i]
                    );
                }

                labelList& send = sendDonors_[proci];
                send.setSize(taken.size());
                forAll(taken, i)
                {
                    send[i] = requestToCell[taken[i]];
                }
            }
        }
    }

    label nLocal = 0;
    for (const label donori : ghostDonors_)
    {
        if (donori >= 0)
        {
            ++nLocal;
        }
    }

    label nRemote = 0;
    for (const label proci : remoteDonorProcs_)
    {
        if (proci >= 0)
        {
            ++nRemote;
        }
    }

    reduce(nLocal, sumOp<label>());
    reduce(nRemote, sumOp<label>());
    const label nTotal = returnReduce(nGhosts, sumOp<label>());

    Info<< "    Ghost cells: " << nTotal << ", with a donor on another "
        << "processor: " << nRemote << ", without a donor fluid cell (the "
        << "body velocity is imposed): " << nTotal - nLocal - nRemote << endl;
}


Foam::label Foam::fv::immersedBoundaryForce::findDonor
(
    const point& wallPoint,
    const vector& normal,
    const scalar width,
    const volScalarField& inside
) const
{
    const vectorField& C = mesh_.C();

    for (label k = 0; k < 4; ++k)
    {
        const point imagePoint = wallPoint + (1 + 0.5*k)*width*normal;

        const label donori = mesh_.findCell(imagePoint);

        if
        (
            donori >= 0
         && inside[donori] < 0.5
         && !penalisedOld_[donori]
         && ((C[donori] - wallPoint) & normal) > 0.25*width
        )
        {
            return donori;
        }
    }

    return -1;
}


void Foam::fv::immersedBoundaryForce::gatherDonors(const volVectorField& U)
{
    const vectorField& C = mesh_.C();
    const label nGhosts = ghostCells_.size();

    ghostDonorU_.setSize(nGhosts, Zero);
    ghostDonorC_.setSize(nGhosts, Zero);

    forAll(ghostCells_, gi)
    {
        const label donori = ghostDonors_[gi];
        if (donori >= 0)
        {
            ghostDonorU_[gi] = U[donori];
            ghostDonorC_[gi] = C[donori];
        }
    }

    if (Pstream::parRun())
    {
        const label nProcs = Pstream::nProcs();
        const label myProc = Pstream::myProcNo();

        PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);
        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci != myProc)
            {
                UOPstream os(proci, pBufs);
                os  << vectorField(U.primitiveField(), sendDonors_[proci])
                    << pointField(C, sendDonors_[proci]);
            }
        }
        pBufs.finishedSends();

        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci != myProc)
            {
                UIPstream is(proci, pBufs);
                vectorField donorU(is);
                pointField donorC(is);

                const labelList& ks = recvGhosts_[proci];
                forAll(ks, i)
                {
                    const label gi = remoteGhosts_[ks[i]];
                    ghostDonorU_[gi] = donorU[i];
                    ghostDonorC_[gi] = donorC[i];
                    // Mark the ghost as having a donor
                    ghostDonors_[gi] = -2;
                }
            }
        }
    }
}


void Foam::fv::immersedBoundaryForce::reconstructTargets
(
    const volVectorField& U
)
{
    // Linear velocity profile along the surface normal, through the body
    // velocity Ub on the surface and the donor cell velocity Ud:
    // U(s) = Ub + (Ud - Ub)*s/sd, where s is the signed distance from the
    // surface (negative inside the body) and sd that of the donor centre
    gatherDonors(U);

    // The reconstructed cells are fully penalised: the target is the fluid
    // velocity at their centre, so it is not also weighted by the occupancy
    const scalar rate = penaltyCoeff_/mesh_.time().deltaTValue();

    vectorField& UiI = Ui_.primitiveFieldRef();
    forAll(ghostCells_, gi)
    {
        if (ghostDonors_[gi] != -1)
        {
            const label celli = ghostCells_[gi];
            const vector& Ub = ghostWallVelocities_[gi];
            const scalar sd =
                (ghostDonorC_[gi] - ghostWallPoints_[gi]) & ghostNormals_[gi];
            const scalar ratio = min(ghostDistances_[gi]/sd, scalar(1));

            // Deeper inside than the donor is outside: the body velocity
            if (ratio >= -1)
            {
                UiI[celli] = Ub + ratio*(ghostDonorU_[gi] - Ub);
            }

            kappa_[celli] = rate;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::immersedBoundaryForce::immersedBoundaryForce
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh),
    bodies_(),
    method_("penalty"),
    penaltyCoeff_(1e3),
    weighting_("volumeFraction"),
    reconstruction_("none"),
    kappa_(mesh.nCells(), Zero),
    couplingCoeff_(0.8),
    surfaceThreshold_(1e-4),
    rho_(1),
    rhoSet_(false),
    writeSurfaces_(true),
    lambda_
    (
        IOobject
        (
            IOobject::scopedName(name, "lambda"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    ),
    Ui_
    (
        IOobject
        (
            IOobject::scopedName(name, "Ui"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero)
    ),
    f_
    (
        IOobject
        (
            IOobject::scopedName(name, "f"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity/dimTime, Zero)
    ),
    f0_(mesh.nCells(), Zero),
    timeIndex_(-1),
    pointsStamp_
    (
        IOobject
        (
            IOobject::scopedName(name, "pointsStamp"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    UEqnPtr_(nullptr),
    UEqnTimeIndex_(-1),
    UEqnCorr_(-1),
    momentum_(),
    momentumTimeIndex_(-1),
    momentum0_(),
    momentum0Valid_(false),
    forceFiles_(),
    force_(),
    torque_(),
    inertia_(),
    forceTimeName_(),
    forcesPending_(false),
    warnedMatrix_(false)
{
    fieldNames_.resize(1, coeffs_.getOrDefault<word>("U", "U"));

    fv::option::resetApplied();

    Info<< "    Immersed bodies:" << endl;

    const dictionary& bodiesDict = coeffs_.subDict("bodies");

    label nBodies = 0;
    for (const entry& e : bodiesDict)
    {
        if (e.isDict())
        {
            ++nBodies;
        }
    }

    if (nBodies == 0)
    {
        FatalIOErrorInFunction(bodiesDict)
            << "No immersed bodies are given" << exit(FatalIOError);
    }

    bodies_.setSize(nBodies);

    label bodyi = 0;
    for (const entry& e : bodiesDict)
    {
        if (e.isDict())
        {
            bodies_.set
            (
                bodyi++,
                new immersedBody(e.keyword(), e.dict(), mesh)
            );
        }
    }

    momentum_.setSize(nBodies, Zero);
    momentum0_.setSize(nBodies, Zero);
    force_.setSize(nBodies, Zero);
    torque_.setSize(nBodies, Zero);
    inertia_.setSize(nBodies, Zero);

    read(dict);

    // Force files
    if (Pstream::master())
    {
        fileName dir(mesh.time().globalPath()/"postProcessing");

        if (mesh.name() != polyMesh::defaultRegion)
        {
            dir = dir/mesh.name();
        }

        dir = dir/name/mesh.time().timeName();

        mkDir(dir);

        forceFiles_.setSize(nBodies);

        forAll(bodies_, bodyi)
        {
            forceFiles_.set
            (
                bodyi,
                new OFstream(dir/(bodies_[bodyi].name() + ".dat"))
            );

            forceFiles_[bodyi].precision(10);

            forceFiles_[bodyi]
                << "# Immersed body " << bodies_[bodyi].name() << nl
                << "# Time" << token::TAB << "force (x y z)" << token::TAB
                << "torque (x y z)" << token::TAB
                << "inertia of the fluid inside (x y z)" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::immersedBoundaryForce::~immersedBoundaryForce()
{
    writeForces();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::immersedBoundaryForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const label timeIndex = mesh_.time().timeIndex();

    bool updated = false;

    if (timeIndex != timeIndex_)
    {
        // The previous time step is complete
        writeForces();

        // New time step: store the forcing of the previous time step (or the
        // forcing read on restart), and move the bodies to the new time
        f0_ = f_.primitiveField();

        bool moving = false;
        for (const immersedBody& body : bodies_)
        {
            moving = moving || body.moving();
        }

        if (timeIndex_ < 0 || moving || !mesh_.upToDatePoints(pointsStamp_))
        {
            updateBodies(true);
            updated = true;
        }

        timeIndex_ = timeIndex;

        if (writeSurfaces_ && mesh_.time().writeTime())
        {
            fileName dir(mesh_.time().globalPath()/"postProcessing");

            if (mesh_.name() != polyMesh::defaultRegion)
            {
                dir = dir/mesh_.name();
            }

            dir = dir/name_/"surfaces"/mesh_.time().timeName();

            for (const immersedBody& body : bodies_)
            {
                if (body.moving())
                {
                    body.writeSurface(dir/(body.name() + ".stl"));
                }
            }
        }
    }
    else if (!mesh_.upToDatePoints(pointsStamp_))
    {
        // The mesh has moved since lambda was calculated
        updateBodies(false);
        updated = true;
    }

    if (method_ == "penalty")
    {
        // Implicit volume penalisation: kappa*(Ui - U)
        const volVectorField& U = eqn.psi();

        setPenalty(U);

        if (reconstruction_ == "linear")
        {
            // Collective: updated is the same on every processor
            if (updated)
            {
                findGhostCells();
            }
            reconstructTargets(U);
        }

        volScalarField::Internal kappa
        (
            IOobject
            (
                IOobject::scopedName(name_, "kappa"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless/dimTime, Zero)
        );
        kappa.field() = kappa_;

        eqn += kappa*Ui_();
        eqn -= fvm::Sp(kappa, U);

        return;
    }

    if (pimple().firstIter())
    {
        // Start the time step, or a repeated solution of the time step, e.g.
        // in a partitioned fluid-solid interaction iteration, from the
        // forcing of the previous time step in the covered cells
        f_.primitiveFieldRef() = f0_*lambda_.primitiveField();
    }
    else if (updated)
    {
        // Remove the forcing from the cells no longer covered
        vectorField& fI = f_.primitiveFieldRef();
        const scalarField& lambdaI = lambda_.primitiveField();

        forAll(fI, celli)
        {
            if (lambdaI[celli] <= SMALL)
            {
                fI[celli] = Zero;
            }
        }
    }

    f_.correctBoundaryConditions();

    eqn += f_;
}


void Foam::fv::immersedBoundaryForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    FatalErrorInFunction
        << "The " << typeName << " option is only implemented for "
        << "incompressible flow" << abort(FatalError);
}


void Foam::fv::immersedBoundaryForce::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    UEqnPtr_ = &eqn;
    UEqnTimeIndex_ = mesh_.time().timeIndex();
    UEqnCorr_ = pimple().corr();
}


void Foam::fv::immersedBoundaryForce::correct(volVectorField& U)
{
    const pimpleControl& pimple = this->pimple();

    // Nothing to do after the momentum predictor
    if (pimple.corrPISO() < 1)
    {
        return;
    }

    if (method_ == "penalty")
    {
        // The forcing exerted by the penalisation, for the forces
        f_.primitiveFieldRef() =
            kappa_*(Ui_.primitiveField() - U.primitiveField());
        f_.correctBoundaryConditions();

        if (pimple.corrPISO() == pimple.nCorrPISO() && pimple.finalIter())
        {
            calcForces(U);
        }

        return;
    }

    // Velocity of the bodies in the covered cells, U elsewhere
    Ui_ = U;

    for (const immersedBody& body : bodies_)
    {
        body.setVelocity(Ui_);
    }

    // Forcing increment, in the cells with a nonzero occupancy
    volVectorField::Internal df
    (
        IOobject
        (
            IOobject::scopedName(name_, "df"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector(f_.dimensions(), Zero)
    );

    {
        const scalar coeff = couplingCoeff_/mesh_.time().deltaTValue();
        const scalarField& lambdaI = lambda_.primitiveField();
        const vectorField& UiI = Ui_.primitiveField();
        const vectorField& UI = U.primitiveField();
        vectorField& dfI = df.field();

        forAll(dfI, celli)
        {
            if (lambdaI[celli] > SMALL)
            {
                dfI[celli] = coeff*(UiI[celli] - UI[celli]);
            }
        }
    }

    f_.primitiveFieldRef() += df.field();
    f_.correctBoundaryConditions();

    if (pimple.corrPISO() < pimple.nCorrPISO())
    {
        // Add the increment to the source of the momentum matrix, so that the
        // next pressure corrector uses the updated forcing. The matrix is only
        // used if it is the one assembled for U in this outer corrector; it
        // is not accessed after the last corrector, as the solver may have
        // cleared it
        const bool matrixValid =
            UEqnPtr_
         && UEqnTimeIndex_ == mesh_.time().timeIndex()
         && UEqnCorr_ == pimple.corr()
         && &UEqnPtr_->psi() == &U;

        if (matrixValid)
        {
            // The forcing is on the right-hand side: eqn == f
            *UEqnPtr_ -= df;
        }
        else if (!warnedMatrix_)
        {
            warnedMatrix_ = true;

            WarningInFunction
                << "The momentum matrix of " << U.name() << " was not "
                << "constrained by the " << typeName << " option " << name_
                << " in this outer corrector, so the updated forcing is only "
                << "used from the next outer corrector" << endl;
        }
    }
    else
    {
        // Last pressure corrector
        UEqnPtr_ = nullptr;

        if (pimple.finalIter())
        {
            calcForces(U);
        }
    }
}


bool Foam::fv::immersedBoundaryForce::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
        coeffs_.readIfPresent("method", method_);

        if (method_ != "incremental" && method_ != "penalty")
        {
            FatalIOErrorInFunction(coeffs_)
                << "Unknown method " << method_ << ": valid methods are "
                << "penalty and incremental" << exit(FatalIOError);
        }

        coeffs_.readCheckIfPresent
        (
            "penaltyCoeff",
            penaltyCoeff_,
            scalarMinMax::ge(1)
        );
        coeffs_.readIfPresent("weighting", weighting_);

        if
        (
            weighting_ != "volumeFraction"
         && weighting_ != "occupancy"
         && weighting_ != "sharp"
        )
        {
            FatalIOErrorInFunction(coeffs_)
                << "Unknown weighting " << weighting_ << ": valid weightings "
                << "are volumeFraction, occupancy and sharp"
                << exit(FatalIOError);
        }

        coeffs_.readIfPresent("reconstruction", reconstruction_);

        if (reconstruction_ != "none" && reconstruction_ != "linear")
        {
            FatalIOErrorInFunction(coeffs_)
                << "Unknown reconstruction " << reconstruction_
                << ": valid reconstructions are none and linear"
                << exit(FatalIOError);
        }

        coeffs_.readCheckIfPresent
        (
            "couplingCoeff",
            couplingCoeff_,
            scalarMinMax::ge(SMALL)
        );

        // An occupancy threshold of 0.5 or more would classify every covered
        // cell as internal
        coeffs_.readCheckIfPresent
        (
            "surfaceThreshold",
            surfaceThreshold_,
            scalarMinMax(0, 0.5 - SMALL)
        );
        coeffs_.readIfPresent("writeSurfaces", writeSurfaces_);

        if (coeffs_.readIfPresent("rho", rho_))
        {
            rhoSet_ = true;
        }

        if (method_ == "penalty")
        {
            Info<< "    Penalty method: coefficient " << penaltyCoeff_
                << ", " << weighting_ << " weighting";
        }
        else
        {
            Info<< "    Incremental method: coupling coefficient "
                << couplingCoeff_;
        }
        Info<< ", surface threshold " << surfaceThreshold_ << endl;

        return true;
    }

    return false;
}


// ************************************************************************* //
