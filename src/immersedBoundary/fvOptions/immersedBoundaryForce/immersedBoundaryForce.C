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
#include "interpolationCellPoint.H"
#include "DynamicField.H"
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


const Foam::label Foam::fv::immersedBoundaryForce::nImageLevels_ = 4;

const Foam::Enum<Foam::fv::immersedBoundaryForce::forcingMethod>
Foam::fv::immersedBoundaryForce::methodNames_
({
    {incremental, "incremental"},
    {penalty, "penalty"},
    {ghostCell, "ghostCell"}
});

const Foam::Enum<Foam::fv::immersedBoundaryForce::penaltyWeighting>
Foam::fv::immersedBoundaryForce::weightingNames_
({
    {volumeFraction, "volumeFraction"},
    {occupancy, "occupancy"}
});

const Foam::Enum<Foam::fv::immersedBoundaryForce::imageInterpolationType>
Foam::fv::immersedBoundaryForce::imageInterpolationNames_
({
    {cellPointImage, "cellPoint"},
    {cellImage, "cell"}
});


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


Foam::fileName Foam::fv::immersedBoundaryForce::outputDir() const
{
    fileName dir(mesh_.time().globalPath()/"postProcessing");

    if (mesh_.name() != polyMesh::defaultRegion)
    {
        dir /= mesh_.name();
    }

    return dir/name_;
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


Foam::scalar Foam::fv::immersedBoundaryForce::nu()
{
    if (!nuSet_)
    {
        const IOdictionary* transportPropertiesPtr =
            mesh_.findObject<IOdictionary>("transportProperties");

        if (!transportPropertiesPtr || !transportPropertiesPtr->found("nu"))
        {
            FatalIOErrorInFunction(coeffs_)
                << "The " << typeName << " option " << name_ << " requires "
                << "the kinematic viscosity for the surface rate: give nu in "
                << "the option, or set surfaceRateCoeff 0" << exit(FatalIOError);
        }

        nu_ =
            dimensionedScalar
            (
                "nu",
                dimViscosity,
                *transportPropertiesPtr
            ).value();
        if (nu_ < 0)
        {
            FatalIOErrorInFunction(*transportPropertiesPtr)
                << "Kinematic viscosity nu must be non-negative, found "
                << nu_ << exit(FatalIOError);
        }
        nuSet_ = true;
    }

    return nu_;
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
    const bool useVolumeFraction = (weighting_ == volumeFraction);
    const scalarField& lambdaI = lambda_.primitiveField();

    Ui_ = U;

    for (const immersedBody& body : bodies_)
    {
        body.setVelocity(Ui_);
    }

    const scalar nu =
        (useVolumeFraction && surfaceRateCoeff_ > 0 ? this->nu() : 0);

    kappa_.primitiveFieldRef() = 0;
    scalarField& kappaI = kappa_.primitiveFieldRef();

    for (const immersedBody& body : bodies_)
    {
        for (const labelList* cellsPtr :
            {&body.internalCells(), &body.surfaceCells()})
        {
            for (const label celli : *cellsPtr)
            {
                const scalar lambda = lambdaI[celli];

                if (useVolumeFraction)
                {
                    // The penalised velocity is lambda*Ui + (1 - lambda)*U
                    kappaI[celli] =
                        min
                        (
                            penaltyCoeff_,
                            lambda/max(1 - lambda, 1/penaltyCoeff_)
                        )/deltaT;

                    // Rate independent of the time step in the partially
                    // covered cells
                    if (surfaceRateCoeff_ > 0 && lambda < 1 - surfaceThreshold_)
                    {
                        const scalar w = body.cellWidths()[celli];

                        kappaI[celli] = min
                        (
                            kappaI[celli],
                            lambda/(1 - lambda)*surfaceRateCoeff_
                           *(nu/sqr(w) + max(mag(Ui_[celli]), mag(U[celli]))/w)
                        );
                    }
                }
                else
                {
                    kappaI[celli] = penaltyCoeff_*lambda/deltaT;
                }
            }
        }
    }
}


void Foam::fv::immersedBoundaryForce::setGhostCellPenalty
(
    const volVectorField& U
)
{
    const scalar deltaT = mesh_.time().deltaTValue();

    Ui_ = U;

    kappa_.primitiveFieldRef() = 0;
    scalarField& kappaI = kappa_.primitiveFieldRef();

    // Sharp penalisation of the cells whose centre is inside a body, towards
    // the body velocity
    for (const immersedBody& body : bodies_)
    {
        body.setVelocity(Ui_);

        for (const label celli : body.insideCells())
        {
            kappaI[celli] = penaltyCoeff_/deltaT;
        }
    }

    // Reconstructed target velocity in the ghost cells
    reconstructTargets(U);
}


Foam::vector Foam::fv::immersedBoundaryForce::solvedSpan
(
    const label celli
) const
{
    vector span
    (
        boundBox(mesh_.points(), mesh_.cellPoints()[celli], false).span()
    );

    const Vector<label>& solD = mesh_.solutionD();
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        if (solD[d] != 1)
        {
            span[d] = 0;
        }
    }

    return span;
}


void Foam::fv::immersedBoundaryForce::findGhostCells(const bool newTimeStep)
{
    // polyMesh::findCell uses the tet base points, whose construction is
    // collective: construct them on every processor before searching
    (void)mesh_.tetBasePtIs();
    (void)mesh_.cellTree();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    // Penalised cells: those whose centre is inside a body
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

    // The body containing each penalised cell (the first, if several)
    labelList cellBody(mesh_.nCells(), -1);

    DynamicList<label> insideCells;
    forAll(bodies_, bodyi)
    {
        for (const label celli : bodies_[bodyi].insideCells())
        {
            if (inside[celli] < 0.5)
            {
                inside[celli] = 1;
                cellBody[celli] = bodyi;
                insideCells.append(celli);
            }
        }
    }
    inside.correctBoundaryConditions();

    // Cells penalised at the end of the previous time step, which are not
    // donors; when the bodies are updated again within a time step (e.g. the
    // mesh moves), these are kept. On startup or restart, insideCells_ holds
    // the cells inside the bodies at the start time (set in addSup). The mesh
    // topology is assumed not to change
    if (newTimeStep || penalisedOld_.size() != mesh_.nCells())
    {
        penalisedOld_.setSize(mesh_.nCells());
        penalisedOld_ = false;
        for (const label celli : insideCells_)
        {
            if (celli < mesh_.nCells())
            {
                penalisedOld_[celli] = true;
            }
        }
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

    // Ghost cells: the penalised cells next to a fluid cell
    DynamicList<label> ghostCells;
    for (const label celli : insideCells_)
    {
        if (isGhost[celli])
        {
            ghostCells.append(celli);
        }
    }
    ghostCells_.transfer(ghostCells);

    const label nGhosts = ghostCells_.size();
    ghostDonors_.setSize(nGhosts);
    ghostDonors_ = -1;
    ghostRemoteDonors_.setSize(nGhosts);
    ghostRemoteDonors_ = false;
    ghostWallPoints_.setSize(nGhosts);
    ghostWallPoints_ = Zero;
    ghostNormals_.setSize(nGhosts);
    ghostNormals_ = Zero;
    ghostWallVelocities_.setSize(nGhosts);
    ghostWallVelocities_ = Zero;
    ghostDistances_.setSize(nGhosts);
    ghostDistances_ = Zero;
    ghostImagePoints_.setSize(nGhosts);
    ghostImagePoints_ = Zero;

    // Cell spans in the solved directions, so that the empty direction of a
    // 2D mesh does not set the search radius or the tolerances
    vectorField spans(nGhosts);
    forAll(ghostCells_, gi)
    {
        spans[gi] = solvedSpan(ghostCells_[gi]);
    }

    // Ghost cells of each body
    List<DynamicList<label>> bodyGhosts(bodies_.size());
    forAll(ghostCells_, gi)
    {
        bodyGhosts[cellBody[ghostCells_[gi]]].append(gi);
    }

    // Nearest surface point of the body that contains each ghost cell
    forAll(bodies_, bodyi)
    {
        const labelList& ghosts = bodyGhosts[bodyi];

        const pointField centres(mesh_.C(), labelList(ghostCells_, ghosts));
        scalarField searchDistSqr(ghosts.size());
        forAll(ghosts, i)
        {
            searchDistSqr[i] = 4*magSqr(spans[ghosts[i]]);
        }

        List<pointIndexHit> hits;
        vectorField normals;
        bodies_[bodyi].nearest(centres, searchDistSqr, hits, normals);

        DynamicList<label> hitGhosts(ghosts.size());
        DynamicField<point> hitPoints(ghosts.size());

        forAll(ghosts, i)
        {
            const label gi = ghosts[i];

            if (!hits[i].hit())
            {
                continue;
            }

            ghostWallPoints_[gi] = hits[i].hitPoint();

            // Outward normal along the line from the cell centre, which is
            // inside the body, to the surface point, so that it varies
            // continuously as the body moves; the normal of the nearest face
            // if the centre is on the surface
            const vector r(centres[i] - hits[i].hitPoint());
            const scalar magR = mag(r);
            if (magR > 1e-3*mag(spans[gi]))
            {
                ghostNormals_[gi] = -r/magR;
            }
            else
            {
                ghostNormals_[gi] = normals[i];
            }

            // Signed distance from the surface, along the outward normal:
            // negative inside the body
            ghostDistances_[gi] = r & ghostNormals_[gi];

            hitGhosts.append(gi);
            hitPoints.append(hits[i].hitPoint());
        }

        const vectorField Ub(bodies_[bodyi].velocity(hitPoints));
        forAll(hitGhosts, i)
        {
            ghostWallVelocities_[hitGhosts[i]] = Ub[i];
        }

        // A ghost cell without a surface point has the body velocity
        const label nMissed =
            returnReduce(ghosts.size() - hitGhosts.size(), sumOp<label>());
        if (nMissed)
        {
            WarningInFunction
                << "Immersed body " << bodies_[bodyi].name() << ": no "
                << "surface point found within two cell sizes of " << nMissed
                << " ghost cells, which have the body velocity" << endl;
        }
    }

    // Donor fluid cell of each ghost cell: the fluid cell containing the
    // first valid image point along the normal (level 0, 1, ...), about
    // imageDistance cell widths beyond the surface. The image points of a
    // level are searched on all the processors before the next level, so
    // that the donor does not depend on the decomposition: a ghost cell
    // without a local donor at level 0 is also requested from the other
    // processors
    scalarField widths(nGhosts, Zero);
    labelList localLevels(nGhosts, nImageLevels_);
    DynamicList<label> remoteGhosts;
    forAll(ghostCells_, gi)
    {
        const vector& n = ghostNormals_[gi];
        widths[gi] = immersedBody::normalWidth(n, spans[gi]);

        if (magSqr(n) > SMALL)
        {
            ghostDonors_[gi] = findDonor
            (
                ghostWallPoints_[gi],
                n,
                widths[gi],
                inside,
                ghostImagePoints_[gi],
                localLevels[gi]
            );

            if (localLevels[gi] > 0 && Pstream::parRun())
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
    sendImagePoints_.clear();
    sendImagePoints_.setSize(nProcs);
    recvGhosts_.clear();
    recvGhosts_.setSize(nProcs);
    remoteDonorProcs_.setSize(remoteGhosts_.size());
    remoteDonorProcs_ = -1;

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
        List<labelList> canServeLevels(nProcs);
        List<labelList> canServeCells(nProcs);
        List<pointField> canServeImages(nProcs);
        forAll(reqPoints, proci)
        {
            if (proci == myProc)
            {
                continue;
            }

            DynamicList<label> reqs;
            DynamicList<label> levels;
            DynamicList<label> cells;
            DynamicList<point> images;
            forAll(reqPoints[proci], k)
            {
                point imagePoint;
                label level;
                const label donori = findDonor
                (
                    reqPoints[proci][k],
                    reqNormals[proci][k],
                    reqWidths[proci][k],
                    inside,
                    imagePoint,
                    level
                );

                if (donori >= 0)
                {
                    reqs.append(k);
                    levels.append(level);
                    cells.append(donori);
                    images.append(imagePoint);
                }
            }
            canServe[proci].transfer(reqs);
            canServeLevels[proci].transfer(levels);
            canServeCells[proci].transfer(cells);
            canServeImages[proci].transfer(images);
        }

        // Tell the requesters which requests can be served
        PstreamBuffers offerBufs(UPstream::commsTypes::nonBlocking);
        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci != myProc)
            {
                UOPstream os(proci, offerBufs);
                os << canServe[proci] << canServeLevels[proci];
            }
        }
        offerBufs.finishedSends();

        // Choose the offer with the lowest level, then the lowest
        // processor, if its level is lower than that of the local donor
        labelList bestLevels(remoteGhosts_.size());
        forAll(remoteGhosts_, k)
        {
            bestLevels[k] = localLevels[remoteGhosts_[k]];
        }

        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci != myProc)
            {
                UIPstream is(proci, offerBufs);
                labelList offers(is);
                labelList offerLevels(is);
                forAll(offers, i)
                {
                    const label k = offers[i];
                    if (offerLevels[i] < bestLevels[k])
                    {
                        bestLevels[k] = offerLevels[i];
                        remoteDonorProcs_[k] = proci;
                    }
                }
            }
        }

        // The ghost cells with a remote donor do not use the local one
        forAll(remoteGhosts_, k)
        {
            if (remoteDonorProcs_[k] >= 0)
            {
                ghostDonors_[remoteGhosts_[k]] = -1;
                ghostRemoteDonors_[remoteGhosts_[k]] = true;
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
                Map<label> requestToIndex(2*canServe[proci].size());
                forAll(canServe[proci], i)
                {
                    requestToIndex.insert(canServe[proci][i], i);
                }

                labelList& send = sendDonors_[proci];
                pointField& sendImages = sendImagePoints_[proci];
                send.setSize(taken.size());
                sendImages.setSize(taken.size());
                forAll(taken, i)
                {
                    const label j = requestToIndex[taken[i]];
                    send[i] = canServeCells[proci][j];
                    sendImages[i] = canServeImages[proci][j];
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
    const volScalarField& inside,
    point& imagePoint,
    label& level
) const
{
    const vectorField& C = mesh_.C();

    for (label k = 0; k < nImageLevels_; ++k)
    {
        imagePoint = wallPoint + (imageDistance_ + 0.5*k)*width*normal;

        const label donori = mesh_.findCell(imagePoint);

        if
        (
            donori >= 0
         && inside[donori] < 0.5
         && !penalisedOld_[donori]
         && ((C[donori] - wallPoint) & normal) > 0.25*width
        )
        {
            level = k;
            return donori;
        }
    }

    level = nImageLevels_;
    return -1;
}


void Foam::fv::immersedBoundaryForce::gatherDonors(const volVectorField& U)
{
    const vectorField& C = mesh_.C();
    const label nGhosts = ghostCells_.size();
    const bool cellPoint = (imageInterpolation_ == cellPointImage);

    // The interpolation is only constructed on the processors with local or
    // requested donors
    bool hasDonors = false;
    for (const label donori : ghostDonors_)
    {
        hasDonors = hasDonors || donori >= 0;
    }
    for (const labelList& donors : sendDonors_)
    {
        hasDonors = hasDonors || donors.size();
    }

    autoPtr<interpolationCellPoint<vector>> interpPtr;
    if (cellPoint && hasDonors)
    {
        interpPtr.reset(new interpolationCellPoint<vector>(U));
    }

    // Velocity and position used for a donor cell and its image point
    auto donorValue = [&](const label donori, const point& imagePoint)
    {
        return
        (
            cellPoint
          ? interpPtr->interpolate(imagePoint, donori)
          : U[donori]
        );
    };
    auto donorPoint = [&](const label donori, const point& imagePoint)
    {
        return (cellPoint ? imagePoint : C[donori]);
    };

    ghostDonorU_.setSize(nGhosts, Zero);
    ghostDonorC_.setSize(nGhosts, Zero);

    forAll(ghostCells_, gi)
    {
        const label donori = ghostDonors_[gi];
        if (donori >= 0)
        {
            ghostDonorU_[gi] = donorValue(donori, ghostImagePoints_[gi]);
            ghostDonorC_[gi] = donorPoint(donori, ghostImagePoints_[gi]);
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
                const labelList& donors = sendDonors_[proci];
                const pointField& images = sendImagePoints_[proci];

                vectorField values(donors.size());
                pointField points(donors.size());
                forAll(donors, i)
                {
                    values[i] = donorValue(donors[i], images[i]);
                    points[i] = donorPoint(donors[i], images[i]);
                }

                UOPstream os(proci, pBufs);
                os  << values << points;
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

    vectorField& UiI = Ui_.primitiveFieldRef();
    forAll(ghostCells_, gi)
    {
        if (ghostDonors_[gi] >= 0 || ghostRemoteDonors_[gi])
        {
            // The ghost cell centre is inside the body (s <= 0) and the
            // donor accepted at more than a quarter cell width outside it
            // (sd > 0), so the ratio is not positive
            const label celli = ghostCells_[gi];
            const vector& Ub = ghostWallVelocities_[gi];
            const scalar sd =
                (ghostDonorC_[gi] - ghostWallPoints_[gi]) & ghostNormals_[gi];
            const scalar ratio = ghostDistances_[gi]/sd;

            // Deeper inside than the donor is outside: the body velocity,
            // as set by setPenalty
            if (ratio >= -1)
            {
                UiI[celli] = Ub + ratio*(ghostDonorU_[gi] - Ub);
            }
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
    method_(penalty),
    penaltyCoeff_(1e3),
    weighting_(volumeFraction),
    surfaceRateCoeff_(3),
    nu_(0),
    nuSet_(false),
    kappa_
    (
        IOobject
        (
            IOobject::scopedName(name, "kappa"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero)
    ),
    imageDistance_(1.5),
    imageInterpolation_(cellPointImage),
    insideCells_(),
    penalisedOld_(),
    ghostCells_(),
    ghostDonors_(),
    ghostRemoteDonors_(),
    ghostWallPoints_(),
    ghostNormals_(),
    ghostWallVelocities_(),
    ghostDistances_(),
    remoteGhosts_(),
    remoteDonorProcs_(),
    ghostImagePoints_(),
    sendDonors_(),
    sendImagePoints_(),
    recvGhosts_(),
    ghostDonorU_(),
    ghostDonorC_(),
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
        const fileName dir(outputDir()/mesh.time().timeName());

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
    const bool newTimeStep = timeIndex != timeIndex_;

    bool updated = false;

    if (newTimeStep)
    {
        // The previous time step is complete
        writeForces();

        // On startup or restart, recover the previous-time fluid momentum
        // from the loaded velocity and the body configuration at startTime.
        // The body is moved to the new time by updateBodies() below.
        if (momentumTimeIndex_ < 0)
        {
            lambda_.primitiveFieldRef() = 0;

            for (immersedBody& body : bodies_)
            {
                body.addOccupancy(lambda_, surfaceThreshold_);
            }

            lambda_.correctBoundaryConditions();

            forAll(bodies_, bodyi)
            {
                momentum_[bodyi] =
                    bodies_[bodyi].momentum(eqn.psi(), lambda_);
            }

            momentum0_ = momentum_;
            momentum0Valid_ = true;
            momentumTimeIndex_ = timeIndex;

            // The cells inside the bodies at startTime, which an
            // uninterrupted run penalised in the previous time step and
            // which are therefore not ghost cell donors in this one
            if (method_ == ghostCell)
            {
                DynamicList<label> insideCells;
                for (const immersedBody& body : bodies_)
                {
                    insideCells.append(body.insideCells());
                }
                insideCells_.transfer(insideCells);
            }
        }

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
            const fileName dir
            (
                outputDir()/"surfaces"/mesh_.time().timeName()
            );

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

    if (method_ == penalty || method_ == ghostCell)
    {
        // Implicit volume penalisation: kappa*(Ui - U)
        const volVectorField& U = eqn.psi();

        if (method_ == ghostCell)
        {
            // Collective: updated is the same on every processor
            if (updated)
            {
                findGhostCells(newTimeStep);
            }
            setGhostCellPenalty(U);
        }
        else
        {
            setPenalty(U);
        }

        eqn += kappa_()*Ui_();
        eqn -= fvm::Sp(kappa_(), U);

        return;
    }

    if (newTimeStep || pimple().firstIter())
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

    if (method_ == penalty || method_ == ghostCell)
    {
        // The forcing exerted by the penalisation, for the forces
        f_.primitiveFieldRef() =
            kappa_.primitiveField()*(Ui_.primitiveField() - U.primitiveField());
        f_.correctBoundaryConditions();

        if
        (
            pimple.corrPISO() == pimple.nCorrPISO()
         && (pimple.finalIter() || pimple.corr() == 0)
        )
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

        // pisoControl does not enter an outer-corrector loop, so corr() stays
        // zero and firstIter()/finalIter() are both false
        if (pimple.finalIter() || pimple.corr() == 0)
        {
            calcForces(U);
        }
    }
}


bool Foam::fv::immersedBoundaryForce::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
        methodNames_.readIfPresent("method", coeffs_, method_);

        coeffs_.readCheckIfPresent
        (
            "penaltyCoeff",
            penaltyCoeff_,
            scalarMinMax::ge(1)
        );
        weightingNames_.readIfPresent("weighting", coeffs_, weighting_);
        coeffs_.readCheckIfPresent
        (
            "surfaceRateCoeff",
            surfaceRateCoeff_,
            scalarMinMax::ge(0)
        );
        if (coeffs_.readCheckIfPresent("nu", nu_, scalarMinMax::ge(0)))
        {
            nuSet_ = true;
        }

        coeffs_.readCheckIfPresent
        (
            "imageDistance",
            imageDistance_,
            scalarMinMax::ge(0.5)
        );
        imageInterpolationNames_.readIfPresent
        (
            "imageInterpolation",
            coeffs_,
            imageInterpolation_
        );

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

        if (method_ == ghostCell)
        {
            Info<< "    Ghost cell method: penalty coefficient "
                << penaltyCoeff_ << ", image distance " << imageDistance_
                << ", " << imageInterpolationNames_[imageInterpolation_]
                << " image interpolation";
        }
        else if (method_ == penalty)
        {
            Info<< "    Penalty method: coefficient " << penaltyCoeff_
                << ", " << weightingNames_[weighting_]
                << " weighting, surface rate "
                << "coefficient " << surfaceRateCoeff_;
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
