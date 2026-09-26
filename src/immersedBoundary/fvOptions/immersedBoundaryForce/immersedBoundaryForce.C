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
#include "cutFaceIso.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcGrad.H"
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
        nuSet_ = true;
    }

    return nu_;
}


Foam::scalar Foam::fv::immersedBoundaryForce::cellWidth
(
    const label celli
) const
{
    const vector span
    (
        boundBox(mesh_.points(), mesh_.cellPoints()[celli], false).span()
    );

    const Vector<label>& solD = mesh_.solutionD();

    scalar w = GREAT;
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        if (solD[d] == 1)
        {
            w = min(w, span[d]);
        }
    }

    return w;
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


Foam::tmp<Foam::vectorField>
Foam::fv::immersedBoundaryForce::ddtU(const volVectorField& U) const
{
    return tmp<vectorField>::New(fvc::ddt(U)().primitiveField());
}


Foam::tmp<Foam::vectorField>
Foam::fv::immersedBoundaryForce::residual(const volVectorField& U)
{
    // The residual ddt(U) + div(phi, U) - laplacian(nu, U) + grad(p) of the
    // final fields, which is the forcing that the fluid receives
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>("phi");

    const dimensionedScalar nu("nu", dimViscosity, this->nu());

    return tmp<vectorField>::New
    (
        (
            fvc::ddt(U)
          + fvc::div(phi, U)
          - fvc::laplacian(nu, U)
          + fvc::grad(p)
        )().primitiveField()
    );
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

    tmp<vectorField> residualPtr;
    tmp<vectorField> ddtUPtr;
    if (method_ == "cutLink")
    {
        residualPtr = residual(U);
        ddtUPtr = ddtU(U);
    }

    forAll(bodies_, bodyi)
    {
        const immersedBody& body = bodies_[bodyi];

        if (method_ == "cutLink")
        {
            // Momentum exchange from the residual of the momentum equation
            // of the final velocity and pressure without the forcing, in the
            // penalised cells and the fluid cells next to them, and the rate
            // of change of the momentum of the penalised cells with the time
            // scheme of the momentum equation, both over the same cells, so
            // that the force does not jump as cells change role
            const vectorField& R = residualPtr();
            const vectorField& dUdt = ddtUPtr();
            const scalarField& V = mesh_.V();
            const vectorField& C = mesh_.C();
            const point CofR(body.CofR());

            vector F(Zero);
            vector T(Zero);
            for (const label celli : body.insideCells())
            {
                const vector fc((R[celli] - dUdt[celli])*V[celli]);
                F -= fc;
                T -= (C[celli] - CofR) ^ fc;
            }
            // A fluid cell next to several bodies is shared in proportion to
            // the rates of its links to each body
            const labelList& lCells = linkCells_[bodyi];
            const scalarField& lRates = linkRates_[bodyi];
            forAll(lCells, i)
            {
                const label celli = lCells[i];
                const scalar share =
                (
                    linkRateSum_[celli] > VSMALL
                  ? lRates[i]/linkRateSum_[celli]
                  : 1
                );
                const vector fc(share*R[celli]*V[celli]);
                F -= fc;
                T -= (C[celli] - CofR) ^ fc;
            }
            reduce(F, sumOp<vector>());
            reduce(T, sumOp<vector>());

            // Force from the surface traction
            if (tractions_.size() != bodies_.size())
            {
                tractions_.setSize(bodies_.size());
            }
            if (!tractions_.set(bodyi))
            {
                tractions_.set
                (
                    bodyi,
                    new immersedSurfaceTraction(body, mesh_)
                );
            }

            vector Fs(Zero);
            vector Ts(Zero);
            tractions_[bodyi].force
            (
                U,
                mesh_.lookupObject<volScalarField>("p"),
                nu(),
                rho,
                CofR,
                Fs,
                Ts
            );

            // The selected estimate is written as the force, and the other
            // after the inertia (which is included in both)
            if (forceEstimator_ == "surfaceTraction")
            {
                force_[bodyi] = Fs;
                torque_[bodyi] = Ts;
                secondaryForce_[bodyi] = rho*F;
            }
            else
            {
                force_[bodyi] = rho*F;
                torque_[bodyi] = rho*T;
                secondaryForce_[bodyi] = Fs;
            }
            momentum_[bodyi] = Zero;
        }
        else
        {
            force_[bodyi] = body.force(f_, rho);
            torque_[bodyi] = body.torque(f_, rho);

            momentum_[bodyi] = body.momentum(U, lambda_);
        }

        inertia_[bodyi] = Zero;
        if (momentum0Valid_ && method_ != "cutLink")
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
                << Fi.x() << token::SPACE << Fi.y() << token::SPACE << Fi.z();

            if (method_ == "cutLink")
            {
                const vector& Ft = secondaryForce_[bodyi];
                os  << token::TAB << Ft.x() << token::SPACE << Ft.y()
                    << token::SPACE << Ft.z();
            }

            os  << endl;
        }
    }
}


Foam::scalar Foam::fv::immersedBoundaryForce::pinnedRate
(
    const label celli,
    const vector& Ub
)
{
    const scalar rate = penaltyCoeff_/mesh_.time().deltaTValue();

    if (pinnedRateCoeff_ <= 0)
    {
        return rate;
    }

    // Rate independent of the time step, relative to the viscous and
    // convective rates of the cell
    const scalar w = cellWidth(celli);
    return min(rate, pinnedRateCoeff_*(nu()/sqr(w) + mag(Ub)/w));
}


void Foam::fv::immersedBoundaryForce::setPenalty(const volVectorField& U)
{
    const scalar deltaT = mesh_.time().deltaTValue();
    const bool volumeFraction = (weighting_ == "volumeFraction");
    const scalarField& lambdaI = lambda_.primitiveField();

    Ui_ = U;

    for (const immersedBody& body : bodies_)
    {
        body.setVelocity(Ui_);
    }

    const scalar nu = (surfaceRateCoeff_ > 0 ? this->nu() : 0);

    kappa_ = 0;

    if (method_ == "cutLink")
    {
        Ui_ = U;
        vectorField& UiI = Ui_.primitiveFieldRef();

        // Penalised cells: the body velocity at their centre
        for (const immersedBody& body : bodies_)
        {
            const labelList& cells = body.insideCells();
            const vectorField Ub
            (
                body.velocity(pointField(mesh_.C(), cells))
            );

            forAll(cells, i)
            {
                kappa_[cells[i]] = pinnedRate(cells[i], Ub[i]);
                UiI[cells[i]] = Ub[i];
            }
        }

        // Fluid cells next to them: the body velocity at the surface,
        // weighted by the rates of the bodies
        forAll(bodies_, bodyi)
        {
            for (const label celli : linkCells_[bodyi])
            {
                UiI[celli] = Zero;
            }
        }
        forAll(bodies_, bodyi)
        {
            const labelList& cells = linkCells_[bodyi];
            const scalarField& rates = linkRates_[bodyi];
            const vectorField Ub
            (
                bodies_[bodyi].velocity(linkWallPoints_[bodyi])
            );

            forAll(cells, i)
            {
                kappa_[cells[i]] += rates[i];
                UiI[cells[i]] += rates[i]*Ub[i];
            }
        }
        // Normalise the target of each fluid cell once, also if it is next
        // to several bodies, and keep the sum of its rates for the forces
        linkRateSum_.setSize(mesh_.nCells());
        forAll(bodies_, bodyi)
        {
            for (const label celli : linkCells_[bodyi])
            {
                linkRateSum_[celli] = kappa_[celli];
            }
        }
        forAll(bodies_, bodyi)
        {
            for (const label celli : linkCells_[bodyi])
            {
                if (kappa_[celli] <= 0)
                {
                    continue;
                }

                // Normalised once: the rate is then marked as negative
                if (kappa_[celli] > VSMALL)
                {
                    UiI[celli] /= kappa_[celli];
                }
                else
                {
                    UiI[celli] = U[celli];
                }
                kappa_[celli] = -kappa_[celli];
            }
        }
        forAll(bodies_, bodyi)
        {
            for (const label celli : linkCells_[bodyi])
            {
                kappa_[celli] = mag(kappa_[celli]);
            }
        }

        // The rate of a fluid cell tends to that of the penalised cells as
        // its centre reaches the surface, where it becomes penalised, so
        // that the coefficients are continuous as the body moves
        forAll(bodies_, bodyi)
        {
            const vectorField& UbL = UiI;
            for (const label celli : linkCells_[bodyi])
            {
                kappa_[celli] =
                    min(kappa_[celli], pinnedRate(celli, UbL[celli]));
            }
        }

        return;
    }

    for (const immersedBody& body : bodies_)
    {
        if (method_ == "ghostCell")
        {
            // Sharp penalisation of the cells whose centre is inside
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

                    // Rate independent of the time step in the partially
                    // covered cells
                    if (surfaceRateCoeff_ > 0 && lambda < 1 - surfaceThreshold_)
                    {
                        const scalar w = cellWidth(celli);

                        kappa_[celli] = min
                        (
                            kappa_[celli],
                            lambda/(1 - lambda)*surfaceRateCoeff_
                           *(nu/sqr(w) + mag(Ui_[celli])/w)
                        );
                    }
                }
                else
                {
                    kappa_[celli] = penaltyCoeff_*lambda/deltaT;
                }
            }
        }
    }
}


void Foam::fv::immersedBoundaryForce::findLinks()
{
    const label nBodies = bodies_.size();
    const scalar nu = this->nu();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& C = mesh_.C();
    const scalarField& V = mesh_.V();
    const surfaceScalarField& deltaCoeffs = mesh_.nonOrthDeltaCoeffs();
    const surfaceScalarField& magSf = mesh_.magSf();

    // Body whose surface contains the centre of each cell (the first), -1
    // if none
    volScalarField insideBody
    (
        IOobject
        (
            IOobject::scopedName(name_, "insideBody"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimless, -1)
    );
    forAll(bodies_, bodyi)
    {
        for (const label celli : bodies_[bodyi].insideCells())
        {
            if (insideBody[celli] < 0)
            {
                insideBody[celli] = bodyi;
            }
        }
    }
    insideBody.correctBoundaryConditions();

    auto bodyOf = [](const scalar b) { return label(std::lround(b)); };

    // Links from a fluid cell P to a penalised cell N: P, the centre of N and
    // the diffusion coefficient of the face divided by the volume of P
    List<DynamicList<label>> linkP(nBodies);
    List<DynamicList<point>> linkEnd(nBodies);
    List<DynamicList<scalar>> linkCoeff(nBodies);
    List<DynamicList<label>> linkN(nBodies);
    List<DynamicList<label>> linkPatch(nBodies);
    List<DynamicList<label>> linkPatchFace(nBodies);

    auto addLink = [&]
    (
        const label P,
        const label bodyi,
        const point& end,
        const scalar coeff,
        const label N,
        const label patchi,
        const label patchFacei
    )
    {
        linkP[bodyi].append(P);
        linkEnd[bodyi].append(end);
        linkCoeff[bodyi].append(coeff);
        linkN[bodyi].append(N);
        linkPatch[bodyi].append(patchi);
        linkPatchFace[bodyi].append(patchFacei);
    };

    forAll(nei, facei)
    {
        const label o = own[facei];
        const label n = nei[facei];
        const label bo = bodyOf(insideBody[o]);
        const label bn = bodyOf(insideBody[n]);
        const scalar coeff = nu*magSf[facei]*deltaCoeffs[facei];

        if (bo < 0 && bn >= 0)
        {
            addLink(o, bn, C[n], coeff/V[o], n, -1, -1);
        }
        else if (bn < 0 && bo >= 0)
        {
            addLink(n, bo, C[o], coeff/V[n], o, -1, -1);
        }
    }

    forAll(insideBody.boundaryField(), patchi)
    {
        const fvPatchScalarField& pf = insideBody.boundaryField()[patchi];

        if (pf.coupled())
        {
            const scalarField nbrBody(pf.patchNeighbourField());
            const labelUList& faceCells = pf.patch().faceCells();
            const vectorField delta(pf.patch().delta());
            const scalarField& pMagSf = magSf.boundaryField()[patchi];
            const scalarField& pDeltaCoeffs =
                deltaCoeffs.boundaryField()[patchi];

            forAll(faceCells, i)
            {
                const label P = faceCells[i];
                const label bn = bodyOf(nbrBody[i]);

                if (bodyOf(insideBody[P]) < 0 && bn >= 0)
                {
                    addLink
                    (
                        P,
                        bn,
                        C[P] + delta[i],
                        nu*pMagSf[i]*pDeltaCoeffs[i]/V[P],
                        -1,
                        patchi,
                        i
                    );
                }
            }
        }
    }

    linkCells_.setSize(nBodies);
    linkRates_.setSize(nBodies);
    linkWallPoints_.setSize(nBodies);
    linkFluidCells_.setSize(nBodies);
    linkPenalisedCells_.setSize(nBodies);
    linkPatches_.setSize(nBodies);
    linkPatchFaces_.setSize(nBodies);
    linkCoeffs_.setSize(nBodies);
    linkPoints_.setSize(nBodies);

    label nLinks = 0;

    forAll(bodies_, bodyi)
    {
        const labelList& P = linkP[bodyi];
        const pointField start(C, P);
        const pointField end(linkEnd[bodyi]);

        List<pointIndexHit> hits;
        bodies_[bodyi].findLine(start, end, hits);

        pointField linkPts(P.size());

        Map<label> cellToLink(2*P.size());
        DynamicList<label> cells(P.size());
        DynamicList<scalar> rates(P.size());
        DynamicList<point> wallPoints(P.size());
        DynamicList<scalar> weights(P.size());

        forAll(P, linki)
        {
            const scalar d = mag(end[linki] - start[linki]);

            // Distance to the surface along the link: half the link if the
            // segment does not cross the surface (round-off)
            scalar phi = 0.5*d;
            point wallPoint = 0.5*(start[linki] + end[linki]);
            if (hits[linki].hit())
            {
                wallPoint = hits[linki].hitPoint();
                phi = mag(wallPoint - start[linki]);
            }
            phi = min(max(phi, 1e-6*d), d);
            linkPts[linki] = wallPoint;

            const scalar rate = linkCoeff[bodyi][linki]*(d/phi - 1);

            // Weight of the surface point: the rate, and a small part of the
            // coefficient for the links whose rate is zero
            const scalar weight = rate + 1e-6*linkCoeff[bodyi][linki];

            const auto iter = cellToLink.cfind(P[linki]);
            if (iter.good())
            {
                rates[iter.val()] += rate;
                wallPoints[iter.val()] += weight*wallPoint;
                weights[iter.val()] += weight;
            }
            else
            {
                cellToLink.insert(P[linki], cells.size());
                cells.append(P[linki]);
                rates.append(rate);
                wallPoints.append(weight*wallPoint);
                weights.append(weight);
            }
        }

        forAll(cells, i)
        {
            wallPoints[i] /= weights[i];
        }

        linkCells_[bodyi].transfer(cells);
        linkRates_[bodyi].transfer(rates);
        linkWallPoints_[bodyi].transfer(wallPoints);

        linkFluidCells_[bodyi] = P;
        linkPenalisedCells_[bodyi] = linkN[bodyi];
        linkPatches_[bodyi] = linkPatch[bodyi];
        linkPatchFaces_[bodyi] = linkPatchFace[bodyi];
        linkCoeffs_[bodyi] = linkCoeff[bodyi];
        linkPoints_[bodyi] = linkPts;

        nLinks += linkCells_[bodyi].size();
    }

    Info<< "    Cut link cells: " << returnReduce(nLinks, sumOp<label>())
        << endl;
}


void Foam::fv::immersedBoundaryForce::updateApertures()
{
    if (!aperturePtr_)
    {
        aperturePtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "immersedBoundaryAperture",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                mesh_,
                dimensionedScalar(dimless, 1)
            )
        );
        wallFluxPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "immersedBoundaryWallFlux",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                mesh_,
                dimensionedScalar(dimVolume/dimTime, Zero)
            )
        );

        // The flux is oriented with the face area vectors
        wallFluxPtr_->setOriented();
    }

    surfaceScalarField& alpha = *aperturePtr_;
    surfaceScalarField& wallFlux = *wallFluxPtr_;

    // Signed distance to the bodies at the points, positive in the fluid
    scalarField f(mesh_.nPoints(), GREAT);
    for (const immersedBody& body : bodies_)
    {
        body.addPointDistances(f);
    }

    cutFaceIso cutFace(mesh_, f);

    // Aperture of a face: the fluid (positive) part of its area
    auto aperture = [&](const label facei, const scalar magS)
    {
        const label status = cutFace.calcSubFace(facei, 0);
        if (status == -1)
        {
            return scalar(1);
        }
        else if (status == 1)
        {
            return scalar(0);
        }
        return min(mag(cutFace.subFaceArea())/max(magS, VSMALL), scalar(1));
    };

    const surfaceScalarField& magSf = mesh_.magSf();
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    // Apertures of the internal and coupled boundary faces; the other
    // boundary faces keep their fluxes
    scalarField& alphaI = alpha.primitiveFieldRef();
    forAll(alphaI, facei)
    {
        alphaI[facei] = aperture(facei, magSf[facei]);
    }
    forAll(alpha.boundaryField(), patchi)
    {
        fvsPatchScalarField& palpha = alpha.boundaryFieldRef()[patchi];
        const fvPatch& patch = mesh_.boundary()[patchi];

        if (patch.coupled())
        {
            forAll(palpha, i)
            {
                palpha[i] = aperture(patch.start() + i, patch.magSf()[i]);
            }
        }
        else
        {
            palpha = 1;
        }
    }

    // Faces with a solid part, and their centres
    DynamicList<label> solidFaces;
    DynamicList<point> solidCentres;
    forAll(alphaI, facei)
    {
        if (alphaI[facei] < 1)
        {
            solidFaces.append(facei);
            solidCentres.append(Cf[facei]);
        }
    }
    const label nInternal = solidFaces.size();
    forAll(alpha.boundaryField(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        const fvsPatchScalarField& palpha = alpha.boundaryField()[patchi];
        forAll(palpha, i)
        {
            if (palpha[i] < 1)
            {
                solidFaces.append(patch.start() + i);
                solidCentres.append(patch.Cf()[i]);
            }
        }
    }

    // Body velocity at these faces: that of the body nearest to the face
    // centre
    const pointField centres(solidCentres);
    vectorField Ub(centres.size(), Zero);
    {
        scalarField bestDist(centres.size(), GREAT);
        forAll(bodies_, bodyi)
        {
            List<pointIndexHit> hits;
            vectorField normals;
            bodies_[bodyi].nearest
            (
                centres,
                scalarField(centres.size(), GREAT),
                hits,
                normals
            );
            const vectorField Ubody(bodies_[bodyi].velocity(centres));
            forAll(hits, i)
            {
                if (hits[i].hit())
                {
                    const scalar d = mag(hits[i].hitPoint() - centres[i]);
                    if (d < bestDist[i])
                    {
                        bestDist[i] = d;
                        Ub[i] = Ubody[i];
                    }
                }
            }
        }
    }

    // Flux of the body velocity through the solid part of the faces
    wallFlux = dimensionedScalar(wallFlux.dimensions(), Zero);
    scalarField& wallFluxI = wallFlux.primitiveFieldRef();
    forAll(solidFaces, i)
    {
        const label facei = solidFaces[i];
        if (i < nInternal)
        {
            wallFluxI[facei] = (1 - alphaI[facei])*(Ub[i] & Sf[facei]);
        }
        else
        {
            const label patchi = mesh_.boundaryMesh().whichPatch(facei);
            const fvPatch& patch = mesh_.boundary()[patchi];
            const label pfacei = facei - patch.start();
            if (patch.coupled())
            {
                wallFlux.boundaryFieldRef()[patchi][pfacei] =
                    (1 - alpha.boundaryField()[patchi][pfacei])
                   *(Ub[i] & patch.Sf()[pfacei]);
            }
        }
    }
}


void Foam::fv::immersedBoundaryForce::findGhostCells(const bool newTimeStep)
{
    // polyMesh::findCell uses the tet base points, whose construction is
    // collective: construct them on every processor before searching
    (void)mesh_.tetBasePtIs();
    (void)mesh_.cellTree();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const pointField& meshPoints = mesh_.points();
    const labelListList& cellPoints = mesh_.cellPoints();

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
    // mesh moves), these are kept. The mesh topology is assumed not to change
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

    const pointField ghostCentres(mesh_.C(), ghostCells_);

    vectorField spans(nGhosts);
    scalarField searchDistSqr(nGhosts);
    forAll(ghostCells_, gi)
    {
        spans[gi] =
            boundBox(meshPoints, cellPoints[ghostCells_[gi]], false).span();
        searchDistSqr[gi] = 4*magSqr(spans[gi]);
    }

    // Nearest surface point of the body that contains each ghost cell
    forAll(bodies_, bodyi)
    {
        List<pointIndexHit> hits;
        vectorField normals;
        bodies_[bodyi].nearest(ghostCentres, searchDistSqr, hits, normals);

        DynamicList<label> hitGhosts(nGhosts);
        DynamicField<point> hitPoints(nGhosts);

        forAll(hits, gi)
        {
            if (cellBody[ghostCells_[gi]] == bodyi && hits[gi].hit())
            {
                ghostWallPoints_[gi] = hits[gi].hitPoint();

                // Outward normal along the line from the cell centre, which
                // is inside the body, to the surface point, so that it varies
                // continuously as the body moves; the normal of the nearest
                // face if the centre is on the surface
                const vector r(ghostCentres[gi] - hits[gi].hitPoint());
                const scalar magR = mag(r);
                if (magR > 1e-3*mag(spans[gi]))
                {
                    ghostNormals_[gi] = -r/magR;
                }
                else
                {
                    ghostNormals_[gi] = normals[gi];
                }

                // Signed distance from the surface, along the outward
                // normal: negative inside the body
                ghostDistances_[gi] = r & ghostNormals_[gi];

                hitGhosts.append(gi);
                hitPoints.append(hits[gi].hitPoint());
            }
        }

        const vectorField Ub(bodies_[bodyi].velocity(hitPoints));
        forAll(hitGhosts, i)
        {
            ghostWallVelocities_[hitGhosts[i]] = Ub[i];
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
        widths[gi] = cmptSum(cmptMultiply(cmptMag(n), spans[gi]));

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
    const bool cellPoint = (imageInterpolation_ == "cellPoint");

    autoPtr<interpolationCellPoint<vector>> interpPtr;
    if (cellPoint)
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
    method_("penalty"),
    penaltyCoeff_(1e3),
    weighting_("volumeFraction"),
    surfaceRateCoeff_(3),
    nu_(0),
    nuSet_(false),
    kappa_(mesh.nCells(), Zero),
    linkCells_(),
    linkRates_(),
    linkWallPoints_(),
    linkFluidCells_(),
    linkPenalisedCells_(),
    linkPatches_(),
    linkPatchFaces_(),
    linkCoeffs_(),
    linkPoints_(),
    linkCorrection_(true),
    linkRateSum_(),
    tractions_(),
    pinnedRateCoeff_(100),
    apertureCoupling_(true),
    aperturePtr_(),
    wallFluxPtr_(),
    forceEstimator_("momentumExchange"),
    imageDistance_(1.5),
    imageInterpolation_("cellPoint"),
    insideCells_(),
    penalisedOld_(),
    ghostCells_(),
    ghostDonors_(),
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
    secondaryForce_(),
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
    secondaryForce_.setSize(nBodies, Zero);

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
                << "inertia of the fluid inside (x y z)";

            if (method_ == "cutLink")
            {
                forceFiles_[bodyi]
                    << token::TAB
                    << (
                           forceEstimator_ == "surfaceTraction"
                         ? "force from the momentum exchange (x y z)"
                         : "force from the surface traction (x y z)"
                       );
            }

            forceFiles_[bodyi] << endl;
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
    const bool newTimeStep = (timeIndex != timeIndex_);

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

    if
    (
        method_ == "penalty"
     || method_ == "ghostCell"
     || method_ == "cutLink"
    )
    {
        // Implicit volume penalisation: kappa*(Ui - U)
        const volVectorField& U = eqn.psi();

        // Collective: updated is the same on every processor
        if (method_ == "cutLink" && (updated || linkCells_.empty()))
        {
            findLinks();

            if (apertureCoupling_)
            {
                updateApertures();
            }
        }

        setPenalty(U);

        if (method_ == "ghostCell")
        {
            // Collective: updated is the same on every processor
            if (updated)
            {
                findGhostCells(newTimeStep);
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

        if (method_ == "cutLink" && linkCorrection_)
        {
            // Explicit correction of the flux from each fluid cell to its
            // penalised neighbours: nu|Sf|deltaCoeff*(Ub - U_N)/V_P, so that
            // the total flux is that of a linear profile through the body
            // velocity on the surface
            volVectorField::Internal correction
            (
                IOobject
                (
                    IOobject::scopedName(name_, "linkCorrection"),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_,
                dimensionedVector(dimVelocity/dimTime, Zero)
            );
            vectorField& cI = correction.field();

            // Velocity of the neighbour cells of the coupled patches
            PtrList<vectorField> nbrU(mesh_.boundary().size());
            forAll(U.boundaryField(), patchi)
            {
                if (U.boundaryField()[patchi].coupled())
                {
                    nbrU.set
                    (
                        patchi,
                        U.boundaryField()[patchi].patchNeighbourField()
                    );
                }
            }

            forAll(bodies_, bodyi)
            {
                const labelList& P = linkFluidCells_[bodyi];
                const labelList& N = linkPenalisedCells_[bodyi];
                const scalarField& coeff = linkCoeffs_[bodyi];
                const vectorField Ub(bodies_[bodyi].velocity(linkPoints_[bodyi]));

                forAll(P, linki)
                {
                    vector UN(Zero);
                    if (N[linki] >= 0)
                    {
                        UN = U[N[linki]];
                    }
                    else
                    {
                        const label patchi = linkPatches_[bodyi][linki];
                        UN = nbrU[patchi][linkPatchFaces_[bodyi][linki]];
                    }
                    cI[P[linki]] += coeff[linki]*(Ub[linki] - UN);
                }
            }

            eqn += correction;
        }

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

    if
    (
        method_ == "penalty"
     || method_ == "ghostCell"
     || method_ == "cutLink"
    )
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

        if
        (
            method_ != "incremental"
         && method_ != "penalty"
         && method_ != "ghostCell"
         && method_ != "cutLink"
        )
        {
            FatalIOErrorInFunction(coeffs_)
                << "Unknown method " << method_ << ": valid methods are "
                << "penalty, cutLink, ghostCell and incremental"
                << exit(FatalIOError);
        }

        coeffs_.readCheckIfPresent
        (
            "penaltyCoeff",
            penaltyCoeff_,
            scalarMinMax::ge(1)
        );
        coeffs_.readIfPresent("weighting", weighting_);
        coeffs_.readIfPresent("surfaceRateCoeff", surfaceRateCoeff_);
        if (coeffs_.readIfPresent("nu", nu_))
        {
            nuSet_ = true;
        }

        if (weighting_ != "volumeFraction" && weighting_ != "occupancy")
        {
            FatalIOErrorInFunction(coeffs_)
                << "Unknown weighting " << weighting_ << ": valid weightings "
                << "are volumeFraction and occupancy" << exit(FatalIOError);
        }

        coeffs_.readCheckIfPresent
        (
            "imageDistance",
            imageDistance_,
            scalarMinMax::ge(0.5)
        );
        coeffs_.readIfPresent("imageInterpolation", imageInterpolation_);
        coeffs_.readIfPresent("apertureCoupling", apertureCoupling_);
        coeffs_.readIfPresent("pinnedRateCoeff", pinnedRateCoeff_);
        coeffs_.readIfPresent("linkCorrection", linkCorrection_);
        coeffs_.readIfPresent("forceEstimator", forceEstimator_);
        if
        (
            forceEstimator_ != "momentumExchange"
         && forceEstimator_ != "surfaceTraction"
        )
        {
            FatalIOErrorInFunction(coeffs_)
                << "Unknown forceEstimator " << forceEstimator_
                << ": valid estimators are momentumExchange and "
                << "surfaceTraction" << exit(FatalIOError);
        }

        if (imageInterpolation_ != "cellPoint" && imageInterpolation_ != "cell")
        {
            FatalIOErrorInFunction(coeffs_)
                << "Unknown imageInterpolation " << imageInterpolation_
                << ": valid options are cellPoint and cell"
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

        if (method_ == "cutLink")
        {
            Info<< "    Cut link method: penalty coefficient "
                << penaltyCoeff_;
        }
        else if (method_ == "ghostCell")
        {
            Info<< "    Ghost cell method: penalty coefficient "
                << penaltyCoeff_ << ", image distance " << imageDistance_
                << ", " << imageInterpolation_ << " image interpolation";
        }
        else if (method_ == "penalty")
        {
            Info<< "    Penalty method: coefficient " << penaltyCoeff_
                << ", " << weighting_ << " weighting, surface rate "
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
