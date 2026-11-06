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

Class
    Foam::thinBodyPockets

Description
    Pressure regions that the closed faces of a thin immersed body (zero
    diffusivity, aperture 0) separate from every fixed-pressure boundary face
    (pockets), used by pimpleFluid with immersedBoundaryForce
    thinBodyBlocking.

    The pressure matrix of a pocket is singular (a constant is a solution),
    and its continuity equation is solvable only if the net volume flux into
    it, from the walls that move, is zero. The class
      - finds the pockets (regionSplit of the mesh across the closed faces,
        also across processors),
      - gives each a reference pressure cell,
      - treats a net flux that is not zero according to pocketTreatment of
        the registered immersedBoundaryThinBody dictionary: flag (the default:
        the net flux is left in the equation, absorbed at the reference
        cell, and reported), reject (a fatal error above pocketTolerance) or
        source (the net flux is removed uniformly from the pocket cells: an
        artificial volume source, NOT physical conservation, explicit only),
      - keeps the volume mean of the pressure of the cells of a pocket
        through the solution, so that its level does not jump when the
        pocket forms or merges (the level does not enter the momentum
        equation; it enters only pressure output and pressure integrals),
      - writes one line per time step to
        postProcessing/thinBodyPockets/pockets.dat (time, pockets, pockets
        without a fixed pressure, largest and summed absolute net volume
        flux that is not zero, largest relative imbalance, volume of the
        pocket with the largest imbalance) and warns above
        pocketTolerance, once per time step.

SourceFiles
    thinBodyPockets.H

\*---------------------------------------------------------------------------*/

#ifndef thinBodyPockets_H
#define thinBodyPockets_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "regionSplit.H"
#include "fvMatrix.H"
#include "fvcDiv.H"
#include "OFstream.H"
#include "IOdictionary.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"

namespace Foam
{

class thinBodyPockets
{
    // Private data

        //- Time index of the last log line
        label logIndex_;

        //- Pockets (region numbers without a fixed pressure) and their data
        labelList regionOf_;
        label nRegions_;
        DynamicList<label> pockets_;
        scalarField vol_;
        scalarField netDiv_;
        scalarField absDiv_;
        scalarField meanP_;
        labelList refCells_;
        boolList pocketFlag_;

        // Values of the current time step, for the report
        scalar pocketNet_, pocketVol_, removed_;
        scalar maxNetStep_, sumNetStep_, maxRelStep_, tolStep_;
        label nWarnStep_;
        word treatmentStep_;

public:

    thinBodyPockets()
    :
        logIndex_(-1),
        regionOf_(),
        nRegions_(0),
        pocketNet_(0), pocketVol_(0), removed_(0),
        maxNetStep_(0), sumNetStep_(0), maxRelStep_(0), tolStep_(0),
        nWarnStep_(0)
    {}

    //- Find the pockets. closed is the face flag (internal and coupled),
    //  needReference/pRefCell the closed-domain reference of the whole
    //  domain, whose region is not a pocket
    void analyse
    (
        const fvMesh& mesh,
        const volScalarField& p,
        const boolList& closed,
        const bool needReference,
        const label pRefCell
    )
    {
        const regionSplit regions(mesh, closed);
        nRegions_ = regions.nRegions();
        regionOf_ = regions;

        boolList hasFixed(nRegions_, false);
        forAll(p.boundaryField(), patchi)
        {
            if (p.boundaryField()[patchi].fixesValue())
            {
                for (const label celli : mesh.boundary()[patchi].faceCells())
                {
                    hasFixed[regionOf_[celli]] = true;
                }
            }
        }
        if (needReference && pRefCell >= 0)
        {
            hasFixed[regionOf_[pRefCell]] = true;
        }
        Pstream::listCombineGather(hasFixed, orEqOp<bool>());
        Pstream::listCombineScatter(hasFixed);

        vol_ = scalarField(nRegions_, 0.0);
        meanP_ = scalarField(nRegions_, 0.0);
        labelList firstCell(nRegions_, -1);
        const scalarField& V = mesh.V();
        forAll(regionOf_, celli)
        {
            const label ri = regionOf_[celli];
            vol_[ri] += V[celli];
            meanP_[ri] += V[celli]*p[celli];
            if (firstCell[ri] < 0)
            {
                firstCell[ri] = celli;
            }
        }
        reduce(vol_, sumOp<scalarField>());
        reduce(meanP_, sumOp<scalarField>());

        // The lowest processor with cells in a region sets its reference
        labelList owner(nRegions_, Pstream::nProcs());
        forAll(firstCell, ri)
        {
            if (firstCell[ri] >= 0)
            {
                owner[ri] = Pstream::myProcNo();
            }
        }
        Pstream::listCombineGather(owner, minEqOp<label>());
        Pstream::listCombineScatter(owner);

        pockets_.clear();
        DynamicList<label> refCells;
        forAll(hasFixed, ri)
        {
            if (!hasFixed[ri] && vol_[ri] > VSMALL)
            {
                pockets_.append(ri);
                if (owner[ri] == Pstream::myProcNo())
                {
                    refCells.append(firstCell[ri]);
                }
            }
        }
        refCells_ = refCells;
        pocketFlag_ = boolList(nRegions_, false);
        for (const label ri : pockets_)
        {
            pocketFlag_[ri] = true;
        }
        forAll(meanP_, ri)
        {
            if (vol_[ri] > VSMALL)
            {
                meanP_[ri] /= vol_[ri];
            }
        }
    }

    label nPockets() const
    {
        return pockets_.size();
    }

    //- Measure the net volume flux of each pocket and apply the policy (see
    //  the class description). div is the divergence of phiHbyA
    void treatSources
    (
        const fvMesh& mesh,
        scalarField& div
    )
    {
        const dictionary& dict =
            mesh.lookupObject<IOdictionary>("immersedBoundaryThinBody");
        const word treatment(dict.get<word>("pocketPolicy"));
        const scalar tol(dict.get<scalar>("pocketTolerance"));

        const scalarField& V = mesh.V();
        netDiv_ = scalarField(nRegions_, 0.0);
        absDiv_ = scalarField(nRegions_, 0.0);
        forAll(div, celli)
        {
            netDiv_[regionOf_[celli]] += div[celli]*V[celli];
            absDiv_[regionOf_[celli]] += mag(div[celli])*V[celli];
        }
        reduce(netDiv_, sumOp<scalarField>());
        reduce(absDiv_, sumOp<scalarField>());

        scalar maxNet = 0, sumNet = 0, maxRel = 0;
        label nWarn = 0;
        for (const label ri : pockets_)
        {
            const scalar rel =
                (absDiv_[ri] > 0) ? mag(netDiv_[ri])/absDiv_[ri] : 0;
            // Round-off of a compatible pocket is far below the tolerance
            const scalar net = (rel > tol) ? mag(netDiv_[ri]) : 0;
            maxNet = max(maxNet, net);
            sumNet += net;
            maxRel = max(maxRel, rel);
            if (rel > tol)
            {
                nWarn++;
            }
        }

        // Kept for the report written after the solution
        pocketNet_ = 0;
        pocketVol_ = 0;
        removed_ = 0;
        maxNetStep_ = maxNet;
        sumNetStep_ = sumNet;
        maxRelStep_ = maxRel;
        nWarnStep_ = nWarn;
        treatmentStep_ = treatment;
        tolStep_ = tol;
        scalar best = -1;
        for (const label ri : pockets_)
        {
            if (mag(netDiv_[ri]) > best)
            {
                best = mag(netDiv_[ri]);
                pocketNet_ = netDiv_[ri];
                pocketVol_ = vol_[ri];
            }
        }
        if (treatment == "source")
        {
            removed_ = sumNet;
        }

        if (treatment == "reject" && nWarn)
        {
            FatalErrorInFunction
                << "A pressure pocket has a net volume flux above the "
                << "relative tolerance " << tol << " (pocketPolicy reject)"
                << exit(FatalError);
        }

        if (treatment == "source")
        {
            forAll(div, celli)
            {
                const label ri = regionOf_[celli];
                if (vol_[ri] > VSMALL && isPocket(ri))
                {
                    div[celli] -= netDiv_[ri]/vol_[ri];
                }
            }
        }
    }

    //- Reference pressure of the pockets, as setReference of a closed domain
    void addReferences(fvScalarMatrix& pEqn, const scalar refValue) const
    {
        forAll(refCells_, i)
        {
            const label celli = refCells_[i];
            scalar& d = pEqn.diag()[celli];
            if (d < VSMALL)
            {
                d = 1;
                pEqn.source()[celli] += refValue;
            }
            else
            {
                pEqn.source()[celli] += d*refValue;
                d += d;
            }
        }
    }

    //- Report, after the solution of the pressure equation of a corrector:
    //  the net volume flux that the closed faces and the moving walls
    //  prescribe on the pocket with the largest one (outflow positive), the
    //  volume of the pocket cells, the source removed from the continuity
    //  equation (policy source only), and the net flux of the corrected flux
    //  phi on the pocket cells (the residual of continuity after the
    //  solution). One line per time step in
    //  postProcessing/thinBodyPockets/pockets.dat
    void report(const fvMesh& mesh, const surfaceScalarField& phi)
    {
        const scalarField divPhi(fvc::div(phi)().primitiveField());
        const scalarField& V = mesh.V();
        scalarField resid(nRegions_, 0.0);
        forAll(divPhi, celli)
        {
            resid[regionOf_[celli]] += divPhi[celli]*V[celli];
        }
        reduce(resid, sumOp<scalarField>());
        scalar maxResid = 0;
        for (const label ri : pockets_)
        {
            maxResid = max(maxResid, mag(resid[ri]));
        }

        const label ti = mesh.time().timeIndex();
        if (ti == logIndex_)
        {
            return;
        }
        logIndex_ = ti;

        if (Pstream::master())
        {
            const fileName dir
            (
                mesh.time().globalPath()/"postProcessing"/"thinBodyPockets"
            );
            mkDir(dir);
            const fileName file(dir/"pockets.dat");
            const bool fresh = !isFile(file);
            OFstream os
            (
                file,
                IOstreamOption(IOstreamOption::ASCII),
                IOstreamOption::APPEND
            );
            if (fresh)
            {
                os  << "# time regions pockets prescribedNetFlux(largest, "
                    << "outflow +) pocketVolume removedSource "
                    << "residualAfterSolve maxRelativeImbalance" << nl;
            }
            os  << mesh.time().value() << ' ' << nRegions_ << ' '
                << pockets_.size() << ' ' << pocketNet_ << ' ' << pocketVol_
                << ' ' << removed_ << ' ' << maxResid << ' ' << maxRelStep_
                << nl;
        }
        if (nWarnStep_)
        {
            WarningInFunction
                << nWarnStep_ << " pressure pocket(s) enclosed by the closed "
                << "faces of a thin body have a net volume flux above the "
                << "relative tolerance " << tolStep_ << " (largest "
                << maxNetStep_ << " m3/s, relative " << maxRelStep_
                << "): a strictly incompressible pocket whose enclosing "
                << "walls change its volume has no solution. pocketPolicy "
                << treatmentStep_ << " (source removes the flux from "
                << "continuity: an artificial volume source, not "
                << "conservation; flag leaves the imbalance at the "
                << "reference cell)" << endl;
        }
    }

    //- Shift the pressure of each pocket so that its volume mean is that
    //  before the solution (no effect on the fluxes or the momentum)
    void keepMeanPressure(const fvMesh& mesh, volScalarField& p) const
    {
        scalarField after(nRegions_, 0.0);
        const scalarField& V = mesh.V();
        forAll(regionOf_, celli)
        {
            after[regionOf_[celli]] += V[celli]*p[celli];
        }
        reduce(after, sumOp<scalarField>());
        scalarField& pI = p.primitiveFieldRef();
        forAll(regionOf_, celli)
        {
            const label ri = regionOf_[celli];
            if (vol_[ri] > VSMALL && isPocket(ri))
            {
                pI[celli] += meanP_[ri] - after[ri]/vol_[ri];
            }
        }
        p.correctBoundaryConditions();
    }

private:

    bool isPocket(const label ri) const
    {
        return pocketFlag_[ri];
    }
};

} // End namespace Foam

#endif
