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

Description
    Estimation of the Robin coefficient rho_s*hs of the elasticWallPressure
    boundary condition: a-priori models (pWaveSpeed, constant,
    thicknessLimited) and the iteration-based secant model, plus the
    per-iteration coupling diagnostics.

\*---------------------------------------------------------------------------*/

#include "elasticWallPressureFvPatchScalarField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidSolidInterface.H"
#include "compatibilityFunctions.H"
#include "backwardDdtScheme.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "processorPolyPatch.H"
#include "labelledTri.H"
#ifdef OPENFOAM_COM
    #include "triSurface.H"
    #include "triSurfaceSearch.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * * //

#ifndef OPENFOAM_COM
// Moller-Trumbore ray-triangle intersection: returns the ray parameter t of
// the hit in [0, 1] along start -> end, or -1 if there is no hit
static scalar rayTriangleHit
(
    const point& start,
    const vector& dir,
    const point& a,
    const point& b,
    const point& c
)
{
    const vector e1(b - a);
    const vector e2(c - a);
    const vector pv(dir ^ e2);
    const scalar det = e1 & pv;

    if (mag(det) < VSMALL)
    {
        return -1;
    }

    const scalar invDet = 1.0/det;
    const vector tv(start - a);
    const scalar u = (tv & pv)*invDet;

    if (u < 0 || u > 1)
    {
        return -1;
    }

    const vector qv(tv ^ e1);
    const scalar v = (dir & qv)*invDet;

    if (v < 0 || u + v > 1)
    {
        return -1;
    }

    const scalar t = (e2 & qv)*invDet;

    return (t >= 0 && t <= 1) ? t : -1;
}
#endif


// Scale factor lambda of the seed coefficient which minimises the largest
// Robin-Neumann contraction factor over the recorded impedance samples
// (solid impedance Ss, fluid impedance Sf, mean seed coefficient m):
//     rate_k(lambda) = (Sf/Ss) |Ss - lambda*m| / (Sf + lambda*m)
// Returns -1 if there are no usable samples
static scalar minimaxScale
(
    const UList<vector>& samples,
    const scalar lambdaMin,
    const scalar lambdaMax,
    scalar& bestRate
)
{
    bestRate = GREAT;
    scalar bestLambda = -1;

    const label nScan = 400;
    const scalar logMin = Foam::log(lambdaMin);
    const scalar logMax = Foam::log(lambdaMax);

    for (label i = 0; i <= nScan; i++)
    {
        const scalar lambda =
            Foam::exp(logMin + (logMax - logMin)*scalar(i)/nScan);

        scalar maxRate = 0;
        forAll(samples, sI)
        {
            const scalar Ss = samples[sI].x();
            const scalar Sf = samples[sI].y();
            const scalar alpha = lambda*samples[sI].z();

            maxRate = max(maxRate, (Sf/Ss)*mag(Ss - alpha)/(Sf + alpha));
        }

        // Strict inequality: ties favour the smaller, safer, coefficient
        if (maxRate < bestRate)
        {
            bestRate = maxRate;
            bestLambda = lambda;
        }
    }

    return bestLambda;
}


// * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * * //

const fluidSolidInterface& elasticWallPressureFvPatchScalarField::fsi() const
{
#ifdef OPENFOAM_NOT_EXTEND
    const fvMesh& mesh = internalField().mesh();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
#endif

    return
        mesh.objectRegistry::parent().lookupObject<fluidSolidInterface>
        (
            "fsiProperties"
        );
}


label elasticWallPressureFvPatchScalarField::interfaceIndex() const
{
    const labelList& fluidPatchIndices = fsi().fluidPatchIndices();

    forAll(fluidPatchIndices, i)
    {
        if (fluidPatchIndices[i] == patch().index())
        {
            return i;
        }
    }

    FatalErrorInFunction
        << "Are you sure this patch is an FSI interface?"
        << abort(FatalError);

    return -1;
}


tmp<scalarField> elasticWallPressureFvPatchScalarField::solidPatchToFluidPatch
(
    const scalarField& solidField
) const
{
    const fluidSolidInterface& fsi = this->fsi();
    const label interfaceID = interfaceIndex();

    // Take references to zones
    const standAlonePatch& fluidZone =
        fsi.fluid().globalPatches()[interfaceID].globalPatch();
    const standAlonePatch& solidZone =
        fsi.solid().globalPatches()[interfaceID].globalPatch();

    // Map the solid patch field to the global zone
    const scalarField solidZoneField
    (
        fsi.solid().globalPatches()[interfaceID].patchFaceToGlobal(solidField)
    );
    scalarField fluidZoneField(fluidZone.size(), 0.0);

    // Transfer the field from the solid interface to the fluid interface
    fsi.interfaceToInterfaceList()[interfaceID].transferFacesZoneToZone
    (
        solidZone,          // from zone
        fluidZone,          // to zone
        solidZoneField,     // from field
        fluidZoneField      // to field
    );

    // Map the zone field to the fluid patch
    return tmp<scalarField>
    (
        new scalarField
        (
            fsi.fluid().globalPatches()
            [
                interfaceID
            ].globalFaceToPatch(fluidZoneField)
        )
    );
}


scalar elasticWallPressureFvPatchScalarField::effectiveDeltaT() const
{
#ifdef OPENFOAM_NOT_EXTEND
    const fvMesh& mesh = internalField().mesh();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
#endif

    scalar deltaT = mesh.time().deltaT().value();

    const word ddtScheme =
#ifdef OPENFOAM_NOT_EXTEND
        word(mesh.ddtScheme("ddt(U)"));
#else
        mesh.schemesDict().ddtScheme("ddt(U)");
#endif

    if
    (
        ddtScheme == fv::backwardDdtScheme<vector>::typeName
     && mesh.time().timeIndex() > 1
    )
    {
        const scalar deltaT0 = mesh.time().deltaT0().value();
        const scalar Cn = 1 + deltaT/(deltaT + deltaT0);

        deltaT /= Cn;
    }

    return deltaT;
}


const scalarField& elasticWallPressureFvPatchScalarField::solidThickness() const
{
    if (solidThicknessPtr_.valid())
    {
        return solidThicknessPtr_();
    }

    const fluidSolidInterface& fsi = this->fsi();
    const fvMesh& sMesh = fsi.solidMesh();
    const polyBoundaryMesh& bm = sMesh.boundaryMesh();
    const pointField& meshPoints = sMesh.points();
    const label solidPatchID = fsi.solidPatchIndices()[interfaceIndex()];

    // Triangulate all physical solid boundary faces; the triangle region is
    // the solid patch index and each triangle also records the index of the
    // boundary face it belongs to
    DynamicList<point> localPoints;
    DynamicList<labelledTri> localTris;
    DynamicList<label> localTriFaces;
    labelList originFace(bm[solidPatchID].size(), -1);
    label nLocalFaces = 0;

    forAll(bm, patchI)
    {
        const polyPatch& pp = bm[patchI];

        if
        (
            pp.coupled()
         || isA<emptyPolyPatch>(pp)
         || isA<wedgePolyPatch>(pp)
        )
        {
            continue;
        }

        forAll(pp, faceI)
        {
            const face& f = pp[faceI];
            const label start = localPoints.size();

            forAll(f, fpI)
            {
                localPoints.append(meshPoints[f[fpI]]);
            }

            for (label fpI = 1; fpI < f.size() - 1; fpI++)
            {
                localTris.append
                (
                    labelledTri(start, start + fpI, start + fpI + 1, patchI)
                );
                localTriFaces.append(nLocalFaces);
            }

            if (patchI == solidPatchID)
            {
                originFace[faceI] = nLocalFaces;
            }

            nLocalFaces++;
        }
    }

    // Gather the boundary triangles from all processors so that rays can
    // cross processor boundaries
    List<pointField> procPoints(Pstream::nProcs());
    List<List<labelledTri>> procTris(Pstream::nProcs());
    List<labelList> procTriFaces(Pstream::nProcs());
    labelList procNFaces(Pstream::nProcs(), 0);
    procPoints[Pstream::myProcNo()].transfer(localPoints);
    procTris[Pstream::myProcNo()].transfer(localTris);
    procTriFaces[Pstream::myProcNo()].transfer(localTriFaces);
    procNFaces[Pstream::myProcNo()] = nLocalFaces;
    Pstream::gatherList(procPoints);
    Pstream::gatherList(procTris);
    Pstream::gatherList(procTriFaces);
    Pstream::gatherList(procNFaces);
    Pstream::scatterList(procPoints);
    Pstream::scatterList(procTris);
    Pstream::scatterList(procTriFaces);
    Pstream::scatterList(procNFaces);

    label nPoints = 0;
    label nTris = 0;
    forAll(procPoints, procI)
    {
        nPoints += procPoints[procI].size();
        nTris += procTris[procI].size();
    }

    pointField allPoints(nPoints);
    List<labelledTri> allTris(nTris);
    labelList allTriFaces(nTris);
    nPoints = 0;
    nTris = 0;
    label faceOffset = 0;
    label myFaceOffset = 0;
    forAll(procPoints, procI)
    {
        const pointField& pts = procPoints[procI];
        const List<labelledTri>& tris = procTris[procI];
        const labelList& triFaces = procTriFaces[procI];

        forAll(pts, pI)
        {
            allPoints[nPoints + pI] = pts[pI];
        }

        forAll(tris, tI)
        {
            const labelledTri& t = tris[tI];
            allTriFaces[nTris] = triFaces[tI] + faceOffset;
            allTris[nTris++] =
                labelledTri
                (
                    t[0] + nPoints, t[1] + nPoints, t[2] + nPoints, t.region()
                );
        }

        if (procI == Pstream::myProcNo())
        {
            myFaceOffset = faceOffset;
        }

        nPoints += pts.size();
        faceOffset += procNFaces[procI];
    }

    forAll(originFace, faceI)
    {
        originFace[faceI] += myFaceOffset;
    }

    // Maximum ray length: the solid bounding-box diagonal
    const boundBox bb(allPoints, false);
    const scalar maxLength = 2.0*mag(bb.max() - bb.min()) + SMALL;

    // Cast a ray from each interface face into the solid, along the inward
    // normal. Hits on the origin face itself, which occur for warped faces,
    // are skipped.
    const polyPatch& solidPatch = bm[solidPatchID];
    const vectorField Sf(solidPatch.faceAreas());
    const vectorField n(Sf/mag(Sf));
    const vectorField Cf(solidPatch.faceCentres());
    const pointField start(Cf + 1e-6*maxLength*n);
    const pointField end(Cf - maxLength*n);

    // Region (patch) and distance of the first hit
    labelList hitPatch(solidPatch.size(), -1);
    scalarField H(solidPatch.size(), maxLength);

#ifdef OPENFOAM_COM
    const triSurface surf(allTris, allPoints);
    const triSurfaceSearch search(surf);

    List<List<pointIndexHit>> hits;
    search.findLineAll(start, end, hits);

    forAll(hits, faceI)
    {
        const List<pointIndexHit>& faceHits = hits[faceI];
        scalar dMin = GREAT;

        forAll(faceHits, hI)
        {
            const label triI = faceHits[hI].index();

            if (allTriFaces[triI] == originFace[faceI])
            {
                continue;
            }

            const scalar d = mag(faceHits[hI].hitPoint() - Cf[faceI]);

            if (d < dMin)
            {
                dMin = d;
                H[faceI] = d;
                hitPatch[faceI] = surf[triI].region();
            }
        }
    }
#else
    // Brute-force search, only used by forks without triSurfaceSearch
    forAll(start, faceI)
    {
        const vector dir(end[faceI] - start[faceI]);
        scalar tMin = GREAT;

        forAll(allTris, tI)
        {
            if (allTriFaces[tI] == originFace[faceI])
            {
                continue;
            }

            const labelledTri& t = allTris[tI];
            const scalar tHit =
                rayTriangleHit
                (
                    start[faceI],
                    dir,
                    allPoints[t[0]],
                    allPoints[t[1]],
                    allPoints[t[2]]
                );

            if (tHit >= 0 && tHit < tMin)
            {
                tMin = tHit;
                hitPatch[faceI] = t.region();
            }
        }

        if (hitPatch[faceI] != -1)
        {
            H[faceI] = mag(start[faceI] + tMin*dir - Cf[faceI]);
        }
    }
#endif

    if (debug)
    {
        // Hit statistics per patch, over the non-processor patches, which
        // are the same on all processors
        label nPatches = bm.size();
        forAll(bm, patchI)
        {
            if (isA<processorPolyPatch>(bm[patchI]))
            {
                nPatches = patchI;
                break;
            }
        }
        scalarField nHits(nPatches, 0.0);
        scalarField sumH(nPatches, 0.0);
        forAll(hitPatch, faceI)
        {
            if (hitPatch[faceI] != -1 && hitPatch[faceI] < nPatches)
            {
                nHits[hitPatch[faceI]]++;
                sumH[hitPatch[faceI]] += H[faceI];
            }
        }
        forAll(nHits, patchI)
        {
            reduce(nHits[patchI], sumOp<scalar>());
            reduce(sumH[patchI], sumOp<scalar>());
        }
        forAll(nHits, patchI)
        {
            if (nHits[patchI] > 0)
            {
                Info<< "    rays exiting through " << bm[patchI].name()
                    << ": " << nHits[patchI] << ", mean distance = "
                    << sumH[patchI]/nHits[patchI] << endl;
            }
        }
    }

    // Walls wetted on both sides: the bending (translation) mode is loaded
    // from both faces, so each face carries half of the wall inertia
    const labelList& solidPatchIndices = fsi.solidPatchIndices();
    label nTwoSided = 0;
    label nMissed = 0;
    forAll(H, faceI)
    {
        if (hitPatch[faceI] == -1)
        {
            nMissed++;
        }
        else if
        (
            hsControls_.twoSidedHalving
         && findIndex(solidPatchIndices, hitPatch[faceI]) != -1
        )
        {
            H[faceI] *= 0.5;
            nTwoSided++;
        }
    }

    Info<< type() << " " << patch().name() << ": solid wall thickness: min = "
        << gMin(H) << ", max = " << gMax(H) << ", mean = " << gAverage(H)
        << "; two-sided faces = " << returnReduce(nTwoSided, sumOp<label>())
        << ", faces without a hit = " << returnReduce(nMissed, sumOp<label>())
        << " of " << returnReduce(H.size(), sumOp<label>()) << endl;

    solidThicknessPtr_.reset(new scalarField(H));

    return solidThicknessPtr_();
}


tmp<scalarField> elasticWallPressureFvPatchScalarField::calcSeedCoeff
(
    const word& model
) const
{
    const fluidSolidInterface& fsi = this->fsi();
    const label solidPatchID = fsi.solidPatchIndices()[interfaceIndex()];

    // Solid density on the solid patch; taken from the solid model because
    // the rho field is only registered once the solid has been solved,
    // e.g. not before couplingStartTime
    const scalarField& rho = fsi.solid().rho().boundaryField()[solidPatchID];

    const scalar deltaT = effectiveDeltaT();

    // Distance travelled by the p-wave in one effective time step
    scalarField ell(rho.size(), GREAT);
    if (model == "pWaveSpeed" || model == "thicknessLimited")
    {
        if (hsControls_.waveSpeed > 0)
        {
            ell = hsControls_.waveSpeed*deltaT;
        }
        else if (fsi.solidMesh().foundObject<volScalarField>("impK"))
        {
            // Solid stiffness (impK for generality)
            const scalarField& impK =
                fsi.solidMesh().lookupObject<volScalarField>
                (
                    "impK"
                ).boundaryField()[solidPatchID];

            ell = sqrt(impK/rho)*deltaT;
        }
        else if (model == "pWaveSpeed")
        {
            FatalErrorInFunction
                << "The solid model does not provide impK: set waveSpeed or "
                << "use hsModel thicknessLimited or constant on patch "
                << patch().name() << abort(FatalError);
        }
        else if (seedCoeffTimeIndex_ == -1)
        {
            WarningInFunction
                << "The solid model does not provide impK: hs on patch "
                << patch().name() << " is limited by the wall thickness only"
                << endl;
        }
    }

    // Solid "virtual thickness"
    scalarField hs(rho.size(), 0.0);
    if (model == "constant")
    {
        hs = constantHs_;
    }
    else if (model == "pWaveSpeed")
    {
        // Calculate a virtual thickness based on the p-wave speed and the
        // effective time step used by the fluid momentum equation
        hs = ell;
    }
    else if (model == "thicknessLimited")
    {
        // Implicit elastodynamic slab with a free back face: the interface
        // impedance is rho*l*tanh(H/l), which tends to the half-space value
        // rho*l for thick walls and the wall inertia rho*H for thin walls
        const scalarField& H = solidThickness();

        if (hsControls_.thicknessBlend == "tanh")
        {
            hs = ell*tanh(H/ell);
        }
        else
        {
            hs = min(ell, H);
        }
    }

    if (debug || seedCoeffTimeIndex_ == -1)
    {
        Info<< type() << " " << patch().name() << ": " << model
            << " hs: min = " << gMin(hs) << ", max = " << gMax(hs)
            << ", mean = " << gAverage(hs) << endl;
    }

    if (model == "constant" && hsUserPtr_.valid())
    {
        // The user hs field is defined on the fluid patch
        return solidPatchToFluidPatch(rho)*hsUserPtr_();
    }

    return solidPatchToFluidPatch(rho*hs);
}


const scalarField& elasticWallPressureFvPatchScalarField::seedCoeff() const
{
    const label timeIndex = db().time().timeIndex();

    const word model =
        hsControls_.model == "secant"
      ? hsControls_.seedModel
      : hsControls_.model;

    if
    (
        seedCoeffPtr_.empty()
     || (model != "constant" && seedCoeffTimeIndex_ != timeIndex)
    )
    {
        seedCoeffPtr_.reset(calcSeedCoeff(model).ptr());
        seedCoeffTimeIndex_ = timeIndex;
    }

    return seedCoeffPtr_();
}


tmp<scalarField>
elasticWallPressureFvPatchScalarField::fluidNormalAcceleration() const
{
    // Normal fluid acceleration implied by the pressure gradient,
    // a_n = -(1/rho_f) dp/dn, consistent with the Robin condition
    tmp<scalarField> tacc(-snGrad());

#ifdef OPENFOAM_NOT_EXTEND
    const fvMesh& mesh = internalField().mesh();
    const word fieldName = internalField().name();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
    const word fieldName = dimensionedInternalField().name();
#endif

    const dimensionSet& pDims =
        mesh.lookupObject<volScalarField>(fieldName).dimensions();

    if (pDims != dimPressure/dimDensity)
    {
        if (mesh.foundObject<volScalarField>("rho"))
        {
            // interFluid (p_rgh): include the gravity and surface-tension
            // source of the wall-normal momentum balance,
            // rho*a_n = gSrc - snGrad(p_rgh), where phig = gSrc*rAUf*|Sf|,
            // if these fields are available (they are local to the pressure
            // equation; at walls with a zero-gradient phase fraction the
            // source vanishes, so this only affects the impedance
            // diagnostics near a contact line)
            if
            (
                mesh.foundObject<surfaceScalarField>("phig")
             && mesh.foundObject<surfaceScalarField>("rAUf")
            )
            {
                const scalarField& phig =
                    patch().lookupPatchField<surfaceScalarField, scalar>
                    (
                        "phig"
                    );
                const scalarField& rAUf =
                    patch().lookupPatchField<surfaceScalarField, scalar>
                    (
                        "rAUf"
                    );
                tmpRef(tacc) += phig/(rAUf*patch().magSf());
            }

            tmpRef(tacc) /=
                patch().lookupPatchField<volScalarField, scalar>("rho");
        }
        else
        {
            const dictionary& transportProperties =
                db().lookupObject<IOdictionary>("transportProperties");

            const dimensionedScalar rhoFluid
            (
                transportProperties.lookup("rho")
            );

            tmpRef(tacc) /= rhoFluid.value();
        }
    }

    return tacc;
}


void elasticWallPressureFvPatchScalarField::smoothPatchField
(
    scalarField& fld,
    const label nSweeps
) const
{
    const labelListList& faceFaces = patch().patch().faceFaces();

    for (label sweepI = 0; sweepI < nSweeps; sweepI++)
    {
        const scalarField fld0(fld);

        forAll(fld, faceI)
        {
            const labelList& nbrs = faceFaces[faceI];

            if (nbrs.size())
            {
                scalar sumNbr = 0;
                forAll(nbrs, nI)
                {
                    sumNbr += fld0[nbrs[nI]];
                }

                fld[faceI] = 0.5*fld0[faceI] + 0.5*sumNbr/nbrs.size();
            }
        }
    }
}


void elasticWallPressureFvPatchScalarField::setSecantScale
(
    const scalarField& newScale
)
{
    forAll(secantScale_, faceI)
    {
        scalar s = newScale[faceI];

        // Limit the change per update
        s = min
        (
            max(s, secantScale_[faceI]/hsControls_.secantMaxChange),
            secantScale_[faceI]*hsControls_.secantMaxIncrease
        );

        // Bound relative to the seed
        secantScale_[faceI] =
            min
            (
                max(s, hsControls_.secantMinFactor),
                hsControls_.secantMaxFactor
            );
    }

    Info<< type() << " " << patch().name() << ": secant scale: min = "
        << gMin(secantScale_) << ", max = " << gMax(secantScale_)
        << ", mean = " << gAverage(secantScale_) << endl;
}


void elasticWallPressureFvPatchScalarField::applySecantUpdate()
{
    scalar lambda = -1;

    if (hsControls_.secantFit == "minimax")
    {
        // Samples of the current and recent time steps: modes that are
        // quiet for the current coefficient may become dominant after it
        // changes, so all recently observed modes are considered
        DynamicList<vector> samples(secantSamples_);
        samples.append(secantSamplesHist_);

        if (samples.size())
        {
            scalar rate = 0;
            lambda =
                minimaxScale
                (
                    samples,
                    hsControls_.secantMinFactor,
                    hsControls_.secantMaxFactor,
                    rate
                );

            if (debug || hsControls_.writeDiagnostics)
            {
                Info<< type() << " " << patch().name() << ": minimax scale = "
                    << lambda << " from " << samples.size()
                    << " samples, predicted worst rate = " << rate << endl;
            }
        }
    }

    if (lambda < 0)
    {
        if (secantLambdas_.empty())
        {
            return;
        }

        // Median of the per-iteration least-squares fits of
        // dp = lambda*alpha0*da, so that no single iteration pair dominates
        scalarList lambdas(secantLambdas_);
        sort(lambdas);
        const label n = lambdas.size();
        lambda =
            n % 2
          ? lambdas[n/2]
          : 0.5*(lambdas[n/2 - 1] + lambdas[n/2]);
    }

    scalarField newScale(size(), lambda);

    if (hsControls_.secantSpatial && hsControls_.secantFit == "median")
    {
        scalarField num(secantNumField_);
        scalarField den(secantDenField_);
        smoothPatchField(num, hsControls_.secantSmoothing);
        smoothPatchField(den, hsControls_.secantSmoothing);

        const scalar denTol = 1e-3*gMax(den);

        forAll(newScale, faceI)
        {
            if (den[faceI] > denTol && num[faceI] > 0)
            {
                newScale[faceI] = num[faceI]/den[faceI];
            }
        }
    }

    setSecantScale(newScale);
}


void elasticWallPressureFvPatchScalarField::leakageFluxes
(
    scalar& leakFlux,
    scalar& netLeakFlux,
    scalar& wallMotionFlux,
    scalar& carryOverFlux
) const
{
    leakFlux = 0;
    netLeakFlux = 0;
    wallMotionFlux = 0;
    carryOverFlux = 0;

#ifdef OPENFOAM_NOT_EXTEND
    const fvMesh& mesh = internalField().mesh();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
#endif

    const label patchI = patch().index();

    // Flux relative to the mesh motion through the interface
    if (mesh.foundObject<surfaceScalarField>("phi"))
    {
        const scalarField& phi =
            mesh.lookupObject<surfaceScalarField>("phi").boundaryField()
            [
                patchI
            ];

        leakFlux = gSum(mag(phi));
        netLeakFlux = gSum(phi);
    }

    if (!mesh.moving())
    {
        return;
    }

    const scalarField& meshPhi = mesh.phi().boundaryField()[patchI];
    wallMotionFlux = gSum(mag(meshPhi));

    // Mismatch between the fluid interface velocity at the start of the time
    // step, from which the Robin interface flux is built, and the old mesh
    // velocity
    if (mesh.foundObject<surfaceVectorField>("Uf"))
    {
        const surfaceVectorField& Uf =
            mesh.lookupObject<surfaceVectorField>("Uf");

        carryOverFlux =
            gSum
            (
                mag
                (
                    (
                        Uf.oldTime().boundaryField()[patchI]
                      & mesh.Sf().boundaryField()[patchI]
                    )
                  - mesh.phi().oldTime().boundaryField()[patchI]
                )
            );
    }
}


OFstream& elasticWallPressureFvPatchScalarField::diagFile()
{
    if (diagFilePtr_.empty())
    {
        fileName dir;
        if (Pstream::parRun())
        {
            dir = db().time().path()/".."/"postProcessing";
        }
        else
        {
            dir = db().time().path()/"postProcessing";
        }

        mkDir(dir);
        diagFilePtr_.reset
        (
            new OFstream(dir/("robinCoefficient_" + patch().name() + ".dat"))
        );

        diagFilePtr_()
            << "# Time iteration alphaMean alphaMin alphaMax seedMean "
            << "secantScaleMean safeguardScale solidImpedance "
            << "fluidImpedance predictedRate observedRate deltaPNorm "
            << "leakFlux netLeakFlux wallMotionFlux carryOverFlux"
            << endl;
    }

    return diagFilePtr_();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void elasticWallPressureFvPatchScalarField::updateRobinCoefficient()
{
    const label timeIndex = db().time().timeIndex();
    const scalarField& magSf = patch().magSf();
    const scalar area = gSum(magSf);
    const bool secant = hsControls_.model == "secant";

    const scalarField seed(seedCoeff());
    const scalarField solidAccN(patch().nf() & prevAcceleration_);
    const scalarField fluidAccN(fluidNormalAcceleration());

    if (secantTimeIndex_ != timeIndex)
    {
        // Rescale using all FSI iterations of the previous time step
        if (secant && hsControls_.secantUpdate == "timeStep")
        {
            applySecantUpdate();
        }

        secantTimeIndex_ = timeIndex;
        secantIter_ = 0;
        secantLambdas_.clear();
        // Move the samples of the previous time step to the history and
        // forget those older than secantMemory time steps
        forAll(secantSamples_, sI)
        {
            secantSamplesHist_.append(secantSamples_[sI]);
            secantSamplesHistStep_.append(timeIndex - 1);
        }
        secantSamples_.clear();

        DynamicList<vector> keptSamples;
        DynamicList<label> keptSteps;
        forAll(secantSamplesHist_, sI)
        {
            if
            (
                secantSamplesHistStep_[sI]
             >= timeIndex - hsControls_.secantMemory
            )
            {
                keptSamples.append(secantSamplesHist_[sI]);
                keptSteps.append(secantSamplesHistStep_[sI]);
            }
        }
        secantSamplesHist_.transfer(keptSamples);
        secantSamplesHistStep_.transfer(keptSteps);
        secantNumField_ = 0.0;
        secantDenField_ = 0.0;
        prevDeltaPNorm_ = -1;
        nDiverging_ = 0;
    }
    else if (iterPressure_.size() == size())
    {
        // Changes between successive FSI iterations; the solid acceleration
        // is the response to the pressure applied in the same iteration
        const scalarField dp(prevPressure_ - iterPressure_);
        const scalarField da(solidAccN - iterSolidAccN_);
        const scalarField daF(fluidAccN - iterFluidAccN_);
        const scalarField alpha0da(seed*da);

        const scalar dpNorm = sqrt(gSum(magSf*sqr(dp))/area);
        const scalar daDen = gSum(magSf*sqr(da));
        const scalar daFDen = gSum(magSf*sqr(daF));
        const scalar aNorm = sqrt(gSum(magSf*sqr(solidAccN)));

        // Area-weighted Rayleigh-quotient estimates of the solid and fluid
        // interface impedances (added mass per unit area)
        const scalar solidImpedance =
            daDen > VSMALL ? gSum(magSf*dp*da)/daDen : 0;
        const scalar fluidImpedance =
            daFDen > VSMALL ? -gSum(magSf*dp*daF)/daFDen : 0;

        const scalar alphaMean = gSum(magSf*rhoSolidHs())/area;

        // Contraction factor predicted by the Robin-Neumann model
        scalar predictedRate = -1;
        if (solidImpedance > SMALL && fluidImpedance + alphaMean > SMALL)
        {
            predictedRate =
                (fluidImpedance/solidImpedance)
               *mag(solidImpedance - alphaMean)
               /(fluidImpedance + alphaMean);
        }

        const scalar observedRate =
            prevDeltaPNorm_ > VSMALL ? dpNorm/prevDeltaPNorm_ : -1;

        // Record the least-squares fit of this pair. The first pair of a time
        // step also contains the physical load change of the new step (e.g.
        // shear tractions, which move the solid without a normal pressure
        // change), so only pairs of subsequent FSI iterations are used.
        // Pairs dominated by round-off are ignored.
        const scalar num = gSum(magSf*dp*alpha0da);
        const scalar den = gSum(magSf*sqr(alpha0da));
        if
        (
            secantIter_ > 1
         && den > VSMALL
         && num > 0
         && sqrt(daDen) > 1e-10*(aNorm + VSMALL)
        )
        {
            secantLambdas_.append(num/den);

            // Impedance sample for the minimax fit, with the mean seed
            // coefficient weighted by the acceleration change of the pair
            if (solidImpedance > SMALL && fluidImpedance > SMALL)
            {
                secantSamples_.append
                (
                    vector
                    (
                        solidImpedance,
                        fluidImpedance,
                        gSum(magSf*seed*sqr(da))/daDen
                    )
                );
            }

            // Normalise so that each pair has the same weight
            secantNumField_ += dp*alpha0da*area/den;
            secantDenField_ += sqr(alpha0da)*area/den;
        }

        if
        (
            secant
         && hsControls_.secantUpdate == "iteration"
         && (secantLambdas_.size() || secantSamplesHist_.size())
        )
        {
            applySecantUpdate();
        }

        // Reduce the coefficient if the iterations diverge: in the
        // Robin-Neumann model this only happens when it is too large.
        // Only used for the a-priori models: the secant fit already reacts
        // to diverging modes, and halving on top of it was found to cause
        // oscillations
        if (hsControls_.divergenceSafeguard && !secant && observedRate > 1)
        {
            if (++nDiverging_ >= 2)
            {
                safeguardScale_ *= 0.5;

                Info<< type() << " " << patch().name()
                    << ": diverging FSI iterations: halving the Robin "
                    << "coefficient (safeguard scale = " << safeguardScale_
                    << ")" << endl;

                nDiverging_ = 0;
            }
        }
        else
        {
            nDiverging_ = 0;
        }

        prevDeltaPNorm_ = dpNorm;

        if (debug || hsControls_.writeDiagnostics)
        {
            Info<< type() << " " << patch().name() << ": iteration "
                << secantIter_ << ": solid impedance = " << solidImpedance
                << ", fluid impedance = " << fluidImpedance
                << ", alpha = " << alphaMean
                << ", predicted rate = " << predictedRate
                << ", observed rate = " << observedRate << endl;
        }

        if (hsControls_.writeDiagnostics)
        {
            const scalarField alpha(rhoSolidHs());
            const scalar alphaMin = gMin(alpha);
            const scalar alphaMax = gMax(alpha);
            const scalar seedMean = gSum(magSf*seed)/area;
            const scalar secantScaleMean = gSum(magSf*secantScale_)/area;

            scalar leakFlux = 0;
            scalar netLeakFlux = 0;
            scalar wallMotionFlux = 0;
            scalar carryOverFlux = 0;
            leakageFluxes
            (
                leakFlux, netLeakFlux, wallMotionFlux, carryOverFlux
            );

            if (Pstream::master())
            {
                diagFile()
                    << db().time().timeName() << " " << secantIter_ << " "
                    << alphaMean << " " << alphaMin << " " << alphaMax << " "
                    << seedMean << " " << secantScaleMean << " "
                    << safeguardScale_ << " " << solidImpedance << " "
                    << fluidImpedance << " " << predictedRate << " "
                    << observedRate << " " << dpNorm << " "
                    << leakFlux << " " << netLeakFlux << " "
                    << wallMotionFlux << " " << carryOverFlux << endl;
            }
        }
    }

    iterPressure_ = prevPressure_;
    iterSolidAccN_ = solidAccN;
    iterFluidAccN_ = fluidAccN;
    secantIter_++;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
