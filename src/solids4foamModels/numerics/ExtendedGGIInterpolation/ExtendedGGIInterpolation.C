/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Interpolation class dealing with transfer of data between two
    primitivePatches

Author
    Hrvoje Jasak, Wikki Ltd.

Contributor:
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "ExtendedGGIInterpolationTemplate.H"
#include "demandDrivenData.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FromPatch, class ToPatch>
void ExtendedGGIInterpolation<FromPatch, ToPatch>::
calcMasterPointAddressing() const
{
    // Find master points addressing
    if (masterPointAddressingPtr_)
    {
        FatalErrorIn
        (
            "void ExtendedGGIInterpolation::"
            "calcMasterPointAddressing() const"
        )
            << "Master points addressing already exists"
                << abort(FatalError);
    }

    masterPointAddressingPtr_ =
        new List<labelPair>
        (
            this->masterPatch().nPoints(),
            labelPair(-1,-1)
        );
    List<labelPair>& masterPointAddr = *masterPointAddressingPtr_;

    masterPointDistancePtr_ =
        new scalarField
        (
            this->masterPatch().nPoints(),
            GREAT
        );
    scalarField& masterPointDist = *masterPointDistancePtr_;

    const labelListList& masterFaceAddr = this->masterAddr();

    const labelListList& pointFaces = this->masterPatch().pointFaces();

    const faceList& slaveFaces = this->slavePatch().localFaces();
    const pointField& slavePoints = this->slavePatch().localPoints();

    const pointField& masterPoints = this->masterPatch().localPoints();

    forAll(masterPointAddr, pointI)
    {
        const point& P = masterPoints[pointI];

        labelHashSet possibleSlaveFacesSet;

        const labelList& curPointFaces = pointFaces[pointI];
        forAll(curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            const labelList& curSlaveFaces = masterFaceAddr[curFace];

            forAll(curSlaveFaces, fI)
            {
                if (!possibleSlaveFacesSet.found(curSlaveFaces[fI]))
                {
                    possibleSlaveFacesSet.insert(curSlaveFaces[fI]);
                }
            }
        }

        labelList possibleSlaveFaces = possibleSlaveFacesSet.toc();

        scalar MinEta = -GREAT;
        labelPair faceTriangle(-1, -1);
        scalar distance = GREAT;

        forAll(possibleSlaveFaces, faceI)
        {
            label curSlaveFace = possibleSlaveFaces[faceI];

            const face& f = slaveFaces[curSlaveFace];

            point ctr = Foam::average(f.points(slavePoints));

            point nextPoint = ctr;

            for (label pI = 0; pI < f.size(); pI++)
            {
                nextPoint = slavePoints[f.nextLabel(pI)];

                triPointRef t
                (
                    slavePoints[f[pI]],
                    nextPoint,
                    ctr
                );

                vector n = t.normal();
                scalar A = mag(n);
                n /= A;

                // Intersection point
                point I = P + n*(n&(t.a() - P));

                // Areal coordinates
                scalarField eta(3, 0);

                eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
                eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
                eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

                scalar minEta = min(eta);

                if (minEta > MinEta)
                {
                    MinEta = minEta;
                    faceTriangle.first() = curSlaveFace;
                    faceTriangle.second() = pI;

                    distance = ((P - I)&n);
                }
            }
        }

        masterPointAddr[pointI] = faceTriangle;
        masterPointDist[pointI] = distance;
    }
}

template<class FromPatch, class ToPatch>
void ExtendedGGIInterpolation<FromPatch, ToPatch>::
calcMasterPointWeights() const
{
    // Find master point weights
    if (masterPointWeightsPtr_)
    {
        FatalErrorIn
        (
            "void ExtendedGGIInterpolation::"
            "calcMasterPointAddressing() const"
        )
            << "Master point weights already exist"
                << abort(FatalError);
    }

    masterPointWeightsPtr_ =
        new FieldField<Field, scalar>(this->masterPatch().nPoints());
    FieldField<Field, scalar>& masterPointWeights = *masterPointWeightsPtr_;

    const faceList& slaveFaces = this->slavePatch().localFaces();
    const pointField& slavePoints = this->slavePatch().localPoints();

    const pointField& masterPoints = this->masterPatch().localPoints();

    const List<labelPair>& addr = this->masterPointAddr();

    forAll(masterPointWeights, pointI)
    {
        if (addr[pointI].first() != -1)
        {
            const point& P = masterPoints[pointI];

            const face& hitFace =
                slaveFaces[addr[pointI].first()];

            point ctr = Foam::average(hitFace.points(slavePoints));

            label pI = addr[pointI].second();

            triPointRef t
            (
                slavePoints[hitFace[pI]],
                slavePoints[hitFace.nextLabel(pI)],
                ctr
            );

            vector n = t.normal();
            n /= mag(n);

            // Intersection point
            point I = P + n*(n&(t.a() - P));

            masterPointWeights.set(pointI, scalarField(3));

            masterPointWeights[pointI][0] = t.Ni(0, I);
            masterPointWeights[pointI][1] = t.Ni(1, I);
            masterPointWeights[pointI][2] = t.Ni(2, I);
        }
        else
        {
            masterPointWeights.set(pointI, scalarField(0));
        }
    }
}

template<class FromPatch, class ToPatch>
void ExtendedGGIInterpolation<FromPatch, ToPatch>::
calcSlavePointAddressing() const
{
    // Find master points addressing
    if (slavePointAddressingPtr_)
    {
        FatalErrorIn
        (
            "void ExtendedGGIInterpolation::"
            "calcSlavePointAddressing() const"
        )
            << "Slave points addressing already exists"
                << abort(FatalError);
    }

    slavePointAddressingPtr_ =
        new List<labelPair>
        (
            this->slavePatch().nPoints(),
            labelPair(-1,-1)
        );
    List<labelPair>& slavePointAddr = *slavePointAddressingPtr_;

    slavePointDistancePtr_ =
        new scalarField
        (
            this->slavePatch().nPoints(),
            GREAT
        );
    scalarField& slavePointDist = *slavePointDistancePtr_;

    const labelListList& slaveFaceAddr = this->slaveAddr();

    const labelListList& pointFaces = this->slavePatch().pointFaces();

    const faceList& masterFaces = this->masterPatch().localFaces();
    const pointField& masterPoints = this->masterPatch().localPoints();

    const pointField& slavePoints = this->slavePatch().localPoints();

    forAll(slavePointAddr, pointI)
    {
        const point& P = slavePoints[pointI];

        labelHashSet possibleMasterFacesSet;

        const labelList& curPointFaces = pointFaces[pointI];
        forAll(curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            const labelList& curMasterFaces = slaveFaceAddr[curFace];

            forAll(curMasterFaces, fI)
            {
                if (!possibleMasterFacesSet.found(curMasterFaces[fI]))
                {
                    possibleMasterFacesSet.insert(curMasterFaces[fI]);
                }
            }
        }

        labelList possibleMasterFaces = possibleMasterFacesSet.toc();

        scalar MinEta = -GREAT;
        labelPair faceTriangle(-1, -1);
        scalar distance = GREAT;

        forAll(possibleMasterFaces, faceI)
        {
            label curMasterFace = possibleMasterFaces[faceI];

            const face& f = masterFaces[curMasterFace];

            point ctr = Foam::average(f.points(masterPoints));

            point nextPoint = ctr;

            for (label pI = 0; pI < f.size(); pI++)
            {
                nextPoint = masterPoints[f.nextLabel(pI)];

                triPointRef t
                (
                    masterPoints[f[pI]],
                    nextPoint,
                    ctr
                );

                vector n = t.normal();
                scalar A = mag(n);
                n /= A;

                // Intersection point
                point I = P + n*(n&(t.a() - P));

                // Areal coordinates
                scalarField eta(3, 0);

                eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
                eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
                eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

                scalar minEta = min(eta);

                if (minEta > MinEta)
                {
                    MinEta = minEta;
                    faceTriangle.first() = curMasterFace;
                    faceTriangle.second() = pI;

                    distance = ((P - I)&n);
                }
            }
        }

//         Info << "MinEta " << MinEta << endl;

        slavePointAddr[pointI] = faceTriangle;
        slavePointDist[pointI] = distance;
    }
}

template<class FromPatch, class ToPatch>
void ExtendedGGIInterpolation<FromPatch, ToPatch>::
calcSlavePointWeights() const
{
    // Find master point weights
    if (slavePointWeightsPtr_)
    {
        FatalErrorIn
        (
            "void ExtendedGGIInterpolation::"
            "calcSlavePointAddressing() const"
        )
            << "Slave point weights already exist"
                << abort(FatalError);
    }

    slavePointWeightsPtr_ =
        new FieldField<Field, scalar>(this->slavePatch().nPoints());
    FieldField<Field, scalar>& slavePointWeights = *slavePointWeightsPtr_;

    const faceList& masterFaces = this->masterPatch().localFaces();
    const pointField& masterPoints = this->masterPatch().localPoints();

    const pointField& slavePoints = this->slavePatch().localPoints();

    const List<labelPair>& addr = this->slavePointAddr();

    forAll(slavePointWeights, pointI)
    {
        if (addr[pointI].first() != -1)
        {
            const point& P = slavePoints[pointI];

            const face& hitFace =
                masterFaces[addr[pointI].first()];

            point ctr = Foam::average(hitFace.points(masterPoints));

            label pI = addr[pointI].second();

            triPointRef t
            (
                masterPoints[hitFace[pI]],
                masterPoints[hitFace.nextLabel(pI)],
                ctr
            );

            vector n = t.normal();
            n /= mag(n);

            // Intersection point
            point I = P + n*(n&(t.a() - P));

            slavePointWeights.set(pointI, scalarField(3));

            slavePointWeights[pointI][0] = t.Ni(0, I);
            slavePointWeights[pointI][1] = t.Ni(1, I);
            slavePointWeights[pointI][2] = t.Ni(2, I);
        }
        else
        {
            slavePointWeights.set(pointI, scalarField(0));
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class MasterPatch, class SlavePatch>
ExtendedGGIInterpolation<MasterPatch, SlavePatch>::ExtendedGGIInterpolation
(
    const MasterPatch& masterPatch,
    const SlavePatch&  slavePatch,
    const tensorField& forwardT,
    const tensorField& reverseT,
    const vectorField& forwardSep,
    const scalar masterNonOverlapFaceTol,
    const scalar slaveNonOverlapFaceTol,
    const bool rescaleGGIWeightingFactors,
    //const GGIInterpolationName::quickReject reject,
    const newGGIInterpolationName::quickReject reject,
    const boundBox& regionOfInterest
)
:
    //GGIInterpolation<MasterPatch, SlavePatch>
    newGGIInterpolation<MasterPatch, SlavePatch>
    (
        masterPatch,
        slavePatch,
        forwardT,
        reverseT,
        forwardSep,
        masterNonOverlapFaceTol,
        slaveNonOverlapFaceTol,
        rescaleGGIWeightingFactors,
        reject,
        regionOfInterest
    ),
    masterPointAddressingPtr_(NULL),
    masterPointWeightsPtr_(NULL),
    masterPointDistancePtr_(NULL),
    slavePointAddressingPtr_(NULL),
    slavePointWeightsPtr_(NULL),
    slavePointDistancePtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
ExtendedGGIInterpolation<MasterPatch, SlavePatch>::~ExtendedGGIInterpolation()
{
    deleteDemandDrivenData(masterPointAddressingPtr_);
    deleteDemandDrivenData(masterPointWeightsPtr_);
    deleteDemandDrivenData(masterPointDistancePtr_);

    deleteDemandDrivenData(slavePointAddressingPtr_);
    deleteDemandDrivenData(slavePointWeightsPtr_);
    deleteDemandDrivenData(slavePointDistancePtr_);

    // not needed as automatically called
    // GGIInterpolation<MasterPatch, SlavePatch>::~GGIInterpolation();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class MasterPatch, class SlavePatch>
const Foam::List<labelPair>&
ExtendedGGIInterpolation<MasterPatch, SlavePatch>::masterPointAddr() const
{
    if (!masterPointAddressingPtr_)
    {
        calcMasterPointAddressing();
    }

    return *masterPointAddressingPtr_;
}

template<class MasterPatch, class SlavePatch>
const Foam::FieldField<Field, scalar>&
ExtendedGGIInterpolation<MasterPatch, SlavePatch>::masterPointWeights() const
{
    if (!masterPointWeightsPtr_)
    {
        calcMasterPointWeights();
    }

    return *masterPointWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarField&
ExtendedGGIInterpolation<MasterPatch, SlavePatch>
::masterPointDistanceToIntersection() const
{
    if (!masterPointDistancePtr_)
    {
        calcMasterPointAddressing();
    }

    return *masterPointDistancePtr_;
}


template<class MasterPatch, class SlavePatch>
const Foam::List<labelPair>&
ExtendedGGIInterpolation<MasterPatch, SlavePatch>::slavePointAddr() const
{
    if (!slavePointAddressingPtr_)
    {
        calcSlavePointAddressing();
    }

    return *slavePointAddressingPtr_;
}

template<class MasterPatch, class SlavePatch>
const Foam::FieldField<Field, scalar>&
ExtendedGGIInterpolation<MasterPatch, SlavePatch>::slavePointWeights() const
{
    if (!slavePointWeightsPtr_)
    {
        calcSlavePointWeights();
    }

    return *slavePointWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarField&
ExtendedGGIInterpolation<MasterPatch, SlavePatch>
::slavePointDistanceToIntersection() const
{
    if (!slavePointDistancePtr_)
    {
        calcSlavePointAddressing();
    }

    return *slavePointDistancePtr_;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> > ExtendedGGIInterpolation<MasterPatch, SlavePatch>::
slaveToMasterPointInterpolate
(
    const Field<Type>& pf
) const
{
    if (pf.size() != this->slavePatch().nPoints())
    {
        FatalErrorIn
        (
            "ExtendedGGIInterpolation::slaveToMasterPointInterpolate"
            "(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << this->slavePatch().nPoints() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            this->masterPatch().nPoints(),
            pTraits<Type>::zero
        )
    );

    // Escape the interpolation if there are no faces in the target patch
    if (this->masterPatch().nPoints() == 0)
    {
        return tresult;
    }

    Field<Type>& result = tresult();

    const List<typename SlavePatch::FaceType>& slaveFaces =
        this->slavePatch().localFaces();

    const List<labelPair>& addr = masterPointAddr();

    const FieldField<Field, scalar>& weights = masterPointWeights();

    forAll (result, pointI)
    {
        if (addr[pointI].first() > -1)
        {
            const face& hitFace =
                slaveFaces[addr[pointI].first()];

            label pI = addr[pointI].second();

            Type ctrF = average(Field<Type>(pf, hitFace));

            result[pointI] =
                weights[pointI][0]*pf[hitFace[pI]]
              + weights[pointI][1]*pf[hitFace.nextLabel(pI)]
              + weights[pointI][2]*ctrF;
        }
    }

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> > ExtendedGGIInterpolation<MasterPatch, SlavePatch>::
masterToSlavePointInterpolate
(
    const Field<Type>& pf
) const
{
    if (pf.size() != this->masterPatch().nPoints())
    {
        FatalErrorIn
        (
            "ExtendedGGIInterpolation::masterToSlavePointInterpolate"
            "(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << this->masterPatch().nPoints() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            this->slavePatch().nPoints(),
            pTraits<Type>::zero
        )
    );

    // Escape the interpolation if there are no faces in the target patch
    if (this->slavePatch().nPoints() == 0)
    {
        return tresult;
    }

    Field<Type>& result = tresult();

    const List<typename SlavePatch::FaceType>& masterFaces =
        this->masterPatch().localFaces();

    const List<labelPair>& addr = slavePointAddr();

    const FieldField<Field, scalar>& weights = slavePointWeights();

    forAll (result, pointI)
    {
        if (addr[pointI].first() > -1)
        {
            const face& hitFace =
                masterFaces[addr[pointI].first()];

            label pI = addr[pointI].second();

            Type ctrF = average(Field<Type>(pf, hitFace));

            result[pointI] =
                weights[pointI][0]*pf[hitFace[pI]]
              + weights[pointI][1]*pf[hitFace.nextLabel(pI)]
              + weights[pointI][2]*ctrF;
        }
    }

    return tresult;
}


template<class MasterPatch, class SlavePatch>
// bool ExtendedGGIInterpolation<MasterPatch, SlavePatch>::movePoints()
bool ExtendedGGIInterpolation<MasterPatch, SlavePatch>::movePoints
(
    const tensorField& forwardT,
    const tensorField& reverseT,
    const vectorField& forwardSep
)
{
    deleteDemandDrivenData(masterPointAddressingPtr_);
    deleteDemandDrivenData(masterPointWeightsPtr_);
    deleteDemandDrivenData(masterPointDistancePtr_);

    deleteDemandDrivenData(slavePointAddressingPtr_);
    deleteDemandDrivenData(slavePointWeightsPtr_);
    deleteDemandDrivenData(slavePointDistancePtr_);

    //return GGIInterpolation<MasterPatch, SlavePatch>::movePoints
    return newGGIInterpolation<MasterPatch, SlavePatch>::movePoints
        (
            forwardT,
            reverseT,
            forwardSep
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
