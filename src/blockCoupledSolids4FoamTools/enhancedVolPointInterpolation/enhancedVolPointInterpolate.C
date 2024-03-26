/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "enhancedVolPointInterpolation.H"
#include "volFields.H"
#include "pointFields.H"
#include "emptyFvPatch.H"
#include "coupledPointPatchField.H"
#include "pointConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::enhancedVolPointInterpolation::pushUntransformedData
(
    List<Type>& pointData
) const
{
    // Transfer onto coupled patch
    const globalMeshData& gmd = mesh().globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    // Combine master data with slave data
    forAll(slaves, i)
    {
        const labelList& slavePoints = slaves[i];

        // Copy master data to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elems[i];
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}


template<class Type>
void Foam::enhancedVolPointInterpolation::addSeparated
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "enhancedVolPointInterpolation::addSeparated" << endl;
    }

    typename GeometricField<Type, pointPatchField, pointMesh>::
        Internal& pfi = pf.ref();

    typename GeometricField<Type, pointPatchField, pointMesh>::
        Boundary& pfbf = pf.boundaryFieldRef();

    const label nReq = Pstream::nRequests();

    forAll(pfbf, patchi)
    {
        if (pfbf[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type>>
                (pfbf[patchi]).initSwapAddSeparated
                (
                    Pstream::commsTypes::nonBlocking,
                    pfi
                );
        }
    }

    // Block for any outstanding requests
    Pstream::waitRequests(nReq);

    forAll(pfbf, patchi)
    {
        if (pfbf[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type>>
                (pfbf[patchi]).swapAddSeparated
                (
                    Pstream::commsTypes::nonBlocking,
                    pfi
                );
        }
    }
}


template<class Type>
void Foam::enhancedVolPointInterpolation::interpolateInternalField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "enhancedVolPointInterpolation::interpolateInternalField("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field " << vf.name()
            << " from cells to points " << pf.name() << endl;
    }

    const labelListList& pointCells = vf.mesh().pointCells();

    // Multiply volField by weighting factor matrix to create pointField
    forAll(pointCells, pointi)
    {
        if (!isPatchPoint_[pointi])
        {
            const scalarList& pw = pointWeights_[pointi];
            const labelList& ppc = pointCells[pointi];

            pf[pointi] = Zero;

            forAll(ppc, pointCelli)
            {
                pf[pointi] += pw[pointCelli]*vf[ppc[pointCelli]];
            }
        }
    }
}


template<class Type>
void Foam::enhancedVolPointInterpolation::interpolateDimensionedInternalField
(
    const DimensionedField<Type, volMesh>& vf,
    DimensionedField<Type, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "enhancedVolPointInterpolation::interpolateDimensionedInternalField("
            << "const DimensionedField<Type, volMesh>&, "
            << "DimensionedField<Type, pointMesh>&) : "
            << "interpolating field " << vf.name() << " from cells to points "
            << pf.name() << endl;
    }

    const fvMesh& mesh = vf.mesh();

    const labelListList& pointCells = mesh.pointCells();
    const pointField& points = mesh.points();
    const vectorField& cellCentres = mesh.cellCentres();

    // Re-do weights and interpolation since normal interpolation
    // pointWeights_ are for non-boundary points only. Not efficient but
    // then saves on space.

    // Multiply volField by weighting factor matrix to create pointField
    scalarField sumW(points.size(), Zero);
    forAll(pointCells, pointi)
    {
        const labelList& ppc = pointCells[pointi];

        pf[pointi] = Type(Zero);

        forAll(ppc, pointCelli)
        {
            label celli = ppc[pointCelli];
            scalar pw = 1.0/mag(points[pointi] - cellCentres[celli]);

            pf[pointi] += pw*vf[celli];
            sumW[pointi] += pw;
        }
    }

    // Sum collocated contributions
    pointConstraints::syncUntransformedData(mesh, sumW, plusEqOp<scalar>());
    pointConstraints::syncUntransformedData(mesh, pf, plusEqOp<Type>());

    // Normalise
    forAll(pf, pointi)
    {
        scalar s = sumW[pointi];
        if (s > ROOTVSMALL)
        {
            pf[pointi] /= s;
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::enhancedVolPointInterpolation::flatBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = vf.mesh();
    const fvBoundaryMesh& bm = mesh.boundary();

    tmp<Field<Type>> tboundaryVals
    (
        new Field<Type>(mesh.nBoundaryFaces())
    );
    Field<Type>& boundaryVals = tboundaryVals.ref();

    forAll(vf.boundaryField(), patchi)
    {
        label bFacei = bm[patchi].patch().start() - mesh.nInternalFaces();

        if
        (
           !isA<emptyFvPatch>(bm[patchi])
        && !vf.boundaryField()[patchi].coupled()
        )
        {
            SubList<Type>
            (
                boundaryVals,
                vf.boundaryField()[patchi].size(),
                bFacei
            ) = vf.boundaryField()[patchi];
        }
        else
        {
            const polyPatch& pp = bm[patchi].patch();

            forAll(pp, i)
            {
                boundaryVals[bFacei++] = Zero;
            }
        }
    }

    return tboundaryVals;
}


template<class Type>
void Foam::enhancedVolPointInterpolation::interpolateBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    const primitivePatch& boundary = boundaryPtr_();

    Field<Type>& pfi = pf.primitiveFieldRef();

    // Get face data in flat list
    tmp<Field<Type>> tboundaryVals(flatBoundaryField(vf));
    const Field<Type>& boundaryVals = tboundaryVals();


    // Do points on 'normal' patches from the surrounding patch faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const labelList& mp = boundary.meshPoints();

    forAll(mp, i)
    {
        label pointi = mp[i];

        if (isPatchPoint_[pointi])
        {
            const labelList& pFaces = boundary.pointFaces()[i];
            const scalarList& pWeights = boundaryPointWeights_[i];

            Type& val = pfi[pointi];

            val = Zero;
            forAll(pFaces, j)
            {
                if (boundaryIsPatchFace_[pFaces[j]])
                {
                    val += pWeights[j]*boundaryVals[pFaces[j]];
                }
            }
        }
    }

    // Sum collocated contributions
    pointConstraints::syncUntransformedData(mesh(), pfi, plusEqOp<Type>());

    // And add separated contributions
    addSeparated(pf);

    // Optionally normalise
    if (normalisationPtr_)
    {
        const scalarField& normalisation = normalisationPtr_();
        forAll(mp, i)
        {
            pfi[mp[i]] *= normalisation[i];
        }
    }


    // Push master data to slaves. It is possible (not sure how often) for
    // a coupled point to have its master on a different patch so
    // to make sure just push master data to slaves.
    pushUntransformedData(pfi);
}


template<class Type>
void Foam::enhancedVolPointInterpolation::interpolateBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf,
    const bool overrideFixedValue
) const
{
    interpolateBoundaryField(vf, pf);

    // Apply constraints
    const pointConstraints& pcs = pointConstraints::New(pf.mesh());

    pcs.constrain(pf, overrideFixedValue);
}


template<class Type>
void Foam::enhancedVolPointInterpolation::interpolate
(
    const DimensionedField<Type, volMesh>& vf,
    const DimensionedField
    <
        typename outerProduct<vector, Type>::type,
        volMesh
    >& gradVf,
    DimensionedField<Type, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "enhancedVolPointInterpolation::interpolateDimensionedInternalField("
            << "const DimensionedField<Type, volMesh>&, "
            << "const DimensionedField<GradType, volMesh>&, "
            << "DimensionedField<Type, pointMesh>&) : "
            << "interpolating field " << vf.name() << " from cells to points "
            << pf.name() << endl;
    }

    const fvMesh& mesh = vf.mesh();

    const labelListList& pointCells = mesh.pointCells();
    const pointField& points = mesh.points();
    const vectorField& cellCentres = mesh.cellCentres();

    // Re-do weights and interpolation since normal interpolation
    // pointWeights_ are for non-boundary points only. Not efficient but
    // then saves on space.

    // Multiply volField by weighting factor matrix to create pointField
    scalarField sumW(points.size(), Zero);
    forAll(pointCells, pointi)
    {
        const labelList& ppc = pointCells[pointi];

        pf[pointi] = Type(Zero);

        forAll(ppc, pointCelli)
        {
            const label celli = ppc[pointCelli];

            // Distance vector
            const vector delta = points[pointi] - cellCentres[celli];

            // Inverse distance weight
            const scalar pw = 1.0/mag(delta);
            // const scalar pw = 1.0/mag(points[pointi] - cellCentres[celli]);

            // Weighted extrapolated value
            pf[pointi] += pw*(vf[celli] + (delta & gradVf[celli]));
            // pf[pointi] += pw*vf[celli];
            sumW[pointi] += pw;
        }
    }

    // Sum collocated contributions
    pointConstraints::syncUntransformedData(mesh, sumW, plusEqOp<scalar>());
    pointConstraints::syncUntransformedData(mesh, pf, plusEqOp<Type>());

    // Normalise
    forAll(pf, pointi)
    {
        scalar s = sumW[pointi];
        if (s > ROOTVSMALL)
        {
            pf[pointi] /= s;
        }
    }
}


template<class Type>
void Foam::enhancedVolPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "enhancedVolPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field " << vf.name() <<  " from cells to points "
            << pf.name() << endl;
    }

    interpolateInternalField(vf, pf);

    // Interpolate to the patches preserving fixed value BCs
    interpolateBoundaryField(vf, pf, false);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::enhancedVolPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const wordList& patchFieldTypes
) const
{
    const pointMesh& pm = pointMesh::New(vf.mesh());

    // Construct tmp<pointField>
    auto tpf = tmp<GeometricField<Type, pointPatchField, pointMesh>>::New
    (
        IOobject
        (
            "enhancedVolPointInterpolate(" + vf.name() + ')',
            vf.instance(),
            pm.thisDb()
        ),
        pm,
        vf.dimensions(),
        patchFieldTypes
    );

    interpolateInternalField(vf, tpf.ref());

    // Interpolate to the patches overriding fixed value BCs
    interpolateBoundaryField(vf, tpf.ref(), true);

    return tpf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::enhancedVolPointInterpolation::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    const wordList& patchFieldTypes
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh>> tpf =
        interpolate(tvf(), patchFieldTypes);
    tvf.clear();
    return tpf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::enhancedVolPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name,
    const bool cache
) const
{
    typedef GeometricField<Type, pointPatchField, pointMesh> PointFieldType;

    const pointMesh& pm = pointMesh::New(vf.mesh());
    const objectRegistry& db = pm.thisDb();

    PointFieldType* pfPtr =
        db.objectRegistry::template getObjectPtr<PointFieldType>(name);

    if (!cache || vf.mesh().changing())
    {
        // Delete any old occurrences to avoid double registration
        if (pfPtr && pfPtr->ownedByRegistry())
        {
            solution::cachePrintMessage("Deleting", name, vf);
            delete pfPtr;
        }

        auto tpf = tmp<GeometricField<Type, pointPatchField, pointMesh>>::New
        (
            IOobject
            (
                name,
                vf.instance(),
                pm.thisDb()
            ),
            pm,
            vf.dimensions()
        );

        interpolate(vf, tpf.ref());

        return tpf;
    }


    if (!pfPtr)
    {
        solution::cachePrintMessage("Calculating and caching", name, vf);

        pfPtr = interpolate(vf, name, false).ptr();
        regIOobject::store(pfPtr);
    }
    else
    {
        PointFieldType& pf = *pfPtr;

        if (pf.upToDate(vf))    //TBD: , vf.mesh().points()))
        {
            solution::cachePrintMessage("Reusing", name, vf);
        }
        else
        {
            solution::cachePrintMessage("Updating", name, vf);
            interpolate(vf, pf);
        }
    }

    return *pfPtr;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::enhancedVolPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return interpolate(vf, "enhancedVolPointInterpolate(" + vf.name() + ')', false);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::enhancedVolPointInterpolation::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh>> tpf =
        interpolate(tvf());
    tvf.clear();
    return tpf;
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::pointMesh>>
Foam::enhancedVolPointInterpolation::interpolate
(
    const DimensionedField<Type, volMesh>& vf,
    const word& name,
    const bool cache
) const
{
    typedef DimensionedField<Type, pointMesh> PointFieldType;

    const pointMesh& pm = pointMesh::New(vf.mesh());
    const objectRegistry& db = pm.thisDb();


    PointFieldType* pfPtr =
        db.objectRegistry::template getObjectPtr<PointFieldType>(name);

    if (!cache || vf.mesh().changing())
    {
        // Delete any old occurrences to avoid double registration
        if (pfPtr && pfPtr->ownedByRegistry())
        {
            solution::cachePrintMessage("Deleting", name, vf);
            delete pfPtr;
        }

        auto tpf = tmp<DimensionedField<Type, pointMesh>>::New
        (
            IOobject
            (
                name,
                vf.instance(),
                pm.thisDb()
            ),
            pm,
            vf.dimensions()
        );

        interpolateDimensionedInternalField(vf, tpf.ref());

        return tpf;
    }


    if (!pfPtr)
    {
        solution::cachePrintMessage("Calculating and caching", name, vf);
        pfPtr = interpolate(vf, name, false).ptr();

        regIOobject::store(pfPtr);
    }
    else
    {
        PointFieldType& pf = *pfPtr;

        if (pf.upToDate(vf))    //TBD: , vf.mesh().points()))
        {
            solution::cachePrintMessage("Reusing", name, vf);
        }
        else
        {
            solution::cachePrintMessage("Updating", name, vf);
            interpolateDimensionedInternalField(vf, pf);
        }
    }

    return *pfPtr;
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::pointMesh>>
Foam::enhancedVolPointInterpolation::interpolate
(
    const DimensionedField<Type, volMesh>& vf
) const
{
    return interpolate(vf, "enhancedVolPointInterpolate(" + vf.name() + ')', false);
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::pointMesh>>
Foam::enhancedVolPointInterpolation::interpolate
(
    const tmp<DimensionedField<Type, volMesh>>& tvf
) const
{
    // Construct tmp<pointField>
    tmp<DimensionedField<Type, pointMesh>> tpf = interpolate(tvf.ref());
    tvf.clear();
    return tpf;
}


// ************************************************************************* //
