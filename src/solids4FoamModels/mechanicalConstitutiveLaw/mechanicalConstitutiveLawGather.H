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
    gatherToIntegrationPoints: fill a law's integration points from a volField,
    as the mechanicalConstitutiveLawManager does for a prescribed state and a
    coupling input; and headerIsA, a file header check that works on every
    fork. Free functions, so that the manager and its helpers can share them.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawGather_H
#define mechanicalConstitutiveLawGather_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceInterpolate.H"
#include "linear.H"
#include "integrationPointTopology.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Whether a file for this IOobject exists and holds this type.
//  foam-extend has no typed header check, so the class name is compared after
//  reading the header, which is what the typed check does anyway
template<class FieldType>
inline bool mechanicalConstitutiveLawHeaderIsA(IOobject& io)
{
#ifdef FOAMEXTEND
    return io.headerOk() && io.headerClassName() == FieldType::typeName;
#else
    return io.typeHeaderOk<FieldType>(true);
#endif
}


//- Fill a law's integration points from a volField on the mesh.
//
//  The value is specified per cell, which is how a user prepares one
//  and how another model produces one. How it reaches a point depends
//  on where the point is, as the topology says through
//  integrationPointFace.
//
//  A point on a mesh face takes the interpolated face value, because a
//  face lies between two cells rather than in one. That holds for a
//  face-centred topology and for a compact face one alike, and on a
//  processor face too: the interpolation there uses the neighbour
//  values the rank already holds, so a decomposed run sees what a
//  serial one does without anything being exchanged.
//
//  A point inside a cell takes that cell's value. A point shared
//  between cells without being a face, as a point-centred one is,
//  takes the average of this law's cells around it on this rank; the
//  pointField overload, whose points these are, refuses a parallel run.
//  The average is independent of the order the cells are visited in up
//  to rounding, and exactly so for two cells.
//
//  Every type is treated alike. A direction, such as a fibre, is
//  therefore interpolated too, and a law reading one normalises it
template<class Type>
inline void gatherToIntegrationPoints
(
    const fvMesh& mesh,
    const labelList& lawCells,
    const GeometricField<Type, fvPatchField, volMesh>& src,
    const integrationPointTopology& topo,
    const labelList& ipIDs,
    Field<Type>& fld
)
{
    // The face each point lies on, or -1 for a point that lies in a cell
    const label nInternal = mesh.nInternalFaces();

    labelList ipFace(ipIDs.size(), -1);
    bool anyFace = false;

    forAll(ipIDs, localI)
    {
        ipFace[localI] = topo.integrationPointFace(ipIDs[localI]);

        if (ipFace[localI] >= 0)
        {
            anyFace = true;
        }
    }

    // A state field is indexed by the law's own integration points,
    // 0 to n-1, while the topology hands out global indices, so the
    // two have to be related before anything can be written
    labelList globalToLocal(topo.nIntegrationPoints(), -1);
    forAll(ipIDs, localI)
    {
        globalToLocal[ipIDs[localI]] = localI;
    }

    // The points in cells: the average of this law's cells around each,
    // which for a point inside one cell is that cell's value
    Field<Type> sum(fld.size(), pTraits<Type>::zero);
    labelList count(fld.size(), 0);

    forAll(lawCells, i)
    {
        const label cellI = lawCells[i];
        const labelUList ips = topo.cellIntegrationPointIDs(cellI);

        forAll(ips, j)
        {
            const label localI = globalToLocal[ips[j]];

            // A point of this cell that is not one of this law's own
            // points cannot happen while materials are assigned per
            // cell, but the guard costs nothing and the alternative is
            // writing outside the field
            if (localI >= 0 && ipFace[localI] < 0)
            {
                sum[localI] += src[cellI];
                ++count[localI];
            }
        }
    }

    forAll(fld, localI)
    {
        if (count[localI] > 0)
        {
            fld[localI] = sum[localI]/scalar(count[localI]);
        }
    }

    if (!anyFace)
    {
        return;
    }

    // The points on faces
    const GeometricField<Type, fvsPatchField, surfaceMesh> srcf
    (
        linearInterpolate(src)
    );

    forAll(ipIDs, localI)
    {
        const label faceI = ipFace[localI];

        if (faceI < 0)
        {
            continue;
        }

        if (faceI < nInternal)
        {
            fld[localI] = srcf[faceI];
            continue;
        }

        const label patchI = mesh.boundaryMesh().whichPatch(faceI);
        const fvsPatchField<Type>& pf = srcf.boundaryField()[patchI];

        if (pf.empty())
        {
            // A face of an empty patch holds no value, and has one cell
            fld[localI] = src[mesh.faceOwner()[faceI]];
        }
        else
        {
            fld[localI] = pf[faceI - mesh.boundaryMesh()[patchI].start()];
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
