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
//  on where the point is. A point inside a cell takes that cell's
//  value. A point that is a mesh face takes the interpolated value,
//  because a face lies between two cells rather than in one.
//
//  A point shared between cells without being a face, as a
//  point-centred one is, takes the average of this law's cells around
//  it, which does not depend on the order they are visited in
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
    if (topo.integrationPointsAreFaces())
    {
        // The integration point index is the face index here, which is
        // what lets the interpolated field be read straight off
        const GeometricField<Type, fvsPatchField, surfaceMesh> srcf
        (
            linearInterpolate(src)
        );

        const label nInternal = mesh.nInternalFaces();

        forAll(ipIDs, localI)
        {
            const label faceI = ipIDs[localI];

            if (faceI < nInternal)
            {
                fld[localI] = srcf[faceI];
            }
            else
            {
                const label patchI =
                    mesh.boundaryMesh().whichPatch(faceI);

                const polyPatch& pp = mesh.boundaryMesh()[patchI];

                fld[localI] =
                    srcf.boundaryField()[patchI][faceI - pp.start()];
            }
        }

        return;
    }

    // A state field is indexed by the law's own integration points,
    // 0 to n-1, while the topology hands out global indices, so the
    // two have to be related before anything can be written
    labelList globalToLocal(topo.nIntegrationPoints(), -1);
    forAll(ipIDs, localI)
    {
        globalToLocal[ipIDs[localI]] = localI;
    }

    const labelList& cells = lawCells;

    // A point shared between cells, as a point-centred one is, takes
    // the average of this law's cells around it, so that the value
    // does not depend on the order the cells are visited in. A point
    // inside one cell has one contribution, and takes that cell's
    Field<Type> sum(fld.size(), pTraits<Type>::zero);
    labelList count(fld.size(), 0);

    forAll(cells, i)
    {
        const label cellI = cells[i];
        const labelUList ips = topo.cellIntegrationPointIDs(cellI);

        forAll(ips, j)
        {
            const label localI = globalToLocal[ips[j]];

            // A point of this cell that is not one of this law's own
            // points cannot happen while materials are assigned per
            // cell, but the guard costs nothing and the alternative is
            // writing outside the field
            if (localI >= 0)
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
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
