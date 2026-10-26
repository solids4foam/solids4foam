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

#ifndef cartesianMeshInfo_H
#define cartesianMeshInfo_H

#include "fvMesh.H"
#include "compatibilityFunctions.H"
#include "boundBox.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cartesianMeshInfo Declaration
\*---------------------------------------------------------------------------*/

struct cartesianMeshInfo
{
    //- Mesh spacings and counts, including the one-cell empty direction
    vector h;
    label N[3];

    //- Discrete cell index in each direction
    labelList index[3];

    //- Construct and validate a uniform, axis-aligned periodic Cartesian mesh
    explicit cartesianMeshInfo(const fvMesh& mesh)
    :
        h(vector::zero)
    {
        const boundBox bounds(mesh.points(), false);
        const vector span(bounds.span());
        const Vector<label> active(mesh.geometricD());
        bool periodic[3] = {false, false, false};

        const auto checkFace = [&](const vector& area, const scalar delta)
        {
            const vector normal(area/mag(area));
            direction axis = 0;
            for (direction d = 1; d < 3; ++d)
            {
                if (mag(normal[d]) > mag(normal[axis])) axis = d;
            }
            if (mag(normal[axis]) < 1 - 1e-10 || delta <= 0)
            {
                FatalErrorInFunction << "Non-Cartesian face or invalid delta"
                    << exit(FatalError);
            }
            const scalar spacing = 1/delta;
            if (h[axis] == 0) h[axis] = spacing;
            if (mag(h[axis] - spacing) > 1e-12*h[axis])
            {
                FatalErrorInFunction << "Non-uniform spacing in axis " << axis
                    << exit(FatalError);
            }
            return axis;
        };

        forAll(primitiveField(mesh.Sf()), facei)
        {
            checkFace(mesh.Sf()[facei], mesh.deltaCoeffs()[facei]);
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& patch = mesh.boundary()[patchi];
            if (patch.type() != "cyclic" && patch.type() != "empty")
            {
                FatalErrorInFunction << "Only cyclic and empty patches supported"
                    << exit(FatalError);
            }
            if (patch.type() == "empty") continue;
            forAll(mesh.Sf().boundaryField()[patchi], facei)
            {
                const direction axis = checkFace
                (
                    mesh.Sf().boundaryField()[patchi][facei],
                    mesh.deltaCoeffs().boundaryField()[patchi][facei]
                );
                periodic[axis] = true;
            }
        }

        label count = 1;
        scalar volume = 1;
        for (direction d = 0; d < 3; ++d)
        {
            if (active[d] == -1)
            {
                h[d] = span[d];
                N[d] = 1;
            }
            else
            {
                if (!periodic[d] || h[d] <= 0)
                {
                    FatalErrorInFunction << "Axis " << d << " is not periodic"
                        << exit(FatalError);
                }
                N[d] = std::lround(span[d]/h[d]);
            }
            if (N[d] < 1 || mag(span[d] - N[d]*h[d]) > 1e-12*span[d])
            {
                FatalErrorInFunction << "Invalid Cartesian extent in axis " << d
                    << exit(FatalError);
            }
            count *= N[d];
            volume *= h[d];
            index[d].setSize(mesh.nCells());
        }
        if (count != mesh.nCells())
        {
            FatalErrorInFunction << "Cell count does not match Cartesian lattice"
                << exit(FatalError);
        }

        boolList seen(count, false);
        forAll(primitiveField(mesh.C()), celli)
        {
            if
            (
                mesh.cells()[celli].size() != 6
             || mag(mesh.V()[celli] - volume) > 1e-12*volume
            )
            {
                FatalErrorInFunction << "Expected uniform hexahedral volumes"
                    << exit(FatalError);
            }
            for (direction d = 0; d < 3; ++d)
            {
                const scalar x = (mesh.C()[celli][d] - bounds.min()[d])/h[d] - 0.5;
                const label i = std::lround(x);
                if (i < 0 || i >= N[d] || mag(x - i) > 1e-10)
                {
                    FatalErrorInFunction << "Cell centre is off the lattice"
                        << exit(FatalError);
                }
                index[d][celli] = i;
            }
            const label flat = index[0][celli]
                + N[0]*(index[1][celli] + N[1]*index[2][celli]);
            if (seen[flat])
            {
                FatalErrorInFunction << "Duplicate Cartesian cell index"
                    << exit(FatalError);
            }
            seen[flat] = true;
        }
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
