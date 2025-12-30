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

#include "cellZoneInterface.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    //- Return a boolList indicating, for each internal face, if it
    //  straddles two cell zones (true) or not (false)
    tmp<Field<bool>> cellZoneInterface(const fvMesh& mesh, const bool debug)
    {
        // Create a field indicating internal faces on bi-material interfaces; that
        // is, they straddle two contiguous cell zones

        // Set local references to mesh data
        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        // Prepare result field
        tmp<Field<bool>> tinterface(new Field<bool>(owner.size(), false));
        Field<bool>& interface = tmpRef(tinterface);

        if (mesh.cellZones().size() > 1)
        {
            labelList cellZoneID(mesh.nCells(), -1);
            forAll(mesh.cellZones(), czI)
            {
                const labelList& curCellZone = mesh.cellZones()[czI];

                forAll(curCellZone, cI)
                {
                    const label cellID = curCellZone[cI];

                    if (cellZoneID[cellID] != -1)
                    {
                        FatalErrorInFunction
                            << "Cell " << cellID << " is in more than on cell zone!"
                            << nl << "It is in cells zones " << cellZoneID[cellID]
                            << " and " << cI << exit(FatalError);
                    }

                    cellZoneID[cellID] = czI;
                }
            }

            // Check all cells are in at least one cell zone
            if (gMin(cellZoneID) == -1)
            {
                FatalErrorInFunction
                    << "There is at least one cell not in a cell zone!"
                    << exit(FatalError);
            }

            // Set the interface field
            label nInterfaceFaces = 0;
            forAll(interface, faceI)
            {
                const label own = owner[faceI];
                const label nei = neighbour[faceI];

                if (cellZoneID[own] != cellZoneID[nei])
                {
                    interface[faceI] = true;
                    nInterfaceFaces++;
                }
            }

            if (debug)
            {
                InfoInFunction
                    << nl << "There are " << nInterfaceFaces << " faces on an interface"
                    << nl << endl;
            }

            if (Pstream::parRun())
            {
                notImplemented
                (
                    "Multiple cell-zone gradient calculation not yet implemented"
                    " for parallel"
                )
            }
        }

        return tinterface;
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
