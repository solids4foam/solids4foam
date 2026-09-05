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
#include "fvMesh.H"
#include "IOdictionary.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    labelList cellMaterialID(const fvMesh& mesh)
    {
        labelList materialID(mesh.nCells(), -1);

        // The materials are the entries of the "mechanical" list in
        // mechanicalProperties, each named after the cell zone it occupies.
        // Looked up rather than read, because the mechanical model registers
        // it and reading it again on every gradient would be wasteful
        if (!mesh.foundObject<IOdictionary>("mechanicalProperties"))
        {
            // Nothing says otherwise, so treat the mesh as one material
            return materialID;
        }

        const IOdictionary& mechProps =
            mesh.lookupObject<IOdictionary>("mechanicalProperties");

        if (!mechProps.found("mechanical"))
        {
            return materialID;
        }

        const PtrList<entry> lawEntries(mechProps.lookup("mechanical"));

        if (lawEntries.size() < 2)
        {
            // A single material fills the mesh whatever zones it carries, so
            // there is no interface anywhere and every cell stays -1
            return materialID;
        }

        forAll(lawEntries, lawI)
        {
            const word& materialName = lawEntries[lawI].keyword();

            const label zoneI = mesh.cellZones().findZoneID(materialName);

            if (zoneI == -1)
            {
                FatalErrorInFunction
                    << "Material '" << materialName << "' in "
                    << "mechanicalProperties has no cell zone of that name."
                    << nl
                    << "With more than one material, each must occupy the "
                    << "cell zone it is named after."
                    << exit(FatalError);
            }

            const labelList& zoneCells = mesh.cellZones()[zoneI];

            forAll(zoneCells, i)
            {
                const label cellI = zoneCells[i];

                if (materialID[cellI] != -1)
                {
                    FatalErrorInFunction
                        << "Cell " << cellI << " is in more than one material "
                        << "zone: " << materialID[cellI] << " and " << lawI
                        << exit(FatalError);
                }

                materialID[cellI] = lawI;
            }
        }

        return materialID;
    }


    tmp<Field<bool>> cellZoneInterface(const fvMesh& mesh, const bool debug)
    {
        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        tmp<Field<bool>> tinterface(new Field<bool>(owner.size(), false));
        Field<bool>& interface = tmpRef(tinterface);

        const labelList materialID(cellMaterialID(mesh));

        label nInterfaceFaces = 0;

        forAll(interface, faceI)
        {
            const label ownID = materialID[owner[faceI]];
            const label neiID = materialID[neighbour[faceI]];

            // -1 on both sides is the single-material case, where nothing is
            // an interface
            if (ownID != neiID)
            {
                interface[faceI] = true;
                nInterfaceFaces++;
            }
        }

        if (debug && nInterfaceFaces > 0)
        {
            InfoInFunction
                << nl << "There are " << nInterfaceFaces
                << " faces on a material interface" << nl << endl;
        }

        return tinterface;
    }


    List<boolList> cellZoneInterfaceCoupled(const fvMesh& mesh)
    {
        List<boolList> interface(mesh.boundary().size());

        forAll(mesh.boundary(), patchI)
        {
            interface[patchI].setSize(mesh.boundary()[patchI].size(), false);
        }

        const labelList materialID(cellMaterialID(mesh));

        // Nothing to do for a single material: every cell is -1, so no face
        // can lie between two different materials
        if (min(materialID) == -1 && max(materialID) == -1)
        {
            return interface;
        }

        // The material of the cell on the other side of each coupled face.
        // Carried across as a field so that the exchange is the one the mesh
        // already knows how to do, rather than a hand-rolled swap
        volScalarField materialIDField
        (
            IOobject
            (
                "materialID",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("minusOne", dimless, -1.0),
            "zeroGradient"
        );

#ifdef OPENFOAM_NOT_EXTEND
        scalarField& materialIDI = materialIDField.primitiveFieldRef();
#else
        scalarField& materialIDI = materialIDField.internalField();
#endif

        forAll(materialIDI, cellI)
        {
            materialIDI[cellI] = materialID[cellI];
        }

        materialIDField.correctBoundaryConditions();

        forAll(mesh.boundary(), patchI)
        {
            if (!mesh.boundary()[patchI].coupled())
            {
                continue;
            }

            const labelUList& faceCells =
                mesh.boundary()[patchI].faceCells();

            const scalarField nbrMaterialID
            (
                materialIDField.boundaryField()[patchI].patchNeighbourField()
            );

            forAll(nbrMaterialID, faceI)
            {
                const label ownID = materialID[faceCells[faceI]];
                const label neiID = label(nbrMaterialID[faceI]);

                if (ownID != neiID)
                {
                    interface[patchI][faceI] = true;
                }
            }
        }

        return interface;
    }

} // End namespace Foam


// ************************************************************************* //
