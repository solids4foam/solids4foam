/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "beamModel.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::beamModel::beamPointData
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf,
    const label bI
) const
{
    const fvMesh& mesh = this->mesh();

    label nCellZones = mesh.cellZones().size();

    if (nCellZones < 2)
    {
        label nPoints = mesh.nCells() + 1;

        tmp<Field<Type> > tPointData
        (
            new Field<Type>(nPoints, pTraits<Type>::zero)
        );
        Field<Type>& pointData = tPointData.ref();

        const Field<Type>& sfI = sf.internalField();

        if (sf.boundaryField()[startPatchIndex()].size())
        {
            pointData[0] = sf.boundaryField()[startPatchIndex()][0];
        }
        else
        {
            // Find left processor patch
            forAll(mesh.boundary(), patchI)
            {
                if (isA<processorFvPatch>(mesh.boundary()[patchI]))
                {
                    const labelList& fc = mesh.boundary()[patchI].faceCells();

                    if (fc[0] == 0)
                    {
                        pointData[0] =
                            sf.boundaryField()[patchI][0];
                    }
                }
            }
        }

        if (sf.boundaryField()[endPatchIndex()].size())
        {
            pointData[nPoints-1] = sf.boundaryField()[endPatchIndex()][0];
        }
        else
        {
            // Find right processor patch
            forAll(mesh.boundary(), patchI)
            {
                if (isA<processorFvPatch>(mesh.boundary()[patchI]))
                {
                    const labelList& fc = mesh.boundary()[patchI].faceCells();

                    if (fc[0] == (mesh.nCells()-1))
                    {
                        pointData[nPoints-1] =
                            sf.boundaryField()[patchI][0];
                    }
                }
            }
        }

        for (label i=0; i<sfI.size(); i++)
        {
            pointData[i+1] = sfI[i];
        }

        return tPointData;
    }
    else
    {
        const cellZone& cz = mesh.cellZones()[bI];

        if (cz.size() == 0)
        {
            return tmp<Field<Type> >
            (
                new Field<Type>(0, pTraits<Type>::zero)
            );
        }

        label nPoints = cz.size() + 1;

        tmp<Field<Type> > tPointData
        (
            new Field<Type>(nPoints, pTraits<Type>::zero)
        );
        Field<Type>& pointData = tPointData.ref();

        const Field<Type>& sfI = sf.internalField();

        const labelList startPatchCells =
            mesh.boundary()[startPatchIndex(bI)].faceCells();

        forAll(sf.boundaryField()[startPatchIndex(bI)], faceI)
        {
            if (cz.whichCell(startPatchCells[faceI]) != -1)
            {
                pointData[0] =
                    sf.boundaryField()[startPatchIndex(bI)][faceI];
            }
        }

        // Find left and right processor patch
        forAll(mesh.boundary(), patchI)
        {
            if (isA<processorFvPatch>(mesh.boundary()[patchI]))
            {
                const labelList& fc =
                    mesh.boundary()[patchI].faceCells();

                label firstCell = min(cz);
                label lastCell = max(cz);

                forAll(fc, cI)
                {
                    if (cz.whichCell(fc[cI]) != -1)
                    {
                        if (fc[cI] == firstCell)
                        {
                            pointData[0] =
                                sf.boundaryField()[patchI][cI];
                        }
                        else if (fc[cI] == lastCell)
                        {
                            label lastPointIndex =
                                pointData.size()-1;

                            pointData[lastPointIndex] =
                                sf.boundaryField()[patchI][cI];
                        }
                    }
                }
            }
        }

        const labelList endPatchCells =
            mesh.boundary()[endPatchIndex(bI)].faceCells();

        forAll(sf.boundaryField()[endPatchIndex(bI)], faceI)
        {
            if (cz.whichCell(endPatchCells[faceI]) != -1)
            {
                pointData[nPoints-1] =
                    sf.boundaryField()[endPatchIndex(bI)][faceI];
            }
        }

        const labelList& nei = mesh.neighbour();

        label index = 1;
        for (label i=0; i<sfI.size(); i++)
        {
            if (cz.whichCell(nei[i]) != -1)
            {
                pointData[index] = sfI[i];
                index++;
            }
        }

        return tPointData;
    }
}

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::beamModel::beamPointData
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf,
    const label bI,
    const label segI
) const
{
    const fvMesh& mesh = this->mesh();

    label nCellZones = mesh.cellZones().size();

    if (nCellZones < 2)
    {
        label nPoints = mesh.nCells() + 1;

        tmp<Field<Type> > tPointData
        (
            new Field<Type>(2, pTraits<Type>::zero)
        );
        Field<Type>& pointData = tPointData();

        const Field<Type>& sfI = sf.internalField();

	if (segI == 0)
	{
	    pointData[0] = sf.boundaryField()[startPatchIndex()][0];
	    pointData[1] = sfI[0];
	}
	else if (segI == (nPoints-1))
	{
	    pointData[0] = sfI[sfI.size()-1];
	    pointData[1] = sf.boundaryField()[endPatchIndex()][0];
	}
	else
	{
	    pointData[0] = sfI[segI-1];
	    pointData[1] = sfI[segI];
	}

        return tPointData;
    }
    else
    {
        const cellZone& cz = mesh.cellZones()[bI];

        label nPoints = cz.size() + 1;

        tmp<Field<Type> > tPointData
        (
            new Field<Type>(2, pTraits<Type>::zero)
        );
        Field<Type>& pointData = tPointData();

        const Field<Type>& sfI = sf.internalField();

	label start = 0;
	label i = 0;
	while (i < bI)
	{
	    start += mesh.cellZones()[i].size();
	    i++;
	}
	label glSegI = start + segI;

	if (segI == 0)
	{
            const labelList startPatchCells =
                mesh.boundary()[startPatchIndex(bI)].faceCells();

	    forAll(sf.boundaryField()[startPatchIndex(bI)], faceI)
	    {
                if (cz.whichCell(startPatchCells[faceI]) != -1)
		{
                    pointData[0] =
                        sf.boundaryField()[startPatchIndex(bI)][faceI];

		    label glIntFace = glSegI - bI;
		    pointData[1] = sfI[glIntFace];
		}
	    }
	}
	else if (segI == (nPoints-2))
	{
            const labelList endPatchCells =
                mesh.boundary()[endPatchIndex(bI)].faceCells();

	    forAll(sf.boundaryField()[endPatchIndex(bI)], faceI)
	    {
                if (cz.whichCell(endPatchCells[faceI]) != -1)
		{
		    label glIntFace = glSegI - bI;
		    pointData[0] = sfI[glIntFace-1];

                    pointData[1] =
                        sf.boundaryField()[endPatchIndex(bI)][faceI];
		}
	    }
	}
	else
	{
	    label glIntFace = glSegI - bI;

	    pointData[0] = sfI[glIntFace-1];
	    pointData[1] = sfI[glIntFace];
	}

        return tPointData;
    }
}

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::beamModel::beamSegmentData
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const label bI,
    const label segI
) const
{
    const fvMesh& mesh = this->mesh();

    label nCellZones = mesh.cellZones().size();

    if (nCellZones < 2)
    {
        // label nPoints = mesh.nCells() + 1;

        tmp<Field<Type> > tSegmentData
        (
            new Field<Type>(1, pTraits<Type>::zero)
        );
        Field<Type>& segmentData = tSegmentData();

        const Field<Type>& vfI = vf.internalField();

        segmentData[0] = vfI[segI];

        return tSegmentData;
    }
    else
    {
        // const cellZone& cz = mesh.cellZones()[bI];

        // label nPoints = cz.size() + 1;

        tmp<Field<Type> > tSegmentData
        (
            new Field<Type>(1, pTraits<Type>::zero)
        );
        Field<Type>& segmentData = tSegmentData();

        const Field<Type>& vfI = vf.internalField();

	label start = 0;
	label i = 0;
	while (i < bI)
	{
	    start += mesh.cellZones()[i].size();
	    i++;
	}
	label glSegI = start + segI;

        segmentData[0] = vfI[glSegI];

        return tSegmentData;
    }
}
// ************************************************************************* //
