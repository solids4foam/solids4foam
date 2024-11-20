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

#include "higherOrderGrad.H"



// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


template<class Type>
void Foam::higherOrderGrad::requestGlobalStencilData
(
    const Field<Type>& localField,
    Map<Type>& globalField
) const
{
    globalField.clear();

    if (Pstream::parRun())
    {
        const labelListList& globalCellStencils = this->globalCellStencils();

        // We will send neiGlobalCellI to the processor that contains
        // neiGlobalCellI. That processor will then send back the localField
        // value for that neiGlobalCellI
        // Is it important to initialise the size of the labelLists?
        Map<labelList> requestedData;

        forAll(globalCellStencils, localOriginCellI)
        {
            const labelList& curStencil = globalCellStencils[localOriginCellI];

            forAll(curStencil, cI)
            {
                const label neiGlobalCellID = curStencil[cI];

                if (!globalCells_.isLocal(neiGlobalCellID))
                {
                    // neiGlobalCellID is not on this processor so we need to
                    // request its field data

                    // Determine which processor owns this cell
                    const label procID =
                        globalCells_.whichProcID(neiGlobalCellID);

                    // Request data for this cell
                    requestedData(procID).append(neiGlobalCellID);
                }
            }
        }

        // Exchange requestedData between processors

        // Prepare data to send to neighboring processors
        Map<labelList> toSend(Pstream::nProcs());
        Map<labelList> toReceive(Pstream::nProcs());

        // Populate toSend lists
        forAllIter(Map<labelList>, requestedData, iter)
        {
            const label procNo = iter.key();
            labelList& sendData = iter();

            toSend(procNo).transfer(sendData);
        }

        // Exchange data with neighboring processors
        Pstream::exchange<labelList, label>
        (
            toSend, toReceive
        );

        // Clear requestedData as it is no longer needed
        requestedData.clear();

        // Process received data
        // This is the requested from other processors for localField data on
        // this procesor

        // Create map to hold to requested field data
        // We will record Tuple2(globalCellI, fieldValue)
        typedef Map<List<Tuple2<label, Type>>> FieldDataMapType;
        FieldDataMapType requestedFieldData;

        forAllConstIter(Map<labelList>, toReceive, iter)
        {
            const label procI = iter.key();
            const labelList& receivedData = iter();

            forAll(receivedData, idx)
            {
                const label globalCellI = receivedData[idx];

                if (!globalCells_.isLocal(globalCellI))
                {
                    FatalErrorInFunction
                        << "Global cell " << globalCellI
                        << " is not on this proc!" << abort(FatalError);
                }

                // Get local ID
                const label localCellI = globalCells_.toLocal(globalCellI);

                // Record local field value
                requestedFieldData(procI).append
                (
                    Tuple2<label, Type>(globalCellI, localField[localCellI])
                );
            }
        }

        // Send the request field data back to the processors who requested it

        // Prepare data to send to neighboring processors
        FieldDataMapType toSendField(Pstream::nProcs());
        FieldDataMapType toReceiveField(Pstream::nProcs());

        // Populate toSend lists
        forAllIter(typename FieldDataMapType, requestedFieldData, iter)
        {
            const label procNo = iter.key();
            List<Tuple2<label, Type>>& sendData = iter();

            toSendField(procNo).transfer(sendData);
        }

        // Exchange data with neighboring processors
        Pstream::exchange<List<Tuple2<label, Type>>, Tuple2<label, Type>>
        (
            toSendField, toReceiveField
        );

        // Clear requestedFieldData as it is no longer needed
        requestedFieldData.clear();

        // Populate the globalField map with the toReceiveField data
        forAllConstIter(typename FieldDataMapType, toReceiveField, iter)
        {
            //const label procI = iter.key();
            const List<Tuple2<label, Type>>& receivedData = iter();

            forAll(receivedData, cI)
            {
                const label globalCellI = receivedData[cI].first();
                const Type& fieldValue = receivedData[cI].second();

                // Record data in the global field map
                globalField(globalCellI) = fieldValue;
            }
        }
    }
}


// ************************************************************************* //
