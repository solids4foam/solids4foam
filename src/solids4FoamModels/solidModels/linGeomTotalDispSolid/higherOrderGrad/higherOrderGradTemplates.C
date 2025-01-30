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
        // typedef Map<List<Tuple2<label, Type>>> FieldDataMapType;
        // FieldDataMapType requestedFieldData;
        Map<labelList> requestedCellIDs;
        Map<List<Type>> requestedCellFieldValues;

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

                // Record the cell ID
                requestedCellIDs(procI).append(globalCellI);

                // Record local field value
                requestedCellFieldValues(procI).append(localField[localCellI]);
            }
        }

        // Send the request field data back to the processors who requested it

        // Prepare data to send to neighboring processors
        // We must send the labels and the values seperately as the exchange
        // function only works with contigious data types
        Map<labelList> toSendLabels(Pstream::nProcs());
        Map<List<Type>> toSendField(Pstream::nProcs());
        Map<labelList> toReceiveLabels(Pstream::nProcs());
        Map<List<Type>> toReceiveField(Pstream::nProcs());

        // Populate toSend lists
        forAllIter(Map<labelList>, requestedCellIDs, iter)
        {
            const label procNo = iter.key();
            labelList& sendData = iter();

            toSendLabels(procNo).transfer(sendData);
        }
        forAllIter(typename Map<List<Type>>, requestedCellFieldValues, iter)
        {
            const label procNo = iter.key();
            List<Type>& sendData = iter();

            toSendField(procNo).transfer(sendData);
        }

        // Exchange data with neighboring processors
        Pstream::exchange<labelList, label>
        (
            toSendLabels, toReceiveLabels
        );
        Pstream::exchange<List<Type>, Type>
        (
            toSendField, toReceiveField
        );

        // Clear requested data as it is no longer needed
        requestedCellIDs.clear();
        requestedCellFieldValues.clear();

        // Populate the globalField map with the toReceiveField data
        forAllConstIter(Map<labelList>, toReceiveLabels, iter)
        {
            const label procI = iter.key();
            const labelList& receivedLabels = iter();
            const List<Type>& receivedField = toReceiveField[procI];

            forAll(receivedLabels, cI)
            {
                const label globalCellI = receivedLabels[cI];
                const Type& fieldValue = receivedField[cI];

                // Record data in the global field map
                globalField(globalCellI) = fieldValue;
            }
        }
    }
}


// ************************************************************************* //
