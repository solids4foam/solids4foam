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

Class
    Foam::functionObjects::beamEnergyData

Description
    FunctionObject reports BEAM internal and kinetic energy data

Author
    Seevani Bali, UCD.  All rights reserved

SourceFiles
    beamEnergyData.C

\*---------------------------------------------------------------------------*/

#ifndef beamEnergyData_H
#define beamEnergyData_H

#include "fvMeshFunctionObject.H"
#include "OFstream.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class beamEnergyData Declaration
\*---------------------------------------------------------------------------*/

class beamEnergyData
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Forces/moment file ptr
        autoPtr<OFstream> historyFilePtr_;

        //- Per-beam energy file ptrs
        PtrList<OFstream> beamHistoryFilePtrs_;

    // Private Member Functions

        //- Write data
        bool writeData();

public:

    //- Runtime type information
    TypeName("beamEnergyData");


    // Constructors

        //- Construct from Time and dictionary
        beamEnergyData
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        beamEnergyData(const beamEnergyData&) = delete;

        //- No copy assignment
        void operator=(const beamEnergyData&) = delete;


    //- Destructor
    virtual ~beamEnergyData() = default;


    // Member Functions

        //- Read the beamEnergyData data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the beamEnergyData
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
