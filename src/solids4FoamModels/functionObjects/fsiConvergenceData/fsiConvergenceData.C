/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "fsiConvergenceData.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "boundBox.H"
#include "fluidSolidInterface.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "IOmanip.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fsiConvergenceData, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fsiConvergenceData,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fsiConvergenceData::writeData()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);
    
    const fluidSolidInterface& fsi =
        mesh.thisDb().parent().lookupObject<fluidSolidInterface>("fsiProperties");
        // mesh.parent().lookupObject<fluidSolidInterface>("fsiProperties");
    
    if (Pstream::master())
    {
        historyFilePtr_() << time_.time().value() << tab
                          << fsi.outerCorr() << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fsiConvergenceData::fsiConvergenceData
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object." << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(mesh.time().startTime().value());


            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"history"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"history"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);



            // Open new file at start up

            OStringStream FileName;
            FileName() << "fsiConvergenceData.dat";

            historyFilePtr_.reset
            (
                new OFstream(historyDir/word(FileName.str()))
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time" << tab << "nFsiCorrectors" << endl;

            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fsiConvergenceData::start()
{
    return false;
}


#if FOAMEXTEND
bool Foam::fsiConvergenceData::execute(const bool forceWrite)
#else
bool Foam::fsiConvergenceData::execute()
#endif
{
    return writeData();
}


bool Foam::fsiConvergenceData::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}


#ifdef OPENFOAMESIORFOUNDATION
bool Foam::fsiConvergenceData::write()
{
    return writeData();
}
#endif

// ************************************************************************* //
