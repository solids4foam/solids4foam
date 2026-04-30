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

\*----------------------------------------------------------------------------*/

#include "solidPressureMinMax.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidPressureMinMax, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidPressureMinMax,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidPressureMinMax::writeData()
{
    // Lookup the solid mesh
    const fvMesh* meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh& mesh = *meshPtr;

    if (mesh.foundObject<volSymmTensorField>("sigma"))
    {
        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>("sigma");

        // tr(sigma) = sigmaXX + sigmaYY + sigmaZZ
        const scalarField trSigma(tr(sigma.internalField()));

        const scalar trMin = gMin(trSigma);
        const scalar trMax = gMax(trSigma);

        // Hydrostatic pressure p = -tr(sigma)/3
        const scalar pMin = -trMax / 3.0;
        const scalar pMax = -trMin / 3.0;

        if (Pstream::master())
        {
            historyFilePtr_()
                << time_.time().value()
                << " " << pMin
                << " " << pMax
                << " " << trMin
                << " " << trMax
                << endl;
        }
    }
    else
    {
        InfoIn(this->name() + " function object writeData()")
            << "volSymmTensorField sigma not found" << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidPressureMinMax::solidPressureMinMax
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object" << endl;

    if (Pstream::master())
    {
        fileName historyDir;

        const fvMesh* meshPtr = NULL;
        if (time_.foundObject<fvMesh>("solid"))
        {
            meshPtr = &(time_.lookupObject<fvMesh>("solid"));
        }
        else if (time_.foundObject<fvMesh>("region0"))
        {
            meshPtr = &(time_.lookupObject<fvMesh>("region0"));
        }

        const word startTimeName =
            (meshPtr != NULL)
          ? time_.timeName(meshPtr->time().startTime().value())
          : "0";

        if (Pstream::parRun())
        {
            historyDir =
                time_.path()/".."/"postProcessing"/startTimeName;
        }
        else
        {
            historyDir = time_.path()/"postProcessing"/startTimeName;
        }

        mkDir(historyDir);

        historyFilePtr_.reset
        (
            new OFstream(historyDir/"solidPressureMinMax.dat")
        );

        if (historyFilePtr_.valid())
        {
            historyFilePtr_()
                << "# Time  p_min  p_max  trSigma_min  trSigma_max" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidPressureMinMax::start()
{
    return false;
}


#if FOAMEXTEND
bool Foam::solidPressureMinMax::execute(const bool forceWrite)
#else
bool Foam::solidPressureMinMax::execute()
#endif
{
    return writeData();
}


bool Foam::solidPressureMinMax::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::solidPressureMinMax::write()
{
    return false;
}
#endif

// ************************************************************************* //
