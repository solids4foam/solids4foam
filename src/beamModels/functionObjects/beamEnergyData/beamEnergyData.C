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

#include "beamEnergyData.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(beamEnergyData, 0);
    addToRunTimeSelectionTable(functionObject, beamEnergyData, dictionary);
}
}

namespace
{
    static const Foam::label energyWidth = 16;

    void writeEnergyHeader(Foam::OFstream& os)
    {
        os  << Foam::setw(energyWidth) << "Time"
            << Foam::setw(energyWidth) << "E_{int}"
            << Foam::setw(energyWidth) << "E_{kin}"
            << Foam::setw(energyWidth) << "E_{kin-lin}"
            << Foam::setw(energyWidth) << "E_{kin-ang}"
            << Foam::setw(energyWidth) << "E_{tot}"
            << Foam::endl;
    }

    void writeEnergyRow
    (
        Foam::OFstream& os,
        const Foam::scalar time,
        const Foam::scalar intE,
        const Foam::scalar kinE,
        const Foam::scalar kinE_lin,
        const Foam::scalar kinE_ang,
        const Foam::scalar totE
    )
    {
        os  << Foam::setw(energyWidth) << time
            << Foam::setw(energyWidth) << intE
            << Foam::setw(energyWidth) << kinE
            << Foam::setw(energyWidth) << kinE_lin
            << Foam::setw(energyWidth) << kinE_ang
            << Foam::setw(energyWidth) << totE
            << Foam::endl;
    }
}

// * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

bool Foam::functionObjects::beamEnergyData::writeData()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const beamModel& beam =
        mesh.objectRegistry::parent().lookupObject<beamModel>("beamProperties");

    // The internal strains and curvatures for calculation of internal energy
    const surfaceVectorField& Gamma =
        mesh.lookupObject<surfaceVectorField>("Gamma");
    const surfaceVectorField& K = mesh.lookupObject<surfaceVectorField>("K");
    const surfaceTensorField& CQ = mesh.lookupObject<surfaceTensorField>("CQ");
    const surfaceTensorField& CM = mesh.lookupObject<surfaceTensorField>("CM");

    // Properties of the beam like length, area, second moment of area, density
    const volScalarField& L = beam.L();

    // Beam linear and angular velocity for calculation of kinetic energy
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    const volVectorField& Omega = mesh.lookupObject<volVectorField>("Omega");

    // Accessing internal fields
    const vectorField& GammaI = Gamma.internalField();
    const tensorField& CQI = CQ.internalField();
    const vectorField& KI = K.internalField();
    const tensorField& CMI = CM.internalField();
    const vectorField& UI = U.internalField();
    const vectorField& OmegaI = Omega.internalField();
    const labelList& own = mesh.owner();

    // Second moment of area scaling factor-if not specified in beamProperties
    // then the value is equal to 1.0
    const scalar kCI = beam.kCI();

    // Energy parameters initialised
    scalarField intE(beam.nBeams(), 0);
    scalarField kinE_lin(beam.nBeams(), 0);
    scalarField kinE_ang(beam.nBeams(), 0);

    // Calculation of internal energy at internal faces
    forAll(GammaI, faceI)
    {
        const label cellI = own[faceI];
        const label bI = beam.whichBeam(beam.globalCellIndex(cellI));

        intE[bI] += 0.5*L[cellI]*
        (
            (GammaI[faceI] & (CQI[faceI] & GammaI[faceI]))
          + (KI[faceI] & (CMI[faceI] & KI[faceI]))
        );
    }

    // Calculation of kinetic energy at cell centres
    // using cell-centre velocity fields
    forAll(UI, cellI)
    {
        const label bI = beam.whichBeam(beam.globalCellIndex(cellI));

        const scalar rho = beam.rho(bI).value();
        const scalar A = beam.A(bI).value();

        const tensor CI
        (
            kCI*beam.J(bI).value(), 0, 0,
            0, kCI*beam.Iyy(bI).value(), 0,
            0, 0, kCI*beam.Izz(bI).value()
        );

        kinE_lin[bI] += 0.5*rho*L[cellI]*
        (
            A*(UI[cellI] & UI[cellI])
        );
        kinE_ang[bI] += 0.5*rho*L[cellI]*
        (
            (OmegaI[cellI] & (CI & OmegaI[cellI]))
        );
    }

    // Adding contribution of boundary fields to the internal energy
    forAll(Gamma.boundaryField(),patchI)
    {
        const vectorField& pGamma = Gamma.boundaryField()[patchI];
        const tensorField& pCQ = CQ.boundaryField()[patchI];
        const vectorField& pK = K.boundaryField()[patchI];
        const tensorField& pCM = CM.boundaryField()[patchI];
        const labelList& faceCells = mesh.boundary()[patchI].faceCells();

        forAll(pGamma,faceI)
        {
            const label cellI = faceCells[faceI];
            const label bI = beam.whichBeam(beam.globalCellIndex(cellI));

            intE[bI] += 0.25*L[cellI]*
            (
                (pGamma[faceI] & (pCQ[faceI] & pGamma[faceI]))
                + (pK[faceI] & (pCM[faceI] & pK[faceI]))
            );
        }
    }

    if (Pstream::master())
    {
        scalar totalIntE(0);
        scalar totalKinE(0);
        scalar totalKinE_lin(0);
        scalar totalKinE_ang(0);
        scalar totalE(0);

        forAll(intE, beamI)
        {
            const scalar beamKinE = kinE_lin[beamI] + kinE_ang[beamI];
            const scalar beamTotE = intE[beamI] + beamKinE;

            totalIntE += intE[beamI];
            totalKinE += beamKinE;
            totalKinE_lin += kinE_lin[beamI];
            totalKinE_ang += kinE_ang[beamI];
            totalE += beamTotE;

            if (beamHistoryFilePtrs_.set(beamI))
            {
                writeEnergyRow
                (
                    beamHistoryFilePtrs_[beamI],
                    time_.time().value(),
                    intE[beamI],
                    beamKinE,
                    kinE_lin[beamI],
                    kinE_ang[beamI],
                    beamTotE
                );
            }
        }

        writeEnergyRow
        (
            historyFilePtr_(),
            time_.time().value(),
            totalIntE,
            totalKinE,
            totalKinE_lin,
            totalKinE_ang,
            totalE
        );
    }

    return true;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::beamEnergyData::beamEnergyData
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    name_(name),
    time_(runTime),
    regionName_(polyMesh::defaultRegion),
    historyFilePtr_(),
    beamHistoryFilePtrs_()
{
    Info<< "Creating " << this->name() << " function object!" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const beamModel& beam =
        mesh.objectRegistry::parent().lookupObject<beamModel>("beamProperties");

    // Create history file if not already created
    if (!historyFilePtr_.valid())
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
                historyDir = time_.path()/".."/"postProcessing"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"postProcessing"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            OStringStream FileName;
            FileName() << "beamEnergyData.dat";

            historyFilePtr_.reset
            (
                new OFstream(historyDir/word(FileName.str()))
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                writeEnergyHeader(historyFilePtr_());
            }

            beamHistoryFilePtrs_.setSize(beam.nBeams());

            for (label beamI = 0; beamI < beam.nBeams(); ++beamI)
            {
                OStringStream beamFileName;
                beamFileName() << "beamEnergyData_beam" << beamI << ".dat";

                beamHistoryFilePtrs_.set
                (
                    beamI,
                    new OFstream(historyDir/word(beamFileName.str()))
                );

                if (beamHistoryFilePtrs_.set(beamI))
                {
                    writeEnergyHeader(beamHistoryFilePtrs_[beamI]);
                }
            }
        }
    }

    // read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::beamEnergyData::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}


bool Foam::functionObjects::beamEnergyData::execute()
{
    return writeData();
}


bool Foam::functionObjects::beamEnergyData::end()
{
    return true;
}


bool Foam::functionObjects::beamEnergyData::write()
{
    return true;
}


// ************************************************************************* //
