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

Application
    Test-fluxCorrectedVelocityRestart

Description
    Checks that fluxCorrectedVelocity patch values in U_0 are restored when a
    backward-time-scheme case is restarted.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "physicsModel.H"
#include "fluxCorrectedVelocityFvPatchVectorField.H"
#include "compatibilityFunctions.H"

using namespace Foam;

namespace
{
    // Minimal concrete regIOobject for reading a field as a dictionary through
    // the active file handler. localIOdictionary is unavailable in foam-extend.
    class fieldFileReader
    :
        public regIOobject
    {
    public:

        fieldFileReader(const IOobject& io)
        :
            regIOobject(io)
        {}

        virtual bool writeData(Ostream&) const
        {
            return false;
        }
    };
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    IOobject oldUIO
    (
        "U_0",
        runTime.timeName(),
        "fluid",
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

#ifdef OPENFOAM_NOT_EXTEND
    if (!oldUIO.typeHeaderOk<volVectorField>(false))
#else
    if (!oldUIO.headerOk())
#endif
    {
        FatalErrorInFunction
            << "Cannot read " << oldUIO.objectPath() << abort(FatalError);
    }

    fieldFileReader oldUReader(oldUIO);
    const dictionary oldUDict
    (
        oldUReader.readStream(volVectorField::typeName)
    );
    oldUReader.close();

    const dictionary& oldUBoundaryDict = oldUDict.subDict("boundaryField");

    autoPtr<physicsModel> physics = physicsModel::New(runTime);
    const fvMesh& mesh = runTime.lookupObject<fvMesh>("fluid");
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>("phi");

    if (U.nOldTimes() < 1)
    {
        FatalErrorInFunction
            << U.name() << " has no old-time level" << abort(FatalError);
    }

    if (phi.dimensions() != dimVolume/dimTime)
    {
        FatalErrorInFunction
            << phi.name() << " is not a volumetric flux" << abort(FatalError);
    }

    label nChecked = 0;

    forAll(U.boundaryField(), patchI)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchI];

        if (!isA<fluxCorrectedVelocityFvPatchVectorField>(Up))
        {
            continue;
        }

        const word& patchName = Up.patch().name();
        const dictionary& oldUPatchDict =
            oldUBoundaryDict.subDict(patchName);
        const vectorField writtenOldU
        (
            "value",
            oldUPatchDict,
            Up.size()
        );
        const vectorField& restoredOldU =
            U.oldTime().boundaryField()[patchI];

        const scalar oldError = gMax(mag(restoredOldU - writtenOldU));

        const fvPatch& patch = Up.patch();
        const vectorField n(patch.nf());
        const vectorField UzeroGrad(Up.patchInternalField());
        const vectorField expectedCurrent
        (
            UzeroGrad
          - n*(n & UzeroGrad)
          + n*phi.boundaryField()[patchI]/patch.magSf()
        );
        const scalar currentError = gMax(mag(Up - expectedCurrent));

        Info<< patchName << ": old-time error = " << oldError
            << ", current-time error = " << currentError << endl;

        if (oldError > SMALL || currentError > SMALL)
        {
            FatalErrorInFunction
                << "fluxCorrectedVelocity restart mismatch on patch "
                << patchName << abort(FatalError);
        }

        nChecked++;
    }

    if (!nChecked)
    {
        FatalErrorInFunction
            << "No fluxCorrectedVelocity patches found" << abort(FatalError);
    }

    Info<< "Checked " << nChecked << " fluxCorrectedVelocity patches" << endl;

    return 0;
}


// ************************************************************************* //
