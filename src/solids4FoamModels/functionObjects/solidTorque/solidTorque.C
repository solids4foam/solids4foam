/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "solidTorque.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidTorque, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidTorque,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidTorque::writeData()
{
    if (patchFound_)
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

        // Check if the stress tensor field is found
        if (mesh.foundObject<volSymmTensorField>(stressName_))
        {
            const symmTensorField& sigma =
                mesh.lookupObject<volSymmTensorField>
                (
                    stressName_
                ).boundaryField()[historyPatchID_];

            // Calculate moment arms
            const vectorField patchC =
                mesh.C().boundaryField()[historyPatchID_];

            vectorField r = patchC - pointOnAxis_;
            r -= axis_*(axis_ & r);

            // Calculate torque
            const vectorField force =
                mesh.Sf().boundaryField()[historyPatchID_] & sigma;

            const vector torque = gSum(r ^ force);

            if (Pstream::master())
            {
                historyFilePtr_()
                    << time_.time().value()
                    << " " << torque.x() << " " << torque.y()
                    << " " << torque.z();

                historyFilePtr_() << endl;
            }
        }
        else
        {
            InfoIn(this->name() + " function object constructor")
                << stressName_ << " not found" << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidTorque::solidTorque
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    historyPatchID_(-1),
    patchFound_(false),
    stressName_
    (
        dict.found("stressName")
      ? word(dict.lookup("stressName"))
      : word("sigma")
    ),
    pointOnAxis_(dict.lookup("pointOnAxis")),        
    axis_(dict.lookup("axisDirection")),
    historyFilePtr_(NULL)
{
    Info<< "Creating " << this->name() << " function object" << endl;

    // Normalise the axis
    const scalar magAxis = mag(axis_);
    if (magAxis < SMALL)
    {
        FatalErrorIn("Foam::solidTorque::solidTorque(...)")
            << "Invalid axis!" << abort(FatalError);
    }
    axis_ /= magAxis;

    // Check the patch is found
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

    const word historyPatchName = word(dict.lookup("historyPatch"));
    historyPatchID_ = mesh.boundaryMesh().findPatchID(historyPatchName);
    if (historyPatchID_ == -1)
    {
        WarningIn(this->name() + " function object constructor")
            << "history patch " << historyPatchName << " not found"
            << endl;
    }
    else
    {
        patchFound_ = true;
    }

    // Create history file if not already created
    if (historyFilePtr_.empty() && patchFound_)
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
            // Use the function object name in the file name to allow multiple
            // objects defined on the same patch
            historyFilePtr_.reset
            (
                new OFstream
                (
                    historyDir/"solidTorque" + historyPatchName + name + ".dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time" << " "
                    << "torqueX" << " " << "torqueY" << " " << "torqueZ"
                    << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidTorque::start()
{
    return writeData();
}

#if FOAMEXTEND > 40
bool Foam::solidTorque::execute(const bool forceWrite)
#else
bool Foam::solidTorque::execute()
#endif
{
    return writeData();
}


bool Foam::solidTorque::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
