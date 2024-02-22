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

#include "solidTorque.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "lookupSolidModel.H"
#include "OSspecific.H"

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
            // Cauchy stress tensor
            const symmTensorField& sigma =
                mesh.lookupObject<volSymmTensorField>
                (
                    stressName_
                ).boundaryField()[historyPatchID_];

            // Patch area vectors
            const vectorField& patchSf =
                mesh.Sf().boundaryField()[historyPatchID_];

            // Lookup solid model
            const solidModel& solMod = lookupSolidModel(mesh);

            scalar torque = 0.0;

            // Check if it is a linear or nonlinear geometry case
            if
            (
                solMod.nonLinGeom() == nonLinearGeometry::LINEAR_GEOMETRY
             || solMod.nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN
            )
            {
                // Calculate moment arms
                const vectorField patchC
                (
                    mesh.C().boundaryField()[historyPatchID_]
                );
                vectorField r(patchC - pointOnAxis_);
                r -= axis_*(axis_ & r);

                // Calculate force
                const vectorField force(patchSf & sigma);

                // Calculate torque
                torque = gSum(axis_ & (r ^ force));
            }
            else if (solMod.nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
            {
                // Lookup the total displacement field
                const vectorField& D =
                    mesh.lookupObject<volVectorField>
                    (
                        "D"
                    ).boundaryField()[historyPatchID_];

                // Calculate deformed patch face centres
                const vectorField patchDeformC
                (
                    mesh.C().boundaryField()[historyPatchID_] + D
                );

                // Calculate moment arms
                vectorField r(patchDeformC - pointOnAxis_);
                r -= axis_*(axis_ & r);

                // Lookup the inverse of the deformation gradient
                const tensorField& Finv =
                    mesh.lookupObject<volTensorField>
                    (
                        "Finv"
                    ).boundaryField()[historyPatchID_];

                // Lookup the Jacobian
                const scalarField& J =
                    mesh.lookupObject<volScalarField>
                    (
                        "J"
                    ).boundaryField()[historyPatchID_];

                // Calculate area vectors in the deformed configuration
                const vectorField patchDeformSf(J*Finv.T() & patchSf);

                // Calculate force
                const vectorField force(patchDeformSf & sigma);

                // Calculate torque
                torque = gSum(axis_ & (r ^ force));
            }
            else
            {
                FatalErrorIn("bool Foam::solidForces::writeData()")
                    << "Unknown solidModel nonLinGeom type = "
                    << solMod.nonLinGeom() << abort(FatalError);
            }

            // Write to file
            if (Pstream::master())
            {
                historyFilePtr_()
                    << time_.time().value() << " " << torque << endl;
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
    historyFilePtr_()
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
                historyDir = time_.path()/".."/"postProcessing"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"postProcessing"/startTimeName;
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
                    << "# Time torque" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidTorque::start()
{
    return writeData();
}

#if FOAMEXTEND
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

#ifdef OPENFOAM_NOT_EXTEND
bool Foam::solidTorque::write()
{
    return false;
}
#endif

// ************************************************************************* //
