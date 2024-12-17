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

#include "solidTractions.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidTractions, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidTractions,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidTractions::writeData()
{
    if (time_.outputTime())
    {
        Info<< name_ << " functionObject: writing traction field" << nl << endl;

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

        // Lookup the stress field
        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>("sigma");

        if
        (
            mesh.foundObject<volVectorField>("DD")
         && mesh.foundObject<volTensorField>("relF")
        )
        {
            // Updated Lagrangian
            // The mesh has been moved to the deformed configuration

            // Create the traction
            volVectorField traction
            (
                IOobject
                (
                    "traction",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimPressure, vector::zero)
            );

            forAll(traction.boundaryField(), patchI)
            {
                if (!traction.boundaryField()[patchI].coupled())
                {
                    // It is assumed that sigma is the true (Cauchy) stress
#ifdef OPENFOAM_NOT_EXTEND
                    traction.boundaryFieldRef()[patchI] =
#else
                    traction.boundaryField()[patchI] =
#endif
                        mesh.boundary()[patchI].nf()
                      & sigma.boundaryField()[patchI];
                }
            }

            traction.write();
        }
        else if (mesh.foundObject<volTensorField>("F"))
        {
            // Total Lagrangian
            // The mesh is in its initial configuration

            // Lookup the inverse deformation gradient
            const volTensorField& Finv =
                mesh.lookupObject<volTensorField>("Finv");

            // Create the traction
            volVectorField traction
            (
                IOobject
                (
                    "traction",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimPressure, vector::zero)
            );

            forAll(traction.boundaryField(), patchI)
            {
                if (!traction.boundaryField()[patchI].coupled())
                {
                    vectorField nCurrent
                    (
                        Finv.boundaryField()[patchI].T()
                      & mesh.boundary()[patchI].nf()
                    );
                    nCurrent /= mag(nCurrent);

                   // It is assumed that sigma is the true (Cauchy) stress
#ifdef OPENFOAM_NOT_EXTEND
                   traction.boundaryFieldRef()[patchI] =
#else
                   traction.boundaryField()[patchI] =
#endif
                      nCurrent & sigma.boundaryField()[patchI];
                }
            }

            traction.write();
        }
        else
        {
            // Small strain approach

            // Create the traction
            volVectorField traction
            (
                IOobject
                (
                    "traction",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimPressure, vector::zero)
            );

            forAll(traction.boundaryField(), patchI)
            {
                if (!traction.boundaryField()[patchI].coupled())
                {
#ifdef OPENFOAM_NOT_EXTEND
                    traction.boundaryFieldRef()[patchI] =
#else
                    traction.boundaryField()[patchI] =
#endif
                        mesh.boundary()[patchI].nf()
                      & sigma.boundaryField()[patchI];
                }
            }

            traction.write();
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidTractions::solidTractions
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t)
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidTractions::start()
{
    return false;
}


#if FOAMEXTEND
bool Foam::solidTractions::execute(const bool forceWrite)
#else
bool Foam::solidTractions::execute()
#endif
{
    return writeData();
}


bool Foam::solidTractions::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::solidTractions::write()
{
    return false;
}
#endif

// ************************************************************************* //
