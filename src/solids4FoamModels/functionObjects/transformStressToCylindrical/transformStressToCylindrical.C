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

#include "transformStressToCylindrical.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(transformStressToCylindrical, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        transformStressToCylindrical,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::transformStressToCylindrical::writeData()
{
    if (time_.outputTime())
    {
        Info<< name_ << " functionObject: transforming sigma field"
            << nl << endl;

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

        volSymmTensorField sigmaTransformed
        (
            IOobject
            (
                "sigma:Transformed",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
            "calculated"
        );

        forAll(sigmaTransformed, cellI)
        {
            const vector C = mesh.C()[cellI];

            coordinateSystem localCS("localCS", origin_, axis_,  C);

#if OPENFOAM_ORG
            const tensor R = localCS.R().R();
#else
            const tensor R = localCS.R();
#endif

            sigmaTransformed[cellI] = Foam::symm(R.T() & (sigma[cellI] & R));
        }

        forAll(sigmaTransformed.boundaryField(), patchI)
        {

            const symmTensorField& sigmaP = sigma.boundaryField()[patchI];

#ifdef OPENFOAM_NOT_EXTEND
            symmTensorField& sigmaTP =
                sigmaTransformed.boundaryFieldRef()[patchI];
#else
            symmTensorField& sigmaTP =
                sigmaTransformed.boundaryField()[patchI];
#endif

            forAll(sigmaTP, faceI)
            {
                const vector& Cf = mesh.boundary()[patchI].Cf()[faceI];

                coordinateSystem localCS("localCS", origin_, axis_, Cf);

#if OPENFOAM_ORG
            const tensor R = localCS.R().R();
#else
            const tensor R = localCS.R();
#endif

                sigmaTP[faceI] = Foam::symm(R.T() & (sigmaP[faceI] & R));
            }
        }

        sigmaTransformed.write();
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::transformStressToCylindrical::transformStressToCylindrical
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    axis_(dict.lookup("axis")),
    origin_(dict.lookup("origin"))
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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::transformStressToCylindrical::start()
{
    return false;
}


#if FOAMEXTEND
bool Foam::transformStressToCylindrical::execute(const bool forceWrite)
#else
bool Foam::transformStressToCylindrical::execute()
#endif
{
    return writeData();
}


bool Foam::transformStressToCylindrical::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::transformStressToCylindrical::write()
{
    return false;
}
#endif

// ************************************************************************* //
