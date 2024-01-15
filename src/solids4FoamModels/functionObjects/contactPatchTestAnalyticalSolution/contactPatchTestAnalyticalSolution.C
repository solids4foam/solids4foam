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

#include "contactPatchTestAnalyticalSolution.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "coordinateSystem.H"
#include "symmTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(contactPatchTestAnalyticalSolution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        contactPatchTestAnalyticalSolution,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::symmTensor
Foam::contactPatchTestAnalyticalSolution::contactPatchTestStress()
{
    symmTensor sigma = symmTensor::zero;

    sigma.xx() = 0.0;

    sigma.yy() = -E_ / (1-Foam::sqr(nu_)) * Foam::mag(displacement_);

    sigma.xy() = 0.0;

    sigma.zz() = nu_ * sigma.yy();

    return sigma;
}

bool Foam::contactPatchTestAnalyticalSolution::writeData()
{
    Info<< "contactPatchTestAnalyticalSolution:" << nl
        << "\tWriting analytical solution fields"
        << endl;

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

    volSymmTensorField analyticalStress
    (
        IOobject
        (
           "analyticalStress",
           time_.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    symmTensorField& sI = analyticalStress;

    forAll(sI, cellI)
    {
        sI[cellI] = contactPatchTestStress();
    }

    forAll(analyticalStress.boundaryField(), patchI)
    {
        if (mesh.boundary()[patchI].type() != "empty")
        {
#ifdef OPENFOAM_NOT_EXTEND
            symmTensorField& sP = analyticalStress.boundaryFieldRef()[patchI];
#else
            symmTensorField& sP = analyticalStress.boundaryField()[patchI];
#endif
            forAll(sP, faceI)
            {
                sP[faceI] = contactPatchTestStress();
            }
        }
    }

    if (mesh.foundObject<volSymmTensorField>("sigma"))
    {
        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>("sigma");

        volScalarField relError
        (
            IOobject
            (
                "relativeError",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );

        scalarField& sI = relError;

        const scalar analyticalSigmayy = contactPatchTestStress().yy();

        forAll(sI, cellI)
        {
            sI[cellI] =
                    (Foam::mag(sigma[cellI].yy() - analyticalSigmayy))
                  / Foam::mag(analyticalSigmayy);

            sI[cellI] *= 100;
        }

        forAll(relError.boundaryField(), patchI)
        {
            if (mesh.boundary()[patchI].type() != "empty")
            {
#ifdef OPENFOAM_NOT_EXTEND
                scalarField& sP = relError.boundaryFieldRef()[patchI];
#else
                scalarField& sP = relError.boundaryField()[patchI];
#endif
                forAll(sP, faceI)
                {
                    sP[faceI] =
                        (
                            Foam::mag
                            (
                                sigma.boundaryField()[patchI][faceI].yy()
                              - analyticalSigmayy
                            )
                        ) / Foam::mag(analyticalSigmayy);

                    sP[faceI] *= 100;
                }
            }
        }

        Info<< "\tWriting sigma_y relative error field (in %)" << endl;
        relError.write();

        Info<< "\tAverage relative error in sigma_y field: "
            << gAverage(relError) << "%" << nl << endl;
    }

    analyticalStress.write();

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::contactPatchTestAnalyticalSolution::contactPatchTestAnalyticalSolution
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    displacement_(readScalar(dict.lookup("displacement"))),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu")))
{
    Info<< "Creating " << this->name() << " function object" << endl;

    if (E_ < SMALL || nu_ < SMALL)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "E and nu should be positive!"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::contactPatchTestAnalyticalSolution::start()
{
    return true;
}


#if FOAMEXTEND
bool Foam::contactPatchTestAnalyticalSolution::execute(const bool forceWrite)
#else
bool Foam::contactPatchTestAnalyticalSolution::execute()
#endif
{
    return writeData();
}


bool Foam::contactPatchTestAnalyticalSolution::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::contactPatchTestAnalyticalSolution::write()
{
    return false;
}
#endif

// ************************************************************************* //
