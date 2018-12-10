/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "thermoLinearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "mechanicalModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoLinearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, thermoLinearElastic, linGeomMechLaw
    );
}


// * * * * * * * * * * * * *  Private Data Members * * * * * * * * * * * * * //

bool Foam::thermoLinearElastic::readTField()
{
    if (debug)
    {
        Info<< nl << "Attempting to read T from time = "
            << mesh().time().timeName() << nl << endl;
    }

    // Only attempt to read the T field from disk once per time-step
    if (curTimeIndex_ != mesh().time().timeIndex())
    {
        curTimeIndex_ = mesh().time().timeIndex();

        IOobject Theader
        (
            "T",
            mesh().time().timeName(),
            mesh(),
            IOobject::MUST_READ
        );

        if (Theader.headerOk())
        {
            Info<< nl << "Reading T field from time = " << mesh().time().timeName()
                << nl << endl;

            TPtr_.clear();
            TPtr_.set(new volScalarField(Theader, mesh()));

            TFieldWasReadFromDisk_ = true;

            // Set the T field to write to disk: this allows a restart where the T
            // field was only specified in the first time-step
            TPtr_().writeOpt() = IOobject::AUTO_WRITE;
        }
    }

    return TFieldWasReadFromDisk_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::thermoLinearElastic::thermoLinearElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    linearElastic(name, mesh, dict, nonLinGeom),
    alpha_(dict.lookup("alpha")),
    T0_(dict.lookup("T0")),
    TPtr_(),
    TFieldWasReadFromDisk_(false),
    curTimeIndex_(-1)
{
    // Check if the T field needs to be read from disk
    readTField();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermoLinearElastic::~thermoLinearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::thermoLinearElastic::correct(volSymmTensorField& sigma)
{
    // Calculate linear elastic stress
    linearElastic::correct(sigma);

    if (mesh().foundObject<volScalarField>("T") && !TFieldWasReadFromDisk_)
    {
        // Lookup the temperature field from the solver
        const volScalarField& T = mesh().lookupObject<volScalarField>("T");

        // Add thermal stress component
        sigma -= 3.0*K()*alpha_*(T - T0_)*symmTensor(I);
    }
    else if (readTField())
    {
        // Add thermal stress component
        sigma -= 3.0*K()*alpha_*(TPtr_() - T0_)*symmTensor(I);
    }
    else
    {
        FatalErrorIn(type() + "::correct(...)")
            << "No T field found in memory or on disk. Make sure you have "
            << "either specified a solidModel that solves for temperature "
            << "or give the T field in at least the starting time "
            << "directory" << abort(FatalError);
    }
}


void Foam::thermoLinearElastic::correct(surfaceSymmTensorField& sigma)
{
    // Calculate linear elastic stress
    linearElastic::correct(sigma);

    if (mesh().foundObject<volScalarField>("T") && !TFieldWasReadFromDisk_)
    {
        // Lookup the temperature field from the solver
        const volScalarField& T = mesh().lookupObject<volScalarField>("T");

        // Interpolate the volField temperature to the faces
        const surfaceScalarField Tf = fvc::interpolate(T);

        // Add thermal stress component
        sigma -= 3.0*K()*alpha_*(Tf - T0_)*symmTensor(I);
    }
    else if (readTField())
    {
        // Interpolate the volField temperature to the faces
        const surfaceScalarField Tf = fvc::interpolate(TPtr_());

        // Add thermal stress component
        sigma -= 3.0*K()*alpha_*(Tf - T0_)*symmTensor(I);
    }
    else
    {
        FatalErrorIn(type() + "::correct(...)")
            << "No T field found in memory or on disk. Make sure you have "
            << "either specified a solidModel that solves for temperature "
            << "or give the T field in at least the starting time "
            << "directory" << abort(FatalError);
    }
}


// ************************************************************************* //
