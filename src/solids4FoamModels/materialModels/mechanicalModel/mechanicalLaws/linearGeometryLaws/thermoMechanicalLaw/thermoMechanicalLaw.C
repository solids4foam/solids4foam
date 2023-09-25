/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "thermoMechanicalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoMechanicalLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, thermoMechanicalLaw, linGeomMechLaw
    );
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::thermoMechanicalLaw::makeTrunTime()
{
    if (TrunTimePtr_.valid())
    {
        FatalErrorIn(type() + "::makeTrunTime()")
            << "Pointer already set!" << abort(FatalError);
    }

    const fileName Tpath = fileName(mesh().time().caseName()/TcaseDir_);

    Info<< "Creating time for T case = " << Tpath << endl;

    TrunTimePtr_.set
    (
        new Time
        (
            Time::controlDictName,
            mesh().time().rootPath(),
            Tpath,
            "system",
            "constant",
            true
#ifdef OPENFOAM_COM
            ,
            true
#endif
        )
    );

    // Set time to be the same as the main case
    TrunTimePtr_().setTime(mesh().time());
}

Foam::Time& Foam::thermoMechanicalLaw::TrunTime()
{
    if (TrunTimePtr_.empty())
    {
        makeTrunTime();
    }

    return TrunTimePtr_();
}


void Foam::thermoMechanicalLaw::makeTmesh()
{
    if (TmeshPtr_.valid())
    {
        FatalErrorIn(type() + "::makeTmesh()")
            << "Pointer already set!" << abort(FatalError);
    }

    Info<< "Read mesh for T case" << endl;

    TmeshPtr_.set
    (
        new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                TrunTime().timeName(),
                TrunTime(),
                IOobject::MUST_READ
            )
        )
    );

    // Check the Tmesh is the same as the mesh (or baseMesh for multi-material
    // cases), in terms of points, faces and cells
    // Note: baseMesh returns mesh in the case of single material cases
    if
    (
        baseMesh().nPoints() != TmeshPtr_().nPoints()
     || baseMesh().nFaces() != TmeshPtr_().nFaces()
     || baseMesh().nCells() != TmeshPtr_().nCells()
    )
    {
        FatalErrorIn(type() + "::readTField()")
            << "The mesh for the T field has a different number of points, "
            << "faces and/or cells than that of the main mesh!"
            << abort(FatalError);
    }
    else if
    (
        gMax(mag(baseMesh().points() - TmeshPtr_().points())) > SMALL
    )
    {
        FatalErrorIn(type() + "::readTField()")
            << "The points in the mesh for the T field are different number "
            << "than those in the main mesh!"
            << abort(FatalError);
    }
}


Foam::fvMesh& Foam::thermoMechanicalLaw::Tmesh()
{
    if (TmeshPtr_.empty())
    {
        makeTmesh();
    }

    return TmeshPtr_();
}


bool Foam::thermoMechanicalLaw::readTField()
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

        // Set TrunTime to be the same as the main case
        TrunTime().setTime(mesh().time());

        IOobject Theader
        (
            "T",
            TrunTime().timeName(),
            Tmesh(),
            IOobject::MUST_READ
        );

#ifdef OPENFOAM_NOT_EXTEND
        if (Theader.typeHeaderOk<volScalarField>())
#else
        if (Theader.headerOk())
#endif
        {
            Info<< nl << "Reading T field from time = "
                << Tmesh().time().timePath() << nl << endl;

            // Read T field from T case
            const volScalarField TField(Theader, Tmesh());

            // Check if a single
            if (mesh() == baseMesh())
            {
                // Copy T field (defined on T mesh) to TPtr_ (defined on mechanical
                // law mesh)
                TPtr_.clear();
                TPtr_.set
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "T",
                            mesh().time().timeName(),
                            mesh(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh(),
                        dimensionedScalar("zero", dimTemperature, 0.0)
                    )
                );

                // Copy internal and boundary fields
#ifdef OPENFOAM_NOT_EXTEND
                TPtr_().primitiveFieldRef() = TField.primitiveField();
#else
                TPtr_().internalField() = TField.internalField();
#endif
                forAll(TField.boundaryField(), patchI)
                {
#ifdef OPENFOAM_NOT_EXTEND
                    TPtr_().boundaryFieldRef()[patchI] =
#else
                    TPtr_().boundaryField()[patchI] =
#endif
                        scalarField(TField.boundaryField()[patchI]);
                }
            }
            else
            {
                // Copy TField to a T field on the baseMesh
                volScalarField TbaseMesh
                (
                    IOobject
                    (
                        "TbaseMesh",
                        baseMesh().time().timeName(),
                        baseMesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    baseMesh(),
                    dimensionedScalar("zero", dimTemperature, 0.0)
                );

                // Copy internal and boundary fields
#ifdef OPENFOAM_NOT_EXTEND
                TbaseMesh.primitiveFieldRef() = TField.primitiveField();
#else
                TbaseMesh.internalField() = TField.internalField();
#endif
                forAll(TField.boundaryField(), patchI)
                {
#ifdef OPENFOAM_NOT_EXTEND
                    TbaseMesh.boundaryFieldRef()[patchI] =
#else
                    TbaseMesh.boundaryField()[patchI] =
#endif
                        scalarField(TField.boundaryField()[patchI]);
                }

                // Map TbaseMesh to T on the current material mesh
                TPtr_.clear();
                TPtr_.set
                (
                    new volScalarField
                    (
                        baseMesh().lookupObject<mechanicalModel>
                        (
                            "mechanicalProperties"
                        ).solSubMeshes().lookupBaseMeshVolField<scalar>
                        (
                            "TbaseMesh", mesh()
                        )()
                    )
                );
            }

            TFieldWasReadFromDisk_ = true;

            // Set the T field to write to disk: this allows a restart where
            // the T field was only specified in the first time-step
            // Also, it is convenient to see the field
            TPtr_().writeOpt() = IOobject::AUTO_WRITE;
        }
    }

    return TFieldWasReadFromDisk_;
}


const Foam::volScalarField&
Foam::thermoMechanicalLaw::lookupTemperatureField()
{
    if (mesh().foundObject<volScalarField>("T") && !TFieldWasReadFromDisk_)
    {
        // Lookup the temperature field from the solver
        return  mesh().lookupObject<volScalarField>("T");
    }
    else if (readTField())
    {
        // Read T from disk
        return TPtr_();
    }
    else
    {
        FatalErrorIn("Foam::poroMechanicalLaw::lookupTemperatureField()")
            << "No T field found in memory or on disk. Make sure you have "
            << "either specified a solidModel that solves for temperature "
            << "or give the T field in at least the starting time "
            << "directory" << abort(FatalError);
    }

    // Keep compiler happy
    return mesh().lookupObject<volScalarField>("null");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::thermoMechanicalLaw::thermoMechanicalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mechLawPtr_
    (
        mechanicalLaw::NewLinGeomMechLaw
        (
            word(dict.subDict("mechanicalLaw").lookup("type")),
            mesh,
            dict.subDict("mechanicalLaw"),
            nonLinGeom
        )
    ),
    alpha_(dict.lookup("alpha")),
    T0_(dict.lookup("T0")),
    TPtr_(),
    TFieldWasReadFromDisk_(false),
    TcaseDir_(mechanicalLaw::dict().lookupOrAddDefault<fileName>("TcaseDirectory", ".")),
    TrunTimePtr_(),
    TmeshPtr_(),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermoMechanicalLaw::~thermoMechanicalLaw()
{}


Foam::tmp<Foam::volScalarField> Foam::thermoMechanicalLaw::impK() const
{
    return mechLawPtr_->impK();
}


void Foam::thermoMechanicalLaw::correct(volSymmTensorField& sigma)
{
    // Calculate stress without thermal stress term
    mechLawPtr_->correct(sigma);

    // Lookup the temperature field
    const volScalarField& T = lookupTemperatureField();

    // Add thermal stress term
    sigma -= 3.0*mechLawPtr_->bulkModulus()*alpha_*(T - T0_)*symmTensor(I);
}


void Foam::thermoMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    // Calculate stress without thermal stress term
    mechLawPtr_->correct(sigma);

    // Lookup the temperature field
    const volScalarField& T = lookupTemperatureField();

    // Interpolate T to the faces
    const surfaceScalarField Tf(fvc::interpolate(T));

    // Interpolate bulk modulus to the faces
    const surfaceScalarField Kf(fvc::interpolate(mechLawPtr_->bulkModulus()));

    // Add thermal stress term
    sigma -= 3.0*Kf*alpha_*(Tf - T0_)*symmTensor(I);
}


// ************************************************************************* //
