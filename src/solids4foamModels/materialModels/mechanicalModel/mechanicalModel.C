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

#include "mechanicalModel.H"
#include "ZoneID.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::mechanicalModel::makeSubMeshes() const
{
    if (!subMeshes_.empty())
    {
        FatalErrorIn("void Foam::mechanicalModel::makeSubMeshes() const")
            << "sub-meshes already exist" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        // Sub-meshes should not be defined if there is only one material
        FatalErrorIn("void Foam::mechanicalModel::makeSubMeshes() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    // Check that each cell is in exactly one cellZone
    checkCellZones();

    labelList region(mesh().nCells(), -1);

    forAll(laws, lawI)
    {
        const ZoneID<cellZone> cellZoneID =
            ZoneID<cellZone>(cellZoneNames_[lawI], mesh().cellZones());

        if (!cellZoneID.active())
        {
            FatalErrorIn("void Foam::mechanicalModel::makeSubMeshes() const")
                << "cellZone not found for material " << laws[lawI].name()
                << abort(FatalError);
        }

        const labelList& curCellZone = mesh().cellZones()[cellZoneID.index()];

        forAll(curCellZone, cI)
        {
            region[curCellZone[cI]] = lawI;
        }
    }

    subMeshes_.setSize(laws.size());

    forAll(subMeshes_, matI)
    {
        subMeshes_.set
        (
            matI,
            new fvMeshSubset
            (
                IOobject
                (
                    Foam::name(Pstream::myProcNo()) + '_' + Foam::name(matI),
                    mesh().time().constant(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh()
            )
        );

        subMeshes_[matI].setLargeCellSubset(region, matI);
    }
}

const Foam::PtrList<Foam::fvMeshSubset>&
Foam::mechanicalModel::subMeshes() const
{
    if (subMeshes_.empty())
    {
        makeSubMeshes();
    }

    return subMeshes_;
}


void Foam::mechanicalModel::checkCellZones() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        // Cell zones need not be defined if there is only one material
        return;
    }

    // We will check that every cell is in exactly one cellZone
    volScalarField nCellZones
    (
        IOobject
        (
            "nCellZones",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    );

    scalarField& nCellZonesI = nCellZones.internalField();

    forAll(laws, lawI)
    {
        const ZoneID<cellZone> cellZoneID =
            ZoneID<cellZone>(cellZoneNames_[lawI], mesh().cellZones());

        if (!cellZoneID.active())
        {
            FatalErrorIn("void Foam::mechanicalModel::checkCellZones()")
                << "cellZone " << cellZoneNames_[lawI]
                << " not found for material " << cellZoneNames_[lawI]
                << abort(FatalError);
        }

        const labelList& curCellZone = mesh().cellZones()[cellZoneID.index()];

        forAll(curCellZone, cI)
        {
            const label cellID = curCellZone[cI];
            nCellZonesI[cellID] = nCellZonesI[cellID] + 1.0;
        }
    }

    if (mag(gMin(nCellZonesI)) < SMALL)
    {
        FatalErrorIn("void Foam::mechanicalModel::checkCellZones()")
            << "There are cells that are not in a material cellZone!"
            << abort(FatalError);
    }

    if (mag(gMax(nCellZonesI) - 1) > SMALL)
    {
        FatalErrorIn("void Foam::mechanicalModel::checkCellZones()")
            << "There are cells that are in more than one material cellZone!"
            << abort(FatalError);
    }
}


void Foam::mechanicalModel::calcSubMeshSigma() const
{
    if (!subMeshSigma_.empty())
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshSigma() const")
            << "pointer list already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        // Sub-meshes should not be defined if there is only one material
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshSigma() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshSigma_.setSize(laws.size());

    const PtrList<fvMeshSubset>& subMeshes = this->subMeshes();

    forAll(laws, lawI)
    {
        subMeshSigma_.set
        (
            lawI,
            new volSymmTensorField
            (
                IOobject
                (
                    "sigma",
                    subMeshes[lawI].subMesh().time().timeName(),
                    subMeshes[lawI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[lawI].subMesh(),
                dimensionedSymmTensor
                (
                    "zero",
                    dimForce/dimArea,
                    symmTensor::zero
                )
            )
        );
    }
}


Foam::PtrList<Foam::volSymmTensorField>& Foam::mechanicalModel::subMeshSigma()
{
    if (subMeshSigma_.empty())
    {
        calcSubMeshSigma();
    }

    return subMeshSigma_;
}


void Foam::mechanicalModel::calcSubMeshSigmaf() const
{
    if (!subMeshSigmaf_.empty())
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshSigmaf() const")
            << "pointer list already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        // Sub-meshes should not be defined if there is only one material
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshSigmaf() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshSigmaf_.setSize(laws.size());

    const PtrList<fvMeshSubset>& subMeshes = this->subMeshes();

    forAll(laws, lawI)
    {
        subMeshSigmaf_.set
        (
            lawI,
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "sigmaf",
                    subMeshes[lawI].subMesh().time().timeName(),
                    subMeshes[lawI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[lawI].subMesh(),
                dimensionedSymmTensor
                (
                    "zero",
                    dimForce/dimArea,
                    symmTensor::zero
                )
            )
        );
    }
}


Foam::PtrList<Foam::surfaceSymmTensorField>&
Foam::mechanicalModel::subMeshSigmaf()
{
    if (subMeshSigmaf_.empty())
    {
        calcSubMeshSigmaf();
    }

    return subMeshSigmaf_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalModel::mechanicalModel(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "mechanicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    PtrList<mechanicalLaw>(),
    mesh_(mesh),
    planeStress_(lookup("planeStress")),
    cellZoneNames_(),
    subMeshes_(),
    subMeshSigma_(),
    subMeshSigmaf_()
{
    Info<< "Creating the mechanicalModel" << endl;

    // Read the mechanical laws
    const PtrList<entry> lawEntries(lookup("mechanical"));

    PtrList<mechanicalLaw>& laws = *this;
    laws.setSize(lawEntries.size());

    if (laws.size() == 1)
    {
        laws.set
        (
            0,
            mechanicalLaw::New
            (
                lawEntries[0].keyword(),
                mesh,
                lawEntries[0].dict()
            )
        );
    }
    else
    {
        // We must create the list of cellZones names before creating the
        // subMeshes as they are used during the construction of the subMeshes
        cellZoneNames_.setSize(laws.size());

        forAll(laws, lawI)
        {
            cellZoneNames_[lawI] = lawEntries[lawI].keyword();
        }

        forAll(laws, lawI)
        {
            laws.set
            (
                lawI,
                mechanicalLaw::New
                (
                    lawEntries[lawI].keyword(),
                    subMeshes()[lawI].subMesh(),
                    lawEntries[lawI].dict()
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mechanicalModel::~mechanicalModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::mechanicalModel::mesh() const
{
    return mesh_;
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::rho() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        return laws[0].rho();
    }
    else
    {
        // Accumulate data for all fields
        tmp<volScalarField> tresult
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimDensity, 0),
                calculatedFvPatchScalarField::typeName
            )
        );
        volScalarField& result = tresult();

        forAll(laws, lawI)
        {
            // Map subMesh field to the base field
            const volScalarField curRho = laws[lawI].rho();
            mapSubMeshVolField<scalar>(lawI, curRho, result);
        }

        return tresult;
    }
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::impK() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        return laws[0].impK();
    }
    else
    {
        // Accumulate data for all fields
        tmp<volScalarField> tresult
        (
            new volScalarField
            (
                IOobject
                (
                    "impK",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimForce/dimArea, 0),
                calculatedFvPatchScalarField::typeName
            )
        );
        volScalarField& result = tresult();

        forAll(laws, lawI)
        {
            // Map subMesh field to the base field
            const volScalarField curImpK = laws[lawI].impK();
            mapSubMeshVolField<scalar>(lawI, curImpK, result);
        }

        return tresult;
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::mechanicalModel::impKf() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        return laws[0].impKf();
    }
    else
    {
        // Accumulate data for all fields
        tmp<surfaceScalarField> tresult
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "impKf",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimForce/dimArea, 0),
                calculatedFvPatchScalarField::typeName
            )
        );
        surfaceScalarField& result = tresult();

        forAll(laws, lawI)
        {
            // Map subMesh field to the base field
            const surfaceScalarField curImpKf = laws[lawI].impKf();
            mapSubMeshSurfaceField<scalar>(lawI, curImpKf, result);
        }

        return tresult;
    }
}


void Foam::mechanicalModel::correct(volSymmTensorField& sigma)
{
    PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        laws[0].correct(sigma);
    }
    else
    {
        // Accumulate data for all fields
        forAll(laws, lawI)
        {
            laws[lawI].correct(subMeshSigma()[lawI]);

            // Map subMesh stress to the base stress field
            mapSubMeshVolField<symmTensor>
            (
                lawI, subMeshSigma()[lawI], sigma
            );
        }
    }
}


void Foam::mechanicalModel::correct(surfaceSymmTensorField& sigma)
{
    PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        laws[0].correct(sigma);
    }
    else
    {
        // Accumulate data for all fields
        forAll(laws, lawI)
        {
            laws[lawI].correct(subMeshSigmaf()[lawI]);

            // Map subMesh stress to the base stress field
            mapSubMeshSurfaceField<symmTensor>
            (
                lawI, subMeshSigmaf()[lawI], sigma
            );
        }
    }
}


Foam::scalar Foam::mechanicalModel::residual()
{
    PtrList<mechanicalLaw>& laws = *this;

    scalar maxResidual = 0.0;

    forAll(laws, lawI)
    {
        maxResidual = max(maxResidual, laws[lawI].residual());
    }

    return maxResidual;
}


void Foam::mechanicalModel::updateTotalFields()
{
    PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        laws[lawI].updateTotalFields();
    }
}


// ************************************************************************* //
