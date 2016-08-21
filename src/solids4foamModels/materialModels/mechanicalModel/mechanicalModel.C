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
#include "fvc.H"
#include "gaussGrad.H"

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


const Foam::PtrList<Foam::volSymmTensorField>&
Foam::mechanicalModel::subMeshSigma() const
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


const Foam::PtrList<Foam::surfaceSymmTensorField>&
Foam::mechanicalModel::subMeshSigmaf() const
{
    if (subMeshSigmaf_.empty())
    {
        calcSubMeshSigmaf();
    }

    return subMeshSigmaf_;
}


void Foam::mechanicalModel::calcSubMeshD() const
{
    if (!subMeshD_.empty())
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshD() const")
            << "pointer list already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        // Sub-meshes should not be defined if there is only one material
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshD() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshD_.setSize(laws.size());

    const PtrList<fvMeshSubset>& subMeshes = this->subMeshes();

    forAll(laws, lawI)
    {
        subMeshD_.set
        (
            lawI,
            new volVectorField
            (
                IOobject
                (
                    "D",
                    subMeshes[lawI].subMesh().time().timeName(),
                    subMeshes[lawI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[lawI].subMesh(),
                dimensionedVector("zero", dimLength, vector::zero)
            )
        );
    }
}


Foam::PtrList<Foam::volVectorField>& Foam::mechanicalModel::subMeshD()
{
    if (subMeshD_.empty())
    {
        calcSubMeshD();
    }

    return subMeshD_;
}


const Foam::PtrList<Foam::volVectorField>&
Foam::mechanicalModel::subMeshD() const
{
    if (subMeshD_.empty())
    {
        calcSubMeshD();
    }

    return subMeshD_;
}


void Foam::mechanicalModel::calcInterfaceShadowIDs() const
{
    if (!interfaceShadowSubMeshID_.empty())
    {
        FatalErrorIn
        (
            "void Foam::mechanicalModel::calcInterfaceShadowIDs() const"
        )   << "pointer already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    interfaceShadowSubMeshID_.setSize(laws.size());
    interfaceShadowPatchID_.setSize(laws.size());
    interfaceShadowFaceID_.setSize(laws.size());

    // Reverse maps for the interface faces
    labelList baseMeshShadSubMeshID = labelList(mesh().nInternalFaces(), -1);
    labelList baseMeshShadPatchID = labelList(mesh().nInternalFaces(), -1);
    labelList baseMeshShadFaceID = labelList(mesh().nInternalFaces(), -1);

    forAll(laws, lawI)
    {
        const fvMesh& subMesh = subMeshes()[lawI].subMesh();

        const labelList& patchMap = subMeshes()[lawI].patchMap();
        const labelList& faceMap = subMeshes()[lawI].faceMap();

        bool uniquePatchFound = false;

        forAll(subMesh.boundaryMesh(), patchI)
        {
            if (patchMap[patchI] == -1)
            {
                if (uniquePatchFound)
                {
                    FatalErrorIn
                    (
                        "void Foam::mechanicalModel::"
                        "calcInterfaceShadowIDs() const"
                    )   << "There are more than one interface patches!"
                        << abort(FatalError);
                }

                uniquePatchFound = true;

                interfaceShadowSubMeshID_.set
                (
                    lawI,
                    new labelList(subMesh.boundaryMesh()[patchI].size(), -1)
                );

                interfaceShadowPatchID_.set
                (
                    lawI,
                    new labelList(subMesh.boundaryMesh()[patchI].size(), -1)
                );

                interfaceShadowFaceID_.set
                (
                    lawI,
                    new labelList(subMesh.boundaryMesh()[patchI].size(), -1)
                );

                labelList& shadSubMeshID = interfaceShadowSubMeshID_[lawI];
                labelList& shadPatchID = interfaceShadowPatchID_[lawI];
                labelList& shadFaceID = interfaceShadowFaceID_[lawI];

                const label start = subMesh.boundaryMesh()[patchI].start();

                forAll(shadSubMeshID, faceI)
                {
                    const label baseMeshFaceID = faceMap[start + faceI];

                    // Check if the face has been set in the baseMesh lists
                    if (baseMeshShadSubMeshID[baseMeshFaceID] == -1)
                    {
                        // Store the local IDs in the baseMesh lists
                        baseMeshShadSubMeshID[baseMeshFaceID] = lawI;
                        baseMeshShadPatchID[baseMeshFaceID] = patchI;
                        baseMeshShadFaceID[baseMeshFaceID] = faceI;
                    }
                    else
                    {
                        // Store the shadow values in the local lists
                        shadSubMeshID[faceI] =
                            baseMeshShadSubMeshID[baseMeshFaceID];
                        shadPatchID[faceI] =
                            baseMeshShadPatchID[baseMeshFaceID];
                        shadFaceID[faceI] =
                            baseMeshShadFaceID[baseMeshFaceID];

                        // Update the shadow lists with the local values
                        const label shadSubMeshID =
                            baseMeshShadSubMeshID[baseMeshFaceID];
                        //const label shadPatchID =
                        //    baseMeshShadPatchID[baseMeshFaceID];
                        const label shadFaceID =
                            baseMeshShadFaceID[baseMeshFaceID];

                        interfaceShadowSubMeshID_[shadSubMeshID][shadFaceID] =
                            lawI;
                        interfaceShadowPatchID_[shadSubMeshID][shadFaceID] =
                            patchI;
                        interfaceShadowFaceID_[shadSubMeshID][shadFaceID] =
                            faceI;
                    }
                }
            }
        }
    }
}


const Foam::PtrList<Foam::labelList>&
Foam::mechanicalModel::interfaceShadowSubMeshID() const
{
    if (interfaceShadowSubMeshID_.empty())
    {
        calcInterfaceShadowIDs();
    }

    return interfaceShadowSubMeshID_;
}


const Foam::PtrList<Foam::labelList>&
Foam::mechanicalModel::interfaceShadowPatchID() const
{
    if (interfaceShadowPatchID_.empty())
    {
        calcInterfaceShadowIDs();
    }

    return interfaceShadowPatchID_;
}


const Foam::PtrList<Foam::labelList>&
Foam::mechanicalModel::interfaceShadowFaceID() const
{
    if (interfaceShadowFaceID_.empty())
    {
        calcInterfaceShadowIDs();
    }

    return interfaceShadowFaceID_;
}


void Foam::mechanicalModel::calcImpKfcorr() const
{
    if (impKfcorrPtr_)
    {
        FatalErrorIn
        (
            "const Foam::volScalarField& "
            "Foam::mechanicalModel::calcImpKfcorr() const"
        )   << "pointer already set" << abort(FatalError);
    }

    impKfcorrPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "impKfcorr",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            impKf()
        );

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() > 1)
    {
        // To disable Rhie-Chow correction on bi-material interface, we will set
        // newImpKf to zero on bi-material interface faces

        scalarField& impKfcorrI = impKfcorrPtr_->internalField();

        forAll(laws, lawI)
        {
            const fvMesh& subMesh = subMeshes()[lawI].subMesh();
            const labelList& patchMap = subMeshes()[lawI].patchMap();
            const labelList& faceMap = subMeshes()[lawI].faceMap();

            forAll(subMesh.boundaryMesh(), patchI)
            {
                if (patchMap[patchI] == -1)
                {
                    const polyPatch& ppatch = subMesh.boundaryMesh()[patchI];
                    const label start = ppatch.start();

                    forAll(ppatch, faceI)
                    {
                        impKfcorrI[faceMap[start + faceI]] = 0.0;
                    }
                }
            }
        }
    }
}


const Foam::surfaceScalarField& Foam::mechanicalModel::impKfcorr() const
{
    if (!impKfcorrPtr_)
    {
        calcImpKfcorr();
    }

    return *impKfcorrPtr_;
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
    subMeshSigmaf_(),
    subMeshD_(),
    interfaceShadowSubMeshID_(),
    interfaceShadowPatchID_(),
    interfaceShadowFaceID_(),
    impKfcorrPtr_()
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
        // Reset sigma before performing the accumulatation as interface values
        // will be added for each material
        // This is not necessary for volFields as they store no value on the
        // interface
        sigma = dimensionedSymmTensor("zero", dimPressure, symmTensor::zero);

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


void Foam::mechanicalModel::grad
(
    volTensorField& gradD,
    const volVectorField& D
)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        gradD = fvc::grad(D);
    }
    else
    {
        // Accumulate data for all fields
        forAll(laws, lawI)
        {
            volVectorField& subMeshD = this->subMeshD()[lawI];

            const fvMesh& subMesh = subMeshes()[lawI].subMesh();
            const labelList& faceMap = subMeshes()[lawI].faceMap();
            const labelList& patchMap = subMeshes()[lawI].patchMap();
            const labelList& cellMap = subMeshes()[lawI].cellMap();

            // Store interface field as it is overwritten with the interpolated
            // value by the interpolate function
            vectorField Dinterface(0);
            forAll(subMeshD.boundaryField(), patchI)
            {
                if (patchMap[patchI] == -1)
                {
                    Dinterface = subMeshD.boundaryField()[patchI];
                }
            }

            // Map the base displacement field to the subMesh; this overwrites
            // the interface with the interpolated values
            subMeshD = subMeshes()[lawI].interpolate(D);

            // Overwrite the interface values with the previous interface values
            forAll(subMeshD.boundaryField(), patchI)
            {
                if (patchMap[patchI] == -1)
                {
                    subMeshD.boundaryField()[patchI] = Dinterface;
                }
            }

            // Calculate the new correction to the interface valaues
            forAll(subMeshD.boundaryField(), patchI)
            {
                if (patchMap[patchI] == -1)
                {
                    // Interface displacement
                    vectorField& Dinterface = subMeshD.boundaryField()[patchI];

                    const label start = subMesh.boundaryMesh()[patchI].start();

                    // Base mesh owner cells
                    const unallocLabelList& baseOwn = mesh().owner();

                    // Base mesh neighbour cells
                    const unallocLabelList& baseNei = mesh().neighbour();

                    // Base mesh face interpolation weights
                    const scalarField& baseWeightsI =
                        mesh().weights().internalField();

                    // Base mesh cell centres
                    const vectorField& baseCI = mesh().C().internalField();

                    // Base mesh face area vectors
                    const vectorField& baseSf = mesh().Sf().internalField();

                    // Base mesh face area vector magnitudes
                    const scalarField& baseMagSf =
                        mesh().magSf().internalField();

                    // Implicit stiffness field
                    const scalarField& KI =
                        mesh().lookupObject<volScalarField>
                        (
                            "impK"
                        ).internalField();

                    // Interface face centres
                    const vectorField& Cf = subMesh.boundary()[patchI].Cf();

                    // Stress in the current subMesh at the interface
                    const symmTensorField& sigmaPatch =
                        subMeshSigma()[lawI].boundaryField()[patchI];

                    const labelList& faceCells =
                        subMesh.boundaryMesh()[patchI].faceCells();

                    // The stress fields in all the subMeshes
                    const PtrList<volSymmTensorField>& subMeshSigma =
                        this->subMeshSigma();

                    const labelList& interfaceShadowSubMeshID =
                        this->interfaceShadowSubMeshID()[lawI];
                    const labelList& interfaceShadowPatchID =
                        this->interfaceShadowPatchID()[lawI];
                    const labelList& interfaceShadowFaceID =
                        this->interfaceShadowFaceID()[lawI];

                    // Calculate the interface displacements
                    forAll(Dinterface, faceI)
                    {
                        // Base mesh face index
                        const label baseFaceID = faceMap[start + faceI];

                        // Base mesh owner cell index
                        const label baseOwnID = baseOwn[baseFaceID];

                        // Base mesh neighbour cell index
                        const label baseNeiID = baseNei[baseFaceID];

                        // Base mesh face interpolation weight
                        const scalar baseW = baseWeightsI[baseFaceID];

                        // Interface unit normal (on base mesh); this may be in
                        // the opposite direction to the subMesh normal
                        const vector n =
                            baseSf[baseFaceID]/baseMagSf[baseFaceID];

                        // Normal distance from the interface to the cell-centre
                        // on side-a
                        const scalar da =
                            mag(n & (Cf[faceI] - baseCI[baseOwnID]));

                        // Normal distance from the interface to the cell-centre
                        // on side-b
                        const scalar db =
                            mag(n & (baseCI[baseNeiID] - Cf[faceI]));

                        // Calculate implicit stiffness for the interface; this
                        // value only affects the convergence and not the result
                        // For now, we ill linearly interpolate the value
                        // Weighted-harmonic interpolation may be a better
                        // candidate
                        const scalar K =
                            baseW*KI[baseOwnID] + (1.0 - baseW)*KI[baseNeiID];

                        // Lookup stress for side-a and side-b
                        symmTensor sigmaa = symmTensor::zero;
                        symmTensor sigmab = symmTensor::zero;

                        // ID of the subMesh on the other side of
                        // interface
                        const label shadowSubMeshID =
                            interfaceShadowSubMeshID[faceI];

                        // ID of the subMesh patch on the other side of
                        // interface
                        const label shadowPatchID =
                            interfaceShadowPatchID[faceI];

                        // ID of the subMesh patch face on the other side of
                        // interface
                        const label shadowFaceID =
                            interfaceShadowFaceID[faceI];

                        if (baseOwnID == cellMap[faceCells[faceI]])
                        {
                            // Stress calculated at the side-a (owner) of the
                            // interface
                            sigmaa = sigmaPatch[faceI];

                            // Stress calculated at the side-b (neighbour) of
                            // the interface
                            sigmab =
                                subMeshSigma
                                [
                                    shadowSubMeshID
                                ].boundaryField()[shadowPatchID][shadowFaceID];
                        }
                        else
                        {
                            // Stress calculated at the side-b (neighbour) of
                            // the interface
                            sigmaa =
                                subMeshSigma
                                [
                                    shadowSubMeshID
                                ].boundaryField()[shadowPatchID][shadowFaceID];

                            // Stress calculated at the side-b (neighbour) of
                            // the interface
                            sigmab = sigmaPatch[faceI];
                        }

                        // Add correction to the interface displacement
                        // This correction goes to zero on convergence
                        Dinterface[faceI] +=
                            ((da*db)/(da + db))*(n & (sigmab - sigmaa)/K);
                    }
                }
            }

            // Calculate gradient on subMesh
            // This will use the values at the interface
            volTensorField subMeshGradD = fvc::grad(subMeshD);

            // Map subMesh gradD to the base gradient field
            mapSubMeshVolField<tensor>
            (
                lawI, subMeshGradD, gradD
            );
        }

        // Correct boundary snGrad
        fv::gaussGrad<vector>
        (
            mesh()
        ).correctBoundaryConditions(D, gradD);
    }
}


Foam::tmp<Foam::volVectorField> Foam::mechanicalModel::RhieChowCorrection
(
    const volVectorField& D,
    const volTensorField& gradD
) const
{
    return
    (
        fvc::laplacian(impKfcorr(), D, "laplacian(DD,D)")
      - fvc::div(impKfcorr()*mesh().Sf() & fvc::interpolate(gradD))
    );
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
