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
#include "fvcGradf.H"
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
        word subMeshName = cellZoneNames_[matI];

        if (Pstream::parRun())
        {
            subMeshName =
                "proc" + Foam::name(Pstream::myProcNo()) + "_"
               + cellZoneNames_[matI];
        }

        subMeshes_.set
        (
            matI,
            new newFvMeshSubset
            (
                IOobject
                (
                    subMeshName,
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

const Foam::PtrList<Foam::newFvMeshSubset>&
Foam::mechanicalModel::subMeshes() const
{
    if (subMeshes_.empty())
    {
        makeSubMeshes();
    }

    return subMeshes_;
}


Foam::PtrList<Foam::newFvMeshSubset>& Foam::mechanicalModel::subMeshes()
{
    if (subMeshes_.empty())
    {
        makeSubMeshes();
    }

    return subMeshes_;
}


void Foam::mechanicalModel::makeVolToPoint() const
{
    if (volToPointPtr_)
    {
        FatalErrorIn
        (
            "void Foam::mechanicalModel::makeVolToPoint() const"
        )   << "pointer already set" << abort(FatalError);
    }

    volToPointPtr_ = new newLeastSquaresVolPointInterpolation(mesh());
}


void Foam::mechanicalModel::makeSubMeshVolToPoint() const
{
    if (!subMeshVolToPoint_.empty())
    {
        FatalErrorIn
        (
            "void Foam::mechanicalModel::makeSubMeshVolToPoint() const"
        )   << "sub-meshes already exist" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        FatalErrorIn
        (
            "void Foam::mechanicalModel::makeSubMeshVolToPoint() const"
        )   << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshVolToPoint_.setSize(laws.size());

    forAll(subMeshVolToPoint_, lawI)
    {
        subMeshVolToPoint_.set
        (
            lawI,
            new newLeastSquaresVolPointInterpolation
            (
                subMeshes()[lawI].subMesh()
            )
        );
    }
}


const Foam::PtrList<Foam::newLeastSquaresVolPointInterpolation>&
Foam::mechanicalModel::subMeshVolToPoint() const
{
    if (subMeshVolToPoint_.empty())
    {
        makeSubMeshVolToPoint();
    }

    return subMeshVolToPoint_;
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

    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

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

    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

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
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshD() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshD_.setSize(laws.size());

    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

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


void Foam::mechanicalModel::calcSubMeshDD() const
{
    if (!subMeshDD_.empty())
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshDD() const")
            << "pointer list already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshDD() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshDD_.setSize(laws.size());

    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

    forAll(laws, lawI)
    {
        subMeshDD_.set
        (
            lawI,
            new volVectorField
            (
                IOobject
                (
                    "DD",
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


Foam::PtrList<Foam::volVectorField>& Foam::mechanicalModel::subMeshDD()
{
    if (subMeshDD_.empty())
    {
        calcSubMeshDD();
    }

    return subMeshDD_;
}


const Foam::PtrList<Foam::volVectorField>&
Foam::mechanicalModel::subMeshDD() const
{
    if (subMeshDD_.empty())
    {
        calcSubMeshDD();
    }

    return subMeshDD_;
}


void Foam::mechanicalModel::calcSubMeshGradD() const
{
    if (!subMeshGradD_.empty())
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshGradD() const")
            << "pointer list already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshGradD() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshGradD_.setSize(laws.size());

    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

    forAll(laws, lawI)
    {
        subMeshGradD_.set
        (
            lawI,
            new volTensorField
            (
                IOobject
                (
                    "grad(D)",
                    subMeshes[lawI].subMesh().time().timeName(),
                    subMeshes[lawI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[lawI].subMesh(),
                dimensionedTensor("zero", dimless, tensor::zero)
            )
        );
    }
}


Foam::PtrList<Foam::volTensorField>& Foam::mechanicalModel::subMeshGradD()
{
    if (subMeshGradD_.empty())
    {
        calcSubMeshGradD();
    }

    return subMeshGradD_;
}


const Foam::PtrList<Foam::volTensorField>&
Foam::mechanicalModel::subMeshGradD() const
{
    if (subMeshGradD_.empty())
    {
        calcSubMeshGradD();
    }

    return subMeshGradD_;
}


void Foam::mechanicalModel::calcSubMeshGradDf() const
{
    if (!subMeshGradDf_.empty())
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshGradDf() const")
            << "pointer list already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshGradDf() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshGradDf_.setSize(laws.size());

    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

    forAll(laws, lawI)
    {
        subMeshGradDf_.set
        (
            lawI,
            new surfaceTensorField
            (
                IOobject
                (
                    "grad(D)f",
                    subMeshes[lawI].subMesh().time().timeName(),
                    subMeshes[lawI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[lawI].subMesh(),
                dimensionedTensor("zero", dimless, tensor::zero)
            )
        );
    }
}


Foam::PtrList<Foam::surfaceTensorField>& Foam::mechanicalModel::subMeshGradDf()
{
    if (subMeshGradDf_.empty())
    {
        calcSubMeshGradDf();
    }

    return subMeshGradDf_;
}


const Foam::PtrList<Foam::surfaceTensorField>&
Foam::mechanicalModel::subMeshGradDf() const
{
    if (subMeshGradDf_.empty())
    {
        calcSubMeshGradDf();
    }

    return subMeshGradDf_;
}




void Foam::mechanicalModel::calcSubMeshGradDD() const
{
    if (!subMeshGradDD_.empty())
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshGradDD() const")
            << "pointer list already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshGradDD() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshGradDD_.setSize(laws.size());

    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

    forAll(laws, lawI)
    {
        subMeshGradDD_.set
        (
            lawI,
            new volTensorField
            (
                IOobject
                (
                    "grad(DD)",
                    subMeshes[lawI].subMesh().time().timeName(),
                    subMeshes[lawI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[lawI].subMesh(),
                dimensionedTensor("zero", dimless, tensor::zero)
            )
        );
    }
}


Foam::PtrList<Foam::volTensorField>& Foam::mechanicalModel::subMeshGradDD()
{
    if (subMeshGradDD_.empty())
    {
        calcSubMeshGradDD();
    }

    return subMeshGradDD_;
}


const Foam::PtrList<Foam::volTensorField>&
Foam::mechanicalModel::subMeshGradDD() const
{
    if (subMeshGradDD_.empty())
    {
        calcSubMeshGradDD();
    }

    return subMeshGradDD_;
}


void Foam::mechanicalModel::calcSubMeshGradDDf() const
{
    if (!subMeshGradDDf_.empty())
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshGradDDf() const")
            << "pointer list already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshGradDDf() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshGradDDf_.setSize(laws.size());

    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

    forAll(laws, lawI)
    {
        subMeshGradDDf_.set
        (
            lawI,
            new surfaceTensorField
            (
                IOobject
                (
                    "grad(DD)f",
                    subMeshes[lawI].subMesh().time().timeName(),
                    subMeshes[lawI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[lawI].subMesh(),
                dimensionedTensor("zero", dimless, tensor::zero)
            )
        );
    }
}


Foam::PtrList<Foam::surfaceTensorField>& Foam::mechanicalModel::subMeshGradDDf()
{
    if (subMeshGradDDf_.empty())
    {
        calcSubMeshGradDDf();
    }

    return subMeshGradDDf_;
}


const Foam::PtrList<Foam::surfaceTensorField>&
Foam::mechanicalModel::subMeshGradDDf() const
{
    if (subMeshGradDDf_.empty())
    {
        calcSubMeshGradDDf();
    }

    return subMeshGradDDf_;
}


void Foam::mechanicalModel::calcSubMeshPointD() const
{
    if (!subMeshPointD_.empty())
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshPointD() const")
            << "pointer list already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        FatalErrorIn("void Foam::mechanicalModel::calcSubMeshPointD() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshPointD_.setSize(laws.size());

    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

    forAll(laws, lawI)
    {
        subMeshPointD_.set
        (
            lawI,
            new pointVectorField
            (
                IOobject
                (
                    "pointD",
                    subMeshes[lawI].subMesh().time().timeName(),
                    subMeshes[lawI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[lawI].subPointMesh(),
                dimensionedVector("zero", dimLength, vector::zero)
            )
        );
    }
}


Foam::PtrList<Foam::pointVectorField>& Foam::mechanicalModel::subMeshPointD()
{
    if (subMeshPointD_.empty())
    {
        calcSubMeshPointD();
    }

    return subMeshPointD_;
}


const Foam::PtrList<Foam::pointVectorField>&
Foam::mechanicalModel::subMeshPointD() const
{
    if (subMeshPointD_.empty())
    {
        calcSubMeshPointD();
    }

    return subMeshPointD_;
}


void Foam::mechanicalModel::calcInterfaceShadowIDs() const
{
    if
    (
        !interfaceShadowSubMeshID_.empty()
     || !interfaceShadowPatchID_.empty()
     || !interfaceShadowFaceID_.empty()
    )
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


void Foam::mechanicalModel::makeInterfaceBaseFaces() const
{
    if (interfaceBaseFacesPtr_)
    {
        FatalErrorIn
        (
            "void Foam::mechanicalModel::makeInterfaceBaseFaces() const"
        )   << "pointer already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() > 1)
    {
        // Create material index field from cellZones

        volScalarField materials
        (
            IOobject
            (
                "materials",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("0", dimless, 0)
        );

        scalarField& materialsI = materials.internalField();

        forAll(laws, lawI)
        {
            const ZoneID<cellZone> cellZoneID =
                ZoneID<cellZone>(cellZoneNames_[lawI], mesh().cellZones());

            if (!cellZoneID.active())
            {
                FatalErrorIn
                (
                    "void Foam::mechanicalModel::makeInterfaceBaseFaces() const"
                )   << "cellZone " << cellZoneNames_[lawI]
                    << " not found for material " << cellZoneNames_[lawI]
                    << abort(FatalError);
            }

            const labelList& curCellZone =
                mesh().cellZones()[cellZoneID.index()];

            forAll(curCellZone, cI)
            {
                const label cellID = curCellZone[cI];
                materialsI[cellID] = lawI;
            }
        }

        const unallocLabelList& owner = mesh().owner();
        const unallocLabelList& neighbour = mesh().neighbour();

        labelHashSet interFacesSet;

        forAll(neighbour, faceI)
        {
            if
            (
                mag(materialsI[neighbour[faceI]] - materialsI[owner[faceI]])
              > SMALL
            )
            {
                interFacesSet.insert(faceI);
            }
        }

        forAll(materials.boundaryField(), patchI)
        {
            if (mesh().boundary()[patchI].type() == processorFvPatch::typeName)
            {
                const scalarField ownMat =
                    materials.boundaryField()[patchI].patchInternalField();

                const scalarField ngbMat =
                    materials.boundaryField()[patchI].patchNeighbourField();

                forAll(ownMat, faceI)
                {
                    if (mag(ownMat[faceI] - ngbMat[faceI]) > SMALL)
                    {
                        const label globalFaceID =
                            mesh().boundaryMesh()[patchI].start() + faceI;

                        interFacesSet.insert(globalFaceID);
                    }
                }
            }
        }

        interfaceBaseFacesPtr_ = new labelList(interFacesSet.toc());
    }
    else
    {
        interfaceBaseFacesPtr_ = new labelList(0);
    }
}


const Foam::labelList& Foam::mechanicalModel::interfaceBaseFaces() const
{
    if (!interfaceBaseFacesPtr_)
    {
        makeInterfaceBaseFaces();
    }

    return *interfaceBaseFacesPtr_;
}


void Foam::mechanicalModel::makePointNumOfMaterials() const
{
    if (pointNumOfMaterialsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::mechanicalModel::makePointNumOfMaterials() const"
        )   << "Pointer already set" << abort(FatalError);
    }

    pointNumOfMaterialsPtr_ = new labelList(mesh().nPoints(), 0);
    labelList& pointNumOfMaterials = *pointNumOfMaterialsPtr_;

    const PtrList<mechanicalLaw>& laws = *this;

    // Create material index field from cellZones

    scalarField materialsI(mesh().nCells(), 0);

    forAll(laws, lawI)
    {
        const ZoneID<cellZone> cellZoneID =
            ZoneID<cellZone>(cellZoneNames_[lawI], mesh().cellZones());

        if (!cellZoneID.active())
        {
            FatalErrorIn
            (
                "void Foam::mechanicalModel::"
                "makeIsolatedInterfacePoints() const"
            )   << "cellZone " << cellZoneNames_[lawI]
                << " not found for material " << cellZoneNames_[lawI]
                << abort(FatalError);
        }

        const labelList& curCellZone = mesh().cellZones()[cellZoneID.index()];

        forAll(curCellZone, cI)
        {
            const label cellID = curCellZone[cI];
            materialsI[cellID] = lawI;
        }
    }

    const labelListList& pointCells = mesh().pointCells();

    forAll(pointNumOfMaterials, pointI)
    {
        // Count the number of unique materials in adjacent cells
        const labelList& curCells = pointCells[pointI];

        labelHashSet matSet;

        forAll(curCells, cellI)
        {
            if (!matSet.found(materialsI[curCells[cellI]]))
            {
                matSet.insert(materialsI[curCells[cellI]]);
            }
        }

        pointNumOfMaterials[pointI] = matSet.toc().size();
    }
}


const Foam::labelList& Foam::mechanicalModel::pointNumOfMaterials() const
{
    if (!pointNumOfMaterialsPtr_)
    {
        makePointNumOfMaterials();
    }

    return *pointNumOfMaterialsPtr_;
}


void Foam::mechanicalModel::makeIsolatedInterfacePoints() const
{
    if (isolatedInterfacePointsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::mechanicalModel::makeIsolatedInterfacePoints() const"
        )   << "pointer already set" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    // Create material index field from cellZones

    scalarField materialsI(mesh().nCells(), 0);

    forAll(laws, lawI)
    {
        const ZoneID<cellZone> cellZoneID =
            ZoneID<cellZone>(cellZoneNames_[lawI], mesh().cellZones());

        if (!cellZoneID.active())
        {
            FatalErrorIn
            (
                "void Foam::mechanicalModel::"
                "makeIsolatedInterfacePoints() const"
            )   << "cellZone " << cellZoneNames_[lawI]
                << " not found for material " << cellZoneNames_[lawI]
                << abort(FatalError);
        }

        const labelList& curCellZone = mesh().cellZones()[cellZoneID.index()];

        forAll(curCellZone, cI)
        {
            const label cellID = curCellZone[cI];
            materialsI[cellID] = lawI;
        }
    }

    pointMesh pMesh(mesh());

    pointScalarField pointMaterials
    (
        IOobject
        (
            "pointMaterials",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh,
        dimensionedScalar("0", dimless, 0)
    );
    scalarField& pointMaterialsI = pointMaterials.internalField();

    const labelListList& pointCells = mesh().pointCells();

    forAll(pointMaterialsI, pointI)
    {
        const labelList& curPointCells = pointCells[pointI];

        forAll(curPointCells, cellI)
        {
            label curCell = curPointCells[cellI];
            pointMaterialsI[pointI] += materialsI[curCell];
        }

        pointMaterialsI[pointI] /= curPointCells.size() + SMALL;
    }

    pointMaterials.correctBoundaryConditions();

    scalarField matInter(pointMaterialsI.size(), 0);

    forAll(pointMaterialsI, pointI)
    {
        const labelList& curPointCells = pointCells[pointI];

        if
        (
            mag(pointMaterialsI[pointI] - materialsI[curPointCells[0]])
          > SMALL
        )
        {
            matInter[pointI] = 1;
        }
    }

    labelHashSet isolatedPointsSet;

    const labelList& noMat = pointNumOfMaterials();

    const labelList& spLabels =
        mesh().globalData().sharedPointLabels();

    const labelListList& pointFaces = mesh().pointFaces();
    forAll(matInter, pointI)
    {
        if (matInter[pointI] && (noMat[pointI] == 1))
        {
            const bool sharedPoint(findIndex(spLabels, pointI) != -1);

            if (!sharedPoint)
            {
                bool hasProcessorFace = false;

                const labelList& curPointFaces = pointFaces[pointI];
                forAll(curPointFaces, faceI)
                {
                    label faceID = curPointFaces[faceI];
                    label patchID =
                        mesh().boundaryMesh().whichPatch(faceID);

                    if (patchID != -1)
                    {
                        if
                        (
                            isA<processorPolyPatch>
                            (
                                mesh().boundaryMesh()[patchID]
                            )
                        )
                        {
                            if (findIndex(interfaceBaseFaces(), faceID) == -1)
                            {
                                hasProcessorFace = true;
                                break;
                            }
                        }
                    }
                }

                if (hasProcessorFace)
                {
                    isolatedPointsSet.insert(pointI);
                }
            }
        }
    }

    isolatedInterfacePointsPtr_ = new labelList(isolatedPointsSet.toc());
}


const Foam::labelList& Foam::mechanicalModel::isolatedInterfacePoints() const
{
    if (!isolatedInterfacePointsPtr_)
    {
        makeIsolatedInterfacePoints();
    }

    return *isolatedInterfacePointsPtr_;
}


void Foam::mechanicalModel::interpolateDtoSubMeshD
(
    const volVectorField& D,
    const bool useVolFieldSigma
)
{
    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

    if (interfaceBaseFaces().size() == 0)
    {
        // No need for any corrections if there are no bi-material interfaces
        forAll(subMeshes, lawI)
        {
            if (D.name() == "D")
            {
                subMeshD()[lawI] = subMeshes[lawI].interpolate(D);
            }
            else
            {
                subMeshDD()[lawI] = subMeshes[lawI].interpolate(D);
            }
        }

        return;
    }

    forAll(subMeshes, lawI)
    {
        volVectorField* subMeshDPtr = NULL;
        if (D.name() == "D")
        {
            subMeshDPtr = &(this->subMeshD()[lawI]);
        }
        else
        {
            subMeshDPtr = &(this->subMeshDD()[lawI]);
        }
        volVectorField& subMeshD = *subMeshDPtr;

        const fvMesh& subMesh = subMeshes[lawI].subMesh();

        const labelList& faceMap = subMeshes[lawI].faceMap();
        const labelList& patchMap = subMeshes[lawI].patchMap();
        const labelList& cellMap = subMeshes[lawI].cellMap();

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
        subMeshD = subMeshes[lawI].interpolate(D);

        // Overwrite the interface values with the previous interface values
        forAll(subMeshD.boundaryField(), patchI)
        {
            if (patchMap[patchI] == -1)
            {
                subMeshD.boundaryField()[patchI] = Dinterface;
            }
        }

        // Check if a large strain procedure is being used, if so we must
        // calculate the deformed normals
        // If the deformation gradient field 'F' is found, we will assume it is
        // a large/finite strain procedure
        // For finite strain procedures, we will look up the deformation
        // gradient: relative deformation gradient for updated Lagrangian
        // approaches and the total deformation gradient for total
        // approaches
        const bool useDeformedNormals = mesh().foundObject<volTensorField>("F");
        const volTensorField* FinvPtr = NULL;
        const volScalarField* JPtr = NULL;
        if (useDeformedNormals)
        {
            if (mesh().foundObject<volTensorField>("relF"))
            {
                // Updated Lagrangian approach: use the inverse of the relative
                // deformation gradient
                FinvPtr = &(mesh().lookupObject<volTensorField>("relFinv"));
                JPtr = &(mesh().lookupObject<volScalarField>("relJ"));
            }
            else
            {
                // Total Lagrangian approach: use the inverse of the total
                // deformation gradient
                FinvPtr = &(mesh().lookupObject<volTensorField>("Finv"));
                JPtr = &(mesh().lookupObject<volScalarField>("J"));
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
                const symmTensorField* sigmaPatchPtr = NULL;
                if (useVolFieldSigma)
                {
                    sigmaPatchPtr =
                        &(subMeshSigma()[lawI].boundaryField()[patchI]);
                }
                else
                {
                    sigmaPatchPtr =
                        &(subMeshSigmaf()[lawI].boundaryField()[patchI]);
                }
                const symmTensorField& sigmaPatch = *sigmaPatchPtr;

                const labelList& faceCells =
                    subMesh.boundaryMesh()[patchI].faceCells();

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
                    vector n = baseSf[baseFaceID]/baseMagSf[baseFaceID];
                    if (useDeformedNormals)
                    {
                        // Interpolate Finv and J to the face
                        const tensorField& FinvI = FinvPtr->internalField();
                        const tensor Finv =
                            baseW*FinvI[baseOwn[baseFaceID]]
                          + (1.0 - baseW)*FinvI[baseOwn[baseFaceID]];

                        const scalarField& JI = JPtr->internalField();
                        const scalar J =
                            baseW*JI[baseOwn[baseFaceID]]
                          + (1.0 - baseW)*JI[baseOwn[baseFaceID]];

                        // Nanson's formula
                        // Note: for updated Lagrangian approach, F is the
                        // relative deformation gradient, whereas for total
                        // Lagrangian approaches, it is the total deformation
                        // gradient
                        n = J*Finv.T() & n;
                        n /= mag(n);
                    }

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
                    // For now, we will linearly interpolate the value
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
                        if (useVolFieldSigma)
                        {
                            sigmab =
                                subMeshSigma()
                                [
                                    shadowSubMeshID
                                ].boundaryField()[shadowPatchID][shadowFaceID];
                        }
                        else
                        {
                            sigmab =
                                subMeshSigmaf()
                                [
                                    shadowSubMeshID
                                ].boundaryField()[shadowPatchID][shadowFaceID];
                        }
                    }
                    else
                    {
                        // Stress calculated at the side-b (neighbour) of
                        // the interface
                        if (useVolFieldSigma)
                        {
                            sigmaa =
                                subMeshSigma()
                                [
                                    shadowSubMeshID
                                ].boundaryField()[shadowPatchID][shadowFaceID];
                        }
                        else
                        {
                            sigmaa =
                                subMeshSigmaf()
                                [
                                    shadowSubMeshID
                                ].boundaryField()[shadowPatchID][shadowFaceID];
                        }

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
    }
}


void Foam::mechanicalModel::correctInterfaceSnGrad
(
    PtrList<volVectorField>& subMeshDList,
    PtrList<volTensorField>& subMeshGradDList
)
{
    if (interfaceBaseFaces().size() == 0)
    {
        // No need for any corrections if there are no bi-material interfaces
        return;
    }

    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

    forAll(subMeshes, lawI)
    {
        const volVectorField& subMeshD = subMeshDList[lawI];
        volTensorField& subMeshGradD = subMeshGradDList[lawI];
        const fvMesh& subMesh = subMeshes[lawI].subMesh();
        const labelList& patchMap = subMeshes[lawI].patchMap();

        forAll(subMeshGradD.boundaryField(), patchI)
        {
            if (patchMap[patchI] == -1)
            {
                tensorField& patchGradD = subMeshGradD.boundaryField()[patchI];
                const tensorField patchGradDif =
                    subMeshGradD.boundaryField()[patchI].patchInternalField();

                const vectorField& patchD = subMeshD.boundaryField()[patchI];
                const vectorField patchDif =
                    subMeshD.boundaryField()[patchI].patchInternalField();

                const vectorField n = subMesh.boundary()[patchI].nf();
                const vectorField delta = subMesh.boundary()[patchI].delta();
                const vectorField k = delta - n*(n & delta);
                const scalarField& deltaCoeffs =
                    subMesh.boundary()[patchI].deltaCoeffs();

                const vectorField correctedSnGrad =
                    (patchD - (patchDif + (k & patchGradDif)))*deltaCoeffs;

                patchGradD += n*(correctedSnGrad - (n & patchGradD));
            }
        }
    }
}


void Foam::mechanicalModel::correctInterfaceSnGradf
(
    PtrList<volVectorField>& subMeshDList,
    PtrList<surfaceTensorField>& subMeshGradDfList,
    PtrList<volTensorField>& subMeshGradDList
)
{
    if (interfaceBaseFaces().size() == 0)
    {
        // No need for any corrections if there are no bi-material interfaces
        return;
    }

    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

    forAll(subMeshes, lawI)
    {
        const volVectorField& subMeshD = subMeshDList[lawI];
        surfaceTensorField& subMeshGradDf = subMeshGradDfList[lawI];
        const volTensorField& subMeshGradD = subMeshGradDList[lawI];
        const fvMesh& subMesh = subMeshes[lawI].subMesh();
        const labelList& patchMap = subMeshes[lawI].patchMap();

        forAll(subMeshGradDf.boundaryField(), patchI)
        {
            if (patchMap[patchI] == -1)
            {
                tensorField& patchGradDf =
                    subMeshGradDf.boundaryField()[patchI];
                const tensorField patchGradDif =
                    subMeshGradD.boundaryField()[patchI].patchInternalField();

                const vectorField& patchD = subMeshD.boundaryField()[patchI];
                const vectorField patchDif =
                    subMeshD.boundaryField()[patchI].patchInternalField();

                const vectorField n = subMesh.boundary()[patchI].nf();
                const vectorField delta = subMesh.boundary()[patchI].delta();
                const vectorField k = delta - n*(n & delta);
                const scalarField& deltaCoeffs =
                    subMesh.boundary()[patchI].deltaCoeffs();

                const vectorField correctedSnGrad =
                    (patchD - (patchDif + (k & patchGradDif)))*deltaCoeffs;

                patchGradDf += n*(correctedSnGrad - (n & patchGradDf));
            }
        }
    }
}


void Foam::mechanicalModel::clearOut() const
{
    deleteDemandDrivenData(volToPointPtr_);
    subMeshVolToPoint_.clear();
    subMeshSigma_.clear();
    subMeshSigmaf_.clear();
    subMeshD_.clear();
    subMeshDD_.clear();
    subMeshGradD_.clear();
    subMeshGradDf_.clear();
    subMeshGradDD_.clear();
    subMeshGradDDf_.clear();
    subMeshPointD_.clear();
    deleteDemandDrivenData(interfaceBaseFacesPtr_);
    interfaceShadowSubMeshID_.clear();
    interfaceShadowPatchID_.clear();
    interfaceShadowFaceID_.clear();
    deleteDemandDrivenData(impKfcorrPtr_);
    deleteDemandDrivenData(pointNumOfMaterialsPtr_);
    deleteDemandDrivenData(isolatedInterfacePointsPtr_);

    // Make sure to clear the subMeshes after (not before) clearing the subMesh
    // fields
    subMeshes_.clear();
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
    volToPointPtr_(),
    subMeshVolToPoint_(),
    subMeshSigma_(),
    subMeshSigmaf_(),
    subMeshD_(),
    subMeshDD_(),
    subMeshGradD_(),
    subMeshGradDf_(),
    subMeshGradDD_(),
    subMeshGradDDf_(),
    subMeshPointD_(),
    interfaceBaseFacesPtr_(NULL),
    interfaceShadowSubMeshID_(),
    interfaceShadowPatchID_(),
    interfaceShadowFaceID_(),
    impKfcorrPtr_(NULL),
    pointNumOfMaterialsPtr_(NULL),
    isolatedInterfacePointsPtr_(NULL)
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

            if (lookupOrDefault<Switch>("writeSubMeshes",  false))
            {
                Info<< "Writing subMeshes "
                    << subMeshes()[lawI].subMesh().name() << endl;
                subMeshes()[lawI].subMesh().setInstance("constant");
                subMeshes()[lawI].subMesh().write();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mechanicalModel::~mechanicalModel()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::mechanicalModel::mesh() const
{
    return mesh_;
}


const Foam::newLeastSquaresVolPointInterpolation&
Foam::mechanicalModel::volToPoint() const
{
    if (!volToPointPtr_)
    {
        makeVolToPoint();
    }

    return *volToPointPtr_;
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

        // Accumulated subMesh fields and then map to the base mesh
        PtrList<volScalarField> rhos(laws.size());

        forAll(laws, lawI)
        {
            rhos.set
            (
                lawI,
                new volScalarField(laws[lawI].rho())
            );
        }

        // Map subMesh fields to the base mesh
        mapSubMeshVolFields<scalar>(rhos, result);

        // Clear subMesh fields
        rhos.clear();

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

        // Accumulated subMesh fields and then map to the base mesh
        PtrList<volScalarField> impKs(laws.size());

        forAll(laws, lawI)
        {
            impKs.set
            (
                lawI,
                new volScalarField(laws[lawI].impK())
            );
        }

        // Map subMesh fields to the base mesh
        mapSubMeshVolFields<scalar>(impKs, result);

        // Clear subMesh fields
        impKs.clear();

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

        // Accumulated subMesh fields and then map to the base mesh
        PtrList<surfaceScalarField> impKfs(laws.size());

        forAll(laws, lawI)
        {
            impKfs.set
            (
                lawI,
                new surfaceScalarField(laws[lawI].impKf())
            );
        }

        // Map subMesh fields to the base mesh
        mapSubMeshSurfaceFields<scalar>(impKfs, result);

        // Clear subMesh fields
        impKfs.clear();

        return tresult;
    }
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::K() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        return laws[0].K();
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
                    "K",
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

        // Accumulated subMesh fields and then map to the base mesh
        PtrList<volScalarField> Ks(laws.size());

        forAll(laws, lawI)
        {
            Ks.set
            (
                lawI,
                new volScalarField(laws[lawI].K())
            );
        }

        // Map subMesh fields to the base mesh
        mapSubMeshVolFields<scalar>(Ks, result);

        // Clear subMesh fields
        Ks.clear();

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
        }

        // Map subMesh fields to the base field
        mapSubMeshVolFields<symmTensor>
        (
            subMeshSigma(), sigma
        );
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
        }

        // Map subMesh fields to the base field
        mapSubMeshSurfaceFields<symmTensor>
        (
            subMeshSigmaf(), sigma
        );
    }
}


void Foam::mechanicalModel::grad
(
    const volVectorField& D,
    volTensorField& gradD
)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (interfaceBaseFaces().size() == 0)
    {
        gradD = fvc::grad(D);
    }
    else
    {
        // Interpolate the base D to the subMesh D
        // If necessary, corrections are applied on bi-material interfaces
        interpolateDtoSubMeshD(D, true);

        // Accumulate data for all fields
        forAll(laws, lawI)
        {
            // Calculate gradient on subMesh
            // This will use the values at the interface
            if (D.name() == "D")
            {
                volTensorField& subMeshGradD = this->subMeshGradD()[lawI];
                subMeshGradD = fvc::grad(subMeshD()[lawI]);
            }
            else
            {
                volTensorField& subMeshGradDD = this->subMeshGradDD()[lawI];
                subMeshGradDD = fvc::grad(subMeshDD()[lawI]);
            }
        }

        // Map subMesh gradD to the base gradD
        if (D.name() == "D")
        {
            // Correct snGrad on the interface patch of subMeshGradD because the
            // default calculated boundaries disable non-orthogonal correction
            correctInterfaceSnGrad(subMeshD(), subMeshGradD());

            mapSubMeshVolFields<tensor>
            (
                subMeshGradD(), gradD
            );
        }
        else
        {
            // Correct snGrad on the interface patch of subMeshGradD because the
            // default calculated boundaries disable non-orthogonal correction
            correctInterfaceSnGrad(subMeshDD(), subMeshGradDD());

            mapSubMeshVolFields<tensor>
            (
                subMeshGradDD(), gradD
            );
        }

        // Correct boundary snGrad
        fv::gaussGrad<vector>
        (
            mesh()
        ).correctBoundaryConditions(D, gradD);
    }
}


void Foam::mechanicalModel::grad
(
    const volVectorField& D,
    const pointVectorField& pointD,
    volTensorField& gradD
)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (interfaceBaseFaces().size() == 0)
    {
        gradD = fvc::grad(D, pointD);
    }
    else
    {
        // Accumulate data for all fields
        forAll(laws, lawI)
        {
            // Map subMesh gradD to the base gradD
            if (D.name() == "D")
            {
                volTensorField& subMeshGradD = this->subMeshGradD()[lawI];
                subMeshGradD =
                    fvc::grad(subMeshD()[lawI], subMeshPointD()[lawI]);
            }
            else
            {
                volTensorField& subMeshGradDD = this->subMeshGradDD()[lawI];
                subMeshGradDD =
                    fvc::grad(subMeshDD()[lawI], subMeshPointD()[lawI]);
            }
        }

        if (D.name() == "D")
        {
            // Correct snGrad on the interface patch of subMeshGradDf because
            // the default calculated boundaries disable non-orthogonal
            // correction
            correctInterfaceSnGrad(subMeshD(), subMeshGradD());

            // Map subMesh gradD fields to the base gradD field
            mapSubMeshVolFields<tensor>
            (
                this->subMeshGradD(), gradD
            );
        }
        else
        {
            // Correct snGrad on the interface patch of subMeshGradDf because
            // the default calculated boundaries disable non-orthogonal
            // correction
            correctInterfaceSnGrad(subMeshDD(), subMeshGradDD());

            // Map subMesh gradD fields to the base gradD field
            mapSubMeshVolFields<tensor>
            (
                this->subMeshGradDD(), gradD
            );
        }

        // Correct boundary snGrad
        fv::gaussGrad<vector>
        (
            mesh()
        ).correctBoundaryConditions(D, gradD);
    }
}


void Foam::mechanicalModel::grad
(
    const volVectorField& D,
    const pointVectorField& pointD,
    surfaceTensorField& gradDf
)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (interfaceBaseFaces().size() == 0)
    {
        gradDf = fvc::fGrad(D, pointD);
    }
    else
    {
        // Accumulate data for all fields
        forAll(laws, lawI)
        {
            if (D.name() == "D")
            {
                surfaceTensorField& subMeshGradDf = this->subMeshGradDf()[lawI];
                subMeshGradDf =
                    fvc::fGrad(subMeshD()[lawI], subMeshPointD()[lawI]);
            }
            else
            {
                surfaceTensorField& subMeshGradDDf =
                    this->subMeshGradDDf()[lawI];
                subMeshGradDDf =
                    fvc::fGrad(subMeshDD()[lawI], subMeshPointD()[lawI]);
            }
        }

        if (D.name() == "D")
        {
            // Correct snGrad on the interface patch of subMeshGradDf because
            // the default calculated boundaries disable non-orthogonal
            // correction
            correctInterfaceSnGradf
            (
                subMeshD(), subMeshGradDf(), subMeshGradD()
            );

            // Map subMesh gradDf fields to the base gradDf field
            mapSubMeshSurfaceFields<tensor>
            (
                subMeshGradDf(), gradDf
            );
        }
        else
        {
            // Correct snGrad on the interface patch of subMeshGradDf because
            // the default calculated boundaries disable non-orthogonal
            // correction
            correctInterfaceSnGradf
            (
                subMeshDD(), subMeshGradDDf(), subMeshGradDD()
            );

            // Map subMesh gradDf fields to the base gradDf field
            mapSubMeshSurfaceFields<tensor>
            (
                subMeshGradDDf(), gradDf
            );
        }

        // Replace normal component
        // If we don't do this then we don't get convergence in many cases
        const surfaceVectorField n = mesh().Sf()/mesh().magSf();
        gradDf += n*fvc::snGrad(D) - (sqr(n) & gradDf);
    }
}


void Foam::mechanicalModel::interpolate
(
    const volVectorField& D,
    pointVectorField& pointD,
    const bool useVolFieldSigma
)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (interfaceBaseFaces().size() == 0)
    {
        volToPoint().interpolate(D, pointD);
    }
    else
    {
        // Interpolate the base D to the subMesh D
        // If necessary, corrections are applied on bi-material interfaces
        interpolateDtoSubMeshD(D, useVolFieldSigma);

        // Accumulate data for all fields
        forAll(laws, lawI)
        {
            if (D.name() == "D")
            {
                // Interpolate the subMeshD to the subMeshPointD
                subMeshVolToPoint()[lawI].interpolate
                (
                    subMeshD()[lawI],
                    subMeshPointD()[lawI]
                );
            }
            else
            {
                // Interpolate the subMeshD to the subMeshPointD
                subMeshVolToPoint()[lawI].interpolate
                (
                    subMeshDD()[lawI],
                    subMeshPointD()[lawI]
                );
            }
        }

        // Map subMesh pointD fields back to the base pointD field
        mapSubMeshPointFields<vector>
        (
            subMeshPointD(), pointD
        );
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
