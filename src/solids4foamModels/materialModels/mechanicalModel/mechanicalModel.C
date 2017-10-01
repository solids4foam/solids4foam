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
#include "twoDPointCorrector.H"

#include "fixedGradientFvPatchFields.H"

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

    // The subMeshD field can represent D or DD
    word Dname = "D";
    if (mesh().foundObject<volVectorField>("DD"))
    {
        Dname = "DD";
    }

    forAll(laws, lawI)
    {
        subMeshD_.set
        (
            lawI,
            new volVectorField
            (
                IOobject
                (
                    Dname,
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

    // The subMeshD field can represent D or DD
    word gradDname = "grad(D)";
    if (incremental_)
    {
        gradDname = "grad(DD)";
    }

    forAll(laws, lawI)
    {
        subMeshGradD_.set
        (
            lawI,
            new volTensorField
            (
                IOobject
                (
                    gradDname,
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

    // The subMeshD field can represent D or DD
    word gradDname = "grad(D)f";
    if (incremental_)
    {
        gradDname = "grad(DD)f";
    }

    forAll(laws, lawI)
    {
        subMeshGradDf_.set
        (
            lawI,
            new surfaceTensorField
            (
                IOobject
                (
                    gradDname,
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

    // The subMeshD field can represent D or DD
    word pointDname = "pointD";
    if (mesh().foundObject<volVectorField>("DD"))
    {
        pointDname = "pointDD";
    }

    forAll(laws, lawI)
    {
        subMeshPointD_.set
        (
            lawI,
            new pointVectorField
            (
                IOobject
                (
                    pointDname,
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
                    const label baseFaceID = faceMap[start + faceI];

                    if (mesh().isInternalFace(baseFaceID))
                    {
                        // Check if the face has been set in the baseMesh lists
                        if (baseMeshShadSubMeshID[baseFaceID] == -1)
                        {
                            // Store the local IDs in the baseMesh lists
                            baseMeshShadSubMeshID[baseFaceID] = lawI;
                            baseMeshShadPatchID[baseFaceID] = patchI;
                            baseMeshShadFaceID[baseFaceID] = faceI;
                        }
                        else
                        {
                            // Store the shadow values in the local lists
                            shadSubMeshID[faceI] =
                                baseMeshShadSubMeshID[baseFaceID];
                            shadPatchID[faceI] =
                                baseMeshShadPatchID[baseFaceID];
                            shadFaceID[faceI] =
                                baseMeshShadFaceID[baseFaceID];

                            // Update the shadow lists with the local values
                            const label shadSubMeshID =
                                baseMeshShadSubMeshID[baseFaceID];
                            //const label shadPatchID =
                            //    baseMeshShadPatchID[baseFaceID];
                            const label shadFaceID =
                                baseMeshShadFaceID[baseFaceID];

                            interfaceShadowSubMeshID_
                            [
                                shadSubMeshID
                            ][shadFaceID] = lawI;
                            interfaceShadowPatchID_[shadSubMeshID][shadFaceID] =
                                patchI;
                            interfaceShadowFaceID_[shadSubMeshID][shadFaceID] =
                                faceI;
                        }
                    }
                    else
                    {
                        // Shadow IDs not set faces that are on a processor
                        // patch in the base mesh
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
        // impKfcorr to zero on bi-material interface faces

        surfaceScalarField& impKfcorr = *impKfcorrPtr_;
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
                        const label baseFaceID = faceMap[start + faceI];

                        if (mesh().isInternalFace(baseFaceID))
                        {
                            impKfcorrI[baseFaceID] = 0.0;
                        }
                        else
                        {
                            // Face is on a coupled patch
                            const label patchID =
                                mesh().boundaryMesh().whichPatch(baseFaceID);

                            const label basePatchStart =
                                mesh().boundaryMesh()[patchID].start();

                            impKfcorr.boundaryField()
                            [
                                patchID
                            ][baseFaceID - basePatchStart] = 0.0;
                        }
                    }
                }
            }
        }

        impKfcorrPtr_->correctBoundaryConditions();
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

        // Sync coupled boundaries
        materials.correctBoundaryConditions();

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

    if (!biMaterialInterfaceActive())
    {
        // No need for any corrections if there are no bi-material interfaces
        forAll(subMeshes, lawI)
        {
            subMeshD()[lawI] = subMeshes[lawI].interpolate(D);
        }

        return;
    }

    const fvMesh& mesh = this->mesh();

    // Update the interface shadow sigma fields
    // This can require parallel communication
    updateInterfaceShadowSigma(useVolFieldSigma);

    forAll(subMeshes, lawI)
    {
        volVectorField& subMeshD = this->subMeshD()[lawI];

        const fvMesh& subMesh = subMeshes[lawI].subMesh();

        const labelList& faceMap = subMeshes[lawI].faceMap();
        const labelList& patchMap = subMeshes[lawI].patchMap();
        const labelList& cellMap = subMeshes[lawI].cellMap();

        // Store interface field as it is overwritten with the interpolated
        // value by the interpolate function
        vectorField DinterfacePrev(0);
        forAll(subMeshD.boundaryField(), patchI)
        {
            if (patchMap[patchI] == -1)
            {
                DinterfacePrev = subMeshD.boundaryField()[patchI];
            }
        }

        // Map the base displacement field to the subMesh; this overwrites
        // the interface with the interpolated values
        subMeshD = subMeshes[lawI].interpolate(D);

        // Check if a large strain procedure is being used, if so we must
        // calculate the deformed normals
        // If the deformation gradient field 'F' is found, we will assume it is
        // a large/finite strain procedure
        // For finite strain procedures, we will look up the deformation
        // gradient: relative deformation gradient for updated Lagrangian
        // approaches and the total deformation gradient for total
        // approaches
        // What about uns approaches? I may need to include surfaceField options
        // here
        const bool useDeformedNormals = mesh.foundObject<volTensorField>("F");
        const volTensorField* FinvPtr = NULL;
        const volScalarField* JPtr = NULL;
        if (useDeformedNormals)
        {
            if (mesh.foundObject<volTensorField>("relF"))
            {
                // Updated Lagrangian approach: use the inverse of the relative
                // deformation gradient
                FinvPtr = &(mesh.lookupObject<volTensorField>("relFinv"));
                JPtr = &(mesh.lookupObject<volScalarField>("relJ"));
            }
            else
            {
                // Total Lagrangian approach: use the inverse of the total
                // deformation gradient
                FinvPtr = &(mesh.lookupObject<volTensorField>("Finv"));
                JPtr = &(mesh.lookupObject<volScalarField>("J"));
            }
        }

        // Calculate the new correction to the interface valaues
        forAll(subMeshD.boundaryField(), patchI)
        {
            if (patchMap[patchI] == -1)
            {
                // Interface displacement
                vectorField& Dinterface = subMeshD.boundaryField()[patchI];

                // Patch start face index
                const label start = subMesh.boundaryMesh()[patchI].start();

                // Base mesh owner cells
                const unallocLabelList& baseOwn = mesh.owner();

                // Base mesh neighbour cells
                const unallocLabelList& baseNei = mesh.neighbour();

                // Base mesh face interpolation weights
                const surfaceScalarField& baseWeights = mesh.weights();
                const scalarField& baseWeightsI = baseWeights.internalField();

                // Base mesh cell centres
                const volVectorField& baseC = mesh.C();
                const vectorField& baseCI = baseC.internalField();

                // Base mesh face area vectors
                const surfaceVectorField& baseSf = mesh.Sf();
                const vectorField& baseSfI = baseSf.internalField();

                // Base mesh face area vector magnitudes
                const surfaceScalarField& baseMagSf = mesh.magSf();
                const scalarField& baseMagSfI = baseMagSf.internalField();

                // Implicit stiffness field
                const volScalarField& K =
                    mesh.lookupObject<volScalarField>("impK");
                const scalarField& KI = K.internalField();

                // Interface face centres
                const vectorField& patchCf = subMesh.boundary()[patchI].Cf();

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

                // Assemble the shadow sigma field: this is the stress
                // calculated from the other side of the interface (in the
                // subMesh on the other side)
                const symmTensorField& interfaceShadSigma =
                    interfaceShadowSigma()[lawI];

                // Calculate the interface displacements
                forAll(Dinterface, faceI)
                {
                    // Base mesh face index
                    const label baseFaceID = faceMap[start + faceI];

                    if (mesh.isInternalFace(baseFaceID))
                    {
                        // Base mesh owner cell index
                        const label baseOwnID = baseOwn[baseFaceID];

                        // Base mesh neighbour cell index
                        const label baseNeiID = baseNei[baseFaceID];

                        // Base mesh face interpolation weight
                        const scalar baseW = baseWeightsI[baseFaceID];

                        // Interface unit normal (on base mesh); this may be in
                        // the opposite direction to the subMesh normal
                        vector n = baseSfI[baseFaceID]/baseMagSfI[baseFaceID];
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
                            // Lagrangian approaches, it is the total
                            // deformation gradient
                            n = J*Finv.T() & n;
                            n /= mag(n);
                        }

                        // In surfaceInterpolation.C the deltaCoeffs are
                        // calculated as:
                        // 1.0/max(unitArea & delta, 0.05*mag(delta));

                        // Normal distance from the interface to the cell-centre
                        // on side-a
                        scalar da =
                            mag(n & (patchCf[faceI] - baseCI[baseOwnID]));

                        // Normal distance from the interface to the cell-centre
                        // on side-b
                        scalar db =
                            mag(n & (baseCI[baseNeiID] - patchCf[faceI]));

                        // Lookup the stiffness wither side of the interface
                        scalar Ka = KI[baseOwnID];
                        scalar Kb = KI[baseNeiID];

                        // The base own cell should be side a
                        if (baseOwnID != cellMap[faceCells[faceI]])
                        {
                            n = -n;
                            Swap(Ka, Kb);
                            Swap(da, db);
                        }

                        // Calculate the traction at side-a and side-b
                        const vector tractiona = n & sigmaPatch[faceI];
                        const vector tractionb = n & interfaceShadSigma[faceI];

                        // Calculate the displacement at the interface
                        Dinterface[faceI] =
                            DinterfacePrev[faceI]
                          + (da*db/(db*Ka + da*Kb))*(tractionb - tractiona);
                    }
                    else
                    {
                        // These are faces that are on a processor patch in the
                        // baseMesh but were placed in the oldInternalFaces
                        // patch in the subMesh: these faces are on a
                        // bi-material interface

                        // nei: sb => stored on the patch
                        // nei: db => also stored on the patch: to check
                        // So if I find that the baseFaceID is on a proc patch,
                        // then I can lookup the patch ID and directly lookup
                        // sb and db: nothing else is required

                        // Base mesh patch ID
                        const label basePatchID =
                            mesh.boundaryMesh().whichPatch(baseFaceID);

                        // Base patch start
                        const label basePatchStart =
                            mesh.boundaryMesh()[basePatchID].start();

                        // Base mesh patch local face ID
                        const label baseLocalFaceID =
                            baseFaceID - basePatchStart;

                        // Base mesh owner cell index
                        //const label baseOwnID = baseOwn[baseFaceID];
                        const label baseOwnID =
                            mesh.boundaryMesh()
                            [
                                basePatchID
                            ].faceCells()[baseLocalFaceID];

                        // Base mesh face interpolation weight
                        //const scalar baseW = baseWeightsI[baseFaceID];
                        const scalar baseW =
                            baseWeights.boundaryField()
                            [
                                basePatchID
                            ][baseLocalFaceID];

                        // Interface unit normal (on base mesh); this may be in
                        // the opposite direction to the subMesh normal
                        //vector n = baseSf[baseFaceID]/baseMagSf[baseFaceID];
                        vector n =
                            baseSf.boundaryField()[basePatchID][baseLocalFaceID]
                           /baseMagSf.boundaryField()
                            [
                                basePatchID
                            ][baseLocalFaceID];

                        if (useDeformedNormals)
                        {
                            // Interpolate Finv and J to the face
                            // to-do
                            const tensorField& FinvI =
                                FinvPtr->internalField();
                            const tensor Finv =
                                baseW*FinvI[baseOwnID]
                              + (1.0 - baseW)
                               *FinvPtr->boundaryField()
                                [
                                    basePatchID
                                ][baseLocalFaceID];

                            const scalarField& JI = JPtr->internalField();
                            const scalar J =
                                baseW*JI[baseOwnID]
                              + (1.0 - baseW)
                               *JPtr->boundaryField()
                                [
                                    basePatchID
                                ][baseLocalFaceID];

                            // Nanson's formula
                            // Note: for updated Lagrangian approach, F is the
                            // relative deformation gradient, whereas for total
                            // Lagrangian approaches, it is the total
                            // deformation gradient
                            n = J*Finv.T() & n;
                            n /= mag(n);
                        }

                        // Normal distance from the interface to the cell-centre
                        // on side-a
                        const scalar da =
                            mag(n & (patchCf[faceI] - baseCI[baseOwnID]));

                        // Normal distance from the interface to the cell-centre
                        // on side-b
                        //const scalar db =
                        //    mag(n & (baseCI[baseNeiID] - Cf[faceI]));
                        // Note: processor patches store the patchNeighbourField
                        // directly on the patch, so the patch value will
                        // correspond to the patchNeighbourField value
                        const scalar db =
                            mag
                            (
                                n
                              & (
                                  baseC.boundaryField()
                                  [
                                      basePatchID
                                  ][baseLocalFaceID]
                                - patchCf[faceI]
                              )
                            );

                        // Lookup the stiffness at either side
                        const scalar Ka = KI[baseOwnID];
                        const scalar Kb =
                            K.boundaryField()[basePatchID][baseLocalFaceID];

                        // Calculate the traction at side-a and side-b
                        const vector tractiona = n & sigmaPatch[faceI];
                        const vector tractionb = n & interfaceShadSigma[faceI];

                        // Calculate the displacement at the interface
                        Dinterface[faceI] =
                            DinterfacePrev[faceI]
                          + (da*db/(db*Ka + da*Kb))*(tractionb - tractiona);
                    }
                }
            }
            else
            {
                // These are other patches including processor patches, but it
                // seems that interface processor patch faces are placed in the
                // oldInternalFaces patch, so we check for them above.
                // For now, no need to do anything here.
            }
        }
    }
}


void Foam::mechanicalModel::correctBoundarySnGrad
(
    PtrList<volVectorField>& subMeshDList,
    PtrList<volTensorField>& subMeshGradDList
)
{
    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

    forAll(subMeshes, lawI)
    {
        const volVectorField& subMeshD = subMeshDList[lawI];
        volTensorField& subMeshGradD = subMeshGradDList[lawI];
        const fvMesh& subMesh = subMeshes[lawI].subMesh();

        forAll(subMeshGradD.boundaryField(), patchI)
        {
            tensorField& patchGradD = subMeshGradD.boundaryField()[patchI];
            const tensorField patchGradDif =
                subMeshGradD.boundaryField()[patchI].patchInternalField();

            const vectorField& patchD = subMeshD.boundaryField()[patchI];
            const vectorField patchDif =
                subMeshD.boundaryField()[patchI].patchInternalField();

            const vectorField n = subMesh.boundary()[patchI].nf();
            const vectorField delta = subMesh.boundary()[patchI].delta();
            const vectorField k = delta - (sqr(n) & delta);
            const scalarField& deltaCoeffs =
                subMesh.boundary()[patchI].deltaCoeffs();

            const vectorField correctedSnGrad =
                (patchD - (patchDif + (k & patchGradDif)))*deltaCoeffs;

            patchGradD += n*(correctedSnGrad - (n & patchGradD));
        }
    }
}


void Foam::mechanicalModel::correctBoundarySnGradf
(
    PtrList<volVectorField>& subMeshDList,
    PtrList<surfaceTensorField>& subMeshGradDfList,
    PtrList<volTensorField>& subMeshGradDList
)
{
    const PtrList<newFvMeshSubset>& subMeshes = this->subMeshes();

    forAll(subMeshes, lawI)
    {
        const volVectorField& subMeshD = subMeshDList[lawI];
        surfaceTensorField& subMeshGradDf = subMeshGradDfList[lawI];
        const volTensorField& subMeshGradD = subMeshGradDList[lawI];
        const fvMesh& subMesh = subMeshes[lawI].subMesh();

        forAll(subMeshGradDf.boundaryField(), patchI)
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
            const vectorField k = delta - (sqr(n) & delta);
            const scalarField& deltaCoeffs =
                subMesh.boundary()[patchI].deltaCoeffs();

            const vectorField correctedSnGrad =
                (patchD - (patchDif + (k & patchGradDif)))*deltaCoeffs;

            patchGradDf += n*(correctedSnGrad - (n & patchGradDf));
        }
    }
}


const Foam::PtrList<Foam::symmTensorField>&
Foam::mechanicalModel::interfaceShadowSigma() const
{
    if (interfaceShadowSigma_.empty())
    {
        makeInterfaceShadowSigma();
    }

    return interfaceShadowSigma_;
}


void Foam::mechanicalModel::makeInterfaceShadowSigma() const
{
    if (!interfaceShadowSigma_.empty())
    {
        FatalErrorIn
        (
            "void Foam::mechanicalModel::makeInterfaceShadowSigma() const"
        ) << "pointer already set" << abort(FatalError);
    }

    interfaceShadowSigma_.setSize(subMeshes().size());

    // Set values for each subMesh
    forAll(subMeshes(), subMeshI)
    {
        const newFvMeshSubset& subsetMesh = subMeshes()[subMeshI];
        const fvMesh& subMesh = subsetMesh.subMesh();
        const labelList& patchMap = subsetMesh.patchMap();

        // Find the interface patch for the current subMesh
        // we should store this!

        label patchID = -1;

        forAll(subMesh.boundaryMesh(), pI)
        {
            if (patchMap[pI] == -1)
            {
                patchID = pI;
                break;
            }
        }

        if (patchID == -1)
        {
            FatalErrorIn
            (
                "void Foam::mechanicalModel::makeInterfaceShadowSigmaGradD()"
                "const"
            )   << "Interface patch not found!" << abort(FatalError);
        }

        interfaceShadowSigma_.set
        (
            subMeshI,
            new symmTensorField
            (
                subMesh.boundaryMesh()[patchID].size(),
                symmTensor::zero
            )
        );
    }
}


void Foam::mechanicalModel::updateInterfaceShadowSigma
(
    const bool useVolFieldSigma
)
{
    if (interfaceShadowSigma_.empty())
    {
        makeInterfaceShadowSigma();
    }

    // Field used for syncing the processor patch values
    volSymmTensorField baseSigmaForSyncing
    (
        IOobject
        (
            "baseSigmaForSyncing",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    // Set values for each subMesh
    forAll(subMeshes(), subMeshI)
    {
        const newFvMeshSubset& subsetMesh = subMeshes()[subMeshI];
        const fvMesh& subMesh = subsetMesh.subMesh();
        const labelList& patchMap = subsetMesh.patchMap();
        const labelList& faceMap = subsetMesh.faceMap();

        // Find the interface patch for the current subMesh
        // we should store this!

        label patchID = -1;

        forAll(subMesh.boundaryMesh(), pI)
        {
            if (patchMap[pI] == -1)
            {
                patchID = pI;
                break;
            }
        }

        if (patchID == -1)
        {
            FatalErrorIn
            (
                "void Foam::updateInterfaceShadowSigma\n"
                "(\n"
                "    const bool useVolFieldSigma\n"
                ")"
            )   << "Interface patch not found!" << abort(FatalError);
        }

        symmTensorField& resultSigma = interfaceShadowSigma_[subMeshI];

        const labelList& interfaceShadowSubMeshID =
            this->interfaceShadowSubMeshID()[subMeshI];
        const labelList& interfaceShadowPatchID =
            this->interfaceShadowPatchID()[subMeshI];
        const labelList& interfaceShadowFaceID =
            this->interfaceShadowFaceID()[subMeshI];

        // Assemble the shadow stress for each face
        forAll(resultSigma, faceI)
        {
            // Check if the face is not on a processor
            if (interfaceShadowSubMeshID[faceI] != -1)
            {
                // ID of the subMesh on the other side of interface
                const label shadowSubMeshID = interfaceShadowSubMeshID[faceI];

                // ID of the subMesh patch on the other side of interface
                const label shadowPatchID = interfaceShadowPatchID[faceI];

                // ID of the subMesh patch face on the other side of interface
                const label shadowFaceID = interfaceShadowFaceID[faceI];

                // Stress calculated at the other side of the interface
                if (useVolFieldSigma)
                {
                    resultSigma[faceI] =
                        subMeshSigma()
                        [
                            shadowSubMeshID
                        ].boundaryField()[shadowPatchID][shadowFaceID];
                }
                else
                {
                    resultSigma[faceI] =
                        subMeshSigmaf()
                        [
                            shadowSubMeshID
                        ].boundaryField()[shadowPatchID][shadowFaceID];
                }
            }
            else
            {
                // Base face is on a processor boundary

                // Local patch start
                const label start = subMesh.boundaryMesh()[patchID].start();

                // Base mesh face ID
                const label baseFaceID = faceMap[start + faceI];

                // Base mesh patch ID
                const label basePatchID =
                    mesh().boundaryMesh().whichPatch(baseFaceID);

                // Base patch start
                const label basePatchStart =
                    mesh().boundaryMesh()[basePatchID].start();

                // Base mesh patch local face ID
                const label baseLocalFaceID = baseFaceID - basePatchStart;

                // Base mesh patch faceCells
                const unallocLabelList& faceCells =
                    mesh().boundaryMesh()[basePatchID].faceCells();

                // Store local stress on the baseMesh proc patch in the patch
                // internal field
                if (useVolFieldSigma)
                {
                    baseSigmaForSyncing.internalField()
                        [
                            faceCells[baseLocalFaceID]
                        ]
                      = subMeshSigma()
                        [
                            subMeshI
                        ].boundaryField()[patchID][faceI];
                }
                else
                {
                    baseSigmaForSyncing.internalField()
                        [
                            faceCells[baseLocalFaceID]
                        ] = subMeshSigmaf()
                        [
                            subMeshI
                        ].boundaryField()[patchID][faceI];
                }
            }
        }
    }

    // Sync base mesh processor patches
    // This will pass the patch internal field and store it on the neighbour
    // patch
    baseSigmaForSyncing.correctBoundaryConditions();

    // Assemble processor values that have been synced
    forAll(subMeshes(), subMeshI)
    {
        const newFvMeshSubset& subsetMesh = subMeshes()[subMeshI];
        const fvMesh& subMesh = subsetMesh.subMesh();
        const labelList& patchMap = subsetMesh.patchMap();
        const labelList& faceMap = subsetMesh.faceMap();

        // Find the interface patch for the current subMesh
        // we should store this!

        label patchID = -1;

        forAll(subMesh.boundaryMesh(), pI)
        {
            if (patchMap[pI] == -1)
            {
                patchID = pI;
                break;
            }
        }

        if (patchID == -1)
        {
            FatalErrorIn
            (
                "void Foam::updateInterfaceShadowSigma\n"
                "(\n"
                "    const bool useVolFieldSigma\n"
                ")"
            )   << "Interface patch not found!" << abort(FatalError);
        }

        symmTensorField& resultSigma = interfaceShadowSigma_[subMeshI];

        const labelList& interfaceShadowSubMeshID =
            this->interfaceShadowSubMeshID()[subMeshI];

        forAll(resultSigma, faceI)
        {
            // Check if the face is on a processor
            if (interfaceShadowSubMeshID[faceI] == -1)
            {
                // The base mesh field will now have the patchNeighbourField
                // values stored on the patch

                // Local patch start
                const label start = subMesh.boundaryMesh()[patchID].start();

                // Base mesh face ID
                const label baseFaceID = faceMap[start + faceI];

                // Base mesh patch ID
                const label basePatchID =
                    mesh().boundaryMesh().whichPatch(baseFaceID);

                // Base patch start
                const label basePatchStart =
                    mesh().boundaryMesh()[basePatchID].start();

                // Base mesh patch local face ID
                const label baseLocalFaceID = baseFaceID - basePatchStart;

                // Copy patch neighbour field values into the result field
                resultSigma[faceI] =
                    baseSigmaForSyncing.boundaryField()
                    [
                        basePatchID
                    ][baseLocalFaceID];
            }
        }
    }
}


bool Foam::mechanicalModel::biMaterialInterfaceActive() const
{
    if (!biMaterialInterfaceActivePtr_)
    {
        calcBiMaterialInterfaceActive();
    }

    return *biMaterialInterfaceActivePtr_;
}


void Foam::mechanicalModel::calcBiMaterialInterfaceActive() const
{
    if (biMaterialInterfaceActivePtr_)
    {
        FatalErrorIn
        (
            "void Foam::mechanicalModel::calcBiMaterialInterfaceActive() const"
        ) << "pointer already set" << abort(FatalError);
    }

    biMaterialInterfaceActivePtr_ =
        new bool(returnReduce(interfaceBaseFaces().size(), maxOp<int>()));
}


void Foam::mechanicalModel::clearOut()
{
    deleteDemandDrivenData(volToPointPtr_);
    subMeshVolToPoint_.clear();
    subMeshSigma_.clear();
    subMeshSigmaf_.clear();
    subMeshD_.clear();
    subMeshGradD_.clear();
    subMeshGradDf_.clear();
    subMeshPointD_.clear();
    deleteDemandDrivenData(biMaterialInterfaceActivePtr_);
    deleteDemandDrivenData(interfaceBaseFacesPtr_);
    interfaceShadowSubMeshID_.clear();
    interfaceShadowPatchID_.clear();
    interfaceShadowFaceID_.clear();
    interfaceShadowSigma_.clear();
    deleteDemandDrivenData(impKfcorrPtr_);
    deleteDemandDrivenData(pointNumOfMaterialsPtr_);
    deleteDemandDrivenData(isolatedInterfacePointsPtr_);

    // Clear the list of mechanical laws
    // Note: we should do this before clearing the subMeshes, as the mechanical
    // laws can store geometricFields that must be deleted before deleting
    // mesh
    PtrList<mechanicalLaw>::clear();

    // Make sure to clear the subMeshes after (not before) clearing the subMesh
    // fields
    subMeshes_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalModel::mechanicalModel
(
    const fvMesh& mesh,
    const nonLinearGeometry::nonLinearType& nonLinGeom,
    const bool incremental
)
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
    incremental_(incremental),
    cellZoneNames_(),
    subMeshes_(),
    volToPointPtr_(),
    subMeshVolToPoint_(),
    subMeshSigma_(),
    subMeshSigmaf_(),
    subMeshD_(),
    subMeshGradD_(),
    subMeshGradDf_(),
    subMeshPointD_(),
    biMaterialInterfaceActivePtr_(NULL),
    interfaceBaseFacesPtr_(NULL),
    interfaceShadowSubMeshID_(),
    interfaceShadowPatchID_(),
    interfaceShadowFaceID_(),
    interfaceShadowSigma_(),
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
        if (nonLinGeom == nonLinearGeometry::LINEAR_GEOMETRY)
        {
            laws.set
            (
                0,
                mechanicalLaw::NewLinGeomMechLaw
                (
                    lawEntries[0].keyword(),
                    mesh,
                    lawEntries[0].dict(),
                    nonLinGeom
                )
            );
        }
        else if
        (
            nonLinGeom == nonLinearGeometry::UPDATED_LAGRANGIAN
         || nonLinGeom == nonLinearGeometry::TOTAL_LAGRANGIAN
        )
        {
            laws.set
            (
                0,
                mechanicalLaw::NewNonLinGeomMechLaw
                (
                    lawEntries[0].keyword(),
                    mesh,
                    lawEntries[0].dict(),
                    nonLinGeom
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "Foam::mechanicalModel::mechanicalModel\n"
                "(\n"
                "    const fvMesh& mesh,\n"
                "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
                ")"
            )   << "It is not clear what type of mechanical law should be "
                << "created for a solidModel with nonLinGeom = " << nonLinGeom
                << abort(FatalError);
        }
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
            if (nonLinGeom == nonLinearGeometry::LINEAR_GEOMETRY)
            {
                laws.set
                (
                    lawI,
                    mechanicalLaw::NewLinGeomMechLaw
                    (
                        lawEntries[lawI].keyword(),
                        subMeshes()[lawI].subMesh(),
                        lawEntries[lawI].dict(),
                        nonLinGeom
                    )
                );
            }
            else if
            (
                nonLinGeom == nonLinearGeometry::UPDATED_LAGRANGIAN
             || nonLinGeom == nonLinearGeometry::TOTAL_LAGRANGIAN
            )
            {
                laws.set
                (
                    lawI,
                    mechanicalLaw::NewNonLinGeomMechLaw
                    (
                        lawEntries[lawI].keyword(),
                        subMeshes()[lawI].subMesh(),
                        lawEntries[lawI].dict(),
                        nonLinGeom
                    )
                );
            }
            else
            {
                FatalErrorIn
                (
                    "Foam::mechanicalModel::mechanicalModel\n"
                    "(\n"
                    "    const fvMesh& mesh,\n"
                    "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
                    ")"
                )   << "It is not clear what type of mechanical law should be "
                    << "created for a solidModel with nonLinGeom = "
                    << nonLinGeom << abort(FatalError);
            }

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
                    "rhoLaw",
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
    // Check the interpolaiton scheme
    // Harmonic interpolation should typically be used unless you know what you
    // are doing!

    const volScalarField impK = this->impK();
    const word interpName = "interpolate(" + impK.name() + ')';

    if (PtrList<mechanicalLaw>::size() > 1)
    {
        if
        (
            word(mesh().schemesDict().interpolationScheme(interpName))
         != "harmonic"
        )
        {
            WarningIn
            (
                "Foam::tmp<Foam::surfaceScalarField> "
                "Foam::mechanicalModel::impKf() const"
            )   << "The interpolation scheme for " << interpName << " is "
                << " currently set to "
                << word(mesh().schemesDict().interpolationScheme(interpName))
                << "; however," << nl
                << "    harmonic typically gives the best convergence,"
                << " particularly for multi-material cases!" << endl;
        }
    }

    return fvc::interpolate(impK, interpName);
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

    if (laws.size() == 1)
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
            volTensorField& subMeshGradD = this->subMeshGradD()[lawI];
            subMeshGradD = fvc::grad(subMeshD()[lawI]);
        }

        // Map subMesh gradD to the base gradD
        correctBoundarySnGrad(subMeshD(), subMeshGradD());

        mapSubMeshVolFields<tensor>(subMeshGradD(), gradD);

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

    if (laws.size() == 1)
    {
        gradD = fvc::grad(D, pointD);
    }
    else
    {
        // Calculate subMesh gradient fields
        forAll(laws, lawI)
        {
            volTensorField& subMeshGradD = this->subMeshGradD()[lawI];
            subMeshGradD = fvc::grad(subMeshD()[lawI], subMeshPointD()[lawI]);
        }

        // Correct snGrad on boundaries because the default calculated
        // boundaries disable non-orthogonal correction
        correctBoundarySnGrad(subMeshD(), subMeshGradD());

        // Map subMesh gradD fields to the base gradD field
        mapSubMeshVolFields<tensor>(subMeshGradD(), gradD);

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

    if (laws.size() == 1)
    {
        gradDf = fvc::fGrad(D, pointD);
    }
    else
    {
        // Calculate subMesh gradient fields
        forAll(laws, lawI)
        {
            surfaceTensorField& subMeshGradDf = this->subMeshGradDf()[lawI];
            subMeshGradDf = fvc::fGrad(subMeshD()[lawI], subMeshPointD()[lawI]);
        }

        // Correct snGrad on boundaries because the default calculated
        // boundaries disable non-orthogonal correction
        correctBoundarySnGradf(subMeshD(), subMeshGradDf(), subMeshGradD());

        // Map subMesh gradDf fields to the base gradDf field
        mapSubMeshSurfaceFields<tensor>(subMeshGradDf(), gradDf);

        // Replace normal component
        // If we don't do this then we don't get convergence in many cases
        const surfaceVectorField n = mesh().Sf()/mesh().magSf();
        gradDf += n*fvc::snGrad(D) - (sqr(n) & gradDf);
    }
}


void Foam::mechanicalModel::grad
(
    const volVectorField& D,
    const pointVectorField& pointD,
    volTensorField& gradD,
    surfaceTensorField& gradDf
)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        gradD = fvc::grad(D, pointD);
        gradDf = fvc::fGrad(D, pointD);
    }
    else
    {
        // Calculate subMesh gradient fields
        forAll(laws, lawI)
        {
            volTensorField& subMeshGradD = this->subMeshGradD()[lawI];
            subMeshGradD = fvc::grad(subMeshD()[lawI], subMeshPointD()[lawI]);

            surfaceTensorField& subMeshGradDf = this->subMeshGradDf()[lawI];
            subMeshGradDf = fvc::fGrad(subMeshD()[lawI], subMeshPointD()[lawI]);
        }

        // Correct snGrad on boundaries because the default calculated
        // boundaries disable non-orthogonal correction
        correctBoundarySnGrad(subMeshD(), subMeshGradD());
        correctBoundarySnGradf(subMeshD(), subMeshGradDf(), subMeshGradD());

        // Map subMesh fields to the base field

        mapSubMeshVolFields<tensor>(subMeshGradD(), gradD);
        mapSubMeshSurfaceFields<tensor>(subMeshGradDf(), gradDf);

        // Correct boundary snGrad of gradD
        fv::gaussGrad<vector>
        (
            mesh()
        ).correctBoundaryConditions(D, gradD);

        // Correct snGrad component of gradDf
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

    if (laws.size() == 1)
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
            // Interpolate the subMeshD to the subMeshPointD
            subMeshVolToPoint()[lawI].interpolate
            (
                subMeshD()[lawI],
                subMeshPointD()[lawI]
            );
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
    const volTensorField& gradD,
    const surfaceScalarField& gamma
) const
{
    // Mathematically "div(grad(phi))" is equivalent to "laplacian(phi)";
    // however, numerically "div(grad(phi))" uses a larger stencil than the
    // "laplacian(phi)"; the difference between these two approximations is
    // a small amount of numerical diffusion that quells oscillations
    //if (D.name() == "DD" || biMaterialInterfaceActive())
    if (true)
    {
        return
        (
            fvc::laplacian
            (
                gamma,
                D,
                "laplacian(D" + D.name() + ',' + D.name() + ')'
            )
          - fvc::div(gamma*mesh().Sf() & fvc::interpolate(gradD))
        );
    }
    else
    {
        // We will calculate this numerical diffusion based on the increment of
        // displacement, as it may become large of we base it on the total
        // displacement
        // Issue: The increment field "D - D.oldTime()" will be incorrect on
        // non-orthogonal meshes as the grad(D - D.oldTime()) field would not be
        // stored... we can/should fix this
        return
        (
            fvc::laplacian
            (
                gamma,
                D - D.oldTime(),
                "laplacian(D" + D.name() + ',' + D.name() + ')'
            )
          - fvc::div
            (
                gamma*mesh().Sf()
              & fvc::interpolate(gradD - gradD.oldTime())
            )
        );
    }
}


Foam::tmp<Foam::volVectorField> Foam::mechanicalModel::RhieChowCorrection
(
    const volVectorField& D,
    const volTensorField& gradD
) const
{
    return RhieChowCorrection(D, gradD, impKfcorr());
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


Foam::scalar Foam::mechanicalModel::newDeltaT()
{
    // Find the minimum time-step of all the mechanical laws
    PtrList<mechanicalLaw>& laws = *this;

    // Initial set deltaT to as large as possible and then check
    // if any mechanical law wants a smaller time-step
    scalar newDeltaT = mesh().time().endTime().value();

    forAll(laws, lawI)
    {
        newDeltaT = min(newDeltaT, laws[lawI].newDeltaT());
    }

    return newDeltaT;
}


void Foam::mechanicalModel::moveSubMeshes()
{
    PtrList<mechanicalLaw>& laws = *this;

    // Sub-meshes only exist when there is more than one material law
    if (laws.size() > 1)
    {
        forAll(subMeshes(), lawI)
        {
            Info<< "    Moving subMesh " << subMeshes()[lawI].subMesh().name()
                << endl;

            twoDPointCorrector twoDCorrector(subMeshes()[lawI].subMesh());
            pointField newPoints =
                subMeshes()[lawI].subMesh().points() + subMeshPointD()[lawI];
            twoDCorrector.correctPoints(newPoints);

            subMeshes()[lawI].subMesh().movePoints(newPoints);
            subMeshes()[lawI].subMesh().V00();
            subMeshes()[lawI].subMesh().moving(false);
            subMeshes()[lawI].subMesh().changing(false);
            subMeshes()[lawI].subMesh().setPhi().writeOpt() =
                IOobject::NO_WRITE;

            if
            (
                mesh().time().outputTime()
             && lookupOrDefault<Switch>("writeSubMeshes",  false)
            )
            {
                Info<< "    Writing subMesh "
                    << subMeshes()[lawI].subMesh().name() << endl;
                subMeshes()[lawI].subMesh().writeOpt() = IOobject::AUTO_WRITE;
                subMeshes()[lawI].subMesh().setInstance
                (
                    mesh().time().timeName()
                );
                subMeshes()[lawI].subMesh().write();
            }
        }
    }
}

// ************************************************************************* //
