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

\*---------------------------------------------------------------------------*/

#include "mechanicalModel.H"
#include "fvc.H"
#include "fvcGradf.H"
#include "gaussGrad.H"
#include "twoDPointCorrector.H"
#include "fixedGradientFvPatchFields.H"
#include "wedgePolyPatch.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "ZoneIDs.H"
#else
    #include "ZoneID.H"
    #include "crackerFvMesh.H"
#endif

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::mechanicalModel::makeSolSubMeshes() const
{
    if (!solSubMeshes_.empty())
    {
        FatalErrorIn("void Foam::mechanicalModel::makeSolidSubMeshes() const")
            << "solid sub-meshes already exist" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        FatalErrorIn("void Foam::mechanicalModel::makeSubMeshes() const")
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    solSubMeshes_.set
    (
        new solidSubMeshes
        (
            mesh_,
            cellZoneNames_,
            incremental_,
            lookupOrDefault<Switch>("writeSubMeshes",  false)
        )
    );
}


const Foam::solidSubMeshes& Foam::mechanicalModel::solSubMeshes() const
{
    if (solSubMeshes_.empty())
    {
        makeSolSubMeshes();
    }

    return solSubMeshes_();
}


Foam::solidSubMeshes& Foam::mechanicalModel::solSubMeshes()
{
    if (solSubMeshes_.empty())
    {
        makeSolSubMeshes();
    }

    return solSubMeshes_();
}


void Foam::mechanicalModel::calcImpKfcorr() const
{
    if (impKfcorrPtr_.valid())
    {
        FatalErrorIn
        (
            "const Foam::volScalarField& "
            "Foam::mechanicalModel::calcImpKfcorr() const"
        )   << "pointer already set" << abort(FatalError);
    }

    impKfcorrPtr_.set
    (
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
        )
    );

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() > 1)
    {
#ifndef OPENFOAM_NOT_EXTEND
        // To disable Rhie-Chow correction on bi-material interface, we will set
        // impKfcorr to zero on bi-material interface faces

        surfaceScalarField& impKfcorr = impKfcorrPtr_();
        scalarField& impKfcorrI = impKfcorrPtr_().internalField();

        forAll(laws, lawI)
        {
            const fvMesh& subMesh = solSubMeshes().subMeshes()[lawI].subMesh();
            const labelList& patchMap =
                solSubMeshes().subMeshes()[lawI].patchMap();
            const labelList& faceMap =
                solSubMeshes().subMeshes()[lawI].faceMap();

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

        impKfcorrPtr_().correctBoundaryConditions();
#else
        FatalErrorIn(type())
            << "Not implemented for this version of OpenFOAM"
            << abort(FatalError);
#endif
    }
}


const Foam::surfaceScalarField& Foam::mechanicalModel::impKfcorr() const
{
    if (impKfcorrPtr_.empty())
    {
        calcImpKfcorr();
    }

    return impKfcorrPtr_();
}


void Foam::mechanicalModel::clearOut()
{
    // Clear the list of mechanical laws
    // Note: we should do this before clearing the subMeshes, as the mechanical
    // laws can store geometricFields that must be deleted before deleting
    // mesh
    PtrList<mechanicalLaw>::clear();

    // Make sure to clear the subMeshes after (not before) clearing the subMesh
    // fields
    solSubMeshes_.clear();
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
    solSubMeshes_(),
    impKfcorrPtr_()
{
    Info<< "Creating the mechanicalModel" << endl;

    // Read the mechanical laws
    const PtrList<entry> lawEntries(lookup("mechanical"));

    PtrList<mechanicalLaw>& laws = *this;
    laws.setSize(lawEntries.size());

    // Create the list of cellZones names: they are used during the construction
    // of the subMeshes
    cellZoneNames_.setSize(laws.size());
    forAll(laws, lawI)
    {
        cellZoneNames_[lawI] = lawEntries[lawI].keyword();
    }

    // Create mechancial laws
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
                        solSubMeshes().subMeshes()[lawI].subMesh(),
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
                        solSubMeshes().subMeshes()[lawI].subMesh(),
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
        }
    }

#ifndef OPENFOAM_NOT_EXTEND
    // Check: currently crackerFvMesh only works with a single material
    // The challenge here is to update the subMesh and subMesh fields when a
    // topo-change (crack) occurs in the babse mesh
    if (isA<crackerFvMesh>(mesh) && laws.size() > 1)
    {
        FatalErrorIn(type() + "::" + type() + "(...)")
            << "Currently the crackerFvMesh can only be used with a single "
            << "material in mechanicalProperties"
            << abort(FatalError);
    }
#endif
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


#ifdef OPENFOAM_NOT_EXTEND
const Foam::volPointInterpolation& Foam::mechanicalModel::volToPoint() const
{
    return volPointInterpolation::New(mesh_);
}
#else
const Foam::newLeastSquaresVolPointInterpolation& Foam::mechanicalModel::volToPoint() const
{
    return newLeastSquaresVolPointInterpolation::New(mesh_);
}
#endif

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

#ifdef OPENFOAM_NOT_EXTEND
        volScalarField& result = tresult.ref();
#else
        volScalarField& result = tresult();
#endif

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
        solSubMeshes().mapSubMeshVolFields<scalar>(rhos, result);

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

#ifdef OPENFOAM_NOT_EXTEND
        volScalarField& result = tresult.ref();
#else
        volScalarField& result = tresult();
#endif

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
        solSubMeshes().mapSubMeshVolFields<scalar>(impKs, result);

        // Clear subMesh fields
        impKs.clear();

        return tresult;
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::mechanicalModel::impKf() const
{
    // Linear interpolation actually seems to give the best convergence
    const volScalarField impK(this->impK());
    const word interpName = "interpolate(" + impK.name() + ')';
    return fvc::interpolate(impK, interpName);
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::bulkModulus() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        return laws[0].bulkModulus();
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
                    "bulkModulusLaw",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimPressure, 0),
                calculatedFvPatchScalarField::typeName
            )
        );

#ifdef OPENFOAM_NOT_EXTEND
        volScalarField& result = tresult.ref();
#else
        volScalarField& result = tresult();
#endif

        // Accumulated subMesh fields and then map to the base mesh
        PtrList<volScalarField> bulkModulii(laws.size());

        forAll(laws, lawI)
        {
            bulkModulii.set
            (
                lawI,
                new volScalarField(laws[lawI].bulkModulus())
            );
        }

        // Map subMesh fields to the base mesh
        solSubMeshes().mapSubMeshVolFields<scalar>(bulkModulii, result);

        // Clear subMesh fields
        bulkModulii.clear();

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
            laws[lawI].correct(solSubMeshes().subMeshSigma()[lawI]);
        }

        // Map subMesh fields to the base field
        solSubMeshes().mapSubMeshVolFields<symmTensor>
        (
            solSubMeshes().subMeshSigma(), sigma
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
            laws[lawI].correct(solSubMeshes().subMeshSigmaf()[lawI]);
        }

        // Map subMesh fields to the base field
        solSubMeshes().mapSubMeshSurfaceFields<symmTensor>
        (
            solSubMeshes().subMeshSigmaf(), sigma
        );
    }
}


void Foam::mechanicalModel::mapGradToSubMeshes(const volTensorField& gradD)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        return;
    }

    forAll(laws, lawI)
    {
        solSubMeshes().subMeshGradD()[lawI] =
            solSubMeshes().subMeshes()[lawI].interpolate(gradD);
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
        solSubMeshes().interpolateDtoSubMeshD(D, true);

        // Accumulate data for all fields
        forAll(laws, lawI)
        {
            // Calculate gradient on subMesh
            // This will use the values at the interface
            volTensorField& subMeshGradD = solSubMeshes().subMeshGradD()[lawI];
            subMeshGradD = fvc::grad(solSubMeshes().subMeshD()[lawI]);
        }

        // Map subMesh gradD to the base gradD
        solSubMeshes().correctBoundarySnGrad
        (
            solSubMeshes().subMeshD(), solSubMeshes().subMeshGradD()
        );

        solSubMeshes().mapSubMeshVolFields<tensor>
        (
            solSubMeshes().subMeshGradD(), gradD
        );

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
            volTensorField& subMeshGradD = solSubMeshes().subMeshGradD()[lawI];
            subMeshGradD = fvc::grad
            (
                solSubMeshes().subMeshD()[lawI],
                solSubMeshes().subMeshPointD()[lawI]
            );
        }

        // Correct snGrad on boundaries because the default calculated
        // boundaries disable non-orthogonal correction
        solSubMeshes().correctBoundarySnGrad
        (
            solSubMeshes().subMeshD(), solSubMeshes().subMeshGradD()
        );

        // Map subMesh gradD fields to the base gradD field
        solSubMeshes().mapSubMeshVolFields<tensor>
        (
            solSubMeshes().subMeshGradD(), gradD
        );

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
            surfaceTensorField& subMeshGradDf =
                solSubMeshes().subMeshGradDf()[lawI];
            subMeshGradDf = fvc::fGrad
            (
                solSubMeshes().subMeshD()[lawI],
                solSubMeshes().subMeshPointD()[lawI]
            );
        }

        // Correct snGrad on boundaries because the default calculated
        // boundaries disable non-orthogonal correction
        solSubMeshes().correctBoundarySnGradf
        (
            solSubMeshes().subMeshD(),
            solSubMeshes().subMeshGradDf(),
            solSubMeshes().subMeshGradD()
        );

        // Map subMesh gradDf fields to the base gradDf field
        solSubMeshes().mapSubMeshSurfaceFields<tensor>
        (
            solSubMeshes().subMeshGradDf(),
            gradDf
        );

        // Replace normal component
        // If we don't do this then we don't get convergence in many cases
        const surfaceVectorField n(mesh().Sf()/mesh().magSf());
        gradDf += n*fvc::snGrad(D) - (n*(n & gradDf));
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
            volTensorField& subMeshGradD = solSubMeshes().subMeshGradD()[lawI];
            subMeshGradD = fvc::grad
            (
                solSubMeshes().subMeshD()[lawI],
                solSubMeshes().subMeshPointD()[lawI]
            );

            surfaceTensorField& subMeshGradDf =
                solSubMeshes().subMeshGradDf()[lawI];
            subMeshGradDf = fvc::fGrad
            (
                solSubMeshes().subMeshD()[lawI],
                solSubMeshes().subMeshPointD()[lawI]
            );
        }

        // Correct snGrad on boundaries because the default calculated
        // boundaries disable non-orthogonal correction
        solSubMeshes().correctBoundarySnGrad
        (
            solSubMeshes().subMeshD(), solSubMeshes().subMeshGradD()
        );
        solSubMeshes().correctBoundarySnGradf
        (
            solSubMeshes().subMeshD(),
            solSubMeshes().subMeshGradDf(),
            solSubMeshes().subMeshGradD()
        );

        // Map subMesh fields to the base field

        solSubMeshes().mapSubMeshVolFields<tensor>
        (
            solSubMeshes().subMeshGradD(), gradD
        );
        solSubMeshes().mapSubMeshSurfaceFields<tensor>
        (
            solSubMeshes().subMeshGradDf(), gradDf
        );

        // Correct boundary snGrad of gradD
        fv::gaussGrad<vector>
        (
            mesh()
        ).correctBoundaryConditions(D, gradD);

        // Correct snGrad component of gradDf
        const surfaceVectorField n(mesh().Sf()/mesh().magSf());
        gradDf += n*fvc::snGrad(D) - (n*(n & gradDf));
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
        solSubMeshes().interpolateDtoSubMeshD(D, useVolFieldSigma);

        // Accumulate data for all fields
        forAll(laws, lawI)
        {
            // Interpolate the subMeshD to the subMeshPointD
            solSubMeshes().subMeshVolToPoint()[lawI].interpolate
            (
                solSubMeshes().subMeshD()[lawI],
                solSubMeshes().subMeshPointD()[lawI]
            );
        }

        // Map subMesh pointD fields back to the base pointD field
        solSubMeshes().mapSubMeshPointFields<vector>
        (
            solSubMeshes().subMeshPointD(), pointD
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
    if (solSubMeshes_.valid())
    {
        solSubMeshes().moveSubMeshes();
    }
}


void Foam::mechanicalModel::setRestart()
{
    PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        laws[lawI].setRestart();
    }
}

void Foam::mechanicalModel::writeDict()
{
    // This has to be done, because law dicts can only be read as IStream and so
    // they are not references to the entries in the 'mechanicalProperties' but
    // seperate entry object and dictionary objects. Because of this they have
    // to be added back to 'mechanicalProperties' before writing to disk.
    PtrList<mechanicalLaw>& laws = *this;

    // Adding all law dictionaries to a single pointer list
    PtrList<primitiveEntry> lawDicts;
    lawDicts.setSize(laws.size());
    forAll(laws,lawI)
    {
        lawDicts.set
        (
            lawI,
            new primitiveEntry
            (
                laws[lawI].name(),
                laws[lawI].dict()
            )
        );
    }

    // Creating an entry from pointer list and replacing 'mechanical' entry in
    // IOdictionary ('materialProperties') with it
    primitiveEntry lawsEntry("mechanical",lawDicts);
    autoPtr<IOdictionary> outputMechLawProps
    (
        new IOdictionary
        (
            IOobject
            (
                "mechanicalProperties.withDefaultValues",
                mesh().time().constant(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        )
    );

#ifdef OPENFOAM_ESI
    outputMechLawProps.ref() = *this;
    outputMechLawProps.ref().IOdictionary::set(lawsEntry);
#else
    outputMechLawProps() = *this;
    outputMechLawProps().IOdictionary::set(lawsEntry);
#endif

    // Writing to disk
    outputMechLawProps().regIOobject::write();
    outputMechLawProps.clear();
}

// ************************************************************************* //
