/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "pisoFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pisoFluid, 0);
addToRunTimeSelectionTable(fluidModel, pisoFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pisoFluid::pisoFluid(const fvMesh& mesh)
:
    fluidModel(this->typeName, mesh),
    U_
    (
        IOobject
        (
            "U",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    pMesh_(mesh),
    pointD_
    (
        IOobject
        (
            "pointD",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    p_
    (
        IOobject
        (
            "p",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    gradp_(fvc::grad(p_)),
    phi_
    (
        IOobject
        (
            "phi",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_) & mesh.Sf()
    ),
    laminarTransport_(U_, phi_),
    turbulence_
    (
        incompressible::turbulenceModel::New
        (
            U_, phi_, laminarTransport_
        )
    ),
    rho_
    (
        IOdictionary
        (
            IOobject
            (
                "transportProperties",
                runTime().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).lookup("rho")
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const volVectorField& pisoFluid::U() const
{
    return U_;
}


const volScalarField& pisoFluid::p() const
{
    return p_;
}


//- Patch viscous force (N/m2)
 tmp<vectorField> pisoFluid::patchViscousForce(const label patchID) const
 {
     tmp<vectorField> tvF
     (
         new vectorField(mesh().boundary()[patchID].size(), vector::zero)
     );

     tvF() =
         rho_.value()
        *(
             mesh().boundary()[patchID].nf()
           & turbulence_->devReff()().boundaryField()[patchID]
         );

 //     vectorField n = mesh().boundary()[patchID].nf();
 //     tvF() -= n*(n&tvF());

     return tvF;
 }

// //- Patch pressure force (N/m2)
 tmp<scalarField> pisoFluid::patchPressureForce(const label patchID) const
 {
     tmp<scalarField> tpF
     (
         new scalarField(mesh().boundary()[patchID].size(), 0)
     );

     tpF() = rho_.value()*p().boundaryField()[patchID];

     return tpF;
 }

// //- Patch viscous force (N/m2)
 tmp<vectorField> pisoFluid::faceZoneViscousForce
 (
     const label zoneID,
     const label patchID
 ) const
 {
     vectorField pVF = patchViscousForce(patchID);

     tmp<vectorField> tvF
     (
         new vectorField(mesh().faceZones()[zoneID].size(), vector::zero)
     );
     vectorField& vF = tvF();

     const label patchStart =
         mesh().boundaryMesh()[patchID].start();

     forAll(pVF, i)
     {
         vF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
             pVF[i];
     }

//     // Parallel data exchange: collect pressure field on all processors
     reduce(vF, sumOp<vectorField>());


     return tvF;
 }

// //- Patch pressure force (N/m2)
 tmp<scalarField> pisoFluid::faceZonePressureForce
 (
     const label zoneID,
     const label patchID
 ) const
 {
     scalarField pPF = patchPressureForce(patchID);

     tmp<scalarField> tpF
     (
         new scalarField(mesh().faceZones()[zoneID].size(), 0)
     );
     scalarField& pF = tpF();

     const label patchStart =
         mesh().boundaryMesh()[patchID].start();

     forAll(pPF, i)
     {
         pF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
             pPF[i];
     }

//     // Parallel data exchange: collect pressure field on all processors
     reduce(pF, sumOp<scalarField>());

     return tpF;
 }

 tmp<scalarField> pisoFluid::faceZoneMuEff
 (
     const label zoneID,
     const label patchID
 ) const
 {
     scalarField pMuEff =
         rho_.value()*turbulence_->nuEff()().boundaryField()[patchID];

     tmp<scalarField> tMuEff
     (
         new scalarField(mesh().faceZones()[zoneID].size(), 0)
     );
     scalarField& muEff = tMuEff();

     const label patchStart =
         mesh().boundaryMesh()[patchID].start();

     forAll(pMuEff, i)
     {
         muEff[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
             pMuEff[i];
     }

//     // Parallel data exchange: collect pressure field on all processors
     reduce(muEff, sumOp<scalarField>());

     return tMuEff;
 }

tmp<vectorField>
 pisoFluid::currentFaceZonePoints(const label zoneID) const
 {
     vectorField pointDisplacement
     (
         mesh().faceZones()[zoneID]().localPoints().size(),
         vector::zero
     );

     const vectorField& pointDI = pointD_.internalField();

     label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

     if (globalZoneIndex != -1)
     {
         // global face zone
         const labelList& curPointMap =
             globalToLocalFaceZonePointMap()[globalZoneIndex];

         const labelList& zoneMeshPoints =
             mesh().faceZones()[zoneID]().meshPoints();

         vectorField zonePointsDisplGlobal
         (
             zoneMeshPoints.size(),
             vector::zero
         );

         //- Inter-proc points are shared by multiple procs
         //  pointNumProc is the number of procs which a point lies on
         scalarField pointNumProcs(zoneMeshPoints.size(), 0);

         forAll(zonePointsDisplGlobal, globalPointI)
         {
             label localPoint = curPointMap[globalPointI];

             if(zoneMeshPoints[localPoint] < mesh().nPoints())
             {
                 label procPoint = zoneMeshPoints[localPoint];

                 zonePointsDisplGlobal[globalPointI] =
                     pointDI[procPoint];

                 pointNumProcs[globalPointI] = 1;
             }
         }

         if (Pstream::parRun())
         {
             reduce(zonePointsDisplGlobal, sumOp<vectorField>());
             reduce(pointNumProcs, sumOp<scalarField>());

             //- now average the displacement between all procs
             zonePointsDisplGlobal /= pointNumProcs;
         }

         forAll(pointDisplacement, globalPointI)
         {
             label localPoint = curPointMap[globalPointI];

             pointDisplacement[localPoint] =
                 zonePointsDisplGlobal[globalPointI];
         }
     }
     else
     {
         pointDisplacement =
             vectorField
             (
                 pointDI,
                 mesh().faceZones()[zoneID]().meshPoints()
             );
     }

     tmp<vectorField> tCurrentPoints
     (
         new vectorField
         (
             mesh().faceZones()[zoneID]().localPoints()
           + pointDisplacement
         )
     );

     return tCurrentPoints;
 }

void pisoFluid::evolve()
{
    Info << "Evolving fluid model" << endl;

    const fvMesh& mesh = fluidModel::mesh();

    const int nCorr(readInt(fluidProperties().lookup("nCorrectors")));

    const int nNonOrthCorr =
        readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"));

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p_, fluidProperties(), pRefCell, pRefValue);

    if (mesh.moving())
    {
        // Make the fluxes relative
        phi_ -= fvc::meshPhi(U_);
    }

    // Calculate CourantNo
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;

        if (mesh.nInternalFaces())
        {
            surfaceScalarField SfUfbyDelta =
                mesh.surfaceInterpolation::deltaCoeffs()*mag(phi_);

            CoNum =
                max(SfUfbyDelta/mesh.magSf()).value()
               *runTime().deltaT().value();

            meanCoNum =
                (sum(SfUfbyDelta)/sum(mesh.magSf())).value()
               *runTime().deltaT().value();

            velMag = max(mag(phi_)/mesh.magSf()).value();
        }

        Info<< "Courant Number mean: " << meanCoNum
            << " max: " << CoNum
            << " velocity magnitude: " << velMag << endl;
    }

    // Construct momentum equation
    // Convection-diffusion matrix
    fvVectorMatrix HUEqn
    (
        fvm::div(phi_, U_)
      + turbulence_->divDevReff(U_)
    );

    // Time derivative matrix
    fvVectorMatrix ddtUEqn(fvm::ddt(U_));

    // Solve momentum equation
    solve(ddtUEqn + HUEqn == -gradp_);

    // --- PISO loop

    volScalarField aU = HUEqn.A();

    for (int corr = 0; corr < nCorr; corr++)
    {
        U_ = HUEqn.H()/aU;
        phi_ = (fvc::interpolate(U_) & mesh.Sf());

        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            // Construct pressure equation
            fvScalarMatrix pEqn
            (
                fvm::laplacian(1/aU, p_) == fvc::div(phi_)
            );

            // Solve pressure equation

            pEqn.setReference(pRefCell, pRefValue);

            if
            (
                corr == nCorr-1
             && nonOrth == nNonOrthCorr
            )
            {
                pEqn.solve(mesh.solver("pFinal"));
            }
            else
            {
                pEqn.solve();
            }

            if (nonOrth == nNonOrthCorr)
            {
                phi_ -= pEqn.flux();
            }
        }

        // Calculate Continuity error
        {
            volScalarField contErr = fvc::div(phi_);

            scalar sumLocalContErr =
                runTime().deltaT().value()
               *mag(contErr)().weightedAverage(mesh.V()).value();

            scalar globalContErr =
                runTime().deltaT().value()
               *contErr.weightedAverage(mesh.V()).value();

            Info<< "time step continuity errors : sum local = "
                << sumLocalContErr << ", global = " << globalContErr << endl;
        }

        gradp_ = fvc::grad(p_);

        U_ = 1.0/(aU + ddtUEqn.A())*
            (
                U_*aU - gradp_ + ddtUEqn.H()
            );
        U_.correctBoundaryConditions();
    }

    turbulence_->correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
