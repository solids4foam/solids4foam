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

#include "newtonIcoFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "CorrectPhi.H"
#include "fvc.H"
#include "fvm.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "findRefCell.H"
#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "elasticWallVelocityFvPatchVectorField.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "movingWallPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(newtonIcoFluid, 0);
addToRunTimeSelectionTable(fluidModel, newtonIcoFluid, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //


void newtonIcoFluid::makeRAUf() const
{
    if (rAUfPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    rAUfPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "rAUf",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("0", dimTime, 0.0)
        )
    );
}


const surfaceScalarField& newtonIcoFluid::rAUf() const
{
    if (rAUfPtr_.empty())
    {
        makeRAUf();
    }

    return rAUfPtr_.ref();
}


surfaceScalarField& newtonIcoFluid::rAUf()
{
    if (rAUfPtr_.empty())
    {
        makeRAUf();
    }

    return rAUfPtr_.ref();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

newtonIcoFluid::newtonIcoFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    foamPetscSnesHelper
    (
        fileName
        (
            fluidProperties().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        mesh().nCells(),
        fluidProperties().lookupOrDefault<Switch>("stopOnPetscError", true),
        true
    ),
    Uf_(),
    rAUfPtr_(),
    pRefCell_(-1),
    pRefValue_(0.0),
    laminarTransport_(U(), phi()),
    turbulence_
    (
        incompressible::turbulenceModel::New
        (
            U(), phi(), laminarTransport_
        )
    ),
    rho_(laminarTransport_.lookup("rho")),
    blockSize_(fluidModel::twoD() ? 3 : 4)
{
    setRefCell(p(), fluidProperties(), pRefCell_, pRefValue_);
    //mesh().setFluxRequired(p().name());
    turbulence_->validate();

    U().oldTime().oldTime();

    // if (mesh().dynamic())
    // {
    //     Info<< "Constructing face velocity Uf\n" << endl;

        // Uf_.reset
        // (
        //     new surfaceVectorField
        //     (
        //         IOobject
        //         (
        //             "Uf",
        //             runTime.timeName(),
        //             mesh(),
        //             IOobject::READ_IF_PRESENT,
        //             IOobject::AUTO_WRITE
        //         ),
        //         fvc::interpolate(U())
        //     )
        // );

        // Uf_().oldTime();

        // if
        // (
        //     word(mesh().ddtScheme("ddt(" + U().name() +')'))
        //  == fv::backwardDdtScheme<vector>::typeName
        // )
        // {
        //     Uf_().oldTime().oldTime();
        // }
    // }

    const fvMesh& mesh = this->mesh();
    const surfaceScalarField& phi = this->phi();
    #include "CourantNo.H"
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> newtonIcoFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF.ref() = rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & (-turbulence_->devReff()().boundaryField()[patchID])
        );

    return tvF;
}


tmp<scalarField> newtonIcoFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF.ref() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


tmp<vectorField> newtonIcoFluid::patchViscousForce
(
    const label patchID, const solidModel& motion
) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    const tensorField Fm(I + motion.gradD().T());
    const scalarField Jm(det(Fm));
    const tensorField invFm(inv(Fm));
    const scalarField nuEff(turbulence_->nuEff()().boundaryField()[patchID]);
    const vectorField& Sf = mesh().boundary()[patchID].Sf();
    const vectorField deformedSf(Jm*invFm.T() & Sf);
    const vectorField deformedNf(deformedSf/mag(deformedSf));
    const tensorField gradU
    (
        fvc::grad(U())().boundaryField()[patchID]
    ); // do I need to calculate the whole field?
    tvF.ref() = rho_.value()*deformedNf & (nuEff*invFm.T() & gradU);

    // Deformed mesh
    // tvF.ref() = rho_.value()
    //    *(
    //         mesh().boundary()[patchID].nf()
    //       & (-turbulence_->devReff()().boundaryField()[patchID])
    //     );

    return tvF;
}


bool newtonIcoFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    // Take references
    // const Time& runTime = fluidModel::runTime();
    dynamicFvMesh& mesh = this->mesh();
    volVectorField& U = this->U();
    volScalarField& p = this->p();
    surfaceScalarField& phi = this->phi();
    // autoPtr<surfaceVectorField>& Uf = Uf_;
    //scalar& cumulativeContErr = cumulativeContErr_;
    //const bool correctPhi = correctPhi_;
    //const bool checkMeshCourantNo = checkMeshCourantNo_;
    //const bool moveMeshOuterCorrectors = moveMeshOuterCorrectors_;

    // Update U boundary conditions
    U.correctBoundaryConditions();

    // Solution predictor
    // if (predictor_ && newTimeStep())
    // {
    //     predict();

    //     // Map the U field to the SNES solution vector
    //     foamPetscSnesHelper::InsertFieldComponents<vector>
    //     (
    //         U.primitiveFieldRef(),
    //         foamPetscSnesHelper::solution(),
    //         fluidModel::twoD() ? 2 : 3, // Block size of x
    //         0                           // Location of first component
    //     );
    //
    //     // Update phi field
    //
    //     // Also predict p field
    //     // Todo
    // }

    // Update the mesh
    mesh.controlledUpdate();

    // Update the flux
    phi = fvc::interpolate(U) & mesh.Sf();

    // If the mesh moved, update the flux and make it relative to the mesh
    // motion
    if (mesh.changing())
    {
        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);
    }

    // Solve the nonlinear system and check the convergence
    Info<< "Solving the fluid for U and p" << endl;
    foamPetscSnesHelper::solve();

    // Access the raw solution data
    const PetscScalar *xx;
    VecGetArrayRead(foamPetscSnesHelper::solution(), &xx);

    // Retrieve the solution
    // Map the PETSc solution to the U field
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        xx,
        U.primitiveFieldRef(),
        0, // Location of U
        blockSize_,
        fluidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    U.correctBoundaryConditions();

    // Map the PETSc solution to the p field
    // p is located in the final component
    foamPetscSnesHelper::ExtractFieldComponents<scalar>
    (
        xx,
        p.primitiveFieldRef(),
        blockSize_ - 1, // Location of p component
        blockSize_
    );

    p.correctBoundaryConditions();

    // Correct Uf if the mesh is moving
    //fvc::correctUf(Uf, U, phi);

    // Update the flux
    //phi = mesh.Sf() & Uf();
    phi = mesh.Sf() & fvc::interpolate(U);

    if (mesh.changing())
    {
        // Enable: needed for inlet/outlet?
        // #include "correctPhi.esi.H"

        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);
    }

    // Correct transport and turbulence models
    laminarTransport_.correct();
    turbulence_->correct();

    return 0;
}


void newtonIcoFluid::clearRAUf()
{
    rAUfPtr_.clear();
}


label newtonIcoFluid::initialiseJacobian(Mat& jac)
{
    // Initialise based on compact stencil fvMesh
    return Foam::initialiseJacobian(jac, mesh(), blockSize_);
}


label newtonIcoFluid::initialiseSolution(Vec& x)
{
    // Initialise based on mesh.nCells()
    return Foam::initialiseSolution(x, mesh(), blockSize_);
}


label newtonIcoFluid::formResidual
(
    PetscScalar *f,
    const PetscScalar *x
)
{
    if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    // Take references
    //const Time& runTime = fluidModel::runTime();
    dynamicFvMesh& mesh = this->mesh();
    volVectorField& U = const_cast<volVectorField&>(this->U());
    volScalarField& p = const_cast<volScalarField&>(this->p());
    surfaceScalarField& phi = const_cast<surfaceScalarField&>(this->phi());
    //autoPtr<surfaceVectorField>& Uf = Uf_;
    //scalar& cumulativeContErr = cumulativeContErr_;
    //const bool correctPhi = correctPhi_;
    // const bool checkMeshCourantNo = checkMeshCourantNo_;
    //const bool moveMeshOuterCorrectors = moveMeshOuterCorrectors_;

    // Copy x into the U field
    vectorField& UI = U;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        UI,
        0, // Location of first UI component
        blockSize_, // Block size of x
        fluidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    // Enforce the boundary conditions
    U.correctBoundaryConditions();

    // Update the flux
    phi = fvc::interpolate(U) & mesh.Sf();

    if (mesh.changing())
    {
        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);
    }

    // TESTING
    // Set the flux to zero on walls, including FSI interfaces
    // makeRelative should do this but mat not work as expected
    forAll(U.boundaryField(), patchI)
    {
        if (mesh.boundaryMesh()[patchI].type() == "wall")
        {
            phi.boundaryFieldRef()[patchI] = 0.0;
        }
    }

    // Copy x into the p field
    scalarField& pI = p;
    foamPetscSnesHelper::ExtractFieldComponents<scalar>
    (
        x, pI, blockSize_ - 1, blockSize_
    );

    // Enforce the boundary conditions
    p.correctBoundaryConditions();

    // Update the pressure BCs to ensure flux consistency
    // constrainPressure(p, U, phiHbyA, rAtU(), MRF);
    // CHECK
    //constrainPressure(p, U, phiHbyA, rAtU());

    // Correct Uf if the mesh is moving
    //fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    //fvc::makeRelative(phi, U);

    // Correct the transport and turbulence models
    laminarTransport_.correct();
    turbulence_->correct();

    // The residual vector is defined as
    // F = div(sigma) - ddt(U) - div(phi*U)
    //   = div(dev(sigma)) - grad(p) - ddt(U) - div(phi*U)
    //   = div(2*nuEff*symm(gradU)) - grad(p) - ddt(U) - div(phi*U)
    //   = laplacian(nuEff,U) + div(nuEff*gradU.T())
    //     - grad(p) - ddt(U) - div(phi*U)
    //
    // Check: do we want to include div(gradU.T).. it makes the stencil
    // larger and should be zero anyway, although it may increase accuracy
    // To be checked ...
    vectorField residual
    (
        fvc::laplacian(turbulence_->nuEff(), U)
        //+ fvc::div((turbulence_->nuEff())*dev2(T(fvc::grad(U))))
      - fvc::grad(p)
      - fvc::ddt(U)
      - fvc::div(phi, U)
    );

    // Make residual extensive as fvc operators are intensive (per unit volume)
    residual *= mesh.V();

    // Copy the residual into the f field
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        residual,
        f,
        0, // Location of first component
        blockSize_, // Block size of x
        fluidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    // Calculate pressure equation residual
    // Fp = stabilisation - div(U)
    //    = stabilisation - tr(grad(U))
    // where stabilisation = laplacian(pD, p) - div(pD*grad(p))
    const dimensionedScalar omega("omega", fluidProperties());
    const scalar localReRef(readScalar(fluidProperties().lookup("localReRef")));
    const scalar omegaExponent
    (
        readScalar(fluidProperties().lookup("omegaExponent"))
    );
    const surfaceScalarField localRe
    (
        fvc::interpolate(mag(U)/turbulence_->nuEff())/mesh.deltaCoeffs()
    );
    scalarField pressureResidual
    (
        fvc::laplacian
        (
            omega*Foam::pow(1.0 + localRe/localReRef, omegaExponent)
           /sqr(mesh.deltaCoeffs()),
            p,
            "laplacian(rAU,p)"
        )
      //   fvc::laplacian(rAUf(), p, "laplacian(rAU,p)")
      // - fvc::div
      //   (
      //       (rAUf()*mesh.Sf()) & fvc::interpolate(fvc::grad(p))
      //   )
      // - tr(fvc::grad(U)) // probably more accurate on a bad grid?
      - fvc::div(phi)
      // - fvc::div(U)
    );

    //Info<< "localRe/localReRef = " << max(mag(localRe/localReRef)) << endl;

    // Make residual extensive
    pressureResidual *= mesh.V();

    // If required, set pressure reference value
    if (pRefCell_ != -1)
    {
        // Info<< "Setting the pressure residual row for cell " << pRefCell_
        //     << " to be zero" << endl;

        // Set the residual to zero for the pressure equation in the pRefCell
        // cell
        pressureResidual[pRefCell_] = 0.0;

        // Set p in the pRefCell
        p[pRefCell_] = pRefValue_;
    }

    // Copy the pressureResidual into the f field as the final equation
    foamPetscSnesHelper::InsertFieldComponents<scalar>
    (
        pressureResidual, f, blockSize_ - 1, blockSize_
    );

    return 0;
}


label newtonIcoFluid::formResidual
(
    PetscScalar *f,
    const PetscScalar *x,
    const solidModel& motion
)
{
    if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    // Take references
    //const Time& runTime = fluidModel::runTime();
    dynamicFvMesh& mesh = this->mesh();
    volVectorField& U = const_cast<volVectorField&>(this->U());
    volScalarField& p = const_cast<volScalarField&>(this->p());
    surfaceScalarField& phi = const_cast<surfaceScalarField&>(this->phi());

    // Lookup the motion gradD and calculate the deformation gradient and its
    // determinant
    const volTensorField Fm(I + motion.gradD().T());
    const volScalarField Jm(det(Fm));
    const volTensorField invFm(inv(Fm));
    const surfaceTensorField Fmf(fvc::interpolate(Fm));
    const surfaceScalarField Jmf(det(Fmf));
    const surfaceTensorField invFmf(inv(Fmf));

    //autoPtr<surfaceVectorField>& Uf = Uf_;
    //scalar& cumulativeContErr = cumulativeContErr_;
    //const bool correctPhi = correctPhi_;
    // const bool checkMeshCourantNo = checkMeshCourantNo_;
    //const bool moveMeshOuterCorrectors = moveMeshOuterCorrectors_;

    // Copy x into the U field
    vectorField& UI = U;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        UI,
        0, // Location of first UI component
        blockSize_, // Block size of x
        fluidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    // Enforce the boundary conditions
    U.correctBoundaryConditions();

    // Calculate the gradient of velocity
    const volTensorField gradU(fvc::grad(U));
    const surfaceTensorField gradUf(fvc::interpolate(gradU));

    // Calculate the deformed area vectors
    const surfaceVectorField deformedSf(Jmf*invFmf.T() & mesh.Sf());

    // Calculate the relative flux
    //phi = fvc::interpolate(U) & mesh.Sf();
    phi = deformedSf & (fvc::interpolate(U - motion.U()));

    // See comment below
    // if (mesh.changing())
    // {
    //     // Make the flux relative to the mesh motion
    //     fvc::makeRelative(phi, U);
    // }

    const volScalarField nuEff(turbulence_->nuEff());
    const surfaceScalarField nuEfff(fvc::interpolate(nuEff));

    // WIP
    // Set the flux to zero on walls, including FSI interfaces
    // makeRelative should do this but mat not work as expected
    // But this solution may not be right for non-FSI cases, where the mesh is
    // moved at the start of the time step
    // Wait, this is not needed since U == motion.U() at FSI interfaces
    // forAll(U.boundaryField(), patchI)
    // {
    //     if (mesh.boundaryMesh()[patchI].type() == "wall")
    //     {
    //         Info<< "Setting the flux to 0 on patch "
    //             << mesh.boundaryMesh()[patchI].name() << endl;
    //         phi.boundaryFieldRef()[patchI] = 0.0;
    //     }
    // }

    // Copy x into the p field
    scalarField& pI = p;
    foamPetscSnesHelper::ExtractFieldComponents<scalar>
    (
        x, pI, blockSize_ - 1, blockSize_
    );

    // Enforce the boundary conditions
    p.correctBoundaryConditions();

    // Update the pressure BCs to ensure flux consistency
    // constrainPressure(p, U, phiHbyA, rAtU(), MRF);
    // CHECK
    //constrainPressure(p, U, phiHbyA, rAtU());

    // Correct Uf if the mesh is moving
    //fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    //fvc::makeRelative(phi, U);

    // Correct the transport and turbulence models
    // NOTE: these assume the fluid mesh is in the deformed configuration and may
    // not be correct in the reference configuration
    // DISABLED for now
    // laminarTransport_.correct();
    // turbulence_->correct();

    // Deformed configuration
    // The residual vector is defined as
    // F = div(sigma) - ddt(U) - div(phi*U)
    //   = div(dev(sigma)) - grad(p) - ddt(U) - div(phi*U)
    //   = div(2*nuEff*symm(gradU)) - grad(p) - ddt(U) - div(phi*U)
    //   = laplacian(nuEff,U) + div(nuEff*gradU.T())
    //     - grad(p) - ddt(U) - div(phi*U)
    //
    // Check: do we want to include div(gradU.T).. it makes the stencil
    // larger and should be zero anyway, although it may increase accuracy
    // To be checked ...
    // ... ignored below
    //
    // Reference configuration
    // The residual vector is defined as
    // F = div((Jmf*invFmf.T() & Sf) & (nuEff*invFmf.T() & gradUf))
    //     - (Jm*invFm.T() & grad(p))
    //     - Jm*ddt(U)
    //     - div(phi*U)
    //     + Du
    //
    // where phi is calculated in terms of the reference configuration
    // quantities
    // Check: div(phi,U) will use the reference mesh weights but we need the
    // deformed configuration
    // Du is the stabilisation term, where
    // Du = alphaU*laplacian(nuEff,U) - alphaU*div(nuEff*gradU)
    //

    // Lookup the stabilisation scale factor
    const scalar alphaU(readScalar(fluidProperties().lookup("alphaU")));

    // Calculate the residual over the reference configuration
    vectorField residual
    (
        fvc::div(deformedSf & (nuEfff*invFmf.T() & gradUf))
      - (Jm*invFm.T() & fvc::grad(p))
      - Jm*fvc::ddt(U)
      - fvc::div(phi, U)
      + alphaU*fvc::laplacian(nuEfff, U)
      - alphaU*fvc::div(mesh.Sf() & nuEfff*gradUf)
    );

    // Make residual extensive as fvc operators are intensive (per unit volume)
    residual *= mesh.V();

    // Copy the residual into the f field
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        residual,
        f,
        0, // Location of first component
        blockSize_, // Block size of x
        fluidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    // Deformed configuration
    // Calculate pressure equation residual
    // Fp = Dp - div(U) // or Dp - tr(grad(U))
    // where the stabilisation Dp = laplacian(pD, p) - div(pD*grad(p))
    //
    // Reference configuration
    // Fp = Dp - div((Jmf*invFmf.T() & Sf) & Uf)
    //

    // Update rAUf
    // {
    //     const scalar pressureSmoothingCoeff
    //     (
    //         readScalar(fluidProperties().lookup("pressureSmoothingCoeff"))
    //     );
    //     rAUf() = pressureSmoothingCoeff*mesh.magSf()/nuEfff;
    // }

    // const surfaceVectorField n(mesh.Sf()/mesh.magSf());
    // const surfaceVectorField gradpf(fvc::interpolate(fvc::grad(p)));
    // const surfaceScalarField snGradp(fvc::snGrad(p));

    const dimensionedScalar omega("omega", fluidProperties());
    const scalar localReRef(readScalar(fluidProperties().lookup("localReRef")));
    const scalar omegaExponent
    (
        readScalar(fluidProperties().lookup("omegaExponent"))
    );
    const surfaceScalarField localRe
    (
        fvc::interpolate(mag(U)/turbulence_->nuEff())/mesh.deltaCoeffs()
    );
    scalarField pressureResidual
    (
        fvc::laplacian
        (
            omega*Foam::pow(1.0 + localRe/localReRef, omegaExponent)
           /sqr(mesh.deltaCoeffs()),
            p,
            "laplacian(rAU,p)"
        )
        // Reference configuration
      //   fvc::laplacian(rAUf(), p, "laplacian(rAU,p)")
      // - fvc::div
      //   (
      //       (rAUf()*mesh.Sf()) & fvc::interpolate(fvc::grad(p))
      //   )
        // // Deformed configuration
        // fvc::div
        // (
        //     deformedSf & (rAUf()*invFmf.T() & (n*snGradp - gradpf))
        // )
        //deformedSf
        // - tr(fvc::grad(U)) // probably more accurate on a bad grid?
      - fvc::div(phi)
      // - fvc::div(U)
    );

    // Make residual extensive
    pressureResidual *= mesh.V();

    // If required, set pressure reference value
    if (pRefCell_ != -1)
    {
        // Set the residual to zero for the pressure equation in the pRefCell
        // cell
        pressureResidual[pRefCell_] = pRefValue_;

        // Set p in the pRefCell
        p[pRefCell_] = pRefValue_;
    }

    // Copy the pressureResidual into the f field as the final equation
    foamPetscSnesHelper::InsertFieldComponents<scalar>
    (
        pressureResidual, f, blockSize_ - 1, blockSize_
    );

    return 0;
}


label newtonIcoFluid::formJacobian
(
    Mat jac,
    const PetscScalar *x
)
{
    if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    const fvMesh& mesh = this->mesh();

    // Copy x into the U field
    volVectorField& U = const_cast<volVectorField&>(this->U());
    vectorField& UI = U;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        UI,
        0, // Location of first component
        blockSize_, // Block size of x
        fluidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    // Enforce the boundary conditions
    U.correctBoundaryConditions();

    // Update the flux
    phi() = fvc::interpolate(U) & mesh.Sf();

    if (mesh.changing())
    {
        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi(), U);
    }

    // Copy x into the p field
    volScalarField& p = const_cast<volScalarField&>(this->p());
    scalarField& pI = p;
    foamPetscSnesHelper::ExtractFieldComponents<scalar>
    (
        x, pI, blockSize_ - 1, blockSize_
    );

    // Enforce the boundary conditions
    p.correctBoundaryConditions();

    // Correct the transport and turbulence models
    laminarTransport_.correct();
    turbulence_->correct();

    // Calculate the segregated approximatoion of momentum equation Jacobian
    // Note: the nonlinear convection term is added separately below
    fvVectorMatrix UEqn
    (
        fvm::laplacian(turbulence_->nuEff(), U)
      - fvm::ddt(U)
    );

    UEqn.relax();

    if (debug)
    {
        Info<< "Inserting U equation in Afluid" << endl;
    }

    // Convert fvMatrix matrix to PETSc matrix
    foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
    (
        UEqn, jac, 0, 0, fluidModel::twoD() ? 2 : 3
    );

    // Insert linearisation of convection term
    // The linearisation assumes an upwind discretisation
    foamPetscSnesHelper::InsertFvmDivPhiUIntoPETScMatrix
    (
        U,
        phi(),
        jac,
        0,                         // row offset
        0,                         // column offset
        fluidModel::twoD() ? 2 : 3 // number of scalar equations to insert
    );

    // Update rAUf
    // const scalar pressureSmoothingCoeff
    // (
    //     readScalar(fluidProperties().lookup("pressureSmoothingCoeff"))
    // );
    // {
    //     const volScalarField nuEff(turbulence_->nuEff());
    //     const surfaceScalarField nuEfff(fvc::interpolate(nuEff));
    //     rAUf() = pressureSmoothingCoeff*mesh.magSf()/nuEfff;
    // }
    const dimensionedScalar omega("omega", fluidProperties());
    const scalar localReRef(readScalar(fluidProperties().lookup("localReRef")));
    const scalar omegaExponent
    (
        readScalar(fluidProperties().lookup("omegaExponent"))
    );
    const surfaceScalarField localRe
    (
        fvc::interpolate(mag(U)/turbulence_->nuEff())/mesh.deltaCoeffs()
    );

    // Calculate pressure equation matrix
    fvScalarMatrix pEqn
    (
        // fvm::laplacian(rAUf(), p, "jacobian-laplacian(rAU,p)")
        fvm::laplacian
        (
            omega*Foam::pow(1.0 + localRe/localReRef, omegaExponent)
           /sqr(mesh.deltaCoeffs()),
            p,
            "jacobian-laplacian(rAU,p)"
        )
    );

    if (debug)
    {
        Info<< "Inserting p equation in Afluid" << endl;
    }

    // If required, set pressure reference value
    if (pRefCell_ != -1)
    {
        // Info<< "Setting the pressure equation row for cell " << pRefCell_
        //     << " to be diagonal" << endl;

        // Set the off-diagonal to zero for cell pRefCell
        pEqn.setValues(labelList(1, pRefCell_), 0.0);

        // Set the diagonal to unity for cell pRefCell
        pEqn.diag()[pRefCell_] = -1.0;
    }

    // Insert the pressure equation
    foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
    (
        pEqn, jac, blockSize_ - 1, blockSize_ - 1, 1
    );

    // Calculate U-in-p equation coeffs coming from tr(grad(U)) == div(U)
    foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
    (
        p,
        U,
        jac,
        blockSize_ - 1,            // row offset
        0,                         // column offset
        fluidModel::twoD() ? 2 : 3 // number of scalar components of U
    );

    // Insert p-in-U term
    // Insert "-grad(p)" (equivalent to "-div(p*I)") into the U equation
    foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
    (
        p,
        jac,
        0,                         // row offset
        blockSize_ - 1,            // column offset
        fluidModel::twoD() ? 2 : 3 // number of scalar equations to insert
    );

    if (debug)
    {
        Info<< "End" << endl;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
