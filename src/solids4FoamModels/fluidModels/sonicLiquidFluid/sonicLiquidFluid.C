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

#include "sonicLiquidFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "adjustPhi.H"
#include "fvc.H"
#include "fvm.H"
#ifdef OPENFOAMESIORFOUNDATION
	#include "CorrectPhi.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sonicLiquidFluid, 0);
addToRunTimeSelectionTable(physicsModel, sonicLiquidFluid, fluid);
addToRunTimeSelectionTable(fluidModel, sonicLiquidFluid, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void sonicLiquidFluid::compressibleContinuityErrs()
{
	scalar sumLocalContErr =
        (sum(mag(rho_ - rho0_ - psi_*(p() - p0_)))/sum(rho_)).value();

    scalar globalContErr = 
		(
			sum(rho_ - rho0_ - psi_*(p() - p0_))/sum(rho_)
		).value();

    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors : sum local = " 
		<< sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_ << endl;
}


void sonicLiquidFluid::solveRhoEqn()
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho_)
      + fvc::div(phi())
#ifdef OPENFOAMESIORFOUNDATION
	  ==
		fvOptions_(rho_)
#endif
    );

#ifdef OPENFOAMESIORFOUNDATION
	fvOptions_.constrain(rhoEqn);

    rhoEqn.solve();

	fvOptions_.correct(rho_);
#else
    rhoEqn.solve();
#endif
}


#ifndef OPENFOAMESIORFOUNDATION
void sonicLiquidFluid::CorrectPhi()
{
	// Based this implementation in the version found in OFv1812 and 
	// the one from the pimpleFluid model 
	Info<< "Correcting flux for moving mesh" << endl;

	// Store div(rhoU) before update
	autoPtr<volScalarField> divrhoU;

	divrhoU.set
	(
		new volScalarField
		(
			"divrhoU",
			fvc::div(phi())
		)
	);

	volScalarField pcorr("pcorr", p());
	pcorr *= 0;

	// Initialise flux with interpolated velocity
	phi() = fvc::interpolate(rho_*U()) & mesh().Sf();

	// Not sure if needed...
//    adjustPhi(phi(), U(), pcorr);

	mesh().schemesDict().setFluxRequired(pcorr.name());

	surfaceScalarField rhorAUf
	(
		"rhorAUf",
		fvc::interpolate(rho_*rAU_)
	);

    while (pimple().correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
			fvm::ddt(psi_, pcorr)
		  + fvc::div(phi())
          - fvm::laplacian(rhorAUf, pcorr)
		 == divrhoU()
		);

        pcorrEqn.solve();

        if (pimple().finalNonOrthogonalIter())
        {
            phi() -= pcorrEqn.flux();
        }
    }
}
#endif

void sonicLiquidFluid::compressibleCourantNo()
{
	scalar CoNum = 0.0;                                                                               
	scalar meanCoNum = 0.0;                                                     
																				
	surfaceScalarField SfUfbyDelta =
		mesh().surfaceInterpolation::deltaCoeffs()
	  * mag(phi())
	  / fvc::interpolate(rho_);

	CoNum = max(SfUfbyDelta/mesh().magSf())
		.value()*runTime().deltaT().value();

	meanCoNum = (sum(SfUfbyDelta)/sum(mesh().magSf()))
		.value()*runTime().deltaT().value();

	Info<< "Region: " << mesh().name()
		<< " Courant Number mean: " << meanCoNum
		<< " max: " << CoNum << endl;
}


void sonicLiquidFluid::solvePEqn
(
	const fvVectorMatrix& UEqn
)
{
    rAU_ = 1.0/UEqn.A();

	surfaceScalarField rhorAUf
	(
		"rhorAUf",
		fvc::interpolate(rho_*rAU_)
	);

	U() = rAU_*UEqn.H();

	surfaceScalarField phid
	(
		"phid",
		psi_
	   *(
#ifdef OPENFOAMESIORFOUNDATION
			fvc::flux(U())
		  + rhorAUf*fvc::ddtCorr(rho_, U(), phi())/fvc::interpolate(rho_)
#else
			(fvc::interpolate(U()) & mesh().Sf())
		  + fvc::ddtPhiCorr(rAU_, rho_, U(), phi())
#endif
		)
	);

	phi() = fvc::interpolate(rhoO_)*phid/psi_;

	fvScalarMatrix pEqn
	(
		fvm::ddt(psi_, p())
	  + fvc::div(phi())
	  + fvm::div(phid, p())
	  - fvm::laplacian(rhorAUf, p())
	  + fvc::ddt(rhoO_) // Essential to account for this term with mesh motion
	);

	pEqn.solve();

	phi() += pEqn.flux();

	solveRhoEqn();

	// State equation errors
	compressibleContinuityErrs();
	
	// Correct velocity
	U() -= rAU_*fvc::grad(p());
	U().correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sonicLiquidFluid::sonicLiquidFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
	mu_(transportProperties_.lookup("mu")),
	K_(transportProperties_.lookup("K")),
    thermodynamicProperties_
    (
        IOobject
        (
            "thermodynamicProperties",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
	rho0_(thermodynamicProperties_.lookup("rho0")),
	p0_(thermodynamicProperties_.lookup("p0")),
    psi_(thermodynamicProperties_.lookup("psi")),
	rhoO_
	(
	 	IOobject
		(
		 	"rhoO",
			runTime.timeName(),
			mesh(),
			IOobject::NO_READ,
            IOobject::NO_WRITE
		),
		mesh(),
		rho0_ - psi_*p0_
	),
	rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoO_ + psi_*p()
    ),
#ifdef OPENFOAMESIORFOUNDATION
	fvOptions_(fv::options::New(mesh())),
#endif
    rAU_
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        runTime.deltaT()/rho0_,
        zeroGradientFvPatchScalarField::typeName
    ),
    cumulativeContErr_(0)
{
    UisRequired();
    pisRequired();

    // Reset phi dimensions: compressible
    Info<< "Resetting the dimensions of phi" << endl;
    phi().dimensions().reset(dimVelocity*dimArea*dimDensity);

    phi() = linearInterpolate(rho_*U()) & mesh().Sf();

#ifdef OPENFOAMESIORFOUNDATION
	mesh().setFluxRequired(p().name());

	if (mesh().dynamic())
	{
		Info<< "Constructing face momentum rhoUf" << endl;

		rhoUf_.set
		(
			new surfaceVectorField
			(
				IOobject
				(
					"rhoUf",
					runTime.timeName(),
					mesh(),
					IOobject::READ_IF_PRESENT,
					IOobject::AUTO_WRITE
				),
				fvc::interpolate(rho_*U())
			)
		);
	}

	// Check if any finite volume option is present
	if (!fvOptions_.optionList::size())
	{
		Info << "No finite volume options present\n" << endl;
	}
#else
    mesh().schemesDict().setFluxRequired(p().name());
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> sonicLiquidFluid::patchViscousForce
(
    const label patchID
) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

#ifdef OPENFOAMESIORFOUNDATION
	tvF.ref() =
#else
	tvF() =
#endif
        mu_.value()*U().boundaryField()[patchID].snGrad();
        
    return tvF;
}


tmp<scalarField> sonicLiquidFluid::patchPressureForce
(
    const label patchID
) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

	// Pressure here is already in Pa
#ifdef OPENFOAMESIORFOUNDATION
	tpF.ref() =
#else
	tpF() =
#endif
   	p().boundaryField()[patchID];

    return tpF;
}


bool sonicLiquidFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    dynamicFvMesh& mesh = this->mesh();

	bool meshChanged = false;

	if (fluidModel::fsiMeshUpdate())
    {
        // The FSI interface is in charge of calling mesh.update()
        meshChanged = fluidModel::fsiMeshUpdateChanged();
    }
    else
    {
        meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());
    }

    if (meshChanged)
    {
        const Time& runTime = fluidModel::runTime();
#       include "volContinuity.H"
    }

    bool correctPhi
    (
        pimple().dict().lookupOrDefault("correctPhi", false)
	);

	if (correctPhi && meshChanged)
    {
#ifdef OPENFOAMESIORFOUNDATION
		// Store divrhoU from the previous mesh so that it can be mapped
		// and used in correctPhi to ensure the corrected phi has the
		// same divergence
		autoPtr<volScalarField> divrhoU;
		if (correctPhi)
		{
			divrhoU.set
			(
				new volScalarField
				(
					"divrhoU",
					fvc::div(phi())
				)
			);
		}

		// Store momentum to set rhoUf for introduced faces.
		autoPtr<volVectorField> rhoU;
		//        if (rhoUf_.valid())
		//        {
		//            rhoU.set(new volVectorField("rhoU", rho_*U()));
		//        }

		// Calculate absolute flux
		// from the mapped surface velocity
		phi() = mesh.Sf() & rhoUf_();

		// Define volScalarField to hold phi
		// to pass to compressible CorrectPhi function
		const volScalarField psi
		(
			 IOobject
			 (
			  	"psi",
				runTime().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			 ),
			 mesh,
			 psi_
		);

		CorrectPhi
		(
			U(),
			phi(),
			p(),
			rho_,
			psi,
			dimensionedScalar("rAUf", dimTime, 1),
			divrhoU(),
			pimple()
#ifndef OPENFOAMESI
			,
			true
#endif
		);
#else
		   CorrectPhi();
#endif
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), rho_, U());

    // Calculate CourantNo
	compressibleCourantNo();

    // Pressure-velocity corrector
    while (pimple().loop())
    {
#ifdef OPENFOAMFOUNDATION
		if (!pimple().simpleRho())
		{
			solveRhoEqn();
		}
#else
		solveRhoEqn();
#endif

        // --- Momentum equation
        // Time-derivative matrix
        fvVectorMatrix ddtUEqn(fvm::ddt(rho_, U()));

        // Convection-diffusion matrix
        fvVectorMatrix HUEqn
		(
		    fvm::div(phi(), U())
		  - fvm::laplacian(mu_, U())
#ifdef OPENFOAMESIORFOUNDATION
		 ==
			fvOptions_(rho_, U())
#endif
		);

#ifdef OPENFOAMESIORFOUNDATION
		tmp<fvVectorMatrix> tUEqn(ddtUEqn + HUEqn);
		fvVectorMatrix& UEqn = tUEqn.ref();

		UEqn.relax();

		fvOptions_.constrain(UEqn);

		if (pimple().momentumPredictor())
		{
			solve(UEqn == -fvc::grad(p()));
		}

        fvOptions_.correct(U());
#else
        fvVectorMatrix UEqn(ddtUEqn + HUEqn);

		// Get under-relaxation factor
		scalar UUrf =
			mesh.solutionDict().equationRelaxationFactor
			(
				U().select(pimple().finalIter())
			);

		if (pimple().momentumPredictor())
		{
			solve
			(
				ddtUEqn
			  + relax(HUEqn, UUrf)
			 ==
			  - fvc::grad(p()),
				mesh.solutionDict().solver((U().select(pimple().finalIter())))
			);
		}
		else
        {
            // Explicit update
            U() = (ddtUEqn.H() + HUEqn.H() - fvc::grad(p()))/(HUEqn.A() + ddtUEqn.A());
            U().correctBoundaryConditions();
        }
#endif

        // --- Pressure corrector loop
        while (pimple().correct())
        {
#ifdef OPENFOAMESIORFOUNDATION
			solvePEqn(tUEqn.ref());
#else
			solvePEqn(UEqn);
#endif
        }

#ifdef OPENFOAMESIORFOUNDATION
		tUEqn.clear();
#endif

        gradU() = fvc::grad(U());

    }

	// Update rho based on equation of state
	rho_ = rhoO_ + psi_*p();

    // Make the fluxes absolute to the mesh motion
    fvc::makeAbsolute(phi(), rho_, U());

	// Print variables for inspection
	// Density variation
	scalar deltaRho = max(rho_).value() - min(rho_).value();
	scalar refDeltaRho = 0.01*rho0_.value();

	// Info variable values
	Info<< nl << "Density: min " << min(rho_).value() 
		<< " max " << max(rho_).value() << endl;

	Info<< "Density variation: " << deltaRho 
		<< " ref: " << refDeltaRho << endl; 

	Info<< "Pressure: min " << min(p()).value()
		<< " max " << max(p()).value() << endl;

	Info<< "Velocity: min " << min(mag(U())).value()
		<< " max " << max(mag(U())).value() << nl << nl;

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
