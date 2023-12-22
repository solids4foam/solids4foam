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

#include "unsLinGeomSolidMMS.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(unsLinGeomSolidMMS, 0);
addToRunTimeSelectionTable(solidModel, unsLinGeomSolidMMS, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsLinGeomSolidMMS::unsLinGeomSolidMMS
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    gradDf_
    (
        IOobject
        (
            "grad(" + D().name() + ")f",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_)
{
    DisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(D());

    // For consistent restarts, we will calculate the gradient field
    D().correctBoundaryConditions();
    D().storePrevIter();
    mechanical().interpolate(D(), pointD(), false);
    mechanical().grad(D(), pointD(), gradD(), gradDf_);

    // Store old times
    gradDf_.oldTime();
    sigmaf_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool unsLinGeomSolidMMS::evolve()
{
    Info << "Evolving solid solver" << endl;

    const scalar lambda_ = 1.1538e11;
    const scalar mu_ = 7.6923e10;
    const scalar ax_ = 2;
    const scalar ay_ = 4;
    const scalar az_ = 6;
    
	int iCorr = 0;
#ifdef OPENFOAMESIORFOUNDATION
	SolverPerformance<vector> solverPerfD;
	SolverPerformance<vector>::debug = 0;
#else
	lduSolverPerformance solverPerfD;
	blockLduMatrix::debug = 0;
#endif

	// Initialise the body force field to zero
	// This length of this list is equal to the number of cells (mesh.size) 
	DimensionedField<vector, volMesh> bodyForces
	(
	    IOobject
	    (
	        "bodyForces",
	        mesh().time().timeName(),
	        mesh(),
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    mesh(),
	    dimensionedVector("zero", dimForce/dimVolume, vector::zero)
	);
	
	// Body force field
	//volVectorField bodyForces(mesh().nCells(), vector::zero);
	
	// List of primary cells
	const cellList& cells = mesh().cells();
	
	// Cell-centre position vectors
	const vectorField& C = mesh().C().internalField();
	
	//pi
	const scalar pi = constant::mathematical::pi; 
 
	// List of primary points
	//const pointField& points = mesh().points();

	// List of primary faces
	//const faceList& faces = mesh().faces();
	
	// List of dual cells
	//const cellList& dualCells = dualMesh().cells();
	
	// List of dual points
	//const pointField& dualPoints = dualMesh().points();

	// List of dual faces
	//const faceList& dualFaces = dualMesh().faces();

	// List of dual cell volumes
	//const scalarField& dualV = dualMesh().V();
	
	// List of cell volumes
	const scalarField& V = mesh().V();  
	
	// Body forces field (list of values)
	vectorField& bodyForcesI = bodyForces.field();  
	
	forAll(bodyForcesI, cellI)
	{
		
		// Take a reference to the current primary volume
		const scalar curV = V[cellI]; 
		
		const scalar x = C[cellI].x();
		const scalar y = C[cellI].y();
		const scalar z = C[cellI].z();
		                                                                                                                                         
		//Set vector in dualCellI for x-equation     	
		bodyForcesI[cellI].x() =
		lambda_*(8*ay_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z) 
		+ 4*az_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y) 
		- 16*ax_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
		+ mu_*(8*ay_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z) 
		+ 4*az_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y) 
		- 5*ax_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
		- 32*ax_*mu_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z);
 

		//Set vector in dualCellI for y-equation  
		bodyForcesI[cellI].y() =
		lambda_*(8*ax_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z) 
		+ 2*az_*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x) 
		- 4*ay_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
		+ mu_*(8*ax_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z) 
		+ 2*az_*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x) 
		- 17*ay_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
		- 8*ay_*mu_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z);
 
 
		//Set vector in dualCellI for -equation  
		bodyForcesI[cellI].z() =
		lambda_*(4*ax_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y) 
		+ 2*ay_*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x) 
		- az_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
		+ mu_*(4*ax_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y) 
		+ 2*ay_*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x) 
		- 20*az_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
		- 2*az_*mu_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z);
	}

    Info<< "Solving the momentum equation for D" << endl;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Linear momentum equation total displacement form
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(mesh().Sf() & sigmaf_)
          + rho()*g()
          + bodyForces
        );

        // Under-relaxation the linear system
        DEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DEqn);

        // Hack to avoid expensive copy of residuals
#ifdef OPENFOAMESI
        const_cast<dictionary&>(mesh().solverPerformanceDict()).clear();
#endif

        // Solve the linear system
        solverPerfD = DEqn.solve();

        // Under-relax the field
        relaxField(D(), iCorr);

        // Update increment of displacement
        //DD() = D() - D().oldTime();

        // Interpolate D to pointD
        mechanical().interpolate(D(), pointD(), false);

        // Update gradient of displacement
        mechanical().grad(D(), pointD(), gradD(), gradDf_);

        // Update gradient of displacement increment
        //gradDD() = gradD() - gradD().oldTime();

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigmaf_);
        mechanical().correct(sigma());
    }
    while
    (
       !converged
        (
            iCorr,
#ifdef OPENFOAMESIORFOUNDATION
            mag(solverPerfD.initialResidual()),
            cmptMax(solverPerfD.nIterations()),
#else
            solverPerfD.initialResidual(),
            solverPerfD.nIterations(),
#endif
            D()
        ) && ++iCorr < nCorr()
    );

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


tmp<vectorField> unsLinGeomSolidMMS::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& gradD = gradDf_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& sigma = sigmaf_.boundaryField()[patchID];

    // Patch unit normals
    const vectorField n(patch.nf());

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & (sigma - impK*gradD))
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
