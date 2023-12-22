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

#include "linGeomSolidMMS.H"
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

defineTypeNameAndDebug(linGeomSolidMMS, 0);
addToRunTimeSelectionTable(solidModel, linGeomSolidMMS, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomSolidMMS::linGeomSolidMMS
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    rhoDdtD_0_
    (
        IOobject
        (
            "rhoDdtD_0",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    )
{
    DDisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(DD());
    fvc::d2dt2(D().oldTime());

    // For consistent restarts, we will calculate the gradient field
    DD().correctBoundaryConditions();
    DD().storePrevIter();
    mechanical().grad(DD(), gradDD());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool linGeomSolidMMS::evolve()
{
    Info<< "Evolving solid solver" << endl;
    
    const scalar lambda_ = 1.1538e11;
    const scalar mu_ = 7.6923e10;
    const scalar ax_ = 2;
    const scalar ay_ = 4;
    const scalar az_ = 6;

    // Mesh update loop
    do
    {
        int iCorr = 0;
#ifdef OPENFOAMESIORFOUNDATION
        SolverPerformance<vector> solverPerfDD;
        SolverPerformance<vector>::debug = 0;
#else
        lduSolverPerformance solverPerfDD;
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

	//	    // Take a reference to the current primary volume
	//	    const scalar curDualV = dualV[vertexI];
	//	    
	//	    //Take a reference to the current dual cell
	//	    const cell& curDualCell = dualCells[vertexI];
	//	    
	//	    //Create a bounding box around the current dual cell
	//	    const boundBox bb(curDualCell.points(dualFaces, dualPoints), false);
	//	    
	//	    //Get the most negative x and y coordinates of the dual cell
	//	    const scalar x1 = bb.min().x();
	//	    const scalar y1 = bb.min().y();
	//	    const scalar z1 = bb.min().z();
	//	   
	//	    //Get the most negative x and y coordinates of the dual cell
	//	    const scalar x2 = bb.max().x();
	//	    const scalar y2 = bb.max().y();
	//	    const scalar z2 = bb.max().z();
			
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

        Info<< "Solving the momentum equation for DD" << endl;

        // Momentum equation loop
        do
        {
            // Store fields for under-relaxation and residual calculation
            DD().storePrevIter();

            // Linear momentum equation total displacement form
            fvVectorMatrix DDEqn
            (
                rho()*fvm::d2dt2(DD())
              + rhoDdtD_0_
             == fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
              - fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
              + fvc::div(sigma(), "div(sigma)")
              + rho()*g()
              + stabilisation().stabilisation(DD(), gradDD(), impK_)
              + bodyForces
            );

            // Under-relaxation the linear system
            DDEqn.relax();

            // Enforce any cell displacements
            solidModel::setCellDisps(DDEqn);

            // Hack to avoid expensive copy of residuals
#ifdef OPENFOAMESI
            const_cast<dictionary&>(mesh().solverPerformanceDict()).clear();
#endif

            // Solve the linear system
            solverPerfDD = DDEqn.solve();

            // Fixed or adaptive field under-relaxation
            relaxField(DD(), iCorr);

            // Update the total displacement
            D() = D().oldTime() + DD();

            // Update gradient of displacement increment
            mechanical().grad(DD(), gradDD());

            // Update gradient of total displacement
            gradD() = gradD().oldTime() + gradDD();

            // Calculate the stress using run-time selectable mechanical law
            const volScalarField DDEqnA("DEqnA", DDEqn.A());
            mechanical().correct(sigma());
        }
        while
        (
            !converged
            (
                iCorr,
#ifdef OPENFOAMESIORFOUNDATION
                mag(solverPerfDD.initialResidual()),
                cmptMax(solverPerfDD.nIterations()),
#else
                solverPerfDD.initialResidual(),
                solverPerfDD.nIterations(),
#endif
                DD()
            ) && ++iCorr < nCorr()
        );

        // Update point displacement increment
        mechanical().interpolate(DD(), pointDD());

        // Update point displacement
        pointD() = pointD().oldTime() + pointDD();

        // Update velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());

    // Store ddt old term
    rhoDdtD_0_ = rho()*fvc::d2dt2(D());

#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


tmp<vectorField> linGeomSolidMMS::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& pImpK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& pRImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradDD = gradDD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField n(patch.nf());

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & (pSigma - pImpK*pGradDD))
            )*pRImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
